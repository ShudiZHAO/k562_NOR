#!/usr/bin/env Rscript
# =============================================================================
# individual_unit_methylation.R
#
# R 版本的 individual_unit_methylation.py
# 每个NOR一行：同一NOR的所有units画在同一条x轴上，x轴数值以该NOR第一个unit的
# 起点为0（相对坐标，单位kb），相邻unit之间用竖虚线分隔；行左侧水平标注"NOR X"，
# 各行标签左对齐（共用同一左缘）。
# 所有行共用同一x轴scale（bp/英寸 = 参考NOR14的units跨度 / 6.72英寸）：每行的
# x轴长度按该NOR自身units跨度伸展（左缘对齐、右端随数据不同），页宽由最长的行
# 决定（不限宽），所有unit完整显示（不裁剪）。
# 整份PDF只有一页（所有NOR行画在同一页上）；无边框、无网格线、无panel标题，
# 仅保留左/下坐标轴。
# 柱状图为不透明实色（不做透明度/alpha处理）；x轴标题为"Consensus position (kb)"，
# 在整页底部居中标注，x轴每10kb一个label；y轴标题为"Methylation (%)"，在整页
# 左侧外边距中竖排、竖直居中标注（整页所有NOR共用该y轴标题，各行的y轴刻度仍逐行自适应）。
#
# 用法: 四个输入文件已内嵌默认路径（相对运行时的工作目录），无参直接运行即可:
#   Rscript individual_unit_methylation.R
# 默认值:
#   --unit_type_file  ../Figure3_Variations/3C_unit_type_map/unit_to_type.xls
#   --array_nor_file  ../name_of_the_arrays.txt
#   --bed_file        ../Figure2_Assembly/2A_circle_profile/named_units.bed
#   --bedgraph_file   stats_Me.bedgraph
# 需要时可在命令行覆盖:
#   Rscript individual_unit_methylation.R \
#     [--unit_type_file unit_to_type.xls] [--array_nor_file name_of_the_arrays.txt] \
#     [--bed_file named_units.bed] [--bedgraph_file stats_Me.bedgraph] \
#     [--bin_size 100] [--output individual_unit_methylation] [--help]
#
# 依赖: data.table
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# 定义type颜色（根据数据中有5种type，定义5种颜色）
TYPE_COLORS <- c(
  "1" = "#5bb8ac",  # Type 1
  # "2" = "#e5e5e5",  # Type 2
  "2" = "#a5a5a5",  # Type 2
  "3" = "#f0c238",  # Type 3
  "4" = "#762A83",  # Type 4
  "5" = "#f07979"   # Type 5
)
DEFAULT_COLOR <- "#808080"  # 未定义type的默认颜色

# ---- 布局常量（单位：英寸）：每个NOR一行，整页一列 ----
ROW_PLOT_H_IN  <- 0.5    # 每行（每个NOR）绘图区高度
ROW_MAI  <- c(0.16, 0.72, 0.06, 0.02)   # 每行 下/左/上/右 边距（英寸，mai；左侧加宽以容纳水平的"NOR X"行标签）
ROW_LABEL_PAD_LEFT_IN <- 0.05           # "NOR X"行标签左对齐时距内区域左缘的空隙（英寸）
# 注意：必须用omi（英寸）而不是oma（oma的单位是文本行数，传英寸值会被误缩小）
PAGE_OMI <- c(0.38, 0.30, 0.10, 0.02)   # 整页 下/左/上/右 外边距（英寸，omi；底部放居中的x轴标题，左侧放竖排的y轴标题）
# ---- 公共x轴scale：所有行与NOR14一致（相同bp/英寸），行长度随各自数据伸展 ----
REF_NOR_NUM       <- 14L    # 参考NOR编号
# NOR14的绘图区宽度（英寸，固定值6.72）：公共scale = NOR14的units跨度 / 该值。
# 6.72 = 参考版面宽度7.5 - 外边距(0.02+0.02) - 行边距(0.72+0.02)，即原脚本的实际
# 绘图区宽度；固定为常数可保证页面外边距变化（如左侧新增竖排y轴标题）时x轴scale
# （bp/英寸）保持与原来完全一致。
REF_PLOT_WIDTH_IN <- 6.72

# ---- 样式常量 ----
LINE_LWD       <- 0.25        # 所有线宽（pt）：轴线/刻度线/unit分隔线
TICK_LEN_PT    <- 1.2         # 刻度线长度（pt，向外）
TICK_LAB_CEX   <- 5 / 8       # 刻度标签5pt（pointsize=8时的cex）
AXIS_TITLE_CEX <- 1           # 轴标题与行标签8pt
UNIT_SEP_COL   <- gray(0.4)   # unit分隔虚线颜色
UNIT_SEP_LTY   <- 2           # unit分隔虚线线型
X_TICK_BY      <- 10000       # x轴刻度间隔（bp）
X_LABEL_BY     <- 10000       # x轴标签间隔（bp，标签数值除以1000以kb显示；10000=每10kb一个label）
Y_TICK_BY      <- 25          # y轴刻度间隔
Y_LABEL_BY     <- 50          # y轴标签间隔

# ---------------------------------------------------------------- 参数解析 --
# 支持 "--flag value" 与 "--flag=value" 两种写法
get_arg <- function(args, flag, default = NULL, required = FALSE, convert = as.character) {
  idx <- match(flag, args)
  value <- NULL
  if (!is.na(idx)) {
    if (idx == length(args)) {
      stop(sprintf("Argument %s needs a value", flag), call. = FALSE)
    }
    value <- args[idx + 1L]
  } else {
    hit <- grep(paste0("^", flag, "="), args)
    if (length(hit) > 0) value <- sub(paste0("^", flag, "="), "", args[hit[1L]])
  }
  if (is.null(value)) {
    if (required) {
      stop(sprintf("Missing required argument: %s", flag), call. = FALSE)
    }
    return(default)
  }
  convert(value)
}

show_help <- function() {
  cat("Plot individual unit methylation distribution by type\n\n")
  cat("Usage:\n")
  cat("  Rscript individual_unit_methylation.R\n")
  cat("    [--unit_type_file  <unit_to_type.xls>]\n")
  cat("        default: ../Figure3_Variations/3C_unit_type_map/unit_to_type.xls\n")
  cat("    [--array_nor_file  <name_of_the_arrays.txt>]\n")
  cat("        default: ../name_of_the_arrays.txt\n")
  cat("    [--bed_file        <named_units.bed>]\n")
  cat("        default: ../Figure2_Assembly/2A_circle_profile/named_units.bed\n")
  cat("    [--bedgraph_file   <stats_Me.bedgraph>]\n")
  cat("        default: stats_Me.bedgraph\n")
  cat("    [--bin_size       <bp, default: 100>]\n")
  cat("    [--output         <file prefix, default: individual_unit_methylation>]\n")
  cat("    [--help]\n")
}

# ------------------------------------------------------------ 数据加载函数 --
load_unit_type_data <- function(unit_type_file) {
  # 加载unit的type和subtype数据
  df <- read.delim(unit_type_file, stringsAsFactors = FALSE, check.names = FALSE)
  # 去掉由行尾多余制表符产生的空列（pandas中的 Unnamed 列，不参与计算）
  keep <- nzchar(names(df)) & !vapply(df, function(x) all(is.na(x)), logical(1))
  df[, keep, drop = FALSE]
}

load_array_to_nor_mapping <- function(array_nor_file) {
  # 加载Array到NOR的映射
  lines <- readLines(array_nor_file, warn = FALSE)
  array_names <- character(0)
  nor_nums <- integer(0)
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line)) next
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) >= 3) {
      array_name <- parts[2]  # Array_X
      nor_name <- parts[3]    # NOR X
      if (nor_name == "Dropped") next
      # 提取NOR编号
      nor_num <- suppressWarnings(
        as.integer(tail(strsplit(trimws(nor_name), "\\s+")[[1]], 1))
      )
      if (is.na(nor_num)) {
        warning(sprintf("Cannot parse NOR number from '%s', skip line: %s",
                        nor_name, line))
        next
      }
      array_names <- c(array_names, array_name)
      nor_nums <- c(nor_nums, nor_num)
    }
  }
  stats::setNames(nor_nums, array_names)
}

load_unit_positions <- function(bed_file) {
  # 加载unit位置信息
  df <- read.delim(bed_file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) < 5) {
    stop("BED file must contain at least 5 columns", call. = FALSE)
  }
  names(df)[1:5] <- c("chrom", "start", "end", "name", "score")
  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  # 过滤：保留开头为"Array"但排除开头为"Array_1_"的行
  df[grepl("^Array", df$name) & !grepl("^Array_1_", df$name), , drop = FALSE]
}

load_bedgraph <- function(bedgraph_file) {
  # 加载bedgraph甲基化数据（chrom, start, end, value）
  # 注意用file=显式传路径，避免fread把含空格的路径当作系统命令
  bg <- fread(file = bedgraph_file, header = FALSE, sep = "\t",
              showProgress = FALSE,
              colClasses = list(character = "V1", integer = c("V2", "V3"),
                                numeric = "V4"))
  names(bg) <- c("chrom", "bg_start", "bg_end", "value")
  # 跳过空行与注释行（# / track / browser）
  bg <- bg[!is.na(bg_start) & !is.na(bg_end) & !is.na(value) &
             nzchar(chrom) & !grepl("^(#|track|browser)", chrom)]
  setkey(bg, chrom, bg_start)
  bg[]
}

# ------------------------------------------------------------ 核心计算函数 --
parse_unit_id <- function(unit_id) {
  # 解析unit ID中的x和y值（如 Array_3_2 -> c(3, 2)）
  m <- regmatches(unit_id, regexec("^Array_(\\d+)_(\\d+)", unit_id))[[1]]
  if (length(m) == 3L) as.integer(m[2:3]) else c(0L, 0L)
}

parse_new_type <- function(new_type_str) {
  # 从new_type字符串中解析type和subtype（如 "1-1" -> c(1, 1)）
  if (is.null(new_type_str)) return(c(0L, 0L))
  if (length(new_type_str) == 0L || is.na(new_type_str)) return(c(0L, 0L))
  parts <- strsplit(trimws(as.character(new_type_str)), "-", fixed = TRUE)[[1]]
  if (length(parts) >= 2) {
    t1 <- suppressWarnings(as.integer(trimws(parts[1])))
    t2 <- suppressWarnings(as.integer(trimws(parts[2])))
    if (!is.na(t1) && !is.na(t2)) return(c(t1, t2))
  }
  c(0L, 0L)
}

extract_methylation_for_region <- function(bg, chr, start, end) {
  # 从bedgraph数据中提取指定区域的甲基化数据，按bp展开为逐bp向量
  # （等价于Python版本：每个interval按重叠长度重复其value）
  # 注意：参数用chr命名，避免与bg的chrom列名冲突被data.table误解析
  if (is.na(chr) || length(chr) == 0L || !chr %in% bg$chrom) {
    return(numeric(0))
  }
  sub <- bg[.(chr), on = "chrom"]                # 该chrom的所有行（已按bg_start排序）
  sub <- sub[bg_end > start & bg_start < end]    # 半开区间重叠判定，与Python一致
  if (nrow(sub) == 0) return(numeric(0))

  ov_start <- pmax(start, sub$bg_start)
  ov_end <- pmin(end, sub$bg_end)
  ov_len <- ov_end - ov_start
  keep <- ov_len > 0
  if (!any(keep)) return(numeric(0))

  # 每个interval按重叠bp数展开（与Python的 values.extend([bg_value]*n) 一致）
  rep(sub$value[keep], times = ov_len[keep])
}

calculate_unit_methylation_profile <- function(chrom, start, end, bg, bin_size) {
  # 计算单个unit的甲基化分布
  unit_length <- end - start

  # 按bp展开该区域的甲基化值（第j个元素对应unit内偏移 j-1 bp）
  meth_values <- extract_methylation_for_region(bg, chrom, start, end)

  # 如果没有甲基化数据，返回空向量
  if (length(meth_values) == 0) {
    return(list(binned = numeric(0), unit_length = unit_length))
  }

  # 按bin_size进行bin化：bin i 覆盖偏移 [i*bin_size, (i+1)*bin_size)
  num_bins <- ceiling(unit_length / bin_size)
  bin_idx <- (seq_along(meth_values) - 1L) %/% bin_size   # 0-based bin编号

  binned <- rep(NA_real_, num_bins)
  sums <- rowsum(meth_values, group = bin_idx, reorder = TRUE)
  counts <- tabulate(bin_idx + 1L, nbins = num_bins)      # bin 0 -> 位置1
  bin_ids <- as.integer(rownames(sums))
  binned[bin_ids + 1L] <- sums[, 1] / counts[bin_ids + 1L]

  list(binned = binned, unit_length = unit_length)
}

# ------------------------------------------------------------ 绘图函数 --
# 整页高度（英寸）：每个NOR一行
nor_page_height <- function(num_nor) {
  num_nor * (ROW_PLOT_H_IN + ROW_MAI[1L] + ROW_MAI[3L]) +
    PAGE_OMI[1L] + PAGE_OMI[3L]
}

# 参考NOR（默认NOR14）的units总跨度（最右end - 最左start）：用于定义所有行共用
# 的x轴scale；找不到参考NOR或跨度无效时，退化为所有NOR中的最大跨度
reference_x_span <- function(nor_data, ref_nor = REF_NOR_NUM) {
  span_of <- function(nd) {
    if (length(nd$units) == 0L) return(NA_real_)
    max(vapply(nd$units, function(x) x$end, numeric(1))) -
      min(vapply(nd$units, function(x) x$start, numeric(1)))
  }
  hit <- Filter(function(nd) identical(nd$nor_num, as.integer(ref_nor)), nor_data)
  if (length(hit) > 0L) {
    s <- span_of(hit[[1L]])
    if (is.finite(s) && s > 0) return(s)
  }
  spans <- vapply(nor_data, span_of, numeric(1))
  spans <- spans[is.finite(spans) & spans > 0]
  if (length(spans) > 0L) max(spans) else 0
}

# 公共x轴scale（bp/英寸）：参考NOR的units跨度 / 其绘图区宽度，所有行一致
common_scale_bp_per_in <- function(nor_data) {
  ref_span <- reference_x_span(nor_data, REF_NOR_NUM)
  if (is.finite(ref_span) && ref_span > 0) {
    ref_span / REF_PLOT_WIDTH_IN
  } else {
    NA_real_
  }
}

# 整页宽度（英寸）：外边距 + 行左/右边距 + 最长行的绘图宽（页宽随最长NOR自适应，
# 不限宽，保证长NOR的所有unit完整显示、各行scale一致）
nor_page_width <- function(nor_data, scale_bp_per_in) {
  span_of <- function(nd) {
    if (length(nd$units) == 0L) return(0)
    max(vapply(nd$units, function(x) x$end, numeric(1))) -
      min(vapply(nd$units, function(x) x$start, numeric(1)))
  }
  max_span <- max(vapply(nor_data, span_of, numeric(1)))
  PAGE_OMI[2L] + PAGE_OMI[4L] + ROW_MAI[2L] + ROW_MAI[4L] +
    max_span / scale_bp_per_in
}

draw_nor_row <- function(units_data, nor_num, bin_size) {
  # 单个NOR一行的面板：该NOR所有unit的甲基化柱画在同一条x轴上，x轴原点(0)取该NOR
  # 第一个unit的起点（相对坐标，刻度数值以kb显示），相邻unit之间用竖虚线分隔；
  # 只保留左/下坐标轴（无边框、无网格、无标题）
  # 注意：绘图帧(plot.new)由draw_all_nor_rows统一创建，本函数不再开帧——
  # 否则每行多翻一帧，PDF会在开头多出一页空白页
  # x轴范围按该NOR自身数据自适应：0=第一个unit的起点，右端=最后一个unit的终点，
  # 所有unit完整显示（不裁剪、各行宽度可以不同），数据结束处即轴结束
  xmin <- min(vapply(units_data, function(x) x$start, numeric(1)))
  xmax <- max(vapply(units_data, function(x) x$end, numeric(1)))

  # y轴范围（该NOR所有unit合并后自动缩放：柱底粘性为0，顶部留5%边距）
  all_vals <- unlist(lapply(units_data, `[[`, "binned"), use.names = FALSE)
  all_vals <- all_vals[!is.na(all_vals)]
  ymax <- if (length(all_vals) > 0) max(all_vals) else 0
  ylim <- if (!is.finite(ymax) || ymax <= 0) c(0, 1) else c(0, ymax * 1.05)

  plot.window(xlim = c(0, xmax - xmin), ylim = ylim, xaxs = "i", yaxs = "i")

  # pt -> 用户坐标换算（绘图区实际尺寸取自par("pin")，单位英寸）
  u <- par("usr")
  pin <- par("pin")
  xu <- function(dx_pt) dx_pt / 72 / pin[1] * (u[2] - u[1])
  yu <- function(dy_pt) dy_pt / 72 / pin[2] * (u[4] - u[3])

  # 柱状图：按unit在基因组上的实际位置摆放（宽度bin_size*0.8，无边框；
  # yaxs="i"使柱底贴住下边框）；颜色为不透明实色（不做alpha/透明度处理）
  for (ud in units_data) {
    binned <- ud$binned
    if (length(binned) == 0L) next
    positions <- (ud$start - xmin) + (seq_along(binned) - 1L) * bin_size + bin_size / 2
    valid <- !is.na(binned)
    if (!any(valid)) next
    rect(positions[valid] - bin_size * 0.4, 0,
         positions[valid] + bin_size * 0.4, binned[valid],
         col = ud$color, border = NA)
  }

  # unit分隔线：相邻unit边界处的竖虚线（unit相连时即共享边界；有间隔时取中点）
  if (length(units_data) > 1L) {
    ends <- vapply(units_data, function(x) x$end, numeric(1))
    starts <- vapply(units_data, function(x) x$start, numeric(1))
    seps <- unique((head(ends, -1L) + tail(starts, -1L)) / 2) - xmin
    if (length(seps) > 0L) {
      segments(seps, u[3], seps, u[4],
               col = UNIT_SEP_COL, lty = UNIT_SEP_LTY, lwd = LINE_LWD)
    }
  }

  # 坐标轴：只画左轴与下轴的轴线（0.25pt），不加边框、不画网格
  segments(u[1], u[3], u[1], u[4], lwd = LINE_LWD)  # y轴线
  segments(u[1], u[3], u[2], u[3], lwd = LINE_LWD)  # x轴线

  # 刻度位置：x轴数值从0开始（相对该NOR第一个unit的起点），每10000bp一个刻度、
  # 每50000bp一个标签（数值以kb显示）；y轴每25一个刻度、每50一个标签
  xticks <- seq(from = 0,
                to = floor((xmax - xmin) / X_TICK_BY) * X_TICK_BY, by = X_TICK_BY)
  xlab_ticks <- xticks[xticks %% X_LABEL_BY == 0L]
  yticks <- seq(from = 0, to = floor(ylim[2] / Y_TICK_BY) * Y_TICK_BY, by = Y_TICK_BY)
  ylab_ticks <- yticks[yticks %% Y_LABEL_BY == 0L]

  # 刻度线：向外1.2pt、线宽0.25pt
  if (length(xticks) > 0L) {
    segments(xticks, u[3], xticks, u[3] - yu(TICK_LEN_PT),
             lwd = LINE_LWD, xpd = NA)
  }
  if (length(yticks) > 0L) {
    segments(u[1] - xu(TICK_LEN_PT), yticks, u[1], yticks,
             lwd = LINE_LWD, xpd = NA)
  }

  # 刻度标签（5pt）：贴近坐标轴——y轴右对齐于轴线左侧4pt（每50一个标签）；
  # x轴居中于轴线下方6pt，缩小标签与轴线的距离
  if (length(ylab_ticks) > 0L) {
    text(u[1] - xu(4), ylab_ticks, as.character(ylab_ticks), adj = c(1, 0.5),
         cex = TICK_LAB_CEX, xpd = NA)
  }
  if (length(xlab_ticks) > 0L) {
    # x轴标签：bp数值除以1000，以kb显示（如150000 -> 150）
    text(xlab_ticks, u[3] - yu(6), as.character(xlab_ticks / 1000),
         adj = c(0.5, 0.5), cex = TICK_LAB_CEX, xpd = NA)
  }

  # 行标签：该行左侧水平标注"NOR X"（8pt），各行共用同一左缘（左对齐）——锚点取
  # 内区域左缘右侧ROW_LABEL_PAD_LEFT_IN英寸处（各行左边距几何一致，故左缘严格
  # 对齐），位于y轴刻度标签左侧，代替原来的ylab与panel标题
  text(u[1] - xu(72 * (ROW_MAI[2L] - ROW_LABEL_PAD_LEFT_IN)), mean(u[3:4]),
       sprintf("NOR %d", nor_num), adj = c(0, 0.5),
       cex = AXIS_TITLE_CEX, xpd = NA)
}

prepare_nor_data <- function(unit_type_df, array_to_nor, unit_positions,
                             bg, bin_size) {
  # 添加NOR列到unit_type_df（由origin_name解析Array编号后查映射，找不到记0）
  parsed <- t(vapply(unit_type_df$origin_name, parse_unit_id, integer(2)))
  arr_num <- parsed[, 1]
  unit_type_df$NOR <- as.integer(array_to_nor[paste0("Array_", arr_num)])
  unit_type_df$NOR[is.na(unit_type_df$NOR)] <- 0L

  # 解析type和subtype：优先用new_type列（如"1-1"）；无该列时直接用type/subtype列
  if ("new_type" %in% names(unit_type_df)) {
    ts <- t(vapply(unit_type_df$new_type, parse_new_type, integer(2)))
    unit_type_df$type_from_new <- ts[, 1]
    unit_type_df$subtype_from_new <- ts[, 2]
  } else {
    unit_type_df$type_from_new <-
      suppressWarnings(as.integer(unit_type_df$type))
    unit_type_df$subtype_from_new <-
      suppressWarnings(as.integer(unit_type_df$subtype))
    unit_type_df$type_from_new[is.na(unit_type_df$type_from_new)] <- 0L
    unit_type_df$subtype_from_new[is.na(unit_type_df$subtype_from_new)] <- 0L
  }

  # 移除NOR为0或type为0的行
  unit_type_df <- unit_type_df[
    unit_type_df$NOR != 0 & unit_type_df$type_from_new != 0, , drop = FALSE]

  # 按NOR排序（1-18），同NOR内按index排序
  unit_type_df <- unit_type_df[order(unit_type_df$NOR, unit_type_df$index), , drop = FALSE]
  rownames(unit_type_df) <- NULL

  if (nrow(unit_type_df) == 0L) return(list())

  # 按NOR分组（每个NOR一行），逐unit收集绘图数据
  nor_groups <- split(unit_type_df, unit_type_df$NOR)
  nor_nums <- as.integer(names(nor_groups))
  lapply(seq_along(nor_groups), function(i) {
    rows <- nor_groups[[i]]
    units_data <- vector("list", nrow(rows))
    n_kept <- 0L
    for (k in seq_len(nrow(rows))) {
      unit_id <- rows$origin_name[k]
      idx <- which(unit_positions$name == unit_id)[1L]
      if (is.na(idx)) {
        message(sprintf("Warning: Unit %s not found in BED file", unit_id))
        next
      }
      unit_row <- unit_positions[idx, ]

      # 获取颜色（柱状图为不透明实色，不再按subtype设置透明度）
      unit_type <- rows$type_from_new[k]
      color <- if (as.character(unit_type) %in% names(TYPE_COLORS)) {
        unname(TYPE_COLORS[[as.character(unit_type)]])
      } else {
        DEFAULT_COLOR  # 默认灰色
      }

      # 计算甲基化分布（无bedgraph覆盖时binned为空，该unit仅占位不画柱子）
      prof <- calculate_unit_methylation_profile(
        unit_row$chrom, unit_row$start, unit_row$end, bg, bin_size)

      n_kept <- n_kept + 1L
      units_data[[n_kept]] <- list(
        start = unit_row$start, end = unit_row$end,
        binned = prof$binned, color = color)
    }
    units_data <- units_data[seq_len(n_kept)]
    # 同一NOR内的unit按基因组位置从左到右排列
    if (n_kept > 0L) {
      units_data <- units_data[order(
        vapply(units_data, function(x) x$start, numeric(1)))]
    }
    list(nor_num = nor_nums[i], units = units_data)
  })
}

draw_all_nor_rows <- function(nor_data, bin_size, scale_bp_per_in) {
  # 在当前设备上绘制所有NOR行（每个NOR一行，全部画在同一页上，PDF仅1页）
  # 布局：不用mfrow（各行绘图区宽度需随数据不同），改用par(plt)手动指定每行的
  # 绘图区（内区域的分数坐标）。所有行共用同一scale（bp/英寸）：每行绘图宽 =
  # 该NOR的units跨度 / scale，左缘全部对齐（y轴同一起点），右端随各自数据伸展；
  # 页宽已按最长行设置，长NOR完整显示（不裁剪）。
  # 每行的绘图帧(plot.new)在此统一创建：首行直接开帧，后续行先par(new=TRUE)
  # 再开帧（不翻页、不清除已画内容）；行内绘制不再单独开帧。
  num_nor <- length(nor_data)
  # mai清零（绘图区位置完全由plt控制），omi留出整页外边距（底部放x轴标题）
  par(mai = c(0, 0, 0, 0), omi = PAGE_OMI, cex = 1)
  fin <- par("fin")  # 内区域尺寸（英寸）= 页面 - 外边距
  block_h <- ROW_PLOT_H_IN + ROW_MAI[1L] + ROW_MAI[3L]

  first_done <- FALSE
  for (i in seq_len(num_nor)) {
    nd <- nor_data[[i]]
    if (length(nd$units) == 0L) next  # 该NOR没有可画的unit时留空行
    span <- max(vapply(nd$units, function(x) x$end, numeric(1))) -
      min(vapply(nd$units, function(x) x$start, numeric(1)))
    w_in <- span / scale_bp_per_in                  # 该行绘图区宽度（英寸）
    y0_in <- (num_nor - i) * block_h + ROW_MAI[1L]  # 该行绘图区下缘（英寸，自内区域底部起）
    par(plt = c(ROW_MAI[2L] / fin[1L],
                (ROW_MAI[2L] + w_in) / fin[1L],
                y0_in / fin[2L],
                (y0_in + ROW_PLOT_H_IN) / fin[2L]))
    if (first_done) par(new = TRUE)  # 后续行在新区域开新帧且不清除已画内容
    plot.new()
    draw_nor_row(nd$units, nd$nor_num, bin_size)
    first_done <- TRUE
  }

  # x轴标题：仅在整页最下方的外边距中居中标注一次（8pt）；
  # outer=TRUE 时 at 以绘图区（内区域）为参照，0.5 = 水平中点 → 标题整页居中
  mtext("Consensus position (kb)", side = 1, outer = TRUE, line = 0.5,
        cex = AXIS_TITLE_CEX, adj = 0.5, at = 0.5)

  # y轴标题：仅在整页最左侧的外边距中竖排标注一次（8pt），整页竖直居中
  # （省略 at，mtext以绘图区竖直中点定位；整页所有NOR共用这一个y轴标题）
  mtext("Methylation (%)", side = 2, outer = TRUE, line = 0.5,
        cex = AXIS_TITLE_CEX)
}

plot_individual_unit_charts <- function(unit_type_df, array_to_nor, unit_positions,
                                        bg, bin_size, output_file) {
  # 收集每个NOR一行的绘图数据并绘制
  nor_data <- prepare_nor_data(unit_type_df, array_to_nor, unit_positions,
                               bg, bin_size)
  num_nor <- length(nor_data)
  if (num_nor == 0L) {
    message("No units to plot")
    return(invisible(NULL))
  }

  # 公共x轴scale：以参考NOR（NOR14）的bp/英寸为准，所有行一致
  scale_bp_per_in <- common_scale_bp_per_in(nor_data)
  if (!is.finite(scale_bp_per_in) || scale_bp_per_in <= 0) {
    stop("Cannot determine reference x-axis scale", call. = FALSE)
  }

  # 页面：每个NOR一行；宽度按最长行自适应（不限宽），所有unit完整显示
  page_w <- nor_page_width(nor_data, scale_bp_per_in)
  pdf(output_file, width = page_w, height = nor_page_height(num_nor),
      onefile = TRUE, pointsize = 8)
  draw_all_nor_rows(nor_data, bin_size, scale_bp_per_in)
  dev.off()
  cat(sprintf(paste0("Generated %s with %d NOR rows (%d units), ",
                     "common scale %.0f bp/inch, page %.2f x %.2f in\n"),
              output_file, num_nor,
              sum(vapply(nor_data, function(x) length(x$units), integer(1))),
              scale_bp_per_in, page_w, nor_page_height(num_nor)))
}

# ------------------------------------------------------------ 主函数 --
main <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% argv || "-h" %in% argv) {
    show_help()
    return(invisible(NULL))
  }

  # 输入文件默认路径（相对运行时的工作目录），命令行可用 --flag 覆盖
  unit_type_file <- get_arg(
    argv, "--unit_type_file",
    default = "../Figure3_Variations/3C_unit_type_map/unit_to_type.xls")
  array_nor_file <- get_arg(
    argv, "--array_nor_file",
    default = "../name_of_the_arrays.txt")
  bed_file <- get_arg(
    argv, "--bed_file",
    default = "../Figure2_Assembly/2A_circle_profile/named_units.bed")
  bedgraph_file <- get_arg(
    argv, "--bedgraph_file",
    default = "stats_Me.bedgraph")
  bin_size <- get_arg(argv, "--bin_size", default = 100L,
                      convert = function(x) as.integer(x))
  output <- get_arg(argv, "--output", default = "individual_unit_methylation")

  if (is.null(bin_size) || is.na(bin_size) || bin_size <= 0) {
    stop("--bin_size must be a positive integer", call. = FALSE)
  }

  cat("Loading unit type data...\n")
  unit_type_df <- load_unit_type_data(unit_type_file)
  cat(sprintf("Loaded %d units\n", nrow(unit_type_df)))

  cat("Loading Array to NOR mapping...\n")
  array_to_nor <- load_array_to_nor_mapping(array_nor_file)
  cat(sprintf("Loaded %d mappings\n", length(array_to_nor)))

  cat("Loading unit positions...\n")
  unit_positions <- load_unit_positions(bed_file)
  cat(sprintf("Loaded %d unit positions\n", nrow(unit_positions)))

  cat("Loading bedgraph methylation data...\n")
  bg <- load_bedgraph(bedgraph_file)
  cat(sprintf("Loaded %d bedgraph intervals\n", nrow(bg)))

  # 绘制所有units的单独图表（按NOR分组）
  cat("Plotting individual unit charts...\n")
  output_file <- sprintf("%s.pdf", output)
  plot_individual_unit_charts(unit_type_df, array_to_nor, unit_positions,
                              bg, bin_size, output_file)

  cat("Done!\n")
}

# 等同于Python的 if __name__ == "__main__"：直接用Rscript执行时才运行main()
if (sys.nframe() == 0L) {
  main()
}
