# 根据simplify_variation_structure_merged.py复刻的R版画图脚本
# 读取diff_changes.txt，按NOR从上到下（NOR1在顶部）堆叠18个panel
# 注：diff_changes.txt的Assembly列已经是NOR名称（NOR1..NOR18，按数值顺序写出），
#     故不再读取name_of_the_arrays.txt，也不再做Assembly -> NOR的编号转换，
#     panel的分组/排序/标签全部直接使用Assembly列的原值
# 内容：黑色线段（Start/End_Query_Pos 对 Start/End_Ref_Index，端点带圆点marker），
#       相邻线段之间用红色线连接；x轴以该NOR最左数据点为0、右端扩1000bp，ylim自动加5%边距
#
# 坐标系：数值缩小1000以kb显示，轴标题带(k)；x/y都只保留10(k)一档的label+tick，
#   不再画2(k)的次刻度；tick向外，刻度标签贴近坐标轴
# 比例尺：panel高固定，x/y采用同一比例尺（pt/bp）=> panel宽 = x_range×比例尺，
#   x轴数值-长度比例与y轴严格一致（各NOR的y_range相同，x从0起算=>panel宽与x_max成正比）；
#   画布宽度由最宽panel（NOR 5）自动推算
#
# 轴标题布局：
#   1) NOR编号放在原y轴标题位置：panel左侧、水平、垂直居中，右缘距panel左缘约12pt
#   2) x轴标题"Query Position (kb)"全图只出现一次：画布最底部、水平居中
#   3) y轴标题"Ref Index (kb)"全图只出现一次：画布最左侧、垂直居中（旋转90°）
#   4) x轴左侧无空白：每个panel先把Query坐标平移到"该NOR最左数据点=0"（xlim左端=0），
#      即第一个数据点（最左边的黑线段起点）正好落在0刻度、也就是panel左缘上
#   实现方式：18个panel画完后，再用一个覆盖整页的"标题层"（plot.new + text，坐标0~1即画布），
#   按pt精确换算位置绘制两个共享标题 —— 比mtext(outer=)/oma的定位更直观、可控
#
# format.txt排版：不加边框只留坐标轴；轴标题 8pt、刻度标签 5pt、NOR编号 8pt；
#   tick长1.2pt、坐标轴/刻度线0.25pt（tcl>0，tick朝外）；
#   行高91pt（行间间隔16 + panel 60 + 下方刻度标签区15），另在画布底部/左侧各留16pt给共享轴标题

# ---- 输入输出 ----
input_file <- "diff_changes.txt"   # Assembly列即NOR名称，无需额外映射文件
output_pdf <- "simplify_variation_structure_merged_R.pdf"
output_png <- "simplify_variation_structure_merged_R_preview.png"  # 预览用

# ---- format.txt 排版参数 ----
axis_title_pt <- 8      # 轴标题 8pt
axis_text_pt  <- 5      # 刻度标签(axis text) 5pt
nor_text_pt   <- 8      # NOR编号 8pt
tick_len_pt   <- 1.2    # tick长 1.2pt
axis_lwd_pt   <- 0.25   # 坐标轴/刻度线 0.25pt
seg_lwd_pt    <- 0.7    # 数据线段粗细（format.txt的0.25pt指坐标轴/刻度线；
                        #   数据墨水取0.6pt保证窄panel内可见，可按需调整）
conn_lwd_pt   <- 0.25    # 红色连接线粗细（Python为lw=1.5）
marker_cex    <- 0       # 端点圆点已去除（原0.35；>0才绘制）

# ---- 坐标轴刻度/间距参数 ----
label_step_bp  <- 10000   # 刻度间隔：label与tick都只保留10(k)一档（不再画2(k)次刻度）
x_start_bp     <- 0       # x轴左端固定从0开始；数据先按"最左数据点=0"平移，故0就是第一个数据点
x_pad_bp       <- 1000    # x轴右端额外留白
mgp_label_line <- 0.13    # 刻度标签距轴：数字右缘距轴约2.5pt（0.3时为5pt，太远）
xlabel_line    <- -0.4    # x轴数字用mtext绘制的行位置（axis()过密时会自动丢标签，故数字改用mtext）
                          # 实测：数字顶距轴约2.3pt、基线约6pt（0.3时顶距12.4pt太远）；与外向tick尖留约1.1pt
nor_label_line <- 0.83    # NOR编号行位置（原y轴标题的位置）：mtext(side=2,las=1,adj=1)的右缘恰好落在
                          #   line×csi = line×14.4pt 处，0.83→右缘距panel左缘约12pt；
                          #   y刻度数字左缘约7.4pt（右缘1.9pt+宽5.6pt），故间隙约4.6pt

# ---- 画布/版面参数（pt）----
# 实测语义（R4.5.1 300dpi标定）：1个margin line = 1.2 x pointsize = 14.4pt；
# axis(mgp=)的标签位置参数不生效，mgp须经par()设置；mtext基线约在轴外(line+0.82)行处
strip_h_pt      <- 16           # 行间间隔（原为NOR条带标签区，标签已移到panel左侧）
axes_h_pt       <- 60           # panel（绘图区）高度
chrome_below_pt <- 15           # panel下方刻度标签区（轴标题已移出panel，只放刻度数字）
row_height_pt   <- strip_h_pt + axes_h_pt + chrome_below_pt   # 每行高度 = 16+60+15 = 91pt
mai_left_pt     <- 46           # 左边距：y刻度数字(≈5.6pt) + 间隙4.6pt + 水平NOR编号(≈28.9pt) + 余量
mai_right_pt    <- 5
title_left_pt   <- 18           # 画布最左侧共享y轴标题区宽度（8pt旋转标题高约11pt，居中→左缘约2.5pt）
title_bottom_pt <- 18           # 画布最底部共享x轴标题区高度（居中→下缘约4.5pt）
panel_left_pt   <- title_left_pt + 5 + mai_left_pt   # panel(绘图区)左缘距画布左缘

# ---------- 1. 读取数据 ----------
# Assembly列已经是NOR名称（NOR1..NOR18），直接使用，不做任何转换
df <- read.delim(input_file, stringsAsFactors = FALSE)
stopifnot(all(c("Assembly", "Start_Query_Pos", "End_Query_Pos",
                "Start_Ref_Index", "End_Ref_Index") %in% names(df)))

# 按Assembly分组：保持其在diff_changes.txt中首次出现的顺序
# （文件按NOR1->NOR18顺序写出，故NOR1在顶部，与Python groupby排序结果一致）
nor_labels <- unique(df$Assembly)
n_panels   <- length(nor_labels)
# 画布高度 = 各行高度之和 + 底部共享x轴标题区
canvas_height_pt <- title_bottom_pt + n_panels * row_height_pt

# lwd换算：R图形设备的lwd单位为1/96英寸，pt需乘96/72
pt2lwd <- function(pt) pt * 96 / 72
# tcl换算：单位为基准行高par("csi")（英寸）。
# 注意：R4.5.1/Windows实测与经典文档相反——tcl正值刻度朝面板内、负值朝外（300dpi标定实验证实），
# 故取负号使刻度向外伸出1.2pt
pt2tcl <- function(pt) -((pt / 72) / par("csi"))

# ---------- 2. 等比布局：x/y同一比例尺 ----------
# 每行先算xlim/ylim，比例尺 = axes_h_pt / y_range（pt/bp），
# panel宽 = x_range x 比例尺（x轴数值-长度比例与y轴一致）；
# 画布宽度由最宽panel自动推算
row_info <- lapply(nor_labels, function(label) {
  g <- df[df$Assembly == label, ]
  g <- g[order(g$Start_Query_Pos), ]
  x0 <- min(g$Start_Query_Pos)                      # 该panel最左数据点（也就是全panel最左的墨水）
  g$Start_Query_Pos <- g$Start_Query_Pos - x0       # 4) 平移到"最左数据点=0"：去掉panel左侧空白，
  g$End_Query_Pos   <- g$End_Query_Pos   - x0       #    第一个数据点正好落在0刻度(=panel左缘)上
  x_max <- max(g$End_Query_Pos)
  xlim <- c(x_start_bp, x_max + x_pad_bp)           # 左端=x_start_bp(=0)，右端扩1000bp
  ymin <- min(g$Start_Ref_Index, g$End_Ref_Index)
  ymax <- max(g$Start_Ref_Index, g$End_Ref_Index)
  span <- max(ymax - ymin, 1e-9)
  ylim <- c(ymin - 0.05 * span, ymax + 0.05 * span) # matplotlib默认5%边距
  list(label = label, g = g, xlim = xlim, ylim = ylim)
})
panel_w_pt     <- vapply(row_info, function(r) diff(r$xlim) * axes_h_pt / diff(r$ylim), numeric(1))
panel_w_max_pt <- max(panel_w_pt)
canvas_width_pt <- panel_left_pt + panel_w_max_pt + mai_right_pt + 8

# ---------- 3. 单panel绘制（完全复刻Python绘图逻辑） ----------
draw_panel <- function(r, panel_w_pt) {
  g    <- r$g
  xlim <- r$xlim; ylim <- r$ylim

  plot.new()
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")
  box(bty = "l", lwd = pt2lwd(axis_lwd_pt))         # 不加边框，只留L形坐标轴

  # 红色线连接跳跃的点（相邻线段之间）
  n <- nrow(g)
  if (n > 1) {
    segments(g$End_Query_Pos[-n], g$End_Ref_Index[-n],
             g$Start_Query_Pos[-1], g$Start_Ref_Index[-1],
             col = "red", lwd = pt2lwd(conn_lwd_pt), lend = 1)
  }
  # 黑色线段
  segments(g$Start_Query_Pos, g$Start_Ref_Index,
           g$End_Query_Pos, g$End_Ref_Index,
           col = "black", lwd = pt2lwd(seg_lwd_pt), lend = 0)
  # 端点圆点marker（Python marker='o', markersize=1）
  if (marker_cex > 0) {
    points(c(g$Start_Query_Pos, g$End_Query_Pos),
           c(g$Start_Ref_Index, g$End_Ref_Index),
           pch = 16, cex = marker_cex, col = "black")
  }

  # x轴：数值缩小1000以kb显示；只在每10(k)处画一档tick（不再有2(k)次刻度），tick向外
  # 注：axis()对放不下的标签会自动省略（panel内间距约12pt时5pt数字被丢一半），
  #     因此刻度线用axis(labels=FALSE)，数字用mtext(at=)逐个绘制，保证每10(k)都显示
  lab_at  <- seq(ceiling(xlim[1] / label_step_bp) * label_step_bp, xlim[2], by = label_step_bp)
  axis(1, at = lab_at, labels = FALSE,
       tcl = pt2tcl(tick_len_pt), lwd = 0, lwd.ticks = pt2lwd(axis_lwd_pt))
  mtext(format(lab_at / 1000, scientific = FALSE, trim = TRUE), side = 1,
        at = lab_at, line = xlabel_line, cex = axis_text_pt / 12)
  # 注：x轴标题已改为全图共享（画在画布最底部、水平居中，见draw_all），此处不再绘制

  # y轴：与x轴同一规则（缩小1000、只保留10(k)一档、向外、贴近）
  ylab_at <- seq(ceiling(ylim[1] / label_step_bp) * label_step_bp, ylim[2], by = label_step_bp)
  axis(2, at = ylab_at, labels = format(ylab_at / 1000, scientific = FALSE, trim = TRUE),
       tcl = pt2tcl(tick_len_pt), las = 1,
       lwd = 0, lwd.ticks = pt2lwd(axis_lwd_pt), cex.axis = axis_text_pt / 12)
  # NOR编号：放在原y轴标题的位置——panel左侧、水平、垂直居中（8pt）
  # 直接用Assembly列的原值（NOR1..NOR18），不做"NOR + 编号"的拼装
  # mtext(side=2, las=1, adj=1)：右缘恰好落在 line×csi(=14.4pt) 处，故取 nor_label_line
  # 注：y轴标题"Ref Index (kb)"已改为全图共享（画在画布最左侧、垂直居中，见draw_all）
  mtext(r$label, side = 2, line = nor_label_line,
        las = 1, adj = 1, cex = nor_text_pt / 12)
}

# ---------- 4. 整页绘制 ----------
draw_all <- function() {
  par(oma = c(0, 0, 0, 0))
  for (k in seq_len(n_panels)) {
    r       <- row_info[[k]]
    panel_w <- panel_w_pt[k]
    # 行k（从上往下）：row_top为该行顶缘的绝对pt坐标
    row_top <- canvas_height_pt - (k - 1) * row_height_pt
    par(fig = c((panel_left_pt - mai_left_pt) / canvas_width_pt,
                (panel_left_pt + panel_w + mai_right_pt) / canvas_width_pt,
                (row_top - row_height_pt) / canvas_height_pt,
                row_top / canvas_height_pt),
        mai = c(chrome_below_pt, mai_left_pt, strip_h_pt, mai_right_pt) / 72,
        las = 1, mgp = c(0, mgp_label_line, 0))
    # 多panel堆叠：第2个起必须new=TRUE，否则plot.new()会清除整张画布
    if (k > 1) par(new = TRUE)
    draw_panel(r, panel_w)
  }

  # ---- 4.1 全图共享的轴标题（各只出现一次）----
  # 覆盖整页的"标题层"：fig=整页 + mar=0 => plot区=整页，用户坐标(0,1)x(0,1)即整张画布
  # 先par(fig=...)再单独par(new=TRUE)（与panel写法一致），保证不清除已画好的18个panel
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), las = 1)
  par(new = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  # x轴标题：整张图最底部、水平居中（y取底部标题区title_bottom_pt的中线）
  text(0.5, (title_bottom_pt / 2) / canvas_height_pt, "Consensus Position (kb)",
       cex = axis_title_pt / 12)
  # y轴标题：整张图最左侧、垂直居中（x取左侧标题区title_left_pt的中线），旋转90°
  text((title_left_pt / 2) / canvas_width_pt, 0.5, "Reference position (kb)",
       cex = axis_title_pt / 12, srt = 90)
}

# ---------- 5. 输出 ----------
pdf(output_pdf, width = canvas_width_pt / 72, height = canvas_height_pt / 72,
    pointsize = 12, useDingbats = FALSE)
draw_all()
dev.off()
