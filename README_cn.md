# NOR 阵列组装、甲基化与变异注释流程

本仓库包含三条流程的脚本：

| 流程 | 目录 | 说明 |
| --- | --- | --- |
| 组装 (Assembly) | `scripts/assembly/` | 从含 rDNA 的 ONT reads 出发，按 barcode 路径组装出 NOR1–NOR18 的完整序列 |
| 甲基化 (Methylation) | `scripts/methylation/` | 生成 query/模板序列，minimap2 比对后生成带有准确CIGAR信息的BAM文件的作为 modkit的输入， 统计 5mC/5hmC |
| 结构变异与注释 (SV & annotation) | `scripts/SV_and_annotation/` | k-mer mapping → 生成 GFF 注释 → 找SV与 R 作图 |

结果数据：`NOR_k562.fa`、`NOR_k562.fa.fai`、`NOR_k562.gff`、`NOR_Methylation.bedgraph`。

## 环境依赖

- Python 3：`biopython`、`pandas`、`numpy`、`matplotlib`、`scikit-learn`、`openpyxl`
- 外部工具：`mafft`、`minimap2`、`samtools`、`tabix`、`modkit`（[下载页](https://github.com/nanoporetech/modkit/releases)）
- R（画图脚本 `simplify_variation_structure_merged_R.R`，仅用 base R）

以下命令中的 `<...>` 为参数占位符，实际运行时替换为真实路径/文件名。

---

## 1. Assembly

输入：含 rDNA 区域的 ONT（Dorado）reads。输出：`NOR_k562.fa`（18 条 NOR 阵列）。

### 1.1 提取含 rDNA 的 reads
用 `grep`，筛出ONT sequence data中带Prime editing spacer序列的 reads。

### 1.2 去重
Dorado 会输出补充比对产生的重复 read，若某条 read 的序列被另一条完全包含则去掉。

```bash
python dedup.py --fasta <input.fasta> --out <output.fasta>
```

### 1.3 提取 barcode

```bash
python get_bc_allow_mismatch.py <input.fasta> <output_barcode_table.xls>
```

### 1.4 人工检查 barcode 分布
检查误加入/遗漏的 barcode。

### 1.5 找 flank（head/tail 未比对上的末端）

```bash
python find_flank.py -i <input.fasta> -r <reference_unit.fasta> -o <output_dir>
```

### 1.6 修正 barcode 频率
结合上下文把频率过低（<=25）的 barcode 人工修正

### 1.7 按 barcode 数量拆分 reads

```bash
python divide_reads_by_barcodes.py --xlsx <barcode_table.xlsx> --outdir <output_dir>
```

输出 barcode 数为 2 与 >2 的两张 read 表，分别供后面的比对与组装使用。

### 1.8 组装（barcode 路径拼接）
用 `bc>2` 的 reads 组装：有重叠的拼接；若出现多条可能路径（ABC/BCD/BCE → ABCD 或 ABCE），
用 barcode 之间的 interval 距离判断选哪条；每完成一次拼接就把被包含的短 path 移除并进入下一轮循环。
结果得到 head–tail 路径，剩余路径再由单 barcode 连接补齐。

```bash
python assemble_reads.py -i <input_reads_gt2bc.csv>
```

输出：barcode path 表 `path.txt`、组装图 `assembly.gfa`、path 图 `path.gfa` 等。

### 1.9 比对到 barcode-only scaffold

```bash
python align.py -i <barcode_path_file.txt>
```

能比对到多个位置的 path，先按 barcode 间 interval 距离做 k-means 分组，再把 interval 序列与已确定的长 path 做 blast，取 identity 最高的位置归属。

### 1.10 检查 interval 长度离群值并剔除

```bash
python testify_interval.py <alignment_table.xlsx>
```

### 1.11 提取 interval 片段并分组

```bash
python get_aligned_fragments.py <alignment_table.xlsx> <input.fasta> <output_dir>
```

按 assembly / fragment 分组把 reads 上的 interval 片段写入输出目录（即甲基化流程使用的 `aligned_fragments/`）。

### 1.12 生成 consensus
对每个 fragment 做多序列比对后生成 gap-aware consensus：

```bash
mafft --thread <threads> --nofft <fragment.fasta> > <aligned_fragment.fasta>
python gap_aware_consensus_with_score.py <aligned_fragment.fasta>
```

### 1.13 连接 consensus 得到最终阵列

```bash
python link_consensus.py <trimmed_consensus_dir> <barcode_path_file.txt>
```

---

## 2. Methylation

输入：组装流程产出的 `aligned_fragments/`、`consensus_trimed/` 目录，barcode path 表，以及组装好的参考基因组 fasta。
输出：50 bp bin 的甲基化比例 bedgraph。

### 2.1 生成 query 与模板序列

```bash
python format_query.py <aligned_fragments_dir> <barcode_path_file.txt>
python get_templete_for_methylation.py <trimmed_consensus_dir> <barcode_path_file.txt>
```

> 模板目录为组装生成的consensus。

### 2.2 minimap2 比对（每个 array × fragment 一条 SAM）

```bash
minimap2 -t <threads> -ax map-pb <template.fasta> <query.fasta> > <output.sam>
```

### 2.3 整理 SAM 有效信息

```bash
awk -v file=<sample_name> -v OFS='\t' '$1!~/^@/ && $2==0 {print file,$1,$4,length($10),$6}' <input.sam>
```

汇总为比对信息表（fragment id / read id / 起点 / 长度 / CIGAR）。

### 2.4 补充比对长度等信息

```bash
python cigar_to_lengths.py <align_info_table.xls> <template_dir> <formatted_align_info.xls>
```

### 2.5 生成合并 CIGAR

```bash
python map_read_cigar.py <formatted_align_info.xls> <output_cigar_table.txt>
```

### 2.6 生成 BAM 并建索引

```bash
samtools view -bS <merged_alignment.sam> > <output.bam>
samtools sort -@ <threads> <output.bam> -o <sorted.bam>
samtools index -@ <threads> <sorted.bam>
```

### 2.7 modkit 统计甲基化

```bash
modkit pileup <sorted.bam> <output_pileup.bed.gz> --modified-bases 5mC 5hmC \
              --reference <reference.fasta> --combine-mods --bgzf --log <log.txt>
tabix -p bed <output_pileup.bed.gz>
modkit stats --regions <bin_region.bed> --out-table <stats_table.bed> <output_pileup.bed.gz>
```

其中 50 bp bin 的 region 由 fai 生成：

```bash
awk -v bin=<bin_size> '{for(i=1;i<=$2;i+=bin){end=i+bin-1; if(end>$2) end=$2; print $1"\t"i-1"\t"end}}' <reference.fasta.fai> > <bin_region.bed>
```

最后取出比例列即得到甲基化 bedgraph。

---

## 3. SV_and_annotation

输入：参考基因组 fasta、各样本组装得到的 query fasta、参考 unit 的 feature 表、定义 unit 划分的 bed。
`generate_gff.py` 会 `import kmer_map`，因此相关脚本需放在同一目录下运行。

### 3.1 k-mer mapping（query → 参考）

```bash
python kmer_map.py <reference.fasta> <query.fasta> <output_map.txt> [-s <query_step>] [-f <flank.bed>] [-k <kmer_size>] [-b <bin_size>]
```

- 参考序列按 k-mer（默认 20；未指定时自动取「无重复 k-mer」的最小 k）建索引，query 按步长（默认 10）扫描，正反链都匹配；
- `-f` 可屏蔽 flank 区域内的匹配点，`-b` 用于去除冗余匹配点；
- 输出三列 `query / query_pos / ref_index` 的 mapping 表（如 `unit_structure_map.txt`，即下一步的默认输入名）。

### 3.2 生成 GFF 注释

```bash
python generate_gff.py <reference.fasta> <features.txt> <query.fasta> <output.gff> [<units.bed>]
```

- `<features.txt>`：参考 unit 上的 feature 表（ID / type / start / end）；
- `<units.bed>`：定义各 assembly 的 unit 划分（默认 `named_units.bed`，5 列：assembly / start / end / name / length）；
- 输出注释 GFF：每个 unit 及其内部 feature（重复 feature 按末端 motif 与拷贝长度拆分，编号为 y.n）。

### 3.3 统计结构差异区间

```bash
python process_diff.py --input <structure_map.txt>
```

计算 `diff = query_pos - ref_index`，按 assembly 记录 diff 保持不变的区段（起止坐标、diff 值、行数），
并过滤掉「与前后 diff 差值均 >3 kb、但前后两者相差 <1 kb」的孤立跳变。固定输出 `diff_changes.txt`。

### 3.4 R 作图

```r
# 修改脚本顶部的 input_file / output_pdf 后运行
Rscript simplify_variation_structure_merged_R.R
```

把 `diff_changes.txt` 按 NOR1–NOR18 从上到下堆叠成 18 个 panel：黑线段表示 query 与参考的对应区段，
相邻线段用红线连接，x 轴为 Query Position、y 轴为 Ref Index，用于展示各 NOR 阵列的单元组成与长度变异。
输出 `simplify_variation_structure_merged_R.pdf`（及 PNG 预览）。

### 3.5 unit 序列的 k-mer PCA

```bash
python pca.py --bed <units.bed> --fasta <reference.fasta> --k <kmer_len> \
              --components <n_components> --clusters <n_clusters> --output <output_prefix>
```

