# NOR Array Assembly, Methylation and Variation Annotation Pipelines

This repository contains the scripts of three pipelines:

| Pipeline | Directory | Description |
| --- | --- | --- |
| Assembly | `scripts/assembly/` | Starts from ONT reads containing rDNA and assembles the full NOR1–NOR18 arrays by barcode-path joining |
| Methylation | `scripts/methylation/` | Builds query/template sequences and, after minimap2 alignment, produces a BAM with accurate CIGAR records as modkit input to call 5mC/5hmC |
| SV & annotation | `scripts/SV_and_annotation/` | k-mer mapping → GFF annotation → SV detection and R plotting |

Result data: `NOR_k562.fa`, `NOR_k562.fa.fai`, `NOR_k562.gff`, `NOR_Methylation.bedgraph`.

## Dependencies

- Python 3: `biopython`, `pandas`, `numpy`, `matplotlib`, `scikit-learn`, `openpyxl`
- External tools: `mafft`, `minimap2`, `samtools`, `tabix`, `modkit` ([releases](https://github.com/nanoporetech/modkit/releases))
- R (only base R, for `simplify_variation_structure_merged_R.R`)

In the commands below, `<...>` is a parameter placeholder and must be replaced with the real path/file name when running.

---

## 1. Assembly

Input: ONT (Dorado) reads containing the rDNA region. Output: `NOR_k562.fa` (18 NOR arrays).

### 1.1 Extract rDNA-containing reads
Use `grep` to select reads carrying the Prime editing spacer sequence from the ONT sequence data.

### 1.2 Deduplication
Dorado also emits duplicate reads from supplementary alignments; if the sequence of a read is fully contained in another read, it is removed.

```bash
python dedup.py --fasta <input.fasta> --out <output.fasta>
```

### 1.3 Barcode extraction

```bash
python get_bc_allow_mismatch.py <input.fasta> <output_barcode_table.xls>
```

### 1.4 Manual inspection of the barcode distribution
Check for wrongly included or missed barcodes.

### 1.5 Flank detection (unmatched head/tail ends)

```bash
python find_flank.py -i <input.fasta> -r <reference_unit.fasta> -o <output_dir>
```

### 1.6 Barcode frequency correction
Manually correct barcodes with too low a frequency (<=25) using the local context.

### 1.7 Split reads by barcode count

```bash
python divide_reads_by_barcodes.py --xlsx <barcode_table.xlsx> --outdir <output_dir>
```

Two read tables are produced (barcode count = 2 and > 2), used later for alignment and assembly respectively.

### 1.8 Assembly (barcode-path joining)
Reads with `bc>2` are assembled: overlapping reads are joined; if several paths are possible (ABC/BCD/BCE → ABCD or ABCE),
the inter-barcode interval distance decides which one is taken; after each join the contained short paths are removed and the
next round of joining starts. This yields the head–tail paths, and the remaining paths are completed by single-barcode links.

```bash
python assemble_reads.py -i <input_reads_gt2bc.csv>
```

Outputs: barcode path table `path.txt`, assembly graph `assembly.gfa`, path graph `path.gfa`, etc.

### 1.9 Align to the barcode-only scaffold

```bash
python align.py -i <barcode_path_file.txt>
```

Paths that align to several loci are first clustered by k-means on inter-barcode interval distances; the interval sequences are
then BLASTed against the corresponding barcodes of the already anchored long paths, and the locus with the highest percent
identity is assigned.

### 1.10 Check and discard interval-length outliers

```bash
python testify_interval.py <alignment_table.xlsx>
```

### 1.11 Extract interval fragments and group them

```bash
python get_aligned_fragments.py <alignment_table.xlsx> <input.fasta> <output_dir>
```

The interval fragments of each read are written to the output directory, grouped by assembly / fragment (this is the
`aligned_fragments/` used by the methylation pipeline).

### 1.12 Consensus generation
A multiple sequence alignment is built per fragment and a gap-aware consensus is called:

```bash
mafft --thread <threads> --nofft <fragment.fasta> > <aligned_fragment.fasta>
python gap_aware_consensus_with_score.py <aligned_fragment.fasta>
```

### 1.13 Link the consensus into the final arrays

```bash
python link_consensus.py <trimmed_consensus_dir> <barcode_path_file.txt>
```

---

## 2. Methylation

Input: the `aligned_fragments/` and `consensus_trimed/` directories produced by the assembly pipeline, the barcode path table,
and the assembled reference genome fasta.
Output: methylation-ratio bedgraph in 50 bp bins.

### 2.1 Generate query and template sequences

```bash
python format_query.py <aligned_fragments_dir> <barcode_path_file.txt>
python get_templete_for_methylation.py <trimmed_consensus_dir> <barcode_path_file.txt>
```

> The template directory holds the consensus produced by the assembly.

### 2.2 minimap2 alignment (one SAM per array × fragment)

```bash
minimap2 -t <threads> -ax map-pb <template.fasta> <query.fasta> > <output.sam>
```

### 2.3 Extract the informative fields of the SAMs

```bash
awk -v file=<sample_name> -v OFS='\t' '$1!~/^@/ && $2==0 {print file,$1,$4,length($10),$6}' <input.sam>
```

Collected into an alignment info table (fragment id / read id / start / length / CIGAR).

### 2.4 Add alignment length and further fields

```bash
python cigar_to_lengths.py <align_info_table.xls> <template_dir> <formatted_align_info.xls>
```

### 2.5 Build the merged CIGARs

```bash
python map_read_cigar.py <formatted_align_info.xls> <output_cigar_table.txt>
```

### 2.6 BAM generation and indexing

```bash
samtools view -bS <merged_alignment.sam> > <output.bam>
samtools sort -@ <threads> <output.bam> -o <sorted.bam>
samtools index -@ <threads> <sorted.bam>
```

### 2.7 Methylation calling with modkit

```bash
modkit pileup <sorted.bam> <output_pileup.bed.gz> --modified-bases 5mC 5hmC \
              --reference <reference.fasta> --combine-mods --bgzf --log <log.txt>
tabix -p bed <output_pileup.bed.gz>
modkit stats --regions <bin_region.bed> --out-table <stats_table.bed> <output_pileup.bed.gz>
```

The 50 bp bin regions are generated from the fai:

```bash
awk -v bin=<bin_size> '{for(i=1;i<=$2;i+=bin){end=i+bin-1; if(end>$2) end=$2; print $1"\t"i-1"\t"end}}' <reference.fasta.fai> > <bin_region.bed>
```

Finally, the ratio column is taken out to give the methylation bedgraph.

---

## 3. SV_and_annotation

Input: reference genome fasta, the query fasta assembled for each sample, the feature table of the reference unit, and the bed
defining the units. `generate_gff.py` does `import kmer_map`, so the scripts have to be run from the same directory.

### 3.1 k-mer mapping (query → reference)

```bash
python kmer_map.py <reference.fasta> <query.fasta> <output_map.txt> [-s <query_step>] [-f <flank.bed>] [-k <kmer_size>] [-b <bin_size>]
```

- The reference is indexed by k-mer (default 20; if not given, the smallest k without duplicated k-mers is chosen automatically);
  the query is scanned with a fixed step (default 10) and both strands are matched;
- `-f` masks matches falling inside flank regions, `-b` removes redundant matches;
- Outputs a three-column mapping table `query / query_pos / ref_index` (e.g. `unit_structure_map.txt`, the default input name of the next step).

### 3.2 Generate the GFF annotation

```bash
python generate_gff.py <reference.fasta> <features.txt> <query.fasta> <output.gff> [<units.bed>]
```

- `<features.txt>`: feature table of the reference unit (ID / type / start / end);
- `<units.bed>`: defines the unit partition of each assembly (default `named_units.bed`, 5 columns: assembly / start / end / name / length);
- Outputs the annotated GFF: every unit together with its internal features (repeated features are split by terminal motif and copy length and numbered y.n).

### 3.3 Summarize the structurally different intervals

```bash
python process_diff.py --input <structure_map.txt>
```

`diff = query_pos - ref_index` is computed and, per assembly, the stretches with an unchanged diff are recorded (start/end
coordinates, diff value, row count); isolated jumps whose diff differs by >3 kb from both neighbours while the two neighbours
differ by <1 kb are filtered out. The output name is fixed to `diff_changes.txt`.

### 3.4 Plotting in R

```r
# edit input_file / output_pdf at the top of the script, then run
Rscript simplify_variation_structure_merged_R.R
```

`diff_changes.txt` is stacked from NOR1 down to NOR18 into 18 panels: the black segments give the query-versus-reference
stretches and adjacent segments are joined by red lines, with Query Position on the x axis and Ref Index on the y axis, so the
unit composition and length variation of each NOR array can be read off directly. Outputs
`simplify_variation_structure_merged_R.pdf` (plus a PNG preview).

### 3.5 k-mer PCA of the unit sequences

```bash
python pca.py --bed <units.bed> --fasta <reference.fasta> --k <kmer_len> \
              --components <n_components> --clusters <n_clusters> --output <output_prefix>
```
