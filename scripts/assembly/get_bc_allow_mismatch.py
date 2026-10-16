import regex
from Bio import SeqIO
import sys

# -------- 参数区（你可以改） --------
if len(sys.argv) != 3:
    print("Usage: python get_bc_new.py <fasta> <out>")
    sys.exit(1)

FASTA = sys.argv[1]
OUT = sys.argv[2]

# -------- 参数区 --------
LEFT = "tggcggagatgggcgcc"
RIGHT = "GCGAGG"

LEFT_ERR = 2
RIGHT_ERR = 1

BARCODE_MIN = 6
BARCODE_MAX = 30
# -----------------------

pattern = (
    f"({LEFT}){{e<={LEFT_ERR}}}"
    f"(.{{{BARCODE_MIN},{BARCODE_MAX}}}?)"
    f"({RIGHT}){{e<={RIGHT_ERR}}}"
)

regex_pat = regex.compile(pattern, regex.IGNORECASE)
# with open(FASTA, "r", encoding="utf-8") as handle:
with open(OUT, "w") as out:
    out.write("read_id\tbarcode_index\tbarcode\tbarcode_start\tbarcode_end\tbarcode_length\tinterval\ttotal_length\n")
    for record in SeqIO.parse(FASTA, "fasta"):
        seq = str(record.seq)
        total_length = len(seq)
        idx = 0
        for m in regex_pat.finditer(seq):                
            idx += 1
            barcode = m.group(2)
            start = m.start(2)
            end = m.end(2)
            if idx == 1:
                interval = 0
            else:
                interval = start - start_last
            start_last = start
            out.write(f"{record.id}\t{idx}\t{barcode}\t{start}\t{end}\t{end-start}\t{interval}\t{total_length}\n")
            print(f"{record.id}\t{idx}\t{barcode}\t{start}\t{end}\t{end-start}\t{interval}\t{total_length}\n")


