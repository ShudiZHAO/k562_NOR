#!/usr/bin/env python3

import sys
from collections import Counter
from Bio import SeqIO

def load_alignment(fasta):
    records = list(SeqIO.parse(fasta, "fasta"))
    seqs = [str(r.seq).upper() for r in records]
    names = [r.id for r in records]
    return names, seqs

def find_effective_range(seq, window=600, min_ratio=0.85):
    """
    输入：
        seq        : str，对齐后的序列（包含 '-'）
        window     : 滑动窗口大小
        min_ratio  : 窗口内 ATCG 占比阈值

    输出：
        (real_start, real_end)
        - 若该 read 没有任何真实对齐区域，返回 (None, None)
        - 坐标为 alignment 坐标（0-based）
    """

    seq = seq.upper()
    n = len(seq)

    def is_valid_base(c):
        return c in "ATCG"

    # ---------- 找 real_start ----------
    real_start = None
    for i in range(0, n - window + 1):
        win = seq[i:i + window]
        valid = sum(is_valid_base(c) for c in win)
        if valid / window >= min_ratio:
            # 找该窗口中第一个 ATCG
            for j in range(i, i + window):
                if is_valid_base(seq[j]):
                    real_start = j
                    break
            break

    # ---------- 找 real_end ----------
    real_end = None
    for i in range(n - window, -1, -1):
        win = seq[i:i + window]
        valid = sum(is_valid_base(c) for c in win)
        if valid / window >= min_ratio:
            # 找该窗口中最后一个 ATCG
            for j in range(i + window - 1, i - 1, -1):
                if is_valid_base(seq[j]):
                    real_end = j
                    break
            break

    if real_start is None or real_end is None or real_start > real_end:
        return None, None

    return real_start, real_end

def gap_aware_consensus(names, seqs):
    n_seq = len(seqs)
    aln_len = len(seqs[0])

    # 计算每条序列的有效区间
    ranges = [find_effective_range(seq) for seq in seqs]

    for i in range(n_seq):
        start, end = ranges[i]
        if start is None:
            sys.stderr.write(f"Warning: {names[i]} has no valid alignment start!")
            continue
        if end is None:
            sys.stderr.write(f"Warning: {names[i]} has no valid alignment end!")
            continue

    consensus = []

    for col in range(aln_len):
        bases = []

        for i in range(n_seq):
            start, end = ranges[i]
            if start is None:
                continue
            if end is None:
                continue
            if col < start or col > end:
                # 前导 / 尾随 gap，不算覆盖
                continue

            base = seqs[i][col]
            if base in "ATCG-":
                bases.append(base)
            # gap 或其他字符直接忽略

        if not bases:
            # coverage = 0，直接跳过，不输出，不补 N
            continue

        # coverage >= 1：ATCG-only 投票
        counter = Counter(bases)
        consensus_base = counter.most_common(1)[0][0]
        if consensus_base != '-':
            consensus.append(consensus_base)

    return "".join(consensus)

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: python gap_aware_consensus.py aln.fasta\n")
        sys.exit(1)

    aln_fasta = sys.argv[1]
    names, seqs = load_alignment(aln_fasta)

    # 简单 sanity check
    aln_len = len(seqs[0])
    for s in seqs:
        if len(s) != aln_len:
            raise ValueError("Alignment sequences have different lengths!")

    consensus = gap_aware_consensus(names, seqs)

    print(">consensus_gap_aware")
    # 每行 60 bp 输出
    for i in range(0, len(consensus), 60):
        print(consensus[i:i+60])

if __name__ == "__main__":
    main()