#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
find_flank.py —— 基于 split_bam.py 的 flank（head/tail）检测逻辑，输入 FASTA
reads（每条 record 视为一个 read），检测 read 两端不匹配参考的 flank 区域，
输出 flank.csv 与 flank.fa

处理流程
--------
1. 读入参考 FASTA（取第一条序列），建立参考全部 k-mer 集合（默认 k=20；
   含出现多次的 k-mer，含 N 的窗口跳过），与 split_bam.py 的 flank 检测一致。
2. 逐条处理输入 FASTA record（先做方向检查 orient_read）：
   a. 正、反两个方向各做一次匹配区定位（b-d）；反向互补方向的真正匹配
      区更宽时判定 read 为反向，将其反转（revcomp），后续步骤均在反转
      后的序列上进行（如 read 804）；两方向均无匹配区或宽度相同（正向
      read 的反向扫描只剩零散误匹配）时保持正向不反转；
   b. 滑窗扫描 read，命中参考 k-mer 的位置形成 hits；
   c. 连续 < _MIN_RUN_HITS（20）个命中的簇（含孤立命中）视为误匹配丢弃；
   d. 某位置起 _DENSE_SPAN（1000）bp 内有 >= _DENSE_HITS（300）个命中才
      认为是真正的匹配区 [p0, q1+k)；
   e. head = read 起点到匹配区起点；tail = 匹配区终点到 read 末尾；整条
      read 均无真正匹配区（如纯 flank read）时不产生任何 flank 记录。
3. 左/右端未匹配 >= --head-tail-len（默认 5000）bp 记为 head/tail，记录到
   flank.csv（readid, type, flank_length, end_barcode, dis_to_end,
   barcode_num）；
   左或右未匹配 >= --flank-fa-min（默认 50000）bp 时把该 flank 序列写入
   flank.fa（> <read名>_head_left / > <read名>_tail_right）。
4. flank.csv 的 barcode 各列（算法同 get_bc_allow_mismatch.py：对整条 read
   扫描模糊匹配 LEFT{e<=2} + 6-30bp 非贪婪任意序列 + RIGHT{e<=1}，忽略
   大小写）：
   * end_barcode——head 取整条 read 的第一个 barcode，tail 取最后一个；
     整条 read 无匹配时为空；
   * dis_to_end——barcode 与 flank 重复区一侧边界的间隔（通常非负，barcode
     横跨边界时为负）：
     head 为 flank 末端（即匹配区起点）与第一个 barcode 的间隔，tail 为
     flank 起点（即匹配区终点）与最后一个 barcode 的间隔（注意 head 不是
     到 read 起始的距离）。barcode 位于 flank 内或匹配区内均按此定义；
     无匹配时为 -1；
   * barcode_num——整条 read 上检出的 barcode 总数（与 end_barcode 同一次
     扫描，head/tail 两行取值相同）；regex 不可用或无匹配时为 0。

依赖：Python >= 3.8；barcode 各列需要第三方模块 regex（缺失时 end_barcode
留空、dis_to_end 为 -1、barcode_num 为 0，其余功能不受影响）

用法示例
--------
python find_flank.py -i reads.fa -r KY962518.1_start_from_feature4.fasta -o flank_out
python find_flank.py -i reads.fa -r ref.fa -o flank_out --head-tail-len 5000

注意：flank 检测参数与 split_bam.py 完全一致；barcode 模式默认与
get_bc_allow_mismatch.py 相同（LEFT=tggcggagatgggcgcc e<=2，
RIGHT=GCGAGG e<=1，barcode 6-30bp），可用 --bc-* 参数调整。
"""

import argparse
import csv
import glob
import gzip
import os
import sys

__version__ = "1.4.0"

try:                    # regex 为可选依赖：flank.csv 的 barcode 两列需要
    import regex
    HAVE_REGEX = True
except ImportError:
    regex = None
    HAVE_REGEX = False

_VALID_BASES = set("ACGTN")

# flank 信息表 flank.csv 的列（前三列与 split_bam.py 的 head_and_tail.csv
# 一致；end_barcode/dis_to_end/barcode_num 为 barcode 各列，见
# find_read_barcodes 与 pick_barcode）
FLANK_FIELDS = ["readid", "type", "flank_length", "end_barcode", "dis_to_end",
                "barcode_num"]

# flank（head/tail）检测参数：连续 >= _MIN_RUN_HITS 个 kmer 命中才算“比对上”
#（连续数不足的命中视为误匹配）；某位置起 _DENSE_SPAN bp 内有 >= _DENSE_HITS
# 个命中才认为是真正的匹配区。
_MIN_RUN_HITS = 20
_DENSE_SPAN = 1000
_DENSE_HITS = 300


def _open_text(path):
    """打开文本文件，自动处理 gzip（*.gz），并兼容 UTF-8 BOM。"""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return open(path, "r", encoding="utf-8-sig")


def load_fasta(path):
    """
    解析参考 FASTA，返回第一条序列（大写）。文件含多条序列时使用第一条并警告。
    """
    records = []
    name, chunks = None, []
    with _open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if chunks:
                    records.append((name, "".join(chunks)))
                    chunks = []
                name = line[1:].split()[0] if len(line) > 1 else "ref"
            elif line:
                chunks.append(line)
    if chunks:
        records.append((name, "".join(chunks)))
    if not records:
        sys.exit("ERROR: no sequence parsed from %s" % path)
    if len(records) > 1:
        print("[WARNING] %s contains %d sequences, using the first one (%s)"
              % (path, len(records), records[0][0]), file=sys.stderr)
    name, seq = records[0]
    seq = seq.upper()
    invalid = set(seq) - _VALID_BASES
    if invalid:
        print("[WARNING] reference contains invalid base character(s) %r, "
              "k-mers covering them are skipped"
              % "".join(sorted(invalid)), file=sys.stderr)
    return name, seq


def load_fasta_records(path):
    """
    解析输入 FASTA，返回全部 [(name, seq), ...]（大写）。每条 record 视为
    一个 read；read 名重复时给出警告（输出记录按原名保留）。
    """
    records = []
    name, chunks = None, []
    with _open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if chunks:
                    records.append((name, "".join(chunks)))
                    chunks = []
                name = line[1:].split()[0] if len(line) > 1 else "read"
            elif line:
                chunks.append(line)
    if chunks:
        records.append((name, "".join(chunks)))
    if not records:
        sys.exit("ERROR: no sequence parsed from %s" % path)
    seen = set()
    for name, _seq in records:
        if name in seen:
            print("[WARNING] %s: duplicate read name %s" % (path, name),
                  file=sys.stderr)
        seen.add(name)
    return [(name, seq.upper()) for name, seq in records]


# barcode 模糊匹配默认参数（与 get_bc_allow_mismatch.py 相同）
BC_LEFT = "tggcggagatgggcgcc"
BC_RIGHT = "GCGAGG"
BC_LEFT_ERR = 2
BC_RIGHT_ERR = 1
BC_MIN = 6
BC_MAX = 30


def build_barcode_pattern(left, right, left_err, right_err, bc_min, bc_max):
    """
    编译 barcode 模糊匹配模式（get_bc_allow_mismatch.py 的写法）：
    LEFT{e<=left_err} + bc_min..bc_max bp 非贪婪任意序列 +
    RIGHT{e<=right_err}，忽略大小写。返回编译后的 regex 模式。
    """
    pattern = ("(%s){e<=%d}(.{%d,%d}?)(%s){e<=%d}"
               % (regex.escape(left), left_err, bc_min, bc_max,
                  regex.escape(right), right_err))
    return regex.compile(pattern, regex.IGNORECASE)


def find_read_barcodes(seq, regex_pat):
    """
    在整条 read 上查找全部 barcode 匹配（与 get_bc_allow_mismatch.py 相同：
    对整条序列扫描，模式见 build_barcode_pattern）。返回按 read 内位置升序
    的 [(start, end, barcode), ...]（start/end 为 barcode 在 read 内的绝对
    位置；regex 不可用或无匹配时为 []）。
    """
    if regex_pat is None:
        return []
    return [(m.start(2), m.end(2), m.group(2))
            for m in regex_pat.finditer(seq)]


def pick_barcode(barcodes, side, boundary):
    """
    从整条 read 的全部 barcode 匹配中按 flank 行取值（get_bc_allow_mismatch.py
    的语义）：side="head" 取第一个匹配（read 起点一侧最靠前的），side="tail"
    取最后一个匹配（read 末端一侧最靠后的）。boundary 为该 flank 的重复区
    一侧边界在 read 内的位置（head 侧 = len(head_seq)，tail 侧 = 匹配区
    终点 = len(seq) - len(tail_seq)）。
    返回 (barcode, dis_to_end)：
      barcode——barcode 序列（整条 read 无匹配时为 ""）；
      dis_to_end——barcode 与 flank 重复区一侧边界之间的间隔（通常非负，barcode 横跨边界时为负）：
        head 为 flank 末端（匹配区起点）与第一个 barcode 的间隔——barcode
        在 flank 内时为 flank 末端到 barcode 末端的距离，barcode 已进入匹配
        区时为 barcode 起始到 flank 末端的距离；
        tail 为 flank 起点（匹配区终点）与最后一个 barcode 的间隔——barcode
        在 flank 内时为 barcode 起始到 flank 起点的距离，barcode 在匹配区内
        时为 flank 起点到 barcode 末端的距离；
      整条 read 无匹配时为 -1。
    """
    if not barcodes:
        return "", -1
    if side == "head":
        start, end, bc = barcodes[0]
        if end <= boundary:
            return bc, boundary - end
        return bc, start - boundary
    start, end, bc = barcodes[-1]
    if start >= boundary:
        return bc, start - boundary
    return bc, boundary - end


def build_kmer_set(ref, k):
    """
    参考全部 k-mer 集合（含出现多次的 k-mer；含 N 的窗口跳过）。
    用于 flank（head/tail）检测：单位内部的重复区段 k-mer 也可命中，
    使匹配区边界不因 unique k-mer 的缺失而提前截断。
    """
    kmers = set()
    for i in range(len(ref) - k + 1):
        km = ref[i:i + k]
        if "N" not in km:
            kmers.add(km)
    return kmers


def _keep_long_runs(hits, min_run=_MIN_RUN_HITS):
    """
    hits 为升序的 read 位置。只保留连续命中簇（相邻位置差 1）长度
    >= min_run 的命中；连续数不足的命中（含孤立命中）视为误匹配丢弃。
    返回过滤后的升序位置列表。
    """
    kept = []
    n = len(hits)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and hits[j + 1] == hits[j] + 1:
            j += 1
        if j - i + 1 >= min_run:
            kept.extend(hits[i:j + 1])
        i = j + 1
    return kept


def _first_dense_pos(hits, span=_DENSE_SPAN, min_hits=_DENSE_HITS):
    """
    hits 为升序的 read 位置（已过滤短连续簇）。返回第一个位置 p 使得
    [p, p+span] 内存在 >= min_hits 个命中（即真正匹配区的起点）；
    命中不足的 flank 区域被跳过、视为未匹配。无则返回 None。
    """
    n = len(hits)
    j = 0
    for i in range(n):
        if j < i:
            j = i
        while j < n and hits[j] <= hits[i] + span:
            j += 1
        if j - i >= min_hits:
            return hits[i]
    return None


def _last_dense_pos(hits, span=_DENSE_SPAN, min_hits=_DENSE_HITS):
    """
    对称版本（从右往左）：返回最后一个位置 p 使得 [p-span, p] 内存在
    >= min_hits 个命中；无则返回 None。
    """
    n = len(hits)
    j = n - 1
    for i in range(n - 1, -1, -1):
        if j > i:
            j = i
        while j >= 0 and hits[j] >= hits[i] - span:
            j -= 1
        if i - j >= min_hits:
            return hits[i]
    return None


def find_match_region(seq, all_kmers, k, span=_DENSE_SPAN,
                      min_hits=_DENSE_HITS, min_run=_MIN_RUN_HITS):
    """
    定位 read 的真正匹配区 [p0, q1 + k)：连续 < min_run 个 kmer 命中的段
    （含孤立命中）视为误匹配、不算比对上（_keep_long_runs）；p0 为第一个
    span bp 内 >= min_hits 个命中的位置，q1 为最后一个。all_kmers 为参考
    全部 k-mer 集合（含重复出现的 k-mer）。整条 read 均无真正匹配区（如
    纯 flank read）时返回 (None, None)，否则返回 (p0, q1)。
    """
    hits = [i for i in range(len(seq) - k + 1)
            if seq[i:i + k] in all_kmers]
    hits = _keep_long_runs(hits, min_run)
    p0 = _first_dense_pos(hits, span, min_hits)
    if p0 is None:
        return None, None
    q1 = _last_dense_pos(hits, span, min_hits)
    return p0, q1


def revcomp(seq):
    """反向互补（非 ACGTN 字符原样保留，位置随序列反转）。"""
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def orient_read(seq, all_kmers, k, span=_DENSE_SPAN,
                min_hits=_DENSE_HITS, min_run=_MIN_RUN_HITS):
    """
    k-mer 匹配时检查 read 方向：正向与反向互补各定位一次真正匹配区
    （find_match_region），反向互补方向的匹配区更宽时判定 read 为反向，
    返回其反向互补序列，后续 head/tail 与 barcode 检测均在反转后的序列
    上进行。正向参考 k-mer 集合下，反向 read 的正向扫描只剩零散误匹配、
    匹配区远窄于反转后的真正匹配区（如 read 804）；正向 read 的反向
    扫描同理不会更宽，故宽度相同（并列）时保持正向。
    两个方向均无真正匹配区（如纯 flank read）、或正向匹配区已覆盖整条
    read（反向最多打平）时保持正向不反转。
    返回 (seq, p0, q1, reversed)：seq 为定向后的序列（可能已反转），
    (p0, q1) 为其真正匹配区，reversed=True 表示已反转。
    """
    p0, q1 = find_match_region(seq, all_kmers, k, span, min_hits, min_run)
    if p0 is not None and p0 == 0 and q1 + k >= len(seq):
        return seq, p0, q1, False      # 正向已覆盖整条 read，反向最多打平
    rc = revcomp(seq)
    r_p0, r_q1 = find_match_region(rc, all_kmers, k, span, min_hits, min_run)
    if r_p0 is None:                   # 反向无匹配区：保持正向
        return seq, p0, q1, False
    if p0 is None:                     # 仅反向有匹配区：read 为反向
        return rc, r_p0, r_q1, True
    if r_q1 + k - r_p0 > q1 + k - p0:  # 反向匹配区更宽：read 为反向
        return rc, r_p0, r_q1, True
    return seq, p0, q1, False          # 正向更宽或并列：保持正向


def find_flanks(seq, all_kmers, k, min_len, span=_DENSE_SPAN,
                min_hits=_DENSE_HITS, min_run=_MIN_RUN_HITS,
                match_region=None):
    """
    检测 read 两端不匹配参考的 flank 区域（head / tail）：
      head = seq 起点到匹配区起点；tail = 匹配区终点到 seq 末尾。
    匹配区由 find_match_region 判定；match_region 提供已算好的 (p0, q1)
    时不再重复扫描。head/tail 长度 >= min_len 时返回对应序列，否则为
    None；整条 read 均无真正匹配区时返回 (None, None)。返回 (head, tail)。
    """
    p0, q1 = (match_region if match_region is not None
              else find_match_region(seq, all_kmers, k, span,
                                     min_hits, min_run))
    if p0 is None:
        return None, None
    head = seq[:p0] if p0 >= min_len else None
    tail = seq[q1 + k:] if len(seq) - q1 - k >= min_len else None
    return head, tail
def collect_paths(inputs):
    """整理输入：支持目录、通配符与具体文件路径，返回去重后的文件列表。"""
    exts = ("*.fa", "*.fasta", "*.fna", "*.fa.gz", "*.fasta.gz", "*.fna.gz")
    paths = []
    for item in inputs:
        if os.path.isdir(item):
            found = []
            for ext in exts:
                found.extend(glob.glob(os.path.join(item, ext)))
            if not found:
                print("[WARNING] no FASTA file found in directory %s" % item,
                      file=sys.stderr)
            paths.extend(sorted(found))
        elif any(ch in item for ch in "*?["):
            found = sorted(glob.glob(item))
            if not found:
                print("[WARNING] wildcard %s matched no file" % item,
                      file=sys.stderr)
            paths.extend(found)
        else:
            paths.append(item)
    unique, seen = [], set()
    for path in paths:
        key = os.path.abspath(path)
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="find_flank.py",
        description="Detect head/tail flanks (unmatched ends) of FASTA "
                    "reads against a reference repeat unit using the same "
                    "k-mer logic as split_bam.py; writes flank.csv and "
                    "flank.fa.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="input FASTA file(s) holding the reads; "
                             "multiple files, a directory or wildcards are "
                             "accepted")
    parser.add_argument("-r", "--reference", required=True,
                        help="reference FASTA holding one repeat unit "
                             "(e.g. KY962518.1_start_from_feature4.fasta)")
    parser.add_argument("-o", "--out-dir", default="find_flank_output",
                        help="output directory")
    parser.add_argument("-k", "--kmer-size", type=int, default=20,
                        help="k-mer size for the reference k-mer set")
    parser.add_argument("--head-tail-len", type=int, default=5000,
                        help="min consecutive unmatched bp at read ends to "
                             "call a head/tail flank")
    parser.add_argument("--flank-fa-min", type=int, default=25000,
                        help="min unmatched bp for a head/tail flank to be "
                             "written to flank.fa")
    parser.add_argument("--bc-left", default=BC_LEFT,
                        help="barcode LEFT flanking sequence "
                             "(fuzzy-matched against reads)")
    parser.add_argument("--bc-right", default=BC_RIGHT,
                        help="barcode RIGHT flanking sequence "
                             "(fuzzy-matched against reads)")
    parser.add_argument("--bc-left-err", type=int, default=BC_LEFT_ERR,
                        help="max errors (edit distance) allowed in the "
                             "LEFT flanking sequence")
    parser.add_argument("--bc-right-err", type=int, default=BC_RIGHT_ERR,
                        help="max errors (edit distance) allowed in the "
                             "RIGHT flanking sequence")
    parser.add_argument("--bc-min", type=int, default=BC_MIN,
                        help="min barcode length between LEFT and RIGHT")
    parser.add_argument("--bc-max", type=int, default=BC_MAX,
                        help="max barcode length between LEFT and RIGHT")
    return parser.parse_args(argv)
def main(argv=None):
    # 输出统一为 UTF-8，避免特殊字符在不同终端编码下报错
    for _stream in (sys.stdout, sys.stderr):
        try:
            _stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass

    args = parse_args(argv)

    if args.kmer_size < 5:
        sys.exit("ERROR: --kmer-size must be >= 5")
    if args.head_tail_len < 1:
        sys.exit("ERROR: --head-tail-len must be >= 1")
    if args.flank_fa_min < 1:
        sys.exit("ERROR: --flank-fa-min must be >= 1")
    if args.bc_left_err < 0 or args.bc_right_err < 0:
        sys.exit("ERROR: --bc-left-err/--bc-right-err must be >= 0")
    if not (1 <= args.bc_min <= args.bc_max):
        sys.exit("ERROR: --bc-min/--bc-max must satisfy 1 <= bc-min <= bc-max")

    fa_paths = collect_paths(args.input)
    if not fa_paths:
        sys.exit("ERROR: no input FASTA file found")
    for path in fa_paths:
        if not os.path.isfile(path):
            sys.exit("ERROR: input file does not exist: %s" % path)
    if not os.path.isfile(args.reference):
        sys.exit("ERROR: reference file does not exist: %s" % args.reference)

    ref_name, ref_seq = load_fasta(args.reference)
    k = args.kmer_size
    if len(ref_seq) < k:
        sys.exit("ERROR: reference (%d bp) is shorter than --kmer-size (%d)"
                 % (len(ref_seq), k))
    all_kmers = build_kmer_set(ref_seq, k)
    print("[INFO] reference %s: %d bp, %d distinct %d-mer(s)"
          % (ref_name, len(ref_seq), len(all_kmers), k))
    os.makedirs(args.out_dir, exist_ok=True)

    if HAVE_REGEX:
        bc_pat = build_barcode_pattern(args.bc_left, args.bc_right,
                                       args.bc_left_err, args.bc_right_err,
                                       args.bc_min, args.bc_max)
        print("[INFO] barcode pattern: %s{e<=%d} + %d-%dbp + %s{e<=%d}"
              % (args.bc_left, args.bc_left_err, args.bc_min, args.bc_max,
                 args.bc_right, args.bc_right_err))
    else:
        bc_pat = None
        print("[WARNING] module 'regex' not installed, barcode detection "
              "disabled: end_barcode left empty and dis_to_end set to -1 "
              "(pip install regex)", file=sys.stderr)

    stats = {"reads": 0, "flank_reads": 0, "no_match": 0,
             "reversed_reads": 0, "barcode_rows": 0}
    flank_info_rows = []
    flank_fa_rows = []

    for path in fa_paths:
        print("[Processing] %s" % path)
        for name, seq in load_fasta_records(path):
            stats["reads"] += 1
            # k-mer 匹配时检查方向：反向 read 反转后再进行下面的步骤
            seq, p0, q1, is_rev = orient_read(seq, all_kmers, k)
            if is_rev:
                stats["reversed_reads"] += 1
                print("[INFO] read %s: reverse orientation detected, read "
                      "reverse-complemented before flank/barcode detection"
                      % name, file=sys.stderr)
            head_seq, tail_seq = find_flanks(
                seq, all_kmers, k, args.head_tail_len,
                match_region=(p0, q1))
            if p0 is None:
                stats["no_match"] += 1
                print("[WARNING] read %s: no match region found against the "
                      "reference, no flank record" % name, file=sys.stderr)
            if head_seq is None and tail_seq is None:
                continue
            stats["flank_reads"] += 1
            # barcode：对整条 read 扫描一次（同 get_bc_allow_mismatch.py），
            # head 行取第一个匹配，tail 行取最后一个匹配；barcode_num 为
            # 整条 read 检出的 barcode 总数
            barcodes = find_read_barcodes(seq, bc_pat)
            bc_num = len(barcodes)
            if head_seq is not None:
                bc, dis = pick_barcode(barcodes, "head", len(head_seq))
                flank_info_rows.append((name, "head", len(head_seq),
                                        bc, dis, bc_num))
                if bc:
                    stats["barcode_rows"] += 1
                if len(head_seq) >= args.flank_fa_min:
                    flank_fa_rows.append(("%s_head_left" % name, head_seq))
            if tail_seq is not None:
                bc, dis = pick_barcode(barcodes, "tail",
                                       len(seq) - len(tail_seq))
                flank_info_rows.append((name, "tail", len(tail_seq),
                                        bc, dis, bc_num))
                if bc:
                    stats["barcode_rows"] += 1
                if len(tail_seq) >= args.flank_fa_min:
                    flank_fa_rows.append(("%s_tail_right" % name, tail_seq))

    # ---- 写出 flank.csv ----
    flank_csv_path = os.path.join(args.out_dir, "flank.csv")
    with open(flank_csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.writer(fh)
        writer.writerow(FLANK_FIELDS)
        writer.writerows(flank_info_rows)

    # ---- 写出 flank.fa ----
    flank_fa_path = os.path.join(args.out_dir, "flank.fa")
    with open(flank_fa_path, "w", encoding="utf-8", newline="") as fh:
        for fname, fseq in flank_fa_rows:
            fh.write(">%s\n%s\n" % (fname, fseq))

    print("=" * 62)
    print("Input FASTA files                   : %d" % len(fa_paths))
    print("Total reads                         : %d" % stats["reads"])
    print("Reads with no match region          : %d" % stats["no_match"])
    print("Reads reversed (reverse orientation): %d" % stats["reversed_reads"])
    print("Reads with head/tail flank          : %d" % stats["flank_reads"])
    print("Flank rows written                  : %d" % len(flank_info_rows))
    print("Flank rows with barcode             : %d" % stats["barcode_rows"])
    print("-" * 62)
    print("Flank info    : %s (%d rows)"
          % (flank_csv_path, len(flank_info_rows)))
    print("Long flanks   : %s (%d seqs)"
          % (flank_fa_path, len(flank_fa_rows)))
    return 0


if __name__ == "__main__":
    main()