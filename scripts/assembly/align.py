#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 barcode 表里的 read 比对到 assembly 的 barcode 路径上。

输入（默认值都在 asm/ 下）：
  reads_2bc.csv     含 2 个 barcode 的 read，每行一个 barcode
  reads_gt2bc.csv   含 2 个以上 barcode 的 read，每行一个 barcode
  correct_list.csv  barcode 纠错表（按 read_id + barcode_index 替换）
  remove.txt        需要剔除的 read_id（每行一个）
  -i/--input        assembly fasta（'>名字' + 每行一条 barcode 路径）

流程：读入两个 barcode 表 -> 按 correct_list.csv 纠错 -> 去掉 remove.txt
的 read -> 逐个 read 在 assembly 路径里找完全匹配 -> 写 alignments.xls
（每个命中一行）和 unaligned.xls（没有命中的 read，每个 read 一行）。

两个输出都是制表符分隔，列固定为 read_id / assembly / start_index /
end_index / bc_num / barcodes / intervals / first_barcode_to_start /
last_barcode_to_end；未比对的行 assembly / start_index / end_index 写 NA，
空（barcode_length=0 的占位）barcode 统一写成 empty。

start_index 是 assembly 里 barcode 路径的 token 下标（0-based），
end_index = start_index + barcode 数（半开区间）；同一条 read 在同一个
assembly 里命中多处时，每处写一行。
"""

import argparse
import os

import pandas as pd

from assemble_reads import (
    DEFAULT_CORRECT_LIST,
    DEFAULT_REMOVE,
    HERE,
    apply_barcode_corrections,
    natural_key,
    read_correct_list,
    read_remove_list,
)

DEFAULT_READS_2BC = os.path.join(HERE, 'asm', 'reads_2bc.csv')
DEFAULT_READS_GT2BC = os.path.join(HERE, 'asm', 'reads_gt2bc.csv')

MIN_BARCODES = 2
DISABLED_VALUES = ('none', 'no', 'off', '-')
REQUIRED_COLUMNS = ('read_id', 'barcode_index', 'barcode', 'barcode_start',
                    'barcode_end', 'interval', 'total_length')


def is_disabled(value):
    """'none' / 'no' / 'off' / '-' 表示关闭这一步"""
    return str(value).strip().lower() in DISABLED_VALUES


def assembly_id(header, index):
    """assembly编号取 '>' 后第一个字段；为空时退回 Assembly_<index>"""
    fields = header.split()
    if fields:
        return fields[0]
    return 'Assembly_%d' % index


def load_assembly(assembly_file):
    """从input文件加载assembly路径，返回 (assembly编号, barcode列表) 的列表"""
    assemblies = []
    if assembly_file:
        try:
            entries = []
            current_header = ""
            current_units = []
            with open(assembly_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_header:
                            entries.append((current_header, current_units))
                        current_header = line[1:]
                        current_units = []
                    elif line:
                        units = line.split()
                        current_units.extend(units)
                if current_header:
                    entries.append((current_header, current_units))
            seen = set()
            for index, (header, units) in enumerate(entries, 1):
                name = assembly_id(header, index)
                if name in seen:
                    print(f"[WARN] duplicate assembly id {name} in {assembly_file}")
                seen.add(name)
                assemblies.append((name, units))
            print(f"[INFO] load {assembly_file} success, {len(assemblies)} assemblies")
        except Exception as e:
            print(f"[ERROR] load {assembly_file} failed: {e}")
    return assemblies


def to_number(value):
    """数值化单元格（'' / nan / 非数字 -> None）"""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if pd.isna(number):
        return None
    return number


def distance_to_read_end(total_length, last_barcode_end):
    """read 长度 - 最后一个 barcode 的结尾；任一项不可用时返回 None"""
    total = to_number(total_length)
    end = to_number(last_barcode_end)
    if total is None or end is None:
        return None
    return total - end


def load_read_table(paths):
    """读入两个 barcode 表（reads_2bc.csv / reads_gt2bc.csv）并合并"""
    frames = []
    for path in paths:
        if is_disabled(path):
            print(f"[INFO] reads      : {path} disabled")
            continue
        if not os.path.isfile(path):
            raise SystemExit(f"[ERROR] read file not found: {path}")
        df = pd.read_csv(path)
        missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
        if missing:
            raise SystemExit(f"[ERROR] {path} lacks column {missing}")
        df = df[df['read_id'].notna()].copy()
        df['read_id'] = df['read_id'].astype(str).str.strip()
        df = df[df['read_id'] != '']
        frames.append(df)
        print(f"[INFO] load {path} success, {len(df)} rows / {df['read_id'].nunique()} reads")
    if not frames:
        raise SystemExit('[ERROR] no read file given (--reads_2bc / --reads_gt2bc)')
    df = pd.concat(frames, ignore_index=True, sort=False)
    duplicated = int(df.duplicated(['read_id', 'barcode_index']).sum())
    if duplicated:
        print(f"[WARN] {duplicated} duplicated read_id/barcode_index rows in the read tables")
    return df


def read_barcode_rows(df):
    """按 barcode_index 排序，把每个 read 的 barcode 整理成

    {read_id: [(barcode, barcode_index, interval, barcode_start,
                barcode_end, total_length), ...]}
    """
    rows = {}
    for rid, group in df.groupby('read_id', sort=False):
        # 按barcode_index排序，非数字索引排到最后
        sort_key = pd.to_numeric(group['barcode_index'], errors='coerce')
        sorted_group = group.assign(
            _sort_index=sort_key.fillna(float('inf'))
        ).sort_values('_sort_index', kind='mergesort')
        items = []
        for _, row in sorted_group.iterrows():
            barcode = row['barcode']
            # 空占位（barcode_length=0）在表里可能是空单元格，统一写成 empty
            if pd.isna(barcode):
                barcode = 'empty'
            items.append((str(barcode), row['barcode_index'], row['interval'],
                          row['barcode_start'], row['barcode_end'],
                          row['total_length']))
        rows[rid] = items
    return rows


def correct_barcodes(rows, correct_list):
    """用 correct_list.csv 纠正 barcode

    复用 assemble_reads.apply_barcode_corrections 的规则：只有 read 里存在该
    barcode_index 且当前 barcode 与 origin 一致时才替换，否则跳过并告警。
    """
    stats = {'entries': 0, 'reads': 0, 'skipped': 0, 'warnings': [],
             'applied': 0, 'missing_index': 0, 'mismatch': 0, 'same': 0,
             'reads_changed': 0, 'missing_read': 0, 'missing_reads': []}
    if is_disabled(correct_list):
        print(f"[INFO] correct    : disabled (--correct_list {correct_list})")
        return rows, stats
    if not os.path.isfile(correct_list):
        print(f"[WARN] correct    : {correct_list} not found, no barcode corrections")
        return rows, stats
    corrections, read_stats = read_correct_list(correct_list)
    stats.update(read_stats)
    if not corrections:
        print(f"[INFO] correct    : {correct_list} holds no usable barcode entry, no corrections")
        return rows, stats
    for rid, entries in corrections.items():
        if rid not in rows:
            stats['missing_read'] += len(entries)
            stats['missing_reads'].append(rid)
            continue
        # apply_barcode_corrections 需要 (barcode, ends, index, interval, start, ...)
        items = [(item[0], 'no', item[1], item[2], item[3], item[4], item[5])
                 for item in rows[rid]]
        items, changed = apply_barcode_corrections(rid, items, entries, stats)
        rows[rid] = [(item[0], item[2], item[3], item[4], item[5], item[6])
                     for item in items]
        if changed:
            stats['reads_changed'] += 1
    unapplied = (stats['skipped'] + stats['missing_index'] + stats['mismatch']
                 + stats['missing_read'])
    print(f"[INFO] correct    : {correct_list} ({stats['entries'] + stats['skipped']} rows / "
          f"{stats['reads']} reads): {stats['applied']} applied to {stats['reads_changed']} reads, "
          f"{unapplied} skipped (index missing/ambiguous {stats['missing_index']}, "
          f"barcode differs {stats['mismatch']}, read not in input {stats['missing_read']}, "
          f"unusable row {stats['skipped']})"
          + (f", {stats['same']} already corrected" if stats['same'] else ''))
    for warning in stats['warnings']:
        print(f"[WARN] {warning}")
    if stats['missing_reads']:
        names = stats['missing_reads']
        print(f"[WARN] {len(names)} correct-list read(s) are not in the input: "
              f"{', '.join(names[:10])}{' ...' if len(names) > 10 else ''}")
    return rows, stats


def build_barcode_info(rows, min_barcodes=MIN_BARCODES):
    """整理成比对/输出需要的 read 信息，丢掉 barcode 少于 min_barcodes 的 read"""
    barcode_info = {}
    dropped = 0
    for rid, items in rows.items():
        if len(items) < min_barcodes:
            dropped += 1
            continue
        barcodes = [item[0] for item in items]
        intervals = [item[2] for item in items]
        start_pos = [item[3] for item in items]
        end_pos = [item[4] for item in items]
        lengths = [to_number(item[5]) for item in items]
        lengths = [length for length in lengths if length is not None]
        total_length = lengths[0] if lengths else None
        # 最后一个 barcode 结尾到 read 末端的距离（缺 total_length 时记 0）
        dis_to_end = distance_to_read_end(total_length, end_pos[-1])
        if dis_to_end is None:
            dis_to_end = 0
        barcode_info[rid] = {
            'barcodes': barcodes,
            'intervals': intervals,
            'start_pos': start_pos,
            'end_pos': end_pos,
            'total_length': total_length,
            'first_barcode_to_start': start_pos[0],
            'last_barcode_to_end': dis_to_end,
            'bc_num': len(items)
        }
    print(f"[INFO] {len(barcode_info)} reads with >= {min_barcodes} barcodes"
          + (f", {dropped} reads skipped" if dropped else ''))
    return barcode_info


def apply_remove_list(barcode_info, remove_file):
    """去掉 remove.txt 里列出的 reads"""
    if is_disabled(remove_file):
        print(f"[INFO] remove     : disabled (--remove {remove_file})")
        return barcode_info
    if not os.path.isfile(remove_file):
        print(f"[WARN] remove     : {remove_file} not found, no reads removed")
        return barcode_info
    remove_ids, stats = read_remove_list(remove_file)
    if not remove_ids:
        print(f"[INFO] remove     : {remove_file} holds no read id, no reads removed")
        return barcode_info
    missing = sorted(remove_ids - set(barcode_info), key=natural_key)
    removed = [rid for rid in barcode_info if rid in remove_ids]
    for rid in removed:
        del barcode_info[rid]
    print(f"[INFO] remove     : {remove_file} ({len(remove_ids)} ids, "
          f"{stats['duplicates']} listed twice): {len(removed)} reads removed, "
          f"{len(barcode_info)} reads kept"
          + (f", {len(missing)} ids not in the input" if missing else ''))
    if missing:
        print(f"[WARN] {len(missing)} remove-list read(s) are not in the input: "
              f"{', '.join(missing[:10])}{' ...' if len(missing) > 10 else ''}")
    if not barcode_info:
        raise SystemExit('[ERROR] every read was removed, nothing to align')
    return barcode_info


def find_path_positions(assembly, path):
    """在 assembly（'|'.join(barcodes)）里找出 path 所有完全匹配的起始下标
    """
    if not path:
        return []
    positions = []
    haystack = '|' + assembly + '|'
    needle = '|' + path + '|'
    start = haystack.find(needle)
    while start != -1:
        positions.append(haystack.count('|', 0, start))
        start = haystack.find(needle, start + 1)
    return positions


def align(assemblies, barcode_info):
    """比对barcode到assembly，允许一个read比对到同一个assembly的多个位置"""
    alignments = []
    str_assemblies = []
    for name, units in assemblies:
        str_assemblies.append((name, '|'.join(units)))
    for rid, info in barcode_info.items():
        path = '|'.join(info['barcodes'])
        for name, assembly in str_assemblies:
            # 不break，允许一个read比对到多个位置
            for index in find_path_positions(assembly, path):
                print(f"[INFO] path {rid} aligned to {name} at {index}")
                alignments.append({
                    'read_id': rid,
                    'assembly': name,
                    'start_index': index,
                    'end_index': index + len(info['barcodes']),
                    'barcodes': info['barcodes'],
                    'intervals': info['intervals'],
                    'first_barcode_to_start': info['first_barcode_to_start'],
                    'last_barcode_to_end': info['last_barcode_to_end'],
                    'bc_num': info['bc_num']
                })
    # 按start_index排序
    alignments = sorted(alignments, key=lambda x: (x['start_index'], -x['bc_num'], -x['first_barcode_to_start']))
    print(f'[INFO] {len(barcode_info)} paths in total, {len(alignments)} alignments align to {len(assemblies)} assemblies')
    return alignments


def collect_unaligned(alignments, read_info):
    """收集没有比对到任何assembly的read，每个read一条记录

    未比对的read同样按 alignments.xls 的列排布，assembly /
    start_index / end_index 三列填 NA。
    """
    mapped = {aln['read_id'] for aln in alignments}
    unaligned = []
    for rid, info in read_info.items():
        if rid in mapped:
            continue
        unaligned.append({
            'read_id': rid,
            'assembly': 'NA',
            'start_index': 'NA',
            'end_index': 'NA',
            'bc_num': info['bc_num'],
            'barcodes': info['barcodes'],
            'intervals': info['intervals'],
            'first_barcode_to_start': info['first_barcode_to_start'],
            'last_barcode_to_end': info['last_barcode_to_end']
        })
    # barcode数多的排前面，barcode数相同时按barcode路径聚合，便于查看
    unaligned = sorted(unaligned,
                       key=lambda x: (-int(to_number(x['bc_num']) or 0),
                                      '|'.join(map(str, x['barcodes'])),
                                      str(x['read_id'])))
    print(f"[INFO] {len(unaligned)} reads not aligned to any assembly")
    return unaligned


def write_alignment_table(records, out_file):
    """按 alignments.xls 的列写表（制表符分隔）"""
    header = ('read_id', 'assembly', 'start_index', 'end_index', 'bc_num',
              'barcodes', 'intervals', 'first_barcode_to_start', 'last_barcode_to_end')
    with open(out_file, "w") as f:
        f.write('\t'.join(header) + '\n')
        for aln in records:
            barcodes_str = '|'.join(map(str, aln['barcodes']))
            intervals_str = '|'.join(map(str, aln['intervals']))
            f.write(f"{aln['read_id']}\t{aln['assembly']}\t{aln['start_index']}\t{aln['end_index']}\t{aln['bc_num']}\t{barcodes_str}\t{intervals_str}\t{aln['first_barcode_to_start']}\t{aln['last_barcode_to_end']}\n")


def write_alignments_file(alignments, out_file):
    write_alignment_table(alignments, out_file)
    print(f"[INFO] {len(alignments)} alignments written to {out_file}")


def write_unaligned_file(unaligned, out_file):
    write_alignment_table(unaligned, out_file)
    print(f"[INFO] {len(unaligned)} unmapped reads written to {out_file}")


def main():
    parser = argparse.ArgumentParser(description="Align barcodes to assembly")
    parser.add_argument("-i", "--input", required=True, help="input assembly fasta file")
    parser.add_argument("--reads_2bc", "--reads-2bc", dest="reads_2bc",
                        default=DEFAULT_READS_2BC,
                        help="reads with 2 barcodes (CSV)")
    parser.add_argument("--reads_gt2bc", "--reads-gt2bc", dest="reads_gt2bc",
                        default=DEFAULT_READS_GT2BC,
                        help="reads with more than 2 barcodes (CSV)")
    parser.add_argument("--correct_list", "--correct-list", dest="correct_list",
                        default=DEFAULT_CORRECT_LIST,
                        help="barcode correction table ('none' to disable)")
    parser.add_argument("--remove", default=DEFAULT_REMOVE,
                        help="reads to drop before alignment ('none' to disable)")
    parser.add_argument("--alignments_out", default="alignments.xls", help="output alignments file")
    parser.add_argument("--unaligned_out", default="unaligned.xls",
                        help="output unmapped reads file (same format as alignments file)")
    args = parser.parse_args()

    # 1) 读入两个 barcode 表（reads_2bc.csv + reads_gt2bc.csv）
    df = load_read_table([args.reads_2bc, args.reads_gt2bc])
    # 2) 按 barcode_index 汇总每个 read 的 barcode
    rows = read_barcode_rows(df)
    # 3) 用 correct_list.csv 纠错
    rows, _corr_stats = correct_barcodes(rows, args.correct_list)
    # 4) 整理成比对信息，并去掉 remove.txt 里的 reads
    barcode_info = build_barcode_info(rows)
    barcode_info = apply_remove_list(barcode_info, args.remove)
    # 5) 解析assembly并比对
    assemblies = load_assembly(args.input)
    alignments = align(assemblies, barcode_info)
    # 6) 整理未比对的read（每个read一行）
    unaligned = collect_unaligned(alignments, barcode_info)
    # 写入比对结果
    write_alignments_file(alignments, args.alignments_out)
    write_unaligned_file(unaligned, args.unaligned_out)
    print("[INFO] Done!")

if __name__ == "__main__":
    main()