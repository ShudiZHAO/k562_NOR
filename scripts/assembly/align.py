#!/usr/bin/env python3

import argparse
import plotly.graph_objects as go
from collections import defaultdict
import pandas as pd
import os

def load_assembly(assembly_file):
    """从input文件加载assembly路径"""
    assemblies = []
    if assembly_file:
        try:
            current_header = ""
            current_units = []
            with open(assembly_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_header:
                            assemblies.append(current_units)
                        current_header = line[1:]
                        current_units = []
                    else:
                        units = line.split()
                        current_units.extend(units)
                if current_header:
                    assemblies.append(current_units)
            print(f"[INFO] load {assembly_file} success, {len(assemblies)} assemblies")
        except Exception as e:
            print(f"[ERROR] load {assembly_file} failed: {e}")
    return assemblies

def load_barcode_info(barcode_file):
    """从xlsx文件加载barcode信息"""
    try:
        df = pd.read_excel(barcode_file)
        # 按read_id分组，收集barcode信息
        barcode_info = {}
        for rid, group in df.groupby('read_id'):
            # 只处理包含2个或以上barcode的read
            if len(group) >= 2:
                barcodes = []
                intervals = []
                start_pos = []
                end_pos = []
                total_length = group.iloc[0]['total_length']
                dis_to_end = group.iloc[-1]['dis_to_end'] if not pd.isna(group.iloc[-1]['dis_to_end']) else 0
                # 按barcode_index排序
                sorted_group = group.sort_values('barcode_index')
                for _, row in sorted_group.iterrows():
                    barcodes.append(row['barcode'])
                    intervals.append(row['interval'])
                    start_pos.append(row['barcode_start'])
                    end_pos.append(row['barcode_end'])
                # 计算第一个barcode到read起始的距离
                first_barcode_start = start_pos[0]
                bc_num = len(group)
                if 'bc_num' in df.columns:
                    bc_num = group.iloc[0]['bc_num']
                barcode_info[rid] = {
                    'barcodes': barcodes,
                    'intervals': intervals,
                    'start_pos': start_pos,
                    'end_pos': end_pos,
                    'total_length': total_length,
                    'first_barcode_to_start': first_barcode_start,
                    'last_barcode_to_end': dis_to_end,
                    'bc_num': bc_num
                }
        print(f"[INFO] load {len(barcode_info)} reads with 2 or more barcodes from {barcode_file}")
        return barcode_info
    except Exception as e:
        print(f"[ERROR] load {barcode_file} failed: {e}")
        return {}

def align(assemblies, barcode_info):
    """比对barcode到assembly，允许一个read比对到不同的位置"""
    alignments = []
    str_assemblies = []
    for assembly in assemblies:
        str_assemblies.append('|'.join(assembly))
    for rid, info in barcode_info.items():
        path = info['barcodes']
        for i, assembly in enumerate(str_assemblies, 1):
            if '|'.join(path) in assembly:
                print(f"[INFO] 路径 {rid} 比对到 Assembly_{i}")
                index = int(assembly.find('|'.join(path))/7)
                alignments.append({
                    'read_id': rid,
                    'assembly': f'Assembly_{i}',
                    'start_index': index,
                    'end_index': index + len(path),
                    'barcodes': path,
                    'intervals': info['intervals'],
                    'first_barcode_to_start': info['first_barcode_to_start'],
                    'last_barcode_to_end': info['last_barcode_to_end'],
                    'bc_num': info['bc_num']
                })
                # 不break，允许一个read比对到多个位置
                # break
    # 按start_index排序
    alignments = sorted(alignments, key=lambda x: (x['start_index'], -x['bc_num'], -x['first_barcode_to_start']))
    print(f'[INFO] {len(barcode_info)} paths in total, {len(alignments)} alignments align to {len(assemblies)} assemblies')
    return alignments

def write_alignments_file(alignments, out_file):
    with open(out_file, "w") as f:
        f.write("read_id\tassembly\tstart_index\tend_index\tbc_num\tbarcodes\tintervals\tfirst_barcode_to_start\tlast_barcode_to_end\n")
        for aln in alignments:
            barcodes_str = '|'.join(aln['barcodes'])
            intervals_str = '|'.join(map(str, aln['intervals']))
            f.write(f"{aln['read_id']}\t{aln['assembly']}\t{aln['start_index']}\t{aln['end_index']}\t{aln['bc_num']}\t{barcodes_str}\t{intervals_str}\t{aln['first_barcode_to_start']}\t{aln['last_barcode_to_end']}\n")
    print(f"[INFO] {len(alignments)} alignments written to {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Align barcodes to assembly")
    parser.add_argument("-i", "--input", required=True, help="input assembly fasta file")
    parser.add_argument("--barcode_file", required=True, help="barcode info")
    parser.add_argument("--alignments_out", default="alignments.xls", help="output alignments file")
    args = parser.parse_args()

    # 加载barcode信息
    barcode_info = load_barcode_info(args.barcode_file)
    # 解析输入文件
    assemblies = load_assembly(args.input)
    # 执行比对
    alignments = align(assemblies, barcode_info)
    # 写入比对结果
    write_alignments_file(alignments, args.alignments_out)
    print("[INFO] Done!")

if __name__ == "__main__":
    main()