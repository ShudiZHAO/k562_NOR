#!/usr/bin/env python3
"""
模拟连接过程，为每个 consensus 序列生成包含 N 的完整长度序列

输入目录结构为 consensus_trimed/NORx/<y>_consensus.fasta（与 trim_consensus_ends.py 输出一致）。
"""

import os
import re
import sys

from Bio import SeqIO

# barcode 文件中的空占位符：占一个 junction 槽位，但长度为 0
EMPTY_BARCODE_TOKENS = ('', 'empty', 'nan', 'none')

# consensus_trimed/NORx/ 下的 fragment 文件名约定：<y>_consensus.fasta
FRAGMENT_PATTERN = re.compile(r'^(\d+)_consensus\.fasta$')


def normalize_barcode(token):
    """把 empty 等占位符转成空字符串，其余 barcode 去空白并转大写。"""
    token = token.strip()
    if token.lower() in EMPTY_BARCODE_TOKENS:
        return ''
    return token.upper()


def read_barcode_info(barcode_file):
    """
    读取 barcode 信息
    返回一个字典：{array_name: [barcode1, barcode2, ...]}
    empty 占位符保留槽位，但转成空字符串
    """
    barcode_info = {}
    current_assembly = None
    barcodes = []

    with open(barcode_file, 'r', encoding='utf-8-sig') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_assembly and barcodes:
                    barcode_info[current_assembly] = barcodes
                current_assembly = line[1:].strip()
                barcodes = []
            elif line:
                # 分割多个 barcode；empty 表示该位置没有 barcode
                line_barcodes = [normalize_barcode(x) for x in line.split()]
                barcodes.extend(line_barcodes)

    # 处理最后一个 array
    if current_assembly and barcodes:
        barcode_info[current_assembly] = barcodes
    return barcode_info

def get_fragment_sequences(consensus_dir, assembly):
    """
    获取 consensus_dir/NORx/ 下的所有 fragment 序列
    返回一个字典：{y: (sequence, filename)}
    """
    sequences = {}
    assembly_dir = os.path.join(consensus_dir, assembly)

    if not os.path.isdir(assembly_dir):
        return sequences

    for filename in os.listdir(assembly_dir):
        match = FRAGMENT_PATTERN.match(filename)
        if not match:
            continue
        y = int(match.group(1))

        # 读取序列
        filepath = os.path.join(assembly_dir, filename)
        for record in SeqIO.parse(filepath, 'fasta'):
            sequences[y] = (str(record.seq), filename)
            break  # 只取第一个序列

    return sequences

def calculate_positions(sequences, barcodes):
    """
    计算每个 fragment 在完整序列中的位置
    返回一个字典：{y: (start, end)}
    """
    sorted_ys = sorted(sequences.keys())
    positions = {}
    current_pos = 0
    
    for i, y in enumerate(sorted_ys):
        seq_len = len(sequences[y][0])
        positions[y] = (current_pos, current_pos + seq_len)
        
        # 更新当前位置
        current_pos += seq_len
        if i < len(sorted_ys) - 1 and i < len(barcodes):
            # junction 只添加 barcode 长度
            current_pos += len(barcodes[i])
    
    return positions, current_pos

def generate_padded_sequence(target_y, sequences, barcodes, total_length, positions):
    """
    为指定的 fragment 生成包含 N 的完整长度序列
    """
    # 初始化结果序列为 N
    result = ['N'] * total_length
    
    # 获取当前 fragment 的信息
    current_seq, _ = sequences[target_y]
    start, end = positions[target_y]
    
    # 填充当前 fragment 的序列
    for i in range(start, end):
        if i < total_length:
            result[i] = current_seq[i - start]
    
    # 非开头的 fragment 在左端只添加 barcode
    sorted_ys = sorted(sequences.keys())
    target_index = sorted_ys.index(target_y)
    if target_index > 0:
        # 计算添加 barcode 的位置
        barcode_index = target_index - 1
        if barcode_index < len(barcodes):
            barcode = barcodes[barcode_index]
            insert_start = start - len(barcode)
            insert_end = start
            if insert_start >= 0:
                for i in range(insert_start, insert_end):
                    result[i] = barcode[i - insert_start]
    
    return ''.join(result)

def main():
    if len(sys.argv) != 3:
        print("Usage: python get_templete_for_methylation.py <consensus_dir> <barcode_file>")
        sys.exit(1)
    
    consensus_dir = sys.argv[1]
    barcode_file = sys.argv[2]
    
    # 读取 barcode 信息
    barcode_info = read_barcode_info(barcode_file)
    
    # 处理每个 array（NORx）
    output_dir = 'templete_for_methylation'
    os.makedirs(output_dir, exist_ok=True)
    for assembly, barcodes in barcode_info.items():
        # 获取 fragment 序列
        sequences = get_fragment_sequences(consensus_dir, assembly)
        
        if not sequences:
            print(f"Warning: No fragments found for {assembly} "
                  f"(expected {assembly}/<y>_consensus.fasta in {consensus_dir})")
            continue
        
        # 计算每个 fragment 的位置和总长度
        positions, total_length = calculate_positions(sequences, barcodes)
        
        # 为每个 fragment 生成包含 N 的完整长度序列
        for y, (_, filename) in sequences.items():
            padded_seq = generate_padded_sequence(y, sequences, barcodes, total_length, positions)
            
            # 输出文件名补上 array 前缀，例如 NOR1_0_consensus.fasta
            output_file = os.path.join(output_dir, f"{assembly}_{filename}")
            with open(output_file, 'w') as f:
                f.write(">consensus_gap_aware\n")
                # 每行 60 个碱基
                for i in range(0, len(padded_seq), 60):
                    f.write(padded_seq[i:i+60] + "\n")
            
            print(f"Created {output_file}")

if __name__ == "__main__":
    main()