#!/usr/bin/env python3
"""
模拟连接过程，为每个consensus_seq中的序列生成包含N的完整长度序列
"""

import os
import sys
from Bio import SeqIO

def read_barcode_info(barcode_file):
    """
    读取 barcode 信息
    返回一个字典：{assembly: [barcode1, barcode2, ...]}
    """
    barcode_info = {}
    current_assembly = None
    barcodes = []
    
    with open(barcode_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_assembly and barcodes:
                    barcode_info[current_assembly] = barcodes
                current_assembly = line[1:]
                barcodes = []
            elif line:
                # 分割多个 barcode
                line_barcodes = line.split()
                barcodes.extend(line_barcodes)
    
    # 处理最后一个 assembly
    if current_assembly and barcodes:
        barcode_info[current_assembly] = barcodes
    return barcode_info

def get_fragment_sequences(consensus_dir, assembly):
    """
    获取指定 assembly 的所有 fragment 序列
    返回一个字典：{y: (sequence, filename)}
    """
    sequences = {}
    prefix = f"{assembly}_"
    
    for filename in os.listdir(consensus_dir):
        if filename.startswith(prefix) and filename.endswith('.fasta'):
            # 提取 y 值
            try:
                y_str = filename[len(prefix):].split('_')[0]
                y = int(y_str)
            except:
                continue
            
            # 读取序列
            filepath = os.path.join(consensus_dir, filename)
            for record in SeqIO.parse(filepath, 'fasta'):
                sequences[y] = (str(record.seq), filename)
                break  # 只取第一个序列
    
    return sequences

def calculate_positions(sequences, barcodes, prefix_seq, suffix_seq):
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
            # 添加 barcode 长度
            current_pos += len(prefix_seq) + len(barcodes[i]) + len(suffix_seq)
    
    return positions, current_pos

def generate_padded_sequence(target_y, sequences, barcodes, prefix_seq, suffix_seq, total_length, positions):
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
    
    # 非结尾的 fragment 在右端添加 prefix_seq
    sorted_ys = sorted(sequences.keys())
    target_index = sorted_ys.index(target_y)
    if target_index < len(sorted_ys) - 1:
        # 计算添加 prefix_seq 的位置
        prefix_start = end
        prefix_end = prefix_start + len(prefix_seq)
        for i in range(prefix_start, prefix_end):
            if i < total_length:
                result[i] = prefix_seq[i - prefix_start]
    
    # 非开头的 fragment 在左端添加 barcode 和 suffix_seq
    if target_index > 0:
        # 计算添加 barcode 和 suffix_seq 的位置
        barcode_index = target_index - 1
        if barcode_index < len(barcodes):
            barcode = barcodes[barcode_index]
            insert_seq = barcode + suffix_seq
            insert_start = start - len(insert_seq)
            insert_end = start
            if insert_start >= 0:
                for i in range(insert_start, insert_end):
                    result[i] = insert_seq[i - insert_start]
    
    return ''.join(result)

def main():
    if len(sys.argv) != 3:
        print("Usage: python get_templete_for_methylation.py <consensus_dir> <barcode_file>")
        sys.exit(1)
    
    consensus_dir = sys.argv[1]
    barcode_file = sys.argv[2]
    
    # 读取 barcode 信息
    barcode_info = read_barcode_info(barcode_file)
    
    # 定义 barcode 前后的序列
    prefix_seq = "TGGCGGAGATGGGCGCC"
    suffix_seq = "GCGAGG"
    
    # 处理每个 assembly
    output_dir = 'templete_for_methylation'
    os.makedirs(output_dir, exist_ok=True)
    for assembly, barcodes in barcode_info.items():
        # 获取 fragment 序列
        sequences = get_fragment_sequences(consensus_dir, assembly)
        
        if not sequences:
            print(f"Warning: No fragments found for Assembly_{assembly}")
            continue
        
        # 计算每个 fragment 的位置和总长度
        positions, total_length = calculate_positions(sequences, barcodes, prefix_seq, suffix_seq)
        
        # 为每个 fragment 生成包含 N 的完整长度序列
        for y, (_, filename) in sequences.items():
            padded_seq = generate_padded_sequence(y, sequences, barcodes, prefix_seq, suffix_seq, total_length, positions)
            
            # 输出结果
            output_file = os.path.join(output_dir, filename)
            with open(output_file, 'w') as f:
                f.write(">consensus_gap_aware\n")
                # 每行 60 个碱基
                for i in range(0, len(padded_seq), 60):
                    f.write(padded_seq[i:i+60] + "\n")
            
            print(f"Created {output_file}")

if __name__ == "__main__":
    main()