#!/usr/bin/env python3
"""
连接 consensus 序列，使用指定的 barcode
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
    # print(barcode_info)
    return barcode_info

def get_fragment_sequences(consensus_dir, assembly):
    """
    获取指定 assembly 的所有 fragment 序列
    返回一个字典：{y: sequence}
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
                sequences[y] = str(record.seq)
                break  # 只取第一个序列
    
    return sequences

def link_sequences(sequences, barcodes, prefix_seq, suffix_seq):
    """
    连接序列，使用处理后的 barcode
    """
    # 按 y 值排序
    sorted_ys = sorted(sequences.keys())
    if not sorted_ys:
        return ""
    
    # 开始连接
    linked_seq = sequences[sorted_ys[0]]
    
    for i in range(1, len(sorted_ys)):
        y = sorted_ys[i]
        # 获取对应的 barcode
        barcode_idx = i - 1
        if barcode_idx < len(barcodes):
            barcode = barcodes[barcode_idx]
            # 添加处理后的 barcode
            linked_seq += prefix_seq + barcode + suffix_seq
        # 添加下一个序列
        linked_seq += sequences[y]
    
    return linked_seq

def main():
    if len(sys.argv) != 3:
        print("Usage: python link_consensus.py <consensus_dir> <barcode_file>")
        sys.exit(1)
    
    consensus_dir = sys.argv[1]
    barcode_file = sys.argv[2]
    
    # 读取 barcode 信息
    barcode_info = read_barcode_info(barcode_file)
    
    # 定义 barcode 前后的序列
    prefix_seq = "TGGCGGAGATGGGCGCC"
    suffix_seq = "GCGAGG"
    
    # 处理每个 assembly
    os.makedirs('linked_consensus', exist_ok=True)
    for assembly, barcodes in barcode_info.items():
        # 获取 fragment 序列
        sequences = get_fragment_sequences(consensus_dir, assembly)
        
        if not sequences:
            print(f"Warning: No fragments found for Assembly_{assembly}")
            continue
        
        # 连接序列
        linked_seq = link_sequences(sequences, barcodes, prefix_seq, suffix_seq)
        
        # 输出结果
        output_file = f"linked_consensus/{assembly}.fasta"
        with open(output_file, 'w') as f:
            f.write(f">{assembly}\n")
            # 每行 60 个碱基
            for i in range(0, len(linked_seq), 60):
                f.write(linked_seq[i:i+60] + "\n")
        
        print(f"Created {output_file}")

if __name__ == "__main__":
    main()