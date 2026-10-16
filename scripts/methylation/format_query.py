#!/usr/bin/env python3
"""
为aligned_fragments中的每个fasta文件处理序列，添加prefix_seq、suffix_seq和barcode
python.exe format_query.py aligned_fragments\ assembly_for_link_consensus.txt
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

def process_assembly_fragments(assembly_dir, assembly_name, barcodes, prefix_seq, suffix_seq, output_dir):
    """
    处理指定Assembly的所有fragment序列
    """
    # 获取所有fasta文件并排序
    fasta_files = []
    for filename in os.listdir(assembly_dir):
        if filename.endswith('.fasta'):
            try:
                # 提取文件名中的数字
                seq_num = int(filename.split('.')[0])
                fasta_files.append((seq_num, filename))
            except:
                continue
    
    # 按数字排序
    fasta_files.sort()
    sorted_files = [filename for _, filename in fasta_files]
    
    # 为当前Assembly创建输出目录
    assembly_output_dir = os.path.join(output_dir, f"{assembly_name}")
    os.makedirs(assembly_output_dir, exist_ok=True)
    
    # 处理每个fasta文件
    for i, filename in enumerate(sorted_files):
        filepath = os.path.join(assembly_dir, filename)
        output_file = os.path.join(assembly_output_dir, f"formatted_{filename}")
        
        # 读取序列
        sequences = list(SeqIO.parse(filepath, 'fasta'))
        if not sequences:
            continue
        
        # 处理每个序列
        formatted_sequences = []
        for record in sequences:
            seq = str(record.seq)
            
            # 如果readid结尾不为-1，在右端添加prefix_seq
            if not record.id.endswith('_-1'):
                seq += prefix_seq
            
            # 如果readid的两个"_"之间的值不为0，在左端添加barcode和suffix_seq
            start_str= record.id.rsplit('_', 2)[1]
            if start_str != '0':
                barcode_index = i - 1
                if barcode_index < len(barcodes):
                    barcode = barcodes[barcode_index]
                    seq = barcode + suffix_seq + seq
            
            # 创建新的记录
            record.seq = record.seq.__class__(seq)
            formatted_sequences.append(record)
        
        # 写入输出文件
        with open(output_file, 'w') as f:
            SeqIO.write(formatted_sequences, f, 'fasta')
        
        print(f"Processed {filepath} -> {output_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python format_query.py <aligned_fragments_dir> <barcode_file>")
        sys.exit(1)
    
    aligned_fragments_dir = sys.argv[1]
    barcode_file = sys.argv[2]
    
    # 读取 barcode 信息
    barcode_info = read_barcode_info(barcode_file)
    
    # 定义 barcode 前后的序列
    prefix_seq = "TGGCGGAGATGGGCGCC"
    suffix_seq = "GCGAGG"
    
    # 创建输出目录
    output_dir = 'query_for_methylation'
    os.makedirs(output_dir, exist_ok=True)
    
    # 处理每个 Assembly 目录
    for assembly_dir in os.listdir(aligned_fragments_dir):
        assembly_path = os.path.join(aligned_fragments_dir, assembly_dir)
        if not os.path.isdir(assembly_path):
            continue
        # 提取 Assembly 编号
        if not assembly_dir.startswith('Assembly_'):
            print(f"Warning: {assembly_dir} does not start with 'Assembly_'")
            continue
        
        assembly_num = assembly_dir
        if assembly_num not in barcode_info:
            print(f"Warning: No barcode information found for {assembly_dir}")
            continue
        
        # 获取该 Assembly 的 barcode 列表
        barcodes = barcode_info[assembly_num]
        
        # 处理该 Assembly 的所有 fragment
        process_assembly_fragments(assembly_path, assembly_num, barcodes, prefix_seq, suffix_seq, output_dir)

if __name__ == "__main__":
    main()