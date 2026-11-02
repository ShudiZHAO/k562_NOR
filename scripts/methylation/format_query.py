#!/usr/bin/env python3
"""
为 aligned_fragments 中的每个 fasta 文件处理序列，只添加 barcode

目录名为 NORx，barcode 文件表头也需要是 >NORx。
python format_query.py aligned_fragments/ assembly_final_with_empty.txt
"""

import os
import sys

from Bio import SeqIO

# 新的 array 目录名前缀
ARRAY_PREFIX = 'NOR'

# barcode 文件中的空占位符：占一个 junction 槽位，但长度为 0
EMPTY_BARCODE_TOKENS = ('', 'empty', 'nan', 'none')


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

def process_assembly_fragments(assembly_dir, assembly_name, barcodes, output_dir):
    """
    处理指定 array（NORx）的所有 fragment 序列
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
    
    # 为当前 array 创建输出目录
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
            
            # 如果readid的两个"_"之间的值不为0，在左端只添加barcode
            start_str = record.id.rsplit('_', 2)[1]
            if start_str != '0':
                barcode_index = i - 1
                if barcode_index < len(barcodes):
                    barcode = barcodes[barcode_index]
                    seq = barcode + seq
            
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
    
    # 创建输出目录
    output_dir = 'query_for_methylation'
    os.makedirs(output_dir, exist_ok=True)
    
    # 处理每个 array 目录（NORx）
    for assembly_dir in sorted(os.listdir(aligned_fragments_dir)):
        assembly_path = os.path.join(aligned_fragments_dir, assembly_dir)
        if not os.path.isdir(assembly_path):
            continue
        # 提取 array 名称
        if not assembly_dir.startswith(ARRAY_PREFIX):
            print(f"Warning: {assembly_dir} does not start with '{ARRAY_PREFIX}'")
            continue
        
        assembly_num = assembly_dir
        if assembly_num not in barcode_info:
            print(f"Warning: No barcode information found for {assembly_dir}")
            continue
        
        # 获取该 array 的 barcode 列表
        barcodes = barcode_info[assembly_num]
        
        # 处理该 array 的所有 fragment
        process_assembly_fragments(assembly_path, assembly_num, barcodes, output_dir)

if __name__ == "__main__":
    main()