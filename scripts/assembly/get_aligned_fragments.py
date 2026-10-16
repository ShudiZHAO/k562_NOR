#!/usr/bin/env python3
"""
根据xlsx文件中的barcode比对信息，提取每个assembly中每个barcode位置的序列
"""

import pandas as pd
import os
from Bio import SeqIO
from collections import defaultdict
import sys

if len(sys.argv) != 4:
    print("Usage: python get_aligned_fragments.py <input.xlsx> <fasta> <output_dir>")
    sys.exit(1)

def parse_xlsx_and_generate_fasta(xlsx_file, fasta_file, output_dir):
    """
    xlsx_file: 包含barcode比对信息的Excel文件
    fasta_file: 包含所有reads序列的fasta文件
    output_dir: 输出目录
    """
    
    # 读取xlsx文件
    df = pd.read_excel(xlsx_file)
    
    # 首先读取所有reads的序列到字典中
    print(f"reading fasta file: {fasta_file}")
    read_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        read_sequences[record.id] = str(record.seq)
    print(f"total {len(read_sequences)} reads sequences")
    
    # 按assembly分组
    assemblies = df.groupby('assembly')
    
    # 为每个assembly创建文件夹
    for assembly_name, assembly_df in assemblies:
        print(f"\nprocessing assembly: {assembly_name}")
        
        # 创建assembly文件夹
        assembly_dir = os.path.join(output_dir, assembly_name)
        os.makedirs(assembly_dir, exist_ok=True)
        
        # 确定这个assembly中最大的barcode索引
        # 根据start_index和end_index，每个read覆盖的barcode范围
        # 需要找出所有可能出现的barcode位置索引
        max_barcode_idx = 0
        for _, row in assembly_df.iterrows():
            max_barcode_idx = max(max_barcode_idx, row['end_index'])
        
        print(f"total barcode index range: 0-{max_barcode_idx}")
        
        # 初始化每个barcode位置的序列列表
        barcode_sequences = defaultdict(list)
        
        # 处理每条read
        for _, row in assembly_df.iterrows():
            read_id = row['read_id']
            start_idx = row['start_index']
            end_idx = row['end_index']
            # barcode_num = row['end_index'] - row['start_index'] + 1  # 实际barcode长度
            barcode_len = [len(x) for x in str(row['barcodes']).split('|')]
            intervals = [int(x) for x in str(row['intervals']).split('|')]
            first_barcode_start = row['first_barcode_to_start']
            
            # 获取该read的完整序列
            if str(read_id) not in read_sequences:
                print(f"warning: read {read_id} not found in fasta file")
                continue
            
            full_seq = read_sequences[str(read_id)]
            # read_len = len(full_seq)
            
            # 计算每个barcode在read上的起始位置
            # 第一个barcode起始位置
            barcode_starts = [first_barcode_start]
            
            # 计算后续barcode的起始位置
            for i in range(len(intervals)-1):
                barcode_starts.append(barcode_starts[-1] + intervals[i+1])
            
            # 根据start_index和end_index确定这个read覆盖的barcode范围
            # 注意：start_index和end_index是barcode在assembly上的索引
            # 例如start_index=0, end_index=3 表示覆盖barcode 0,1,2,3
            fragment_start = 0
            for i in range(start_idx, end_idx + 1):
                # 计算这个barcode在read上的索引
                # 第一个barcode在read上的索引是0
                read_barcode_idx = i - start_idx
                
                if read_barcode_idx < len(barcode_starts):
                    # 该barcode在read上的起始位置
                    fragment_end = barcode_starts[read_barcode_idx]
                else:
                    # 最后一个barcode，结束位置是read结尾减去last_barcode_to_end
                    fragment_end = -1
                    
                    # 提取该barcode对应的序列
                    # 注意：barcode_start是read上的位置，barcode_end也是read上的位置
                    # 序列从bc_start到bc_end（不包括bc_end，保持python的切片习惯）
                if fragment_start == 0 and fragment_end != -1:
                    fragment_seq = full_seq[fragment_start:fragment_end-17]
                elif fragment_start != 0 and fragment_end == -1:
                    fragment_seq = full_seq[fragment_start+6:fragment_end]
                elif fragment_start != 0 and fragment_end != -1:
                    fragment_seq = full_seq[fragment_start+6:fragment_end-17]
                else:
                    print(f"warning: unknown case for read {read_id}, barcode {i}")

                # 记录这个序列
                barcode_sequences[i].append({
                    'read_id': read_id,
                    'sequence': fragment_seq,
                    'start': fragment_start,
                    'end': fragment_end
                })
                if fragment_end == -1:
                    break
                fragment_start = fragment_end + barcode_len[read_barcode_idx]
        
        # 为每个barcode位置生成fasta文件
        for barcode_idx, sequences in barcode_sequences.items():
            output_file = os.path.join(assembly_dir, f"{barcode_idx}.fasta")
            
            with open(output_file, 'w') as f:
                for seq_info in sequences:
                    # fasta格式：>read_id_start_end
                    header = f">{seq_info['read_id']}_{seq_info['start']}_{seq_info['end']}"
                    f.write(header + '\n')
                    f.write(seq_info['sequence'] + '\n')
            
            print(f"{len(sequences)} sequences generated in {output_file}")
        
        # 检查是否有缺失的barcode索引
        all_possible_indices = set(range(max_barcode_idx + 1))
        found_indices = set(barcode_sequences.keys())
        missing_indices = all_possible_indices - found_indices
        if missing_indices:
            print(f"warning: these barcode index has no sequences: {sorted(missing_indices)}")

def main():
    # 配置参数
    xlsx_file = sys.argv[1]  # 替换为你的xlsx文件路径
    fasta_file = sys.argv[2]  # 你的fasta文件路径
    output_dir = sys.argv[3]  # 输出目录
    
    # 运行处理
    parse_xlsx_and_generate_fasta(xlsx_file, fasta_file, output_dir)
    print("\nDone!")

if __name__ == "__main__":
    main()