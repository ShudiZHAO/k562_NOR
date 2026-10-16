#!/usr/bin/env python3
"""
处理fragment_align_info_formatted.xlsx文件，生成连接的cigar字段
"""

import pandas as pd
import re
from cigar_to_lengths import parse_cigar

def process_cigar(cigar, start, end):
    """
    处理cigar字符串：
    - 若start不为0，其cigar起始的“nS”应替换为“nI”
    - 若end不为-1，其cigar结尾的“nS”应替换为“nI”
    """
    # 处理起始的nS
    if start != 0:
        cigar = re.sub(r'^(\d+)S', r'\1I', cigar)
    
    # 处理结尾的nS
    if end != -1:
        cigar = re.sub(r'(\d+)S$', r'\1I', cigar)
    
    return cigar

def merge_consecutive_cigar(cigar):
    """
    合并连续的相同cigar操作，例如12M13M -> 25M
    """
    if not cigar:
        return cigar
    
    # 提取所有cigar操作和数字
    parts = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    if not parts:
        return cigar
    
    merged = []
    current_count = int(parts[0][0])
    current_op = parts[0][1]
    
    for i in range(1, len(parts)):
        count, op = parts[i]
        if op == current_op:
            current_count += int(count)
        else:
            merged.append(f"{current_count}{current_op}")
            current_count = int(count)
            current_op = op
    
    # 添加最后一个操作
    merged.append(f"{current_count}{current_op}")
    
    return ''.join(merged)

def extract_assembly(fragment_id):
    """
    从fragment_id中提取Assembly_x信息
    """
    # 假设fragment_id格式类似 "Assembly_1_0"
    match = re.match(r'Assembly_(\d+)_', fragment_id)
    if match:
        return f"Assembly_{match.group(1)}"
    return None

def main():
    # 读取xlsx文件
    input_file = 'fragment_align_info_formatted.xlsx'
    output_file = 'read_cigar_mapped.txt'
    
    try:
        df = pd.read_excel(input_file)
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    # 确保必要的列存在
    required_columns = ['fragment_id', 'read_id', 'start', 'end', 'cigar', 'aligned_pos', 'aligned_end']
    for col in required_columns:
        if col not in df.columns:
            print(f"Error: Column '{col}' not found in input file")
            return
    
    # 按read_id和start排序
    df = df.sort_values(by=['read_id', 'start'])
    
    # 分组处理每个read
    result = []
    for read_id, group in df.groupby('read_id'):
        # 确保按start排序
        group = group.sort_values(by='start')
        
        # 处理每个fragment的cigar
        processed_cigars = []
        prev_end = None
        min_start_aligned_pos = None
        assemblies = set()
        
        for _, row in group.iterrows():
            start = row['start']
            end = row['end']
            cigar = row['cigar']
            aligned_pos = row['aligned_pos']
            aligned_end = row['aligned_end']
            fragment_id = row['fragment_id']
            
            # 记录start最小的aligned_pos
            if min_start_aligned_pos is None:
                min_start_aligned_pos = aligned_pos

            # 提取并收集Assembly信息
            assembly = extract_assembly(fragment_id)
            if assembly:
                assemblies.add(assembly)

            # 处理cigar
            processed_cigar = process_cigar(cigar, start, end)
            processed_cigars.append(processed_cigar)
            
            # 检查是否需要插入D
            if prev_end is not None:
                gap = aligned_pos - prev_end -1
                if gap > 0:
                    processed_cigars.append(f"{gap}D")
                    print(f"Inserted gap of {gap} at read {read_id} between {prev_end} and {aligned_pos}")
            
            prev_end = aligned_end
        
        # 连接所有cigar
        combined_cigar = ''.join(processed_cigars)

        # 将开头和结尾的第一个nI替换为nS
        combined_cigar = re.sub(r'^(\d+)I', r'\1S', combined_cigar)  
        combined_cigar = re.sub(r'(\d+)I$', r'\1S', combined_cigar)

        # 合并连续的相同cigar操作
        combined_cigar = merge_consecutive_cigar(combined_cigar)

        # 调用parse_cigar函数计算query_length_from_cigar和ref_consumed
        query_length_from_cigar, ref_consumed = parse_cigar(combined_cigar)

        # 处理Assembly信息
        if len(assemblies) == 1:
            assembly = assemblies.pop()
        else:
            assembly = "ERROR: Multiple assemblies found"
            print(f"Warning: Read {read_id} has multiple assemblies: {assemblies}")

        result.append({
            'read_id': read_id,
            'assembly': assembly,
            'aligned_pos': min_start_aligned_pos,
            'query_length_from_cigar': query_length_from_cigar,
            'ref_consumed': ref_consumed,
            'combined_cigar': combined_cigar
        })
    
    # 生成输出文件
    output_df = pd.DataFrame(result)
    output_df.to_csv(output_file, sep='\t', index=False, header=False)
    
    print(f"Processed {len(result)} reads")
    print(f"Output written to {output_file}")

if __name__ == "__main__":
    main()