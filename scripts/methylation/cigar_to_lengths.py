#!/usr/bin/env python3

import sys
import re
import os
from Bio import SeqIO

def parse_cigar(cigar):
    """
    根据 CIGAR 字符串计算比对的长度（所有 M 和 I 的碱基数之和）
    """
    length = 0
    ref_consumed = 0
    # 使用正则表达式匹配数字+字母的组合
    pattern = re.compile(r'(\d+)([MIDS])')
    for match in pattern.finditer(cigar):
        count = int(match.group(1))
        op = match.group(2)
        if op in ('M', 'I', 'S'):  # M 、S和 I 计入查询序列长度
            length += count
        if op in ('M', 'D'):
            ref_consumed += count
    return length, ref_consumed

def get_first_non_n_position(template_dir):
    """
    读取templete_for_methylation中的文件，返回一个字典{'NORx_y': 第一个不是N的字符出现的位置}
    当遇到NORx下y的最大值时，记录NORx_ymax的序列长度作为key NORx_ymax+1的值
    """
    result = {}
    array_info = {}  # 用于存储每个array（NORx）的y值和序列长度

    # 模板文件名格式：NOR1_0_consensus.fasta -> NOR1_0
    template_pattern = re.compile(r'^([^_]+)_(\d+)_consensus\.fasta$')

    # 第一遍遍历：收集所有array的信息
    for filename in os.listdir(template_dir):
        match = template_pattern.match(filename)
        if not match:
            continue

        # 从文件名中提取array信息和fragment编号
        # 文件名格式：NOR1_0_consensus.fasta -> NOR1_0
        print(f'Reading {filename} to get expected ref pos')
        array_name = match.group(1)
        fragment_num = int(match.group(2))
        key = f'{array_name}_{fragment_num}'

        # 读取fasta文件
        filepath = os.path.join(template_dir, filename)
        for record in SeqIO.parse(filepath, 'fasta'):
            seq = str(record.seq)
            # 找到第一个不是N的字符位置
            for i, base in enumerate(seq):
                if base != 'N':
                    result[key] = i
                    break

            # 存储array信息
            if array_name not in array_info:
                array_info[array_name] = []
            array_info[array_name].append((fragment_num, len(seq)))
            break  # 只处理第一个序列

    # 第二遍遍历：处理每个array的最大值情况
    for array_name, fragments in array_info.items():
        if fragments:
            # 找到y的最大值
            max_y, max_seq_length = max(fragments, key=lambda x: x[0])
            # 生成key NORx_ymax+1
            new_key = f'{array_name}_{max_y + 1}'
            # 将序列长度作为值
            result[new_key] = max_seq_length

    return result

def process_file(input_file, template_dir, output_file):
    """
    处理输入文件，计算比对长度和位置，写入输出文件。
    """
    ref_positions = get_first_non_n_position(template_dir)

    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        # 写入表头（可选）
        fout.write("fragment_id\tread_id\tstart\tend\tlength\taligned_length\texpected_ref_start\taligned_pos\texpected_ref_end\taligned_end\tcigar\n")

        for line in fin:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            fragment_id = parts[0]
            read_part = parts[1]
            aligned_pos = parts[2]
            seq_length = parts[3]
            cigar = parts[4]

            # 解析 read_id_start_end
            try:
                read_id, start_str, end_str = read_part.rsplit('_', 2)
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                # 格式不符合预期，跳过或报错
                continue

            # 计算 aligned_length
            aligned_length, ref_consumed = parse_cigar(cigar)

            # 计算 比对在参考序列上的起止位置
            # fragment_id 格式：NORx_<y>，例如 NOR10_0
            try:
                array_name, fragment_num_str = fragment_id.rsplit('_', 1)
                fragment_num = int(fragment_num_str)
            except ValueError:
                continue
            first_non_n_pos = ref_positions[f'{array_name}_{fragment_num}']
            expected_ref_start = first_non_n_pos + 1
            expected_ref_end = ref_positions[f'{array_name}_{fragment_num + 1}']
            aligned_end = int(aligned_pos) + ref_consumed -1
            if start == 0:
                expected_ref_start = expected_ref_end - ref_consumed + 1
            elif end == -1:
                expected_ref_end = expected_ref_start + ref_consumed - 1
                
                
            fout.write(f"{fragment_id}\t{read_id}\t{start}\t{end}\t{seq_length}\t{aligned_length}\t{expected_ref_start}\t{aligned_pos}\t{expected_ref_end}\t{aligned_end}\t{cigar}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python cigar_to_lengths.py input.txt template_dir output.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    template_dir = sys.argv[2]
    output_file = sys.argv[3]
    process_file(input_file, template_dir, output_file)
