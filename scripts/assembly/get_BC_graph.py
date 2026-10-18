#!/usr/bin/env python3

import argparse
from collections import Counter, defaultdict
from Bio import SeqIO

def load_barcode_replacements(replace_file):
    """从replace.txt文件加载barcode替换规则"""
    replacements = {}
    if replace_file:
        try:
            with open(replace_file, 'r') as f:
                for line in f:
                    # 跳过空行和注释行
                    if line.strip() == '' or line.startswith('#'):
                        continue
                    # 分割行，取前两列
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        original = parts[0]
                        replacement = parts[1]
                        replacements[original] = replacement
            print(f"[INFO] 从{replace_file}加载了{len(replacements)}条barcode替换规则")
        except Exception as e:
            print(f"[ERROR] 读取{replace_file}文件失败: {e}")
    return replacements

def parse_fastq(fq_file, left_flank, right_flank, replacements={}):
    """解析FASTQ，提取每条read的barcode路径"""
    read_barcodes = {}
    barcodes_pos = {}
    all_barcodes = []
    print("[INFO] 开始解析reads...")
    for i, record in enumerate(SeqIO.parse(fq_file, "fastq"), 1):
        seq = str(record.seq)
        barcodes = []
        poses = []
        start = 0
        while True:
            left_pos = seq.find(left_flank, start)
            if left_pos == -1:
                break
            right_pos = seq.find(right_flank, left_pos + len(left_flank))
            if right_pos == -1:
                break
            pos = left_pos + len(left_flank)
            barcode = seq[pos:right_pos]
            if barcode:
                if len(barcode) > 8:
                    barcode = barcode[:8]
                if barcode in replacements:
                    barcode = replacements[barcode]
                barcodes.append(barcode)
                poses.append(pos)
                all_barcodes.append(barcode)
            start = right_pos + len(right_flank)
        if barcodes:
            read_barcodes[record.id] = barcodes
            barcodes_pos[record.id] = poses
        if i % 50 == 0:
            print(f"[INFO] 已处理 {i} 条reads...")
    print(f"[INFO] 总共解析到 {len(read_barcodes)} 条含barcode的reads, 总barcode数={len(all_barcodes)}")
    return read_barcodes, all_barcodes, barcodes_pos

def filter_barcodes(read_barcodes, all_barcodes, min_count):
    """过滤低频barcode并根据上下文替换"""
    barcode_counts = Counter(all_barcodes)    
    # 创建所有read的barcode序列的列表，用于上下文搜索
    all_paths = list(read_barcodes.values())   
    # 存储需要替换的barcode及其可能的替换选项
    barcode_replacements = {}   
    # 找出低频barcode
    low_freq_barcodes = [b for b, c in barcode_counts.items() if (c <= min_count or len(b) != 6)]
    # 为每个低频barcode寻找可能的替换选项
    replaced = 0
    for bc in low_freq_barcodes:
        candidates = []        
        # 遍历所有read的barcode路径
        for path in all_paths:
            if bc in path:
                idx = path.index(bc)                
                # 1. 搜索该barcode左边的两个barcode顺序，得到右边的barcode
                if idx >= 2:
                    # 左边的两个barcode顺序
                    left_pair = (path[idx-2], path[idx-1])
                    # 在所有路径中搜索这个顺序并收集右边的barcode
                    for p in all_paths:
                        for i in range(len(p)-2):
                            if p[i] == left_pair[0] and p[i+1] == left_pair[1] and i+2 < len(p):
                                candidates.append(p[i+2])                
                # 2. 搜索该barcode右边的两个barcode顺序，得到左边的barcode
                if idx <= len(path)-3:
                    # 右边的两个barcode顺序
                    right_pair = (path[idx+1], path[idx+2])
                    # 在所有路径中搜索这个顺序并收集左边的barcode
                    for p in all_paths:
                        for i in range(1, len(p)-1):
                            if p[i] == right_pair[0] and p[i+1] == right_pair[1] and i-1 >= 0:
                                candidates.append(p[i-1])               
                # 3. 搜索该barcode左边和右边的barcode，将其中间的barcode放入临时列表
                if idx >= 1 and idx <= len(path)-2:
                    # 左边和右边的barcode
                    context = (path[idx-1], path[idx+1])
                    # 在所有路径中搜索这个上下文并收集中间的barcode
                    for p in all_paths:
                        for i in range(1, len(p)-1):
                            if p[i-1] == context[0] and p[i+1] == context[1]:
                                candidates.append(p[i])        
        # 如果有候选替换barcode，选择出现频率最高的
        if candidates:
            candidate_counts = Counter(candidates)
            # 找出出现频率最高的barcode
            most_common = candidate_counts.most_common(1)[0][0]
            print(candidate_counts.most_common)
            barcode_replacements[bc] = most_common
            if most_common != bc:
                replaced += 1
                print(f"[INFO] 替换低频barcode {bc} (出现{barcode_counts[bc]}次) 为 {most_common}")    
        else:
            print(f"[INFO] 未找到 {bc} 的替换选项，将被直接过滤")
    # 执行替换并过滤
    filtered_reads = {}
    for rid, path in read_barcodes.items():
        filtered = []
        for bc in path:
            if bc in barcode_replacements:
                # 替换低频barcode
                filtered.append(barcode_replacements[bc])
            elif barcode_counts[bc] > min_count:
                # 保留高频barcode
                filtered.append(bc)
            # 注意：没有足够上下文信息的低频barcode会被直接过滤掉   
        if filtered:
            filtered_reads[rid] = filtered   
    # 计算最终保留的barcode集合
    final_barcodes = []
    for path in filtered_reads.values():
        for bc in path:
            final_barcodes.append(bc)
    print(f"[INFO] 总barcode种类数={len(barcode_counts)}, 低频barcode个数={len(low_freq_barcodes)}, \
          替换的barcode个数={replaced}, 过滤后barcode个数={len(final_barcodes)}")   
    return filtered_reads, final_barcodes

def build_overlap_assemblies(read_paths):
    """用路径重叠方式合并barcode路径，避免缠绕"""
    # 按长度降序处理，优先长read
    assemblies = []
    sorted_paths = sorted(read_paths, key=lambda x: len(x), reverse=True)
    length_path = len(sorted_paths)
    print(f"[INFO] 开始路径拼接，总条数={len(sorted_paths)}")
    i = 0
    while i < length_path:
        max_overlap = 0
        best_link = None
        j = 0
        while j < length_path:
            if i == j:
                j += 1
                continue
            # 找最长 overlap
            overlap = 0
            for k in range(1, min(len(sorted_paths[i]), len(sorted_paths[j])) + 1):
                if sorted_paths[i][-k:] == sorted_paths[j][:k]:
                    overlap = k
            if overlap > max_overlap:
                max_overlap = overlap
                if max_overlap > 1:
                    best_link = sorted_paths[j]
            j += 1
        if best_link is not None:
            assembly = sorted_paths[i] + best_link[max_overlap:]
            assemblies.append(assembly)
            need_remove = sorted_paths[i]
            while need_remove in sorted_paths:
                sorted_paths.remove(need_remove)
            while best_link in sorted_paths:
                sorted_paths.remove(best_link)
            length_path = len(sorted_paths)
        i += 1
    print(f"[INFO] 拼接完成，共拼接 {len(assemblies)} 条assembly路径，剩余 {len(sorted_paths)} 条路径")
    assemblies += sorted_paths
    return assemblies

def remove_contained_paths(assemblies):
    # 添加功能：移除被其他路径包含的子路径
    print("[INFO] 开始检测并移除被包含的子路径...")
    to_keep = [True] * len(assemblies)
    # 检查每条路径是否被其他路径包含
    for i in range(len(assemblies)):
        if not to_keep[i]:  # 已经标记为删除，跳过
            continue        
        current_path = assemblies[i]
        current_len = len(current_path)        
        # 检查是否被其他更长的路径包含
        for j in range(len(assemblies)):
            if i == j or not to_keep[j]:  # 跳过自己和已经标记为删除的路径
                continue            
            other_path = assemblies[j]
            other_len = len(other_path)           
            # 如果当前路径比其他路径长，则不可能被包含，跳过
            if current_len >= other_len:
                continue           
            # 检查当前路径是否是其他路径的子序列
            # 将两个路径转换为字符串，使用特殊分隔符连接，避免部分匹配
            current_str = '|'.join(current_path)
            other_str = '|'.join(other_path)           
            if current_str in other_str:
                # 当前路径被其他路径包含，标记为删除
                to_keep[i] = False
                print(f"[INFO] 路径 {current_path} 被路径 {other_path} 包含，将被移除")
                break   
    # 根据标记保留路径
    filtered_assemblies = [assemblies[i] for i in range(len(assemblies)) if to_keep[i]]   
    # 再次按长度排序（可选）
    filtered_assemblies.sort(key=lambda x: len(x), reverse=True) 
    print(f"[INFO] 移除被包含路径后，剩余 {len(filtered_assemblies)} 条assembly路径")
    return filtered_assemblies

def align(filtered_assemblies, filtered_reads):
    alignments = {}
    str_assemblies = []
    for assembly in filtered_assemblies:
        str_assemblies.append('|'.join(assembly))
    for rid, path in filtered_reads.items():
        for i, assembly in enumerate(str_assemblies, 1):
            if '|'.join(path) in assembly:
                print(f"[INFO] 路径 {rid} 比对到 Assembly_{i}")
                index = int(assembly.find('|'.join(path))/7)
                alignments[rid] = [f'Assembly_{i}', index, index + len(path)]
                break
    print(f'[INFO] 共有 {len(filtered_reads)} 条路径，其中 {len(alignments)} 条比对到assembly中')
    return alignments

def write_origin(read_barcodes, out_file):
    with open(out_file, "w") as f:
        for rid, barcode in read_barcodes.items():
            f.write(f">{rid}\n{barcode}\n")
    print(f"[INFO] barcode原始文件输出: {out_file}")

def write_pos(barcodes_pos, out_file):
    with open(out_file, "w") as f:
        for rid, poses in barcodes_pos.items():
            f.write(f">{rid}\n")
            f.write(" ".join(map(str, poses)) + "\n")
    print(f"[INFO] barcode位置文件输出: {out_file}")

def write_paths_file(assemblies, out_file):
    with open(out_file, "w") as f:
        for i, path in enumerate(assemblies, 1):
            f.write(f">Assembly_{i}\n")
            f.write(" ".join(path) + "\n")
    print(f"[INFO] barcode路径文件输出: {out_file}")

def build_gfa(assemblies, out_gfa):
    """输出简单GFA图"""
    print(f"[INFO] 输出GFA文件: {out_gfa}")
    nodes = set()
    edges = set()
    for path in assemblies:
        nodes.update(path)
        for i in range(len(path)-1):
            edges.add((path[i], path[i+1]))
    with open(out_gfa, "w") as f:
        for node in nodes:
            f.write(f"S\t{node}\t*\n")
        for u, v in edges:
            f.write(f"L\t{u}\t+\t{v}\t+\t0M\n")

def write_alignments_file(alignments, out_file):
    with open(out_file, "w") as f:
        for rid, (asm, index, end_index) in alignments.items():
            f.write(f"{rid}\t{asm}\t{index}\t{end_index}\n")
    print(f"[INFO] 比对结果文件输出: {out_file}")

def get_assemble_poses(alignments, barcodes_pos, filtered_assemblies):
    """根据barcode路径拼接真实序列"""
    # 建立barcode -> 序列索引
    assemble_index = {}
    assemble_poses = {}
    assemblies = {}
    for i, path in enumerate(filtered_assemblies, 1):
        assemblies[f'Assembly_{i}'] = path
    for rid, (asm, index, end_index) in alignments.items():
        print(f"[INFO] reads {rid} 映射到 {asm} 位置 {index} 到 {end_index}")
        if asm not in assemble_index:
            assemble_index[asm] = []
        assemble_index[asm].append([rid, index, end_index])
    for asm in assemble_index:
    # 使用lambda函数按列表中的第二个元素（索引为1）排序
        assemble_index[asm].sort(key=lambda x: (x[1], -x[2]))
    for asm, aligned_reads in assemble_index.items():
        #去掉比对的开头相同的路径
        if len(aligned_reads) > 1:
            i = 1
            while i < len(aligned_reads):
                if aligned_reads[i][1] == aligned_reads[i-1][1]:
                    aligned_reads.remove(aligned_reads[i])
                else:
                    i += 1
        print(f"[INFO] {asm} 共有 {len(aligned_reads)} 条不同比对起始的路径")
        if len(aligned_reads) == 1:
            if aligned_reads[0][1] != 0 or aligned_reads[0][2] != len(assemblies[asm]):
                print(f'reads比对未能覆盖单路径 {asm} ，停止生成fasta')
                # return
                continue
            assemble_poses[asm] = [[aligned_reads[0][0], 0, 'end']]
        elif len(aligned_reads) > 1:
            if aligned_reads[0][1] != 0:
                print(f'reads比对未能覆盖 {asm} 开头，停止生成fasta')
                # return
                continue
            no_tail = True
            #第一个路径的end位置为第二个路径的start位置（aligned_reads[0][0]：当前路径rid，aligned_reads[1][1]：下一路径start的barcode index）
            assemble_poses[asm] = [[aligned_reads[0][0], 0, barcodes_pos[aligned_reads[0][0]][aligned_reads[1][1]]]]
            for i in range(1, len(aligned_reads)):
                if aligned_reads[i][2] == len(assemblies[asm]):
                    no_tail = False
                    assemble_poses[asm].append([aligned_reads[i][0], barcodes_pos[aligned_reads[i][0]][0], 'end'])
                    break
                if i == len(aligned_reads)-1:
                    print(f'reads比对未能覆盖 {asm} 结尾，停止生成fasta')
                    # return
                    continue
                if aligned_reads[i][2] < aligned_reads[i+1][1]:
                    print(f'{asm} 路径存在gap，停止生成fasta')
                    # return
                    continue
                assemble_poses[asm].append([aligned_reads[i][0], barcodes_pos[aligned_reads[i][0]][0], barcodes_pos[aligned_reads[i][0]][aligned_reads[i+1][1]-aligned_reads[i][1]]])
            if no_tail:
                print(f'reads比对未能覆盖 {asm} 结尾，停止生成fasta')
                # return
                continue
    # print(assemble_poses)
    return assemble_poses

def get_assemblied_fasta(fastq_file, assemble_poses, out_fasta):
    print(f"[INFO] 开始生成FASTA: {out_fasta}")
    fasta_input = {}
    fasta_output = {}
    for record in SeqIO.parse(fastq_file, "fastq"):
        rid = record.id
        seq = str(record.seq)
        fasta_input[rid] = seq
    for asm, reads_poses in assemble_poses.items():
        asm_seq = ''
        for read in reads_poses:
            if read[2] == 'end':
                asm_seq += fasta_input[read[0]][read[1]:]
            else:
                asm_seq += fasta_input[read[0]][read[1]:read[2]]
        fasta_output[asm] = asm_seq
    with open(out_fasta, "w") as f:
        for asm, asm_seq in fasta_output.items():
            f.write(f">{asm}\n")
            f.write(f"{asm_seq}\n")
    print(f"[INFO] FASTA输出完成: {out_fasta}")

def write_barcode_file(read_barcodes, out_file):
    sorted_reads = sorted(read_barcodes.items(), key=lambda item: len(item[1]), reverse=True)
    with open(out_file, "w") as f:
        for key,seq in sorted_reads:
            f.write(f">{key}\n")
            f.write(" ".join(seq) + "\n")
    print(f"[INFO] barcode路径文件输出: {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Nanopore barcode-guided assembly (overlap merge)")
    parser.add_argument("-i", "--input", required=True, help="输入FASTQ文件")
    parser.add_argument("--left", required=True, help="barcode左侧固定序列")
    parser.add_argument("--right", required=True, help="barcode右侧固定序列")
    parser.add_argument("--min_count", type=int, default=3, help="过滤低频barcode")
    parser.add_argument("--paths_out", default="assembly_paths.txt", help="输出barcode路径文件")
    parser.add_argument("--barcodes_out", default="barcode_paths.txt", help="输出barcode文件")
    parser.add_argument("--gfa_out", default="assembly.gfa", help="输出GFA文件")
    parser.add_argument("--fasta_out", default="assembly.fasta", help="输出FASTA文件")
    parser.add_argument("--replace_file", help="barcode替换规则文件路径(replace.txt)")
    parser.add_argument("--alignments_out", default="alignments.xls", help="输出比对结果文件")
    args = parser.parse_args()

    replacements = load_barcode_replacements(args.replace_file)
    read_barcodes, all_barcodes, barcodes_pos = parse_fastq(args.input, args.left, args.right, replacements)
    filtered_reads, valid_barcodes = filter_barcodes(read_barcodes, all_barcodes, args.min_count)
    write_barcode_file(filtered_reads, args.barcodes_out)
    write_pos(barcodes_pos, args.barcodes_out.replace('.txt', '_pos.txt'))
    while valid_barcodes != all_barcodes:
        print('[INFO] barcode列表已更新，重新校正和过滤')
        all_barcodes = valid_barcodes
        read_barcodes = filtered_reads
        filtered_reads, valid_barcodes = filter_barcodes(read_barcodes, all_barcodes, args.min_count)
    assembly_input = filtered_reads.values()
    assemblies = build_overlap_assemblies(assembly_input)
    while len(assemblies) != len(assembly_input):
        assembly_input = assemblies
        assemblies = build_overlap_assemblies(assembly_input)
    assemblies = remove_contained_paths(assemblies)
    alignments = align(assemblies, filtered_reads)
    write_paths_file(assemblies, args.paths_out)
    build_gfa(assemblies, args.gfa_out)
    write_alignments_file(alignments, args.alignments_out)
    assemble_poses = get_assemble_poses(alignments, barcodes_pos, assemblies)
    get_assemblied_fasta(args.input, assemble_poses, args.fasta_out)
    print("[INFO] 任务完成！")

if __name__ == "__main__":
    main()