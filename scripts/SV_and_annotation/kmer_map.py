#!/usr/bin/env python3

from Bio import SeqIO
import argparse

def reverse_complement(seq):
    """Generate reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

def build_ref_index(ref_seq, K = 20, REF_STEP = 1):
    """build kmer -> index dictionary"""
    index = {}

    for i in range(0, len(ref_seq) - K + 1, REF_STEP):
        kmer = ref_seq[i:i+K]
        index[kmer] = i + 1   # 1-based index

    return index

def calculate_optimal_k(ref_seq, min_k=20):
    """
    计算最佳K值，从min_k开始增加，直到build_ref_index得到的index没有重复的kmer
    """
    K = min_k
    while True:
        index = build_ref_index(ref_seq, K)
        # 计算理论上的kmer数量
        expected_kmers = len(ref_seq) - K + 1
        # 如果索引长度等于理论kmer数量，说明没有重复
        if len(index) == expected_kmers:
            return K
        # 否则增加K值
        K += 1
        # 防止无限循环
        if K > len(ref_seq):
            return min_k  # 如果所有K都有重复，返回最小K

def read_flank_regions(flank_file):
    """
    读取flank.bed文件，返回一个字典，键是seq_id，值是flank区域列表
    """
    flank_regions = {}
    with open(flank_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            seq_id = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            
            if seq_id not in flank_regions:
                flank_regions[seq_id] = []
            flank_regions[seq_id].append((start, end))
    return flank_regions

def is_in_flank(qpos, flank_regions, seq_id):
    """
    检查qpos是否在flank区域中
    """
    if seq_id not in flank_regions:
        return False
    
    for start, end in flank_regions[seq_id]:
        if start <= qpos <= end:
            return True
    return False

def query_sequence(seq, ref_index, QUERY_STEP = 10, K = 20):
    """scan query sequence and find matching reference index"""
    results = []

    for i in range(0, len(seq) - K + 1, QUERY_STEP):
        kmer = seq[i:i+K]

        if kmer in ref_index:
            results.append((i+1, ref_index[kmer]))
        
        else:
            rev_comp_kmer = reverse_complement(kmer)
            if rev_comp_kmer in ref_index:
                results.append((i+1, ref_index[rev_comp_kmer]))

    return results

def remove_redundant_matches(matches, bin_size):
    """
    移除冗余的匹配点，当某个点距离上一个保留的点小于bin且后边bin范围内有点，才去掉该点
    例如：1 2 3 4 5 6 100 这种情况下，保留1、6和100
    """
    if not matches:
        return []
    
    # 按qpos排序
    sorted_matches = sorted(matches, key=lambda x: x[0])
    
    non_redundant = [sorted_matches[0]]
    
    i = 1
    while True:
        current = sorted_matches[i]
        last = non_redundant[-1]
        
        # 检查当前点是否距离上一个保留的点小于bin且后面有在bin范围内的点
        if i == len(sorted_matches) - 1:
            non_redundant.append(current)
            break
        if current[0] - last[0] < bin_size and sorted_matches[i+1][0] - current[0] < bin_size:
            i += 1
            continue
        
        # 否则保留当前点
        non_redundant.append(current)
        i += 1
    
    return non_redundant


def main(args):
    ref_fasta = args.ref_fasta
    query_fasta = args.query_fasta
    out_file = args.out_file
    q_step = args.query_step
    flank_file = args.flank    
    K = args.kmer_size
    bin_size = args.bin_size

    ref_record = next(SeqIO.parse(ref_fasta, "fasta"))
    ref_seq = str(ref_record.seq).upper()

    if K is None:
        print("Calculating optimal K value...")
        k_value = calculate_optimal_k(ref_seq)
        print(f"Optimal K value: {k_value}")
        K = k_value

    REF_STEP = 1

    print("Building reference index...")
    ref_index = build_ref_index(ref_seq, K, REF_STEP)

    print("Reference fragments:", len(ref_index))

    flank_regions = {}
    if flank_file:
        print(f"Reading flank regions from {flank_file}...")
        flank_regions = read_flank_regions(flank_file)

    with open(out_file, "w") as out:

        out.write("query\tquery_pos\tref_index\n")

        for rec in SeqIO.parse(query_fasta, "fasta"):

            seq = str(rec.seq).upper()
            print(f"Processing {rec.id}...")
            matches = query_sequence(seq, ref_index, q_step, K)

            # 过滤flank区域内的点
            filtered_matches = []
            for qpos, ref_idx in matches:
                if not is_in_flank(qpos, flank_regions, rec.id):
                    filtered_matches.append((qpos, ref_idx))

            #如果bin_size不是1，移除冗余匹配点
            if bin_size and bin_size > 1:
                print(f"Removing redundant matches with bin size {bin_size}...")
                filtered_matches = remove_redundant_matches(filtered_matches, bin_size)

            # 写入结果
            for qpos, ref_idx in filtered_matches:
                out.write(f"{rec.id}\t{qpos}\t{ref_idx}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map query sequences to reference using k-mers")
    parser.add_argument("ref_fasta", help="Reference FASTA file")
    parser.add_argument("query_fasta", help="Query FASTA file")
    parser.add_argument("out_file", help="Output file")
    parser.add_argument("-s", "--query_step", type=int, default=10, help="Query step size (default: 10)")
    parser.add_argument("-f", "--flank", help="Flank regions BED file")
    parser.add_argument("-k", "--kmer_size", type=int, help="K-mer size (if not set, the kmer size will be calculated)")
    parser.add_argument("-b", "--bin_size", type=int, help="Bin size for removing redundant matches (default: 1)")

    args = parser.parse_args()
    main(args)
