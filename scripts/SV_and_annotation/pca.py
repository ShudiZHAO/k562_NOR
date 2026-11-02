#!/usr/bin/env python3
"""
对BED文件中的unit取序列，用k-mer频率做PCA，输出PCA结果CSV（不画图）。

BED第四列不是 NORx_y 命名的单元（如 Right_fragment、Left_fragment）不参与分析，直接跳过。

用法:
    python pca.py --bed named_units.bed --fasta NOR_k562.fa --k 6 \
                  --components 2 --clusters 7 --output pca_result

输入:
    --bed     BED文件，至少4列: chrom / start / end / name（0-based半开区间）
    --fasta   基因组FASTA，序列名与BED第1列一致

输出:
    <output>_components.csv      PC1..PCn、ID、Cluster
    <output>_kmer_frequency.csv  各k-mer频率（加 --save_kmer 时）
    pca_model.pkl               StandardScaler + PCA模型（可用 --model_output 改）
"""

import pickle
import re
import sys
from itertools import product

import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# unit名 NORx_y，例如 NOR11_1 -> NOR编号11、unit序号1
NOR_UNIT_RE = re.compile(r'^NOR(\d+)_(\d+)$', re.IGNORECASE)


def is_nor_unit(name):
    """BED第四列是否为 NORx_y 命名（只有这类unit参与分析）"""
    return NOR_UNIT_RE.match(name) is not None


def read_bed_file(bed_path):
    """
    读取BED文件，返回unit信息（chrom/start/end/name），第4列之后的内容忽略
    BED坐标为0-based半开区间 [start, end)
    """
    units = []
    with open(bed_path, 'r') as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith(('#', 'track', 'browser')):
                continue
            parts = line.split()
            if len(parts) < 4:
                print(f"Warning: line {line_no} skipped (need >=4 columns): {line}")
                continue
            units.append({
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'name': parts[3],
            })
    return units


def read_fasta_file(fasta_path):
    """
    读取FASTA文件，返回 {序列名: 大写序列}
    使用Bio.SeqIO读取
    """
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences


def get_unit_sequences(units, fasta_dict):
    """
    根据unit信息从fasta字典中获取对应的序列
    BED文件中的位置是0-based，半开区间 [start, end)
    """
    unit_names = []
    unit_sequences = []

    for unit in units:
        chrom = unit['chrom']
        start = unit['start']
        end = unit['end']
        name = unit['name']

        if chrom not in fasta_dict:
            print(f"Warning: chromosome {chrom} not found in FASTA file (unit {name})")
            continue

        chrom_seq = fasta_dict[chrom]
        if start < 0 or end > len(chrom_seq) or start >= end:
            print(f"Warning: unit {name} range out of bounds for {chrom} "
                  f"(length {len(chrom_seq)}, range [{start}, {end}))")
            continue

        unit_sequences.append(chrom_seq[start:end])
        unit_names.append(name)

    return unit_names, unit_sequences


def kmer_encode(sequences, k=6):
    """
    对DNA序列进行k-mer频率编码
    计算每个序列中所有可能的k-mer（4^k维）的频率，行内归一化到和为1
    """
    kmers = [''.join(p) for p in product('ACGT', repeat=k)]
    kmer_to_idx = {kmer: i for i, kmer in enumerate(kmers)}

    encoded = np.zeros((len(sequences), len(kmers)))
    for i, seq in enumerate(sequences):
        for j in range(len(seq) - k + 1):
            kmer = seq[j:j + k]
            if kmer in kmer_to_idx:
                encoded[i, kmer_to_idx[kmer]] += 1
        # 归一化
        if np.sum(encoded[i]) > 0:
            encoded[i] /= np.sum(encoded[i])

    return encoded, kmers


def perform_pca(features, n_components=2):
    """
    标准化后执行PCA分析
    """
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)

    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(scaled_features)

    # 解释方差比例
    explained_variance = pca.explained_variance_ratio_

    return principal_components, explained_variance, scaler, pca


def save_pca_model(scaler, pca, params, model_path='pca_model.pkl'):
    """
    保存PCA模型和scaler到文件
    """
    model = {
        'scaler': scaler,
        'pca': pca,
        'params': params  # 保存k值等参数
    }
    with open(model_path, 'wb') as f:
        pickle.dump(model, f)
    print(f"PCA model saved to {model_path}")


def main():
    """
    主函数
    """
    import argparse

    parser = argparse.ArgumentParser(
        description='PCA analysis of DNA units using k-mer frequency encoding')
    parser.add_argument('--bed', required=True, help='Path to BED file containing unit information')
    parser.add_argument('--fasta', required=True, help='Path to FASTA file containing DNA sequences')
    parser.add_argument('--k', type=int, default=6, help='k-mer length (default: 6)')
    parser.add_argument('--components', type=int, default=2, help='Number of PCA components')
    parser.add_argument('--clusters', type=int, default=3, help='Number of clusters for K-means')
    parser.add_argument('--output', default='pca_result', help='Output file prefix for PCA result csv')
    parser.add_argument('--save_kmer', action='store_true', help='Save kmer frequency data')
    parser.add_argument('--model_output', default='pca_model.pkl', help='Path to save PCA model')
    args = parser.parse_args()

    if args.k < 1:
        parser.error('--k must be >= 1')
    if args.components < 1:
        parser.error('--components must be >= 1')
    if args.clusters < 1:
        parser.error('--clusters must be >= 1')

    # 读取BED文件
    units = read_bed_file(args.bed)
    print(f"Read {len(units)} units from BED file")
    if not units:
        sys.exit(f"Error: no unit found in {args.bed}, check the BED format (need >=4 columns)")

    # 只分析BED第四列以NOR开头的unit（Right_fragment / Left_fragment 等跳过）
    skipped = [u['name'] for u in units if not is_nor_unit(u['name'])]
    if skipped:
        skipped_counts = pd.Series(skipped).value_counts()
        print(f"Skipped {len(skipped)} units not named NORx_y: "
              + ', '.join(f'{n}({c})' for n, c in skipped_counts.items()))
    units = [u for u in units if is_nor_unit(u['name'])]
    print(f"Kept {len(units)} units named NORx_y")
    if not units:
        sys.exit('Error: no unit named NORx_y found in the BED file (expected names like NOR11_1)')

    # 读取FASTA文件（使用Bio.SeqIO）
    fasta_dict = read_fasta_file(args.fasta)
    print(f"Read {len(fasta_dict)} sequences from FASTA file")

    # 获取unit对应的序列
    unit_names, unit_sequences = get_unit_sequences(units, fasta_dict)
    print(f"Found {len(unit_sequences)} sequences for units")
    if not unit_sequences:
        sys.exit('Error: no unit sequence extracted, '
                 'check that column 1 of the BED matches the FASTA sequence names')
    if args.clusters > len(unit_sequences):
        parser.error(f'--clusters ({args.clusters}) cannot exceed '
                     f'the number of units ({len(unit_sequences)})')

    # k-mer频率编码
    features, kmers = kmer_encode(unit_sequences, k=args.k)
    print(f"k-mer encoding (k={args.k}): {features.shape[1]} features")
    if args.components > min(features.shape):
        parser.error(f'--components ({args.components}) cannot exceed '
                     f'min(n_units, n_features) = {min(features.shape)}')

    if args.save_kmer:
        kmer_df = pd.DataFrame(features, columns=kmers)
        kmer_df.insert(0, 'unit_index', unit_names)
        kmer_df.to_csv(f'{args.output}_kmer_frequency.csv', index=False)
        print(f"K-mer frequency data saved to {args.output}_kmer_frequency.csv")

    # 执行PCA
    principal_components, explained_variance, scaler, pca = perform_pca(features, args.components)
    print(f"PCA completed with {args.components} components")
    print(f"Explained variance: {explained_variance}")

    # 保存PCA模型
    save_pca_model(scaler, pca, {'k': args.k}, args.model_output)

    # 保存PCA结果和聚类信息
    out_csv = f'{args.output}_components.csv'
    pca_df = pd.DataFrame(principal_components,
                          columns=[f'PC{i+1}' for i in range(args.components)])
    pca_df['ID'] = unit_names

    # 执行K-means聚类以保存聚类信息
    kmeans = KMeans(n_clusters=args.clusters, random_state=42)
    pca_df['Cluster'] = kmeans.fit_predict(principal_components)

    pca_df.to_csv(out_csv, index=False)
    print(f"PCA components and cluster information saved to {out_csv}")


if __name__ == "__main__":
    main()

