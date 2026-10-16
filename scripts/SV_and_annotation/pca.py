import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from Bio import SeqIO
import pickle

def read_bed_file(bed_path):
    """
    读取BED文件，返回unit信息
    """
    units = []
    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3]
                # 只保留以"Array"开头的unit
                if name.startswith('Array') and not name.startswith('Array_1_'):
                    units.append({'chrom': chrom, 'start': start, 'end': end, 'name': name})
    return units

def read_fasta_file(fasta_path):
    """
    读取FASTA文件，返回序列字典
    使用Bio.SeqIO读取
    """
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences

def read_array_nor_mapping(mapping_path):
    """
    读取name_of_the_arrays.txt文件，返回Array到NOR的映射字典
    """
    array_to_nor = {}
    with open(mapping_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                array_name = parts[1]  # Array_1, Array_2, etc.
                nor_name = parts[2]    # NOR 1, NOR 10, etc.
                # 移除空格，转换为NOR1格式
                nor_name = nor_name.replace(' ', '')
                array_to_nor[array_name] = nor_name
    return array_to_nor

def get_unit_sequences(units, fasta_dict):
    """
    根据unit信息从fasta字典中获取对应的序列
    BED文件中的位置是0-based，半开区间 [start, end)
    """
    unit_sequences = []
    unit_names = []
    
    for unit in units:
        chrom = unit['chrom']
        start = unit['start']
        end = unit['end']
        name = unit['name']
        
        if chrom in fasta_dict:
            chrom_seq = fasta_dict[chrom]
            # 确保位置在染色体长度范围内
            if start < len(chrom_seq) and end <= len(chrom_seq):
                # 提取子序列（BED文件是0-based，Python字符串也是0-based）
                unit_seq = chrom_seq[start:end]
                unit_sequences.append(unit_seq)
                unit_names.append(name)
            else:
                print(f"Warning: Unit {name} position out of range for chromosome {chrom}")
                print(f"Chromosome length: {len(chrom_seq)}, Unit range: [{start}, {end})")
        else:
            print(f"Warning: Chromosome {chrom} not found in FASTA file")
    
    return unit_names, unit_sequences

def one_hot_encode(sequences):
    """
    对DNA序列进行独热编码
    每个核苷酸被编码为4维向量：A=[1,0,0,0], C=[0,1,0,0], G=[0,0,1,0], T=[0,0,0,1]
    对于不同长度的序列，将使用最长序列的长度， shorter sequences will be padded with zeros
    """
    # 确定最长序列长度
    max_len = max(len(seq) for seq in sequences)
    
    # 初始化编码矩阵
    encoded = np.zeros((len(sequences), max_len * 4))
    
    # 核苷酸到索引的映射
    nuc_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    for i, seq in enumerate(sequences):
        for j, nuc in enumerate(seq):
            if nuc in nuc_map:
                encoded[i, j*4 + nuc_map[nuc]] = 1
    
    return encoded

def kmer_encode(sequences, k=3):
    """
    对DNA序列进行k-mer频率编码
    计算每个序列中所有可能的k-mer的频率
    """
    # 生成所有可能的k-mer
    from itertools import product
    kmers = [''.join(p) for p in product('ACGT', repeat=k)]
    kmer_to_idx = {kmer: i for i, kmer in enumerate(kmers)}
    
    # 初始化编码矩阵
    encoded = np.zeros((len(sequences), len(kmers)))
    
    for i, seq in enumerate(sequences):
        for j in range(len(seq) - k + 1):
            kmer = seq[j:j+k]
            if kmer in kmer_to_idx:
                encoded[i, kmer_to_idx[kmer]] += 1
        
        # 归一化
        if np.sum(encoded[i]) > 0:
            encoded[i] /= np.sum(encoded[i])
    
    return encoded, kmers

def align_encode(sequences, ref_path, k=20):
    """
    对DNA序列进行基于k-mer比对的编码
    使用kmer_map.py中的函数构建索引并查询序列
    """
    # 读取参考序列
    ref_record = next(SeqIO.parse(ref_path, "fasta"))
    ref_seq = str(ref_record.seq).upper()
    
    # 构建参考索引
    ref_index = kmer_map.build_ref_index(ref_seq, K=k)
    print(f"Built reference index with {len(ref_index)} kmers")
    
    # 对每个序列进行查询
    max_matches = 0
    all_matches = []
    
    for seq in sequences:
        matches = kmer_map.query_sequence(seq, ref_index, QUERY_STEP=1, K=k)
        all_matches.append(matches)
        if len(matches) > max_matches:
            max_matches = len(matches)
    
    # 构建编码矩阵
    # 对于每个序列，我们记录其匹配的参考位置
    # 为了统一长度，我们只取前max_matches个匹配，不足的用0填充
    encoded = np.zeros((len(sequences), max_matches))
    
    for i, matches in enumerate(all_matches):
        for j, (_, ref_pos) in enumerate(matches[:max_matches]):
            encoded[i, j] = ref_pos
    
    return encoded

def nucleotide_composition(sequences):
    """
    基于核苷酸组成的编码
    计算每个序列中A、C、G、T的比例
    """
    encoded = np.zeros((len(sequences), 4))
    nuc_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    for i, seq in enumerate(sequences):
        counts = np.zeros(4)
        for nuc in seq:
            if nuc in nuc_map:
                counts[nuc_map[nuc]] += 1
        
        # 归一化
        if len(seq) > 0:
            counts /= len(seq)
        
        encoded[i] = counts
    
    return encoded

def perform_pca(features, n_components=2):
    """
    执行PCA分析
    """
    # 标准化特征
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)
    
    # 执行PCA
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(scaled_features)
    
    # 计算解释方差比例
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

def plot_pca_by_cluster(principal_components, ids, explained_variance, n_clusters=3):
    """
    绘制PCA结果，按cluster配色，不添加文字和网格线
    """
    plt.figure(figsize=(12, 10))
    
    # 执行K-means聚类
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(principal_components)
    
    # 为每个聚类绘制点和置信圈
    colors = ['#DE77AE','#7FBC41','#9970AB','#D6604D','#4393C3','#8C510A','#35978F']
    for i in range(n_clusters):
        # 获取当前聚类的点
        cluster_points = principal_components[clusters == i]
        
        if len(cluster_points) > 0:
            # 计算聚类中心
            center = np.mean(cluster_points, axis=0)
            
            # 计算协方差矩阵
            cov = np.cov(cluster_points.T)
            
            # 计算置信椭圆
            from matplotlib.patches import Ellipse
            import math
            
            # 计算特征值和特征向量
            eigenvalues, eigenvectors = np.linalg.eigh(cov)
            
            # 确保特征值非负（处理浮点精度问题）
            eigenvalues = np.maximum(eigenvalues, 0)
            
            # 找到最大和最小特征值的索引
            order = eigenvalues.argsort()[::-1]
            eigenvalues = eigenvalues[order]
            eigenvectors = eigenvectors[:, order]
            
            # 计算椭圆的角度
            angle = math.degrees(math.atan2(eigenvectors[1, 0], eigenvectors[0, 0]))
            
            # 计算椭圆的宽度和高度（2倍标准差）
            width = 2 * math.sqrt(eigenvalues[0]) * 2 if eigenvalues[0] > 0 else 0
            height = 2 * math.sqrt(eigenvalues[1]) * 2 if eigenvalues[1] > 0 else 0
            
            # 只有当宽度和高度都大于0时才绘制椭圆
            if width > 0 and height > 0:
                # 创建椭圆
                ellipse = Ellipse(xy=center, width=width, height=height, angle=angle,
                                edgecolor=colors[i % len(colors)], facecolor='none', 
                                alpha=0.5, linewidth=2)
                
                # 添加椭圆到图中
                plt.gca().add_patch(ellipse)
            
            # 绘制聚类中的点
            plt.scatter(cluster_points[:, 0], cluster_points[:, 1], 
                        c=colors[i % len(colors)], label=f'Cluster {i+1}', alpha=0.8)
    
    # 设置标题和轴标签
    plt.title('PCA Analysis of DNA Units with Clusters')
    plt.xlabel(f'PC1 ({explained_variance[0]:.2%} variance)')
    plt.ylabel(f'PC2 ({explained_variance[1]:.2%} variance)')
    
    # 添加图例
    plt.legend()

    # 获取当前图的范围
    xlim = plt.gca().get_xlim()
    ylim = plt.gca().get_ylim()

    return plt, xlim, ylim

def plot_pca_by_array(principal_components, ids, explained_variance, xlim, ylim, array_to_nor=None):
    """
    绘制PCA结果，按array号配色，不添加文字和网格线
    """
    plt.figure(figsize=(12, 10))

    nor_colors = {
        'NOR1': '#8C510A',
        'NOR2': '#4D4D4D',
        'NOR3': '#35978F',
        'NOR4': '#80CDC1',
        'NOR5': '#C51B7D',
        'NOR6': '#DE77AE',
        'NOR7': '#F1B6DA',
        'NOR8': '#4D9221',
        'NOR9': '#7FBC41',
        'NOR10': '#B8E186',
        'NOR11': '#762A83',
        'NOR12': '#9970AB',
        'NOR13': '#D6604D',
        'NOR14': '#F4A582',
        'NOR15': '#2166AC',
        'NOR16': '#4393C3',
        'NOR17': '#6BACD1',
        'NOR18': '#92C5DE'
    }

    # 从ID中提取array号，并映射到NOR
    nor_ids = []
    for id in ids:
        # 提取完整的Array名称，例如Array_3_1 -> Array_3
        parts = id.split('_')
        if len(parts) >= 2 and parts[0] == 'Array':
            array_name = f'Array_{parts[1]}'
            # 如果有映射文件，使用NOR名称
            if array_to_nor and array_name in array_to_nor:
                nor_ids.append(array_to_nor[array_name])
            else:
                nor_ids.append(f'Array_{parts[1]}')
        else:
            nor_ids.append('Unknown')
    
    # 获取唯一的NOR号并排序
    unique_nors = sorted(list(set(nor_ids)), key=lambda x: (int(x.replace('NOR', '').replace('Array_', '')) if x.replace('NOR', '').replace('Array_', '').isdigit() else 999))
    
    # 为每个NOR绘制点，确保每个NOR只添加一个图例项
    for nor_id in unique_nors:
        # 获取当前NOR的所有点
        nor_points = principal_components[np.array(nor_ids) == nor_id]
        if len(nor_points) > 0:
            # 使用用户指定的颜色，如果没有则使用默认颜色
            color = nor_colors.get(nor_id, '#808080')  # 灰色作为默认
            plt.scatter(nor_points[:, 0], nor_points[:, 1], 
                        c=[color], label=f'{nor_id}', alpha=0.6)
    
    # 设置与cluster图相同的范围
    plt.xlim(xlim)
    plt.ylim(ylim)
    
    # 设置标题和轴标签
    plt.title('PCA Analysis of DNA Units by NOR')
    plt.xlabel(f'PC1 ({explained_variance[0]:.2%} variance)')
    plt.ylabel(f'PC2 ({explained_variance[1]:.2%} variance)')
    
    # 添加图例
    plt.legend()
    
    return plt

def main():
    """
    主函数
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='PCA analysis of DNA sequences from BED and FASTA files')
    parser.add_argument('--bed', required=True, help='Path to BED file containing unit information')
    parser.add_argument('--fasta', required=True, help='Path to FASTA file containing DNA sequences')
    parser.add_argument('--method', choices=['onehot', 'kmer', 'composition', 'align'], 
                        default='kmer', help='Encoding method for DNA sequences')
    parser.add_argument('--k', type=int, help='k-mer length for kmer or align encoding')
    parser.add_argument('--ref', help='Reference FASTA file path for align encoding')
    parser.add_argument('--components', type=int, default=2, help='Number of PCA components')
    parser.add_argument('--clusters', type=int, default=3, help='Number of clusters for K-means')
    parser.add_argument('--output', default='pca_result', help='Output file prefix for PCA plots')
    parser.add_argument('--save_kmer', action='store_true', help='Save kmer frequency data')
    parser.add_argument('--model_output', default='pca_model.pkl', help='Path to save PCA model')
    parser.add_argument('--mapping', default='name_of_the_arrays.txt', help='Path to name_of_the_arrays.txt mapping file')
    args = parser.parse_args()
    
    # 检查参数
    if args.method == 'kmer' and args.k is None:
        parser.error('--k is required when using kmer method')
    if args.method == 'align':
        if args.k is None:
            parser.error('--k is required when using align method')
        if args.ref is None:
            parser.error('--ref is required when using align method')
    
    # 读取BED文件（只保留Array开头的unit）
    units = read_bed_file(args.bed)
    print(f"Read {len(units)} Array units from BED file")
    
    # 读取FASTA文件（使用Bio.SeqIO）
    fasta_dict = read_fasta_file(args.fasta)
    print(f"Read {len(fasta_dict)} sequences from FASTA file")
    
    # 获取unit对应的序列
    unit_names, unit_sequences = get_unit_sequences(units, fasta_dict)
    print(f"Found {len(unit_sequences)} sequences for Array units")
    
    # 编码序列
    if args.method == 'onehot':
        features = one_hot_encode(unit_sequences)
        print(f"One-hot encoding: {features.shape[1]} features")
    elif args.method == 'kmer':
        features, kmers = kmer_encode(unit_sequences, k=args.k)
        print(f"k-mer encoding (k={args.k}): {features.shape[1]} features")
        # 保存kmer频率数据
        if args.save_kmer:
            # 创建DataFrame，第一列是unit的index，后面是kmer频率
            kmer_df = pd.DataFrame(features, columns=kmers)
            kmer_df.insert(0, 'unit_index', unit_names)
            kmer_df.to_csv(f'{args.output}_kmer_frequency.csv', index=False)
            print(f"K-mer frequency data saved to {args.output}_kmer_frequency.csv")
    elif args.method == 'composition':
        features = nucleotide_composition(unit_sequences)
        print(f"Nucleotide composition: {features.shape[1]} features")
    elif args.method == 'align':
        import kmer_map
        features = align_encode(unit_sequences, args.ref, k=args.k)
        print(f"Align encoding (k={args.k}): {features.shape[1]} features")
    
    # 执行PCA
    principal_components, explained_variance, scaler, pca = perform_pca(features, args.components)
    print(f"PCA completed with {args.components} components")
    print(f"Explained variance: {explained_variance}")
    # 保存PCA模型（移动到这里）
    if args.method == 'kmer':
        save_pca_model(scaler, pca, {'k': args.k}, args.model_output)
    else:
        save_pca_model(scaler, pca, {'method': args.method}, args.model_output)
    
    # 绘制结果（按cluster配色）
    plt_cluster, xlim, ylim = plot_pca_by_cluster(principal_components, unit_names, explained_variance, n_clusters=args.clusters)
    plt_cluster.savefig(f'{args.output}_cluster.pdf', dpi=300, bbox_inches='tight')
    print(f"PCA plot by cluster saved to {args.output}_cluster.pdf")
    
    # 绘制结果（按array号配色）
    array_to_nor = read_array_nor_mapping(args.mapping)
    print(f"Read {len(array_to_nor)} Array-NOR mappings")
    plt_array = plot_pca_by_array(principal_components, unit_names, explained_variance, xlim, ylim, array_to_nor)
    plt_array.savefig(f'{args.output}_array.pdf', dpi=300, bbox_inches='tight')
    print(f"PCA plot by array saved to {args.output}_array.pdf")
    
    # 保存PCA结果和聚类信息
    out_csv = f'{args.output}_components.csv'
    pca_df = pd.DataFrame(principal_components, columns=[f'PC{i+1}' for i in range(args.components)])
    pca_df['ID'] = unit_names
    
    # 执行K-means聚类以保存聚类信息
    kmeans = KMeans(n_clusters=args.clusters, random_state=42)
    clusters = kmeans.fit_predict(principal_components)
    pca_df['Cluster'] = clusters
    
    pca_df.to_csv(out_csv, index=False)
    print(f"PCA components and cluster information saved to {out_csv}")

if __name__ == "__main__":
    main()