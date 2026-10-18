import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# 定义type颜色（根据数据中有5种type，定义5种颜色）
TYPE_COLORS = {
    1: '#5bb8ac',  # Type 1 - 粉色
    2: '#e5e5e5',  # Type 2 - 蓝色
    3: '#f0c238',  # Type 3 - 绿色
    4: '#762A83',  # Type 4 - 黄色
    5: '#f07979',  # Type 5 - 紫色
}

def load_unit_type_data(unit_type_file):
    """加载unit的type和subtype数据"""
    df = pd.read_csv(unit_type_file, sep='\t')
    return df

def load_array_to_nor_mapping(array_nor_file):
    """加载Array到NOR的映射"""
    mapping = {}
    with open(array_nor_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                array_name = parts[1]  # Array_X
                nor_name = parts[2]    # NOR X
                if nor_name == 'Dropped':
                    continue
                # 提取NOR编号
                nor_num = int(nor_name.split()[-1])
                mapping[array_name] = nor_num
    return mapping

def load_unit_positions(bed_file):
    """加载unit位置信息"""
    df = pd.read_csv(bed_file, sep='\t', header=None, 
                     names=['chrom', 'start', 'end', 'name', 'score'])
    # 过滤：保留开头为"Array"但排除开头为"Array_1_"的行
    df = df[df['name'].str.startswith('Array') & ~df['name'].str.startswith('Array_1_')]
    return df

def extract_methylation_for_region(bedgraph_file, chrom, start, end):
    """从bedgraph文件中提取指定区域的甲基化数据"""
    values = []
    with open(bedgraph_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            bg_chrom = parts[0]
            bg_start = int(parts[1])
            bg_end = int(parts[2])
            bg_value = float(parts[3])
            
            if bg_chrom == chrom and not (bg_end <= start or bg_start >= end):
                overlap_start = max(start, bg_start)
                overlap_end = min(end, bg_end)
                values.extend([bg_value] * (overlap_end - overlap_start))
    
    return values

def calculate_unit_methylation_profile(unit_info, bedgraph_file, bin_size):
    """计算单个unit的甲基化分布"""
    chrom = unit_info['chrom']
    start = unit_info['start']
    end = unit_info['end']
    unit_length = end - start
    
    meth_values = extract_methylation_for_region(bedgraph_file, chrom, start, end)
    
    # 如果没有甲基化数据，返回空列表和长度
    if len(meth_values) == 0:
        return [], unit_length
    
    # 按bin_size进行bin化
    num_bins = int(np.ceil(unit_length / bin_size))
    binned_values = []
    
    for i in range(num_bins):
        start_idx = i * bin_size
        end_idx = min((i + 1) * bin_size, len(meth_values))
        if start_idx < len(meth_values):
            bin_values = meth_values[start_idx:end_idx]
            binned_values.append(np.mean(bin_values) if bin_values else np.nan)
        else:
            binned_values.append(np.nan)
    
    return binned_values, unit_length

def parse_unit_id(unit_id):
    """解析unit ID中的x和y值"""
    match = re.match(r'Array_(\d+)_(\d+)', unit_id)
    if match:
        return int(match.group(1)), int(match.group(2))
    return 0, 0

def parse_new_type(new_type_str):
    """从new_type字符串中解析type和subtype"""
    if pd.isna(new_type_str) or new_type_str is None:
        return 0, 0
    parts = str(new_type_str).strip().split('-')
    if len(parts) >= 2:
        try:
            return int(parts[0]), int(parts[1])
        except ValueError:
            return 0, 0
    return 0, 0

def plot_individual_unit_charts(unit_type_df, array_to_nor, unit_positions, bedgraph_file, bin_size, output_file):
    # 添加NOR列到unit_type_df
    unit_type_df['NOR'] = unit_type_df['origin_name'].apply(
        lambda x: array_to_nor.get(f'Array_{parse_unit_id(x)[0]}', 0)
    )
    
    # 从new_type列解析type和subtype
    unit_type_df[['type_from_new', 'subtype_from_new']] = unit_type_df['new_type'].apply(
        lambda x: pd.Series(parse_new_type(x))
    )
    
    # 移除NOR为0或type为0的行
    unit_type_df = unit_type_df[(unit_type_df['NOR'] != 0) & (unit_type_df['type_from_new'] != 0)]
    
    # 按NOR排序（1-18），同NOR内按index排序
    unit_type_df = unit_type_df.sort_values(by=['NOR', 'index'])
    
    # 获取排序后的units
    units = unit_type_df['origin_name'].tolist()
    num_units = len(units)
    
    if num_units == 0:
        print("No units to plot")
        return
    
    # 计算布局：每行18个unit
    num_cols = 18
    
    # 计算每个NOR有多少个units
    nor_counts = unit_type_df['NOR'].value_counts().sort_index()
    
    # 计算总行数（每个NOR可能占用多行，如果unit数量超过18）
    num_rows = sum((count + num_cols - 1) // num_cols for _, count in nor_counts.items())
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(2 * num_cols, 2 * num_rows))
    
    # 处理单行情况
    if num_rows == 1:
        axes = axes.reshape(1, -1)
    
    # 初始化所有子图为完全空白
    for r in range(num_rows):
        for c in range(num_cols):
            axes[r, c].axis('off')
            axes[r, c].set_frame_on(False)
    
    # 为每个type的subtype计算透明度映射（从new_type解析）
    transparency_map = {}
    for type_val in unit_type_df['type_from_new'].unique():
        subtypes = sorted(unit_type_df[unit_type_df['type_from_new'] == type_val]['subtype_from_new'].unique())
        max_subtype = len(subtypes) if len(subtypes) > 0 else 1
        for idx, subtype in enumerate(subtypes):
            # 透明度从0.3到0.9，根据subtype数量分配
            transparency_map[(type_val, subtype)] = 0.3 + (idx / max_subtype) * 0.6
    
    # 按NOR分组填充子图
    current_row = 0
    current_col = 0
    
    for nor_num in sorted(nor_counts.index):
        # 获取当前NOR的所有units
        nor_units = unit_type_df[unit_type_df['NOR'] == nor_num]['origin_name'].tolist()
        
        for unit_id in nor_units:
            # 如果当前行已满，换到下一行
            if current_col >= num_cols:
                current_row += 1
                current_col = 0
            
            ax = axes[current_row, current_col]
            
            # 获取unit信息
            unit_info = unit_positions[unit_positions['name'] == unit_id]
            if len(unit_info) == 0:
                print(f"Warning: Unit {unit_id} not found in BED file")
                current_col += 1
                continue
            
            # 获取unit的type和subtype（从new_type解析）
            type_info = unit_type_df[unit_type_df['origin_name'] == unit_id]
            if len(type_info) == 0:
                print(f"Warning: Unit {unit_id} not found in unit_to_type.xls")
                current_col += 1
                continue
            
            unit_type = type_info['type_from_new'].iloc[0]
            unit_subtype = type_info['subtype_from_new'].iloc[0]
            unit_index = type_info['index'].iloc[0]
            
            # 获取颜色和透明度
            color = TYPE_COLORS.get(unit_type, '#808080')  # 默认灰色
            alpha = transparency_map.get((unit_type, unit_subtype), 0.8)
            
            # 计算甲基化分布
            binned_values, unit_length = calculate_unit_methylation_profile(
                unit_info.iloc[0], bedgraph_file, bin_size
            )
            
            if len(binned_values) == 0:
                current_col += 1
                continue
            
            # 重新启用坐标轴和边框
            ax.axis('on')
            ax.set_frame_on(True)
            
            # 绘制柱状图
            positions = np.arange(len(binned_values)) * bin_size + bin_size / 2
            valid_mask = ~np.isnan(binned_values)
            valid_positions = positions[valid_mask]
            valid_values = np.array(binned_values)[valid_mask]
            
            ax.bar(valid_positions, valid_values, width=bin_size*0.8, 
                   color=color, alpha=alpha)
            
            # 设置标题：NOR x y (Type a-b)
            ax.set_title(f'NOR {nor_num} {unit_index}\n(Type {unit_type}-{unit_subtype})', fontsize=8)
            
            # 设置x轴范围（每个unit的上限是unit长度）
            ax.set_xlim(0, unit_length)
            
            # 设置x轴刻度：每5000一个tick，每20000一个label
            max_tick = ((unit_length // 5000) + 1) * 5000
            xticks = np.arange(0, max_tick + 1, 5000)
            ax.set_xticks(xticks)
            ax.set_xticklabels([str(int(tick)) if tick % 20000 == 0 else '' for tick in xticks], 
                              fontsize=5)
            
            ax.set_xlabel(f'Position (bp)', fontsize=6)
            ax.set_ylabel('Methylation', fontsize=6)
            ax.tick_params(axis='y', labelsize=5)
            ax.grid(True, alpha=0.3)
            
            current_col += 1
        
        # NOR结束后换到下一行（即使当前行未满）
        current_row += 1
        current_col = 0
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Generated {output_file} with {num_units} unit charts")

def main():
    parser = argparse.ArgumentParser(description='Plot individual unit methylation distribution by type')
    parser.add_argument('--unit_type_file', required=True, help='Path to unit_to_type.xls file')
    parser.add_argument('--array_nor_file', required=True, help='Path to name_of_the_arrays.txt file')
    parser.add_argument('--bed_file', required=True, help='Path to named_units.bed file')
    parser.add_argument('--bedgraph_file', required=True, help='Path to stats_Me.bedgraph file')
    parser.add_argument('--bin_size', type=int, default=100, help='Bin size in bp (default: 100)')
    parser.add_argument('--output', default='individual_unit_methylation', help='Output file prefix')
    
    args = parser.parse_args()
    
    print("Loading unit type data...")
    unit_type_df = load_unit_type_data(args.unit_type_file)
    print(f"Loaded {len(unit_type_df)} units")
    
    print("Loading Array to NOR mapping...")
    array_to_nor = load_array_to_nor_mapping(args.array_nor_file)
    print(f"Loaded {len(array_to_nor)} mappings")
    
    print("Loading unit positions...")
    unit_positions = load_unit_positions(args.bed_file)
    print(f"Loaded {len(unit_positions)} unit positions")
    
    # 绘制所有units的单独图表（按NOR分组）
    print("Plotting individual unit charts...")
    plot_individual_unit_charts(
        unit_type_df,
        array_to_nor,
        unit_positions, 
        args.bedgraph_file, 
        args.bin_size,
        f'{args.output}.pdf'
    )
    
    print("Done!")

if __name__ == "__main__":
    main()