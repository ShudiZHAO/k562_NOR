import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse

def process_assembly_data(input_file):
    """
    处理unit_structure_map.txt文件，按Assembly分组并计算diff值
    """
    # 读取文件，使用pandas分块处理以避免内存问题
    chunksize = 100000  # 每次处理10万行
    chunk_iter = pd.read_csv(input_file, sep='\t', chunksize=chunksize)
    
    # 用于存储每个Assembly的数据
    assembly_data = {}
    
    # 处理每个数据块
    for chunk in chunk_iter:
        # 计算差值列
        chunk['diff'] = chunk['query_pos'] - chunk['ref_index']
        
        # 按query列（Assembly）分组
        for assembly, group in chunk.groupby('query'):
            if assembly not in assembly_data:
                assembly_data[assembly] = []
            assembly_data[assembly].append(group)
    
    # 合并每个Assembly的数据
    for assembly in assembly_data:
        assembly_data[assembly] = pd.concat(assembly_data[assembly])
    
    return assembly_data

def merge_diff_changes(assembly_data):
    """
    记录每个Assembly中diff值的变化，记录第一个diff值对应的行，
    当发现下一个diff与上一个不同时，记录上一个diff行和下一个diff行，
    并计算该diff值保持的行数
    """
    results = {}
    
    for assembly, data in assembly_data.items():
        # 确保数据按query_pos排序
        data = data.sort_values('query_pos').reset_index(drop=True)
        
        # 存储结果
        assembly_results = []
        
        if len(data) == 0:
            results[assembly] = assembly_results
            continue
        
        # 初始化变量
        current_diff = data['diff'].iloc[0]
        start_row = data.iloc[0]
        start_index = 0
        
        # 遍历数据
        for i in range(1, len(data)):
            row = data.iloc[i]
            if row['diff'] != current_diff:
                # 记录变化
                end_row = data.iloc[i-1]
                row_count = i - start_index
                
                assembly_results.append({
                    'start_row': start_row,
                    'end_row': end_row,
                    'diff_value': current_diff,
                    'row_count': row_count
                })
                
                # 更新当前diff和起始位置
                current_diff = row['diff']
                start_row = row
                start_index = i
        
        # 处理最后一组
        if len(data) > start_index:
            end_row = data.iloc[-1]
            row_count = len(data) - start_index
            
            assembly_results.append({
                'start_row': start_row,
                'end_row': end_row,
                'diff_value': current_diff,
                'row_count': row_count
            })
        
        results[assembly] = assembly_results
    
    return results

def filter_diff_changes(results):
    """
    过滤满足条件的行：若某行Diff_Value的上一行与下一行差值小于1000，
    该行与上一行和下一行的差值均大于3000，则去掉该行
    """
    filtered_results = {}
    
    for assembly, assembly_results in results.items():
        filtered = []
        
        for i in range(len(assembly_results)):
            # 跳过首行和末行（没有上一行或下一行）
            if i == 0 or i == len(assembly_results) - 1:
                filtered.append(assembly_results[i])
                continue
            
            current_diff = assembly_results[i]['diff_value']
            prev_diff = assembly_results[i-1]['diff_value']
            next_diff = assembly_results[i+1]['diff_value']
            
            # 检查条件
            prev_next_diff = abs(prev_diff - next_diff)
            curr_prev_diff = abs(current_diff - prev_diff)
            curr_next_diff = abs(current_diff - next_diff)
            
            if prev_next_diff < 1000 and curr_prev_diff > 3000 and curr_next_diff > 3000:
                # 满足条件，跳过该行（不添加到filtered）
                # print(f"Jump {i}: {assembly_results[i]}")
                continue
            
            filtered.append(assembly_results[i])
        
        filtered_results[assembly] = filtered
    
    return filtered_results

def save_merge_results(results, output_file):
    """
    保存合并结果到文件
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        # 写入表头
        f.write('Assembly\tStart_Query_Pos\tStart_Ref_Index\tStart_Diff\tEnd_Query_Pos\tEnd_Ref_Index\tEnd_Diff\tDiff_Value\tRow_Count\n')
        
        # 写入数据
        for assembly, assembly_results in results.items():
            for result in assembly_results:
                start_row = result['start_row']
                end_row = result['end_row']
                
                f.write(f"{assembly}\t{start_row['query_pos']}\t{start_row['ref_index']}\t{start_row['diff']}\t{end_row['query_pos']}\t{end_row['ref_index']}\t{end_row['diff']}\t{result['diff_value']}\t{result['row_count']}\n")

def plot_from_diff_changes(diff_changes_file, output_dir, file_format='pdf'):
    """
    从diff_changes.txt文件生成图像
    """
    # 读取diff_changes.txt文件
    df = pd.read_csv(diff_changes_file, sep='\t')
    
    # 按Assembly分组
    assembly_groups = df.groupby('Assembly')
    
    # 创建输出文件夹
    os.makedirs(output_dir, exist_ok=True)
    
    for assembly, group in assembly_groups:
        # 创建图形
        sorted_group = group.sort_values('Start_Query_Pos').reset_index(drop=True)

        # 计算x的总范围
        x_min = sorted_group['Start_Query_Pos'].min()
        x_max = sorted_group['End_Query_Pos'].max()
        x_range = x_max - x_min

        plt.figure(figsize=(x_range/5000, 8))

        # 设置x轴刻度
        plt.xlim(x_min-1000, x_max+1000)
        plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(10000))
        plt.gca().xaxis.set_minor_locator(ticker.MultipleLocator(2000))
        plt.ticklabel_format(style='plain', axis='both')
        
        # 用红色线连接跳跃的点
        if len(sorted_group) > 1:
            for i in range(1, len(sorted_group)):
                prev_row = sorted_group.iloc[i-1]
                curr_row = sorted_group.iloc[i]
                
                # 连接前一段的结束点和当前段的起始点
                plt.plot([prev_row['End_Query_Pos'], curr_row['Start_Query_Pos']], 
                         [prev_row['End_Ref_Index'], curr_row['Start_Ref_Index']], 
                         color='red', linestyle='-', linewidth=0.7)

        # 绘制线段
        for _, row in group.iterrows():
            plt.plot([row['Start_Query_Pos'], row['End_Query_Pos']], 
                     [row['Start_Ref_Index'], row['End_Ref_Index']], color='black',
                     marker='o', markersize=1, linestyle='-', linewidth=1.5)

        # 添加标题和标签
        plt.title(f'{assembly} Structure')
        plt.xlabel('Query Position')
        plt.ylabel('Ref Index')
        
        # 保存图形
        output_file = os.path.join(output_dir, f'{assembly}.{file_format}')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        
        # 关闭图形以释放内存
        plt.close()
        
        print(f"Saved plot for {assembly} to {output_file}")

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='Process unit structure map and generate plots')
    parser.add_argument('--input', type=str, default='unit_structure_map.txt', help='Input unit structure map file')
    parser.add_argument('--diff-changes', type=str, help='Directly use diff_changes.txt to generate plots (skips processing)')
    parser.add_argument('--output-dir', type=str, default='simplify_continuous_structure', help='Output directory for plots')
    parser.add_argument('--format', type=str, default='pdf', choices=['pdf', 'png'], help='Output file format')
    args = parser.parse_args()

    if args.diff_changes:
        # 直接从diff_changes.txt生成图像
        print(f"Generating plots from {args.diff_changes}...")
        plot_from_diff_changes(args.diff_changes, args.output_dir, args.format)
        print(f"All plots saved to {args.output_dir} directory")
    else:
        # 处理unit_structure_map.txt文件
        print(f"Processing {args.input}...")
        assembly_data = process_assembly_data(args.input)
        print(f"Processed {len(assembly_data)} assemblies")
        
        # 合并diff变化
        print("Merging diff changes...")
        merge_results = merge_diff_changes(assembly_data)

        # 过滤满足条件的行
        print("Filtering diff changes...")
        merge_results = filter_diff_changes(merge_results)
        
        # 保存结果
        output_file = 'diff_changes.txt'
        save_merge_results(merge_results, output_file)
        print(f"Merge results saved to {output_file}")
        
        # 绘制每个assembly的结构
        print("Plotting assembly structures...")
        plot_from_diff_changes(output_file, args.output_dir, args.format)
        print(f"All plots saved to {args.output_dir} directory")

if __name__ == "__main__":
    main()