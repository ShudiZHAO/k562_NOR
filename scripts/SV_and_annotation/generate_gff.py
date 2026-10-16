#!/usr/bin/env python3

from Bio import SeqIO
import sys
import os
import re
from collections import defaultdict

# Import functions from kmer_map.py
from kmer_map import build_ref_index, query_sequence

def get_base_feature_name(feature_id):
    """去掉feature id末尾的参考编号，如 'R Repeat 1' -> 'R Repeat'，'ITS1' -> 'ITS1'"""
    return re.sub(r'\s+\d+$', '', feature_id)

def parse_features_file(features_file):
    """Parse Features.txt file and return features"""
    features = []
    with open(features_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 4:
                # 第一列：ID，第二列：type，第三列：start，第四列：end
                feature_id = parts[0]
                type_ = parts[1]
                start = int(parts[2])
                end = int(parts[3])
                features.append({
                    'id': feature_id,
                    # 去掉参考feature编号后的基础名，用于unit内编号，如 'R Repeat 1' -> 'R Repeat'
                    'base': get_base_feature_name(feature_id),
                    'type': type_,
                    'start': start,
                    'end': end,
                    # 每行一个唯一编号：仅用于区分参考行；合并按基础名（老版本方法）进行，
                    # Butterfly/Long repeat则按末端motif+拷贝长度精确拆分
                    'index': len(features)
                })
    return features

def parse_named_units(units_bed):
    """
    解析named_units.bed，确定每个assembly的unit划分
    返回 {assembly: [unit_dict, ...]}，unit_dict包含name/start/end/y
    bed共5列：assembly, start, end, name, length
    坐标约定与simplify_variation_structure_by_cluster.py的load_named_units一致：
    start为1-based inclusive，end为exclusive，即1-based区间为[start, end-1]
    y的取值：Array_x_y -> y；Right_fragment -> R；Left_fragment -> L；其他 -> 名字首字母大写（如Dropped -> D）
    """
    named_units = {}
    with open(units_bed, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            assembly = parts[0]
            start = int(parts[1])
            end = int(parts[2]) - 1  # end为exclusive，转换为inclusive
            name = parts[3]

            m = re.match(r'^Array_\d+_(\d+)$', name)
            if m:
                y = m.group(1)
            elif name == 'Right_fragment':
                y = 'R'
            elif name == 'Left_fragment':
                y = 'L'
            else:
                y = name[0].upper()
                print(f"Warning: unknown unit name '{name}' in {assembly}, using y='{y}'")

            named_units.setdefault(assembly, []).append({
                'name': name,
                'start': start,
                'end': end,
                'y': y
            })
    return named_units

def get_unit_for_position(pos, units):
    """根据1-based位置找到对应的unit；若不在任何unit内（如边界1bp间隙），归入距离最近的unit"""
    for unit in units:
        if unit['start'] <= pos <= unit['end']:
            return unit
    nearest = min(units, key=lambda u: min(abs(pos - u['start']), abs(pos - u['end'])))
    print(f"Warning: position {pos} not in any unit, assigned to nearest unit {nearest['name']}")
    return nearest

def build_ref_pos_to_feature_map(ref_features):
    """Build a map from reference position to feature"""
    # Sort features by start position
    sorted_features = sorted(ref_features, key=lambda x: x['start'])
    
    def get_feature_for_pos(pos):
        """Get the feature that contains the given position"""
        for feature in sorted_features:
            if feature['start'] <= pos <= feature['end']:
                return feature
        return None
    
    return get_feature_for_pos

def merge_features(features, interval_distance = 50):
    """Merge adjacent or overlapping features"""
    if not features:
        return []
    
    # 按起始位置排序
    sorted_features = sorted(features, key=lambda x: x['start'])
    
    merged_intervals = []
    current_interval = sorted_features[0].copy()
    
    for interval in sorted_features[1:]:
        # 检查是否是同一个feature（按基础名）且相邻或重叠
        # 注意：同一基础名的多行参考feature（如3个R Repeat行）按老版本方法合并为一组；
        # Butterfly/Long repeat不在此处理（走motif精确拆分路径）
        if (interval['base'] == current_interval['base'] and
            interval['start'] <= current_interval['end'] + interval_distance):
            # 合并区间
            current_interval['end'] = max(current_interval['end'], interval['end'])
        else:
            # 添加当前区间到列表
            merged_intervals.append(current_interval)
            # 开始新的区间
            current_interval = interval.copy()
    
    # 添加最后一个区间
    merged_intervals.append(current_interval)
    
    return merged_intervals

def is_in_unit_region(qpos, unit_regions):
    """
    检查qpos是否位于任一unit区间内（unit_regions为[(start, end), ...]，1-based inclusive）
    有了named_units.bed之后，unit区间以外的匹配点（原flank区域）直接被过滤，不再需要flank.bed
    """
    for start, end in unit_regions:
        if start <= qpos <= end:
            return True
    return False

# ---------------- Butterfly/Long repeat 的精确分隔 ----------------
# Butterfly/Long repeat拷贝通常以该30bp motif结尾（测序错误/突变允许少量mismatch），
# 且拷贝长度一般在4000-5000bp。据此把每个unit内butterfly的匹配点切成单个拷贝。
BUTTERFLY_END_MOTIF = 'CCCGATAAGCTCCTCGACTTCAACATAAAC'
BUTTERFLY_MAX_MISMATCH = 3        # motif允许的最大mismatch数
# <=3个mismatch时motif的4段种子中至少一段完全匹配（鸽笼原理），先用种子粗筛再验证，速度快
BUTTERFLY_MOTIF_SEEDS = [(0, 8), (8, 15), (15, 22), (22, 30)]
BUTTERFLY_GAP_MERGE = 3500        # 同一拷贝内部未匹配间隙的聚簇阈值（实测最大~1.8kb）
BUTTERFLY_MIN_COPY_LEN = 3000     # 接受一个motif末端的最短拷贝长度（拷贝一般4000-5000bp）
BUTTERFLY_MIN_FEATURE_LEN = 80    # 短于该长度的拷贝/尾巴视为噪声丢弃
BUTTERFLY_TAIL_MARGIN = 100       # motif可超出最后一个匹配点的大致余量
# 真实butterfly簇的k-mer匹配密度约0.66-0.86；共享k-mer造成的零星假簇密度仅~0.01，
# 另有span仅1-2bp的碎片簇，按下面阈值过滤
BUTTERFLY_MIN_CLUSTER_DENSITY = 0.25
BUTTERFLY_MIN_CLUSTER_SPAN = 300

def is_butterfly_feature(feature):
    """判断是否为Butterfly/Long repeat类feature"""
    return 'Butterfly' in feature['base']

def find_motif_ends(seq, lo, hi, motif=BUTTERFLY_END_MOTIF,
                    max_mismatch=BUTTERFLY_MAX_MISMATCH):
    """
    在seq的[lo, hi]（1-based，含端点，hi可超出序列）范围内查找motif，
    返回[(end_pos, mismatch), ...]，end_pos为motif最后一个碱基的1-based位置（即拷贝末端），升序。
    """
    L = len(motif)
    # 窗口起点w(0-based)：motif占1-based [w+1, w+L]，要求 lo <= w+L <= hi
    w_lo = max(0, lo - L)
    w_hi = min(len(seq) - L, hi - L)
    candidates = set()
    for a, b in BUTTERFLY_MOTIF_SEEDS:
        seed = motif[a:b]
        i = max(0, w_lo + a)
        limit = w_hi + a
        while i <= limit:
            i = seq.find(seed, i, limit + b)
            if i == -1:
                break
            candidates.add(i - a)
            i += 1
    ends = []
    for w in sorted(candidates):
        d = sum(1 for x, y in zip(seq[w:w + L], motif) if x != y)
        if d <= max_mismatch:
            ends.append((w + L, d))
    return ends

def split_butterfly_points(points, seq):
    """
    把一个unit内butterfly的匹配点位置（已排序的qpos列表）按末端motif精确分隔成拷贝。
    返回[(start, end), ...]（1-based inclusive，升序）。
    步骤：先按间隙<=BUTTERFLY_GAP_MERGE聚成连续簇；每簇内查找motif末端，
    用「与前一个边界的距离>=BUTTERFLY_MIN_COPY_LEN」贪心确认拷贝末端，
    最后一motif之后的剩余尾巴>=BUTTERFLY_MIN_FEATURE_LEN时作为(可能不完整的)最后一个拷贝保留。
    """
    if not points:
        return []
    # 聚簇：相邻点间隙<=BUTTERFLY_GAP_MERGE视为同一（或连续拷贝的）区域
    clusters = []
    current = [points[0]]
    for q in points[1:]:
        if q - current[-1] <= BUTTERFLY_GAP_MERGE:
            current.append(q)
        else:
            clusters.append(current)
            current = [q]
    clusters.append(current)

    copies = []
    for cluster in clusters:
        c_start, c_end = cluster[0], cluster[-1]
        span = c_end - c_start + 1
        density = len(cluster) / span
        if span < BUTTERFLY_MIN_CLUSTER_SPAN or density < BUTTERFLY_MIN_CLUSTER_DENSITY:
            print(f"Warning: dropped sparse/short butterfly point cluster "
                  f"[{c_start}, {c_end}] (span={span}, points={len(cluster)}, density={density:.3f})")
            continue
        motif_ends = find_motif_ends(seq, c_start, c_end + BUTTERFLY_TAIL_MARGIN)
        cluster_copies = []
        prev_end = c_start - 1
        for m_end, _ in motif_ends:
            if m_end - prev_end >= BUTTERFLY_MIN_COPY_LEN:
                cluster_copies.append((prev_end + 1, m_end))
                prev_end = m_end
            else:
                print(f"Warning: butterfly motif end {m_end} too close to previous "
                      f"boundary {prev_end} (copy would be <{BUTTERFLY_MIN_COPY_LEN}bp), skipped")
        # 最后一个motif之后的尾巴：可能是被截断/未匹配完的最后一个拷贝
        tail_start = prev_end + 1
        if c_end - tail_start + 1 >= BUTTERFLY_MIN_FEATURE_LEN:
            cluster_copies.append((tail_start, c_end))
        if not motif_ends:
            print(f"Warning: no butterfly end motif found in region "
                  f"[{c_start}, {c_end}], kept as a single feature")
        copies.extend(cluster_copies)
    return copies

def map_and_merge_features(query_matches, get_feature_for_pos, unit_regions, query_seq=None):
    """
    Map reference features to query sequence and merge adjacent features.
    Butterfly/Long repeat与R Repeat的点不按基础名合并：Butterfly按末端motif+拷贝长度、
    R Repeat按首尾motif+拷贝长度精确分隔成单个拷贝；
    其余feature按基础名（老版本方法）合并。
    """
    if not query_matches:
        return []

    # Sort matches by query position
    sorted_matches = sorted(query_matches, key=lambda x: x[0])

    # Map each query position to its corresponding feature
    # butterfly点/R repeat点与其他feature点分开处理
    butterfly_points = []          # butterfly匹配点的qpos（按参考行无关，统一处理）
    butterfly_template = None      # 用于拷贝feature的id/base/type模板
    rrepeat_points = []            # R repeat匹配点的qpos（按参考行无关，统一处理）
    rrepeat_template = None        # 用于拷贝feature的id/base/type模板
    query_pos_feature_pairs = []   # 其他feature的(qpos, feature)
    for qpos, ref_pos in sorted_matches:
        # 只保留落在unit区间内的匹配点
        if is_in_unit_region(qpos, unit_regions):
            feature = get_feature_for_pos(ref_pos)
            if feature:
                if query_seq is not None and is_butterfly_feature(feature):
                    if butterfly_template is None:
                        butterfly_template = feature
                    butterfly_points.append(qpos)
                elif query_seq is not None and is_rrepeat_feature(feature):
                    if rrepeat_template is None:
                        rrepeat_template = feature
                    rrepeat_points.append(qpos)
                else:
                    query_pos_feature_pairs.append((qpos, feature))

    # ---------- 其他feature：按参考行合并 ----------
    filtered_features = []
    if query_pos_feature_pairs:
        # 从匹配的点创建初始区间，每个点视为一个长度为1的区间
        initial_intervals = []
        for qpos, feature in query_pos_feature_pairs:
            initial_intervals.append({
                'id': feature['id'],
                'base': feature['base'],
                'type': feature['type'],
                'index': feature['index'],
                'start': qpos,
                'end': qpos
            })

        # 合并相邻或重叠的同一参考feature行区间
        merged_intervals = merge_features(initial_intervals)

        # 过滤掉长度小于80bp的feature
        for feature in merged_intervals:
            length = feature['end'] - feature['start'] + 1
            if feature['type'] == 'Target':
                filtered_features.append(feature)
            elif length >= 80:
                filtered_features.append(feature)

        # 重新merge一遍，合并距离为100bp
        filtered_features = merge_features(filtered_features, 100)

    # ---------- Butterfly/Long repeat：按末端motif精确分隔拷贝 ----------
    if butterfly_points:
        template = butterfly_template
        for c_start, c_end in split_butterfly_points(sorted(butterfly_points), query_seq):
            filtered_features.append({
                'id': template['id'],
                'base': template['base'],
                'type': template['type'],
                'index': template['index'],
                'start': c_start,
                'end': c_end
            })

    # ---------- R Repeat：按首尾motif+拷贝长度精确分隔拷贝 ----------
    # 与butterfly一样在最终re-merge之后追加，避免相邻拷贝（间隔仅1bp）被重新合并
    if rrepeat_points:
        template = rrepeat_template
        for c_start, c_end in split_rrepeat_points(sorted(rrepeat_points), query_seq):
            filtered_features.append({
                'id': template['id'],
                'base': template['base'],
                'type': template['type'],
                'index': template['index'],
                'start': c_start,
                'end': c_end
            })

    filtered_features.sort(key=lambda x: x['start'])
    return filtered_features

def is_rrepeat_feature(feature):
    """判断是否为R Repeat类feature"""
    return feature['base'] == 'R Repeat'

# ---------------- R Repeat 的精确分隔 ----------------
# R Repeat拷贝以TGAGACT(N)...起始、以CCTTGGTGCCTTCTCGGCTC终止（20bp，其后常跟
# 1个T+下一拷贝起点），相邻拷贝紧密相连；完整拷贝长度约550-800bp。
# 据此把每个unit内R的匹配点切成单个拷贝；feature严格限制在实际比对上的区域内。
RREPEAT_START_MOTIF = 'TGAGACTTAGCCGCTGTCTCGCCGTG'  # 26bp（参考3个拷贝变体两两差异<=3bp）
RREPEAT_START_MAX_MM = 5        # 实测真实起点motif mm 0-4（部分unit有插入变异），拷贝内部
                                # 的退化相似位点mm 7-9，阈值5恰好分开
RREPEAT_END_MOTIF = 'CCTTGGTGCCTTCTCGGCTC'  # 20bp终止序列（参考各拷贝末端完全一致）
RREPEAT_END_MAX_MM = 2          # 实测query终点motif mm 0-1
RREPEAT_END_ADJ_TOL = 15        # 同一连接点的首尾证据相距<=该值视为同一边界
RREPEAT_GAP_MERGE = 1000        # 同一R位点内部未匹配间隔的聚簇阈值（实测~100bp）
RREPEAT_MIN_CLUSTER_DENSITY = 0.25  # 簇密度过滤阈值：真实R位点匹配点极密(~0.8+)，
                                    # 退化的R同源交叉匹配点稀疏(实测~0.02-0.12)
RREPEAT_MIN_COPY_LEN = 300      # 相邻连接点最短间距（完整拷贝最短长度）
RREPEAT_MAX_COPY_LEN = 900      # 完整拷贝最长长度
RREPEAT_MIN_FEATURE_LEN = 80    # 短于该长度的部分拷贝/尾巴视为噪声丢弃
RREPEAT_TAIL_MARGIN = 100       # 终点motif搜索可超出最后一个匹配点的余量
RREPEAT_START_MARGIN = 30       # 起点motif搜索可早于第一个匹配点的余量
RREPEAT_MAX_INTERNAL_GAP = 250  # 拷贝内部匹配点最大间隙（实测真实拷贝<=200），超过则按
                                # 马赛克式拼贴切开分别成段
RREPEAT_MIN_COPY_DENSITY = 0.6  # 单拷贝匹配点密度下限（实测真实拷贝>=0.67），
                                # 低于该值的碎片视为马赛克式噪声丢弃

def _find_motif_hits(seq, lo, hi, motif, max_mm):
    """在seq的[lo,hi]（1-based含端点，可越界）内查找motif所有出现（mismatch<=max_mm），
    返回[(motif起始位置1-based, mismatch), ...]升序。R区域短，直接扫描即可。"""
    L = len(motif)
    hits = []
    w_lo = max(0, lo - 1)
    w_hi = min(len(seq) - L, hi - 1)
    for w in range(w_lo, w_hi + 1):
        d = 0
        for x, y in zip(seq[w:w + L], motif):
            if x != y:
                d += 1
                if d > max_mm:
                    break
        if d <= max_mm:
            hits.append((w + 1, d))
    return hits

def split_rrepeat_points(points, seq):
    """
    把一个unit内R Repeat的匹配点位置（已排序的qpos列表）按首尾motif+拷贝长度
    精确分隔成单个拷贝。返回[(start, end), ...]（1-based inclusive，升序）。
    结构：每个拷贝以TGAGACT(N)...起始、以CCTTGGTGCCTTCTCGGCTC终止（20bp，
    其后常跟1个T+下一拷贝起点），完整拷贝长度约550-800bp。
    步骤：
    1. 按间隔<=RREPEAT_GAP_MERGE聚成R位点簇，密度<RREPEAT_MIN_CLUSTER_DENSITY的
       稀疏簇（IGS中退化的R同源区，马赛克式散点比对）直接丢弃；
    2. 簇内用「起点motif位置」和「终点motif末位+2(=下一拷贝起点)」作为拷贝间
       连接点证据（首尾任一保守即可锚定，相距<=RREPEAT_END_ADJ_TOL的证据合并为
       同一连接点，优先取起点motif位置）；相邻连接点距离<最短拷贝长的视为假阳性丢弃；
    3. 相邻连接点之间的匹配点为一个拷贝。feature严格限制在实际比对上的区域内：
       起点/终点取区间内第一个/最后一个匹配点，motif只用于定界，不把比对不上
       R repeat的序列算进feature（GFF输出终点+19恰好对应末位k-mer窗口末端）；
    4. 拷贝内部匹配点间隙>RREPEAT_MAX_INTERNAL_GAP的（马赛克式拼贴）切开分别成段；
       每段长度>=RREPEAT_MIN_FEATURE_LEN且匹配点密度>=RREPEAT_MIN_COPY_DENSITY才保留。
    """
    if not points:
        return []
    # 聚簇：相邻点间隔<=RREPEAT_GAP_MERGE视为同一R位点
    clusters = []
    current = [points[0]]
    for q in points[1:]:
        if q - current[-1] <= RREPEAT_GAP_MERGE:
            current.append(q)
        else:
            clusters.append(current)
            current = [q]
    clusters.append(current)

    end_L = len(RREPEAT_END_MOTIF)
    copies = []
    for cluster in clusters:
        c_start, c_end = cluster[0], cluster[-1]
        # 簇密度过滤：真实R位点匹配点极密，稀疏簇为退化R同源区的交叉匹配噪声
        span = c_end - c_start + 1
        density = len(cluster) / span
        if density < RREPEAT_MIN_CLUSTER_DENSITY:
            print(f"Warning: dropped sparse/short R repeat point cluster "
                  f"[{c_start}, {c_end}] (span={span}, points={len(cluster)}, density={density:.3f})")
            continue
        start_hits = _find_motif_hits(seq, c_start - RREPEAT_START_MARGIN, c_end + 60,
                                      RREPEAT_START_MOTIF, RREPEAT_START_MAX_MM)
        end_hits = _find_motif_hits(seq, c_start - end_L + 1, c_end + RREPEAT_TAIL_MARGIN,
                                    RREPEAT_END_MOTIF, RREPEAT_END_MAX_MM)
        start_pos = set(p for p, _ in start_hits)
        # 连接点证据：起点motif位置；或终点motif末位+2（参考结构为终点motif+1个T
        # +下一拷贝起点，即终点motif末位+2=下一拷贝起点）
        cands = sorted(start_pos | set(p + end_L + 1 for p, _ in end_hits))
        # 相距<=RREPEAT_END_ADJ_TOL的证据合并为同一连接点，优先用起点motif位置
        junctions = []
        i = 0
        while i < len(cands):
            j = i
            grp = [cands[i]]
            while j + 1 < len(cands) and cands[j + 1] - grp[-1] <= RREPEAT_END_ADJ_TOL:
                j += 1
                grp.append(cands[j])
            in_grp = [p for p in grp if p in start_pos]
            junctions.append(in_grp[0] if in_grp else grp[0])
            i = j + 1

        if not junctions:
            if span >= RREPEAT_MIN_FEATURE_LEN:
                print(f"Warning: no R repeat start/end motif found in region "
                      f"[{c_start}, {c_end}], kept as a single feature")
                copies.append((c_start, c_end))
            else:
                print(f"Warning: dropped short R repeat cluster [{c_start}, {c_end}] "
                      f"(len={span}, no R repeat motif)")
            continue

        # 距离过滤：相邻连接点<最短拷贝长的视为假阳性丢弃
        kept = [junctions[0]]
        for j in junctions[1:]:
            if j - kept[-1] >= RREPEAT_MIN_COPY_LEN:
                kept.append(j)
        junctions = kept

        # 连接点切出的拷贝区间：簇开头到第一个连接点为前部不完整拷贝，
        # 相邻连接点之间为一个拷贝，最后一个连接点之后为末拷贝（可能不完整）
        seg_bounds = []
        if junctions[0] > c_start:
            seg_bounds.append((c_start, junctions[0] - 1))
        prev = junctions[0]
        for j in junctions[1:]:
            seg_bounds.append((prev, min(j - 1, c_end)))
            prev = j
            if j > c_end:
                break
        if prev <= c_end:
            seg_bounds.append((prev, c_end))

        # 每个区间只取实际比对上的匹配点：起点/终点取区间内第一/最后一个匹配点，
        # 内部大间隙（马赛克式拼贴）切开，逐段校验长度与密度后保留
        for lo, hi in seg_bounds:
            seg_pts = [p for p in cluster if lo <= p <= hi]
            if not seg_pts:
                continue
            pieces = [[seg_pts[0]]]
            for q in seg_pts[1:]:
                if q - pieces[-1][-1] > RREPEAT_MAX_INTERNAL_GAP:
                    pieces.append([q])
                else:
                    pieces[-1].append(q)
            for pc in pieces:
                s, e = pc[0], pc[-1]
                length = e - s + 1
                if length < RREPEAT_MIN_FEATURE_LEN:
                    continue
                dens = len(pc) / length
                if dens < RREPEAT_MIN_COPY_DENSITY:
                    print(f"Warning: dropped mosaic-like R repeat copy [{s}, {e}] "
                          f"(len={length}, points={len(pc)}, density={dens:.3f})")
                    continue
                if length > RREPEAT_MAX_COPY_LEN:
                    print(f"Warning: R repeat copy [{s}, {e}] length {length} "
                          f"> {RREPEAT_MAX_COPY_LEN}, kept as a single feature")
                copies.append((s, e))

    return copies

def main(ref_fasta, features_file, query_fasta, output_gff, units_bed='named_units.bed'):
    """Main function"""
    # Read reference sequence
    ref_record = next(SeqIO.parse(ref_fasta, "fasta"))
    ref_seq = str(ref_record.seq).upper()
    
    # Build reference index
    print("Building reference index...")
    ref_index = build_ref_index(ref_seq)
    print(f"Reference fragments: {len(ref_index)}")
    
    # Read Features.txt file
    print("Parsing Features.txt file...")
    ref_features = parse_features_file(features_file)
    print(f"Found {len(ref_features)} features in Features.txt")
    
    # Read named units (named_units.bed)
    print(f"Reading named units from {units_bed}...")
    named_units = parse_named_units(units_bed)
    print(f"Found {sum(len(v) for v in named_units.values())} units in {len(named_units)} assemblies")
    
    # Build position to feature map
    get_feature_for_pos = build_ref_pos_to_feature_map(ref_features)
    
    # Process each sequence in query fasta
    print("Processing query sequences...")
    with open(output_gff, 'w') as out:
        out.write('##gff-version 3\n')
        
        for rec in SeqIO.parse(query_fasta, "fasta"):
            print(f"Processing {rec.id}...")
            seq = str(rec.seq).upper()
            seq_length = len(seq)
            
            # 该assembly在named_units.bed中的unit信息
            units = named_units.get(rec.id, [])
            
            # unit区间（1-based inclusive），用于过滤unit以外的匹配点（替代flank.bed）
            unit_regions = [(u['start'], u['end']) for u in units]
            
            # Get matches
            matches = query_sequence(seq, ref_index, 1)
            
            # Map and merge features
            query_features = map_and_merge_features(matches, get_feature_for_pos, unit_regions, seq)
            
            if not units:
                # named_units.bed中没有该assembly的unit信息时，退回旧方式：整条序列一个gene，编号全局递增
                print(f"Warning: no units found for {rec.id} in {units_bed}, falling back to whole-sequence gene")
                
                # Add sequence region information
                out.write(f'{rec.id}\tself_anno\tgene\t1\t{seq_length}\t.\t+\t.\tID={rec.id}\n')
                
                # Track occurrence order for each ID in this assembly
                id_counters = defaultdict(int)
                
                # Write features for this sequence
                for feature in query_features:
                    id_counters[feature['id']] += 1
                    attributes = f"ID={feature['id']} {id_counters[feature['id']]};Parent={rec.id}"
                    parts = [
                        rec.id,
                        "self_anno",
                        feature['type'],
                        str(feature['start']),
                        str(feature['end']+19),
                        ".",
                        "+",
                        ".",
                        attributes
                    ]
                    out.write('\t'.join(parts) + '\n')
                continue
            
            # 每个unit写一行gene，unit名字为Array_x_y或Left/Right_fragment
            for unit in units:
                out.write(f'{rec.id}\tself_anno\tgene\t{unit["start"]}\t{unit["end"]}\t.\t+\t.\tID={unit["name"]}\n')
            
            # 按feature起始位置将feature分配到对应的unit
            unit_features = {unit['name']: [] for unit in units}
            for feature in query_features:
                unit = get_unit_for_position(feature['start'], units)
                unit_features[unit['name']].append(feature)
            
            # 依次写每个unit内的feature，编号为 xxx y；unit内重复出现的同一feature记为 y.n
            for unit in units:
                feats = sorted(unit_features[unit['name']], key=lambda x: x['start'])
                
                # 统计该unit内每个feature基础名出现的次数
                base_counts = defaultdict(int)
                for feature in feats:
                    base_counts[feature['base']] += 1
                
                base_seen = defaultdict(int)
                for feature in feats:
                    base_seen[feature['base']] += 1
                    n = base_seen[feature['base']]
                    if base_counts[feature['base']] > 1:
                        # 同一个unit内重复的feature，编号为 y.n，如 R Repeat 3.1、R Repeat 3.2
                        feature_id = f"{feature['base']} {unit['y']}.{n}"
                    else:
                        feature_id = f"{feature['base']} {unit['y']}"
                    
                    attributes = f"ID={feature_id};Parent={unit['name']}"
                    parts = [
                        rec.id,
                        "self_anno",
                        feature['type'],
                        str(feature['start']),
                        str(feature['end']+19),
                        ".",
                        "+",
                        ".",
                        attributes
                    ]
                    out.write('\t'.join(parts) + '\n')
    
    print(f"GFF file written to {output_gff}")

if __name__ == "__main__":
    if len(sys.argv) < 5 or len(sys.argv) > 6:
        print("Usage: python generate_gff.py <ref_fasta> <features_file> <query_fasta> <output_gff> [units_bed]")
        print("  units_bed: named_units.bed defining the units (default: named_units.bed in current directory)")
        print("  注意：有了named_units.bed之后不再需要flank.bed，unit区间以外的匹配点会被自动过滤")
        sys.exit(1)
    
    ref_fasta = sys.argv[1]
    features_file = sys.argv[2]
    query_fasta = sys.argv[3]
    output_gff = sys.argv[4]
    units_bed = sys.argv[5] if len(sys.argv) == 6 else 'named_units.bed'
    
    if not os.path.isfile(units_bed):
        print(f"Error: units bed file not found: {units_bed}")
        sys.exit(1)
    
    main(ref_fasta, features_file, query_fasta, output_gff, units_bed)