import pandas as pd
import numpy as np
from collections import defaultdict
import sys

if len(sys.argv) != 2:
    sys.stderr.write("Usage: python testify_interval.py input.tsv\n")
    sys.exit(1)

THRESHOLD = 0.025  # 5%

# ============================
# 读取数据
# ============================
input_file = sys.argv[1]
df = pd.read_excel(input_file)

df["interval_list"] = df["intervals"].apply(
    lambda x: [int(i) for i in x.split("|")]
)

# ============================
# Step 1: 统计每个 assembly 的 interval 均值
# ============================
assembly_interval_values = defaultdict(lambda: defaultdict(list))

for _, row in df.iterrows():
    assembly = row["assembly"]
    start = row["start_index"]
    intervals = row["interval_list"]

    # intervals[0] 是占位 0
    for i in range(1, len(intervals)):
        assembly_pos = start + i - 1
        assembly_interval_values[assembly][assembly_pos].append(intervals[i])

assembly_interval_mean = {
    assembly: {
        pos: np.mean(vals)
        for pos, vals in pos_dict.items()
    }
    for assembly, pos_dict in assembly_interval_values.items()
}

# print(assembly_interval_mean)

# ============================
# Step 2: 计算每个 assembly 的 barcode 数
# ============================
assembly_length = (
    df.groupby("assembly")["end_index"].max()
).to_dict()

# ============================
# Step 3: 异常检测
# ============================
all_abnormal = []

for _, row in df.iterrows():
    read_id = row["read_id"]
    assembly = row["assembly"]
    start = row["start_index"]
    end = row["end_index"]
    intervals = row["interval_list"]

    mean_dict = assembly_interval_mean[assembly]
    asm_len = assembly_length[assembly]

    # ---- (1) barcode interval 异常 ----
    for i in range(1, len(intervals)):
        assembly_pos = start + i - 1
        mean_val = mean_dict[assembly_pos]

        if abs(intervals[i] - mean_val) > mean_val * THRESHOLD:
            all_abnormal.append({
                "type": "interval",
                "read_id": read_id,
                "assembly": assembly,
                "assembly_interval_index": assembly_pos,
                "interval": intervals[i],
                "mean_interval": mean_val
            })

    # ---- (2) first_barcode_to_start ----
    # print(f'mean_dict: {mean_dict}')
    # print(f'assembly: {assembly}')
    if start > 0:
        # print(f'start: {start}')
        prev_mean = mean_dict[start - 1]
        if row["first_barcode_to_start"] > prev_mean * (1 + THRESHOLD):
            all_abnormal.append({
                "type": "first_barcode_to_start",
                "read_id": read_id,
                "assembly": assembly,
                "first_barcode_to_start": row["first_barcode_to_start"],
                "threshold": prev_mean * (1 + THRESHOLD)
            })

    # ---- (3) last_barcode_to_end ----
    if end < asm_len:
        # print(f'end: {end-1}')
        # print(f'assembly: {assembly}')
        next_mean = mean_dict[end-1]
        if row["last_barcode_to_end"] > next_mean * (1 + THRESHOLD):
            all_abnormal.append({
                "type": "last_barcode_to_end",
                "read_id": read_id,
                "assembly": assembly,
                "last_barcode_to_end": row["last_barcode_to_end"],
                "threshold": next_mean * (1 + THRESHOLD)
            })

# ============================
# 输出
# ============================
pd.DataFrame(all_abnormal).to_csv(
    "abnormal.xls", sep="\t", index=False
)

print("All done (multi-assembly).")