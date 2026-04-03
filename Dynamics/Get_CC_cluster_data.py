# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

# ===== 配置 =====
CLASSIFY_ROOT = Path("/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/Classify")
OUT_ROOT      = Path("/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/Classify_out")

ABS_VAL = False          # 是否取绝对值
CHUNK_N = 1_000_000      # 分块写出行数

# 大类目录
BUCKETS = ["anm_top3", "anm_5_per", "anm_5_20_per", "anm_20_50_per", "anm_greate_60_per"]
# 子类目录
SUBFOLDERS = ["Interface", "Non_Interface", "Cross_Interface", "noncluster"]

def read_matrix(csv_path: Path) -> np.ndarray:
    """读取矩阵并返回 numpy 数组"""
    df = pd.read_csv(csv_path, sep=r"[,\t]+", engine="python", header=0, index_col=0)
    df = df.loc[:, ~df.columns.astype(str).str.contains(r"^Unnamed", na=False)]
    df = df.dropna(how="all").dropna(how="all", axis=1)
    df.index   = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    common = df.index.intersection(df.columns)
    if len(common) > 0:
        df = df.loc[common, common]
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all").dropna(how="all", axis=1)
    return df.to_numpy(dtype=float, copy=False)

def write_all_values(files, out_csv, abs_val=False, chunk_n=1_000_000):
    """遍历所有文件，抽出数值并汇总到单列 CSV"""
    if out_csv.exists():
        out_csv.unlink()
    header_written = False

    for fp in tqdm(files, desc=f"Processing {out_csv.stem}", unit="file"):
        try:
            mat = read_matrix(fp)
            vals = mat.ravel()
            vals = vals[~np.isnan(vals)]
            if abs_val:
                vals = np.abs(vals)
            if vals.size == 0:
                continue

            start = 0
            while start < vals.size:
                end = min(start + chunk_n, vals.size)
                pd.DataFrame({"value": vals[start:end]}).to_csv(
                    out_csv, mode="a", index=False, header=not header_written
                )
                header_written = True
                start = end
        except Exception as e:
            print(f"[警告] {fp.name} 处理失败：{e}")

def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    for bucket in BUCKETS:
        for sub in SUBFOLDERS:
            in_dir = CLASSIFY_ROOT / bucket / sub
            out_csv = OUT_ROOT / f"{bucket}_{sub.lower()}.csv"

            files = sorted(in_dir.glob("*.csv"))
            if not files:
                print(f"[提示] 目录为空：{in_dir}")
                continue

            write_all_values(files, out_csv, abs_val=ABS_VAL, chunk_n=CHUNK_N)
            print(f"[完成] {bucket}/{sub} -> {out_csv}")

if __name__ == "__main__":
    main()
