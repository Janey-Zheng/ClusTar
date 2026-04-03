# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

# ========= 配置 =========
# 中文：singleAA 的根目录（其下每个子文件夹都要做）
BASE_DIR = Path("/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/Classify/singleAA")
# 中文：输出目录（所有结果集中放到这里）
OUT_DIR  = BASE_DIR / "aggregated"

# 中文：要汇总的列
TARGET_COLS = ["ANM_sq", "ANM_effectiveness", "ANM_sensitivity", "ANM_stiffness"]
# 中文：块写入的行数，用于降低内存占用
CHUNK_N = 1_000_000


# ========= 工具函数 =========
def iter_files(root: Path):
    """中文：遍历目录下所有 csv/xlsx 文件"""
    for pat in ("*.csv", "*.xlsx"):
        for fp in sorted(root.glob(pat)):
            if fp.is_file():
                yield fp

def read_needed(fp: Path, cols: list[str]) -> pd.DataFrame:
    """中文：只读取并返回需要的列；自动处理分隔符/Excel"""
    if fp.suffix.lower() == ".csv":
        df = pd.read_csv(fp, sep=r"[,\t]+", engine="python")
    else:
        df = pd.read_excel(fp)
    df.columns = df.columns.astype(str).str.strip()
    keep = [c for c in cols if c in df.columns]
    if not keep:
        return pd.DataFrame(columns=cols)
    out = df[keep].apply(pd.to_numeric, errors="coerce")
    for c in cols:
        if c not in out.columns:
            out[c] = np.nan
    return out[cols]

def ensure_parent(path: Path):
    """中文：确保输出文件父目录存在"""
    path.parent.mkdir(parents=True, exist_ok=True)

def append_series_csv(values: np.ndarray, out_csv: Path, colname: str, header_written: dict):
    """中文：把一维数组分块追加写入 CSV；首块写表头"""
    ensure_parent(out_csv)
    start, n = 0, values.size
    while start < n:
        end = min(start + CHUNK_N, n)
        pd.DataFrame({colname: values[start:end]}).to_csv(
            out_csv,
            mode="a",
            index=False,
            header=(not header_written.get(out_csv, False))
        )
        header_written[out_csv] = True
        start = end

def process_one_dir(in_dir: Path):
    """中文：处理一个子目录，输出四个单列 CSV"""
    tag = in_dir.name  # 中文：用子目录名作为文件前缀
    out_map = {
        "ANM_sq":            OUT_DIR / f"{tag}_ANM_sq.csv",
        "ANM_effectiveness": OUT_DIR / f"{tag}_ANM_effectiveness.csv",
        "ANM_sensitivity":   OUT_DIR / f"{tag}_ANM_sensitivity.csv",
        "ANM_stiffness":     OUT_DIR / f"{tag}_ANM_stiffness.csv",
    }
    # 中文：清理旧结果
    for p in out_map.values():
        if p.exists():
            p.unlink()

    header_written = {}
    files = list(iter_files(in_dir))
    if not files:
        print(f"[提示] 目录为空：{in_dir}")
        return

    for fp in tqdm(files, desc=f"Aggregating {tag}", unit="file"):
        try:
            df = read_needed(fp, TARGET_COLS)
            if df.empty:
                continue
            for col, out_csv in out_map.items():
                vals = df[col].to_numpy(dtype=float, copy=False)
                vals = vals[~np.isnan(vals)]
                if vals.size:
                    append_series_csv(vals, out_csv, col, header_written)
        except Exception as e:
            print(f"[警告] 处理失败 {fp.name}: {e}")

    # 中文：简要统计
    for col, out_csv in out_map.items():
        if out_csv.exists():
            n = max(0, sum(1 for _ in open(out_csv, "r")) - 1)
            print(f"[完成] {tag}: {col} -> {out_csv}  共 {n} 条")

def main():
    """中文：枚举 BASE_DIR 下的所有一级子目录，逐个处理"""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    subdirs = [p for p in sorted(BASE_DIR.iterdir()) if p.is_dir()]
    if not subdirs:
        print(f"[提示] 未找到子目录：{BASE_DIR}")
        return
    for d in subdirs:
        process_one_dir(d)

if __name__ == "__main__":
    main()
