# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np

# ======== 路径（按需改） ========
base_table_path = "/home/zjliang/users/zhengjiani/dynamics/data/cluster_md_results.xlsx"
out_path        = "/home/zjliang/users/zhengjiani/dynamics/data/cluster_md_results_final.xlsx"

# ANM
anm_singleAA_cluster_dir   = "/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/singleAA_cluster"
anm_singleAA_noncluster_dir= "/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/singleAA_noncluster"
anm_CC_cluster_dir         = "/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/CC_cluster"
anm_CC_noncluster_dir      = "/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/CC_noncluster"

# GNM
gnm_singleAA_cluster_dir   = "/home/zjliang/users/zhengjiani/dynamics/data/results/gnm_cluster/singleAA_cluster"
gnm_singleAA_noncluster_dir= "/home/zjliang/users/zhengjiani/dynamics/data/results/gnm_cluster/singleAA_noncluster"
gnm_CC_cluster_dir         = "/home/zjliang/users/zhengjiani/dynamics/data/results/gnm_cluster/CC_cluster"
gnm_CC_noncluster_dir      = "/home/zjliang/users/zhengjiani/dynamics/data/results/gnm_cluster/CC_noncluster"

# ======== 工具函数 ========
def read_sheet_any(path):
    """根据扩展名自动读取 CSV/Excel"""
    if path is None or not os.path.exists(path):
        return None
    ext = os.path.splitext(path)[1].lower()
    try:
        if ext == ".csv":
            for enc in ("utf-8","gbk","latin1"):
                try:
                    return pd.read_csv(path, encoding=enc, engine="c", low_memory=False)
                except Exception:
                    continue
            return None
        else:
            return pd.read_excel(path)
    except Exception as e:
        print(f"[WARN] 读文件失败 {path}: {e}")
        return None

def pick_col(df, candidates):
    """返回 df 中第一个存在的候选列名"""
    if df is None:
        return None
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None

def max_singleAA_metrics(file_path, model="ANM"):
    """
    singleAA 文件取最大值：
    - ANM: sq, sensor, effector, stiffness
    - GNM: sq, sensor, effector
    """
    df = read_sheet_any(file_path)
    if df is None or df.empty:
        return (np.nan, np.nan, np.nan, np.nan) if model=="ANM" else (np.nan, np.nan, np.nan)

    if model.upper() == "ANM":
        sq_col        = pick_col(df, ["ANM_sq","sq"])
        sensor_col    = pick_col(df, ["ANM_sensitivity","sensitivity","sensor"])
        effector_col  = pick_col(df, ["ANM_effectiveness","effectiveness","effector"])
        stiffness_col = pick_col(df, ["ANM_stiffness","stiffness"])
    else:  # GNM
        sq_col        = pick_col(df, ["GNM_sq","sq"])
        sensor_col    = pick_col(df, ["GNM_sensitivity","sensitivity","sensor"])
        effector_col  = pick_col(df, ["GNM_effectiveness","effectiveness","effector"])
        stiffness_col = None  # GNM 没有 stiffness

    def safe_max(col):
        if col is None: return np.nan
        s = pd.to_numeric(df[col], errors="coerce")
        return np.nanmax(s.values) if np.isfinite(s).any() else np.nan

    if model.upper() == "ANM":
        return safe_max(sq_col), safe_max(sensor_col), safe_max(effector_col), safe_max(stiffness_col)
    else:
        return safe_max(sq_col), safe_max(sensor_col), safe_max(effector_col)

def max_cc(path):
    """读取 CC 矩阵并取上三角(不含对角)最大值；处理标签列、重复标签、非方阵等情况。"""
    if path is None or not os.path.exists(path):
        return np.nan

    # 读取（优先 CSV，多编码兜底）
    df = None
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            df = pd.read_csv(path, encoding=enc, engine="c", low_memory=False)
            break
        except Exception:
            try:
                df = pd.read_csv(path, encoding=enc, engine="python", low_memory=False, on_bad_lines="skip")
                break
            except Exception:
                continue
    if df is None or df.empty:
        return np.nan

    # 若首列像标签（多数非数字），将其设为索引
    first_col = df.columns[0]
    if df[first_col].dtype == object:
        t = pd.to_numeric(df[first_col], errors="coerce")
        if t.isna().mean() > 0.5:
            df = df.set_index(first_col)

    # 去重行/列标签，避免 loc 后仍非方阵
    if df.index.has_duplicates:
        df = df.loc[~df.index.duplicated(keep="first")]
    if pd.Index(df.columns).has_duplicates:
        df = df.loc[:, ~pd.Index(df.columns).duplicated(keep="first")]

    # 用公共标签对齐为（尽量）方阵
    idx = df.index.astype(str)
    cols = df.columns.astype(str)
    common = idx.intersection(cols)

    if len(common) > 0:
        df = df.loc[common, common]
    else:
        # 没有公共标签：退化为数值子矩阵的左上角方阵
        df = df.apply(pd.to_numeric, errors="coerce")
        if df.size == 0:
            return np.nan
        r, c = df.shape
        n = min(r, c)
        df = df.iloc[:n, :n]

    # 数值化
    df = df.apply(pd.to_numeric, errors="coerce")
    arr = df.to_numpy(dtype=float, copy=False)

    # 仍然防守：若偶发非方阵，裁成左上角方阵
    r, c = arr.shape
    n = min(r, c)
    if n == 0:
        return np.nan
    if r != n or c != n:
        arr = arr[:n, :n]

    # 上三角（不含对角）
    iu = np.triu_indices(n, k=1)
    sub = arr[iu]
    sub = sub[np.isfinite(sub)]
    return float(sub.max()) if sub.size else np.nan



# ======== 主流程 ========
base = pd.read_excel(base_table_path)

new_cols = [
    "anm_cluster_sq","anm_cluster_sensor","anm_cluster_effector","anm_cluster_stiffness","anm_cluster_top3cc","anm_cluster_5percc",
    "gnm_cluster_sq","gnm_cluster_sensor","gnm_cluster_effector","gnm_cluster_top3cc","gnm_cluster_5percc",
    "anm_noncluster_sq","anm_noncluster_sensor","anm_noncluster_effector","anm_noncluster_stiffness","anm_noncluster_top3cc","anm_noncluster_5percc",
    "gnm_noncluster_sq","gnm_noncluster_sensor","gnm_noncluster_effector","gnm_noncluster_top3cc","gnm_noncluster_5percc"
]
for c in new_cols:
    base[c] = np.nan

for i, row in base.iterrows():
    pdb = str(row["PDB"]).strip()
    cluster = str(row["Cluster_New"]).strip()

    # --- ANM singleAA cluster ---
    f = os.path.join(anm_singleAA_cluster_dir, f"{cluster}.csv")
    sq, sensor, eff, stiff = max_singleAA_metrics(f, "ANM")
    base.loc[i, ["anm_cluster_sq","anm_cluster_sensor","anm_cluster_effector","anm_cluster_stiffness"]] = [sq, sensor, eff, stiff]

    # --- ANM singleAA noncluster ---
    f = os.path.join(anm_singleAA_noncluster_dir, f"{pdb}.csv")
    sq, sensor, eff, stiff = max_singleAA_metrics(f, "ANM")
    base.loc[i, ["anm_noncluster_sq","anm_noncluster_sensor","anm_noncluster_effector","anm_noncluster_stiffness"]] = [sq, sensor, eff, stiff]

    # --- ANM CC cluster ---
    f_top3 = os.path.join(anm_CC_cluster_dir, f"{cluster}_top3_cluster_cc.csv")
    f_5per = os.path.join(anm_CC_cluster_dir, f"{cluster}_5_per_cluster_cc.csv")
    base.loc[i, "anm_cluster_top3cc"] = max_cc(f_top3)
    base.loc[i, "anm_cluster_5percc"] = max_cc(f_5per)

    # --- ANM CC noncluster ---
    f_top3 = os.path.join(anm_CC_noncluster_dir, f"{pdb}_top3_noncluster_cc.csv")
    f_5per = os.path.join(anm_CC_noncluster_dir, f"{pdb}_5_per_noncluster_cc.csv")
    base.loc[i, "anm_noncluster_top3cc"] = max_cc(f_top3)
    base.loc[i, "anm_noncluster_5percc"] = max_cc(f_5per)

    # --- GNM singleAA cluster ---
    f = os.path.join(gnm_singleAA_cluster_dir, f"{cluster}.csv")
    sq, sensor, eff = max_singleAA_metrics(f, "GNM")
    base.loc[i, ["gnm_cluster_sq","gnm_cluster_sensor","gnm_cluster_effector"]] = [sq, sensor, eff]

    # --- GNM singleAA noncluster ---
    f = os.path.join(gnm_singleAA_noncluster_dir, f"{pdb}.csv")
    sq, sensor, eff = max_singleAA_metrics(f, "GNM")
    base.loc[i, ["gnm_noncluster_sq","gnm_noncluster_sensor","gnm_noncluster_effector"]] = [sq, sensor, eff]

    # --- GNM CC cluster ---
    f_top3 = os.path.join(gnm_CC_cluster_dir, f"{cluster}_top3_cluster_cc.csv")
    f_5per = os.path.join(gnm_CC_cluster_dir, f"{cluster}_5_per_cluster_cc.csv")
    base.loc[i, "gnm_cluster_top3cc"] = max_cc(f_top3)
    base.loc[i, "gnm_cluster_5percc"] = max_cc(f_5per)

    # --- GNM CC noncluster ---
    f_top3 = os.path.join(gnm_CC_noncluster_dir, f"{pdb}_top3_noncluster_cc.csv")
    f_5per = os.path.join(gnm_CC_noncluster_dir, f"{pdb}_5_per_noncluster_cc.csv")
    base.loc[i, "gnm_noncluster_top3cc"] = max_cc(f_top3)
    base.loc[i, "gnm_noncluster_5percc"] = max_cc(f_5per)

# 导出
base.to_excel(out_path, index=False)
print("Done ->", out_path)
