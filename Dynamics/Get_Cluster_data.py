# -*- coding: utf-8 -*-

import os
import re
import glob
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

# ===== EXCEL_Path =====
cluster_sites_inf_path = r"/home/zjliang/users/zhengjiani/dynamics/data/cluster_sites_for_MD.xlsx"

# ANM
anm_singleAA_path = r"/home/zjliang/users/zhengjiani/dynamics/data/results/anm/singleAA_Data"
anm_CC_path       = r"/home/zjliang/users/zhengjiani/dynamics/data/results/anm/CC_Matrix"
anm_out_single_cluster = r"/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/singleAA_cluster"
anm_out_single_nonclu = r"/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/singleAA_noncluster"
anm_out_cc_cluster    = r"/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/CC_cluster"
anm_out_cc_nonclu     = r"/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/CC_noncluster"

# GNM
gnm_singleAA_path = r"/home/zjliang/users/zhengjiani/dynamics/data/results/gnm/singleAA_Data"
gnm_CC_path       = r"/home/zjliang/users/zhengjiani/dynamics/data/results/gnm/CC_Matrix"
gnm_out_single_cluster = r"/home/zjliang/users/zhengjiani/dynamics/data/results/gnm_cluster/singleAA_cluster"
gnm_out_single_nonclu = r"/home/zjliang/users/zhengjiani/dynamics/data/results/gnm_cluster/singleAA_noncluster"
gnm_out_cc_cluster    = r"/home/zjliang/users/zhengjiani/dynamics/data/results/gnm_cluster/CC_cluster"
gnm_out_cc_nonclu     = r"/home/zjliang/users/zhengjiani/dynamics/data/results/gnm_cluster/CC_noncluster"

for p in [
    anm_out_single_cluster, anm_out_single_nonclu, anm_out_cc_cluster, anm_out_cc_nonclu,
    gnm_out_single_cluster, gnm_out_single_nonclu, gnm_out_cc_cluster, gnm_out_cc_nonclu
]:
    os.makedirs(p, exist_ok=True)

# ===== 参数 =====
WRITE_EXCEL = False           # True: 写 xlsx；False: 写 csv（更快）
EXCEL_ENGINE = "xlsxwriter"   # 写 Excel 的引擎
N_WORKERS = max(os.cpu_count() - 1, 1)  # 并行进程数；设为 1 则串行
KEEP_NONCLUSTER_CC = True     # 导出 noncluster(行) vs interface(列) 的 CC 子矩阵

CC_BANDS = ["top3", "5_per", "5_20_per", "20_50_per", "greate_60_per"]

# ===== Functions =====
def split_sites(site_str):
    """'A, B , C' -> {'A','B','C'}；空值安全"""
    if pd.isna(site_str) or not str(site_str).strip():
        return set()
    return {s.strip() for s in str(site_str).split(',') if s.strip()}

def read_singleAA_file(base_dir, pdb):
    """读取 singleAA csv；找不到返回 None"""
    pattern = os.path.join(base_dir, f"{pdb}_*_singleAA_data.csv")
    files = glob.glob(pattern)
    if not files:
        return None
    df = pd.read_csv(files[0])
    if "Res_Info" not in df.columns:
        return None
    df["Res_Info"] = df["Res_Info"].astype(str).str.strip()
    return df

def _read_csv_guess_sep(path):
    """优先逗号读取；若只有 1 列则回退为制表符"""
    df = pd.read_csv(path, index_col=0)
    if df.shape[1] <= 1:
        df = pd.read_csv(path, sep="\t", index_col=0, engine="c")
    return df

def read_cc_files(base_dir, pdb):

    out = {}
    # 匹配：..._{band}_cc.csv 结尾
    band_re = re.compile(r"_(%s)_cc\.csv$" % "|".join(map(re.escape, CC_BANDS)))
    for fp in sorted(glob.glob(os.path.join(base_dir, f"{pdb}_*_cc.csv"))):
        m = band_re.search(os.path.basename(fp))
        if not m:
            continue
        band = m.group(1)
        if band in out:           # 已读取过则跳过（通常不会重复）
            continue
        df = _read_csv_guess_sep(fp)
        df.index = df.index.astype(str)
        df.columns = df.columns.astype(str)
        out[band] = df
    return out  # 可能为空 {}

def write_table(df, out_path):
    if WRITE_EXCEL:
        df.to_excel(out_path, index=False if isinstance(df, pd.DataFrame) and df.index.name is None else True,
                    engine=EXCEL_ENGINE)
    else:
        df.to_csv(out_path.replace(".xlsx", ".csv"), index=False)

def write_matrix(df, out_path):
    if WRITE_EXCEL:
        df.to_excel(out_path, engine=EXCEL_ENGINE)
    else:
        df.to_csv(out_path.replace(".xlsx", ".csv"))

# —— 位点规范化（仅 CC 用）：把 “..._A_A117” -> “..._A_X117”
_site_pat = re.compile(r"^(?P<pfx>.+?)_(?P<chain>[^_]+)_(?P<aa>[A-Za-z])?(?P<num>\d+)$")
def to_X_form(site: str) -> str:
    if not isinstance(site, str):
        return site
    site = site.strip()
    m = _site_pat.match(site)
    if m:
        return f"{m.group('pfx')}_{m.group('chain')}_X{m.group('num')}"
    m2 = re.search(r"(\d+)$", site)
    if m2:
        num = m2.group(1)
        parts = site.split("_")
        if len(parts) >= 3:
            return "_".join([parts[0], parts[1], f"X{num}"])
    return site

def normalize_sites(sites):
    """批量 X 化 + 去重"""
    return {to_X_form(s) for s in sites if isinstance(s, str) and s.strip()}

def normalize_index_and_columns(M: pd.DataFrame) -> pd.DataFrame:
    """矩阵行列统一 X 化"""
    M = M.copy()
    M.index = [to_X_form(str(i)) for i in M.index]
    M.columns = [to_X_form(str(c)) for c in M.columns]
    return M

def safe_submatrix(M: pd.DataFrame, rows, cols):
    """只取存在的行/列；空则返回 None"""
    r = sorted(set(rows) & set(M.index))
    c = sorted(set(cols) & set(M.columns))
    if not r or not c:
        return None
    return M.loc[r, c]

# ===== 核心处理 =====
def process_one_pdb(pdb, clusters,
                    all_sites_orig,        # 原始位点并集（singleAA 用）
                    all_sites_X,           # X 形式位点并集（判定 noncluster 行）
                    all_interface_X,       # X 形式接口位点并集（判定 noncluster 列）
                    single_dir, cc_dir,
                    out_single_cluster, out_single_nonclu,
                    out_cc_cluster, out_cc_nonclu,
                    tag="ANM"):
    # —— singleAA：按原始位点名过滤（不做 X 化）
    df_single = read_singleAA_file(single_dir, pdb)
    if df_single is not None:
        # noncluster（singleAA）
        mask_non = ~df_single["Res_Info"].isin(all_sites_orig)
        non_df = df_single.loc[mask_non]
        if not non_df.empty:
            write_table(non_df, os.path.join(out_single_nonclu, f"{pdb}.xlsx"))
        # 每个 cluster 的 singleAA
        for cluster_id, cset_orig, _cset_X, _ifset_X in clusters:
            if not cset_orig:
                continue
            sub = df_single.loc[df_single["Res_Info"].isin(cset_orig)]
            if not sub.empty:
                write_table(sub, os.path.join(out_single_cluster, f"{cluster_id}.xlsx"))

    # —— CC：对每个 band（top3 / 5% / 5~20% / 20~50% / ≥60%）都导出
    cc_dict = read_cc_files(cc_dir, pdb)  # {band: df}
    if not cc_dict:
        return f"[{tag}] {pdb} done (no CC)"

    # X 化
    for b in list(cc_dict.keys()):
        cc_dict[b] = normalize_index_and_columns(cc_dict[b])

    # noncluster：行 = 非 cluster 的 X 位点；列 = 接口 X 位点
    if KEEP_NONCLUSTER_CC:
        for band, M in cc_dict.items():
            non_rows = sorted(set(M.index) - set(all_sites_X))
            non_cols = sorted(set(all_interface_X) & set(M.columns))
            sub = safe_submatrix(M, non_rows, non_cols)
            if sub is not None and not sub.empty:
                write_matrix(sub, os.path.join(
                    out_cc_nonclu, f"{pdb}_{band}_noncluster_vs_interface_cc.xlsx"
                ))

    # 每个 cluster：行 = 该 cluster 的 X 位点；列 = 该 cluster 的接口 X 位点
    for cluster_id, _cset_orig, cset_X, ifset_X in clusters:
        if not cset_X or not ifset_X:
            continue
        for band, M in cc_dict.items():
            sub = safe_submatrix(M, cset_X, ifset_X)
            if sub is not None and not sub.empty:
                write_matrix(sub, os.path.join(
                    out_cc_cluster, f"{cluster_id}_{band}_cluster_vs_interface_cc.xlsx"
                ))

    return f"[{tag}] {pdb} done"

# ===== 主流程：读基表 → 构造映射 → 并行处理 =====
# 需要的列：PDB, Cluster_New, Site_New, Interface_Site
cluster_sites_inf = pd.read_excel(cluster_sites_inf_path)

# 构造：{pdb: [(cluster_id, cluster_sites_orig, cluster_sites_X, interface_sites_X), ...]}
clusters_by_pdb = {}
all_sites_by_pdb_orig = {}
all_sites_by_pdb_X = {}
all_interface_by_pdb_X = {}

for pdb, dfp in cluster_sites_inf.groupby("PDB", sort=False):
    lst = []
    union_sites_orig = set()
    union_sites_X = set()
    union_interface_X = set()

    for _, row in dfp.iterrows():
        cluster_id = str(row["Cluster_New"])
        cluster_sites_orig = split_sites(row["Site_New"])
        interface_sites_orig = split_sites(row.get("Interface_Site", None))
        cluster_sites_X   = normalize_sites(cluster_sites_orig)
        interface_sites_X = normalize_sites(interface_sites_orig)

        if cluster_sites_orig:
            lst.append((cluster_id, cluster_sites_orig, cluster_sites_X, interface_sites_X))
            union_sites_orig |= cluster_sites_orig
            union_sites_X    |= cluster_sites_X
            union_interface_X |= interface_sites_X

    clusters_by_pdb[pdb]        = lst
    all_sites_by_pdb_orig[pdb]  = union_sites_orig
    all_sites_by_pdb_X[pdb]     = union_sites_X
    all_interface_by_pdb_X[pdb] = union_interface_X

pdbs = list(clusters_by_pdb.keys())

def run_anm_for(pdb):
    return process_one_pdb(
        pdb,
        clusters_by_pdb[pdb],
        all_sites_by_pdb_orig[pdb],
        all_sites_by_pdb_X[pdb],
        all_interface_by_pdb_X[pdb],
        anm_singleAA_path,
        anm_CC_path,
        anm_out_single_cluster,
        anm_out_single_nonclu,
        anm_out_cc_cluster,
        anm_out_cc_nonclu,
        tag="ANM"
    )

def run_gnm_for(pdb):
    return process_one_pdb(
        pdb,
        clusters_by_pdb[pdb],
        all_sites_by_pdb_orig[pdb],
        all_sites_by_pdb_X[pdb],
        all_interface_by_pdb_X[pdb],
        gnm_singleAA_path,
        gnm_CC_path,
        gnm_out_single_cluster,
        gnm_out_single_nonclu,
        gnm_out_cc_cluster,
        gnm_out_cc_nonclu,
        tag="GNM"
    )

# —— 并行（设 N_WORKERS=1 可改为串行）
if N_WORKERS > 1:
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = [ex.submit(run_anm_for, pdb) for pdb in pdbs] + \
                  [ex.submit(run_gnm_for, pdb) for pdb in pdbs]
        for f in as_completed(futures):
            try:
                print(f.result())
            except Exception as e:
                print("Error:", e)
else:
    for pdb in pdbs:
        print(run_anm_for(pdb))
    for pdb in pdbs:
        print(run_gnm_for(pdb))
