# -*- coding: utf-8 -*-
import re
import shutil
import pandas as pd
from pathlib import Path

# ========= 路径配置（按你的目录改即可） =========
CLUSTER_IN_DIR    = Path("/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/CC_cluster")
NONCLUSTER_IN_DIR = Path("/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/CC_noncluster")
CLASS_ROOT        = Path("/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/Classify")
META_XLSX         = Path("/home/zjliang/users/zhengjiani/dynamics/data/cluster_summary_cc30_0825.xlsx")

# ========= 窗口类型（桶）的识别规则 =========
BUCKET_PATTERNS = [
    ("anm_top3"          , re.compile(r"_top3_")),
    ("anm_5_per"         , re.compile(r"_5_per_")),
    ("anm_5_20_per"      , re.compile(r"_5_20_per_")),
    ("anm_20_50_per"     , re.compile(r"_20_50_per_")),
    ("anm_greate_60_per" , re.compile(r"_greate_60_per_")),
]

# ========= cluster 细分到 Pos_Type 的映射（输出子目录名） =========
POS_DIR = {
    "Interface":      "Interface",
    "Non-Interface":  "Non_Interface",
    "Cross-Interface":"Cross_Interface",
}

# 从文件名抓 cluster id，如：Cluster.29.1_5_per_cluster_vs_interface_cc.csv -> Cluster.29.1
CLUSTER_ID_RE = re.compile(r"^(Cluster\.\d+(?:\.\d+)?)_")

def bucket_by_name(fname: str) -> str | None:
    """根据文件名匹配窗口类型桶名；匹配失败返回 None"""
    for bucket, pat in BUCKET_PATTERNS:
        if pat.search(fname):
            return bucket
    return None

def load_pos_map(xlsx_path: Path) -> dict:
    """Excel → {Cluster_New: Pos_Type}"""
    df = pd.read_excel(xlsx_path, usecols=["Cluster_New", "Pos_Type"]).dropna(subset=["Cluster_New"])
    df["Cluster_New"] = df["Cluster_New"].astype(str).str.strip()
    df["Pos_Type"]    = df["Pos_Type"].astype(str).str.strip()
    return dict(zip(df["Cluster_New"], df["Pos_Type"]))

def ensure_layout():
    """创建输出目录框架：每个桶下含 Interface/Non_Interface/Cross_Interface/noncluster 四类"""
    CLASS_ROOT.mkdir(parents=True, exist_ok=True)
    for bucket, _ in BUCKET_PATTERNS:
        base = CLASS_ROOT / bucket
        (base / "Interface").mkdir(parents=True, exist_ok=True)
        (base / "Non_Interface").mkdir(parents=True, exist_ok=True)
        (base / "Cross_Interface").mkdir(parents=True, exist_ok=True)
        (base / "noncluster").mkdir(parents=True, exist_ok=True)

def copy_cluster_files(pos_map: dict):
    """处理 cluster：按桶 + Pos_Type 复制"""
    n = 0
    for fp in sorted(CLUSTER_IN_DIR.glob("*.csv")):
        bucket = bucket_by_name(fp.name)
        if not bucket:
            continue                               # 约定：不会出现 Unknown；不匹配则跳过
        m = CLUSTER_ID_RE.match(fp.name)
        if not m:
            continue                               # 没抓到 Cluster id 则跳过
        cid = m.group(1)
        pos = pos_map.get(cid)
        if not pos:
            continue                               # 映射里找不到则跳过
        subdir = POS_DIR.get(pos)
        if not subdir:
            continue                               # 不支持的 Pos_Type 跳过

        dst = CLASS_ROOT / bucket / subdir / fp.name
        if not dst.exists():
            shutil.copy2(fp, dst)                 # 复制（保留元数据）
            n += 1
    print(f"[cluster] 已复制 {n} 个文件")

def copy_noncluster_files():
    """处理 noncluster：仅按桶分类，全部放到桶下的 noncluster/ 目录"""
    n = 0
    for fp in sorted(NONCLUSTER_IN_DIR.glob("*.csv")):
        bucket = bucket_by_name(fp.name)
        if not bucket:
            continue                               # 约定：不会出现 Unknown；不匹配则跳过
        dst = CLASS_ROOT / bucket / "noncluster" / fp.name
        if not dst.exists():
            shutil.copy2(fp, dst)                 # 复制
            n += 1
    print(f"[noncluster] 已复制 {n} 个文件")

def main():
    ensure_layout()                       # 准备输出目录结构
    pos_map = load_pos_map(META_XLSX)     # 读取 {Cluster: Pos_Type}
    copy_cluster_files(pos_map)           # 复制 cluster 文件
    copy_noncluster_files()               # 复制 noncluster 文件
    print(f"[完成] 输出目录：{CLASS_ROOT}")

if __name__ == "__main__":
    main()