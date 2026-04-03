# -*- coding: utf-8 -*-
import re
import shutil
import pandas as pd
from pathlib import Path

# ========= 路径配置 =========
SINGLEAA_CLUSTER_DIR = Path("/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/singleAA_cluster")
SINGLEAA_NONCLU_DIR  = Path("/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/singleAA_noncluster")
CLASS_ROOT           = Path("/home/zjliang/users/zhengjiani/dynamics/data/results/anm_cluster/Classify")
META_XLSX            = Path("/home/zjliang/users/zhengjiani/dynamics/data/cluster_summary_cc30_0825.xlsx")

# ========= Pos_Type → 输出子目录 =========
POS_DIR = {
    "Interface":       "Interface",
    "Non-Interface":   "Non_Interface",
    "Cross-Interface": "Cross_Interface",
}

# 匹配 Cluster ID，例如 Cluster.29.1.xlsx
CLUSTER_ID_RE = re.compile(r"(Cluster\.\d+(?:\.\d+)?)")

def load_pos_map(xlsx_path: Path) -> dict:
    """读取 Excel，返回 {Cluster_New: Pos_Type}"""
    df = pd.read_excel(xlsx_path, usecols=["Cluster_New", "Pos_Type"]).dropna(subset=["Cluster_New"])
    df["Cluster_New"] = df["Cluster_New"].astype(str).str.strip()
    df["Pos_Type"]    = df["Pos_Type"].astype(str).str.strip()
    return dict(zip(df["Cluster_New"], df["Pos_Type"]))

def ensure_layout():
    """创建 singleAA 分类目录"""
    saa = CLASS_ROOT / "singleAA"
    for sub in ["Interface", "Non_Interface", "Cross_Interface", "noncluster"]:
        path = saa / sub
        path.mkdir(parents=True, exist_ok=True)
        print(f"[DIR] ensured {path}")

def iter_singleaa_files(root: Path):
    """遍历 singleAA 目录下的文件（支持 csv/xlsx）"""
    for ext in ("*.xlsx", "*.csv"):
        for fp in root.glob(ext):
            if fp.is_file():
                yield fp

def copy_singleaa_cluster_files(pos_map: dict):
    """处理 singleAA cluster：按 Pos_Type 分类"""
    n = 0
    out_root = CLASS_ROOT / "singleAA"
    for fp in sorted(iter_singleaa_files(SINGLEAA_CLUSTER_DIR)):
        m = CLUSTER_ID_RE.search(fp.name)
        if not m:
            continue
        cid = m.group(1)
        pos = pos_map.get(cid)
        if not pos:
            continue
        subdir = POS_DIR.get(pos)
        if not subdir:
            continue

        dst = out_root / subdir / fp.name
        if not dst.exists():
            shutil.copy2(fp, dst)
            n += 1
    print(f"[singleAA cluster] copied {n} files")

def copy_singleaa_noncluster_files():
    """处理 singleAA noncluster：全部放到 noncluster 目录"""
    n = 0
    out_root = CLASS_ROOT / "singleAA" / "noncluster"
    for fp in sorted(iter_singleaa_files(SINGLEAA_NONCLU_DIR)):
        dst = out_root / fp.name
        if not dst.exists():
            shutil.copy2(fp, dst)
            n += 1
    print(f"[singleAA noncluster] copied {n} files")

def main():
    ensure_layout()                       # 确保目录一定创建
    pos_map = load_pos_map(META_XLSX)     # 读 Pos_Type 映射
    copy_singleaa_cluster_files(pos_map)  # 复制 cluster
    copy_singleaa_noncluster_files()      # 复制 noncluster
    print(f"[DONE] output root: {CLASS_ROOT/'singleAA'}")

if __name__ == "__main__":
    main()
