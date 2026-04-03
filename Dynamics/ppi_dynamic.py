
# -*- coding: utf-8 -*-
"""
ppi_pipeline.py
---------------------------------
Unified entry point: Integrates ANM, splitting, classification, aggregation, and final table into a single script with CLI options.

Usage examples:
1) Run ANM calculation (reads pdb/cif with chain pairs, outputs singleAA + CC/PRS)
   python ppi_pipeline.py anm --pdb-dir /path/to/pdb --ppi-info /path/to/dynamic_calc_info.txt \
       --out-singleAA /path/to/results/anm/singleAA_Data \
       --out-cc /path/to/results/anm/CC_Matrix \
       --out-prs /path/to/results/anm/PRS_Matrix

2) Split singleAA/CC into cluster vs noncluster according to cluster/interface info (supports ANM and GNM)
   python ppi_pipeline.py split --cluster-sites-xlsx /path/to/cluster_sites_for_MD.xlsx \
       --anm-singleAA /path/to/results/anm/singleAA_Data \
       --anm-cc /path/to/results/anm/CC_Matrix \
       --gnm-singleAA /path/to/results/gnm/singleAA_Data \
       --gnm-cc /path/to/results/gnm/CC_Matrix \
       --anm-out-root /path/to/results/anm_cluster \
       --gnm-out-root /path/to/results/gnm_cluster

3) Classification (organize CC/noncluster files by mode bucket and Pos_Type; also classify singleAA)
   python ppi_pipeline.py classify-cc --cc-cluster /path/to/results/anm_cluster/CC_cluster \
       --cc-noncluster /path/to/results/anm_cluster/CC_noncluster \
       --class-root /path/to/results/anm_cluster/Classify \
       --meta-xlsx /path/to/cluster_summary_cc30_0825.xlsx
   python ppi_pipeline.py classify-singleaa --saa-cluster /path/to/results/anm_cluster/singleAA_cluster \
       --saa-noncluster /path/to/results/anm_cluster/singleAA_noncluster \
       --class-root /path/to/results/anm_cluster/Classify \
       --meta-xlsx /path/to/cluster_summary_cc30_0825.xlsx

4) Aggregation:
   - CC： CSV
     python ppi_pipeline.py aggregate-cc --classify-root /path/to/results/anm_cluster/Classify \
         --out-root /path/to/results/anm_cluster/Classify_out --abs-val

   - singleAA： CSV
     python ppi_pipeline.py aggregate-singleaa --saa-classify /path/to/results/anm_cluster/Classify/singleAA

5) Generate final table (merge aggregated singleAA and CC statistics into base table):
   python ppi_pipeline.py final --base-xlsx /path/to/cluster_md_results.xlsx \
       --out-xlsx /path/to/cluster_md_results_final.xlsx \
       --anm-root /path/to/results/anm_cluster \
       --gnm-root /path/to/results/gnm_cluster

6) Run the entire pipeline with one command (existing outputs are skipped based on each step's logic):
   python ppi_pipeline.py all --()
"""
import os
import re
import math
import glob
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

# ----------------------------- Common utilities -----------------------------
def ensure_dir(p: Path | str):
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p

# ============================= 1) ANM  =============================
# Based on MD_ANM.py core implementation (supports PDB and mmCIF; mmCIF uses segname for chain selection)
# ： MD_ANM.py

try:
    from prody import ANM, parsePDB, parseMMCIF, calcSqFlucts, calcMechStiff, calcCrossCorr, calcPerturbResponse
except Exception as _e:
    ANM = parsePDB = parseMMCIF = calcSqFlucts = calcMechStiff = calcCrossCorr = calcPerturbResponse = None

three_to_one = {'ALA': 'A','ARG': 'R','ASN': 'N','ASP': 'D','CYS': 'C','GLY': 'G','GLN': 'Q','GLU': 'E','HIS': 'H',
                'ILE': 'I','LEU': 'L','LYS': 'K','MET': 'M','PRO': 'P','PHE': 'F','SER': 'S','THR': 'T','TYR': 'Y','TRP': 'W','VAL': 'V'}
three_amino_acids = list(three_to_one.keys())

def _make_matrix_to_df(matrix, res_info):
    return pd.DataFrame(matrix, index=res_info, columns=res_info)

def _anm_core_build(structure, ca_selector: str):
    """ (calphas, res_info)"""
    calphas = structure.select(ca_selector)
    if calphas is None or calphas.numAtoms() == 0:
        raise ValueError(f"No CA atoms found by selector: {ca_selector}")
    res = calphas.getResnames()
    res = [three_to_one.get(r, r) for r in res]
    res_num = pd.DataFrame(calphas.getResnums())
    return calphas, res, res_num

def _anm_calc(pdb, calphas, res_info, out_cc_dir: Path, out_prs_dir: Path, out_single_dir: Path,
              uniprot1, uniprot2, chain1, chain2):
    anm = ANM(pdb)
    anm.buildHessian(calphas, cutoff=15.0, gamma=1)
    anm.calcModes(n_modes=None)
    anm_sq = calcSqFlucts(anm)
    anm_stiffness = np.mean(calcMechStiff(anm, calphas), axis=1)

    mode_num = anm._n_modes
    if mode_num < 5:
        return "Residue number is less than 5"

    anm_top3_cc       = calcCrossCorr(anm[0:3])
    anm_5_per_cc      = calcCrossCorr(anm[0:math.ceil(mode_num*0.05)])
    anm_5_20_per_cc   = calcCrossCorr(anm[math.ceil(mode_num*0.05):math.ceil(mode_num*0.2)])
    anm_20_50_per_cc  = calcCrossCorr(anm[math.ceil(mode_num*0.2):math.ceil(mode_num*0.5)])
    anm_greate_60_cc  = calcCrossCorr(anm[math.ceil(mode_num*0.6):])

    anm_prs, anm_eff, anm_sens = calcPerturbResponse(anm)

    #  DF 
    cc_list = [
        (_make_matrix_to_df(anm_top3_cc,      res_info), "anm_top3_cc"),
        (_make_matrix_to_df(anm_5_per_cc,     res_info), "anm_5_per_cc"),
        (_make_matrix_to_df(anm_5_20_per_cc,  res_info), "anm_5_20_per_cc"),
        (_make_matrix_to_df(anm_20_50_per_cc, res_info), "anm_20_50_per_cc"),
        (_make_matrix_to_df(anm_greate_60_cc, res_info), "anm_greate_60_per_cc"),
    ]
    for df, suffix in cc_list:
        fname = f"{pdb}_{uniprot1}_{uniprot2}_{chain1}_{chain2}_{suffix}.csv"
        df.to_csv(out_cc_dir / fname, index=True)

    prs_df = _make_matrix_to_df(anm_prs, res_info)
    prs_df.to_csv(out_prs_dir / f"{pdb}_{uniprot1}_{uniprot2}_{chain1}_{chain2}_anm_prs.csv", index=True)

    df = pd.DataFrame(res_info.copy())
    df["ANM_sq"] = anm_sq
    df["ANM_effectiveness"] = anm_eff
    df["ANM_sensitivity"]   = anm_sens
    df["ANM_stiffness"]     = anm_stiffness
    df.to_csv(out_single_dir / f"{pdb}_{chain1}_{chain2}_singleAA_data.csv", index=False)

def run_anm(args):
    if ANM is None:
        raise RuntimeError("ProDy not installed or failed to import. Please install with: pip install prody")

    pdb_dir = Path(args.pdb_dir)
    ppi_file = Path(args.ppi_info)
    out_single = ensure_dir(args.out_singleAA)
    out_cc     = ensure_dir(args.out_cc)
    out_prs    = ensure_dir(args.out_prs)

    pdb_files = {p.name: p for p in pdb_dir.iterdir() if p.suffix.lower() in (".pdb", ".cif")}
    with open(ppi_file, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    for line in lines:
        try:
            uniprot1, uniprot2, pdb, chain1, chain2 = line.split("\t")
            pdb_file = pdb + ".pdb" if pdb + ".pdb" in pdb_files else pdb + ".cif"
            if pdb_file not in pdb_files:
                print(f"[ANM] Skip {pdb}: no PDB/CIF found")
                continue
            path = pdb_files[pdb_file]
            if path.suffix.lower() == ".cif":
                structure = parseMMCIF(str(path))
                calphas, res3, res_num = _anm_core_build(structure, f"not ion and name CA and protein and segname {chain1} {chain2}")
                cif_seg1 = calphas[0].getSegname()
                cif_seg2 = calphas[-1].getSegname()
                if chain1 != cif_seg1 and chain2 != cif_seg2:
                    chain1, chain2 = chain2, chain1
                    uniprot1, uniprot2 = uniprot2, uniprot1
            else:
                structure = parsePDB(str(path))
                calphas, res3, res_num = _anm_core_build(structure, f"not ion and name CA and protein and chain {chain1} {chain2}")
                pdb_chain1 = calphas[0].getChid()
                pdb_chain2 = calphas[-1].getChid()
                if chain1 != pdb_chain1 and chain2 != pdb_chain2:
                    chain1, chain2 = chain2, chain1
                    uniprot1, uniprot2 = uniprot2, uniprot1

            res_num["chain"] = [chain1]*len(structure.select(f"not ion and name CA and protein and chain {chain1}" if path.suffix.lower()==".pdb" else f"not ion and name CA and protein and segname {chain1}")) + \
                               [chain2]*len(structure.select(f"not ion and name CA and protein and chain {chain2}" if path.suffix.lower()==".pdb" else f"not ion and name CA and protein and segname {chain2}"))
            res_num["uniprot"] = [uniprot1]*len(structure.select(f"not ion and name CA and protein and chain {chain1}" if path.suffix.lower()==".pdb" else f"not ion and name CA and protein and segname {chain1}")) + \
                                 [uniprot2]*len(structure.select(f"not ion and name CA and protein and chain {chain2}" if path.suffix.lower()==".pdb" else f"not ion and name CA and protein and segname {chain2}"))
            # ： Res_Info
            # Note: additional select here is only to get sequence length. For complex cases, users should ensure consistent chain labeling.

            #  Res_Info
            #  res3 ，（ calphas.getResnums()）
            res_names = [three_to_one.get(r, r) for r in res3]
            res_num["res"] = res_names
            res_num["Res_Info"] = res_num.apply(lambda row: '_'.join([row["uniprot"], row["chain"], row["res"]]) + str(row[0]), axis=1)
            res_info = res_num["Res_Info"]

            _anm_calc(pdb, calphas, res_info, Path(out_cc), Path(out_prs), Path(out_single), uniprot1, uniprot2, chain1, chain2)
            print(f"[ANM] Done: {pdb}")
        except Exception as e:
            print(f"[ANM] Error at line: {line} -> {e}")

# ============================= 2) Splitting/export (cluster vs noncluster) =============================
#  Get_Cluster_data.py ，/。
# For simplicity, reusing main functions with CLI-provided arguments.

def _split_sites(site_str):
    if pd.isna(site_str) or not str(site_str).strip():
        return set()
    return {s.strip() for s in str(site_str).split(',') if s.strip()}

_site_pat = re.compile(r"^(?P<pfx>.+?)_(?P<chain>[^_]+)_(?P<aa>[A-Za-z])?(?P<num>\d+)$")
def _to_X_form(site: str) -> str:
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

def _normalize_sites(sites):
    return {_to_X_form(s) for s in sites if isinstance(s, str) and s.strip()}

def _read_csv_guess_sep(path):
    df = pd.read_csv(path, index_col=0)
    if df.shape[1] <= 1:
        df = pd.read_csv(path, sep="\t", index_col=0, engine="c")
    return df

def _read_cc_files(base_dir, pdb, bands):
    out = {}
    band_re = re.compile(r"_(%s)_cc\.csv$" % "|".join(map(re.escape, bands)))
    for fp in sorted(glob.glob(os.path.join(base_dir, f"{pdb}_*_cc.csv"))):
        m = band_re.search(os.path.basename(fp))
        if not m:
            continue
        band = m.group(1)
        if band in out:
            continue
        df = _read_csv_guess_sep(fp)
        df.index = df.index.astype(str)
        df.columns = df.columns.astype(str)
        out[band] = df
    return out

def _normalize_index_and_columns(M: pd.DataFrame) -> pd.DataFrame:
    M = M.copy()
    M.index = [_to_X_form(str(i)) for i in M.index]
    M.columns = [_to_X_form(str(c)) for c in M.columns]
    return M

def _safe_submatrix(M: pd.DataFrame, rows, cols):
    r = sorted(set(rows) & set(M.index))
    c = sorted(set(cols) & set(M.columns))
    if not r or not c:
        return None
    return M.loc[r, c]

def run_split(args):
    # 
    xlsx_path = Path(args.cluster_sites_xlsx)
    bands = ["top3", "5_per", "5_20_per", "20_50_per", "greate_60_per"]
    # 
    anm_singleAA = Path(args.anm_singleAA)
    anm_cc       = Path(args.anm_cc)
    gnm_singleAA = Path(args.gnm_singleAA)
    gnm_cc       = Path(args.gnm_cc)
    # 
    anm_root = Path(args.anm_out_root); ensure_dir(anm_root)
    gnm_root = Path(args.gnm_out_root); ensure_dir(gnm_root)
    # 
    anm_out = {
        "singleAA_cluster": ensure_dir(anm_root / "singleAA_cluster"),
        "singleAA_noncluster": ensure_dir(anm_root / "singleAA_noncluster"),
        "CC_cluster": ensure_dir(anm_root / "CC_cluster"),
        "CC_noncluster": ensure_dir(anm_root / "CC_noncluster"),
    }
    gnm_out = {
        "singleAA_cluster": ensure_dir(gnm_root / "singleAA_cluster"),
        "singleAA_noncluster": ensure_dir(gnm_root / "singleAA_noncluster"),
        "CC_cluster": ensure_dir(gnm_root / "CC_cluster"),
        "CC_noncluster": ensure_dir(gnm_root / "CC_noncluster"),
    }

    base = pd.read_excel(xlsx_path)
    clusters_by_pdb = {}
    all_sites_by_pdb_orig = {}
    all_sites_by_pdb_X = {}
    all_interface_by_pdb_X = {}

    for pdb, dfp in base.groupby("PDB", sort=False):
        lst = []
        union_sites_orig = set()
        union_sites_X    = set()
        union_if_X       = set()
        for _, row in dfp.iterrows():
            cluster_id = str(row["Cluster_New"])
            cluster_sites_orig = _split_sites(row["Site_New"])
            interface_sites_orig = _split_sites(row.get("Interface_Site", None))
            cluster_sites_X   = _normalize_sites(cluster_sites_orig)
            interface_sites_X = _normalize_sites(interface_sites_orig)
            if cluster_sites_orig:
                lst.append((cluster_id, cluster_sites_orig, cluster_sites_X, interface_sites_X))
                union_sites_orig |= cluster_sites_orig
                union_sites_X    |= cluster_sites_X
                union_if_X       |= interface_sites_X
        clusters_by_pdb[pdb]        = lst
        all_sites_by_pdb_orig[pdb]  = union_sites_orig
        all_sites_by_pdb_X[pdb]     = union_sites_X
        all_interface_by_pdb_X[pdb] = union_if_X

    def _process_one(pdb, single_dir, cc_dir, out_map, tag="ANM", keep_noncluster_cc=True):
        df_single = None
        # singleAA：noncluster/cluster
        # singleAA ：{pdb}_*_singleAA_data.csv
        files = glob.glob(os.path.join(single_dir, f"{pdb}_*_singleAA_data.csv"))
        if files:
            df_single = pd.read_csv(files[0])
            if "Res_Info" in df_single.columns:
                # noncluster
                mask_non = ~df_single["Res_Info"].isin(all_sites_by_pdb_orig[pdb])
                non_df = df_single.loc[mask_non]
                if not non_df.empty:
                    non_df.to_csv(out_map["singleAA_noncluster"] / f"{pdb}.csv", index=False)
                #  cluster
                for cluster_id, cset_orig, _cset_X, _ifset_X in clusters_by_pdb[pdb]:
                    if not cset_orig:
                        continue
                    sub = df_single.loc[df_single["Res_Info"].isin(cset_orig)]
                    if not sub.empty:
                        sub.to_csv(out_map["singleAA_cluster"] / f"{cluster_id}.csv", index=False)

        # CC： band 
        cc_dict = _read_cc_files(cc_dir, pdb, bands)
        if not cc_dict:
            return f"[{tag}] {pdb} done (no CC)"
        for b in list(cc_dict.keys()):
            cc_dict[b] = _normalize_index_and_columns(cc_dict[b])

        if keep_noncluster_cc:
            for band, M in cc_dict.items():
                non_rows = sorted(set(M.index) - set(all_sites_by_pdb_X[pdb]))
                non_cols = sorted(set(all_interface_by_pdb_X[pdb]) & set(M.columns))
                sub = _safe_submatrix(M, non_rows, non_cols)
                if sub is not None and not sub.empty:
                    sub.to_csv(out_map["CC_noncluster"] / f"{pdb}_{band}_noncluster_vs_interface_cc.csv")

        for cluster_id, _cset_orig, cset_X, ifset_X in clusters_by_pdb[pdb]:
            if not cset_X or not ifset_X:
                continue
            for band, M in cc_dict.items():
                sub = _safe_submatrix(M, cset_X, ifset_X)
                if sub is not None and not sub.empty:
                    sub.to_csv(out_map["CC_cluster"] / f"{cluster_id}_{band}_cluster_vs_interface_cc.csv")

        return f"[{tag}] {pdb} done"

    pdbs = list(clusters_by_pdb.keys())
    # ANM
    for pdb in pdbs:
        print(_process_one(pdb, str(anm_singleAA), str(anm_cc), anm_out, "ANM"))
    # GNM
    for pdb in pdbs:
        print(_process_one(pdb, str(gnm_singleAA), str(gnm_cc), gnm_out, "GNM"))

# ============================= 3) Classification (CC/SingleAA) =============================
# 3.1 CC （ Classify_files_CC.py）
def bucket_by_name(fname: str):
    BUCKET_PATTERNS = [
        ("anm_top3"          , re.compile(r"_top3_")),
        ("anm_5_per"         , re.compile(r"_5_per_")),
        ("anm_5_20_per"      , re.compile(r"_5_20_per_")),
        ("anm_20_50_per"     , re.compile(r"_20_50_per_")),
        ("anm_greate_60_per" , re.compile(r"_greate_60_per_")),
    ]
    for bucket, pat in BUCKET_PATTERNS:
        if pat.search(fname):
            return bucket
    return None

CLUSTER_ID_RE = re.compile(r"^(Cluster\.\d+(?:\.\d+)?)_")

def _load_pos_map(xlsx_path: Path) -> dict:
    df = pd.read_excel(xlsx_path, usecols=["Cluster_New", "Pos_Type"]).dropna(subset=["Cluster_New"])
    df["Cluster_New"] = df["Cluster_New"].astype(str).str.strip()
    df["Pos_Type"]    = df["Pos_Type"].astype(str).str.strip()
    return dict(zip(df["Cluster_New"], df["Pos_Type"]))

def run_classify_cc(args):
    CLUSTER_IN_DIR    = Path(args.cc_cluster)
    NONCLUSTER_IN_DIR = Path(args.cc_noncluster)
    CLASS_ROOT        = Path(args.class_root)
    META_XLSX         = Path(args.meta_xlsx)

    POS_DIR = {"Interface":"Interface","Non-Interface":"Non_Interface","Cross-Interface":"Cross_Interface"}

    # 
    for bucket in ["anm_top3","anm_5_per","anm_5_20_per","anm_20_50_per","anm_greate_60_per"]:
        for sub in ["Interface","Non_Interface","Cross_Interface","noncluster"]:
            ensure_dir(CLASS_ROOT / bucket / sub)

    pos_map = _load_pos_map(META_XLSX)

    n = 0
    for fp in sorted(CLUSTER_IN_DIR.glob("*.csv")):
        bucket = bucket_by_name(fp.name)
        if not bucket: continue
        m = CLUSTER_ID_RE.match(fp.name)
        if not m: continue
        cid = m.group(1)
        pos = pos_map.get(cid)
        subdir = {"Interface":"Interface","Non-Interface":"Non_Interface","Cross-Interface":"Cross_Interface"}.get(pos)
        if not subdir: continue
        dst = CLASS_ROOT / bucket / subdir / fp.name
        if not dst.exists():
            dst.write_bytes(fp.read_bytes()); n += 1
    print(f"[classify CC cluster] copied {n} files")

    n = 0
    for fp in sorted(NONCLUSTER_IN_DIR.glob("*.csv")):
        bucket = bucket_by_name(fp.name)
        if not bucket: continue
        dst = CLASS_ROOT / bucket / "noncluster" / fp.name
        if not dst.exists():
            dst.write_bytes(fp.read_bytes()); n += 1
    print(f"[classify CC noncluster] copied {n} files")
    print(f"[DONE] output root: {CLASS_ROOT}")

# 3.2 singleAA （ Classify_files_singleAA.py）
def run_classify_singleaa(args):
    SINGLEAA_CLUSTER_DIR = Path(args.saa_cluster)
    SINGLEAA_NONCLU_DIR  = Path(args.saa_noncluster)
    CLASS_ROOT           = Path(args.class_root)
    META_XLSX            = Path(args.meta_xlsx)

    POS_DIR = {"Interface":"Interface","Non-Interface":"Non_Interface","Cross-Interface":"Cross_Interface"}
    # 
    for sub in ["Interface","Non_Interface","Cross_Interface","noncluster"]:
        ensure_dir(CLASS_ROOT / "singleAA" / sub)

    pos_map = _load_pos_map(META_XLSX)
    CLUSTER_ID_RE2 = re.compile(r"(Cluster\.\d+(?:\.\d+)?)")

    # cluster
    n = 0
    for ext in ("*.xlsx","*.csv"):
        for fp in sorted(SINGLEAA_CLUSTER_DIR.glob(ext)):
            m = CLUSTER_ID_RE2.search(fp.name)
            if not m: continue
            cid = m.group(1)
            pos = pos_map.get(cid)
            subdir = POS_DIR.get(pos)
            if not subdir: continue
            dst = CLASS_ROOT / "singleAA" / subdir / fp.name
            if not dst.exists():
                dst.write_bytes(fp.read_bytes()); n += 1
    print(f"[classify singleAA cluster] copied {n} files")

    # noncluster
    n = 0
    for ext in ("*.xlsx","*.csv"):
        for fp in sorted(SINGLEAA_NONCLU_DIR.glob(ext)):
            dst = CLASS_ROOT / "singleAA" / "noncluster" / fp.name
            if not dst.exists():
                dst.write_bytes(fp.read_bytes()); n += 1
    print(f"[classify singleAA noncluster] copied {n} files")
    print(f"[DONE] output root: {CLASS_ROOT/'singleAA'}")

# ============================= 4) Aggregation (CC / singleAA) =============================
# 4.1 CC （Get_CC_cluster_data.py ：）
def _read_matrix(csv_path: Path) -> np.ndarray:
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

def _append_values_to_csv(values: np.ndarray, out_csv: Path, chunk_n=1_000_000, header_written=False):
    if out_csv.exists() and not header_written:
        out_csv.unlink()
    start, n = 0, values.size
    while start < n:
        end = min(start + chunk_n, n)
        pd.DataFrame({"value": values[start:end]}).to_csv(
            out_csv, mode="a", index=False, header=not header_written
        )
        header_written = True
        start = end

def run_aggregate_cc(args):
    CLASSIFY_ROOT = Path(args.classify_root)
    OUT_ROOT      = Path(args.out_root)
    ABS_VAL       = bool(args.abs_val)
    CHUNK_N       = int(args.chunk_n)

    BUCKETS = ["anm_top3", "anm_5_per", "anm_5_20_per", "anm_20_50_per", "anm_greate_60_per"]
    SUBFOLDERS = ["Interface", "Non_Interface", "Cross_Interface", "noncluster"]

    ensure_dir(OUT_ROOT)

    for bucket in BUCKETS:
        for sub in SUBFOLDERS:
            in_dir = CLASSIFY_ROOT / bucket / sub
            out_csv = OUT_ROOT / f"{bucket}_{sub.lower()}.csv"
            files = sorted(in_dir.glob("*.csv"))
            if not files:
                print(f"[] ：{in_dir}")
                continue
            if out_csv.exists():
                out_csv.unlink()
            header_written = False
            for fp in files:
                try:
                    mat = _read_matrix(fp)
                    vals = mat.ravel()
                    vals = vals[~np.isnan(vals)]
                    if ABS_VAL:
                        vals = np.abs(vals)
                    if vals.size == 0:
                        continue
                    _append_values_to_csv(vals, out_csv, CHUNK_N, header_written)
                    header_written = True
                except Exception as e:
                    print(f"[] {fp.name} ：{e}")
            print(f"[] {bucket}/{sub} -> {out_csv}")

# 4.2 singleAA （Get_singleAA_cluster_data.py ）
def run_aggregate_singleaa(args):
    BASE_DIR = Path(args.saa_classify)
    OUT_DIR  = ensure_dir(BASE_DIR / "aggregated")
    TARGET_COLS = ["ANM_sq", "ANM_effectiveness", "ANM_sensitivity", "ANM_stiffness"]
    CHUNK_N = int(args.chunk_n)

    def iter_files(root: Path):
        for pat in ("*.csv", "*.xlsx"):
            for fp in sorted(root.glob(pat)):
                if fp.is_file():
                    yield fp

    def read_needed(fp: Path, cols: list[str]) -> pd.DataFrame:
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

    def append_series_csv(values: np.ndarray, out_csv: Path, colname: str, header_written: dict):
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        start, n = 0, values.size
        while start < n:
            end = min(start + CHUNK_N, n)
            pd.DataFrame({colname: values[start:end]}).to_csv(
                out_csv, mode="a", index=False, header=(not header_written.get(out_csv, False))
            )
            header_written[out_csv] = True
            start = end

    subdirs = [p for p in sorted(BASE_DIR.iterdir()) if p.is_dir() and p.name != "aggregated"]
    if not subdirs:
        print(f"[] ：{BASE_DIR}")
        return

    for in_dir in subdirs:
        tag = in_dir.name
        out_map = {
            "ANM_sq":            OUT_DIR / f"{tag}_ANM_sq.csv",
            "ANM_effectiveness": OUT_DIR / f"{tag}_ANM_effectiveness.csv",
            "ANM_sensitivity":   OUT_DIR / f"{tag}_ANM_sensitivity.csv",
            "ANM_stiffness":     OUT_DIR / f"{tag}_ANM_stiffness.csv",
        }
        for p in out_map.values():
            if p.exists():
                p.unlink()
        header_written = {}
        files = list(iter_files(in_dir))
        if not files:
            print(f"[] ：{in_dir}"); continue

        for fp in files:
            try:
                df = read_needed(fp, TARGET_COLS)
                if df.empty: continue
                for col, out_csv in out_map.items():
                    vals = df[col].to_numpy(dtype=float, copy=False)
                    vals = vals[~np.isnan(vals)]
                    if vals.size:
                        append_series_csv(vals, out_csv, col, header_written)
            except Exception as e:
                print(f"[]  {fp.name}: {e}")

        for col, out_csv in out_map.items():
            if out_csv.exists():
                n = max(0, sum(1 for _ in open(out_csv, "r")) - 1)
                print(f"[] {tag}: {col} -> {out_csv}   {n} ")

# ============================= 5) Final table generation (merge statistics) =============================
#  Get_Final_data.py （ CSV/Excel；）
def run_final(args):
    base_table_path = args.base_xlsx
    out_path        = args.out_xlsx

    # “”
    anm_root = Path(args.anm_root)
    gnm_root = Path(args.gnm_root)

    paths = {
        "anm_singleAA_cluster_dir":   str(anm_root / "singleAA_cluster"),
        "anm_singleAA_noncluster_dir":str(anm_root / "singleAA_noncluster"),
        "anm_CC_cluster_dir":         str(anm_root / "CC_cluster"),
        "anm_CC_noncluster_dir":      str(anm_root / "CC_noncluster"),
        "gnm_singleAA_cluster_dir":   str(gnm_root / "singleAA_cluster"),
        "gnm_singleAA_noncluster_dir":str(gnm_root / "singleAA_noncluster"),
        "gnm_CC_cluster_dir":         str(gnm_root / "CC_cluster"),
        "gnm_CC_noncluster_dir":      str(gnm_root / "CC_noncluster"),
    }

    def read_sheet_any(path):
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
            print(f"[WARN]  {path}: {e}")
            return None

    def pick_col(df, candidates):
        if df is None: return None
        lower_map = {c.lower(): c for c in df.columns}
        for cand in candidates:
            if cand.lower() in lower_map:
                return lower_map[cand.lower()]
        return None

    def max_singleAA_metrics(file_path, model="ANM"):
        df = read_sheet_any(file_path)
        if df is None or df.empty:
            return (np.nan, np.nan, np.nan, np.nan) if model=="ANM" else (np.nan, np.nan, np.nan)
        if model.upper() == "ANM":
            sq_col        = pick_col(df, ["ANM_sq","sq"])
            sensor_col    = pick_col(df, ["ANM_sensitivity","sensitivity","sensor"])
            effector_col  = pick_col(df, ["ANM_effectiveness","effectiveness","effector"])
            stiffness_col = pick_col(df, ["ANM_stiffness","stiffness"])
        else:
            sq_col        = pick_col(df, ["GNM_sq","sq"])
            sensor_col    = pick_col(df, ["GNM_sensitivity","sensitivity","sensor"])
            effector_col  = pick_col(df, ["GNM_effectiveness","effectiveness","effector"])
            stiffness_col = None
        def safe_max(col):
            if col is None: return np.nan
            s = pd.to_numeric(df[col], errors="coerce")
            return np.nanmax(s.values) if np.isfinite(s).any() else np.nan
        if model.upper() == "ANM":
            return safe_max(sq_col), safe_max(sensor_col), safe_max(effector_col), safe_max(stiffness_col)
        else:
            return safe_max(sq_col), safe_max(sensor_col), safe_max(effector_col)

    def max_cc(path):
        if path is None or not os.path.exists(path):
            return np.nan
        df = None
        for enc in ("utf-8", "gbk", "latin1"):
            try:
                df = pd.read_csv(path, encoding=enc, engine="c", low_memory=False); break
            except Exception:
                try:
                    df = pd.read_csv(path, encoding=enc, engine="python", low_memory=False, on_bad_lines="skip"); break
                except Exception:
                    continue
        if df is None or df.empty:
            return np.nan
        first_col = df.columns[0]
        if df[first_col].dtype == object:
            t = pd.to_numeric(df[first_col], errors="coerce")
            if t.isna().mean() > 0.5:
                df = df.set_index(first_col)
        if df.index.has_duplicates:
            df = df.loc[~df.index.duplicated(keep="first")]
        if pd.Index(df.columns).has_duplicates:
            df = df.loc[:, ~pd.Index(df.columns).duplicated(keep="first")]
        idx = df.index.astype(str); cols = df.columns.astype(str)
        common = idx.intersection(cols)
        if len(common) > 0:
            df = df.loc[common, common]
        else:
            df = df.apply(pd.to_numeric, errors="coerce")
            if df.size == 0: return np.nan
            r, c = df.shape; n = min(r, c); df = df.iloc[:n, :n]
        df = df.apply(pd.to_numeric, errors="coerce")
        arr = df.to_numpy(dtype=float, copy=False)
        r, c = arr.shape; n = min(r, c)
        if n == 0: return np.nan
        if r != n or c != n: arr = arr[:n, :n]
        iu = np.triu_indices(n, k=1)
        sub = arr[iu]; sub = sub[np.isfinite(sub)]
        return float(sub.max()) if sub.size else np.nan

    base = pd.read_excel(base_table_path)
    new_cols = [
        "anm_cluster_sq","anm_cluster_sensor","anm_cluster_effector","anm_cluster_stiffness","anm_cluster_top3cc","anm_cluster_5percc",
        "gnm_cluster_sq","gnm_cluster_sensor","gnm_cluster_effector","gnm_cluster_top3cc","gnm_cluster_5percc",
        "anm_noncluster_sq","anm_noncluster_sensor","anm_noncluster_effector","anm_noncluster_stiffness","anm_noncluster_top3cc","anm_noncluster_5percc",
        "gnm_noncluster_sq","gnm_noncluster_sensor","gnm_noncluster_effector","gnm_noncluster_top3cc","gnm_noncluster_5percc"
    ]
    for c in new_cols: base[c] = np.nan

    for i, row in base.iterrows():
        pdb = str(row["PDB"]).strip()
        cluster = str(row["Cluster_New"]).strip()

        # ANM singleAA cluster/noncluster
        f = os.path.join(paths["anm_singleAA_cluster_dir"], f"{cluster}.csv")
        sq, sensor, eff, stiff = max_singleAA_metrics(f, "ANM")
        base.loc[i, ["anm_cluster_sq","anm_cluster_sensor","anm_cluster_effector","anm_cluster_stiffness"]] = [sq, sensor, eff, stiff]

        f = os.path.join(paths["anm_singleAA_noncluster_dir"], f"{pdb}.csv")
        sq, sensor, eff, stiff = max_singleAA_metrics(f, "ANM")
        base.loc[i, ["anm_noncluster_sq","anm_noncluster_sensor","anm_noncluster_effector","anm_noncluster_stiffness"]] = [sq, sensor, eff, stiff]

        # ANM CC
        f_top3 = os.path.join(paths["anm_CC_cluster_dir"], f"{cluster}_top3_cluster_vs_interface_cc.csv")
        f_5per = os.path.join(paths["anm_CC_cluster_dir"], f"{cluster}_5_per_cluster_vs_interface_cc.csv")
        if not os.path.exists(f_top3): f_top3 = os.path.join(paths["anm_CC_cluster_dir"], f"{cluster}_top3_cluster_cc.csv")
        if not os.path.exists(f_5per): f_5per = os.path.join(paths["anm_CC_cluster_dir"], f"{cluster}_5_per_cluster_cc.csv")
        base.loc[i, "anm_cluster_top3cc"] = max_cc(f_top3)
        base.loc[i, "anm_cluster_5percc"] = max_cc(f_5per)

        f_top3 = os.path.join(paths["anm_CC_noncluster_dir"], f"{pdb}_top3_noncluster_vs_interface_cc.csv")
        f_5per = os.path.join(paths["anm_CC_noncluster_dir"], f"{pdb}_5_per_noncluster_vs_interface_cc.csv")
        if not os.path.exists(f_top3): f_top3 = os.path.join(paths["anm_CC_noncluster_dir"], f"{pdb}_top3_noncluster_cc.csv")
        if not os.path.exists(f_5per): f_5per = os.path.join(paths["anm_CC_noncluster_dir"], f"{pdb}_5_per_noncluster_cc.csv")
        base.loc[i, "anm_noncluster_top3cc"] = max_cc(f_top3)
        base.loc[i, "anm_noncluster_5percc"] = max_cc(f_5per)

        # GNM singleAA
        f = os.path.join(paths["gnm_singleAA_cluster_dir"], f"{cluster}.csv")
        sq, sensor, eff = max_singleAA_metrics(f, "GNM")
        base.loc[i, ["gnm_cluster_sq","gnm_cluster_sensor","gnm_cluster_effector"]] = [sq, sensor, eff]
        f = os.path.join(paths["gnm_singleAA_noncluster_dir"], f"{pdb}.csv")
        sq, sensor, eff = max_singleAA_metrics(f, "GNM")
        base.loc[i, ["gnm_noncluster_sq","gnm_noncluster_sensor","gnm_noncluster_effector"]] = [sq, sensor, eff]

        # GNM CC
        f_top3 = os.path.join(paths["gnm_CC_cluster_dir"], f"{cluster}_top3_cluster_vs_interface_cc.csv")
        f_5per = os.path.join(paths["gnm_CC_cluster_dir"], f"{cluster}_5_per_cluster_vs_interface_cc.csv")
        if not os.path.exists(f_top3): f_top3 = os.path.join(paths["gnm_CC_cluster_dir"], f"{cluster}_top3_cluster_cc.csv")
        if not os.path.exists(f_5per): f_5per = os.path.join(paths["gnm_CC_cluster_dir"], f"{cluster}_5_per_cluster_cc.csv")
        base.loc[i, "gnm_cluster_top3cc"] = max_cc(f_top3)
        base.loc[i, "gnm_cluster_5percc"] = max_cc(f_5per)

        f_top3 = os.path.join(paths["gnm_CC_noncluster_dir"], f"{pdb}_top3_noncluster_vs_interface_cc.csv")
        f_5per = os.path.join(paths["gnm_CC_noncluster_dir"], f"{pdb}_5_per_noncluster_vs_interface_cc.csv")
        if not os.path.exists(f_top3): f_top3 = os.path.join(paths["gnm_CC_noncluster_dir"], f"{pdb}_top3_noncluster_cc.csv")
        if not os.path.exists(f_5per): f_5per = os.path.join(paths["gnm_CC_noncluster_dir"], f"{pdb}_5_per_noncluster_cc.csv")
        base.loc[i, "gnm_noncluster_top3cc"] = max_cc(f_top3)
        base.loc[i, "gnm_noncluster_5percc"] = max_cc(f_5per)

    base.to_excel(out_path, index=False)
    print("[FINAL] Done ->", out_path)

# ============================= 6) All-in-one pipeline =============================
def run_all(args):
    # May only provide partial parameters; each subcommand validates or defaults internally.
    if args.do_anm:
        print("== Step 1: ANM ==")
        run_anm(args)
    if args.do_split:
        print("== Step 2: Split ==")
        run_split(args)
    if args.do_classify_cc:
        print("== Step 3a: Classify CC ==")
        run_classify_cc(args)
    if args.do_classify_saa:
        print("== Step 3b: Classify singleAA ==")
        run_classify_singleaa(args)
    if args.do_agg_cc:
        print("== Step 4a: Aggregate CC ==")
        run_aggregate_cc(args)
    if args.do_agg_saa:
        print("== Step 4b: Aggregate singleAA ==")
        run_aggregate_singleaa(args)
    if args.do_final:
        print("== Step 5: Final Table ==")
        run_final(args)

# ============================= CLI =============================
def build_parser():
    p = argparse.ArgumentParser(description="PPI dynamics analysis integrated script")
    sub = p.add_subparsers(dest="cmd", required=True)

    # anm
    a = sub.add_parser("anm", help="Run ANM calculation (PDB/CIF -> singleAA/CC/PRS)")
    a.add_argument("--pdb-dir", required=True)
    a.add_argument("--ppi-info", required=True, help=" (uniprot1, uniprot2, pdb, chain1, chain2)  TSV ")
    a.add_argument("--out-singleAA", required=True)
    a.add_argument("--out-cc", required=True)
    a.add_argument("--out-prs", required=True)
    a.set_defaults(func=run_anm)

    # split
    s = sub.add_parser("split", help="Split singleAA/CC into cluster/noncluster based on cluster_sites table (supports ANM and GNM)")
    s.add_argument("--cluster-sites-xlsx", required=True)
    s.add_argument("--anm-singleAA", required=True)
    s.add_argument("--anm-cc", required=True)
    s.add_argument("--gnm-singleAA", required=True)
    s.add_argument("--gnm-cc", required=True)
    s.add_argument("--anm-out-root", required=True)
    s.add_argument("--gnm-out-root", required=True)
    s.set_defaults(func=run_split)

    # classify-cc
    c1 = sub.add_parser("classify-cc", help="Classify CC files into buckets and Pos_Type categories")
    c1.add_argument("--cc-cluster", required=True)
    c1.add_argument("--cc-noncluster", required=True)
    c1.add_argument("--class-root", required=True)
    c1.add_argument("--meta-xlsx", required=True)
    c1.set_defaults(func=run_classify_cc)

    # classify-singleaa
    c2 = sub.add_parser("classify-singleaa", help="Classify singleAA files into Pos_Type categories")
    c2.add_argument("--saa-cluster", required=True)
    c2.add_argument("--saa-noncluster", required=True)
    c2.add_argument("--class-root", required=True)
    c2.add_argument("--meta-xlsx", required=True)
    c2.set_defaults(func=run_classify_singleaa)

    # aggregate-cc
    g1 = sub.add_parser("aggregate-cc", help="Flatten CC matrix values into single-column CSV (one per directory)")
    g1.add_argument("--classify-root", required=True)
    g1.add_argument("--out-root", required=True)
    g1.add_argument("--abs-val", action="store_true", help="")
    g1.add_argument("--chunk-n", type=int, default=1_000_000, help="")
    g1.set_defaults(func=run_aggregate_cc)

    # aggregate-singleaa
    g2 = sub.add_parser("aggregate-singleaa", help="Aggregate four singleAA metrics into four single-column CSVs (per directory)")
    g2.add_argument("--saa-classify", required=True, help="Classify/singleAA ")
    g2.add_argument("--chunk-n", type=int, default=1_000_000)
    g2.set_defaults(func=run_aggregate_singleaa)

    # final
    f = sub.add_parser("final", help="Merge aggregated metrics into base table to produce final results")
    f.add_argument("--base-xlsx", required=True)
    f.add_argument("--out-xlsx", required=True)
    f.add_argument("--anm-root", required=True, help="ANM  cluster/noncluster ")
    f.add_argument("--gnm-root", required=True, help="GNM  cluster/noncluster ")
    f.set_defaults(func=run_final)

    # all (pipeline)
    allp = sub.add_parser("all", help="Execute multiple steps in order based on flags")
    # ，；
    # anm
    allp.add_argument("--pdb-dir")
    allp.add_argument("--ppi-info")
    allp.add_argument("--out-singleAA")
    allp.add_argument("--out-cc")
    allp.add_argument("--out-prs")
    # split
    allp.add_argument("--cluster-sites-xlsx")
    allp.add_argument("--anm-singleAA")
    allp.add_argument("--anm-cc")
    allp.add_argument("--gnm-singleAA")
    allp.add_argument("--gnm-cc")
    allp.add_argument("--anm-out-root")
    allp.add_argument("--gnm-out-root")
    # classify
    allp.add_argument("--cc-cluster")
    allp.add_argument("--cc-noncluster")
    allp.add_argument("--class-root")
    allp.add_argument("--meta-xlsx")
    allp.add_argument("--saa-cluster")
    allp.add_argument("--saa-noncluster")
    # aggregate
    allp.add_argument("--classify-root")
    allp.add_argument("--out-root")
    allp.add_argument("--abs-val", action="store_true")
    allp.add_argument("--chunk-n", type=int, default=1_000_000)
    allp.add_argument("--saa-classify")
    # final
    allp.add_argument("--base-xlsx")
    allp.add_argument("--out-xlsx")
    allp.add_argument("--anm-root")
    allp.add_argument("--gnm-root")
    # （）
    allp.add_argument("--do-anm", action="store_true")
    allp.add_argument("--do-split", action="store_true")
    allp.add_argument("--do-classify-cc", action="store_true")
    allp.add_argument("--do-classify-saa", action="store_true")
    allp.add_argument("--do-agg-cc", action="store_true")
    allp.add_argument("--do-agg-saa", action="store_true")
    allp.add_argument("--do-final", action="store_true")
    allp.set_defaults(func=run_all)

    return p

def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
