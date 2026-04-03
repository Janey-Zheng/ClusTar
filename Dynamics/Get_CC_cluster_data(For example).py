# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import argparse

# Buckets (main folders)
BUCKETS = ["anm_top3", "anm_5_per", "anm_5_20_per", "anm_20_50_per", "anm_greate_60_per"]
# Subfolders
SUBFOLDERS = ["Interface", "Non_Interface", "Cross_Interface", "noncluster"]

def read_matrix(csv_path: Path) -> np.ndarray:
    """Read matrix from CSV file and return as numpy array"""
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
    """Extract values from matrices and save into a single-column CSV"""
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
            print(f"[Warning] Failed to process {fp.name}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Extract matrix values from CSV files into single-column CSVs.")
    parser.add_argument("--in_dir",  type=str, required=True, help="Input root directory (Classify)")
    parser.add_argument("--out_dir", type=str, required=True, help="Output root directory (Classify_out)")
    parser.add_argument("--abs",     action="store_true", help="Take absolute values")
    parser.add_argument("--chunk",   type=int, default=1_000_000, help="Chunk size for writing (default: 1,000,000)")

    args = parser.parse_args()

    classify_root = Path(args.in_dir)
    out_root      = Path(args.out_dir)
    abs_val       = args.abs
    chunk_n       = args.chunk

    out_root.mkdir(parents=True, exist_ok=True)

    for bucket in BUCKETS:
        for sub in SUBFOLDERS:
            in_dir = classify_root / bucket / sub
            out_csv = out_root / f"{bucket}_{sub.lower()}.csv"

            files = sorted(in_dir.glob("*.csv"))
            if not files:
                print(f"[Info] Directory is empty: {in_dir}")
                continue

            write_all_values(files, out_csv, abs_val=abs_val, chunk_n=chunk_n)
            print(f"[Done] {bucket}/{sub} -> {out_csv}")

if __name__ == "__main__":
    main()
