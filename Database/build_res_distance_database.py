# -*- coding: gbk -*-

import pandas as pd
import psycopg2
from Bio.PDB import PDBParser, MMCIFParser, PDBExceptions
from Bio.SeqUtils import IUPACData
import numpy as np
import glob
import os
import argparse
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

def calculate_min_distance(res1, res2):
    coords1 = np.array([atom.get_coord() for atom in res1.get_atoms()])
    coords2 = np.array([atom.get_coord() for atom in res2.get_atoms()])
    if coords1.size == 0 or coords2.size == 0:
        return float('inf')
    distances = np.linalg.norm(coords1[:, None, :] - coords2[None, :, :], axis=-1)
    return distances.min()

def create_table_name(pdb_id, uniprot1, uniprot2, chain1, chain2):
    return f'"{pdb_id}_{uniprot1}_{uniprot2}_{chain1}_{chain2}"'

def process_pdb(pdb_path, chain_pairs, conn_info, processed_files_log):
    pdb_id = os.path.basename(pdb_path).split('.')[0]
    parser = PDBParser(QUIET=True) if pdb_path.endswith(".pdb") else MMCIFParser(QUIET=True)

    try:
        structure = parser.get_structure(pdb_id, pdb_path)
    except (PDBExceptions.PDBConstructionException, KeyError) as e:
        print(f"Error parsing {pdb_id}: {e}")
        return

    conn = psycopg2.connect(conn_info)
    cursor = conn.cursor()

    if pdb_id not in chain_pairs:
        conn.close()
        return

    for uniprot1, uniprot2, chain1, chain2 in chain_pairs[pdb_id]:
        model = next(iter(structure))
        ch1 = model[chain1] if chain1 in model else None
        ch2 = model[chain2] if chain2 in model else None

        if ch1 is None or ch2 is None:
            print(f"Skipping {pdb_id}: Chain {chain1} or {chain2} not found")
            continue

        table_name = create_table_name(pdb_id, uniprot1, uniprot2, chain1, chain2)
        cursor.execute(f'''
            CREATE TABLE IF NOT EXISTS {table_name} (
                id SERIAL PRIMARY KEY,
                Chain1 TEXT, Residue1_Name TEXT, Residue1_Num INTEGER,
                Chain2 TEXT, Residue2_Name TEXT, Residue2_Num INTEGER,
                Distance FLOAT
            )
        ''')

        data_to_insert = []
        for res1 in ch1:
            for res2 in ch2:
                if res1.id[0] == " " and res2.id[0] == " ":
                    if res1.id[1] == res2.id[1] and chain1 == chain2:
                        continue
                    distance = round(calculate_min_distance(res1, res2), 2)
                    res1name = IUPACData.protein_letters_3to1.get(res1.get_resname().capitalize(), 'X')
                    res2name = IUPACData.protein_letters_3to1.get(res2.get_resname().capitalize(), 'X')
                    data_to_insert.append((chain1, res1name, res1.id[1],
                                           chain2, res2name, res2.id[1], float(distance)))
        for res3 in ch1:
            for res4 in ch1:
                if res3.id[0] == " " and res4.id[0] == " ":
                    if res3.id[1] == res4.id[1]:
                        continue
                    distance_chain1 = round(calculate_min_distance(res3, res4), 2)
                    res3name = IUPACData.protein_letters_3to1.get(res3.get_resname().capitalize(), 'X')
                    res4name = IUPACData.protein_letters_3to1.get(res4.get_resname().capitalize(), 'X')
                    data_to_insert.append((chain1, res3name, res3.id[1],
                                           chain1, res4name, res4.id[1], float(distance_chain1)))
        for res5 in ch2:
            for res6 in ch2:
                if res5.id[0] == " " and res6.id[0] == " ":
                    if res5.id[1] == res6.id[1]:
                        continue
                    distance_chain2 = round(calculate_min_distance(res5, res6), 2)
                    res5name = IUPACData.protein_letters_3to1.get(res5.get_resname().capitalize(), 'X')
                    res6name = IUPACData.protein_letters_3to1.get(res6.get_resname().capitalize(), 'X')
                    data_to_insert.append((chain2, res5name, res5.id[1],
                                           chain2, res6name, res6.id[1], float(distance_chain2)))

        if data_to_insert:
            cursor.executemany(f'''
                INSERT INTO {table_name} (Chain1, Residue1_Name, Residue1_Num, Chain2, Residue2_Name, Residue2_Num, Distance)
                VALUES (%s, %s, %s, %s, %s, %s, %s)
            ''', data_to_insert)

    conn.commit()
    cursor.close()
    conn.close()

    with open(processed_files_log, 'a') as log_file:
        log_file.write(f"{pdb_id}\n")

def get_chain_pairs(csv_path):
    df = pd.read_csv(csv_path, sep="\t")
    chain_pairs = {}
    for _, row in df.iterrows():
        pdb_id = row['PDB_ID']
        chain_pairs.setdefault(pdb_id, []).append((row['Uniprot1'], row['Uniprot2'], row['Chain1'], row['Chain2']))
    return chain_pairs

def get_or_create_processed_log(processed_files_log):
    if not os.path.exists(processed_files_log):
        open(processed_files_log, 'w').close()
    with open(processed_files_log, 'r') as log_file:
        processed_files = set(log_file.read().splitlines())
    return processed_files

def get_pdb_ids_from_csv(csv_path):
    df = pd.read_csv(csv_path, sep="\t")
    return set(df['PDB_ID'].unique())

def filter_pdb_files(input_directory, pdb_ids):
    pdb_files = glob.glob(os.path.join(input_directory, '*.pdb')) + glob.glob(os.path.join(input_directory, '*.cif'))
    return [p for p in pdb_files if os.path.basename(p).split('.')[0] in pdb_ids]

def main():
    parser = argparse.ArgumentParser(description="Process PDB/MMCIF files and store residue distances into PostgreSQL.")
    parser.add_argument('--input_dir', required=True, help='Directory containing PDB or CIF files')
    parser.add_argument('--csv', required=True, help='Path to the TSV file containing chain pair information')
    parser.add_argument('--db', required=True, help='PostgreSQL connection string (e.g. "dbname=... user=... password=... host=... port=...")')
    parser.add_argument('--log', default='processed_files.log', help='Path to the processed files log (default: processed_files.log)')
    parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers (default: 4)')
    
    args = parser.parse_args()

    chain_pairs = get_chain_pairs(args.csv)
    pdb_ids = get_pdb_ids_from_csv(args.csv)
    filtered_pdbs = filter_pdb_files(args.input_dir, pdb_ids)
    processed_files = get_or_create_processed_log(args.log)

    to_process = [p for p in filtered_pdbs if os.path.basename(p).split('.')[0] not in processed_files]

    if to_process:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(process_pdb, pdb_file, chain_pairs, args.db, args.log): pdb_file
                for pdb_file in to_process
            }
            for future in tqdm(as_completed(futures), total=len(to_process), desc="Processing PDBs"):
                try:
                    future.result()
                except Exception as e:
                    print(f"Error in processing: {e}")
    else:
        print("No PDB files to process.")

if __name__ == "__main__":
    main()
