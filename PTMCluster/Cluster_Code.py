# -*- coding: utf-8 -*-
import cProfile, pstats
import io
import os
import pandas as pd
# from prody import *
import numpy as np
import random
import csv
from Bio.PDB import PDBParser, PPBuilder ,PDBList
import math
import networkx as nx
from scipy.sparse.csgraph import floyd_warshall
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from Bio.PDB.Polypeptide import is_aa
from scipy.spatial import KDTree
from multiprocessing import Pool
from sklearn.cluster import AgglomerativeClustering
from collections import defaultdict
import psycopg2
from Bio.SeqUtils import IUPACData
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import shutil
from shutil import rmtree
from concurrent.futures import TimeoutError
from numpy import sqrt
from typing import Tuple, List


# 处理单个文件，提取所有站点对及其距离和p值
def process_file(df, db_conn_info_mon, db_conn_info_ppi):
    table = df["Database"].iat[0]
    class_val = df["Class"].iat[0]
    db_info   = db_conn_info_mon if class_val == 'Monomer' else db_conn_info_ppi
    conn = psycopg2.connect(db_info)
    cur  = conn.cursor()
    cur.execute(f'''
        SELECT residue1_num, residue2_num, chain1, chain2, distance
        FROM "{table}"
    ''')
    rows = cur.fetchall()  # List of (r1, r2, c1, c2, d)
    conn.close()
    dists = np.array([r[4] for r in rows], dtype=float)
    all_distances = np.round(np.sort(dists), 3)
    dist_dict = {}
    for r1, r2, c1, c2, d in rows:
        d3 = round(d, 3)
        dist_dict[(r1, c1, r2, c2)] = d3
        dist_dict[(r2, c2, r1, c1)] = d3
    for r1, r2, c1, c2, _ in rows:
        dist_dict[(r1, c1, r1, c1)] = 0.000

    site_pairs = []
    recs = df.to_dict("records")
    n = len(recs)
    for i in range(n):
        r1, c1 = recs[i]["Site_Num"], recs[i]["Chain"]
        for j in range(i+1, n):
            r2, c2 = recs[j]["Site_Num"], recs[j]["Chain"]
            d = dist_dict.get((r1, c1, r2, c2))
            if d is None:
                continue
            cnt  = np.searchsorted(all_distances, d, side="left")
            pval = cnt / all_distances.size
            site_pairs.append(((r1, c1), (r2, c2), d, pval))

    return site_pairs, dist_dict




        
# 定义筛选出显著性位点对的函数
def filter_significant(site_pairs: list) -> list:
    significant_pairs = []
    for (r1, c1), (r2, c2), d, pval in site_pairs:
        if pval >= 0.05:
            continue
        if c1 == c2:
            if d < 10 and abs(r1 - r2) > 20:
                significant_pairs.append(((r1, c1), (r2, c2), d, pval))
        else:
            if d < 20:
                significant_pairs.append(((r1, c1), (r2, c2), d, pval))
    return significant_pairs


# 使用显著性位点对构建无向图
def build_graph(significant_pairs: list) -> nx.Graph:
    G = nx.Graph()
    for (r1, c1), (r2, c2), d, pval in significant_pairs:
        G.add_node((r1, c1))
        G.add_node((r2, c2))
        G.add_edge((r1, c1), (r2, c2), weight=d)
    return G


# 使用Floyd_Warshall最短路径算法更新出距离矩阵
def calculate_shortest_paths(G, dist_dict):
    # List nodes and create index mapping
    nodes = list(G.nodes())
    index = {node: i for i, node in enumerate(nodes)}
    n = len(nodes)
    # Initialize adjacency matrix with infinities
    adj = np.full((n, n), np.inf, dtype=float)
    np.fill_diagonal(adj, 0)
    # Fill in weights from graph edges
    for i in range(n):
        u = nodes[i]
        for j in range(i+1, n):
            v = nodes[j]
            d = dist_dict.get((u[0], u[1], v[0], v[1]))
            if d is not None:
                adj[i, j] = d
                adj[j, i] = d
    # Compute shortest paths
    dist_matrix = floyd_warshall(adj, overwrite=True)
    return dist_matrix, nodes



# 计算每个残基的C(vi)
def compute_closeness(dist_matrix, nodes, occurrence_counts):
    n = dist_matrix.shape[0]
    centrality = {}
    for i, node in enumerate(nodes):
        # compute base closeness
        row = dist_matrix[i]
        valid = [row[j] for j in range(n) if j != i and not np.isinf(row[j])]
        c_val = sum(1.0 / (2.0 ** d) for d in valid)
        # adjust for occurrences
        if occurrence_counts and node in occurrence_counts:
            k = occurrence_counts[node]
            if k > 1:
                c_val += (k - 1)
        centrality[node] = c_val
    return centrality



# 划分超簇
def assign_subclusters(nodes, centrality, dist_dict):
    clusters = []
    remaining = set(nodes)
    sorted_nodes = sorted(nodes, key=lambda n: centrality.get(n, 0), reverse=True)
    cluster_idx = 0

    for center in sorted_nodes:
        if center not in remaining:
            continue
        cchain = center[1]
        members = [center]
        for node in remaining:
            if node is center:
                continue
            d = dist_dict.get((center[0], cchain, node[0], node[1]), np.inf)
            # 同链≤10 or 异链≤20
            if (cchain == node[1] and d <= 10) or (cchain != node[1] and d <= 20):
                members.append(node)
        if len(members) > 1:
            Cc = sum(centrality[m] for m in members)
            clusters.append({
                'cluster_name': f"Cluster{cluster_idx}",
                'cluster_sites': members,
                'C(vi)': {m: centrality[m] for m in members},
                'Cc': Cc
            })
            cluster_idx += 1
            
        for m in members:
            remaining.discard(m)

        if not remaining:
            break

    return clusters





def main(pdb_txt_path, db_conn_info_mon, db_conn_info_ppi, output_dir):
    # pr = cProfile.Profile()
    # pr.enable()
    file_paths = [os.path.join(pdb_txt_path, filename) for filename in os.listdir(pdb_txt_path)]
    for file_path in tqdm(file_paths):
        df = pd.read_csv(file_path, sep='\t')
        occ_counts = defaultdict(int)
        for rec in df.to_dict('records'):
            key = (rec['Site_Num'], rec['Chain'])
            occ_counts[key] += 1
        site_pairs, dist_dict = process_file(df, db_conn_info_mon, db_conn_info_ppi)
        sig_pairs = filter_significant(site_pairs)
        G = build_graph(sig_pairs)
        shortest_paths_matrix, nodes = calculate_shortest_paths(G, dist_dict)
        centrality = compute_closeness(shortest_paths_matrix, nodes, occurrence_counts=occ_counts)
        clusters = assign_subclusters(nodes, centrality, dist_dict)
        # print(clusters)
        
        base = os.path.splitext(os.path.basename(file_path))[0]
        out_txt = os.path.join(output_dir, f"{base}_clusters.txt")
        os.makedirs(output_dir, exist_ok=True)

        with open(out_txt, 'w') as fw:
            fw.write("PDB\tSite_Inf\tChain\tSite_Num\tCluster\tC(vi)\tCc\n")
            # 对每个小簇
            for c in clusters:
                name   = c['cluster_name']
                sites  = c['cluster_sites']
                civs   = c['C(vi)']
                Cc_val = c['Cc']
                for node in sites:
                    r, chain = node
                    civ = civs[node]
                    # 找原始文件中对应的行，保留 PDB 和 Site_Inf
                    subdf = df[(df['Site_Num']==r) & (df['Chain']==chain)]
                    for _, row in subdf.iterrows():
                        pdb     = row['PDB']
                        siteinf = row.get('Site_Inf', '')
                        fw.write(
                            f"{pdb}\t{siteinf}\t{chain}\t{r}\t"
                            f"{name}\t{civ:.4f}\t{Cc_val:.4f}\n"
                        )
        # delete processed file
        os.remove(file_path)
    # pr.disable()
    # s = io.StringIO()
    # ps = pstats.Stats(pr, stream=s).sort_stats('cumtime')
    # ps.print_stats(5)
    # print(s.getvalue())

if __name__ == "__main__":
    pdb_txt_path = r'/home/zjliang/users/zhengjiani/PPI_Cluster/sites_txt_new'
    # cif_txt_path = r''
    db_conn_info_mon = "dbname='zjn_ppi_db' user='postgres' password='123456' host='localhost' port='5432'"
    db_conn_info_ppi = "dbname='zjn_ppi_db' user='postgres' password='123456' host='localhost' port='5432'"
    output_dir = r'/home/zjliang/users/zhengjiani/PPI_Cluster/cluster_result_new'
    main(pdb_txt_path, db_conn_info_mon, db_conn_info_ppi, output_dir)
    
