from prody import *
import numpy as np
import pandas as pd 
import os
import math
import concurrent.futures
os.chdir("/home/zjliang/users/zhengjiani/dynamics/code")


three_to_one = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
	 'GLY': 'G', 'GLN': 'Q', 'GLU': 'E', 'HIS': 'H', 'ILE': 'I',
	 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PRO': 'P', 'PHE': 'F',
	 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'TRP': 'W', 'VAL': 'V'}
three_amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLY', 'GLN', 'GLU', 'HIS', 'ILE',
						'LEU', 'LYS', 'MET', 'PRO', 'PHE', 'SER', 'THR', 'TRP', 'TYR', 'VAL']


def make_matrix_to_df(matrix, res_info):
    return pd.DataFrame(matrix, index=res_info, columns=res_info)
 
def get_anm_result_pdb(pdb, structure, chain1, chain2, uniprot1, uniprot2, pairAA_CC_Matrix_dir, singleAA_sq_dir, pairAA_PRS_Matrix_dir):
    # 选择非离子原子和CA原子，这些原子需要属于蛋白质并且在指定的两条链上
    choose = 'not ion and name CA and protein and chain ' + chain1 + ' ' + chain2
    calphas = structure.select(choose) # 返回所有CA原子组成的AtomGroup
    res = calphas.getResnames()   # 提取每个CA原子对应残基的三字母列表
    for i in range(len(res)):
        if res[i] in three_amino_acids:
            res[i] = three_to_one[res[i]]
    res_num = calphas.getResnums()
    res_num = pd.DataFrame(res_num)
    
    # 确定链在PDB文件中的排序
    pdb_chain1 = calphas[0].getChid()    
    pdb_chain2 = calphas[len(calphas)-1].getChid()
    if chain1 != pdb_chain1 and chain2 != pdb_chain2:
        chain1, chain2 = chain2, chain1
        uniprot1, uniprot2 = uniprot2, uniprot1
    
    p1 = structure.select("not ion and name CA and protein and chain " + chain1)  # 选出单链集合
    p2 = structure.select("not ion and name CA and protein and chain " + chain2)
    ppi_seq = [uniprot1] * len(p1.getResindices()) + [uniprot2] * len(p2.getResindices())
    ppi_chain = [chain1] * len(p1.getResindices()) + [chain2] * len(p2.getResindices())
    res_num["chain"] = ppi_chain
    res_num["uniprot"] = ppi_seq
    res_num["res"] = res
    res_num["Res_Info"] = res_num.apply(lambda row: '_'.join([row["uniprot"], row["chain"], row["res"]]) + str(row[0]), axis=1)
    res_info = res_num["Res_Info"]
    
    # ----- MD -----
    set_modes_num = None # 计算模式为None，计算全部模式
    anm = ANM(pdb) #新建ANM对象
    # 基于选出的CA原子构建Hessian矩阵，原子间距离小于或等于15Å的对视作用弹簧相连，弹簧常数默认为1
    # 生成并存储一个大小为(3N×3N)的Hessian矩阵
    anm.buildHessian(calphas, cutoff=15.0, gamma=1) 
    # 计算模式（特征值／特征向量）
    anm.calcModes(n_modes=set_modes_num)
    # 计算MSF
    anm_sq = calcSqFlucts(anm)
    # 计算stiffness(对称矩阵的每一行均值)
    anm_stiffness = np.mean(calcMechStiff(anm, calphas), axis=1)
    
    mode_num = anm._n_modes
    if mode_num < 5:
        return "Residue number is less than 5"
    # 计算所有残基之间的交叉相关矩阵(N×N)
    print('cal anm_cc...')
    anm_top3_cc = calcCrossCorr(anm[0:3])
    anm_5_per_cc = calcCrossCorr(anm[0:math.ceil(mode_num*0.05)])
    anm_5_20_per_cc = calcCrossCorr(anm[math.ceil(mode_num*0.05):math.ceil(mode_num*0.2)])
    anm_20_50_per_cc = calcCrossCorr(anm[math.ceil(mode_num*0.2):math.ceil(mode_num*0.5)])
    anm_greate_60_per_cc = calcCrossCorr(anm[math.ceil(mode_num*0.6):])
    
    # 计算PRS
    print('cal anm_prs...')
    anm_prs, anm_effectiveness, anm_sensitivity = calcPerturbResponse(anm)
    
    # 将CC矩阵封装为DataFrame
    anm_top3_cc, anm_5_per_cc, anm_5_20_per_cc, anm_20_50_per_cc, anm_greate_60_per_cc = [
    make_matrix_to_df(mat, res_info) for mat in (anm_top3_cc, anm_5_per_cc, anm_5_20_per_cc, anm_20_50_per_cc, anm_greate_60_per_cc)]
    # 将PRS矩阵封装为DataFrame
    anm_prs = make_matrix_to_df(anm_prs, res_info)
    
    
    # 保存CC矩阵和PRS矩阵
    cc_list = [
        (anm_top3_cc, "anm_top3_cc"),
        (anm_5_per_cc, "anm_5_per_cc"),
        (anm_5_20_per_cc, "anm_5_20_per_cc"),
        (anm_20_50_per_cc, "anm_20_50_per_cc"),
        (anm_greate_60_per_cc, "anm_greate_60_per_cc"),
    ]
    for df, suffix in cc_list:
        fname = f"{pdb}_{uniprot1}_{uniprot2}_{chain1}_{chain2}_{suffix}.csv"
        df.to_csv(os.path.join(pairAA_CC_Matrix_dir, fname), index=True)
    prs_list = [
        (anm_prs, "anm_prs")
    ]
    for df, suffix in prs_list:
        fname = f"{pdb}_{uniprot1}_{uniprot2}_{chain1}_{chain2}_{suffix}.csv"
        df.to_csv(os.path.join(pairAA_PRS_Matrix_dir, fname), index=True)
    
        
    # 保存MSF(sq)/effectiveness/sensitivity/stiffness文件
    # Build file path
    csv_name = f"{pdb}_{chain1}_{chain2}_singleAA_data.csv"
    csv_path = os.path.join(singleAA_sq_dir, csv_name)
    # 把残基信息和sq数据拼成一个df
    res_num_with_sq = pd.DataFrame(res_info.copy())
    res_num_with_sq["ANM_sq"] = anm_sq
    res_num_with_sq["ANM_effectiveness"] = anm_effectiveness
    res_num_with_sq["ANM_sensitivity"] = anm_sensitivity
    res_num_with_sq["ANM_stiffness"] = anm_stiffness

    if os.path.exists(csv_path):
        # Load existing, add/overwrite column
        df = pd.read_csv(csv_path)
        if len(anm_sq) != len(df):
            raise ValueError(f"Length of data ({len(anm_sq)}) does not match "
                             f"existing rows ({len(df)}).")
        df["ANM_sq"] = anm_sq
        df["ANM_effectiveness"] = anm_effectiveness
        df["ANM_sensitivity"] = anm_sensitivity
        df["ANM_stiffness"] = anm_stiffness
    else:
        # Create fresh DataFrame
        df = res_num_with_sq

    # Save out
    df.to_csv(csv_path, index=False)


def get_anm_result_cif(pdb, structure, chain1, chain2, uniprot1, uniprot2, pairAA_CC_Matrix_dir, singleAA_sq_dir, pairAA_PRS_Matrix_dir):
    # 选两条链的CA（使用segname）
    selector = f"not ion and name CA and protein and segname {chain1} {chain2}"
    calphas = structure.select(selector)
    if calphas is None or calphas.numAtoms() == 0:
        raise ValueError(f"[CIF] No CA atoms found for chains {chain1}, {chain2} using selector: {selector}")
    res = calphas.getResnames()
    for i in range(len(res)):
        if res[i] in three_amino_acids:
            res[i] = three_to_one[res[i]]
    res_num = pd.DataFrame(calphas.getResnums())

    # 按segname判断顺序（不要用 getChid，这里 chid=A/B是label_asym_id）
    cif_seg1 = calphas[0].getSegname()
    cif_seg2 = calphas[-1].getSegname()
    if chain1 != cif_seg1 and chain2 != cif_seg2:
        chain1, chain2 = chain2, chain1
        uniprot1, uniprot2 = uniprot2, uniprot1

    # 单链也用segname
    p1 = structure.select(f"not ion and name CA and protein and segname {chain1}")
    p2 = structure.select(f"not ion and name CA and protein and segname {chain2}")
    if p1 is None or p2 is None:
        raise ValueError(f"[CIF] Failed to select CA atoms for chain {chain1} or {chain2}")
    ppi_seq = [uniprot1] * len(p1.getResindices()) + [uniprot2] * len(p2.getResindices())
    ppi_chain = [chain1] * len(p1.getResindices()) + [chain2] * len(p2.getResindices())
    res_num["chain"] = ppi_chain
    res_num["uniprot"] = ppi_seq
    res_num["res"] = res
    res_num["Res_Info"] = res_num.apply(lambda row: '_'.join([row["uniprot"], row["chain"], row["res"]]) + str(row[0]), axis=1)
    res_info = res_num["Res_Info"]
    # ----------- ANM ----------- #
    anm = ANM(pdb)
    anm.buildHessian(calphas, cutoff=15.0, gamma=1)
    anm.calcModes(n_modes=None)
    anm_sq = calcSqFlucts(anm)
    anm_stiffness = np.mean(calcMechStiff(anm, calphas), axis=1)
    mode_num = anm._n_modes
    if mode_num < 5:
        return "Residue number is less than 5"

    print('cal anm_cc...')
    anm_top3_cc = calcCrossCorr(anm[0:3])
    anm_5_per_cc = calcCrossCorr(anm[0:math.ceil(mode_num*0.05)])
    anm_5_20_per_cc = calcCrossCorr(anm[math.ceil(mode_num*0.05):math.ceil(mode_num*0.2)])
    anm_20_50_per_cc = calcCrossCorr(anm[math.ceil(mode_num*0.2):math.ceil(mode_num*0.5)])
    anm_greate_60_per_cc = calcCrossCorr(anm[math.ceil(mode_num*0.6):])

    print('cal anm_prs...')
    anm_prs, anm_effectiveness, anm_sensitivity = calcPerturbResponse(anm)

    # dataframe 封装
    anm_top3_cc, anm_5_per_cc, anm_5_20_per_cc, anm_20_50_per_cc, anm_greate_60_per_cc = [
        make_matrix_to_df(mat, res_info) for mat in (
            anm_top3_cc, anm_5_per_cc, anm_5_20_per_cc, anm_20_50_per_cc, anm_greate_60_per_cc)
    ]
    anm_prs = make_matrix_to_df(anm_prs, res_info)

    # 保存 CC 和 PRS
    cc_list = [
        (anm_top3_cc, "anm_top3_cc"),
        (anm_5_per_cc, "anm_5_per_cc"),
        (anm_5_20_per_cc, "anm_5_20_per_cc"),
        (anm_20_50_per_cc, "anm_20_50_per_cc"),
        (anm_greate_60_per_cc, "anm_greate_60_per_cc"),
    ]
    for df, suffix in cc_list:
        fname = f"{pdb}_{uniprot1}_{uniprot2}_{chain1}_{chain2}_{suffix}.csv"
        df.to_csv(os.path.join(pairAA_CC_Matrix_dir, fname), index=True)
    
    prs_list = [(anm_prs, "anm_prs")]
    for df, suffix in prs_list:
        fname = f"{pdb}_{uniprot1}_{uniprot2}_{chain1}_{chain2}_{suffix}.csv"
        df.to_csv(os.path.join(pairAA_PRS_Matrix_dir, fname), index=True)

    # 保存单链数据
    csv_name = f"{pdb}_{chain1}_{chain2}_singleAA_data.csv"
    csv_path = os.path.join(singleAA_sq_dir, csv_name)
    res_num_with_sq = pd.DataFrame(res_info.copy())
    res_num_with_sq["ANM_sq"] = anm_sq
    res_num_with_sq["ANM_effectiveness"] = anm_effectiveness
    res_num_with_sq["ANM_sensitivity"] = anm_sensitivity
    res_num_with_sq["ANM_stiffness"] = anm_stiffness

    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        if len(anm_sq) != len(df):
            raise ValueError(f"Length mismatch: {len(anm_sq)} vs existing {len(df)}.")
        df["ANM_sq"] = anm_sq
        df["ANM_effectiveness"] = anm_effectiveness
        df["ANM_sensitivity"] = anm_sensitivity
        df["ANM_stiffness"] = anm_stiffness
    else:
        df = res_num_with_sq
    df.to_csv(csv_path, index=False)


        

def main():
    pdb_path = r"/home/zjliang/users/zhengjiani/dynamics/data/pdb"
    PPI_str_info = "/home/zjliang/users/zhengjiani/dynamics/data/dynamic_calc_info.txt"
    singleAA_results_dir = r"/home/zjliang/users/zhengjiani/dynamics/data/results/anm/singleAA_Data"
    pairAA_CC_Matrix_dir = r"/home/zjliang/users/zhengjiani/dynamics/data/results/anm/CC_Matrix"
    pairAA_PRS_Matrix_dir = r"/home/zjliang/users/zhengjiani/dynamics/data/results/anm/PRS_Matrix"
    os.makedirs(singleAA_results_dir, exist_ok=True)
    os.makedirs(pairAA_CC_Matrix_dir, exist_ok=True)
    os.makedirs(pairAA_PRS_Matrix_dir, exist_ok=True)
    
    pdb_file_list = os.listdir(pdb_path)
    singleAA_processed_list = [f for f in os.listdir(singleAA_results_dir) if f.endswith("_singleAA_data.csv")]
    singleAA_processed_ppi = [os.path.basename(s.split('_')[0]) for s in singleAA_processed_list]
    
    with open(PPI_str_info, "r") as f:
        lines = f.readlines()
    for i in lines:
        try:
            uniprot1, uniprot2, pdb, chain1, chain2 = i.strip().split("\t")
            if pdb + ".pdb" in pdb_file_list:
                pdb_file = os.path.join(pdb_path, pdb+".pdb")
            elif pdb + ".cif" in pdb_file_list:
                pdb_file = os.path.join(pdb_path, pdb+".cif")
            else:
                print(f"No PDB or CIF file found for {pdb}")
                continue
            print(f"Processing {pdb_file} with chains {chain1}, {chain2}")
            
            if pdb not in singleAA_processed_ppi:
                if pdb_file.endswith('.cif'):
                    structure = parseMMCIF(pdb_file)
                    get_anm_result_cif(pdb, structure, chain1, chain2, uniprot1, uniprot2, pairAA_CC_Matrix_dir, singleAA_results_dir, pairAA_PRS_Matrix_dir)
                else:
                    structure = parsePDB(pdb_file)
                    get_anm_result_pdb(pdb, structure, chain1, chain2, uniprot1, uniprot2, pairAA_CC_Matrix_dir, singleAA_results_dir, pairAA_PRS_Matrix_dir)
            else:
                print(f"{pdb} has processed in result directory")
                continue
            
        except ValueError as ve:
            print(f"Skipping line due to ValueError: {i.strip()}. Error: {ve}")
        except Exception as e:
            print(f"Error processing line: {i.strip()}. Error: {e}")
            
        
        

if __name__ == '__main__':
    main()







