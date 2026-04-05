# ClusTar

A comprehensive toolkit for identifying, clustering, annotating, and visualizing protein-protein interaction (PPI) sites and residue clusters.

<img width="500" height="500" alt="f31bfa3e2085cb119ea31ecb3d055a77" src="https://github.com/user-attachments/assets/6bba4598-17cf-43cf-920a-a6c28f4effee" />


## Features
- Build a residue distance database from PDB/mmCIF structures  
- Detect and cluster protein interaction sites  
- Integrate multi-source functional annotations (mutation, PTM, pocket info, drug targets, etc.)  
- Visualize cluster composition, distribution, and dynamics  
- Analyze mutation frequency, PTM levels, and structural dynamics  

---

## Requirements
- **Python 3.x**
- **PostgreSQL** (for residue distance database)
- **R** (for downstream analysis and visualization)

---

## Workflow

### 1. Build residue distance database
```bash
python build_res_distance_database.py \
  --input-directory /path/to/database_pdb \
  --csv-path /path/to/pdb_inf.txt \
  --conn-string "dbname='xxx' user='xxx' password='xxx' host='localhost' port='5432'" \
  --threshold 8.0 \
  --workers 8
```
Input: PDB/CIF files
Output: Residue distance database

---

### 2. Cluster calculation
```bash
python cluster_sites.py \
  --input-dir /path/to/sites_txt_new \
  --output-dir /path/to/cluster_result_new \
  --db-mon-conn "dbname='xxx' ..." \
  --db-ppi-conn "dbname='xxx' ..."
```
Output example: Cluster assignment, C(vi), Cc

---

### 3. Data processing & annotation
```bash
Process_PPI_cluster_results.R → outputs Cluster_Raw.xlsx, Single_Site_Data_Raw.xlsx
Merge_cluster_with_annotations.R → integrates PDB info, BioLiP pockets, Fpocket, PTM, mutation, COSMIC, DrugBank, ASD, NACCESS
Plot_Cc_density.R → plots Cc distribution density
Plot_cluster_composition.R → pie charts of cluster composition/distribution
Plot_heatmap_mutation_PTM.R → mutation frequency & PTM heatmap
```

---

### 4. Dynamics analysis
- ANM calculation
```bash
python ppi_dynamic.py anm \
  --pdb-dir /path/to/pdb \
  --ppi-info /path/to/dynamic_calc_info.txt \
  --out-singleAA /path/to/results/anm/singleAA_Data \
  --out-cc /path/to/results/anm/CC_Matrix \
  --out-prs /path/to/results/anm/PRS_Matrix
```

---

### 5. Visualization
```bash
Visualize_MD_results.R → visualization of ANM/GNM dynamics results
```

---
