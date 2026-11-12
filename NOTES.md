# quantum-proj: Notes

**Goal:** pick a small set of genes that separates cancer from normal cells, written as a **QUBO**, then solved with brute force → CPLEX → QAOA.

---

## Setup

- Project folder: `~/Desktop/QuantumOptimization/quantum-proj`
- Environment: **qproj** at `~/venvs/qproj` (outside iCloud, so it starts fast)
  - activate: `source ~/venvs/qproj/bin/activate`
  - installed: `scanpy qiskit ipykernel`
- If the kernel doesn't appear in VS Code:
  ⌘⇧P → `>Python: Select Interpreter` → **Enter interpreter path** → `/Users/dariastrait/venvs/qproj/bin/python`
- Lesson learned: environments (venvs) don't belong on the Desktop (iCloud sync makes them slow and breaks kernels)

---

## The data: melanoma (Tirosh et al. 2016, GEO **GSE72056**)

- File: `data/GSE72056_melanoma_single_cell_revised_v2.txt.gz` (don't unzip it)
- **4,645 cells**, **23,686 genes**, **19 patients** (all have melanoma)
- Values are already normalized (log-scaled), not raw counts
- Genes are **rows**, cells are **columns** in the file (we flipped it)
- The first 3 rows are the **answer key**:
  - `tumor`: patient ID (labels, not amounts; IDs skip numbers)
  - `malignant`: 1 = normal, 2 = cancer, 0 = unresolved
  - `cell_type`: 1 = T, 2 = B, 3 = Macrophage, 4 = Endothelial, 5 = CAF, 6 = NK, 0 = cancer/unknown

**Counts:** T cells 2,068 · cancer/unknown 1,758 · B 515 · Macrophage 126 · Endothelial 65 · CAF 61 · NK 52
Patients contribute very unequal numbers of cells (63 to 896).

---

## What we learned

- **0 in the table** = no RNA messages caught for that gene in that cell. Sometimes the gene is on but was missed (a **dropout**).
- **mRNA** = a copy of a gene's recipe (letters A, U, G, C), used to build a protein. We measure **how many copies** per gene per cell.
- More normal cells than cancer cells, because many cells were **sorted for CD45** (immune marker). Lesson: *how the data was collected shapes what you see.*
- **UMAP:** each dot is a cell; closeness means similar genes; the axes mean nothing.
  - Normal cells from different patients **mix**
  - Cancer cells form **islands by patient** (each person's cancer is unique)
- **"Highly variable" ≠ "useful for my question."** PTPRC separates normal from cancer cells 25× (5.23 vs 0.21), yet ranked **17,374th** in variability. The method never looks at labels.
- Without labels, scientists find cancer cells with clues like **broken chromosomes (inferCNV)**, **patient-specific clusters**, and **marker genes**.

---

## Code so far (`explore.ipynb`)

```python
import pandas as pd, numpy as np
import anndata as ad, scanpy as sc

# 1. Load and split the answer key from the genes
df = pd.read_csv("data/GSE72056_melanoma_single_cell_revised_v2.txt.gz", sep="\t", index_col=0)
labels = df.iloc[:3]
genes = df.iloc[3:]
key = labels.T
key.columns = ["tumor", "malignant", "cell_type"]

# 2. AnnData (cells × genes)
adata = ad.AnnData(X=genes.T.values.astype("float32"))
adata.obs_names = genes.columns
adata.var_names = genes.index
adata.var_names_make_unique()
for col in ["tumor", "malignant", "cell_type"]:
    adata.obs[col] = key[col].values

# 3. Highly variable genes + known markers
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
markers = ["PTPRC", "CD3E", "CD19", "CD14", "NKG7", "PECAM1", "COL1A1", "MLANA", "PMEL"]
markers = [g for g in markers if g in adata.var_names]
adata.var.loc[markers, "highly_variable"] = True        # total: 2,007

# 4. PCA + UMAP
adata_hv = adata[:, adata.var["highly_variable"]].copy()
sc.pp.pca(adata_hv, n_comps=50)
sc.pp.neighbors(adata_hv)
sc.tl.umap(adata_hv)
adata_hv.obs["malignant"] = adata_hv.obs["malignant"].astype(int).astype(str)
adata_hv.obs["tumor"] = adata_hv.obs["tumor"].astype(int).astype(str)

# 5. Train/test split (never grade on the cells you studied)
from sklearn.model_selection import train_test_split
train_idx, test_idx = train_test_split(adata_hv.obs_names, test_size=0.5,
                                       random_state=42, stratify=adata_hv.obs["cell_type"])
# 2,322 train, 2,323 test

# 6. Candidate genes from the TRAINING half only (unresolved cells dropped)
train = adata_hv[train_idx].copy()
train = train[train.obs["malignant"] != "0"].copy()
sc.tl.rank_genes_groups(train, groupby="malignant", method="wilcoxon")
top_cancer = sc.get.rank_genes_groups_df(train, group="2").head(8)
top_normal = sc.get.rank_genes_groups_df(train, group="1").head(8)
candidates = list(top_cancer["names"]) + list(top_normal["names"])

# 7. Redundancy (pair penalty) and relevance (reward)
expr = train[:, candidates].X
corr = np.corrcoef(expr.T)                       # 16 × 16
y = (train.obs["malignant"] == "2").astype(int).values
r = np.array([abs(np.corrcoef(expr[:, i], y)[0, 1]) for i in range(len(candidates))])
```

---

## The 16 candidate genes (= 16 qubits)

| Cancer-high | r | Normal-high | r |
|---|---|---|---|
| APOD | 0.82 | PTPRC | 0.76 |
| SPARC | 0.78 | CD74 | 0.63 |
| PMEL | 0.78 | CD3E | 0.50 |
| MLANA | 0.72 | CABP4 | 0.38 |
| TIMP1 | 0.67 | UGDH-AS1 | 0.37 |
| CALD1 | 0.66 | LYZ | 0.36 |
| CCND1 | 0.65 | MAB21L3 | 0.34 |
| CTSL1 | 0.64 | CCL22 | 0.32 |

- PMEL ↔ MLANA correlation = **0.8**, so they're redundant witnesses
- The correlation heatmap shows **two warm blocks** (each team) and **two blue blocks** (cancer vs. normal)
- CD3E is only 0.50 because it marks **T cells only**, while PTPRC marks **all** immune cells

---

## The QUBO

$$\min_x \; -\sum_i r_i x_i \;+\; \alpha \sum_{i<j} |R_{ij}|\, x_i x_j \;+\; \lambda\Big(\sum_i x_i - k\Big)^2$$

- **−relevance:** we want it big, so we negate it
- **+redundancy:** a cost for picking similar pairs
- **+λ(Σx − k)²:** a penalty for not picking exactly k genes

---

## Next time

1. Choose **k** (how many genes, e.g. 4) and the weights **α**, **λ**
2. Build the **Q matrix** (expand (Σx − k)² like in Homework 1; mind the factor of 2 off the diagonal)
3. Solve: **brute force** (2¹⁶ = 65,536 options) → **CPLEX** → **QAOA**
4. **Test:** train a simple classifier on the chosen genes (training half), score it on the **test half**
5. Compare the QUBO's choice with simply taking the top-k genes by relevance

---

*Note: the Alignment method belongs to a professor, so ask them before using it in this project.*
