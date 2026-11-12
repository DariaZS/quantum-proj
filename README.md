# quantum-proj

Picking a small gene panel that separates cancer from normal cells in melanoma, three ways:
**classical** (top genes), **QUBO + brute force** (exact), and **QAOA** (Qiskit simulator).

Course project for COMS 6998 Quantum Optimization (Columbia, Fall 2026). Dasha Strait.

## Data

Melanoma single-cell RNA-seq, Tirosh et al. 2016, GEO [GSE72056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056).
4,645 cells, 23,686 genes, 19 patients.

The data is not in this repo. Download `GSE72056_melanoma_single_cell_revised_v2.txt.gz`
from GEO and put it in `data/` (don't unzip it).

## Setup

```bash
python -m venv ~/venvs/qproj          # keep the venv outside iCloud folders
source ~/venvs/qproj/bin/activate
pip install scanpy qiskit scikit-learn scipy ipykernel
```

Then open `explore.ipynb` and pick the `qproj` kernel.

## What's in the notebook

1. Load data, split off the labels, build an AnnData object
2. Highly variable genes + known markers, PCA / UMAP
3. 50/50 train/test split (stratified by cell type)
4. 16 candidate genes (top 8 cancer + top 8 normal, Wilcoxon, train half only)
5. QUBO: reward relevance, penalize redundant pairs, panel size k = 2; brute force over 2^16 panels
6. Positive-only redundancy (opposite genes are not redundant)
7. 10-split fairness check
8. Cost Hamiltonian from the QUBO (x = (1 − Z)/2), checked against brute force
9. QAOA p = 1 (angle grid) and p = 2 (COBYLA)

## Results so far

| Panel | How | Balanced accuracy (10 splits) |
|---|---|---|
| APOD + PTPRC | QUBO, positive-only redundancy | 0.931 ± 0.004 |
| SPARC + MLANA | QUBO, \|R\| redundancy | 0.911 ± 0.008 |
| APOD + SPARC | top 2 by relevance | 0.910 ± 0.005 |

| | Avg energy | P(best panel) | P(2-gene panel) |
|---|---|---|---|
| Random | 81.80 | 0.00002 | 0.0018 |
| QAOA p = 1 | 55.96 | 0.00016 | 0.0185 |
| QAOA p = 2 | 52.09 | 0.00174 | 0.0909 |
| Exact | −1.58 | 1 | 1 |

## Next

QAOA p = 3, smaller λ, shot-based sampling, k = 4, a dropout-noise term, and immune cell types.

## References

- Tirosh et al. (2016), Science. GEO GSE72056.
- Romero et al. (2024), *Quantum Annealing for Enhanced Feature Selection in scRNA-seq.* arXiv:2408.08867
- Farhi, Goldstone, Gutmann (2014), *A Quantum Approximate Optimization Algorithm.* arXiv:1411.4028