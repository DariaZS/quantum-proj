# Gene panel selection with QUBO and QAOA

Selecting a 2-gene panel that separates cancer from normal cells in single-cell RNA-seq data, comparing
classical feature selection, an exactly solved QUBO, and QAOA, across three tumor types.

COMS 6998 Quantum Optimization, Columbia University, Fall 2026. Dasha Strait.

## Repository

| | |
|---|---|
| `compare_all.ipynb` | Main analysis: each dataset, then a side-by-side comparison |
| `h.py` | Data loading, gene selection, QUBO, solvers, evaluation, QAOA simulator |
| `NOTES.md` | Findings so far |
| `planning.md` | Project plan |
| `archive/` | Earlier per-dataset notebooks |

## Data

| Dataset | GEO | Cells | Patients |
|---|---|---|---|
| Melanoma (Tirosh et al. 2016) | [GSE72056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056) | 4,645 | 19 |
| Head and neck (Puram et al. 2017) | [GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322) | 5,902 | 18 |
| Breast (Chung et al. 2017) | [GSE75688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688) | 515 | 11 |

Data files are not included. Download the supplementary `.txt.gz` files from GEO into `data/`.

## Setup

Python 3.13.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Method

All selection steps use the training cells only (80/20 split).

1. **Candidate pool (16 genes = 16 qubits)**, chosen four ways: highly variable genes plus known markers,
   Wilcoxon test on all genes, lasso, and stability selection.
2. **Panel selection (2 genes):**
   - Classical: top 2 genes by relevance; lasso on all genes.
   - QUBO: relevance reward, redundancy penalty (positive correlation), and a size penalty, solved exactly.
     The redundancy weight is set by cross-validation; the size penalty is the smallest value that guarantees 2 genes.
   - QAOA: the same QUBO, statevector-simulated, with the standard X mixer and a Hamming-weight-preserving XY mixer.
3. **Evaluation:** class-balanced logistic regression, balanced accuracy on the test set and over 10 random splits.

## Results

Balanced accuracy, 10 random 80/20 splits. Best panel per candidate pool:

| | Hand-picked pool | Automatic pools |
|---|---|---|
| Melanoma | APOD, PTPRC: 0.951 (top 2: 0.927) | S100B, TMSB4X: 0.982 |
| Head and neck | KRT5, KRT17: 0.993 | KRT5, KRT17: 0.993 |
| Breast | EPCAM, KRT8: 0.898 | PTPRF, TFAP2C: 0.919 |

- Candidate pool choice matters more than the selection method; without hand-picking, all methods converge.
- The QUBO beats top-2 selection only when the pool contains redundant genes (melanoma, hand-picked pool).
- QAOA rarely concentrates on the exact optimum: the XY mixer starts at P(optimal) = 1/120 and its layers change
  this by 0.4× to 1.7× across datasets, likely because the optimizer minimizes average energy, not P(optimal).

See `NOTES.md` for details.

## References

- Tirosh et al. (2016) *Science* 352:189. Puram et al. (2017) *Cell* 171:1611. Chung et al. (2017) *Nat. Commun.* 8:15081.
- Farhi, Goldstone, Gutmann (2014). A Quantum Approximate Optimization Algorithm. arXiv:1411.4028
- Hadfield et al. (2019). From the QAOA to a Quantum Alternating Operator Ansatz. *Algorithms* 12:34
- Zhou et al. (2020). QAOA: Performance, Mechanism, and Implementation on Near-Term Devices. *PRX* 10:021067
- Meinshausen & Bühlmann (2010). Stability selection. *JRSS B* 72:417
- Peng, Long, Ding (2005). Feature selection based on mutual information (mRMR). *IEEE TPAMI* 27:1226
- Romero et al. (2024). Quantum annealing for feature selection in scRNA-seq. arXiv:2408.08867
