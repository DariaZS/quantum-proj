# Findings

Comparison of melanoma, head and neck, and breast cancer. Methods: see [README](README.md).
All numbers: 80/20 split, class-balanced logistic regression, balanced accuracy over 10 random splits.

## Datasets

| | Melanoma | Head and neck | Breast |
|---|---|---|---|
| Source | Tirosh 2016 (GSE72056) | Puram 2017 (GSE103322) | Chung 2017 (GSE75688) |
| Cells / patients | 4,645 / 19 | 5,902 / 18 | 515 / 11 |
| Unresolved cells | 132 | 324 | 0 |
| Normal cells mostly | immune (CD45-sorted) | fibroblasts, T cells | B cells, T cells, myeloid |

## 1. Cell structure

- UMAP, all three: cancer cells form one cluster per patient; normal cells from different patients mix.
- Unresolved cells lie at cluster borders (melanoma, head and neck).
- Isomap (breast): a near-linear cancer/normal split; one patient (BC05) far from all others.

## 2. Top genes reflect tissue of origin

| | Cancer side | Normal side |
|---|---|---|
| Melanoma | S100B, PRAME, SOX10, PMEL, MLANA | PTPRC, LAPTM5, TMSB4X |
| Head and neck | KRT5, KRT17, PERP, DSP, SFN | VIM, TMSB4X, B2M |
| Breast | PTPRF, TFAP2C, ERBB3, EPCAM, KRT8 | RGS1, CD53, PTPRC |

Without hand-picking, the selection recovers clinical melanoma markers (S100B, PRAME).

## 3. Candidate pool matters most

| | Hand-picked + HVG | Automatic (Wilcoxon / lasso / stability) |
|---|---|---|
| Melanoma | APOD, PTPRC: 0.951 | S100B, TMSB4X: 0.982 |
| Head and neck | KRT5, KRT17: 0.993 | KRT5, KRT17: 0.993 |
| Breast | EPCAM, KRT8: 0.898 | PTPRF, TFAP2C: 0.919 |

- The highly-variable-gene filter is label-blind and can exclude useful genes (PTPRC ranks 17,374th in melanoma).
- Wilcoxon and stability selection give nearly identical pools; the lasso overlaps partly.
- Head and neck is near ceiling: KRT5 alone separates the classes almost perfectly.

## 4. QUBO vs classical selection

- With automatic pools, cross-validation sets the redundancy weight α to 0 in all three datasets, so the QUBO
  equals top-2 selection. The lasso (choosing from all genes) matches the best automatic pair in melanoma and
  breast, and is within 0.002 in head and neck.
- The QUBO improves on top-2 only when candidates are redundant: melanoma, hand-picked pool, 0.951 vs 0.927.
- Redundancy should use positive correlation only; |correlation| penalizes complementary cancer/immune pairs.
- Correlated genes can still add information (breast: two cancer genes beat cancer + immune), so α is chosen
  per dataset by cross-validation within the training set.
- The class-balanced classifier changed breast's α from 0.05 to 0; an earlier QUBO advantage there
  (0.913 vs 0.887) did not hold.

## 5. QAOA (16 qubits, statevector simulation, stability-selection pool)

P(optimal panel):

| | X mixer, p = 1 / 2 / 3 | XY mixer, p = 1 / 2 | XY start (1/120) |
|---|---|---|---|
| Melanoma | 0.0030 / 0.0037 / 0.0041 | 0.0107 / 0.0142 | 0.0083 |
| Head and neck | 0.0032 / 0.0040 / 0.0047 | 0.0025 / 0.0033 | 0.0083 |
| Breast | 0.0029 / 0.0035 / 0.0038 | 0.0096 / 0.0112 | 0.0083 |

Random guess over all 65,536 bit strings: 0.000015.

- The cost Hamiltonian reproduces the brute-force optimum in every dataset.
- **The XY mixer's advantage comes mostly from its starting state** (uniform over the 120 valid panels).
  Its layers change P(optimal) by 1.7× (melanoma), 1.35× (breast), and 0.4× (head and neck).
- **The optimizer targets average energy, not P(optimal).** With many near-tied panels, lower average energy
  can coexist with a lower chance of the exact optimum.
- **Angle optimization dominated early results.** Energy rescaling, multi-start, and INTERP initialization raised
  P(optimal) about 20× over the first runs (together with the new candidate pool; effects not separated).
- **The size penalty λ had little effect** once angles were well optimized (melanoma: λ = 2 gave 0.0033 / 0.0040
  at p = 1 / 2, the smallest safe λ gave 0.0030 / 0.0037).

## Open questions

- Does a CVaR objective (optimizing the best shots rather than the mean) concentrate QAOA on the optimum?
- Do larger panels (k = 4, 10) generalize better, and is QAOA more competitive there?
- Does the method extend to harder tasks, such as separating immune cell types?
