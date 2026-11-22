# Findings: comparing three cancers

What the comparison of melanoma, head and neck, and breast cancer has shown so far.
Methods are described in the [README](README.md); code is in `h.py` and `compare_all.ipynb`.

## Datasets at a glance

| | Melanoma | Head and neck | Breast |
|---|---|---|---|
| Source | Tirosh 2016 (GSE72056) | Puram 2017 (GSE103322) | Chung 2017 (GSE75688) |
| Cells / patients | 4,645 / 19 | 5,902 / 18 | 515 / 11 |
| Unresolved cells | 132 | 324 | 0 |
| Normal cells are mostly | immune (CD45-sorted) | fibroblasts and T cells | B cells, T cells, myeloid |
| Best single gene (relevance) | APOD, 0.82 | KRT5, 0.96 | KRT18, 0.78 |

## 1. The same cell structure in all three cancers

- In UMAP, cancer cells form **one island per patient**, while normal cells from different patients **mix**.
  Each patient's tumor is unique; normal cell types are shared.
- Unresolved cells sit at the borders between groups (melanoma, head and neck).
- Isomap (breast) shows a nearly linear split between cancer and normal cells, and pulls the patient islands into one
  cancer region, except one patient (BC05) whose tumor is far from all others.

## 2. Top genes reflect each tumor's tissue of origin

| | Cancer-side genes | Normal-side genes |
|---|---|---|
| Melanoma | S100B, PRAME, PMEL, MLANA (pigment cell lineage) | PTPRC, LAPTM5, CXCR4 |
| Head and neck | KRT5, KRT14, KRT17 (squamous keratins) | PTPRC, LYZ, CCL22 |
| Breast | KRT8, KRT18, KRT19, EPCAM, ERBB2 (glandular epithelium) | RGS1, PTPRC, CD53 |

PTPRC, the pan-immune marker, is a top normal-side gene in all three datasets.

## 3. How candidates are chosen matters more than how the final pair is picked

Melanoma, 80/20 split, 10 repeated splits:

| Candidate pool | Best pair | Balanced accuracy |
|---|---|---|
| Hand-picked markers + highly variable genes | APOD, PTPRC | 0.931 ± 0.008 |
| Wilcoxon on all genes / lasso / stability selection | S100B, TMSB4X | 0.976 ± 0.007 |

- Filtering to highly variable genes is label-blind and removed useful genes (PTPRC ranks 17,374th in variability
  in melanoma). Hand-picked markers only partly compensated.
- Without hand-picking, the selection recovered clinically used melanoma markers (S100B, PRAME).
- In breast, the automatic pool shared only 3 of 16 candidates with the hand-picked pool, and its best pair
  (TFAP2C, PTPRF, 0.924) was also chosen independently by the lasso from all genes.

## 4. QUBO vs classical selection

First runs (hand-picked pool, 50/50 split, α = 1, λ = 2):

| | QUBO | Top 2 by relevance | Outcome |
|---|---|---|---|
| Melanoma | APOD, PTPRC: 0.931 | APOD, SPARC: 0.910 | QUBO better |
| Head and neck | KRT5, PTPRC: 0.989 | KRT5, KRT17: 0.995 | tie (one gene already gives 0.989) |
| Breast | KRT18, RGS1: 0.902 | KRT18, EPCAM: 0.905 | tie |

- **Redundancy should count only positive correlation.** Penalizing |correlation| treated a cancer gene and an
  immune gene as redundant, although together they separate the groups best.
- **Correlation is not the same as redundant information.** In breast, two correlated cancer genes outperformed a
  cancer + immune pair (0.924 vs 0.895); the QUBO switches between them at α ≈ 0.10–0.15.
- **The best redundancy weight depends on the dataset** (melanoma needs it, breast does not), so α is now chosen by
  cross-validation within the training cells. With well-chosen candidates, cross-validation selected α = 0 and all
  methods agreed.
- Summary: the QUBO's redundancy term helps most when the candidate pool contains redundant genes.

## 5. QAOA

- The cost Hamiltonian (136 Pauli terms for 16 genes) reproduces the brute-force optimum in every dataset.
- **A large size penalty dominates the energy landscape.** With λ = 2, melanoma and head and neck had identical
  optimal p = 1 angles (β = 0.327, γ = 0.167) and nearly identical probabilities: QAOA was mostly learning the
  size constraint rather than the gene choice. The smallest provably safe λ is about half as large (0.85–0.96).
- **Average energy can hide progress.** From p = 1 to p = 2, P(optimal panel) rose about 11× in melanoma while the
  average energy changed little.
- **Angle optimization can get stuck.** In breast, p = 2 reached a much deeper minimum (average energy 0.59,
  P(optimal) 0.004, 46% of shots with exactly 2 genes) than in the other datasets, suggesting earlier runs stopped in
  local minima. Multi-start optimization with INTERP initialization is now used.

## Open questions

- How much do the smaller λ and the XY mixer (which keeps exactly k genes) improve QAOA?
- Do larger panels (k = 4, 10) generalize better across splits and datasets?
- Can a different redundancy measure (e.g. mutual information) make the QUBO track accuracy more closely?
- Does the method extend to harder tasks, such as separating immune cell types?
