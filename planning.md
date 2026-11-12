# Project Plan: Noise-Aware Gene Selection as a QUBO

*COMS 6998 Quantum Optimization · Dasha Strait*

---

## Research question

**Can QUBO-based gene selection that accounts for single-cell noise produce small gene panels that preserve cell-type structure, including rare cell types, better than standard variance-based selection? And how close do QAOA solutions come to the exact optimum?**

---

## Motivation (what I found in the data)

- Standard practice keeps the **most variable genes**. In the melanoma data, **PTPRC** (CD45) separates normal from cancer cells by about **25×** (5.23 vs 0.21), yet ranks **17,374th of 23,686** in variability. **"Variable" is not the same as "useful."**
- Single-cell data is full of **dropouts**: a gene can be on, but its RNA isn't caught, so the table shows 0. Genes with heavy dropout are unreliable witnesses.
- Tumors contain **rare cell types** (here NK: 52 cells, CAF: 61, endothelial: 65) that are easy to lose in filtering and dimensionality reduction.
- Picking genes is a **subset selection** problem: reward useful genes, penalize redundant pairs. That's the same shape as a QUBO.

---

## The QUBO

Variables: **xᵢ = 1** if gene *i* is in the panel, **0** if not.

$$\min_x \; -\sum_i r_i x_i \;+\; \alpha\sum_{i<j}|R_{ij}|\,x_ix_j \;+\; \beta\sum_i d_i x_i \;+\; \lambda\Big(\sum_i x_i - k\Big)^2$$

| Term | Meaning | Sign |
|---|---|---|
| **rᵢ** | **Relevance**: how strongly gene *i* follows the cell groups (correlation with the label) | − (reward) |
| **Rᵢⱼ** | **Redundancy**: correlation between genes *i* and *j* | + (cost) |
| **dᵢ** | **Noise**: dropout rate of gene *i*, the fraction of zeros among cells of the group where the gene is supposed to be on | + (cost) ← new |
| **(Σx − k)²** | Panel must have exactly **k** genes | + (penalty) |

Weights: **α** (redundancy), **β** (noise), **λ** (size). λ must be large enough that breaking the size rule never pays off (same reasoning as choosing M in Homework 1).

---

## Sub-questions

1. **Noise term:** does adding β·dᵢ change *which* genes get picked? Does it improve results?
2. **Structure:** do the selected panels give embeddings (PCA / UMAP) where cell types, **including rare ones**, stay separate?
3. **Quantum:** how do **QAOA** solutions (p = 1, 2, 3) compare with **brute force** and **CPLEX** as the problem grows (8 → 16 → 20 qubits)?

---

## Data

- **Melanoma scRNA-seq**, Tirosh et al. 2016, GEO **GSE72056** (free, NIH)
- 4,645 cells · 23,686 genes · 19 patients
- Answer key: tumor (patient), malignant (1 = normal, 2 = cancer, 0 = unresolved), cell type (T, B, macrophage, endothelial, CAF, NK)
- Values already log-normalized

---

## Method

1. **Split** cells 50/50 into train / test (stratified by cell type, `random_state=42`). Choose everything on train; grade only on test.
2. **Candidate pool** (train only): top genes per group by `rank_genes_groups` (Wilcoxon), giving 8–20 candidates (= qubits).
3. **Compute** rᵢ, Rᵢⱼ, dᵢ on train.
4. **Build Q** (careful with the factor of 2 off the diagonal; add the constant from expanding (Σx − k)²).
5. **Solve:**
   - brute force (2ⁿ panels; fine up to about 20 genes)
   - CPLEX (exact)
   - QAOA in Qiskit (simulator), p = 1, 2, 3
6. **Evaluate** each panel on the **test half**.

### Tasks (increasing difficulty)
- **A.** Cancer vs. normal (warm-up) ← *started*
- **B.** Immune cell types: T, B, macrophage, NK (includes a rare type)

---

## Baselines

- **Top-k most variable genes** (standard practice)
- **Top-k by relevance alone** (no redundancy or noise terms)
- **QUBO without noise term** (β = 0), to isolate what the noise term adds

---

## Metrics

| What | How |
|---|---|
| Classification | Accuracy / balanced accuracy of a simple classifier (logistic regression) on the test half |
| Rare types | Per-class recall for NK, CAF, endothelial |
| Embedding quality | Silhouette score by cell type in PCA / UMAP of the panel |
| Quantum | **Approximation ratio** = QAOA energy / optimal energy; probability of sampling the optimal panel |

---

## Honest limits

- No quantum advantage is claimed. At 16–20 qubits, brute force and CPLEX are instant. The question is **how close QAOA gets** and how that changes with p and size.
- Results come from one dataset with labels that were themselves computed (from inferCNV, clustering, and markers).
- Dropout rate is a simple noise model. Richer models exist (see below).

---

## Progress

- [x] Environment `qproj` (`~/venvs/qproj`): scanpy, qiskit
- [x] Data loaded, answer key separated, AnnData built
- [x] Highly variable genes + known markers (2,007 genes), UMAP
- [x] Train/test split (2,322 / 2,323)
- [x] 16 candidate genes for task A; relevance **r** and redundancy **R** computed
- [ ] Noise term **d**
- [ ] Build Q, choose k, α, β, λ
- [ ] Brute force → CPLEX → QAOA
- [ ] Evaluate on test half, compare with baselines
- [ ] Task B (immune types)
- [ ] Write-up

---

## Later / side quests

- **Anomaly detection** (Stolz et al. 2020, geometric anomaly detection): are the "unresolved" (malignant = 0) cells geometric anomalies sitting between groups? Start with ~1,000 cells, ~10 PCA dimensions.
- **Compare noise-aware dimensionality reduction** on the melanoma data with its answer key: PCA vs. **scVI**, **GLM-PCA**, **ZINB-WaVE**, **PHATE**. Look for where they all fail; that gap is where a new method could live.
- **Long-term idea:** a dimensionality reduction that models noise directly, plus a way to judge which embedding to trust.

> The Alignment method belongs to a professor. Ask before using or extending it in this project.

---

## References

- Tirosh et al. (2016). *Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq.* Science. Data: GEO GSE72056.
- Romero et al. *Quantum Annealing for Enhanced Feature Selection in Single-Cell RNA Sequencing Data Analysis.* [arXiv:2408.08867](https://arxiv.org/abs/2408.08867)
- Stolz et al. (2020). *Geometric anomaly detection in data.* PNAS.
- Farhi, Goldstone, Gutmann (2014). *A Quantum Approximate Optimization Algorithm.* arXiv:1411.4028
