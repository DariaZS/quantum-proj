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
- Also: head and neck (Puram et al. 2017, GSE103322) and breast (Chung et al. 2017, GSE75688)
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

- [x] Environment: scanpy, qiskit
- [x] Three datasets loaded with the same label format (melanoma, head and neck, breast)
- [x] Train/test split; 16 candidate genes; relevance and redundancy
- [x] Build Q, brute force, evaluate on the test half, compare with top-k
- [x] Positive-only redundancy; 10-split fairness check
- [x] Cost Hamiltonian checked against brute force; QAOA p = 1 (grid) and p = 2 (COBYLA)
- [x] Candidate pools without hand-picking (Wilcoxon on all genes, lasso, stability selection)
- [x] `h.py`: safe automatic λ, α by cross-validation, QAOA simulator with INTERP, multi-start, XY mixer
- [ ] `compare_all.ipynb`: all three datasets side by side (melanoma done)
- [ ] QAOA with the new λ and the XY mixer
- [ ] Noise term d (dropouts)
- [ ] Larger panels (k = 4, 10)
- [ ] Task B (immune types)
- [ ] Write-up and presentation

---

## Later / side quests

- **Anomaly detection** (Stolz et al. 2020, geometric anomaly detection): are the "unresolved" (malignant = 0) cells geometric anomalies sitting between groups? Start with ~1,000 cells, ~10 PCA dimensions.
  - **Observation (PCA / Isomap / UMAP comparison):** the uncertain cells mostly live **at borders and in gaps**.
    - **Normal, unknown type** (malignant = 1, cell type = 0; **416 cells**): in UMAP they form a **bridge between B cells and macrophages** and sit around the edges of the T-cell island; in PCA they're scattered in the middle.
    - **Unresolved** (malignant = 0; **132 cells**): in UMAP some sit **inside cancer islands** and some **inside the T-cell island**; in PCA a line of them runs along the **T-cell / cancer border**. (Possibly doublets: two cells captured together.)
    - The answer key isn't perfectly tidy: 28 cells labeled "T" are also marked unresolved or cancer.
  - **Testable question:** do geometric anomaly scores rank the unknown-type and unresolved cells higher than well-labeled cells?
  - Label code used:
    ```python
    mal = adata_hv.obs["malignant"].astype(int)
    ct = adata_hv.obs["cell_type"].astype(int).map(names)
    adata_hv.obs["label"] = np.where(mal == 2, "cancer", np.where(mal == 0, "unresolved", ct))
    ```
- **Embedding comparison (done once, k = 5 for Isomap):** PCA keeps global distances but squashes T and B together; Isomap sits in between; UMAP gives every rare type its own island but exaggerates gaps (distances between islands mean little). UMAP answers "which groups exist," Isomap "how they connect."
- **Compare noise-aware dimensionality reduction** on the melanoma data with its answer key: PCA vs. **scVI**, **GLM-PCA**, **ZINB-WaVE**, **PHATE**. Look for where they all fail; that gap is where a new method could live.
- **Long-term idea:** a dimensionality reduction that models noise directly, plus a way to judge which embedding to trust.

> The Alignment method belongs to a professor. Ask before using or extending it in this project.

---

## References

- Tirosh et al. (2016). *Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq.* Science. Data: GEO GSE72056.
- Romero et al. *Quantum Annealing for Enhanced Feature Selection in Single-Cell RNA Sequencing Data Analysis.* [arXiv:2408.08867](https://arxiv.org/abs/2408.08867)
- Stolz et al. (2020). *Geometric anomaly detection in data.* PNAS.
- Farhi, Goldstone, Gutmann (2014). *A Quantum Approximate Optimization Algorithm.* arXiv:1411.4028
