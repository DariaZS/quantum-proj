"""
h.py  --  helper functions for my project:
classical vs QUBO vs QAOA for choosing small gene panels in single-cell RNA data.

Dasha Strait, COMS 6998 Quantum Optimization, Columbia, Fall 2026.

Typical use in a notebook:

    import h
    adata = h.load_breast("data/GSE75688_final_sample_information.txt.gz",
                          "data/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt.gz")
    adata = h.prepare(adata)                          # drop unresolved cells
    train, test = h.split(adata)                      # same split every time (seed 42)
    cands = h.candidates(adata, train, method="stable")
    alpha, cv_table = h.choose_alpha_cv(adata, train, cands, k=2)
    Q, const, info = h.build_qubo(adata, train, cands, k=2, alpha=alpha, lam="auto")
    panel = h.best_panel(Q, const, cands, k=2)
    h.evaluate(adata, panel, train, test), h.repeated_splits(adata, panel)
    E = h.all_energies(Q, const)
    table = h.qaoa_ladder(E, p_max=4, k=2)

Sections
  1. Loading the three datasets (all come out in the same format)
  2. Preparing and splitting
  3. Choosing candidate genes (no hand-picking needed)
  4. The QUBO (relevance, redundancy, size penalty with a safe automatic lambda)
  5. Exact and classical solvers (brute force, k-subsets, simulated annealing)
  6. Grading panels (test half, repeated splits, cross-validated alpha)
  7. Quantum: Hamiltonian, fast QAOA simulator, X and XY mixers, multi-start, INTERP
  8. One-call pipeline

Conventions
  * Bit i of the integer b is gene i (same order as Qiskit), so index b of the
    energy vector is the panel {i : (b >> i) & 1}.
  * malignant: "2" = cancer, "1" = normal, "0" = unresolved.
  * Every choice (genes, alpha, lambda) uses the TRAINING cells only.
"""

from __future__ import annotations

import itertools
import math
import time

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler

DEFAULT_SEED = 42

DEFAULT_MARKERS = [
    "PTPRC", "CD3E", "CD19", "CD14", "NKG7", "PECAM1", "COL1A1",   # immune / stromal
    "MLANA", "PMEL",                                               # melanoma
    "KRT5", "KRT14", "KRT17",                                      # squamous (head and neck)
    "EPCAM", "KRT8", "KRT18", "KRT19", "ESR1", "ERBB2",            # glandular (breast)
]


# ---------------------------------------------------------------------------
# 1. Loading
# ---------------------------------------------------------------------------

def _to_anndata(X, cells, genes, obs: dict):
    import anndata as ad
    adata = ad.AnnData(X=np.asarray(X, dtype="float32"))
    adata.obs_names = [str(c) for c in cells]
    adata.var_names = [str(g) for g in genes]
    adata.var_names_make_unique()
    for key, val in obs.items():
        adata.obs[key] = np.asarray(val)
    return adata


def load_melanoma(path="data/GSE72056_melanoma_single_cell_revised_v2.txt.gz"):
    """Tirosh et al. 2016 (GSE72056). Values are already log scale."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    key = df.iloc[:3].T
    key.columns = ["tumor", "malignant", "cell_type"]
    genes = df.iloc[3:]
    names = {1: "T", 2: "B", 3: "Macrophage", 4: "Endothelial", 5: "CAF", 6: "NK", 0: "none"}
    return _to_anndata(
        genes.T.values, genes.columns, genes.index,
        {"malignant": key["malignant"].astype(float).astype(int).astype(str).values,
         "cell_type": key["cell_type"].astype(float).astype(int).map(names).values,
         "tumor": key["tumor"].astype(float).astype(int).astype(str).values},
    )


def load_hnscc(path="data/GSE103322_HNSCC_all_data.txt.gz"):
    """Puram et al. 2017 (GSE103322). Values are already log scale.
    Patient IDs written as HN28 and HNSCC28 are merged (18 patients, as in the paper)."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    lab = df.iloc[:5].T
    lab.columns = ["maxima", "lymph_node", "cancer", "normal", "cell_type"]
    genes = df.iloc[5:].astype(float)
    genes.index = genes.index.str.strip("'")
    c = lab["cancer"].astype(float)
    nm = lab["normal"].astype(float)
    malignant = np.where(c == 1, "2", np.where(nm == 1, "1", "0"))
    tumor = (pd.Series(genes.columns).str.split("_").str[0]
             .str.replace(r"^HN(\d)", r"HNSCC\1", regex=True))
    cell_type = lab["cell_type"].astype(str).replace({"0.0": "none", "0": "none"})
    return _to_anndata(
        genes.T.values, genes.columns, genes.index,
        {"malignant": malignant, "cell_type": cell_type.values, "tumor": tumor.values,
         "lymph_node": lab["lymph_node"].astype(float).astype(int).values},
    )


def load_breast(info_path="data/GSE75688_final_sample_information.txt.gz",
                tpm_path="data/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt.gz"):
    """Chung et al. 2017 (GSE75688). Values are TPM, so log1p is applied here.
    Only single cells are kept (the file also has bulk and pooled samples)."""
    info = pd.read_csv(info_path, sep="\t")
    info = info.set_index(info.columns[0])
    info = info[info["type"].astype(str).str.upper() == "SC"]
    lab = info["index"].astype(str).str.lower()
    malignant = np.where(lab.str.startswith("non"), "1", np.where(lab.str.contains("tumor"), "2", "0"))
    info = info.assign(malignant=malignant)
    tpm = pd.read_csv(tpm_path, sep="\t").set_index("gene_name")
    cells = [c for c in info.index if c in tpm.columns]
    patient = pd.Series(cells).str.split("_").str[0]
    X = np.log1p(tpm[cells].T.values.astype("float32"))
    type_col = "index3" if "index3" in info.columns else "index2"
    return _to_anndata(
        X, cells, tpm.index,
        {"malignant": info.loc[cells, "malignant"].values,
         "cell_type": info.loc[cells, type_col].astype(str).values,
         "tumor": patient.str.replace("LN", "", regex=False).values,
         "lymph_node": patient.str.contains("LN").astype(int).values},
    )


# ---------------------------------------------------------------------------
# 2. Preparing and splitting
# ---------------------------------------------------------------------------

def prepare(adata, drop_unresolved=True):
    """Copy of the data without unresolved cells (malignant == "0")."""
    if drop_unresolved:
        adata = adata[adata.obs["malignant"] != "0"]
    return adata.copy()


def labels(adata, cells=None):
    """1 = cancer, 0 = normal."""
    obs = adata.obs if cells is None else adata[cells].obs
    return (obs["malignant"].astype(str) == "2").astype(int).values


def split(adata, seed=DEFAULT_SEED, test_size=0.5, stratify="malignant"):
    """Train / test cell names. Same seed -> same split every time."""
    strat = adata.obs[stratify] if stratify else None
    tr, te = train_test_split(adata.obs_names, test_size=test_size, random_state=seed, stratify=strat)
    return pd.Index(tr), pd.Index(te)


def expr(adata, cells, genes):
    """Dense expression matrix (cells x genes)."""
    X = adata[cells, list(genes)].X
    return np.asarray(X.toarray() if hasattr(X, "toarray") else X, dtype=float)


def detected_genes(adata, cells, min_cells=10):
    """Genes seen (value > 0) in at least min_cells of the given cells."""
    X = adata[cells].X
    seen = np.asarray((X > 0).sum(axis=0)).ravel()
    return adata.var_names[seen >= min_cells]


# ---------------------------------------------------------------------------
# 3. Choosing candidate genes
# ---------------------------------------------------------------------------

def auc_scores(X, y):
    """Vectorized Mann-Whitney / Wilcoxon AUC for every column of X.
    AUC near 1: gene is higher in cancer. Near 0: higher in normal. 0.5: no difference."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).astype(bool)
    n1, n0 = y.sum(), (~y).sum()
    ranks = rankdata(X, axis=0)
    U = ranks[y].sum(axis=0) - n1 * (n1 + 1) / 2
    return U / (n1 * n0)


def _top_by_auc(auc, genes, n_per_side):
    genes = np.asarray(genes)
    return list(genes[np.argsort(-auc)[:n_per_side]]) + list(genes[np.argsort(auc)[:n_per_side]])


def candidates(adata, train, method="stable", n_per_side=8, min_cells=10,
               markers=None, n_hvg=2000, n_boot=100, sub_frac=0.5, seed=DEFAULT_SEED,
               return_scores=False):
    """Choose 2 * n_per_side candidate genes using ONLY the training cells.

    method:
      "markers_hvg": highly variable genes + hand-picked markers, then Wilcoxon (my first way)
      "wilcoxon":    Wilcoxon AUC on all detected genes (no hand-picking)
      "lasso":       L1 logistic regression on all detected genes
      "stable":      stability selection: Wilcoxon AUC on many random half-samples of the
                     training cells; keep the genes that make the top list most often
                     (Meinshausen & Buhlmann 2010). Most reliable without hand-picking.
    Returns the list: first n_per_side are cancer-side genes, then normal-side genes.
    """
    y = labels(adata, train)
    genes = detected_genes(adata, train, min_cells)

    if method == "markers_hvg":
        import scanpy as sc
        sub = adata[train].copy()
        sc.pp.highly_variable_genes(sub, n_top_genes=n_hvg)
        mk = [g for g in (markers or DEFAULT_MARKERS) if g in sub.var_names]
        sub.var.loc[mk, "highly_variable"] = True
        pool = [g for g in sub.var_names[sub.var["highly_variable"]] if g in set(genes)]
        auc = auc_scores(expr(adata, train, pool), y)
        out = _top_by_auc(auc, pool, n_per_side)
        scores = pd.Series(auc, index=pool)

    elif method == "wilcoxon":
        auc = auc_scores(expr(adata, train, genes), y)
        out = _top_by_auc(auc, genes, n_per_side)
        scores = pd.Series(auc, index=genes)

    elif method == "lasso":
        X = StandardScaler().fit_transform(expr(adata, train, genes))
        coef = None
        for C in np.logspace(-3, 1, 40):
            m = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=5000).fit(X, y)
            coef = m.coef_[0]
            if (coef > 0).sum() >= n_per_side and (coef < 0).sum() >= n_per_side:
                break
        g = np.asarray(genes)
        out = list(g[np.argsort(-coef)[:n_per_side]]) + list(g[np.argsort(coef)[:n_per_side]])
        scores = pd.Series(coef, index=genes)

    elif method == "stable":
        X = expr(adata, train, genes)
        rng = np.random.default_rng(seed)
        idx1, idx0 = np.where(y == 1)[0], np.where(y == 0)[0]
        up = np.zeros(len(genes)); down = np.zeros(len(genes)); auc_sum = np.zeros(len(genes))
        for _ in range(n_boot):
            s = np.concatenate([rng.choice(idx1, max(2, int(sub_frac * len(idx1))), replace=False),
                                rng.choice(idx0, max(2, int(sub_frac * len(idx0))), replace=False)])
            auc = auc_scores(X[s], y[s])
            auc_sum += auc
            up[np.argsort(-auc)[:n_per_side]] += 1
            down[np.argsort(auc)[:n_per_side]] += 1
        mean_auc = auc_sum / n_boot
        g = np.asarray(genes)
        # most frequently selected first; ties broken by mean AUC
        o_up = np.lexsort((-mean_auc, -up))[:n_per_side]
        o_dn = np.lexsort((mean_auc, -down))[:n_per_side]
        out = list(g[o_up]) + list(g[o_dn])
        scores = pd.DataFrame({"freq_cancer": up / n_boot, "freq_normal": down / n_boot,
                               "mean_auc": mean_auc}, index=genes)
    else:
        raise ValueError(f"unknown method {method!r}")

    return (out, scores) if return_scores else out


# ---------------------------------------------------------------------------
# 4. The QUBO
# ---------------------------------------------------------------------------

def relevance(X, y, kind="corr", seed=DEFAULT_SEED):
    """How strongly each gene follows the label.
    "corr": |correlation with the label| (what I used so far)
    "auc":  |AUC - 0.5| * 2  (rank-based, robust to outliers and zeros)
    "mi":   mutual information (also catches nonlinear relationships)"""
    if kind == "corr":
        Xc = X - X.mean(0)
        yc = y - y.mean()
        denom = np.sqrt((Xc**2).sum(0) * (yc**2).sum()) + 1e-12
        return np.abs(Xc.T @ yc) / denom
    if kind == "auc":
        return np.abs(auc_scores(X, y) - 0.5) * 2
    if kind == "mi":
        from sklearn.feature_selection import mutual_info_classif
        return mutual_info_classif(X, y, random_state=seed)
    raise ValueError(f"unknown relevance {kind!r}")


def redundancy(X, kind="pos_corr"):
    """Pair penalty matrix (zero diagonal).
    "pos_corr": only positive correlation counts (opposite genes are complementary)
    "abs_corr": |correlation| (my first version)
    "none":     no redundancy term"""
    n = X.shape[1]
    if kind == "none":
        return np.zeros((n, n))
    C = np.nan_to_num(np.corrcoef(X.T))
    R = np.clip(C, 0, None) if kind == "pos_corr" else np.abs(C)
    np.fill_diagonal(R, 0)
    return R


def safe_lambda(r, R, k, alpha, margin=1.05):
    """Smallest size penalty that provably keeps the optimum at exactly k genes.

    Adding d extra genes can gain at most d * max(r) (redundancy only adds cost).
    Removing d genes can gain at most d * alpha * (k - 1) * max(R).
    The penalty grows like lambda * d^2 >= lambda * d, so any lambda above both
    bounds makes every wrong-size panel worse than some right-size panel.
    A smaller lambda keeps the penalty from drowning out the gene choice in QAOA."""
    bound = max(float(np.max(r)), alpha * (k - 1) * float(np.max(R)) if R.size else 0.0)
    return margin * bound + 1e-9


def qubo_from_scores(r, R, k=2, alpha=1.0, lam="auto"):
    """Q and const for:  - sum r_i x_i + alpha * sum_{i<j} R_ij x_i x_j + lam * (sum x_i - k)^2
    (symmetric convention: each pair term is split between Q[i,j] and Q[j,i])."""
    n = len(r)
    lam = safe_lambda(r, R, k, alpha) if lam == "auto" else float(lam)
    Q = alpha * R / 2 + lam * (np.ones((n, n)) - np.eye(n))
    Q[np.diag_indices(n)] += -np.asarray(r) + lam * (1 - 2 * k)
    return Q, lam * k**2, lam


def build_qubo(adata, train, cands, k=2, alpha=1.0, lam="auto",
               relevance_kind="corr", redundancy_kind="pos_corr"):
    """QUBO for the candidate genes, from the training cells.
    Returns Q, const, and an info dict (r, R, lam)."""
    X = expr(adata, train, cands)
    y = labels(adata, train)
    r = relevance(X, y, relevance_kind)
    R = redundancy(X, redundancy_kind)
    Q, const, lam_used = qubo_from_scores(r, R, k, alpha, lam)
    return Q, const, {"r": r, "R": R, "lam": lam_used, "alpha": alpha, "k": k}


# ---------------------------------------------------------------------------
# 5. Exact and classical solvers
# ---------------------------------------------------------------------------

def bitstrings(n):
    """All 2^n bit strings; row b is the binary of b with bit i = gene i."""
    return ((np.arange(2**n)[:, None] >> np.arange(n)) & 1).astype(np.int8)


def all_energies(Q, const, chunk=1 << 18):
    """QUBO score x^T Q x + const for every bit string (index b = panel of b)."""
    n = Q.shape[0]
    out = np.empty(2**n)
    for start in range(0, 2**n, chunk):
        b = np.arange(start, min(start + chunk, 2**n))
        X = ((b[:, None] >> np.arange(n)) & 1).astype(float)
        out[start:start + len(b)] = np.einsum("bi,ij,bj->b", X, Q, X) + const
    return out


def panel_of(b, cands):
    return [cands[i] for i in range(len(cands)) if (b >> i) & 1]


def index_of(panel, cands):
    return int(sum(1 << cands.index(g) for g in panel))


def brute_force(Q, const, cands):
    """Exact optimum over all 2^n panels. Returns (panel, score, index)."""
    E = all_energies(Q, const)
    b = int(np.argmin(E))
    return panel_of(b, cands), float(E[b]), b


def best_panel(Q, const, cands, k, max_subsets=5_000_000):
    """Exact best panel of exactly k genes, by enumerating all k-subsets
    (much cheaper than 2^n when k is small: 16 choose 2 = 120)."""
    n = len(cands)
    if math.comb(n, k) > max_subsets:
        raise ValueError("too many subsets; use simulated_annealing")
    best, best_e = None, np.inf
    combos = itertools.combinations(range(n), k)
    while True:
        block = np.array(list(itertools.islice(combos, 200_000)))
        if len(block) == 0:
            break
        e = Q[block[:, :, None], block[:, None, :]].sum(axis=(1, 2)) + const
        j = int(np.argmin(e))
        if e[j] < best_e:
            best_e, best = e[j], block[j]
    return [cands[i] for i in best]


def simulated_annealing(Q, const, cands, k, n_steps=20000, n_restarts=10,
                        t_start=None, t_end=1e-3, seed=DEFAULT_SEED):
    """Classical heuristic for large problems: simulated annealing that always keeps
    exactly k genes (each move swaps one gene in for one gene out).
    Kirkpatrick, Gelatt & Vecchi 1983. Returns (panel, score)."""
    rng = np.random.default_rng(seed)
    n = len(cands)
    Qs = (Q + Q.T) / 2
    if t_start is None:
        t_start = float(np.abs(Qs).mean() * k) + 1e-9
    best_x, best_e = None, np.inf
    for _ in range(n_restarts):
        x = np.zeros(n, dtype=bool)
        x[rng.choice(n, k, replace=False)] = True
        e = float(x @ Qs @ x) + const
        Qx = Qs @ x
        for step in range(n_steps):
            t = t_start * (t_end / t_start) ** (step / max(1, n_steps - 1))
            i = rng.choice(np.where(x)[0])        # leaves
            j = rng.choice(np.where(~x)[0])       # joins
            d = (-2 * Qx[i] + Qs[i, i]) + (2 * (Qx[j] - Qs[j, i]) + Qs[j, j])
            if d <= 0 or rng.random() < math.exp(-d / t):
                x[i], x[j] = False, True
                Qx += Qs[:, j] - Qs[:, i]
                e += d
                if e < best_e:
                    best_e, best_x = e, x.copy()
        if e < best_e:
            best_e, best_x = e, x.copy()
    return [cands[i] for i in np.where(best_x)[0]], float(best_e)


# ---------------------------------------------------------------------------
# 6. Grading panels
# ---------------------------------------------------------------------------

def _fit_score(Xtr, ytr, Xte, yte):
    # class_weight="balanced": both classes count equally during training,
    # matching the balanced-accuracy score used for grading
    model = LogisticRegression(max_iter=1000, class_weight="balanced").fit(Xtr, ytr)
    return balanced_accuracy_score(yte, model.predict(Xte))


def evaluate(adata, panel, train, test):
    """Balanced accuracy of class-balanced logistic regression on the panel, trained on train, graded on test."""
    return _fit_score(expr(adata, train, panel), labels(adata, train),
                      expr(adata, test, panel), labels(adata, test))


def repeated_splits(adata, panel, n_splits=10, test_size=0.5):
    """Mean and std of balanced accuracy over n_splits random splits (seeds 0..n-1)."""
    y = labels(adata)
    X = expr(adata, adata.obs_names, panel)
    scores = []
    for seed in range(n_splits):
        tr, te = train_test_split(np.arange(len(y)), test_size=test_size,
                                  random_state=seed, stratify=y)
        scores.append(_fit_score(X[tr], y[tr], X[te], y[te]))
    return float(np.mean(scores)), float(np.std(scores))


def choose_alpha_cv(adata, train, cands, k=2, alphas=(0, 0.05, 0.1, 0.15, 0.25, 0.5, 1.0),
                    folds=5, relevance_kind="corr", redundancy_kind="pos_corr", seed=DEFAULT_SEED):
    """Pick the redundancy weight alpha by cross-validation INSIDE the training cells.
    For each alpha and each fold: build the QUBO on the fold's training part, take its
    best k-gene panel, grade it on the fold's held-out part. The test half is never used.
    Returns (best_alpha, table)."""
    Xall = expr(adata, train, cands)
    yall = labels(adata, train)
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    rows = []
    for a in alphas:
        scores, panels = [], []
        for tr, va in skf.split(Xall, yall):
            r = relevance(Xall[tr], yall[tr], relevance_kind)
            R = redundancy(Xall[tr], redundancy_kind)
            Q, const, _ = qubo_from_scores(r, R, k, a, "auto")
            panel = best_panel(Q, const, list(range(len(cands))), k)
            panels.append(tuple(cands[i] for i in panel))
            scores.append(_fit_score(Xall[tr][:, panel], yall[tr], Xall[va][:, panel], yall[va]))
        most = max(set(panels), key=panels.count)
        rows.append({"alpha": a, "cv balanced accuracy": np.mean(scores),
                     "std": np.std(scores), "most common panel": ", ".join(most)})
    table = pd.DataFrame(rows)
    best = float(table.loc[table["cv balanced accuracy"].idxmax(), "alpha"])
    return best, table


def lasso_panel(adata, train, k=2, min_cells=10):
    """Classical ML rival: the lasso picks k genes by itself from all detected genes."""
    genes = detected_genes(adata, train, min_cells)
    X = StandardScaler().fit_transform(expr(adata, train, genes))
    y = labels(adata, train)
    w = None
    for C in np.logspace(-4, 0, 80):
        w = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=5000).fit(X, y).coef_[0]
        if (w != 0).sum() >= k:
            break
    return list(np.asarray(genes)[np.argsort(-np.abs(w))[:k]])


# ---------------------------------------------------------------------------
# 7. Quantum
# ---------------------------------------------------------------------------

def hamiltonian(Q, const):
    """Cost Hamiltonian (Qiskit SparsePauliOp) with x = (1 - Z)/2, plus the constant offset.
    Its diagonal + offset equals all_energies(Q, const)."""
    from qiskit.quantum_info import SparsePauliOp
    n = Q.shape[0]
    terms, offset = [], const
    for i in range(n):
        terms.append(("Z", [i], -Q[i, i] / 2))
        offset += Q[i, i] / 2
        for j in range(n):
            if i != j:
                terms.append(("ZZ", [i, j], Q[i, j] / 4))
                terms.append(("Z", [i], -Q[i, j] / 4))
                terms.append(("Z", [j], -Q[i, j] / 4))
                offset += Q[i, j] / 4
    return SparsePauliOp.from_sparse_list(terms, num_qubits=n).simplify(), offset


class QAOASim:
    """Fast exact statevector simulator for QAOA on a diagonal cost (NumPy only).

    cost layer:  psi <- exp(-i * gamma * E) * psi          (E = energy of each bit string)
    X mixer:     exp(-i * beta * sum_i X_i)  = Rx(2 beta) on every qubit (standard QAOA)
    XY mixer:    exp(-i * beta * (X_i X_j + Y_i Y_j) / 2) on a ring of pairs
                 (even pairs, then odd pairs). It swaps 01 <-> 10, so the number of 1s
                 never changes: start with exactly k genes, stay with exactly k genes.
                 This is the "quantum alternating operator ansatz" (Hadfield et al. 2019).
    init:        "plus"  = uniform superposition over all bit strings
                 "dicke" = uniform superposition over bit strings with exactly k ones
    Matches Qiskit's QAOAAnsatz for the X mixer (see qiskit_check)."""

    def __init__(self, energies, mixer="x", init="plus", k=None):
        self.E = np.asarray(energies, dtype=float)
        self.N = len(self.E)
        self.n = int(round(math.log2(self.N)))
        self.mixer, self.init_kind, self.k = mixer, init, k
        self.weights = np.array([bin(b).count("1") for b in range(self.N)])
        if mixer == "xy":
            n = self.n
            pairs = [(i, i + 1) for i in range(0, n - 1, 2)] + [(i, i + 1) for i in range(1, n - 1, 2)]
            if n > 2:
                pairs.append((n - 1, 0))
            b = np.arange(self.N)
            self.xy = []
            for i, j in pairs:
                idx01 = b[(((b >> i) & 1) == 0) & (((b >> j) & 1) == 1)]
                self.xy.append((idx01, idx01 ^ (1 << i) ^ (1 << j)))

    def initial_state(self):
        if self.init_kind == "plus":
            return np.full(self.N, 1 / math.sqrt(self.N), dtype=complex)
        if self.init_kind == "dicke":
            psi = np.zeros(self.N, dtype=complex)
            mask = self.weights == self.k
            psi[mask] = 1 / math.sqrt(mask.sum())
            return psi
        raise ValueError(self.init_kind)

    def _mix_x(self, psi, beta):
        c, s = math.cos(beta), -1j * math.sin(beta)
        for q in range(self.n):
            v = psi.reshape(-1, 2, 2**q)
            a0 = v[:, 0, :].copy()
            v[:, 0, :] *= c
            v[:, 0, :] += s * v[:, 1, :]
            v[:, 1, :] *= c
            v[:, 1, :] += s * a0
        return psi

    def _mix_xy(self, psi, beta):
        c, s = math.cos(beta), -1j * math.sin(beta)
        for i01, i10 in self.xy:
            a, b = psi[i01].copy(), psi[i10].copy()
            psi[i01] = c * a + s * b
            psi[i10] = s * a + c * b
        return psi

    def state(self, gammas, betas):
        psi = self.initial_state()
        for g, b in zip(gammas, betas):
            psi *= np.exp(-1j * g * self.E)
            psi = self._mix_x(psi, b) if self.mixer == "x" else self._mix_xy(psi, b)
        return psi

    def probs(self, gammas, betas):
        return np.abs(self.state(gammas, betas)) ** 2

    def expectation(self, gammas, betas):
        return float(self.probs(gammas, betas) @ self.E)


def interp_angles(gammas, betas):
    """INTERP initialization (Zhou, Wang, Choi, Pichler & Lukin, PRX 2020):
    turn good p-layer angles into a smooth starting guess for p + 1 layers."""
    def up(a):
        p = len(a)
        a = np.concatenate([[0.0], a, [0.0]])
        return np.array([(i / p) * a[i] + ((p - i) / p) * a[i + 1] for i in range(p + 1)])
    return up(np.asarray(gammas)), up(np.asarray(betas))


def qaoa_metrics(probs, energies, k=None):
    """avg energy, approximation ratio, P(best), P(right size), shots needed for 99% chance
    of seeing the best panel at least once.
    Approximation ratio here = (E_max - <E>) / (E_max - E_min): 1 = perfect, 0 = worst."""
    E = np.asarray(energies)
    best = int(np.argmin(E))
    avg = float(probs @ E)
    out = {"avg energy": avg,
           "approx ratio": (E.max() - avg) / (E.max() - E.min()),
           "P(best panel)": float(probs[best])}
    if k is not None:
        w = np.array([bin(b).count("1") for b in range(len(E))])
        out["P(k genes)"] = float(probs[w == k].sum())
    p = out["P(best panel)"]
    out["shots for 99%"] = math.ceil(math.log(0.01) / math.log(1 - p)) if 0 < p < 1 else (1 if p >= 1 else np.inf)
    return out


def optimize_qaoa(energies, p, mixer="x", init="plus", k=None, starts=6, init_angles=None,
                  maxiter=500, seed=DEFAULT_SEED, normalize=True):
    """Find good QAOA angles: COBYLA from several starting points, keep the best.

    The energies are divided by their spread first (normalize=True), so the same angle
    ranges work for every dataset and every lambda; the returned gammas are converted
    back to the original energy units. init_angles = (gammas, betas) in original units
    (for example from interp_angles) is tried first."""
    E = np.asarray(energies, dtype=float)
    scale = float(E.std()) if normalize and E.std() > 0 else 1.0
    sim = QAOASim(E / scale, mixer, init, k)
    rng = np.random.default_rng(seed)

    def f(theta):
        return sim.expectation(theta[:p], theta[p:])

    trials = []
    if init_angles is not None:
        g0, b0 = init_angles
        trials.append(np.concatenate([np.asarray(g0) * scale, b0]))
    while len(trials) < starts:
        trials.append(np.concatenate([rng.uniform(0, np.pi, p), rng.uniform(0, np.pi / 2, p)]))

    best = None
    for t0 in trials:
        res = minimize(f, t0, method="COBYLA", options={"maxiter": maxiter})
        if best is None or res.fun < best.fun:
            best = res
    gam, bet = best.x[:p] / scale, best.x[p:]
    probs = QAOASim(E, mixer, init, k).probs(gam, bet)
    return {"p": p, "gammas": gam, "betas": bet, "probs": probs, **qaoa_metrics(probs, E, k)}


def qaoa_ladder(energies, p_max=4, mixer="x", init="plus", k=None, starts=6, seed=DEFAULT_SEED):
    """Run p = 1 .. p_max. Each p starts from INTERP of the previous best angles plus
    random restarts. Returns (table, results list)."""
    results, prev = [], None
    for p in range(1, p_max + 1):
        t = time.time()
        init_angles = interp_angles(prev["gammas"], prev["betas"]) if prev else None
        res = optimize_qaoa(energies, p, mixer, init, k, starts, init_angles, seed=seed + p)
        res["seconds"] = time.time() - t
        results.append(res)
        prev = res
    cols = ["p", "avg energy", "approx ratio", "P(best panel)", "P(k genes)", "shots for 99%", "seconds"]
    table = pd.DataFrame([{c: r.get(c) for c in cols} for r in results]).set_index("p")
    return table, results


def baseline_row(energies, k=None, init="plus"):
    """Metrics of the starting state alone (no QAOA layers): the 'random guess' row."""
    sim = QAOASim(energies, "x", init, k)
    return qaoa_metrics(np.abs(sim.initial_state()) ** 2, energies, k)


def p1_landscape(energies, gammas, betas, mixer="x", init="plus", k=None):
    """Average energy on a grid of (gamma, beta) for p = 1 (rows = gammas)."""
    sim = QAOASim(energies, mixer, init, k)
    return np.array([[sim.expectation([g], [b]) for b in betas] for g in gammas])


def sample(probs, shots=1000, seed=DEFAULT_SEED):
    """Simulate measuring the circuit `shots` times. Returns counts {index: count}."""
    rng = np.random.default_rng(seed)
    draws = rng.choice(len(probs), size=shots, p=probs / probs.sum())
    idx, cnt = np.unique(draws, return_counts=True)
    return dict(zip(idx.tolist(), cnt.tolist()))


def qiskit_check(Q, const, gamma=0.3, beta=0.4):
    """Compare this simulator with Qiskit's QAOAAnsatz (p = 1, X mixer).
    The two numbers should agree to many decimals."""
    from qiskit.circuit.library import QAOAAnsatz
    from qiskit.quantum_info import Statevector
    H, offset = hamiltonian(Q, const)
    qc = QAOAAnsatz(H, reps=1)
    e_qiskit = float(np.real(Statevector(qc.assign_parameters([beta, gamma])).expectation_value(H))) + offset
    e_numpy = QAOASim(all_energies(Q, const)).expectation([gamma], [beta])
    return e_qiskit, e_numpy


# ---------------------------------------------------------------------------
# 8. One-call pipeline
# ---------------------------------------------------------------------------

def run_classical(adata, method="stable", k=2, alpha="cv", lam="auto",
                  relevance_kind="corr", redundancy_kind="pos_corr", seed=DEFAULT_SEED,
                  n_per_side=8):
    """Whole classical part for one dataset. Returns a dict with candidates, Q, panels
    and a results table (seed-42 test score and 10-split mean +- std)."""
    adata = prepare(adata)
    train, test = split(adata, seed)
    cands = candidates(adata, train, method, n_per_side=n_per_side, seed=seed)
    cv_table = None
    if alpha == "cv":
        alpha, cv_table = choose_alpha_cv(adata, train, cands, k,
                                          relevance_kind=relevance_kind,
                                          redundancy_kind=redundancy_kind, seed=seed)
    Q, const, info = build_qubo(adata, train, cands, k, alpha, lam, relevance_kind, redundancy_kind)
    order = np.argsort(-info["r"])
    panels = {
        "QUBO": best_panel(Q, const, cands, k),
        f"top {k} by relevance": [cands[i] for i in order[:k]],
        "best single gene": [cands[order[0]]],
        f"lasso ({k} genes)": lasso_panel(adata, train, k),
    }
    rows = []
    for name, p in panels.items():
        m, s = repeated_splits(adata, p)
        rows.append({"method": name, "panel": ", ".join(p),
                     "test (seed 42)": round(evaluate(adata, p, train, test), 3),
                     "10 splits": f"{m:.3f} ± {s:.3f}"})
    return {"adata": adata, "train": train, "test": test, "candidates": cands,
            "alpha": alpha, "cv_table": cv_table, "Q": Q, "const": const, "info": info,
            "panels": panels, "table": pd.DataFrame(rows).set_index("method")}
