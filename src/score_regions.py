"""Gene-set permutation test: is a gene set more highly expressed in a region than random gene sets?"""
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from src import config


def gene_set_null(expression, k, n_perm=config.N_PERM_GENES, seed=config.SEED, chunk=500):
    """Mean expression of n_perm random k-gene sets in every region -> (n_perm, n_regions)."""
    rng = np.random.default_rng(seed)
    X = expression.to_numpy(dtype=np.float32)          # regions x genes
    n_genes = X.shape[1]
    null = np.empty((n_perm, X.shape[0]), dtype=np.float32)
    for start in range(0, n_perm, chunk):
        stop = min(start + chunk, n_perm)
        # k distinct genes per random set, encoded as a 0/1 mask so the means are one matrix product
        picks = rng.random((stop - start, n_genes)).argpartition(k, axis=1)[:, :k]
        mask = np.zeros((stop - start, n_genes), dtype=np.float32)
        np.put_along_axis(mask, picks, 1.0, axis=1)
        null[start:stop] = mask @ X.T / k
    return null


def score(expression, genes, n_perm=config.N_PERM_GENES, seed=config.SEED):
    """Per-region observed mean, z, one-sided permutation p and BH-FDR q."""
    expression = expression.dropna(axis=0, how="all").dropna(axis=1, how="any")
    available = [g for g in genes if g in expression.columns]

    observed = expression[available].mean(axis=1).to_numpy()
    null = gene_set_null(expression, len(available), n_perm, seed)

    z = (observed - null.mean(axis=0)) / null.std(axis=0)
    p = (1 + (null >= observed).sum(axis=0)) / (n_perm + 1)
    q = multipletests(p, method="fdr_bh")[1]

    results = pd.DataFrame({"observed": observed, "z": z, "p": p, "q": q}, index=expression.index)
    results.index.name = "id"
    return results, null, available


def without_complex_loci(expression, genes):
    """Re-score without the MHC / MAPT / LRRK2 genes."""
    reduced = [g for g in genes if g not in config.COMPLEX_LOCI_GENES]
    return score(expression, reduced)[0]


def leave_one_locus_out(expression, gene_table, n_perm=1_000):
    """z for every region with each GWAS locus removed in turn -> (regions x loci).

    z only needs the null mean and SD, which 1,000 random sets estimate well.
    """
    table = gene_table[(gene_table["locus_group"] == "protein-coding gene")
                       & gene_table["expression_symbol"].notna()].drop_duplicates("expression_symbol")
    genes = table["expression_symbol"].tolist()
    runs = {}
    for locus, dropped in table.groupby("source_locus")["expression_symbol"]:
        kept = [g for g in genes if g not in set(dropped)]
        runs[locus] = score(expression, kept, n_perm=n_perm)[0]["z"]
    return pd.DataFrame(runs)
