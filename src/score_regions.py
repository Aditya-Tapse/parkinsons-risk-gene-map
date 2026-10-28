"""Step 3: is the PD gene set unusually expressed in each region, relative to random gene sets?"""
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from src import config


def gene_set_null(expression, k, n_perm=config.N_PERM_GENES, seed=config.SEED, chunk=500):
    """Mean expression of n_perm random k-gene sets in every region -> (n_perm, n_regions).

    Same idea as the PDF's loop, but done as a matrix product in chunks so 10k
    permutations take seconds instead of minutes.
    """
    rng = np.random.default_rng(seed)
    X = expression.to_numpy(dtype=np.float32)          # regions x genes
    n_genes = X.shape[1]
    null = np.empty((n_perm, X.shape[0]), dtype=np.float32)
    for start in range(0, n_perm, chunk):
        stop = min(start + chunk, n_perm)
        # argpartition of uniform noise = k distinct genes per row, without replacement
        picks = rng.random((stop - start, n_genes)).argpartition(k, axis=1)[:, :k]
        mask = np.zeros((stop - start, n_genes), dtype=np.float32)
        np.put_along_axis(mask, picks, 1.0, axis=1)
        null[start:stop] = mask @ X.T / k
    return null


def score(expression, genes, n_perm=config.N_PERM_GENES, seed=config.SEED, verbose=True):
    expression = expression.dropna(axis=0, how="all").dropna(axis=1, how="any")
    available = [g for g in genes if g in expression.columns]
    if verbose:
        print(f"{len(available)} of {len(genes)} genes found in expression data")

    observed = expression[available].mean(axis=1).to_numpy()
    null = gene_set_null(expression, len(available), n_perm, seed)

    z = (observed - null.mean(axis=0)) / null.std(axis=0)
    p = (1 + (null >= observed).sum(axis=0)) / (n_perm + 1)
    _, q, _, _ = multipletests(p, method="fdr_bh")

    results = pd.DataFrame(
        {"observed": observed, "z": z, "p": p, "q": q}, index=expression.index
    )
    results.index.name = "id"
    return results, null, available


def without_complex_loci(expression, genes):
    """Re-score without the MHC / MAPT / LRRK2 genes, which bring many correlated genes per signal."""
    reduced = [g for g in genes if g not in config.COMPLEX_LOCI_GENES]
    results, _, _ = score(expression, reduced, verbose=False)
    print(f"complex loci removed: {len(genes) - len(reduced)} genes dropped, {len(reduced)} left")
    return results


def leave_one_locus_out(expression, gene_table, n_perm=1_000):
    """z for every region with each GWAS locus removed in turn -> (regions x loci).

    Fewer permutations than the main run: only the null's mean and SD are
    needed for z, and 1,000 sets pin those down well.
    """
    table = gene_table.dropna(subset=["expression_symbol"]).drop_duplicates("expression_symbol")
    table = table[table["locus_group"] == "protein-coding gene"]
    genes = table["expression_symbol"].tolist()
    runs = {}
    for locus, dropped in table.groupby("source_locus")["expression_symbol"]:
        kept = [g for g in genes if g not in set(dropped)]
        results, _, _ = score(expression, kept, n_perm=n_perm, verbose=False)
        runs[locus] = results["z"]
    print(f"leave-one-locus-out: {len(runs)} loci")
    return pd.DataFrame(runs)


if __name__ == "__main__":
    from src.fetch_expression import fetch_expression, region_info
    from src.fetch_gwas import analysis_genes

    expression = fetch_expression()
    info = region_info()
    gene_table = pd.read_csv(config.GENES_CSV)
    genes = analysis_genes(gene_table)
    results, null, _ = score(expression, genes)
    results = info.join(results, how="right")
    np.save(config.NULL_NPY, null)

    results["z_without_complex_loci"] = without_complex_loci(expression, genes)["z"]
    loo = leave_one_locus_out(expression, gene_table)
    results["z_loo_min"], results["z_loo_max"] = loo.min(axis=1), loo.max(axis=1)
    results.sort_values("z", ascending=False).to_csv(config.REGION_SCORES)
    loo.to_csv(config.TABLES / "leave_one_locus_out_z.csv")
    print(results.sort_values("z", ascending=False).head(10)
          [["name", "z", "q", "z_without_complex_loci", "z_loo_min", "z_loo_max"]].round(2))
