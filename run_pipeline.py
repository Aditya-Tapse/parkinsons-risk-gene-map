"""Reproduce every result and figure in the README. Cached steps are skipped on re-runs.

Steps 6 (Braak staging) and 8 (atrophy + spin test) are drafted in
src/braak.py and src/spatial_test.py but not yet part of the pipeline.
"""
import numpy as np

from src import config, plots
from src.fetch_expression import fetch_expression, region_info
from src.fetch_gwas import analysis_genes, load_genes
from src.score_regions import leave_one_locus_out, score, without_complex_loci

FOOTER = f"Allen Human Brain Atlas (6 donors) · Desikan-Killiany atlas · GWAS Catalog {config.TRAIT_ID}"


def main():
    # Step 2: expression data (83 regions x ~15.6k genes)
    expression = fetch_expression()
    info = region_info()

    # Step 4: positive control. If striatal markers don't light up the striatum, stop.
    control, _, control_genes = score(expression, config.POSITIVE_CONTROL)
    control = info.join(control, how="right")
    top = control.sort_values("z", ascending=False)["label"].head(4).tolist()
    print(f"positive control top regions: {top}")
    assert {"caudate", "putamen"} <= set(top), "pipeline failed the striatal positive control"
    control.sort_values("z", ascending=False).to_csv(config.TABLES / "positive_control_scores.csv")

    # Step 3: PD risk genes
    gene_table = load_genes(expression.columns)
    genes = analysis_genes(gene_table)

    # Step 5: score every region against 10,000 random gene sets
    results, null, available = score(expression, genes)
    results = info.join(results, how="right")
    np.save(config.NULL_NPY, null)

    # Sensitivity: drop the complex MHC / MAPT / LRRK2 loci, then each locus in turn.
    results["z_without_complex_loci"] = without_complex_loci(expression, genes)["z"]
    loo = leave_one_locus_out(expression, gene_table)
    results["z_loo_min"], results["z_loo_max"] = loo.min(axis=1), loo.max(axis=1)
    loo.to_csv(config.TABLES / "leave_one_locus_out_z.csv")
    results.sort_values("z", ascending=False).to_csv(config.REGION_SCORES)
    print(f"regions passing FDR q < {config.FDR_ALPHA}: {(results['q'] < config.FDR_ALPHA).sum()} of {len(results)}")

    # Figures
    plots.hero(
        results,
        title="Parkinson's risk genes are not more active where the disease strikes",
        subtitle=f"Mean expression of {len(available)} GWAS risk genes in each brain region, vs. 10,000 random gene sets",
        footer=FOOTER,
    )
    plots.hero(
        control,
        title="Positive control: the method finds the striatum when it should",
        subtitle=f"Mean expression of {len(control_genes)} known striatal marker genes "
                 f"({', '.join(control_genes)}), vs. 10,000 random gene sets",
        footer="Same data and method as the Parkinson's analysis",
        out_name="positive_control", vmax=6.0,
    )
    best = results["z"].idxmax()
    plots.null_histogram(
        null[:, results.index.get_loc(best)], results.loc[best, "observed"],
        f"{plots.pretty(results.loc[best, 'label'])} ({results.loc[best, 'hemisphere']})",
    )


if __name__ == "__main__":
    main()
