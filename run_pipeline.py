"""Reproduce every result table and figure."""
import pandas as pd

from src import braak, config, plots, spatial_test
from src.fetch_expression import fetch_expression, region_info
from src.fetch_gwas import analysis_genes, load_genes
from src.score_regions import leave_one_locus_out, score, without_complex_loci

FOOTER = f"Allen Human Brain Atlas (6 donors) · Desikan-Killiany atlas · GWAS Catalog {config.TRAIT_ID}"


def main():
    expression = fetch_expression()
    info = region_info()

    # Positive control: striatal markers must rank caudate and putamen at the top.
    control, _, control_genes = score(expression, config.POSITIVE_CONTROL)
    control = info.join(control, how="right")
    top = control.sort_values("z", ascending=False)["label"].head(4).tolist()
    assert {"caudate", "putamen"} <= set(top), f"positive control failed: top regions {top}"
    control.sort_values("z", ascending=False).to_csv(config.TABLES / "positive_control_scores.csv")

    # PD risk genes vs. 10,000 random gene sets in every region
    gene_table = load_genes(expression.columns)
    genes = analysis_genes(gene_table)
    results, null, available = score(expression, genes)
    results = info.join(results, how="right")

    # Sensitivity: without the complex loci, and leaving out each locus in turn
    results["z_without_complex_loci"] = without_complex_loci(expression, genes)["z"]
    loo = leave_one_locus_out(expression, gene_table)
    results["z_loo_min"], results["z_loo_max"] = loo.min(axis=1), loo.max(axis=1)
    loo.to_csv(config.TABLES / "leave_one_locus_out_z.csv")
    results.sort_values("z", ascending=False).to_csv(config.REGION_SCORES)
    print(f"regions with FDR q < {config.FDR_ALPHA}: {(results['q'] < config.FDR_ALPHA).sum()} of {len(results)}")

    # Braak staging trend
    trend = braak.braak_trend(results, null, expression.index)
    trend["table"].to_csv(config.TABLES / "braak_stage_regions.csv")
    pd.DataFrame([{"rho": trend["rho"], "p_decreasing": trend["p"], "n_regions": len(trend["table"])}]) \
        .to_csv(config.TABLES / "braak_trend.csv", index=False)
    print(f"Braak trend: rho = {trend['rho']:.3f}, p = {trend['p']:.4f}")

    # Expression map vs. ENIGMA-PD cortical thinning, naive and spin nulls
    atrophy = spatial_test.load_atrophy(info)
    spatial_summary, spatial = spatial_test.compare(results["z"], atrophy)
    spatial_summary.to_csv(config.TABLES / "spatial_test.csv", index=False)
    print(spatial_summary.to_string(index=False))

    # Figures
    plots.braak_plot(trend)
    plots.spatial_null_plot(spatial, spatial_summary)
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
