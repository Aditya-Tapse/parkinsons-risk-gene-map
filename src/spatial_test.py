"""Map-to-map test: does the PD-gene expression map correlate with PD cortical thinning?

Neighbouring regions have similar values (spatial autocorrelation), so shuffling
region labels gives p-values that are too small. The spin test rotates the map
on the sphere instead, which preserves its spatial structure.
"""
import numpy as np
import pandas as pd
from neuromaps import nulls
from nilearn import datasets
from scipy.stats import spearmanr

from src import config
from src.fetch_expression import load_atlas


def load_atrophy(info):
    """ENIGMA-PD cortical thickness Cohen's d (Laansma et al. 2021, Suppl. Table 2a), by region id."""
    atrophy = pd.read_csv(config.ATROPHY_CSV)
    atrophy["name"] = atrophy["label"] + "_" + atrophy["hemisphere"]
    cortex = info[info["structure"] == "cortex"]
    return cortex.reset_index().merge(atrophy[["name", "d"]], on="name").set_index("id")["d"]


def naive_null(x, y, n_perm=config.N_PERM_GENES, seed=config.SEED):
    rng = np.random.default_rng(seed)
    return np.array([spearmanr(rng.permutation(x), y).statistic for _ in range(n_perm)])


def spin_null(x, n_perm=config.N_PERM_SPATIAL, seed=config.SEED):
    """Alexander-Bloch spins of a DK cortical map -> (n_regions, n_perm).

    x must be sorted by region id: neuromaps relabels the parcels (lh 1-34,
    rh 42-75) in sorted order.
    """
    fsavg = datasets.fetch_surf_fsaverage("fsaverage5")
    return nulls.alexander_bloch(
        x, atlas="fsaverage", density="10k", parcellation=load_atlas(surface=True)["image"],
        n_perm=n_perm, seed=seed, surfaces=(fsavg["sphere_left"], fsavg["sphere_right"]),
    )


def compare(expression_z, atrophy_d):
    """Spearman rho between the two maps with naive and spin two-sided p-values."""
    ids = expression_z.index.intersection(atrophy_d.index).sort_values()
    x, y = expression_z.loc[ids].to_numpy(), atrophy_d.loc[ids].to_numpy()
    rho = spearmanr(x, y).statistic

    naive = naive_null(x, y)
    spins = spin_null(x)
    spun = np.array([spearmanr(spins[:, i], y, nan_policy="omit").statistic for i in range(spins.shape[1])])

    def two_sided(null):
        return (1 + (np.abs(null) >= abs(rho)).sum()) / (len(null) + 1)

    summary = pd.DataFrame({
        "null_model": ["naive (shuffle regions)", "spatial (spin test)"],
        "rho": [rho, rho],
        "p": [two_sided(naive), two_sided(spun)],
        "n_perm": [len(naive), len(spun)],
    })
    return summary, {"x": x, "y": y, "naive": naive, "spin": spun}
