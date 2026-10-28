"""Step 5: does the PD-gene expression map line up with where PD thins the cortex?

Neighbouring regions have similar values (spatial autocorrelation), so shuffling
region labels gives p-values that are too small. A spin test rotates the map
on the sphere instead, which keeps its spatial structure intact.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from src import config


def load_atrophy(info):
    """ENIGMA-PD cortical thickness effect sizes (Laansma et al. 2021), one per DK cortical region.

    Expected CSV columns: label, hemisphere (L/R), d. Values must be copied from
    the paper's supplementary table; nothing here is filled in automatically.
    """
    if not config.ATROPHY_CSV.exists():
        raise FileNotFoundError(f"add the ENIGMA-PD effect sizes to {config.ATROPHY_CSV}")
    atrophy = pd.read_csv(config.ATROPHY_CSV)
    atrophy["name"] = atrophy["label"] + "_" + atrophy["hemisphere"]
    cortex = info[info["structure"] == "cortex"]
    return cortex.reset_index().merge(atrophy[["name", "d"]], on="name").set_index("id")["d"]


def naive_null(x, y, n_perm=config.N_PERM_GENES, seed=config.SEED):
    rng = np.random.default_rng(seed)
    return np.array([spearmanr(rng.permutation(x), y).statistic for _ in range(n_perm)])


def spin_null(x, n_perm=config.N_PERM_SPATIAL, seed=config.SEED):
    """Alexander-Bloch spins of a DK cortical map; returns (n_regions, n_perm)."""
    from neuromaps import nulls
    from src.fetch_expression import load_atlas

    parcellation = load_atlas(surface=True)["image"]   # (lh, rh) fsaverage5 label GIFTIs
    return nulls.alexander_bloch(
        x, atlas="fsaverage", density="10k", parcellation=parcellation,
        n_perm=n_perm, seed=seed,
    )


def compare(expression_z, atrophy_d):
    """expression_z, atrophy_d: Series indexed by cortical region id, same regions."""
    ids = expression_z.index.intersection(atrophy_d.index).sort_values()
    x, y = expression_z.loc[ids].to_numpy(), atrophy_d.loc[ids].to_numpy()
    rho = spearmanr(x, y).statistic

    naive = naive_null(x, y)
    spins = spin_null(x)
    spun = np.array([spearmanr(spins[:, i], y, nan_policy="omit").statistic for i in range(spins.shape[1])])

    def two_sided(null):
        return (1 + (np.abs(null) >= abs(rho)).sum()) / (len(null) + 1)

    summary = pd.DataFrame(
        {"null_model": ["naive (shuffle regions)", "spatial (spin test)"],
         "rho": [rho, rho],
         "p": [two_sided(naive), two_sided(spun)],
         "n_perm": [len(naive), len(spun)]}
    )
    return summary, {"x": x, "y": y, "naive": naive, "spin": spun}
