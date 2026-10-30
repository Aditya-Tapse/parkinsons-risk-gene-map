"""Trend test: does PD-gene expression decrease along the Braak staging sequence?

Braak et al. 2003: Lewy pathology starts in the lower brainstem (stages 1-2),
reaches the substantia nigra and amygdala (3), the anteromedial temporal
mesocortex (4), high-order association cortex (5) and finally primary
sensory/motor cortex (6). Desikan-Killiany has a single brainstem region and
no substantia nigra. Basal ganglia and thalamus are not part of the cortical
sequence and are excluded.
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata

BRAAK_STAGE = {
    "brainstem": 1,
    "amygdala": 3,
    "entorhinal": 4, "parahippocampal": 4, "hippocampus": 4, "temporalpole": 4,
    "precentral": 6, "postcentral": 6, "paracentral": 6, "pericalcarine": 6,
    "transversetemporal": 6, "cuneus": 6, "lingual": 6,
}
EXCLUDED = {"caudate", "putamen", "pallidum", "accumbensarea", "thalamusproper", "thalamus"}


def stage_of(label):
    if label in BRAAK_STAGE:
        return BRAAK_STAGE[label]
    if label in EXCLUDED:
        return np.nan
    return 5  # remaining neocortex = association cortex


def _spearman_columns(x, Y):
    """Spearman correlation of vector x with every column of Y."""
    rx = rankdata(x)
    rY = rankdata(Y, axis=0)
    rx = (rx - rx.mean()) / rx.std()
    rY = (rY - rY.mean(axis=0)) / rY.std(axis=0)
    return (rx[:, None] * rY).mean(axis=0)


def braak_trend(results, null, null_ids):
    """Spearman(stage, expression), tested against the same statistic for every random gene set.

    results:  region table with label, observed, z
    null:     (n_perm, n_regions) random-set means, columns in `null_ids` order
    Left and right hemispheres are averaged first, since after mirroring they
    are near copies. A negative rho means expression is highest in early stages.
    """
    df = results.loc[list(null_ids), ["label", "observed", "z"]].rename_axis("id")
    null_df = pd.DataFrame(null.T, index=df["label"].to_numpy())

    by_label = df.groupby("label")[["observed", "z"]].mean()
    by_label["stage"] = by_label.index.map(stage_of)
    by_label = by_label.dropna(subset=["stage"]).sort_values("stage")
    null_by_label = null_df.groupby(level=0).mean().loc[by_label.index].to_numpy()

    stage = by_label["stage"].to_numpy()
    rho_obs = _spearman_columns(stage, by_label[["observed"]].to_numpy())[0]
    rho_null = _spearman_columns(stage, null_by_label)
    p = (1 + (rho_null <= rho_obs).sum()) / (len(rho_null) + 1)
    return {"rho": rho_obs, "p": p, "table": by_label}
