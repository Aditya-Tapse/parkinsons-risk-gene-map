"""Step 4: does PD-gene expression follow the Braak staging sequence?

Braak et al. 2003: Lewy pathology starts in the lower brainstem (stages 1-2),
reaches the nigra/amygdala (3), the anteromedial temporal mesocortex (4), then
association neocortex (5) and finally primary sensory/motor cortex (6).
Desikan-Killiany has a single brainstem region and no substantia nigra, so the
mapping below is coarse. Basal ganglia and thalamus are not part of the
cortical Braak sequence and are left out of the trend test.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

BRAAK_STAGE = {
    "brainstem": 1,
    "amygdala": 3,
    "entorhinal": 4, "parahippocampal": 4, "hippocampus": 4, "temporalpole": 4,
    "precentral": 6, "postcentral": 6, "paracentral": 6, "pericalcarine": 6,
    "transversetemporal": 6, "cuneus": 6, "lingual": 6,
}
SUBCORTICAL_EXCLUDED = {"caudate", "putamen", "pallidum", "accumbensarea", "thalamusproper", "thalamus"}


def stage_of(label):
    if label in BRAAK_STAGE:
        return BRAAK_STAGE[label]
    if label in SUBCORTICAL_EXCLUDED:
        return np.nan
    return 5  # remaining neocortex = high-order association cortex


def braak_trend(results, null):
    """Spearman(stage, score) vs. the same statistic on every random gene set.

    A negative rho means expression is highest in early-affected regions.
    """
    stages = results["label"].map(stage_of)
    keep = stages.notna().to_numpy()
    s = stages[keep].to_numpy()

    rho_obs = spearmanr(s, results["observed"].to_numpy()[keep]).statistic
    rho_null = np.array([spearmanr(s, row[keep]).statistic for row in null])
    p = (1 + (rho_null <= rho_obs).sum()) / (len(rho_null) + 1)
    table = pd.DataFrame({"stage": stages, "z": results["z"], "name": results["name"]})
    return {"rho": rho_obs, "p": p, "rho_null": rho_null, "table": table.dropna()}
