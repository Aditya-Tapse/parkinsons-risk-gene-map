"""Figures. z-scores use a diverging blue-red colormap with a neutral grey midpoint (z = 0)."""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.gridspec import GridSpec

from src import config

INK = "#0b0b0b"
INK_2 = "#52514e"
INK_3 = "#8a8983"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
ACCENT = "#e34948"
NULL_FILL = "#b7d3f6"

Z_CMAP = LinearSegmentedColormap.from_list(
    "z_div", ["#104281", "#3987e5", "#f0efec", "#e66767", "#9b1c1c"]
)

PRETTY = {
    "bankssts": "Banks of STS", "caudalanteriorcingulate": "Caudal anterior cingulate",
    "caudalmiddlefrontal": "Caudal middle frontal", "cuneus": "Cuneus",
    "entorhinal": "Entorhinal", "fusiform": "Fusiform", "inferiorparietal": "Inferior parietal",
    "inferiortemporal": "Inferior temporal", "isthmuscingulate": "Isthmus cingulate",
    "lateraloccipital": "Lateral occipital", "lateralorbitofrontal": "Lateral orbitofrontal",
    "lingual": "Lingual", "medialorbitofrontal": "Medial orbitofrontal",
    "middletemporal": "Middle temporal", "parahippocampal": "Parahippocampal",
    "paracentral": "Paracentral", "parsopercularis": "Pars opercularis",
    "parsorbitalis": "Pars orbitalis", "parstriangularis": "Pars triangularis",
    "pericalcarine": "Pericalcarine", "postcentral": "Postcentral",
    "posteriorcingulate": "Posterior cingulate", "precentral": "Precentral",
    "precuneus": "Precuneus", "rostralanteriorcingulate": "Rostral anterior cingulate",
    "rostralmiddlefrontal": "Rostral middle frontal", "superiorfrontal": "Superior frontal",
    "superiorparietal": "Superior parietal", "superiortemporal": "Superior temporal",
    "supramarginal": "Supramarginal", "frontalpole": "Frontal pole",
    "temporalpole": "Temporal pole", "transversetemporal": "Transverse temporal",
    "insula": "Insula", "caudate": "Caudate", "putamen": "Putamen", "pallidum": "Pallidum",
    "accumbensarea": "Nucleus accumbens", "thalamusproper": "Thalamus", "thalamus": "Thalamus",
    "hippocampus": "Hippocampus", "amygdala": "Amygdala", "brainstem": "Brainstem",
}


def pretty(label):
    return PRETTY.get(label, label)


def _style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "Arial", "DejaVu Sans"],
        "font.size": 11, "text.color": INK, "axes.labelcolor": INK_2,
        "axes.edgecolor": GRID, "xtick.color": INK_2, "ytick.color": INK_2,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    })


def _save(fig, name, out_dir=config.FIGURES):
    fig.savefig(out_dir / f"{name}.png", dpi=200)
    plt.close(fig)
    print(f"saved {out_dir / name}.png")


def _vertex_values(label_gifti, values_by_id):
    """Paint parcel values onto fsaverage5 vertices; unlabelled vertices -> 0 (hidden)."""
    import nibabel as nib
    labels = np.asarray(nib.load(label_gifti).agg_data()).astype(int)
    out = np.zeros(labels.shape, dtype=float)
    for rid, val in values_by_id.items():
        out[labels == rid] = val
    return out


def hero(results, title, subtitle, footer, out_name="hero", out_dir=config.FIGURES, vmax=3.0):
    """Cortical z-map in four views + ranked subcortex/brainstem bars.

    results: DataFrame indexed by region id with columns label, hemisphere,
    structure, z, q.
    """
    from nilearn import datasets, plotting
    from src.fetch_expression import load_atlas

    _style()
    fsavg = datasets.fetch_surf_fsaverage("fsaverage5")
    lh_gii, rh_gii = load_atlas(surface=True)["image"]
    z_by_id = results["z"].to_dict()

    fig = plt.figure(figsize=(16, 9))
    gs = GridSpec(
        3, 3, figure=fig, width_ratios=[1, 1, 0.95], height_ratios=[1, 1, 0.22],
        left=0.0, right=0.97, top=0.84, bottom=0.06, wspace=0.0, hspace=0.0,
    )
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    views = [("left", "lateral", lh_gii, 0, 0), ("right", "lateral", rh_gii, 0, 1),
             ("left", "medial", lh_gii, 1, 0), ("right", "medial", rh_gii, 1, 1)]
    for hemi, view, gii, r, c in views:
        ax = fig.add_subplot(gs[r, c], projection="3d")
        plotting.plot_surf_stat_map(
            fsavg[f"infl_{hemi}"], _vertex_values(gii, z_by_id), hemi=hemi, view=view,
            bg_map=fsavg[f"sulc_{hemi}"], bg_on_data=True,
            cmap=Z_CMAP, vmax=vmax, symmetric_cbar=True, threshold=1e-6,
            colorbar=False, axes=ax, figure=fig,
        )
        ax.set_facecolor(SURFACE)
        ax.set_box_aspect(None, zoom=1.45)
        ax.text2D(0.5, 0.0, f"{hemi.title()} · {view}", transform=ax.transAxes,
                  ha="center", color=INK_3, fontsize=10)

    # Subcortex + brainstem, averaged across hemispheres
    deep = results[results["structure"] != "cortex"]
    deep = deep.groupby("label").agg(z=("z", "mean"), q=("q", "min")).sort_values("z")
    ax = fig.add_subplot(gs[0:2, 2])
    y = np.arange(len(deep))
    ax.barh(y, deep["z"], height=0.55, color=Z_CMAP(norm(deep["z"].clip(-vmax, vmax))))
    ax.axvline(0, color=INK_3, lw=1)
    for yi, (lab, row) in zip(y, deep.iterrows()):
        star = " *" if row["q"] < config.FDR_ALPHA else ""
        off = 0.08 if row["z"] >= 0 else -0.08
        label = f"{row['z']:.1f}".replace("-0.0", "0.0")
        ax.text(row["z"] + off, yi, label + star, va="center",
                ha="left" if row["z"] >= 0 else "right", color=INK, fontsize=10)
    ax.set_yticks(y, [pretty(l) for l in deep.index], fontsize=11, color=INK)
    lim = max(vmax, np.abs(deep["z"]).max() * 1.25)
    ax.set_xlim(-lim, lim)
    ax.set_xlabel("z vs. 10,000 random gene sets")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_title("Subcortex & brainstem", loc="left", fontsize=12, color=INK_2, pad=10)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)

    cax = fig.add_axes([0.16, 0.095, 0.30, 0.022])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=Z_CMAP), cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.set_ticks([-vmax, 0, vmax], labels=[f"−{vmax:g}  lower than random", "0", f"+{vmax:g}  higher than random"])
    cb.ax.tick_params(length=0, labelsize=10, colors=INK_2)

    fig.text(0.02, 0.945, title, fontsize=22, fontweight="bold", color=INK)
    fig.text(0.02, 0.895, subtitle, fontsize=13, color=INK_2)
    fig.text(0.97, 0.03, footer, fontsize=9.5, color=INK_3, ha="right")
    if (results["q"] < config.FDR_ALPHA).any():
        fig.text(0.02, 0.03, "* FDR q < 0.05", fontsize=9.5, color=INK_3)
    _save(fig, out_name, out_dir)


def null_histogram(null_col, observed, region_name, out_dir=config.FIGURES):
    _style()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.hist(null_col, bins=60, color=NULL_FILL, edgecolor=SURFACE, linewidth=0.5)
    ax.axvline(observed, color=ACCENT, lw=2)
    ax.text(observed, ax.get_ylim()[1] * 0.92, "  PD risk genes", color=INK, va="top")
    ax.set_xlabel("Mean expression of gene set")
    ax.set_ylabel("Random gene sets")
    ax.set_title(f"{region_name}: real genes vs. 10,000 random sets", loc="left", fontsize=14)
    fig.tight_layout()
    _save(fig, "null_histogram", out_dir)


def braak_plot(trend, out_dir=config.FIGURES):
    _style()
    t = trend["table"]
    rng = np.random.default_rng(config.SEED)
    fig, ax = plt.subplots(figsize=(8, 5))
    groups = [("1–2\nbrainstem", [1, 2]), ("3\namygdala", [3]), ("4\nmesocortex", [4]),
              ("5\nassociation\ncortex", [5]), ("6\nprimary\ncortex", [6])]
    for i, (_, st) in enumerate(groups):
        vals = t.loc[t["stage"].isin(st), "z"].to_numpy()
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(vals)), vals, s=36, color="#2a78d6",
                   edgecolor=SURFACE, linewidth=1.5, zorder=3)
        if len(vals):
            ax.hlines(np.median(vals), i - 0.25, i + 0.25, color=INK, lw=2, zorder=4)
    ax.set_xticks(range(len(groups)), [g[0] for g in groups])
    ax.set_xlabel("Braak stage (order in which Parkinson's pathology reaches each region)")
    ax.axhline(0, color=INK_3, lw=1)
    ax.set_ylabel("z vs. random gene sets")
    fig.suptitle("Risk-gene activity does not follow the spread of Parkinson's",
                 x=0.02, ha="left", fontsize=14, color=INK)
    ax.set_title(f"Prediction: highest in early stages, falling with each stage.  "
                 f"Observed ρ = {trend['rho']:.2f}, p = {trend['p']:.2f}  ·  bar = median",
                 loc="left", fontsize=10, color=INK_2)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    _save(fig, "braak_trend", out_dir)


def spatial_null_plot(spatial, summary, out_dir=config.FIGURES):
    """Scatter of expression vs atrophy, plus naive vs spin null distributions."""
    _style()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [1, 1.2]})
    a1.scatter(spatial["x"], spatial["y"], s=40, color="#2a78d6", edgecolor=SURFACE, linewidth=1.5)
    a1.axhline(0, color=INK_3, lw=1)
    a1.set_xlabel("PD risk-gene expression (z vs. random gene sets)")
    a1.set_ylabel("Cortical thickness in PD vs. controls (Cohen's d)\n← thinner in PD")
    rho = summary["rho"].iloc[0]
    a1.set_title(f"68 cortical regions  ·  Spearman ρ = {rho:.2f}", loc="left", fontsize=12, color=INK_2)
    a1.grid(color=GRID, lw=0.8)
    a1.set_axisbelow(True)

    bins = np.linspace(-0.6, 0.6, 49)
    a2.hist(spatial["naive"], bins=bins, density=True, color=INK_3, alpha=0.45, label="Naive null (shuffled regions)")
    a2.hist(spatial["spin"], bins=bins, density=True, color="#2a78d6", alpha=0.55, label="Spatial null (spin test)")
    a2.axvline(rho, color=ACCENT, lw=2)
    a2.text(rho, a2.get_ylim()[1] * 0.97, "  observed", color=INK, va="top")
    p_naive, p_spin = summary["p"].iloc[0], summary["p"].iloc[1]
    a2.set_title(f"Naive p = {p_naive:.2f}  ·  Spin p = {p_spin:.2f}  (spin null is wider: neighbours look alike)",
                 loc="left", fontsize=12, color=INK_2)
    a2.set_xlabel("Correlation expected by chance")
    a2.set_yticks([])
    a2.spines["left"].set_visible(False)
    a2.legend(frameon=False, loc="upper right", labelcolor=INK_2)
    fig.suptitle("Risk-gene activity does not match where patients' cortex thins",
                 x=0.02, ha="left", fontsize=14, color=INK)
    fig.tight_layout()
    _save(fig, "spatial_null", out_dir)
