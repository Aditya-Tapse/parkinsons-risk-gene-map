"""Step 1: Allen Human Brain Atlas expression, parcellated into Desikan-Killiany (83 regions)."""
import pandas as pd

from src import config


def load_atlas(surface=False):
    import abagen
    return abagen.fetch_desikan_killiany(surface=surface, data_dir=str(config.ABAGEN_DIR))


def region_info():
    """Region id -> label / hemisphere / structure, indexed by id."""
    info = pd.read_csv(load_atlas()["info"]).set_index("id")
    info["name"] = info["label"] + "_" + info["hemisphere"]
    return info


def fetch_expression(force=False):
    if config.EXPRESSION_CSV.exists() and not force:
        return pd.read_csv(config.EXPRESSION_CSV, index_col=0)

    import abagen
    atlas = load_atlas()
    # ~4 GB download on first run; data_dir makes it resumable and cached.
    expression = abagen.get_expression_data(
        atlas["image"], atlas["info"],
        lr_mirror=config.LR_MIRROR,
        data_dir=str(config.ABAGEN_DIR),
        verbose=1,
    )
    expression.to_csv(config.EXPRESSION_CSV)

    n_regions, n_genes = expression.shape
    print(f"expression matrix: {n_regions} regions x {n_genes} genes")
    assert 75 <= n_regions <= 83 and n_genes > 10_000, "unexpected expression shape"
    empty = expression.index[expression.isna().all(axis=1)].tolist()
    if empty:
        print(f"regions with no samples: {empty}")
    return expression


if __name__ == "__main__":
    fetch_expression()
