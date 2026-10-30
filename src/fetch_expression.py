"""Allen Human Brain Atlas expression, parcellated into the Desikan-Killiany atlas (83 regions)."""
import abagen
import pandas as pd

from src import config


def load_atlas(surface=False):
    return abagen.fetch_desikan_killiany(surface=surface, data_dir=str(config.ABAGEN_DIR))


def region_info():
    """Region id -> label, hemisphere, structure."""
    info = pd.read_csv(load_atlas()["info"]).set_index("id")
    info["name"] = info["label"] + "_" + info["hemisphere"]
    return info


def fetch_expression():
    """Regions x genes expression matrix. Downloads AHBA (~4 GB) on first call, then reads the cache."""
    if config.EXPRESSION_CSV.exists():
        return pd.read_csv(config.EXPRESSION_CSV, index_col=0)

    atlas = load_atlas()
    expression = abagen.get_expression_data(
        atlas["image"], atlas["info"],
        lr_mirror=config.LR_MIRROR,
        data_dir=str(config.ABAGEN_DIR),
        verbose=1,
    )
    expression.to_csv(config.EXPRESSION_CSV)
    print(f"expression matrix: {expression.shape[0]} regions x {expression.shape[1]} genes")
    return expression
