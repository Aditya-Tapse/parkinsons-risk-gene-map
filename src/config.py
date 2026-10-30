"""Paths, random seed and analysis parameters."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"
DATA_EXTERNAL = ROOT / "data" / "external"
ABAGEN_DIR = DATA_RAW / "abagen"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "figures"

EXPRESSION_CSV = DATA_PROCESSED / "ahba_dk.csv"
GWAS_RAW = DATA_RAW / "gwas_pd_associations.tsv"
GENES_CSV = TABLES / "pd_genes.csv"
HGNC_TSV = DATA_RAW / "hgnc" / "hgnc_complete_set.txt"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
REGION_SCORES = TABLES / "region_scores.csv"
ATROPHY_CSV = DATA_EXTERNAL / "pd_atrophy_dk.csv"

# Only 2 of the 6 AHBA donors have right-hemisphere samples, so samples are
# mirrored across hemispheres.
LR_MIRROR = "bidirectional"

# GWAS Catalog trait (Parkinson disease). Only associations whose mapped trait is
# PD risk itself are kept; the trait page also lists progression, age at onset,
# levodopa response and a grey-matter pleiotropy study.
TRAIT_ID = "MONDO_0005180"
GWAS_MAPPED_TRAIT = "Parkinson disease"
GWAS_P_THRESHOLD = 5e-8

SEED = 0
N_PERM_GENES = 10_000
N_PERM_SPATIAL = 1_000
FDR_ALPHA = 0.05

# Loci with long-range LD and many mapped genes (MHC, 17q21.31/MAPT, LRRK2).
COMPLEX_LOCI_GENES = {
    "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1", "BTNL2",
    "MAPT", "KANSL1", "CRHR1", "LRRC37A", "LRRC37A2", "ARL17A", "ARL17B",
    "NSF", "WNT3", "LRRK2",
}

# Striatal medium spiny neuron markers; should score highest in caudate/putamen.
POSITIVE_CONTROL = ["DRD1", "DRD2", "PPP1R1B", "PENK", "TAC1", "ADORA2A", "GPR88", "RGS9"]

for _d in (DATA_RAW, DATA_PROCESSED, DATA_EXTERNAL, ABAGEN_DIR, TABLES, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)
