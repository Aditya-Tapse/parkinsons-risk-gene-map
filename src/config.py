"""Paths, seeds and analysis choices. Every decision the write-up mentions lives here."""
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
NULL_NPY = DATA_PROCESSED / "gene_set_null.npy"
ATROPHY_CSV = DATA_EXTERNAL / "pd_atrophy_dk.csv"

# Only 2 of 6 AHBA donors have right-hemisphere samples; mirror samples across
# hemispheres so every region gets data. Stated in the README limitations.
LR_MIRROR = "bidirectional"

DISEASE = "Parkinson's disease"
# GWAS Catalog trait ID. The catalog moved disease traits to MONDO; this was EFO_0002508.
TRAIT_ID = "MONDO_0005180"
# Keep only rows mapped to PD risk itself. The trait page also lists progression,
# age at onset, levodopa response and a grey-matter-atrophy pleiotropy study;
# the last would make the atrophy comparison circular.
GWAS_MAPPED_TRAIT = "Parkinson disease"
GWAS_P_THRESHOLD = 5e-8
GENE_LIST_RANGE = (20, 500)

SEED = 0
N_PERM_GENES = 10_000
N_PERM_SPATIAL = 1_000
FDR_ALPHA = 0.05

# Loci with long-range LD / many mapped genes; dropped in a sensitivity check.
COMPLEX_LOCI_GENES = {
    "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1", "BTNL2",
    "MAPT", "KANSL1", "CRHR1", "LRRC37A", "LRRC37A2", "ARL17A", "ARL17B",
    "NSF", "WNT3", "LRRK2",
}

# Positive control: striatal medium spiny neuron markers should top caudate/putamen.
POSITIVE_CONTROL = ["DRD1", "DRD2", "PPP1R1B", "PENK", "TAC1", "ADORA2A", "GPR88", "RGS9"]

for _d in (DATA_RAW, DATA_PROCESSED, DATA_EXTERNAL, ABAGEN_DIR, TABLES, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)
