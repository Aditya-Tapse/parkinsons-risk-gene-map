"""Genome-wide significant Parkinson's risk genes from the GWAS Catalog."""
import re

import pandas as pd
import requests

from src import config


def _split_genes(cell):
    """'SNCA - MMRN1', 'GBA1, GBAP1', 'LRRK2; MUC19' -> individual symbols."""
    parts = [g.strip() for g in re.split(r"\s+-\s+|\s*[,;]\s*|\s+x\s+", str(cell))]
    return [g for g in parts if g and g.upper() not in {"NA", "NR", "INTERGENIC"} and not g.startswith("LOC")]


def load_genes(expression_columns):
    """All mapped genes at p < 5e-8, annotated with HGNC locus group and expression-data symbol."""
    gwas = pd.read_csv(config.GWAS_RAW, sep="\t", low_memory=False)
    risk = gwas[gwas["MAPPED_TRAIT"].str.strip() == config.GWAS_MAPPED_TRAIT].copy()
    risk["P-VALUE"] = pd.to_numeric(risk["P-VALUE"], errors="coerce")
    hits = risk[risk["P-VALUE"] < config.GWAS_P_THRESHOLD].dropna(subset=["MAPPED_GENE"])

    rows = [{"gene": g, "source_locus": r["REGION"], "p": r["P-VALUE"]}
            for _, r in hits.iterrows() for g in _split_genes(r["MAPPED_GENE"])]
    genes = pd.DataFrame(rows).sort_values("p").drop_duplicates("gene").reset_index(drop=True)

    genes = annotate(genes, expression_columns)
    genes.to_csv(config.GENES_CSV, index=False)
    print(f"{len(hits)} associations -> {len(genes)} mapped genes -> {len(analysis_genes(genes))} analysed")
    return genes


def load_hgnc():
    if not config.HGNC_TSV.exists():
        config.HGNC_TSV.parent.mkdir(parents=True, exist_ok=True)
        config.HGNC_TSV.write_bytes(requests.get(config.HGNC_URL, timeout=300).content)
    return pd.read_csv(config.HGNC_TSV, sep="\t", low_memory=False,
                       usecols=["symbol", "locus_group", "prev_symbol"])


def annotate(genes, expression_columns):
    """Add HGNC locus_group and the symbol each gene has in the expression data.

    The AHBA microarray uses older gene names (e.g. GBA rather than GBA1), so a
    gene is matched by its current symbol first, then by its previous symbols.
    Aliases are not used because they are shared between genes (BRIP1's alias
    BACH1 is a different gene), and neither is a previous symbol that is now
    another gene's approved symbol.
    """
    hgnc = load_hgnc()
    by_symbol = hgnc.set_index("symbol")
    approved = set(by_symbol.index)
    retired = {}
    for _, r in hgnc.dropna(subset=["prev_symbol"]).iterrows():
        for old in r["prev_symbol"].split("|"):
            retired.setdefault(old, r["symbol"])

    available = set(expression_columns)
    groups, expr_symbols = [], []
    for g in genes["gene"]:
        current = g if g in approved else retired.get(g)
        groups.append(by_symbol.at[current, "locus_group"] if current else "not in HGNC")

        candidates = [g]
        if current:
            candidates.append(current)
            prev = by_symbol.at[current, "prev_symbol"]
            if isinstance(prev, str):
                candidates += [p for p in prev.split("|") if p not in approved]
        expr_symbols.append(next((c for c in candidates if c in available), None))

    return genes.assign(locus_group=groups, expression_symbol=expr_symbols)


def analysis_genes(genes):
    """Protein-coding genes that are measured in the expression data."""
    keep = genes[(genes["locus_group"] == "protein-coding gene") & genes["expression_symbol"].notna()]
    return keep.drop_duplicates("expression_symbol")["expression_symbol"].tolist()
