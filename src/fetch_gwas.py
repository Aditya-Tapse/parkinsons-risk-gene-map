"""Step 2: genome-wide significant Parkinson's risk genes from the GWAS Catalog."""
import re

import pandas as pd

from src import config

# The catalog has shipped both classic and camelCase column names.
P_COLS = ["P-VALUE", "pValue", "p_value"]
GENE_COLS = ["MAPPED_GENE", "mappedGenes", "mapped_gene"]
LOCUS_COLS = ["REGION", "locations", "CHR_ID"]
TRAIT_COLS = ["MAPPED_TRAIT", "mappedTraits", "mapped_trait"]


def _pick(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"none of {candidates} in GWAS file columns: {list(df.columns)[:15]}...")


def _split_genes(cell):
    # "SNCA - MMRN1", "GBA1, GBAP1", "LRRK2; MUC19" -> individual symbols
    parts = re.split(r"\s+-\s+|\s*[,;]\s*|\s+x\s+", str(cell))
    keep = []
    for g in parts:
        g = g.strip()
        if not g or g.upper() in {"NA", "NR", "INTERGENIC"} or g.startswith("LOC"):
            continue
        keep.append(g)
    return keep


def load_genes(expression_columns=None):
    if not config.GWAS_RAW.exists():
        raise FileNotFoundError(
            f"Download the PD associations ({config.TRAIT_ID}) from "
            f"https://www.ebi.ac.uk/gwas/efotraits/{config.TRAIT_ID} "
            f"and save the TSV as {config.GWAS_RAW}"
        )
    gwas = pd.read_csv(config.GWAS_RAW, sep="\t", low_memory=False)
    p_col, gene_col = _pick(gwas, P_COLS), _pick(gwas, GENE_COLS)
    locus_col = next((c for c in LOCUS_COLS if c in gwas.columns), None)

    trait_col = _pick(gwas, TRAIT_COLS)
    risk = gwas[gwas[trait_col].str.strip() == config.GWAS_MAPPED_TRAIT].copy()
    print(f"{len(risk)} of {len(gwas)} rows are PD risk (mapped trait '{config.GWAS_MAPPED_TRAIT}')")

    risk[p_col] = pd.to_numeric(risk[p_col], errors="coerce")
    hits = risk[risk[p_col] < config.GWAS_P_THRESHOLD].dropna(subset=[gene_col])

    rows = []
    for _, r in hits.iterrows():
        for g in _split_genes(r[gene_col]):
            rows.append({"gene": g, "source_locus": r[locus_col] if locus_col else "", "p": r[p_col]})
    genes = (
        pd.DataFrame(rows)
        .sort_values("p")
        .drop_duplicates("gene")
        .reset_index(drop=True)
    )
    print(f"{len(hits)} associations at p < {config.GWAS_P_THRESHOLD:g} -> {len(genes)} unique mapped genes")

    genes = annotate(genes, expression_columns)
    genes.to_csv(config.GENES_CSV, index=False)
    print(genes["locus_group"].value_counts().to_string())

    coding = genes[genes["locus_group"] == "protein-coding gene"]
    matched = coding["expression_symbol"].notna()
    renamed = coding[matched & (coding["expression_symbol"] != coding["gene"])]
    print(f"\nprotein-coding: {len(coding)}, matched to expression data: {matched.sum()} ({matched.mean():.0%})")
    if len(renamed):
        print("matched via previous symbol:",
              ", ".join(f"{a}->{b}" for a, b in zip(renamed["gene"], renamed["expression_symbol"])))
    print("protein-coding but not on the array:", ", ".join(coding.loc[~matched, "gene"]))

    final = analysis_genes(genes)
    lo, hi = config.GENE_LIST_RANGE
    print(f"\nfinal gene set for analysis: {len(final)} genes")
    if not lo <= len(final) <= hi:
        print(f"WARNING: gene list outside the {lo}-{hi} range; adjust the threshold")
    return genes


def load_hgnc():
    if not config.HGNC_TSV.exists():
        import requests
        config.HGNC_TSV.parent.mkdir(parents=True, exist_ok=True)
        config.HGNC_TSV.write_bytes(requests.get(config.HGNC_URL, timeout=300).content)
    hgnc = pd.read_csv(config.HGNC_TSV, sep="\t", low_memory=False,
                       usecols=["symbol", "locus_group", "prev_symbol"])
    return hgnc


def annotate(genes, expression_columns=None):
    """Add HGNC locus_group and the symbol this gene has in the expression data.

    The 2013 microarray uses older gene names (e.g. GBA rather than GBA1), so a
    gene is matched by its current symbol first, then by previous symbols.
    Aliases are not used: they get reused across genes (BRIP1's alias BACH1 is
    a different gene), and neither is a previous symbol that is now another
    gene's approved symbol.
    """
    hgnc = load_hgnc()
    by_symbol = hgnc.set_index("symbol")
    approved = set(by_symbol.index)
    # symbols the GWAS catalog used that HGNC has since retired
    retired = {}
    for _, r in hgnc.dropna(subset=["prev_symbol"]).iterrows():
        for old in r["prev_symbol"].split("|"):
            retired.setdefault(old, r["symbol"])

    available = set(expression_columns) if expression_columns is not None else set()
    groups, expr_symbols = [], []
    for g in genes["gene"]:
        current = g if g in by_symbol.index else retired.get(g)
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
    """Protein-coding GWAS genes that are measured in the expression data, deduplicated."""
    keep = genes[(genes["locus_group"] == "protein-coding gene") & genes["expression_symbol"].notna()]
    return keep.drop_duplicates("expression_symbol")["expression_symbol"].tolist()


if __name__ == "__main__":
    from src.fetch_expression import fetch_expression

    load_genes(fetch_expression().columns)
