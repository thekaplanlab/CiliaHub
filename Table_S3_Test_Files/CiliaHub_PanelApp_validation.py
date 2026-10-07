#!/usr/bin/env python3
"""
CiliaHub validation against independent clinical gene panels (Supplementary Table S3).

Compares the CiliaHub v9.1 gene catalogue (2,786 genes) and its ciliopathy module
with PanelApp Australia gene panels, by PanelApp rating (green, amber, red).

Inputs
  --ciliahub  ciliahub_data.json (CiliaHub release v9.1)
  --hgnc      CiliaHub_PubMed_Search_Symbols.tsv (HGNC current/previous/alias symbols,
              columns: Gene, Ensembl_Gene_ID, Search_Symbol, Symbol_Type)
  --panel     one or more PanelApp Australia panel TSV downloads, as NAME=FILE
              (e.g. "Ciliopathies=Ciliopathies.tsv")
Outputs (in --out directory)
  panel_genes.tsv    one row per panel gene with the matched CiliaHub record
  panel_summary.tsv  recall, fold enrichment and one-sided Fisher exact P values

Usage
  python3 CiliaHub_PanelApp_validation.py \
      --ciliahub ciliahub_data.json --hgnc CiliaHub_PubMed_Search_Symbols.tsv \
      --panel "Ciliopathies=Ciliopathies.tsv" \
      --panel "Ciliary Dyskinesia=Ciliary_Dyskinesia.tsv" \
      --panel "Dilated Cardiomyopathy=Dilated_Cardiomyopathy.tsv" \
      --out results/

Requires Python 3.9+, pandas, scipy.
"""
import argparse, ast, json, os, re
import pandas as pd
from scipy.stats import fisher_exact

UNIVERSE = 19214  # protein-coding gene universe (Ensembl 116 x HGNC; Supplementary Table S1)
RATING = {"3": "Green", "2": "Amber", "1": "Red"}


def as_list(x):
    """CiliaHub list fields may be Python lists or their string representation."""
    if isinstance(x, list):
        return [str(i) for i in x if str(i).strip()]
    try:
        v = ast.literal_eval(x) if str(x).startswith("[") else []
    except (ValueError, SyntaxError):
        v = []
    return [str(i) for i in v if str(i).strip()]


def load_ciliahub(path):
    genes = json.load(open(path))["genes"]
    category = {g: ("Gold-standard" if v["evidence_tier"].startswith("Gold") else "Cilia-associated")
                for g, v in genes.items()}
    # Ciliopathy module: genes linked to >=1 ciliopathy (primary, motile or tissue-restricted);
    # genes linked only to secondary (non-ciliary) disorders are excluded.
    module = {g for g, v in genes.items()
              if any(c != "Secondary Diseases" for c in as_list(v.get("ciliopathy_classification", "")))}
    return genes, category, module


def build_matcher(genes, hgnc_path):
    ens2g, sym2g, syn2g = {}, {g.upper(): g for g in genes}, {}
    for g, v in genes.items():
        for e in re.findall(r"ENSG\d{11}", str(v.get("ensembl_id", ""))):
            ens2g.setdefault(e, g)
        for s in as_list(v.get("synonyms", "")) or re.split(r"[;,|]\s*", str(v.get("synonyms", ""))):
            s = s.strip().upper()
            if s and s not in sym2g:
                syn2g.setdefault(s, set()).add(g)
    H = pd.read_csv(hgnc_path, sep="\t", dtype=str)
    cur = H[H.Symbol_Type == "Current"]
    hgnc_ens = dict(zip(cur.Ensembl_Gene_ID, cur.Gene))
    hgnc_current = set(cur.Gene.str.upper())

    def unique_map(kind):  # previous/alias symbols that point to exactly one gene
        s = H[H.Symbol_Type == kind].groupby(H.Search_Symbol.str.upper()).Gene.agg(set)
        return {k: next(iter(v)) for k, v in s.items() if len(v) == 1}
    hgnc_prev, hgnc_alias = unique_map("Previous"), unique_map("Alias")

    def match(ensembl, symbol):
        """Return (CiliaHub gene or None, match method)."""
        ids = re.findall(r"ENSG\d{11}", str(ensembl or ""))
        for e in ids:                                          # 1. Ensembl ID
            if e in ens2g:
                return ens2g[e], "Ensembl ID"
        for e in ids:                                          # 2. Ensembl -> current HGNC symbol
            c = hgnc_ens.get(e)
            if c and c.upper() in sym2g:
                return sym2g[c.upper()], "Ensembl ID via HGNC current symbol"
        if not isinstance(symbol, str) or not symbol.strip():
            return None, "Not found"
        s = symbol.strip().upper()
        if s in sym2g:                                         # 3. exact current symbol
            return sym2g[s], "Symbol"
        if s in hgnc_current:                                  # a current symbol of a non-CiliaHub gene
            return None, "Not found"
        if s in hgnc_prev and hgnc_prev[s].upper() in sym2g:   # 4a. HGNC previous symbol
            return sym2g[hgnc_prev[s].upper()], "HGNC previous symbol"
        if s in syn2g and len(syn2g[s]) == 1:                  # 4b. CiliaHub synonym
            return next(iter(syn2g[s])), "CiliaHub synonym"
        if s in hgnc_alias and hgnc_alias[s].upper() in sym2g: # 4c. HGNC alias symbol
            return sym2g[hgnc_alias[s].upper()], "HGNC alias symbol"
        return None, "Not found"
    return match


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ciliahub", required=True)
    ap.add_argument("--hgnc", required=True)
    ap.add_argument("--panel", action="append", required=True, help="NAME=FILE")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    genes, category, module = load_ciliahub(a.ciliahub)
    match = build_matcher(genes, a.hgnc)
    n_hub, n_mod = len(genes), len(module)
    print(f"CiliaHub genes: {n_hub}; ciliopathy module genes: {n_mod}; universe: {UNIVERSE}")

    rows = []
    for spec in a.panel:
        name, path = spec.split("=", 1)
        p = pd.read_csv(path, sep="\t", dtype=str)
        p = p[p["Entity type"] == "gene"]            # copy-number regions are excluded
        for r in p.to_dict("records"):
            g, how = match(r.get("EnsemblId(GRch38)"), r.get("Gene Symbol"))
            v = genes.get(g, {})
            rows.append({
                "Panel": name, "Panel version": r.get("version"),
                "Gene symbol (panel)": r.get("Gene Symbol"),
                "PanelApp rating": RATING.get(str(r.get("GEL_Status")), "Unrated"),
                "Ensembl ID GRCh38 (panel)": r.get("EnsemblId(GRch38)"),
                "In CiliaHub": "Yes" if g else "No", "CiliaHub gene": g or "",
                "Match method": how, "CiliaHub category": category.get(g, "Not in CiliaHub"),
                "In ciliopathy module": "Yes" if g in module else "No",
                "CiliaHub ciliopathies": "; ".join(as_list(v.get("ciliopathies", ""))),
            })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.out, "panel_genes.tsv"), sep="\t", index=False)

    def enrichment(k, n, set_size):
        """Fold enrichment and one-sided Fisher exact P of k/n panel genes in a set of set_size genes."""
        table = [[k, n - k], [set_size - k, UNIVERSE - set_size - (n - k)]]
        return (k / n) / (set_size / UNIVERSE), fisher_exact(table, alternative="greater")[1]

    summ = []
    for name, sub in df.groupby("Panel", sort=False):
        for rating in ["Green", "Amber", "Red", "All"]:
            s = sub if rating == "All" else sub[sub["PanelApp rating"] == rating]
            n = len(s)
            if n == 0:
                continue
            k = (s["In CiliaHub"] == "Yes").sum()
            km = (s["In ciliopathy module"] == "Yes").sum()
            fh, ph = enrichment(k, n, n_hub)
            fm, pm = enrichment(km, n, n_mod)
            summ.append({
                "Panel": name, "PanelApp rating": rating, "Genes": n,
                "In CiliaHub": k, "% in CiliaHub": round(100 * k / n, 1),
                "Gold-standard": (s["CiliaHub category"] == "Gold-standard").sum(),
                "Cilia-associated": (s["CiliaHub category"] == "Cilia-associated").sum(),
                "Fold enrichment (CiliaHub)": round(fh, 1), "Fisher P (CiliaHub)": f"{ph:.3g}",
                "In ciliopathy module": km, "% in ciliopathy module": round(100 * km / n, 1),
                "Fold enrichment (module)": round(fm, 1), "Fisher P (module)": f"{pm:.3g}",
                "Not in CiliaHub": ", ".join(s.loc[s["In CiliaHub"] == "No", "Gene symbol (panel)"]),
            })
    out = pd.DataFrame(summ)
    out.to_csv(os.path.join(a.out, "panel_summary.tsv"), sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 20, "display.max_colwidth", 60):
        print(out.drop(columns=["Not in CiliaHub"]).to_string(index=False))


if __name__ == "__main__":
    main()
