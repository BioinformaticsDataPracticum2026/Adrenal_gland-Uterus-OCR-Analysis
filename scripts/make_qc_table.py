#!/usr/bin/env python3
"""
Extract QC metrics from ENCODE ATAC-seq pipeline qc.json files
and generate a summary table for the project report.

Usage:
    python make_qc_table.py \
        --human-adrenal /path/to/human_adrenal/qc.json \
        --mouse-adrenal /path/to/mouse_adrenal/qc.json \
        --human-uterus  /path/to/human_uterus/qc.json \
        --mouse-uterus  /path/to/mouse_uterus/qc.json \
        --out-tsv   results/qc_summary_table.tsv \
        --out-md    results/qc_summary_table.md
"""

import json
import argparse
import csv
import sys

DATASETS = [
    ("Human Adrenal", "human_adrenal", True),   # (label, arg_key, selected)
    ("Mouse Adrenal", "mouse_adrenal", True),
    ("Human Uterus",  "human_uterus",  False),
    ("Mouse Uterus",  "mouse_uterus",  False),
]

def safe_get(d, *keys, default="N/A", digits=3):
    """Navigate nested dict safely; round floats if found."""
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    if isinstance(d, float):
        return round(d, digits)
    return d if d is not None else default

def parse_qc(path):
    """Parse an ENCODE ATAC-seq qc.json and return a flat metrics dict per replicate."""
    with open(path) as f:
        q = json.load(f)

    results = {}
    for rep in ("rep1", "rep2"):
        r = {}

        # --- Mapped reads (from alignment flagstat) ---
        r["mapped_reads"] = safe_get(q, "align", rep, "samstat", "mapped",
                                     default=safe_get(q, "align", rep, "flagstat", "mapped"))

        # --- Distinct fragments (from library complexity) ---
        r["distinct_frags"] = safe_get(q, "lib_complexity", rep, "picard_est_lib_size")

        # --- NRF (Non-Redundant Fraction) ---
        r["NRF"] = safe_get(q, "lib_complexity", rep, "NRF")

        # --- TSS enrichment ---
        r["TSS_enrichment"] = safe_get(q, "tss_enrich", rep, "tss_enrich")

        # --- FRiP ---
        r["FRiP"] = safe_get(q, "frac_reads_in_peaks", rep, "frip")

        # --- Duplication rate ---
        r["pct_dup"] = safe_get(q, "dup", rep, "pct_duplicate_reads")

        # --- Mitochondrial reads ---
        r["pct_mito"] = safe_get(q, "align", rep, "samstat", "pct_mapped_reads_mito",
                                  default=safe_get(q, "align", rep, "pct_mito"))

        results[rep] = r

    # IDR peak counts (one per pair, not per replicate)
    pair_key = next((k for k in q.get("idr", {}) if "rep" in k), None)
    if pair_key:
        results["opt_peaks"]  = safe_get(q, "idr", pair_key, "opt_pks")
        results["cons_peaks"] = safe_get(q, "idr", pair_key, "consv_pks")
    else:
        results["opt_peaks"]  = "N/A"
        results["cons_peaks"] = "N/A"

    return results

def build_rows(paths):
    """Return list of row dicts for TSV/markdown output."""
    metrics = [
        ("TSS Enrichment (Rep1)",     lambda m: m["rep1"]["TSS_enrichment"]),
        ("TSS Enrichment (Rep2)",     lambda m: m["rep2"]["TSS_enrichment"]),
        ("FRiP (Rep1)",               lambda m: m["rep1"]["FRiP"]),
        ("FRiP (Rep2)",               lambda m: m["rep2"]["FRiP"]),
        ("Mapped Reads (Rep1)",       lambda m: m["rep1"]["mapped_reads"]),
        ("Mapped Reads (Rep2)",       lambda m: m["rep2"]["mapped_reads"]),
        ("Distinct Fragments (Rep1)", lambda m: m["rep1"]["distinct_frags"]),
        ("Distinct Fragments (Rep2)", lambda m: m["rep2"]["distinct_frags"]),
        ("NRF (Rep1)",                lambda m: m["rep1"]["NRF"]),
        ("NRF (Rep2)",                lambda m: m["rep2"]["NRF"]),
        ("% Duplication (Rep1)",      lambda m: m["rep1"]["pct_dup"]),
        ("% Duplication (Rep2)",      lambda m: m["rep2"]["pct_dup"]),
        ("% Mitochondrial (Rep1)",    lambda m: m["rep1"]["pct_mito"]),
        ("% Mitochondrial (Rep2)",    lambda m: m["rep2"]["pct_mito"]),
        ("IDR Optimal Peaks",         lambda m: m["opt_peaks"]),
        ("IDR Conservative Peaks",    lambda m: m["cons_peaks"]),
        ("Selected for Analysis",     None),   # filled manually below
    ]

    parsed = {}
    for label, arg_key, selected in DATASETS:
        path = paths.get(arg_key)
        if path:
            try:
                parsed[(label, selected)] = parse_qc(path)
            except Exception as e:
                print(f"Warning: could not parse {path}: {e}", file=sys.stderr)
                parsed[(label, selected)] = None
        else:
            parsed[(label, selected)] = None

    rows = []
    for metric_label, extractor in metrics:
        row = {"Metric": metric_label}
        for (ds_label, selected), m in parsed.items():
            if metric_label == "Selected for Analysis":
                row[ds_label] = "Yes ✓" if selected else "No"
            elif m is not None and extractor:
                row[ds_label] = extractor(m)
            else:
                row[ds_label] = "N/A"
        rows.append(row)

    return rows, [label for label, _, _ in DATASETS]

def write_tsv(rows, col_order, path):
    fieldnames = ["Metric"] + col_order
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"TSV written to {path}")

def write_markdown(rows, col_order, path):
    cols = ["Metric"] + col_order
    header = "| " + " | ".join(cols) + " |"
    sep    = "| " + " | ".join(["---"] * len(cols)) + " |"
    lines  = [header, sep]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(c, "N/A")) for c in cols) + " |")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Markdown table written to {path}")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--human-adrenal", dest="human_adrenal", metavar="QC_JSON")
    parser.add_argument("--mouse-adrenal", dest="mouse_adrenal", metavar="QC_JSON")
    parser.add_argument("--human-uterus",  dest="human_uterus",  metavar="QC_JSON")
    parser.add_argument("--mouse-uterus",  dest="mouse_uterus",  metavar="QC_JSON")
    parser.add_argument("--out-tsv", default="results/qc_summary_table.tsv")
    parser.add_argument("--out-md",  default="results/qc_summary_table.md")
    args = parser.parse_args()

    paths = {
        "human_adrenal": args.human_adrenal,
        "mouse_adrenal": args.mouse_adrenal,
        "human_uterus":  args.human_uterus,
        "mouse_uterus":  args.mouse_uterus,
    }

    rows, col_order = build_rows(paths)
    write_tsv(rows, col_order, args.out_tsv)
    write_markdown(rows, col_order, args.out_md)

if __name__ == "__main__":
    main()
