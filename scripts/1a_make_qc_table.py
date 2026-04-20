#!/usr/bin/env python3
"""
Extract QC metrics from ENCODE ATAC-seq pipeline qc.json files
and generate a summary table for the project report.

Usage:
    python scripts/make_qc_table.py \
        --human-adrenal /ocean/projects/bio200034p/ikaplow/HumanDNase/AdrenalGlandFemaleGoodReps/atac/qc/qc.json \
        --mouse-adrenal /ocean/projects/bio200034p/ikaplow/MouseDNase/AdrenalGland2025/atac/qc/qc.json \
        --human-uterus  /ocean/projects/bio200034p/ikaplow/HumanDNase/UterusFemale/atac/qc/qc.json \
        --mouse-uterus  /ocean/projects/bio200034p/ikaplow/MouseDNase/Uterus2025/atac/qc/qc.json \
        --out-tsv results/qc_summary_table.tsv \
        --out-md  results/qc_summary_table.md
"""

import json
import argparse
import csv
import sys
from pathlib import Path

DATASETS = [
    ("Human Adrenal", "human_adrenal", True),
    ("Mouse Adrenal", "mouse_adrenal", True),
    ("Human Uterus",  "human_uterus",  False),
    ("Mouse Uterus",  "mouse_uterus",  False),
]

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

def safe_get(d, *keys, default="N/A", digits=4):
    """Navigate nested dict safely; round floats if found."""
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    if isinstance(d, float):
        return round(d, digits)
    return d if d is not None else default

def fmt_int(val):
    """Format large integers with commas."""
    if isinstance(val, int):
        return f"{val:,}"
    return val

def parse_qc(path):
    """Parse an ENCODE ATAC-seq qc.json and return a flat metrics dict."""
    with open(path) as f:
        q = json.load(f)

    result = {}

    for rep in ("rep1", "rep2"):
        r = {}

        # Mapped reads (post-alignment, pre-dedup)
        r["mapped_reads"] = fmt_int(safe_get(q, "align", "samstat", rep, "mapped_reads"))

        # Distinct fragments (library complexity)
        r["distinct_frags"] = fmt_int(safe_get(q, "lib_complexity", "lib_complexity", rep, "distinct_fragments"))

        # NRF (Non-Redundant Fraction) - closer to 1 is better; >0.9 is ideal
        r["NRF"] = safe_get(q, "lib_complexity", "lib_complexity", rep, "NRF")

        # TSS enrichment - measure of signal-to-noise; >6 is acceptable, >10 is good
        r["TSS_enrichment"] = safe_get(q, "align_enrich", "tss_enrich", rep, "tss_enrich", digits=2)

        # Duplication rate
        r["pct_dup"] = safe_get(q, "align", "dup", rep, "pct_duplicate_reads")

        # Mitochondrial fraction
        r["frac_mito"] = safe_get(q, "align", "frac_mito", rep, "frac_mito_reads")

        result[rep] = r

    # IDR peak counts (dataset-level, not per replicate)
    result["opt_peaks"]  = fmt_int(safe_get(q, "replication", "reproducibility", "idr", "N_opt"))
    result["cons_peaks"] = fmt_int(safe_get(q, "replication", "reproducibility", "idr", "N_consv"))
    result["reproducibility"] = safe_get(q, "replication", "reproducibility", "idr", "reproducibility")

    # FRiP against IDR peak sets
    # pooled-pr1_vs_pooled-pr2 = optimal set; rep1_vs_rep2 = conservative set
    result["frip_optimal"]       = safe_get(q, "peak_enrich", "frac_reads_in_peaks", "idr", "pooled-pr1_vs_pooled-pr2", "frip")
    result["frip_conservative"]  = safe_get(q, "peak_enrich", "frac_reads_in_peaks", "idr", "rep1_vs_rep2", "frip")

    return result

def build_rows(paths):
    """Return list of row dicts for TSV/markdown output."""
    metrics = [
        ("Mapped Reads (Rep1)",           lambda m: m["rep1"]["mapped_reads"]),
        ("Mapped Reads (Rep2)",           lambda m: m["rep2"]["mapped_reads"]),
        ("Distinct Fragments (Rep1)",     lambda m: m["rep1"]["distinct_frags"]),
        ("Distinct Fragments (Rep2)",     lambda m: m["rep2"]["distinct_frags"]),
        ("NRF (Rep1)",                    lambda m: m["rep1"]["NRF"]),
        ("NRF (Rep2)",                    lambda m: m["rep2"]["NRF"]),
        ("TSS Enrichment (Rep1)",         lambda m: m["rep1"]["TSS_enrichment"]),
        ("TSS Enrichment (Rep2)",         lambda m: m["rep2"]["TSS_enrichment"]),
        ("% Duplication (Rep1)",          lambda m: m["rep1"]["pct_dup"]),
        ("% Duplication (Rep2)",          lambda m: m["rep2"]["pct_dup"]),
        ("Mito Fraction (Rep1)",          lambda m: m["rep1"]["frac_mito"]),
        ("Mito Fraction (Rep2)",          lambda m: m["rep2"]["frac_mito"]),
        ("FRiP - Optimal Set",            lambda m: m["frip_optimal"]),
        ("FRiP - Conservative Set",       lambda m: m["frip_conservative"]),
        ("IDR Optimal Peaks",             lambda m: m["opt_peaks"]),
        ("IDR Conservative Peaks",        lambda m: m["cons_peaks"]),
        ("IDR Reproducibility",           lambda m: m["reproducibility"]),
        ("Selected for Analysis",         None),
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
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["Metric"] + col_order
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"TSV written to {path}")

def write_markdown(rows, col_order, path):
    path.parent.mkdir(parents=True, exist_ok=True)
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
    parser.add_argument("--out-tsv", type=Path, default=REPO_ROOT / "results" / "qc_summary_table.tsv")
    parser.add_argument("--out-md",  type=Path, default=REPO_ROOT / "results" / "qc_summary_table.md")
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
