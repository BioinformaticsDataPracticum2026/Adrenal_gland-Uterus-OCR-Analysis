#!/usr/bin/env python3
"""
Extract QC metrics from ENCODE ATAC-seq QC reports
and generate a summary table for the project report.

By default this script reads the HTML reports under data/qc_html
using paths resolved relative to the repository root, so it can be
run directly from the scripts directory.
"""

import json
import argparse
import csv
import sys
from pathlib import Path
from html.parser import HTMLParser

DATASETS = [
    ("Human Adrenal", "human_adrenal", True),
    ("Mouse Adrenal", "mouse_adrenal", True),
    ("Human Uterus",  "human_uterus",  False),
    ("Mouse Uterus",  "mouse_uterus",  False),
]

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
QC_HTML_DIR = REPO_ROOT / "data" / "qc_html"
QC_RESULTS_DIR = REPO_ROOT / "results" / "qc"

DEFAULT_INPUTS = {
    "human_adrenal": QC_HTML_DIR / "Human_AdrenalGlandFemale_qc.html",
    "mouse_adrenal": QC_HTML_DIR / "Mouse_AdrenalGlandFemale_qc.html",
    "human_uterus": QC_HTML_DIR / "Human_UterusFemale_qc.html",
    "mouse_uterus": QC_HTML_DIR / "Mouse_UterusFemale_qc.html",
}

SUMMARY_LABELS = {
    "human_adrenal": "Human_Adrenal",
    "human_uterus": "Human_Uterus",
    "mouse_adrenal": "Mouse_Adrenal",
    "mouse_uterus": "Mouse_Uterus",
}


class QCHtmlParser(HTMLParser):
    """Extract section headings and table cells from the ENCODE QC HTML report."""

    def __init__(self):
        super().__init__()
        self.current_heading = None
        self.in_heading = False
        self.heading_text = []
        self.in_table = False
        self.current_table = None
        self.current_row = None
        self.current_cell = None
        self.tables = {}

    def handle_starttag(self, tag, attrs):
        if tag == "h2":
            self.in_heading = True
            self.heading_text = []
        elif tag == "table":
            self.in_table = True
            self.current_table = []
        elif tag == "tr" and self.in_table:
            self.current_row = []
        elif tag in {"th", "td"} and self.in_table and self.current_row is not None:
            self.current_cell = []

    def handle_data(self, data):
        if self.in_heading:
            self.heading_text.append(data)
        elif self.current_cell is not None:
            self.current_cell.append(data)

    def handle_endtag(self, tag):
        if tag == "h2":
            self.current_heading = "".join(self.heading_text).strip()
            self.in_heading = False
            self.heading_text = []
        elif tag in {"th", "td"} and self.current_cell is not None:
            cell_text = " ".join(part.strip() for part in self.current_cell if part.strip()).strip()
            self.current_row.append(cell_text)
            self.current_cell = None
        elif tag == "tr" and self.in_table and self.current_row is not None:
            if any(cell for cell in self.current_row):
                self.current_table.append(self.current_row)
            self.current_row = None
        elif tag == "table" and self.in_table:
            if self.current_heading and self.current_table:
                self.tables[self.current_heading] = self.current_table
            self.in_table = False
            self.current_table = None


def parse_numeric(value):
    """Convert HTML cell text to int/float when possible."""
    if value in {"", "N/A"}:
        return "N/A"
    try:
        number = float(value)
    except ValueError:
        return value
    if number.is_integer():
        return int(number)
    return number

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

def parse_qc_json(path):
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
        r["PBC1"] = safe_get(q, "lib_complexity", "lib_complexity", rep, "PBC1")
        r["PBC2"] = safe_get(q, "lib_complexity", "lib_complexity", rep, "PBC2")

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


def table_to_dict(table):
    """Convert an HTML table to a row-label -> column-label mapping."""
    if not table or len(table[0]) < 2:
        return {}
    headers = table[0][1:]
    rows = {}
    for row in table[1:]:
        if not row:
            continue
        label = row[0]
        values = {}
        for idx, header in enumerate(headers, start=1):
            values[header] = parse_numeric(row[idx]) if idx < len(row) else "N/A"
        rows[label] = values
    return rows


def safe_table_get(table_map, row_label, column_label, default="N/A", digits=4):
    value = table_map.get(row_label, {}).get(column_label, default)
    if isinstance(value, float):
        return round(value, digits)
    return value


def parse_qc_html(path):
    """Parse an ENCODE QC HTML report and return the same flat metrics dict as the JSON parser."""
    parser = QCHtmlParser()
    parser.feed(path.read_text(encoding="utf-8"))

    raw_unfiltered = table_to_dict(parser.tables.get("SAMstat (raw unfiltered BAM)", []))
    dedup_metrics = table_to_dict(parser.tables.get("Marking duplicates (filtered BAM)", []))
    mito_metrics = table_to_dict(parser.tables.get("Fraction of mitochondrial reads (unfiltered BAM)", []))
    lib_complexity = table_to_dict(parser.tables.get("Library complexity (filtered non-mito BAM)", []))
    tss_metrics = table_to_dict(parser.tables.get("TSS enrichment (filtered/deduped BAM)", []))
    reproducibility = table_to_dict(parser.tables.get("Reproducibility QC and peak detection statistics", []))
    frip_metrics = table_to_dict(parser.tables.get("Fraction of reads in peaks (FRiP)", []))

    result = {}
    for rep in ("rep1", "rep2"):
        result[rep] = {
            "mapped_reads": fmt_int(safe_table_get(raw_unfiltered, "Mapped Reads", rep)),
            "distinct_frags": fmt_int(safe_table_get(lib_complexity, "Distinct Fragments", rep)),
            "NRF": safe_table_get(lib_complexity, "NRF = Distinct/Total", rep),
            "PBC1": safe_table_get(lib_complexity, "PBC1 = OneRead/Distinct", rep),
            "PBC2": safe_table_get(lib_complexity, "PBC2 = OneRead/TwoRead", rep),
            "TSS_enrichment": safe_table_get(tss_metrics, "TSS enrichment", rep, digits=2),
            "pct_dup": safe_table_get(dedup_metrics, "% Duplicate Reads", rep),
            "frac_mito": safe_table_get(
                mito_metrics,
                "Rm/(Rn+Rm) = Frac. of mitochondrial reads",
                rep,
            ),
        }

    result["opt_peaks"] = fmt_int(safe_table_get(reproducibility, "N optimal", "idr"))
    result["cons_peaks"] = fmt_int(safe_table_get(reproducibility, "N conservative", "idr"))
    result["reproducibility"] = safe_table_get(reproducibility, "Rescue Ratio", "idr")
    result["frip_optimal"] = safe_table_get(frip_metrics, "Fraction of Reads in Peaks", "pooled-pr1_vs_pooled-pr2")
    result["frip_conservative"] = safe_table_get(frip_metrics, "Fraction of Reads in Peaks", "rep1_vs_rep2")

    return result


def parse_qc(path):
    """Parse a QC report from either HTML or JSON input."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".json":
        return parse_qc_json(path)
    if suffix in {".html", ".htm"}:
        return parse_qc_html(path)
    raise ValueError(f"Unsupported QC input format: {path}")

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


def print_concise_summary(paths):
    """Print a quick terminal summary similar to the old fast QC script."""
    for arg_key, sample_name in SUMMARY_LABELS.items():
        path = paths.get(arg_key)
        print(f"\n=== {sample_name} ===")
        if not path:
            print("  Missing input path")
            continue
        if not Path(path).exists():
            print(f"  Missing file: {path}")
            continue
        try:
            metrics = parse_qc(path)
        except Exception as exc:
            print(f"  Failed to parse {path}: {exc}")
            continue

        print(f"  TSS Rep1: {metrics['rep1'].get('TSS_enrichment', 'N/A')}")
        print(f"  TSS Rep2: {metrics['rep2'].get('TSS_enrichment', 'N/A')}")
        print(f"  FRiP (pooled): {metrics.get('frip_optimal', 'N/A')}")
        print(f"  IDR Optimal Peaks: {metrics.get('opt_peaks', 'N/A')}")
        print(f"  IDR Conservative Peaks: {metrics.get('cons_peaks', 'N/A')}")
        for rep in ("rep1", "rep2"):
            rep_metrics = metrics.get(rep, {})
            print(
                f"  {rep}: NRF={rep_metrics.get('NRF', 'N/A')} "
                f"PBC1={rep_metrics.get('PBC1', 'N/A')} "
                f"PBC2={rep_metrics.get('PBC2', 'N/A')}"
            )

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
    parser.add_argument("--human-adrenal", dest="human_adrenal", type=Path, metavar="QC_REPORT", default=DEFAULT_INPUTS["human_adrenal"])
    parser.add_argument("--mouse-adrenal", dest="mouse_adrenal", type=Path, metavar="QC_REPORT", default=DEFAULT_INPUTS["mouse_adrenal"])
    parser.add_argument("--human-uterus",  dest="human_uterus",  type=Path, metavar="QC_REPORT", default=DEFAULT_INPUTS["human_uterus"])
    parser.add_argument("--mouse-uterus",  dest="mouse_uterus",  type=Path, metavar="QC_REPORT", default=DEFAULT_INPUTS["mouse_uterus"])
    parser.add_argument("--out-tsv", type=Path, default=QC_RESULTS_DIR / "qc_summary_table.tsv")
    parser.add_argument("--out-md",  type=Path, default=QC_RESULTS_DIR / "qc_summary_table.md")
    parser.add_argument(
        "--print-summary",
        action="store_true",
        help="Print a concise per-sample QC summary to the terminal.",
    )
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
    if args.print_summary:
        print_concise_summary(paths)

if __name__ == "__main__":
    main()
