#!/usr/bin/env python3
"""Enhance ENCODE-style ATAC QC HTML reports with threshold-based coloring.

This script reads one or more QC HTML files, evaluates selected metrics against
thresholds, and writes out annotated HTML files with pass/fail colors.

Threshold config format (JSON):
{
  "metric name": {"op": "min", "value": 0.4, "label": ">= 0.4"},
  "another metric": {"op": "max", "value": 0.1},
  "boolean metric": {"op": "bool", "value": true}
}
"""

from __future__ import annotations

import argparse
import html
import json
import re
from pathlib import Path
from typing import Any
from urllib.parse import quote


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
QC_HTML_DIR = REPO_ROOT / "data" / "qc_html"
QC_RESULTS_DIR = REPO_ROOT / "results" / "qc"


DEFAULT_THRESHOLDS: dict[str, dict[str, Any]] = {
    "% Mapped Reads": {"op": "min", "value": 90.0, "label": ">= 90"},
    "% Properly Paired Reads": {"op": "min", "value": 80.0, "label": ">= 80"},
    "% Duplicate Reads": {"op": "max", "value": 30.0, "label": "<= 30"},
    "Rm/(Rn+Rm) = Frac. of mitochondrial reads": {
        "op": "max",
        "value": 0.1,
        "label": "<= 0.10",
    },
    "Fraction of reads in NFR": {"op": "min", "value": 0.4, "label": ">= 0.40"},
    "NFR / mono-nuc reads": {"op": "min", "value": 2.5, "label": ">= 2.50"},
    "Presence of NFR peak": {"op": "bool", "value": True, "label": "True"},
    "Presence of Mono-Nuc peak": {"op": "bool", "value": True, "label": "True"},
    "Presence of Di-Nuc peak": {"op": "bool", "value": True, "label": "True"},
    "Fraction of Reads in universal DHS regions": {
        "op": "min",
        "value": 0.3,
        "label": ">= 0.30",
    },
    "Fraction of Reads in blacklist regions": {
        "op": "max",
        "value": 0.01,
        "label": "<= 0.01",
    },
    "Fraction of Reads in promoter regions": {
        "op": "min",
        "value": 0.15,
        "label": ">= 0.15",
    },
    "Fraction of Reads in enhancer regions": {
        "op": "min",
        "value": 0.2,
        "label": ">= 0.20",
    },
    "TSS enrichment": {"op": "min", "value": 6.0, "label": ">= 6.0"},
}


ROW_RE = re.compile(r"<tr\b[^>]*>.*?</tr>", re.IGNORECASE | re.DOTALL)
CELL_RE = re.compile(r"<(th|td)\b([^>]*)>(.*?)</\1>", re.IGNORECASE | re.DOTALL)
TAG_RE = re.compile(r"<[^>]+>", re.DOTALL)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "inputs",
        nargs="*",
        help="Input QC HTML file(s) or directories.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="Optional JSON threshold config. Overrides defaults by metric name.",
    )
    parser.add_argument(
        "--suffix",
        default="_checked",
        help="Suffix added before .html for generated files.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Recursively scan input directories for *.html files.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=QC_RESULTS_DIR,
        help="Directory for generated HTML files.",
    )
    return parser.parse_args()


def load_thresholds(config_path: Path | None) -> dict[str, dict[str, Any]]:
    thresholds = dict(DEFAULT_THRESHOLDS)
    if not config_path:
        return thresholds
    with config_path.open("r", encoding="utf-8") as handle:
        custom = json.load(handle)
    for metric, spec in custom.items():
        thresholds[metric] = spec
    return thresholds


def strip_tags(text: str) -> str:
    text = TAG_RE.sub("", text)
    return html.unescape(text).strip()


def parse_scalar(value: str) -> Any:
    raw = value.strip()
    low = raw.lower()
    if low == "true":
        return True
    if low == "false":
        return False
    try:
        return float(raw)
    except ValueError:
        return raw


def evaluate(metric: str, value: str, thresholds: dict[str, dict[str, Any]]) -> tuple[bool, str] | None:
    spec = thresholds.get(metric)
    if not spec and metric.endswith("(QC pass)"):
        parsed = parse_scalar(value)
        if isinstance(parsed, bool):
            return parsed, "existing QC pass flag"
        return None
    if not spec:
        return None

    parsed = parse_scalar(value)
    op = spec["op"]
    target = spec["value"]
    label = spec.get("label", f"{op} {target}")

    if op == "min":
        passed = isinstance(parsed, (int, float)) and parsed >= float(target)
    elif op == "max":
        passed = isinstance(parsed, (int, float)) and parsed <= float(target)
    elif op == "bool":
        passed = parsed is bool(target)
    else:
        raise ValueError(f"Unsupported op for metric {metric!r}: {op!r}")
    return passed, label


def style_attr(existing: str, passed: bool, note: str) -> str:
    clean = existing.strip()
    if clean and not clean.startswith(" "):
        clean = " " + clean
    color = "#e8f1fb" if passed else "#fff4db"
    border = "#1d4ed8" if passed else "#b45309"
    title = html.escape(f"Threshold: {note}")
    return (
        f'{clean} class="qc-threshold {"qc-pass" if passed else "qc-fail"}" '
        f'style="background:{color};box-shadow: inset 4px 0 0 {border};font-weight:600;" '
        f'title="{title}"'
    )


def row_anchor(metric: str) -> str:
    return "qc-metric-" + quote(metric, safe="").lower()


def annotate_row(row_html: str, thresholds: dict[str, dict[str, Any]], stats: dict[str, Any], headers: list[str]) -> str:
    cells = list(CELL_RE.finditer(row_html))
    if len(cells) < 2:
        return row_html

    metric = strip_tags(cells[0].group(3))
    replacements: list[tuple[int, int, str]] = []

    for index, cell in enumerate(cells[1:], start=1):
        tag, attrs, inner = cell.group(1), cell.group(2), cell.group(3)
        if tag.lower() != "td":
            continue
        result = evaluate(metric, strip_tags(inner), thresholds)
        if result is None:
            continue
        passed, note = result
        stats["pass" if passed else "fail"] += 1
        if passed:
            stats["pass_terms"].add(metric)
        else:
            column_name = headers[index - 1] if index - 1 < len(headers) else f"col{index}"
            stats["fail_terms"].setdefault(metric, set()).add(column_name)
        new_cell = f"<{tag}{style_attr(attrs, passed, note)}>{inner}</{tag}>"
        replacements.append((cell.start(), cell.end(), new_cell))

    if not replacements:
        return row_html

    parts: list[str] = []
    cursor = 0
    for start, end, replacement in replacements:
        parts.append(row_html[cursor:start])
        parts.append(replacement)
        cursor = end
    parts.append(row_html[cursor:])
    updated = "".join(parts)
    anchor = row_anchor(metric)
    return updated.replace("<tr", f'<tr id="{anchor}"', 1)


def build_fail_list(title: str, items: dict[str, set[str]], css_class: str) -> str:
    if not items:
        body = "<li>None</li>"
    else:
        rows = []
        for metric in sorted(items):
            failures = ", ".join(f"{column} failed" for column in sorted(items[metric]))
            rows.append(
                f'<li><a href="#{row_anchor(metric)}"><code>{html.escape(metric)}</code></a>: {html.escape(failures)}</li>'
            )
        body = "".join(rows)
    return f"""
  <details class="{css_class}">
    <summary>{html.escape(title)}</summary>
    <ul>{body}</ul>
  </details>
"""


def build_metric_list(title: str, items: set[str], css_class: str) -> str:
    if not items:
        body = "<li>None</li>"
    else:
        body = "".join(
            f'<li><a href="#{row_anchor(item)}"><code>{html.escape(item)}</code></a></li>'
            for item in sorted(items)
        )
    return f"""
  <details class="{css_class}">
    <summary>{html.escape(title)}</summary>
    <ul>{body}</ul>
  </details>
"""


def build_header(stats: dict[str, Any], thresholds: dict[str, dict[str, Any]]) -> str:
    items = "".join(
        f"<li><code>{html.escape(metric)}</code>: {html.escape(str(spec.get('label', spec['value'])))}</li>"
        for metric, spec in thresholds.items()
    )
    pass_terms = build_metric_list("Metrics with passing results", stats["pass_terms"], "qc-pass-list")
    fail_terms = build_fail_list("Metrics with failing results", stats["fail_terms"], "qc-fail-list")
    return f"""
<style>
.qc-enhancer {{
  margin: 16px 0;
  padding: 14px 16px;
  border: 1px solid #cfd8dc;
  background: #f8fafc;
  font-family: Arial, sans-serif;
}}
.qc-enhancer h2 {{
  margin: 0 0 8px 0;
  font-size: 18px;
}}
.qc-enhancer p {{
  margin: 6px 0;
}}
.qc-enhancer .legend {{
  display: inline-block;
  margin-right: 12px;
  padding: 4px 8px;
  border-radius: 6px;
}}
.qc-enhancer .legend-pass {{
  background: #e8f1fb;
  color: #1e3a8a;
}}
.qc-enhancer .legend-fail {{
  background: #fff4db;
  color: #92400e;
}}
.qc-enhancer details {{
  margin-top: 10px;
}}
.qc-enhancer li {{
  margin: 4px 0;
}}
.qc-enhancer a {{
  color: inherit;
  text-decoration: underline;
}}
.qc-enhancer a:hover {{
  text-decoration-thickness: 2px;
}}
tr:target td,
tr:target th {{
  outline: 3px solid #0f172a;
  outline-offset: -3px;
}}
</style>
<div class="qc-enhancer">
  <h2>Threshold Review</h2>
  <p>
    <span class="legend legend-pass">Pass: {stats["pass"]}</span>
    <span class="legend legend-fail">Fail: {stats["fail"]}</span>
  </p>
  <p>Known metrics are color-coded directly in the tables below.</p>
{pass_terms}
{fail_terms}
  <details>
    <summary>Configured thresholds</summary>
    <ul>{items}</ul>
  </details>
</div>
"""


def inject_header(html_text: str, header: str) -> str:
    marker = "<h1>QC Report</h1>"
    if marker in html_text:
        return html_text.replace(marker, marker + header, 1)
    return header + html_text


def enhance_report(path: Path, thresholds: dict[str, dict[str, Any]], suffix: str, outdir: Path) -> Path:
    source = path.read_text(encoding="utf-8")
    stats = {"pass": 0, "fail": 0, "pass_terms": set(), "fail_terms": {}}

    state = {"headers": []}

    def repl(match: re.Match[str]) -> str:
        row_html = match.group(0)
        cells = list(CELL_RE.finditer(row_html))
        if cells and all(cell.group(1).lower() == "th" for cell in cells):
            labels = [strip_tags(cell.group(3)) for cell in cells[1:]]
            state["headers"] = [label for label in labels if label]
            return row_html
        return annotate_row(row_html, thresholds, stats, state["headers"])

    annotated = ROW_RE.sub(repl, source)
    annotated = inject_header(annotated, build_header(stats, thresholds))

    outdir.mkdir(parents=True, exist_ok=True)
    output = outdir / f"{path.stem}{suffix}{path.suffix}"
    output.write_text(annotated, encoding="utf-8")
    return output


def collect_files(inputs: list[str], recursive: bool, suffix: str) -> list[Path]:
    files: list[Path] = []
    seen: set[Path] = set()
    for item in inputs:
        path = Path(item)
        if path.is_dir():
            iterator = path.rglob("*.html") if recursive else path.glob("*.html")
            for candidate in iterator:
                if candidate.stem.endswith(suffix):
                    continue
                resolved = candidate.resolve()
                if resolved not in seen:
                    seen.add(resolved)
                    files.append(candidate)
        elif path.is_file():
            if path.stem.endswith(suffix):
                continue
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                files.append(path)
    return sorted(files)


def main() -> int:
    args = parse_args()
    thresholds = load_thresholds(args.config)
    inputs = args.inputs or [str(QC_HTML_DIR)]
    files = collect_files(inputs, args.recursive, args.suffix)
    if not files:
        raise SystemExit("No input HTML files found.")

    for file_path in files:
        output = enhance_report(file_path, thresholds, args.suffix, args.outdir)
        print(f"Enhanced: {file_path} -> {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
