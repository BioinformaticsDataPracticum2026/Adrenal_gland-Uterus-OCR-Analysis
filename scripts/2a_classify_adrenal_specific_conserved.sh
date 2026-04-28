#!/bin/bash
# Classify Mouse->Human HALPER adrenal OCRs as conserved or specific by
# overlap with native human adrenal OCRs.
#
# Definitions:
#   Conserved : lifted mouse peak overlaps a human adrenal peak
#   Specific  : lifted mouse peak does not overlap a human adrenal peak
#
# Usage: bash scripts/2a_classify_adrenal_specific_conserved.sh
# Requires: bedtools

set -euo pipefail

if command -v module >/dev/null 2>&1; then
    module load bedtools >/dev/null 2>&1 || true
fi

if ! command -v bedtools >/dev/null 2>&1; then
    echo "Error: bedtools is not available on PATH." >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

INPUT_DIR="$REPO_ROOT/data/idr_Optimal_Peaks"
OUTPUT_DIR="$REPO_ROOT/results/Specific_and_Conserved"
SUMMARY_TXT="$OUTPUT_DIR/adrenal_specific_conserved_summary.txt"

MOUSE_M2H="$INPUT_DIR/Mouse_AdrenalGland_idr.optimal_peak.MouseToHuman.HALPER.narrowPeak"
HUMAN_PEAKS="$INPUT_DIR/Human_AdrenalGland_idr.optimal_peak.narrowPeak"

CONSERVED_BED="$OUTPUT_DIR/Mouse_AdrenalGland_M2H_conserved_with_Human_AdrenalGland.bed"
SPECIFIC_BED="$OUTPUT_DIR/Mouse_AdrenalGland_M2H_specific_no_Human_AdrenalGland_overlap.bed"

for f in "$MOUSE_M2H" "$HUMAN_PEAKS"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: missing required file: $f" >&2
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"

bedtools intersect \
    -a "$MOUSE_M2H" \
    -b "$HUMAN_PEAKS" \
    -u > "$CONSERVED_BED"

bedtools intersect \
    -a "$MOUSE_M2H" \
    -b "$HUMAN_PEAKS" \
    -v > "$SPECIFIC_BED"

TOTAL_MOUSE_M2H=$(wc -l < "$MOUSE_M2H")
TOTAL_HUMAN=$(wc -l < "$HUMAN_PEAKS")
N_CONSERVED=$(wc -l < "$CONSERVED_BED")
N_SPECIFIC=$(wc -l < "$SPECIFIC_BED")

awk -v total="$TOTAL_MOUSE_M2H" \
    -v conserved="$N_CONSERVED" \
    -v specific="$N_SPECIFIC" \
    -v human="$TOTAL_HUMAN" \
    -v m2h="$MOUSE_M2H" \
    -v human_file="$HUMAN_PEAKS" \
    -v conserved_bed="$CONSERVED_BED" \
    -v specific_bed="$SPECIFIC_BED" \
    'BEGIN {
        conserved_pct = (total > 0) ? conserved / total * 100 : 0;
        specific_pct = (total > 0) ? specific / total * 100 : 0;
        print "Mouse adrenal Mouse->Human HALPER vs human adrenal OCR summary";
        print "================================================================";
        print "";
        print "Input files:";
        print "  Mouse M2H HALPER: " m2h;
        print "  Human peaks:      " human_file;
        print "";
        printf "Total Mouse->Human HALPER peaks: %d\n", total;
        printf "Total human adrenal peaks:       %d\n", human;
        printf "Conserved peaks:                 %d (%.2f%%)\n", conserved, conserved_pct;
        printf "Specific peaks:                  %d (%.2f%%)\n", specific, specific_pct;
        print "";
        print "Output files:";
        print "  Conserved: " conserved_bed;
        print "  Specific:  " specific_bed;
    }' > "$SUMMARY_TXT"

cat "$SUMMARY_TXT"
