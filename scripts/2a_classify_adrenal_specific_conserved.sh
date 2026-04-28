#!/bin/bash
# Classify adrenal OCRs as conserved or species-specific using one-directional
# Mouse->Human HALPER liftover overlap with native human adrenal OCRs.
#
# Definitions:
#   Mouse conserved            : source mouse peak has any lifted interval that
#                                overlaps a human adrenal peak; output in mm10
#   Mouse specific             : source mouse peak has no lifted interval that
#                                overlaps a human adrenal peak; output in mm10
#   Human with mouse ortholog  : human peak overlaps a lifted mouse peak
#   Candidate human specific   : human peak does not overlap any lifted mouse peak
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
MOUSE_PEAKS="$INPUT_DIR/Mouse_AdrenalGland_idr.optimal_peak.narrowPeak"
HUMAN_PEAKS="$INPUT_DIR/Human_AdrenalGland_idr.optimal_peak.narrowPeak"

MOUSE_CONSERVED_BED="$OUTPUT_DIR/mouse_coord_mm10.Mouse_AdrenalGland.conserved_by_MouseToHuman_overlap_human_adrenal.bed"
MOUSE_SPECIFIC_BED="$OUTPUT_DIR/mouse_coord_mm10.Mouse_AdrenalGland.specific_no_MouseToHuman_overlap_human_adrenal.bed"
HUMAN_CONSERVED_BED="$OUTPUT_DIR/human_coord_hg38.Human_AdrenalGland.conserved_overlap_mouse_MouseToHuman.bed"
HUMAN_SPECIFIC_BED="$OUTPUT_DIR/human_coord_hg38.Human_AdrenalGland.specific_no_overlap_mouse_MouseToHuman.bed"
TMP_DIR="$OUTPUT_DIR/tmp"

for f in "$MOUSE_M2H" "$MOUSE_PEAKS" "$HUMAN_PEAKS"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: missing required file: $f" >&2
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

bedtools intersect \
    -a "$MOUSE_M2H" \
    -b "$HUMAN_PEAKS" \
    -u > "$TMP_DIR/mouse_M2H_overlapping_human.bed"

awk 'OFS="\t" {$4=$1":"$2"-"$3; print}' "$MOUSE_PEAKS" \
    > "$TMP_DIR/mouse_named.bed"

awk 'OFS="\t" {sub(/:[0-9]+$/, "", $4); print $4}' "$TMP_DIR/mouse_M2H_overlapping_human.bed" \
    | sort -u > "$TMP_DIR/mouse_conserved_ids.txt"

awk 'FILENAME == ARGV[1] {k[$1]=1; next} FILENAME == ARGV[2] && ($4 in k)' \
    "$TMP_DIR/mouse_conserved_ids.txt" \
    "$TMP_DIR/mouse_named.bed" \
    > "$MOUSE_CONSERVED_BED"

awk 'FILENAME == ARGV[1] {k[$1]=1; next} FILENAME == ARGV[2] && !($4 in k)' \
    "$TMP_DIR/mouse_conserved_ids.txt" \
    "$TMP_DIR/mouse_named.bed" \
    > "$MOUSE_SPECIFIC_BED"

bedtools intersect \
    -a "$HUMAN_PEAKS" \
    -b "$MOUSE_M2H" \
    -u > "$HUMAN_CONSERVED_BED"

bedtools intersect \
    -a "$HUMAN_PEAKS" \
    -b "$MOUSE_M2H" \
    -v > "$HUMAN_SPECIFIC_BED"

TOTAL_MOUSE=$(wc -l < "$MOUSE_PEAKS")
TOTAL_MOUSE_M2H=$(wc -l < "$MOUSE_M2H")
TOTAL_HUMAN=$(wc -l < "$HUMAN_PEAKS")
N_MOUSE_CONSERVED=$(wc -l < "$MOUSE_CONSERVED_BED")
N_MOUSE_SPECIFIC=$(wc -l < "$MOUSE_SPECIFIC_BED")
N_HUMAN_CONSERVED=$(wc -l < "$HUMAN_CONSERVED_BED")
N_HUMAN_SPECIFIC=$(wc -l < "$HUMAN_SPECIFIC_BED")

awk -v total_mouse="$TOTAL_MOUSE" \
    -v total_m2h="$TOTAL_MOUSE_M2H" \
    -v mouse_conserved="$N_MOUSE_CONSERVED" \
    -v mouse_specific="$N_MOUSE_SPECIFIC" \
    -v human="$TOTAL_HUMAN" \
    -v human_conserved="$N_HUMAN_CONSERVED" \
    -v human_specific="$N_HUMAN_SPECIFIC" \
    -v m2h="$MOUSE_M2H" \
    -v mouse_file="$MOUSE_PEAKS" \
    -v human_file="$HUMAN_PEAKS" \
    -v mouse_conserved_bed="$MOUSE_CONSERVED_BED" \
    -v mouse_specific_bed="$MOUSE_SPECIFIC_BED" \
    -v human_conserved_bed="$HUMAN_CONSERVED_BED" \
    -v human_specific_bed="$HUMAN_SPECIFIC_BED" \
    'BEGIN {
        mouse_conserved_pct = (total_mouse > 0) ? mouse_conserved / total_mouse * 100 : 0;
        mouse_specific_pct = (total_mouse > 0) ? mouse_specific / total_mouse * 100 : 0;
        human_conserved_pct = (human > 0) ? human_conserved / human * 100 : 0;
        human_specific_pct = (human > 0) ? human_specific / human * 100 : 0;
        print "Mouse adrenal Mouse->Human HALPER vs human adrenal OCR summary";
        print "================================================================";
        print "";
        print "Input files:";
        print "  Mouse M2H HALPER: " m2h;
        print "  Mouse peaks:      " mouse_file;
        print "  Human peaks:      " human_file;
        print "";
        printf "Total mouse adrenal peaks:       %d\n", total_mouse;
        printf "Total Mouse->Human HALPER rows:  %d\n", total_m2h;
        printf "Total human adrenal peaks:       %d\n", human;
        print "";
        print "Mouse source peaks in mm10 coordinates:";
        printf "  Conserved:                     %d (%.2f%%)\n", mouse_conserved, mouse_conserved_pct;
        printf "  Mouse specific:                %d (%.2f%%)\n", mouse_specific, mouse_specific_pct;
        print "";
        print "Human peaks:";
        printf "  With mouse ortholog:           %d (%.2f%%)\n", human_conserved, human_conserved_pct;
        printf "  Candidate human specific:      %d (%.2f%%)\n", human_specific, human_specific_pct;
        print "";
        print "Output files:";
        print "  Mouse conserved:       " mouse_conserved_bed;
        print "  Mouse specific:        " mouse_specific_bed;
        print "  Human with ortholog:   " human_conserved_bed;
        print "  Human candidate specific: " human_specific_bed;
    }' > "$SUMMARY_TXT"

rm -rf "$TMP_DIR"
cat "$SUMMARY_TXT"
