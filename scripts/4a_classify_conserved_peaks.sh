#!/bin/bash
# Classifies ATAC-seq peaks as conserved or species-specific using
# one-directional Mouse→Human HALPER liftover overlap.
#
# Definitions:
#   Mouse-conserved            : mouse peak lifts to human genome AND overlaps
#                                a native human ATAC-seq peak in the same tissue
#   Mouse-specific             : mouse peak that does not lift over, OR lifts over
#                                but finds no overlapping human peak
#   Human w/ mouse ortholog    : human peak overlapped by any M→H HALPER entry
#   Candidate human-specific   : human peak with no M→H overlap
#                                (weaker call; noted as a limitation without H→M data)
#
# Usage: bash scripts/4a_classify_conserved_peaks.sh
# Requires: bedtools >= 2.30

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

DATA_DIR="$REPO_ROOT/data/Promoters_and_Enhancers"
RESULTS_DIR="$REPO_ROOT/results/Enhancer_and_Promoters"
CONSERVED_DIR="$RESULTS_DIR/Conserved"
SPECIFIC_DIR="$RESULTS_DIR/Specific"
TMP_DIR="$RESULTS_DIR/tmp"
SUMMARY_TXT="$RESULTS_DIR/conserved_specific_summary.txt"

mkdir -p "$CONSERVED_DIR" "$SPECIFIC_DIR" "$TMP_DIR"
rm -f "$SUMMARY_TXT"

exec > >(tee "$SUMMARY_TXT") 2>&1

for FEATURE in enhancer promoter; do
    echo "========================================"
    echo "  Feature: adrenal $FEATURE"
    echo "========================================"

    M2H="$DATA_DIR/mouse_adrenal_${FEATURE}_peaks.MouseToHuman.HALPER.narrowPeak"
    HUMAN_PEAKS="$DATA_DIR/human_adrenal_${FEATURE}_peaks.bed"
    MOUSE_PEAKS="$DATA_DIR/mouse_adrenal_${FEATURE}_peaks.bed"

    for f in "$M2H" "$HUMAN_PEAKS" "$MOUSE_PEAKS"; do
        if [[ ! -f "$f" ]]; then
            echo "  [SKIP] Missing file: $f"
            continue 2
        fi
    done

    # ------------------------------------------------------------------ #
    # Step 1: Assign coordinate-based IDs to peak sets                    #
    #         Format: chr:start-end                                        #
    # ------------------------------------------------------------------ #
    awk 'OFS="\t" {$4=$1":"$2"-"$3; print}' "$HUMAN_PEAKS" \
        > "$TMP_DIR/${FEATURE}_human_named.bed"
    awk 'OFS="\t" {$4=$1":"$2"-"$3; print}' "$MOUSE_PEAKS" \
        > "$TMP_DIR/${FEATURE}_mouse_named.bed"

    # Strip summit offset from HALPER name field
    # e.g. "chr3:45000000-45001000:532" -> "chr3:45000000-45001000"
    awk 'OFS="\t" {sub(/:[0-9]+$/, "", $4); print}' "$M2H" \
        > "$TMP_DIR/${FEATURE}_M2H_named.bed"

    # ------------------------------------------------------------------ #
    # Step 2: M→H overlap                                                 #
    #         Mouse peaks lifted to human coords vs. human native peaks    #
    #         $4 = mouse source peak ID  |  $14 = human native peak ID    #
    # ------------------------------------------------------------------ #
    bedtools intersect \
        -a "$TMP_DIR/${FEATURE}_M2H_named.bed" \
        -b "$TMP_DIR/${FEATURE}_human_named.bed" \
        -wa -wb > "$TMP_DIR/${FEATURE}_M2H_overlaps.bed"

    # ------------------------------------------------------------------ #
    # Step 3: Extract peak IDs from overlaps                              #
    # ------------------------------------------------------------------ #
    cut -f4  "$TMP_DIR/${FEATURE}_M2H_overlaps.bed" | sort -u \
        > "$TMP_DIR/${FEATURE}_conserved_mouse_ids.txt"
    cut -f14 "$TMP_DIR/${FEATURE}_M2H_overlaps.bed" | sort -u \
        > "$TMP_DIR/${FEATURE}_human_ortholog_ids.txt"

    # ------------------------------------------------------------------ #
    # Step 4: Write output files                                          #
    # ------------------------------------------------------------------ #
    awk 'NR==FNR {k[$1]=1; next} ($4 in k)' \
        "$TMP_DIR/${FEATURE}_conserved_mouse_ids.txt" \
        "$TMP_DIR/${FEATURE}_mouse_named.bed" \
        > "$CONSERVED_DIR/mouse_adrenal_${FEATURE}_conserved.bed"

    awk 'NR==FNR {k[$1]=1; next} !($4 in k)' \
        "$TMP_DIR/${FEATURE}_conserved_mouse_ids.txt" \
        "$TMP_DIR/${FEATURE}_mouse_named.bed" \
        > "$SPECIFIC_DIR/mouse_adrenal_${FEATURE}_specific.bed"

    awk 'NR==FNR {k[$1]=1; next} ($4 in k)' \
        "$TMP_DIR/${FEATURE}_human_ortholog_ids.txt" \
        "$TMP_DIR/${FEATURE}_human_named.bed" \
        > "$CONSERVED_DIR/human_adrenal_${FEATURE}_with_mouse_ortholog.bed"

    awk 'NR==FNR {k[$1]=1; next} !($4 in k)' \
        "$TMP_DIR/${FEATURE}_human_ortholog_ids.txt" \
        "$TMP_DIR/${FEATURE}_human_named.bed" \
        > "$SPECIFIC_DIR/human_adrenal_${FEATURE}_candidate_specific.bed"

    # ------------------------------------------------------------------ #
    # Step 5: Summary                                                      #
    # ------------------------------------------------------------------ #
    TOTAL_M=$(wc -l < "$MOUSE_PEAKS")
    TOTAL_H=$(wc -l < "$HUMAN_PEAKS")
    N_M2H=$(wc -l < "$TMP_DIR/${FEATURE}_M2H_named.bed")
    N_CONS_M=$(wc -l < "$CONSERVED_DIR/mouse_adrenal_${FEATURE}_conserved.bed")
    N_SPEC_M=$(wc -l < "$SPECIFIC_DIR/mouse_adrenal_${FEATURE}_specific.bed")
    N_H_ORTH=$(wc -l < "$CONSERVED_DIR/human_adrenal_${FEATURE}_with_mouse_ortholog.bed")
    N_H_SPEC=$(wc -l < "$SPECIFIC_DIR/human_adrenal_${FEATURE}_candidate_specific.bed")

    echo ""
    echo "  --- Mouse adrenal $FEATURE peaks ---"
    printf "  %-40s %6d\n" "Total:" "$TOTAL_M"
    printf "  %-40s %6d\n" "  Lifted over (M→H):" "$N_M2H"
    printf "  %-40s %6d  (%.1f%%)\n" "  Conserved (overlap human peak):" \
        "$N_CONS_M" "$(awk "BEGIN{printf \"%.1f\", $N_CONS_M/$TOTAL_M*100}")"
    printf "  %-40s %6d  (%.1f%%)\n" "  Mouse-specific:" \
        "$N_SPEC_M" "$(awk "BEGIN{printf \"%.1f\", $N_SPEC_M/$TOTAL_M*100}")"

    echo ""
    echo "  --- Human adrenal $FEATURE peaks (inferred from M→H direction) ---"
    printf "  %-40s %6d\n" "Total:" "$TOTAL_H"
    printf "  %-40s %6d  (%.1f%%)\n" "  With mouse ortholog:" \
        "$N_H_ORTH" "$(awk "BEGIN{printf \"%.1f\", $N_H_ORTH/$TOTAL_H*100}")"
    printf "  %-40s %6d  (%.1f%%)\n" "  Candidate human-specific:" \
        "$N_H_SPEC" "$(awk "BEGIN{printf \"%.1f\", $N_H_SPEC/$TOTAL_H*100}")"
    echo "  (Note: human classification is inferential without H→M liftover)"
    echo ""
done

echo "Conserved output written to: $CONSERVED_DIR"
echo "Specific output written to: $SPECIFIC_DIR"
rm -rf "$TMP_DIR"
echo "Temporary directory removed: $TMP_DIR"
echo "Summary written to: $SUMMARY_TXT"
echo "Done."
