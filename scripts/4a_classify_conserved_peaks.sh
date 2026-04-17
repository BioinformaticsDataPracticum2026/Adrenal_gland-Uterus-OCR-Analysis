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

HALPER_DIR="data/halper_peaks"
IDR_DIR="data/idr_Conservative_Peaks"
OUT_DIR="results/peak_classification"
TMP_DIR="$OUT_DIR/tmp"

mkdir -p "$OUT_DIR" "$TMP_DIR"

for TISSUE in AdrenalGland Uterus; do
    echo "========================================"
    echo "  Tissue: $TISSUE"
    echo "========================================"

    M2H="$HALPER_DIR/Mouse_${TISSUE}_idr.conservative_peak.MouseToHuman.HALPER.narrowPeak"
    HUMAN_IDR="$IDR_DIR/Human_${TISSUE}_idr.conservative_peak.narrowPeak"
    MOUSE_IDR="$IDR_DIR/Mouse_${TISSUE}_idr.conservative_peak.narrowPeak"

    for f in "$M2H" "$HUMAN_IDR" "$MOUSE_IDR"; do
        if [[ ! -f "$f" ]]; then
            echo "  [SKIP] Missing file: $f"
            continue 2
        fi
    done

    # ------------------------------------------------------------------ #
    # Step 1: Assign coordinate-based IDs to IDR peaks (col 4 = ".")     #
    #         Format: chr:start-end                                        #
    # ------------------------------------------------------------------ #
    awk 'OFS="\t" {$4=$1":"$2"-"$3; print}' "$HUMAN_IDR" \
        > "$TMP_DIR/${TISSUE}_human_named.narrowPeak"
    awk 'OFS="\t" {$4=$1":"$2"-"$3; print}' "$MOUSE_IDR" \
        > "$TMP_DIR/${TISSUE}_mouse_named.narrowPeak"

    # Strip summit offset from HALPER name field
    # e.g. "chr3:45000000-45001000:532" -> "chr3:45000000-45001000"
    awk 'OFS="\t" {sub(/:[0-9]+$/, "", $4); print}' "$M2H" \
        > "$TMP_DIR/${TISSUE}_M2H_named.narrowPeak"

    # ------------------------------------------------------------------ #
    # Step 2: M→H overlap                                                 #
    #         Mouse peaks in human coords vs. human native peaks           #
    #         $4 = mouse source peak ID  |  $14 = human native peak ID    #
    # ------------------------------------------------------------------ #
    bedtools intersect \
        -a "$TMP_DIR/${TISSUE}_M2H_named.narrowPeak" \
        -b "$TMP_DIR/${TISSUE}_human_named.narrowPeak" \
        -wa -wb > "$TMP_DIR/${TISSUE}_M2H_overlaps.bed"

    # ------------------------------------------------------------------ #
    # Step 3: Extract peak IDs from overlaps                              #
    # ------------------------------------------------------------------ #
    cut -f4  "$TMP_DIR/${TISSUE}_M2H_overlaps.bed" | sort -u \
        > "$TMP_DIR/${TISSUE}_conserved_mouse_ids.txt"
    cut -f14 "$TMP_DIR/${TISSUE}_M2H_overlaps.bed" | sort -u \
        > "$TMP_DIR/${TISSUE}_human_ortholog_ids.txt"

    # ------------------------------------------------------------------ #
    # Step 4: Write output narrowPeak files                               #
    # ------------------------------------------------------------------ #
    awk 'NR==FNR {k[$1]=1; next} ($4 in k)' \
        "$TMP_DIR/${TISSUE}_conserved_mouse_ids.txt" \
        "$TMP_DIR/${TISSUE}_mouse_named.narrowPeak" \
        > "$OUT_DIR/${TISSUE}_mouse_conserved.narrowPeak"

    awk 'NR==FNR {k[$1]=1; next} !($4 in k)' \
        "$TMP_DIR/${TISSUE}_conserved_mouse_ids.txt" \
        "$TMP_DIR/${TISSUE}_mouse_named.narrowPeak" \
        > "$OUT_DIR/${TISSUE}_mouse_specific.narrowPeak"

    awk 'NR==FNR {k[$1]=1; next} ($4 in k)' \
        "$TMP_DIR/${TISSUE}_human_ortholog_ids.txt" \
        "$TMP_DIR/${TISSUE}_human_named.narrowPeak" \
        > "$OUT_DIR/${TISSUE}_human_with_mouse_ortholog.narrowPeak"

    awk 'NR==FNR {k[$1]=1; next} !($4 in k)' \
        "$TMP_DIR/${TISSUE}_human_ortholog_ids.txt" \
        "$TMP_DIR/${TISSUE}_human_named.narrowPeak" \
        > "$OUT_DIR/${TISSUE}_human_candidate_specific.narrowPeak"

    # ------------------------------------------------------------------ #
    # Step 5: Summary                                                      #
    # ------------------------------------------------------------------ #
    TOTAL_M=$(wc -l < "$MOUSE_IDR")
    TOTAL_H=$(wc -l < "$HUMAN_IDR")
    N_M2H=$(wc -l < "$TMP_DIR/${TISSUE}_M2H_named.narrowPeak")
    N_CONS_M=$(wc -l < "$OUT_DIR/${TISSUE}_mouse_conserved.narrowPeak")
    N_SPEC_M=$(wc -l < "$OUT_DIR/${TISSUE}_mouse_specific.narrowPeak")
    N_H_ORTH=$(wc -l < "$OUT_DIR/${TISSUE}_human_with_mouse_ortholog.narrowPeak")
    N_H_SPEC=$(wc -l < "$OUT_DIR/${TISSUE}_human_candidate_specific.narrowPeak")

    echo ""
    echo "  --- Mouse peaks ---"
    printf "  %-40s %6d\n" "Total:" "$TOTAL_M"
    printf "  %-40s %6d\n" "  Lifted over (M→H):" "$N_M2H"
    printf "  %-40s %6d  (%.1f%%)\n" "  Conserved (overlap human peak):" \
        "$N_CONS_M" "$(awk "BEGIN{printf \"%.1f\", $N_CONS_M/$TOTAL_M*100}")"
    printf "  %-40s %6d  (%.1f%%)\n" "  Mouse-specific:" \
        "$N_SPEC_M" "$(awk "BEGIN{printf \"%.1f\", $N_SPEC_M/$TOTAL_M*100}")"

    echo ""
    echo "  --- Human peaks (inferred from M→H direction) ---"
    printf "  %-40s %6d\n" "Total:" "$TOTAL_H"
    printf "  %-40s %6d  (%.1f%%)\n" "  With mouse ortholog:" \
        "$N_H_ORTH" "$(awk "BEGIN{printf \"%.1f\", $N_H_ORTH/$TOTAL_H*100}")"
    printf "  %-40s %6d  (%.1f%%)\n" "  Candidate human-specific:" \
        "$N_H_SPEC" "$(awk "BEGIN{printf \"%.1f\", $N_H_SPEC/$TOTAL_H*100}")"
    echo "  (Note: human classification is inferential without H→M liftover)"
    echo ""
done

echo "Output written to: $OUT_DIR"
echo "Done."
