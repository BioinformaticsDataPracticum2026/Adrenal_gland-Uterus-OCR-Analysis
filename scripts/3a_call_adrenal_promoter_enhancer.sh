#!/bin/bash

set -euo pipefail

module load bedtools

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

INPUT_DIR="$REPO_ROOT/data/idr_Optimal_Peaks"
OUTPUT_DIR="$REPO_ROOT/data/Promoters_and_Enhancers"

HUMAN_INPUT="$INPUT_DIR/Human_AdrenalGland_idr.optimal_peak.narrowPeak"
MOUSE_INPUT="$INPUT_DIR/Mouse_AdrenalGland_idr.optimal_peak.narrowPeak"

HUMAN_TSS="${HUMAN_TSS:-/ocean/projects/bio230007p/ikaplow/HumanGenomeInfo/gencode.v27.annotation.protTranscript.TSSsWithStrand_sorted.bed}"
HUMAN_GENOME_INDEX="${HUMAN_GENOME_INDEX:-/ocean/projects/bio230007p/xli51/data/HumanGenomeInfo/hg38.fa.fai}"
MOUSE_TSS="${MOUSE_TSS:-/ocean/projects/bio230007p/ikaplow/MouseGenomeInfo/gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed}"
MOUSE_GENOME_INDEX="${MOUSE_GENOME_INDEX:-/ocean/projects/bio230007p/xli51/data/MouseGenomeInfo/mm10.fa.fai}"

mkdir -p "$OUTPUT_DIR"

for f in "$HUMAN_INPUT" "$MOUSE_INPUT" "$HUMAN_TSS" "$HUMAN_GENOME_INDEX" "$MOUSE_TSS" "$MOUSE_GENOME_INDEX"; do
  if [[ ! -f "$f" ]]; then
    echo "Missing required file: $f" >&2
    exit 1
  fi
done

# HUMAN ADRENAL
echo "Processing HUMAN adrenal..."

cut -f1-3 "$HUMAN_INPUT" > "$OUTPUT_DIR/human_adrenal_optimal_peaks.bed"

bedtools slop \
  -i "$HUMAN_TSS" \
  -g "$HUMAN_GENOME_INDEX" \
  -b 2000 \
  > "$OUTPUT_DIR/human_promoters_2kb.bed"

bedtools intersect -u \
  -a "$OUTPUT_DIR/human_adrenal_optimal_peaks.bed" \
  -b "$OUTPUT_DIR/human_promoters_2kb.bed" \
  > "$OUTPUT_DIR/human_adrenal_promoter_peaks.bed"

bedtools intersect -v \
  -a "$OUTPUT_DIR/human_adrenal_optimal_peaks.bed" \
  -b "$OUTPUT_DIR/human_promoters_2kb.bed" \
  > "$OUTPUT_DIR/human_adrenal_enhancer_peaks.bed"

echo "Human results:"
wc -l "$OUTPUT_DIR/human_adrenal_optimal_peaks.bed"
wc -l "$OUTPUT_DIR/human_adrenal_promoter_peaks.bed"
wc -l "$OUTPUT_DIR/human_adrenal_enhancer_peaks.bed"

# MOUSE ADRENAL
echo "Processing MOUSE adrenal..."

cut -f1-3 "$MOUSE_INPUT" > "$OUTPUT_DIR/mouse_adrenal_optimal_peaks.bed"

bedtools slop \
  -i "$MOUSE_TSS" \
  -g "$MOUSE_GENOME_INDEX" \
  -b 2000 \
  > "$OUTPUT_DIR/mouse_promoters_2kb.bed"

bedtools intersect -u \
  -a "$OUTPUT_DIR/mouse_adrenal_optimal_peaks.bed" \
  -b "$OUTPUT_DIR/mouse_promoters_2kb.bed" \
  > "$OUTPUT_DIR/mouse_adrenal_promoter_peaks.bed"

bedtools intersect -v \
  -a "$OUTPUT_DIR/mouse_adrenal_optimal_peaks.bed" \
  -b "$OUTPUT_DIR/mouse_promoters_2kb.bed" \
  > "$OUTPUT_DIR/mouse_adrenal_enhancer_peaks.bed"

echo "Mouse results:"
wc -l "$OUTPUT_DIR/mouse_adrenal_optimal_peaks.bed"
wc -l "$OUTPUT_DIR/mouse_adrenal_promoter_peaks.bed"
wc -l "$OUTPUT_DIR/mouse_adrenal_enhancer_peaks.bed"

echo "Done!"
