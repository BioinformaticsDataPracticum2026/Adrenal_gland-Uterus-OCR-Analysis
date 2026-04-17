#!/bin/bash

module load bedtools

# HUMAN ADRENAL
echo "Processing HUMAN adrenal..."

cd /ocean/projects/bio230007p/ychen120/HumanAtac/AdrenalGland/idr_reproducibility

zcat idr.optimal_peak.narrowPeak.gz | cut -f1-3 > human_adrenal_optimal_peaks.bed

bedtools slop \
-i /ocean/projects/bio230007p/ikaplow/HumanGenomeInfo/gencode.v27.annotation.protTranscript.TSSsWithStrand_sorted.bed \
-g /ocean/projects/bio230007p/ychen120/genome_files/hg38.fa.fai \
-b 2000 \
> human_promoters_2kb.bed

bedtools intersect -u \
-a human_adrenal_optimal_peaks.bed \
-b human_promoters_2kb.bed \
> human_adrenal_promoter_peaks.bed

bedtools intersect -v \
-a human_adrenal_optimal_peaks.bed \
-b human_promoters_2kb.bed \
> human_adrenal_enhancer_peaks.bed

echo "Human results:"
wc -l human_adrenal_optimal_peaks.bed
wc -l human_adrenal_promoter_peaks.bed
wc -l human_adrenal_enhancer_peaks.bed

# MOUSE ADRENAL
echo "Processing MOUSE adrenal..."

cd /ocean/projects/bio230007p/ychen120/MouseAtac/AdrenalGland/peak/idr_reproducibility

zcat idr.optimal_peak.narrowPeak.gz | cut -f1-3 > mouse_adrenal_optimal_peaks.bed

bedtools slop \
-i /ocean/projects/bio230007p/ikaplow/MouseGenomeInfo/gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed \
-g /ocean/projects/bio230007p/ychen120/genome_files/mm10.fa.fai \
-b 2000 \
> mouse_promoters_2kb.bed

bedtools intersect -u \
-a mouse_adrenal_optimal_peaks.bed \
-b mouse_promoters_2kb.bed \
> mouse_adrenal_promoter_peaks.bed

bedtools intersect -v \
-a mouse_adrenal_optimal_peaks.bed \
-b mouse_promoters_2kb.bed \
> mouse_adrenal_enhancer_peaks.bed

echo "Mouse results:"
wc -l mouse_adrenal_optimal_peaks.bed
wc -l mouse_adrenal_promoter_peaks.bed
wc -l mouse_adrenal_enhancer_peaks.bed

echo "Done!"
