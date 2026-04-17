#!/usr/bin/env bash

mkdir -p logs

HALPER_SCRIPT=/jet/home/xli51/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh
HAL=/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal

DATADIR=/ocean/projects/bio230007p/xli51/repo/Adrenal_gland-Uterus-OCR-Analysis/data/Promoters_and_Enhancers
OUTDIR=$DATADIR

########################################
# HUMAN → MOUSE
########################################

# enhancer
sbatch -A bio230007p -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b ${DATADIR}/human_adrenal_enhancer_peaks.bed \
  -o ${OUTDIR}/human_adrenal_enhancer_to_mouse \
  -s Human \
  -t Mouse \
  -c $HAL

# promoter
sbatch -A bio230007p -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b ${DATADIR}/human_adrenal_promoter_peaks.bed \
  -o ${OUTDIR}/human_adrenal_promoter_to_mouse \
  -s Human \
  -t Mouse \
  -c $HAL

########################################
# MOUSE → HUMAN
########################################

# enhancer
sbatch -A bio230007p -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b ${DATADIR}/mouse_adrenal_enhancer_peaks.bed \
  -o ${OUTDIR}/mouse_adrenal_enhancer_to_human \
  -s Mouse \
  -t Human \
  -c $HAL

# promoter
sbatch -A bio230007p -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b ${DATADIR}/mouse_adrenal_promoter_peaks.bed \
  -o ${OUTDIR}/mouse_adrenal_promoter_to_human \
  -s Mouse \
  -t Human \
  -c $HAL

echo "All jobs submitted."
