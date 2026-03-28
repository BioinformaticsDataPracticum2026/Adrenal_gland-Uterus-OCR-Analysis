#!/bin/bash

mkdir -p logs

HALPER_SCRIPT=/jet/home/xli51/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh
HAL=/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal

# Adrenal Gland
# Human → Mouse
sbatch -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b /ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Human_AdrenalGland_idr.conservative_peak.narrowPeak \
  -o AdrenalGland_halper_human_to_mouse \
  -s Human \
  -t Mouse \
  -c $HAL

# Mouse → Human
sbatch -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b /ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Mouse_AdrenalGland_idr.conservative_peak.narrowPeak \
  -o AdrenalGland_halper_mouse_to_human \
  -s Mouse \
  -t Human \
  -c $HAL

# Uterus
# Human → Mouse
sbatch -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b /ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Human_Uterus_idr.conservative_peak.narrowPeak \
  -o Uterus_halper_human_to_mouse \
  -s Human \
  -t Mouse \
  -c $HAL

# Mouse → Human
sbatch -p RM-shared --mem=2000 -t 08:00:00 $HALPER_SCRIPT \
  -b /ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Mouse_Uterus_idr.conservative_peak.narrowPeak \
  -o Uterus_halper_mouse_to_human \
  -s Mouse \
  -t Human \
  -c $HAL

echo "All jobs submitted."