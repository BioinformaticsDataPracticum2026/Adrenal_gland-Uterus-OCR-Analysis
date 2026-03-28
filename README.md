# Adrenal_gland-Uterus-OCR-Analysis

03-713: Bioinformatics Data Integration Practicum, Spring 2026

Team members: Xinyi L., Jimmy L., Yu Hsin C.

## Project Overview

todo

## Evaluate data quality



## Map open chromatin regions across species

Multi-species alignment: `/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal`

IDR Conservative peaks:

```bash
# Human
/ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Human_AdrenalGland_idr.conservative_peak.narrowPeak
/ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Human_Uterus_idr.conservative_peak.narrowPeak

# Mouse
/ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Mouse_AdrenalGland_idr.conservative_peak.narrowPeak
/ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Mouse_Uterus_idr.conservative_peak.narrowPeak
```

Using https://github.com/pfenninglab/halLiftover-postprocessing/tree/master

run_hal.sh: 

```bash
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
```

Output data: `/ocean/projects/bio230007p/xli51/data/hal_Mapper_Peaks`

## Compare open chromatin between species



## Find biological processes that are likely to be regulated by open chromatin regions



## Compare candidate enhancers to candidate promoters



## Find transcription factors that tend to bind open chromatin regions

