# Adrenal_gland-Uterus-OCR-Analysis

03-713: Bioinformatics Data Integration Practicum, Spring 2026

Team members: Xinyi L., Jimmy L., Yu Hsin C.

## Project Overview

todo

## Evaluate data quality



## Map open chromatin regions across species

Multi-species alignment: `/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal`

ATAC-seq:

```bash
# Human
/ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Human_AdrenalGland_idr.conservative_peak.narrowPeak
/ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Human_Uterus_idr.conservative_peak.narrowPeak

# Mouse
/ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Mouse_AdrenalGland_idr.conservative_peak.narrowPeak
/ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Mouse_Uterus_idr.conservative_peak.narrowPeak
```

Using https://github.com/pfenninglab/halLiftover-postprocessing/tree/master

```bash
# AG: Human to mouse
bash /jet/home/xli51/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Human_AdrenalGland_idr.conservative_peak.narrowPeak \
  -o AdrenalGland_halper_human_to_mouse \
  -s Human \
  -t Mouse \
  -c /ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal
    
# AG: Mouse to human
bash /jet/home/xli51/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Mouse_AdrenalGland_idr.conservative_peak.narrowPeak \
  -o AdrenalGland_halper_mouse_to_human \
  -s Mouse \
  -t Human \
  -c /ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal

# Uterus: Human to mouse
bash /jet/home/xli51/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Human_Uterus_idr.conservative_peak.narrowPeak \
  -o Uterus_halper_human_to_mouse \
  -s Human \
  -t Mouse \
  -c /ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal
  
# Uterus: Mouse to human
bash /jet/home/xli51/repos/halLiftover-postprocessing/halper_map_peak_orthologs.sh \
  -b /ocean/projects/bio230007p/xli51/data/idr_Conservative_Peaks/Mouse_Uterus_idr.conservative_peak.narrowPeak \
  -o Uterus_halper_mouse_to_human \
  -s Mouse \
  -t Human \
  -c /ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal
```

## Compare open chromatin between species



## Find biological processes that are likely to be regulated by open chromatin regions



## Compare candidate enhancers to candidate promoters



## Find transcription factors that tend to bind open chromatin regions

