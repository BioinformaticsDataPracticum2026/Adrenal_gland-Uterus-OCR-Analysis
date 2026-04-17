# Adrenal_gland-Uterus-OCR-Analysis

03-713: Bioinformatics Data Integration Practicum, Spring 2026

Team members: Xinyi L., Jimmy L., Yu Hsin C.

## Project Overview

This project compares open chromatin regions (OCRs) between adrenal gland and uterus in human and mouse ATAC-seq datasets. The current repository includes:

- QC summaries for all four datasets
- IDR conservative and optimal peak sets
- Cross-species OCR mapping with HAL/HALPER
- Local `rGREAT` functional enrichment analysis on optimal OCRs
- Promoter/enhancer partitioning for adrenal gland OCRs

Repository structure:

```text
data/
  qc/                         ENCODE-style QC HTML reports
  idr_Conservative_Peaks/     conservative peak sets used for initial mapping
  idr_Optimal_Peaks/          optimal peak sets used for downstream analyses
  hal_Mapper_Peaks/           HALPER-mapped OCR peak files
  Promoters_and_Enhancers/    promoter/enhancer split peak sets and lifted results
scripts/                      analysis scripts and job submission helpers
results/
  qc/                         QC tables and interpretation
  rGREAT/                     GO enrichment tables and dot plots
  promoters_enhancers.txt     adrenal promoter/enhancer summary
```

## Evaluate data quality

QC summaries are stored in `results/qc/`:

- `qc_summary_table.tsv` and `qc_summary_table.md`: combined metrics table
- `qc_summary.txt`: extracted TSS enrichment and pooled FRiP values
- `qc_interpretation.md`: short project-level interpretation

Relevant scripts:

- `scripts/extract_qc_metrics.py`: extracts key QC metrics from ENCODE ATAC-seq outputs
- `scripts/make_qc_table.py`: builds the final QC summary table
- `scripts/count_peaks.sh`: counts reproducible IDR conservative peaks

Summary of the current QC conclusion:

- Human adrenal: strongest overall signal, with TSS enrichment 19.19 and 20.39, pooled FRiP 0.565, and 150,470 conservative peaks
- Mouse adrenal: strong TSS enrichment 18.37 and 19.00, pooled FRiP 0.287, and 47,792 conservative peaks
- Human uterus: mixed replicate quality, pooled FRiP 0.184, and 79,438 conservative peaks
- Mouse uterus: weakest signal, TSS enrichment 5.70 and 6.24, pooled FRiP 0.153, and 42,412 conservative peaks

Conclusion: adrenal gland datasets show better and more consistent QC than uterus in both species, so adrenal gland was prioritized for downstream cross-species interpretation.

## Map open chromatin regions across species

Multi-species alignment:

`/ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal`

Current local inputs:

- Conservative peaks: `data/idr_Conservative_Peaks/`
- Optimal peaks: `data/idr_Optimal_Peaks/`
- HALPER outputs: `data/hal_Mapper_Peaks/`

Relevant scripts:

- `scripts/run_hal.sh`: maps conservative peaks between human and mouse
- `scripts/run_hal_optimal.sh`: maps optimal peaks between human and mouse

The mapping was performed with HALPER:

https://github.com/pfenninglab/halLiftover-postprocessing/tree/master

Example output files:

```text
data/hal_Mapper_Peaks/Human_AdrenalGland_idr.conservative_peak.HumanToMouse.HALPER.narrowPeak
data/hal_Mapper_Peaks/Mouse_AdrenalGland_idr.conservative_peak.MouseToHuman.HALPER.narrowPeak
data/hal_Mapper_Peaks/Human_Uterus_idr.conservative_peak.HumanToMouse.HALPER.narrowPeak
data/hal_Mapper_Peaks/Mouse_Uterus_idr.conservative_peak.MouseToHuman.HALPER.narrowPeak
```

## Compare open chromatin between species

Cross-species mapped OCR files are available locally in `data/hal_Mapper_Peaks/` for both adrenal gland and uterus. This section is still in progress, but the repository already contains the lifted peak sets needed for overlap, conservation, and category-level comparisons.

At the moment, the main completed downstream analyses built on these OCR sets are:

- functional enrichment with `rGREAT`
- promoter/enhancer partitioning for adrenal gland optimal peaks

## Find biological processes that are likely to be regulated by open chromatin regions

We used local `rGREAT` analysis on the IDR optimal peaks in `data/idr_Optimal_Peaks`, and generated GO enrichment tables and dot plots for all four datasets.

Relevant scripts:

- `scripts/rGREAT.R`: runs local `rGREAT` on each optimal peak set for `GO:BP`, `GO:CC`, and `GO:MF`
- `scripts/rGREAT_plot.R`: creates dot plots from each enrichment table

Results are stored in `results/rGREAT/`, with one directory per sample:

```text
results/rGREAT/Human_AdrenalGland_idr.optimal_peak/
results/rGREAT/Human_Uterus_idr.optimal_peak/
results/rGREAT/Mouse_AdrenalGland_idr.optimal_peak/
results/rGREAT/Mouse_Uterus_idr.optimal_peak/
```

Each sample directory contains:

- `rGREAT_BP.tsv`, `rGREAT_CC.tsv`, `rGREAT_MF.tsv`
- corresponding dot plot PNG files

Current interpretation:

These GO terms are standard Gene Ontology categories enriched from genes linked to OCR regions. Most of the enriched terms are broad functions such as gene expression, biosynthetic and metabolic processes, and general regulation. The results currently lack strong tissue-specific signals, suggesting that enrichment is dominated by housekeeping functions rather than adrenal- or uterus-specific pathways. The relatively modest fold enrichment also suggests limited specificity, so more refined peak-to-gene mapping or more targeted OCR subsets would likely improve biological interpretability.

## Compare candidate enhancers to candidate promoters

Adrenal gland optimal OCRs were partitioned into promoter-proximal and enhancer-like sets using a 2 kb window around annotated TSSs.

Relevant scripts:

- `scripts/adrenal_promoter_enhancer.sh`: creates promoter and enhancer peak sets with `bedtools`
- `scripts/run_hal_promoter_enhancer.sh`: maps adrenal promoter and enhancer sets across species with HALPER

Intermediate and output files are stored in `data/Promoters_and_Enhancers/`.

Current summary in `results/promoters_enhancers.txt`:

- Human adrenal: 206,765 total peaks, 72,666 promoter peaks (35.15%), 134,099 enhancer peaks (64.85%)
- Mouse adrenal: 48,263 total peaks, 27,105 promoter peaks (56.17%), 21,158 enhancer peaks (43.83%)

Mapped promoter/enhancer liftovers are already available for:

- mouse adrenal promoter to human
- mouse adrenal enhancer to human

## Find transcription factors that tend to bind open chromatin regions

todo
