# Adrenal_gland-Uterus-OCR-Analysis

03-713: Bioinformatics Data Integration Practicum, Spring 2026

Team members: Xinyi L., Jimmy L., Yu Hsin C.

## Project Overview

This project compares open chromatin regions (OCRs) between adrenal gland and uterus in human and mouse ATAC-seq datasets. The current repository includes:

- QC summaries for all four datasets
- IDR conservative and optimal peak sets
- Cross-species OCR mapping with HAL/HALPER
- Species-conserved and species-specific peak classification (Mouse→Human)
- Local `rGREAT` functional enrichment analysis on optimal OCRs
- Promoter/enhancer partitioning for adrenal gland OCRs

Repository structure:

```text
data/
  qc_html/                    ENCODE-style QC HTML reports
  idr_Conservative_Peaks/     conservative peak sets used for cross-species mapping
  idr_Optimal_Peaks/          optimal peak sets used for downstream analyses
  halper_peaks/               Mouse→Human HALPER liftover outputs
  Promoters_and_Enhancers/    promoter/enhancer split peak sets and lifted results
scripts/
  1a-1c                       QC summary and peak counting scripts
  2a-2b                       rGREAT enrichment and plotting scripts
  3a-3b                       adrenal promoter/enhancer calling and HALPER mapping
  4a                          species-conserved vs. species-specific peak classification
results/
  qc/                         QC tables and interpretation
  rGREAT/                     GO enrichment tables and dot plots
  Enhancer_and_Promoters/     adrenal promoter/enhancer summary
  peak_classification/        conserved and species-specific peak calls
```

## Dependencies

The analyses in this repository depend on the following software and packages:

- `Python 3` for QC summary scripts
- `R` for functional enrichment and plotting
- R packages: `rGREAT`, `GenomicRanges`, `IRanges`, `ggplot2`
- `bedtools` for promoter/enhancer partitioning and peak classification
- `HALPER` and `halLiftover` for cross-species OCR mapping
- access to the multi-species HAL alignment file `10plusway-master.hal`
- a SLURM environment for submitting the HAL/HALPER jobs in `scripts/3b_run_hal_promoter_enhancer.sh`

## 1. Evaluate data quality

QC summaries are stored in `results/qc/`:

- `qc_summary_table.tsv` and `qc_summary_table.md`: combined metrics table
- `qc_summary.txt`: extracted TSS enrichment and pooled FRiP values
- `qc_interpretation.md`: short project-level interpretation

Relevant scripts:

- `scripts/1a_make_qc_table.py`: builds the final QC summary table
- `scripts/1b_extract_qc_metrics.py`: extracts key QC metrics from ENCODE ATAC-seq outputs
- `scripts/1c_count_peaks.sh`: counts reproducible IDR conservative peaks

Summary of the current QC conclusion:

- Human adrenal: strongest overall signal, with TSS enrichment 19.19 and 20.39, pooled FRiP 0.565, and 150,470 conservative peaks
- Mouse adrenal: strong TSS enrichment 18.37 and 19.00, pooled FRiP 0.287, and 47,792 conservative peaks
- Human uterus: mixed replicate quality, pooled FRiP 0.184, and 79,438 conservative peaks
- Mouse uterus: weakest signal, TSS enrichment 5.70 and 6.24, pooled FRiP 0.153, and 42,412 conservative peaks

Conclusion: adrenal gland datasets show better and more consistent QC than uterus in both species, so adrenal gland was prioritized for downstream cross-species interpretation.

## 2. Find biological processes that are likely to be regulated by open chromatin regions

We used local `rGREAT` analysis on the IDR optimal peaks in `data/idr_Optimal_Peaks`, and generated GO enrichment tables and dot plots for all four datasets.

Relevant scripts:

- `scripts/2a_rGREAT.R`: runs local `rGREAT` on each optimal peak set for `GO:BP`, `GO:CC`, and `GO:MF`
- `scripts/2b_rGREAT_plot.R`: creates dot plots from each enrichment table

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

## 3. Compare candidate enhancers to candidate promoters

Adrenal gland optimal OCRs were partitioned into promoter-proximal and enhancer-like sets using a 2 kb window around annotated TSSs.

Relevant scripts:

- `scripts/3a_call_adrenal_promoter_enhancer.sh`: creates promoter and enhancer peak sets with `bedtools`
- `scripts/3b_run_hal_promoter_enhancer.sh`: maps adrenal promoter and enhancer sets across species with HALPER

Intermediate and output files are stored in `data/Promoters_and_Enhancers/`.

Current summary in `results/Enhancer_and_Promoters/promoters_enhancers.txt`:

- Human adrenal: 206,765 total peaks, 72,666 promoter peaks (35.15%), 134,099 enhancer peaks (64.85%)
- Mouse adrenal: 48,263 total peaks, 27,105 promoter peaks (56.17%), 21,158 enhancer peaks (43.83%)

Mapped promoter/enhancer liftovers are already available for:

- mouse adrenal promoter to human
- mouse adrenal enhancer to human

## 4. Compare open chromatin between species

Mouse IDR conservative peaks were mapped to the human genome using [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing) (Mouse→Human direction) via the Cactus whole-genome alignment `10plusway-master.hal`. Liftover outputs are stored in `data/halper_peaks/`.

Peaks were then classified using one-directional Mouse→Human overlap:

| Class | Definition | Confidence |
|---|---|---|
| Mouse-conserved | Mouse peak lifts to human genome AND overlaps a native human peak | High |
| Mouse-specific | Mouse peak does not lift over, OR lifts over but misses all human peaks | High |
| Human peak with mouse ortholog | Human peak overlapped by an M→H HALPER entry | Moderate |
| Candidate human-specific | Human peak with no M→H overlap | Moderate |

> Human classifications are inferred from the mouse perspective. Without Human→Mouse liftover data, human-specific calls may also reflect alignment gaps rather than true regulatory divergence — this is acknowledged as a methodological limitation.

Relevant script:

- `scripts/4a_classify_conserved_peaks.sh`: runs `bedtools intersect` on HALPER outputs and IDR conservative peaks; writes four output files per tissue to `results/peak_classification/`

To run (requires `bedtools` ≥ 2.30):

```bash
bash scripts/4a_classify_conserved_peaks.sh
```

Output files in `results/peak_classification/`:

```text
{Tissue}_mouse_conserved.narrowPeak
{Tissue}_mouse_specific.narrowPeak
{Tissue}_human_with_mouse_ortholog.narrowPeak
{Tissue}_human_candidate_specific.narrowPeak
```

## 5. Find transcription factors that tend to bind open chromatin regions

todo

## Citations

1. Gu Z, Hubschmann D. rGREAT: an R/Bioconductor package for functional enrichment on genomic regions. *Bioinformatics*. 2023;39(1):btac745. https://doi.org/10.1093/bioinformatics/btac745
2. Wickham H. *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York; 2016. https://ggplot2.tidyverse.org
3. Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010;26(6):841-842. https://doi.org/10.1093/bioinformatics/btq033
4. Zhang X, Kaplow I, Wirthlin M, Park T, Pfenning A. HALPER facilitates the identification of regulatory element orthologs across species. *Bioinformatics*. 2020;36(15):4339-4340. https://doi.org/10.1093/bioinformatics/btaa378
5. Hickey G, Paten B, Earl D, Zerbino D, Haussler D. HAL: a hierarchical format for storing and analyzing multiple genome alignments. *Bioinformatics*. 2013;29(10):1341-1342. https://doi.org/10.1093/bioinformatics/btt128
6. Paten B, Earl D, Nguyen N, Diekhans M, Zerbino D, Haussler D. Cactus: algorithms for genome multiple sequence alignment. *Genome Research*. 2011;21(9):1512-1528. https://doi.org/10.1101/gr.123356.111

Useful package pages for the Bioconductor dependencies used in `scripts/rGREAT.R`:

- GenomicRanges: https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
- IRanges: https://bioconductor.org/packages/release/bioc/html/IRanges.html
