# Adrenal_gland-Uterus-OCR-Analysis

03-713: Bioinformatics Data Integration Practicum, Spring 2026

Team members: Xinyi Li, Jimmy Lee, Yu Hsin Chen

## Project Overview

This repository analyzes open chromatin regions (OCRs) in adrenal gland ATAC-seq datasets from human and mouse, with uterus QC retained as part of the original dataset screening. The main downstream analysis focuses on adrenal gland OCRs, including functional enrichment, cross-species conservation, promoter/enhancer partitioning, and motif enrichment.

The repository currently includes:

- QC summaries for human adrenal gland, mouse adrenal gland, human uterus, and mouse uterus datasets
- Human and mouse adrenal IDR optimal peak sets
- Mouse-to-human HALPER-mapped adrenal OCRs
- Whole-adrenal conserved and species-specific OCR classification
- rGREAT enrichment results and dot plots for adrenal OCRs
- Promoter-proximal and enhancer-like adrenal OCR partitioning
- Promoter/enhancer conserved and species-specific OCR classification
- HOMER motif enrichment results for specific vs. conserved promoter/enhancer OCR subsets

Pipeline design:

![Pipeline_Design](Pipeline_Design.png)

## Repository Structure

```text
.
├── data/
│   ├── qc_html/                    input QC HTML reports
│   ├── idr_Optimal_Peaks/          adrenal optimal peaks and Mouse->Human HALPER file
│   └── Promoters_and_Enhancers/    adrenal promoter/enhancer peak sets and HALPER outputs
├── results/
│   ├── qc/                         checked QC HTML reports and summary tables
│   ├── Specific_and_Conserved/     whole-adrenal conserved/specific OCR calls
│   ├── rGREAT/                     rGREAT TSV outputs and dot plots
│   ├── Enhancer_and_Promoters/     promoter/enhancer conserved/specific OCR calls
│   └── HOMER/                      motif enrichment results
├── scripts/
│   ├── 1a_qc_html_report.py
│   ├── 1b_make_qc_table.py
│   ├── 2a_classify_adrenal_specific_conserved.sh
│   ├── 2b_rGREAT.R
│   ├── 2c_rGREAT_plot.R
│   ├── 3a_call_adrenal_promoter_enhancer.sh
│   ├── 3b_run_hal_promoter_enhancer.sh
│   ├── 4a_classify_conserved_peaks.sh
│   └── 5a_run_HOMER.sh
└── README.md
```

## Dependencies

Core dependencies:

- Python 3.11 or newer for QC parsing/report scripts
- R with `rGREAT`, `GenomicRanges`, `IRanges`, and `ggplot2`
- `bedtools` for interval operations
- HAL / `halLiftover` and HALPER for cross-species OCR mapping
- HOMER for motif enrichment
- Human and mouse genome FASTA files, such as `hg38.fa` and `mm10.fa`
- Human and mouse genome index files, such as `hg38.fa.fai` and `mm10.fa.fai`
- Human and mouse TSS BED annotation files
- A multi-species `.hal` alignment file

The HALPER and HOMER steps are written for a SLURM-based Linux cluster. Other steps use repository-relative paths and can be run independently when the required files are available.

## Quick Start

The instructions below assume you are running in a Linux environment.

### 1. Clone the repository

```bash
git clone https://github.com/BioinformaticsDataPracticum2026/Adrenal_gland-Uterus-OCR-Analysis.git
cd Adrenal_gland-Uterus-OCR-Analysis
```

### 2. Create and activate a conda environment

The Python scripts in this repository currently use only the Python standard library, so there is no required `requirements.txt` at the moment. We recommend using a conda environment to keep the Python, R, and command-line dependencies together.

```bash
conda create -n adrenal-ocr-analysis python=3.11 -y
conda activate adrenal-ocr-analysis
conda install -c conda-forge -c bioconda pip r-base bedtools -y                        
```

### 3. Install analysis dependencies

Install the required R packages:

```bash
Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'
Rscript -e 'BiocManager::install("rGREAT", ask=FALSE, update=FALSE, INSTALL_opts="--no-multiarch")'
```

Install HALPER and the post-processing helper script:

```bash
git clone https://github.com/pfenninglab/halLiftover-postprocessing.git
```

Install HOMER:

```bash
mkdir homer
cd homer
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install
cd ..
```

You can also refer to the official project instructions:

- HAL / halLiftover: https://github.com/ComparativeGenomicsToolkit/hal
- HALPER post-processing: https://github.com/pfenninglab/halLiftover-postprocessing
- HOMER: http://homer.ucsd.edu/homer/
- bedtools: https://bedtools.readthedocs.io/
- Bioconductor `rGREAT`: https://bioconductor.org/packages/rGREAT/
- Bioconductor `GenomicRanges`: https://bioconductor.org/packages/GenomicRanges/
- Bioconductor `IRanges`: https://bioconductor.org/packages/IRanges/
- `ggplot2`: https://ggplot2.tidyverse.org/

### 4. Configure local paths for external resources

Some scripts use repository-relative input/output paths,  but still require external resources that you must provide in your  environment. The current repository already includes part of the initial input data and several intermediate results, so some scripts can be run directly without reproducing every upstream step from scratch.

Update these variables before running the corresponding scripts:

- ```bash
  scripts/3a_call_adrenal_promoter_enhancer.sh
  ```

  - `HUMAN_TSS`
  - `HUMAN_GENOME_INDEX`
  - `MOUSE_TSS`
  - `MOUSE_GENOME_INDEX`

- ```bash
  scripts/3b_run_hal_promoter_enhancer.sh
  ```

  - `HALPER_SCRIPT`
  - `HAL`

- ```bash
  scripts/5a_run_HOMER.sh
  ```

  - `HUMAN_FA`
  - `MOUSE_FA`
  - the `export PATH=.../HOMER/bin:$PATH` line

### 5. Run the pipeline

From the repository root:

```bash
# QC reports and summary table
python scripts/1a_qc_html_report.py
python scripts/1b_make_qc_table.py --print-summary

# Whole-adrenal conserved/specific classification
bash scripts/2a_classify_adrenal_specific_conserved.sh

# rGREAT enrichment and plots
Rscript scripts/2b_rGREAT.R
Rscript scripts/2c_rGREAT_plot.R

# Adrenal promoter/enhancer partitioning
bash scripts/3a_call_adrenal_promoter_enhancer.sh

# Submit HALPER promoter/enhancer liftover jobs
bash scripts/3b_run_hal_promoter_enhancer.sh

# Promoter/enhancer conserved/specific classification
bash scripts/4a_classify_conserved_peaks.sh

# HOMER motif enrichment
sbatch scripts/5a_run_HOMER.sh
```

### 6. Check outputs

Main outputs are written to:

- `results/qc/`
- `results/rGREAT/`
- `results/Enhancer_and_Promoters/`
- `results/HOMER/`

## 1. Evaluate data quality

QC summaries are stored in `results/qc/`:

- `*_checked.html`: threshold-annotated QC HTML reports generated from the raw ENCODE-style reports
- `qc_summary_table.tsv` and `qc_summary_table.md`: combined metrics table across all four datasets

Relevant scripts:

- `scripts/1a_qc_html_report.py`: enhances the ENCODE-style QC HTML reports with threshold annotations
- `scripts/1b_make_qc_table.py`: parses the QC HTML  reports, builds the final QC summary table, and can print a concise  terminal summary including TSS enrichment, pooled FRiP, NRF/PBC metrics, and IDR peak counts

To run:

```bash
# Enhance QC HTML reports and write annotated HTML files to results/qc
python scripts/1a_qc_html_report.py

# Build combined QC summary table from the HTML reports in data/qc_html
# and print a concise QC summary to the terminal
python scripts/1b_make_qc_table.py --print-summary

# Run both QC scripts in one command
# PowerShell:
python scripts/1a_qc_html_report.py; python scripts/1b_make_qc_table.py --print-summary

# bash:
python scripts/1a_qc_html_report.py && python scripts/1b_make_qc_table.py --print-summary 
```

`1a_qc_html_report.py` reads QC HTML files from `data/qc_html/` by default and writes annotated `*_checked.html` reports to `results/qc/`. `1b_make_qc_table.py` reads the same HTML inputs by default, writes `qc_summary_table.tsv` and `qc_summary_table.md` to `results/qc/`, and optionally prints a concise per-sample summary with:

- `TSS Rep1` / `TSS Rep2`
- `FRiP (pooled)`
- `IDR Optimal Peaks` / `IDR Conservative Peaks`
- `NRF`, `PBC1`, and `PBC2` for each replicate

## 2. Conserved and Species-Specific OCRs

`scripts/2a_classify_adrenal_specific_conserved.sh` classifies adrenal OCRs using one-directional Mouse->Human HALPER overlap against native human adrenal OCRs.

Definitions:

- Mouse conserved: source mouse peak has a lifted Mouse->Human interval overlapping a human adrenal peak
- Mouse specific: source mouse peak has no lifted interval overlapping a human adrenal peak
- Human with mouse ortholog: human peak overlaps a lifted mouse peak
- Candidate human specific: human peak has no Mouse->Human overlap

Run:

```bash
bash scripts/2a_classify_adrenal_specific_conserved.sh
```

Current outputs are in `results/Specific_and_Conserved/`:

- `mouse_coord_mm10.Mouse_AdrenalGland.conserved_by_MouseToHuman_overlap_human_adrenal.bed`
- `mouse_coord_mm10.Mouse_AdrenalGland.specific_no_MouseToHuman_overlap_human_adrenal.bed`
- `human_coord_hg38.Human_AdrenalGland.conserved_overlap_mouse_MouseToHuman.bed`
- `human_coord_hg38.Human_AdrenalGland.specific_no_overlap_mouse_MouseToHuman.bed`
- `adrenal_specific_conserved_summary.txt`

## 3. rGREAT Functional Enrichment

`scripts/2b_rGREAT.R` runs local rGREAT enrichment for GO biological process, cellular component, and molecular function terms. It is configured for:

- Human adrenal optimal OCRs
- Mouse adrenal optimal OCRs
- Human conserved adrenal OCRs
- Human-specific adrenal OCRs
- Mouse-specific adrenal OCRs

`scripts/2c_rGREAT_plot.R` reads rGREAT TSV files and generates dot plots for the top GO terms.

Run:

```bash
Rscript scripts/2b_rGREAT.R
Rscript scripts/2c_rGREAT_plot.R
```

Current checked-in rGREAT result directories:

```bash
results/rGREAT/Human_AdrenalGland_idr.optimal_peak/
results/rGREAT/Mouse_AdrenalGland_idr.optimal_peak/
```

Each directory contains:

- `rGREAT_BP.tsv`
- `rGREAT_CC.tsv`
- `rGREAT_MF.tsv`
- corresponding `*_dotplot.png` files

## 4. Compare candidate enhancers to candidate promoters

Adrenal gland optimal OCRs were partitioned into  promoter-proximal and enhancer-like sets using a 2 kb window around  annotated TSSs.

Relevant scripts:

- `scripts/3a_call_adrenal_promoter_enhancer.sh`: creates promoter and enhancer peak sets with `bedtools`
- `scripts/3b_run_hal_promoter_enhancer.sh`: maps adrenal promoter and enhancer sets across species with HALPER

To run:

```bash
# Partition optimal peaks into promoter-proximal and enhancer-like sets
bash scripts/3a_call_adrenal_promoter_enhancer.sh

# Submit HALPER liftover job to the cluster
sbatch scripts/3b_run_hal_promoter_enhancer.sh           
```

`3a_call_adrenal_promoter_enhancer.sh` now reads adrenal optimal peak files from `data/idr_Optimal_Peaks/` using repo-relative paths and writes promoter/enhancer BED files to `data/Promoters_and_Enhancers/`.

Before running `3a_call_adrenal_promoter_enhancer.sh`, you still need to provide the external annotation and genome index files by setting:

- `HUMAN_TSS`
- `HUMAN_GENOME_INDEX`
- `MOUSE_TSS`
- `MOUSE_GENOME_INDEX`

`3b_run_hal_promoter_enhancer.sh` uses the repository-local peak files in `data/Promoters_and_Enhancers/`, but `HALPER` and the `.hal` alignment file are still external dependencies. Before running it, update:

- `HALPER_SCRIPT` to your local `halper_map_peak_orthologs.sh`
- `HAL` to your local multi-species `.hal` alignment file

Intermediate and output files are stored in `data/Promoters_and_Enhancers/`.

Current summary from `data/Promoters_and_Enhancers/`:

- Human adrenal: 206,765 total peaks, 72,666 promoter peaks (35.15%), 134,099 enhancer peaks (64.85%)
- Mouse adrenal: 48,263 total peaks, 27,105 promoter peaks (56.17%), 21,158 enhancer peaks (43.83%)

Mapped promoter/enhancer liftovers are already available for:

- mouse adrenal promoter to human
- mouse adrenal enhancer to human

## 5. Compare open chromatin between species

Adrenal promoter and enhancer peak sets were compared  between mouse and human using the available Mouse→Human HALPER outputs  together with `bedtools intersect`. Conserved peaks were  defined as lifted mouse peaks that overlap the corresponding human peak  set, and the remaining peaks were retained as species-specific  candidates.

Relevant script:

- `scripts/4a_classify_conserved_peaks.sh`: compares the adrenal promoter/enhancer peak sets in `data/Promoters_and_Enhancers/` and writes conserved and specific peak files to `results/Enhancer_and_Promoters/`

To run (requires `bedtools` ≥ 2.30):

```bash
bash scripts/4a_classify_conserved_peaks.sh
```

Results are stored in `results/Enhancer_and_Promoters/`:

```bash
results/Enhancer_and_Promoters/Conserved/
results/Enhancer_and_Promoters/Specific/
results/Enhancer_and_Promoters/conserved_specific_summary.txt
```

## 6. Find transcription factors that tend to bind open chromatin regions using HOMER

We used `HOMER` motif enrichment to compare  species-specific adrenal OCR subsets against conserved OCR subsets as  background. This analysis was run separately for promoter and enhancer  peaks in human and mouse.

Foreground and background sets:

- human adrenal enhancer specific vs. human adrenal enhancer with mouse ortholog
- human adrenal promoter specific vs. human adrenal promoter with mouse ortholog
- mouse adrenal enhancer specific vs. mouse adrenal enhancer conserved
- mouse adrenal promoter specific vs. mouse adrenal promoter conserved

The script is configured to use repository-relative BED  inputs together with local genome FASTA files as the HOMER genome  argument:

- `data/Promoters_and_Enhancers/Specific/`
- `data/Promoters_and_Enhancers/Conserved/`
- local FASTA files such as `hg38.fa` and `mm10.fa`
- `-size given`
- `-mask`
- the matched conserved set as `-bg`
- `SLURM_CPUS_PER_TASK` (default `8`) for parallel threads

To run on the cluster:

```bash
sbatch scripts/5a_run_HOMER.sh
```

**Note:** Before running in a different environment, update the following variables in `scripts/5a_run_HOMER.sh`:

- `HUMAN_FA` / `MOUSE_FA` — paths to the `hg38.fa` and `mm10.fa` genome FASTA files
- the `export PATH` line — path to the HOMER `bin/` directory on your system

`5a_run_HOMER.sh` expects you to install HOMER yourself. It uses local genome FASTA file paths as the second argument to `findMotifsGenome.pl`, rather than HOMER's built-in short genome names.

Results are written to:

```bash
results/HOMER/human_enhancer_specific_vs_conserved/
results/HOMER/human_promoter_specific_vs_conserved/
results/HOMER/mouse_enhancer_specific_vs_conserved/
results/HOMER/mouse_promoter_specific_vs_conserved/  
```

Each output directory contains the standard HOMER motif  enrichment reports for known and de novo motifs, along with the  corresponding log files in `logs/`.

## Authors

Xinyi Li [xinyili3@andrew.cmu.edu](mailto:xinyili3@andrew.cmu.edu), Jimmy Lee [jimmylee@andrew.cmu.edu](mailto:jimmylee@andrew.cmu.edu), Yu Hsin Chen [heather4@andrew.cmu.edu](mailto:heather4@andrew.cmu.edu)

## To Cite this Repository

Xinyi Li, Jimmy Lee, Yu Hsin Chen (2026).  Adrenal_gland-Uterus-OCR-Analysis. 03-713: Bioinformatics  Data  Integration Practicum, Carnegie Mellon Univeristy.

## References

1. Gu Z, Hubschmann D. rGREAT: an R/Bioconductor package for functional enrichment on genomic regions. *Bioinformatics*. 2023;39(1):btac745. https://doi.org/10.1093/bioinformatics/btac745
2. Wickham H. *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York; 2016. https://ggplot2.tidyverse.org
3. Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010;26(6):841-842. https://doi.org/10.1093/bioinformatics/btq033
4. Zhang X, Kaplow I, Wirthlin M, Park T, Pfenning A. HALPER  facilitates the identification of regulatory element orthologs across  species. *Bioinformatics*. 2020;36(15):4339-4340. https://doi.org/10.1093/bioinformatics/btaa378
5. Hickey G, Paten B, Earl D, Zerbino D, Haussler D. HAL: a  hierarchical format for storing and analyzing multiple genome  alignments. *Bioinformatics*. 2013;29(10):1341-1342. https://doi.org/10.1093/bioinformatics/btt128
6. Paten B, Earl D, Nguyen N, Diekhans M, Zerbino D, Haussler D. Cactus: algorithms for genome multiple sequence alignment. *Genome Research*. 2011;21(9):1512-1528. https://doi.org/10.1101/gr.123356.111
7. Heinz S, Benner C, Spann N, et al. Simple combinations of  lineage-determining transcription factors prime cis-regulatory elements  required for macrophage and B cell identities. *Molecular Cell*. 2010;38(4):576-589. https://doi.org/10.1016/j.molcel.2010.05.004
8. Liu C, Wang M, Wei X, et al. An ATAC-seq atlas of chromatin accessibility in mouse tissues. *Scientific Data*. 2019;6:65. https://doi.org/10.1038/s41597-019-0071-0

Useful package pages for the Bioconductor dependencies used in `scripts/rGREAT.R`:

- GenomicRanges: https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
- IRanges: https://bioconductor.org/packages/release/bioc/html/IRanges.html
