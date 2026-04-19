# Comparative Analysis of Open Chromatin Regions in Human and Mouse Adrenal Gland and Uterus

## Introduction

Noncoding regions comprise the majority of the genome and harbor most regulatory elements that control when, where, and to what extent genes are expressed. Because these regulatory elements are often embedded within complex chromatin landscapes, identifying them directly remains challenging. Open chromatin regions (OCRs), which correspond to genomic loci accessible to transcription factors and other regulatory proteins, therefore provide a practical proxy for locating candidate cis-regulatory elements. The development of genome-wide chromatin accessibility assays, particularly ATAC-seq, has made it possible to profile OCRs rapidly and sensitively with low input material. Using such approaches, large-scale efforts such as ENCODE have generated comprehensive accessibility maps across diverse tissues and cell types, revealing both conserved and species-specific regulatory architectures across different species.

From an evolutionary perspective, cis-regulatory divergence has been proposed as a major driver of phenotypic differences between species. Comparative analyses of mammalian regulatory landscapes have shown that different classes of regulatory elements exhibit distinct evolutionary dynamics. Promoters tend to be relatively conserved, whereas enhancers often undergo more rapid turnover, even when overall gene expression programs remain broadly similar. These observations make cross-species OCR analysis a useful framework for studying regulatory conservation and divergence, while also highlighting a central technical challenge: regulatory orthology is often less straightforward to define than orthology for protein-coding genes.

Here, we designed a multi-step pipeline to characterize and compare OCRs between human and mouse. We first assessed the quality of ATAC-seq data from adrenal gland and uterus samples in both species to ensure the reliability of downstream analyses, and performed ontology enrichment analysis on OCRs. We then focused on adrenal OCRs, which exhibited higher data quality, and partitioned them into promoter-like and enhancer-like subsets. We performed cross-species comparisons using Mouse-to-Human HALPER-mapped adrenal promoter and enhancer peaks, classifying them based on their overlap with corresponding human adrenal OCRs. Finally, we used HOMER to identify enriched sequence motifs within these regulatory regions.

## Methods

### Data

TODO (check data resources)

### Quality Control Assessment

Because chromatin accessibility analyses depend strongly on signal quality and reproducible peak calling, the first stage of the project focused on QC. Key metrics were extracted from the QC reports with `scripts/1b_extract_qc_metrics.py`, combined into project-level tables with `scripts/1a_make_qc_table.py`, and supplemented by conservative peak counts generated with `scripts/1c_count_peaks.sh`. These summaries were used to compare TSS enrichment, pooled FRiP, and conservative peak counts across all four datasets. This step guides downstream decisions about which tissues were sufficiently reliable for more detailed analysis.

### Functional Enrichment Analysis

To assess whether OCRs were associated with interpretable biological functions, we ran local `rGREAT` enrichment on IDR optimal peak sets using `scripts/2a_rGREAT.R`, followed by visualization with `scripts/2b_rGREAT_plot.R`. Enrichment was evaluated across Gene Ontology biological process, cellular component, and molecular function categories. Optimal peak sets were used for this step to retain broader OCR coverage, with the understanding that this choice favors sensitivity over maximal reproducibility.

### Promoter and Enhancer Partitioning

For adrenal gland OCRs, `scripts/3a_call_adrenal_promoter_enhancer.sh` separated peaks into promoter-proximal and enhancer-like groups using a 2 kb window around annotated transcription start sites. This stage was limited to adrenal gland because the QC results indicated that adrenal was the stronger tissue in both species. The resulting files, together with later cross-species outputs, are stored in `data/Promoters_and_Enhancers/` and summarized in `results/Enhancer_and_Promoters/promoters_enhancers.txt`.

### Cross-Species Comparison of Adrenal OCRs

The cross-species analysis used the available Mouse-to-Human HALPER outputs together with `bedtools intersect`. HALPER was designed to reconstruct contiguous orthologous regulatory elements from `halLiftover` outputs and is particularly useful for comparative studies of noncoding regulatory DNA. In the current workflow, mouse adrenal promoter and enhancer peak sets had already been lifted into human coordinates, producing `mouse_adrenal_enhancer_peaks.MouseToHuman.HALPER.narrowPeak` and `mouse_adrenal_promoter_peaks.MouseToHuman.HALPER.narrowPeak`. These files were compared against `human_adrenal_enhancer_peaks.bed` and `human_adrenal_promoter_peaks.bed` using `scripts/4a_classify_conserved_peaks.sh`.

The script first assigned coordinate-based identifiers to the native mouse and human peak sets, then removed the summit suffix from the HALPER name field so that lifted intervals remained linked to the original mouse peaks. Lifted mouse peaks were intersected against the corresponding human adrenal peak set in human coordinates. Overlapping peaks were labeled as conserved. Human peaks overlapped by lifted mouse peaks were labeled as peaks with mouse ortholog support, whereas human peaks without overlap were retained as candidate human-specific peaks. Conserved and non-conserved mouse peak IDs were then mapped back to native mouse coordinates and written to `results/Enhancer_and_Promoters/Conserved/` and `results/Enhancer_and_Promoters/Specific/`. A run summary was also written to `results/Enhancer_and_Promoters/conserved_specific_summary.txt`.

### Motif Enrichment Analysis with HOMER

To examine whether conserved and species-specific adrenal OCR subsets were associated with distinct transcription factor binding signatures, we performed motif enrichment analysis using HOMER with `scripts/5_run_HOMER.sh`. This step was run separately for promoter and enhancer OCRs in both species. For each comparison, the species-specific OCR subset was used as the foreground set and the matched conserved subset was used as the background set. Specifically, the four comparisons were human enhancer candidate-specific versus human enhancer peaks with mouse ortholog support, human promoter candidate-specific versus human promoter peaks with mouse ortholog support, mouse enhancer-specific versus mouse enhancer-conserved, and mouse promoter-specific versus mouse promoter-conserved.

The script runs `findMotifsGenome.pl` in a SLURM environment and uses species-matched genome FASTA files as the genome reference argument (`hg38.fa` for human and `mm10.fa` for mouse). Motif discovery was performed with peak-width windows (`-size given`), repeat masking (`-mask`), explicit background peak sets (`-bg`), and parallel threads determined by `SLURM_CPUS_PER_TASK` with a default of 8. Output directories were written under `results/HOMER/` for each of the four specific-versus-conserved comparisons, and standard HOMER reports for known and de novo motifs were generated together with job log files.

### Analytical Constraints

Several decisions affect interpretation of the results. Conservative peak sets were favored when robustness was the main priority, whereas optimal peaks were used for enrichment analysis to retain broader OCR representation. Adrenal gland was prioritized over uterus for promoter/enhancer and cross-species analysis because of stronger QC. The conservation analysis was also one-directional. Because the repository contains Mouse-to-Human HALPER outputs but not reciprocal Human-to-Mouse results, peaks classified as candidate human-specific should be interpreted cautiously. Lack of overlap in this framework may reflect true biological divergence, but it may also reflect incomplete alignability, liftover failure, or the absence of reciprocal mapping.

## Results

### Quality Control Patterns Across Datasets

The QC analysis showed a clear difference between adrenal gland and uterus datasets. Human adrenal had the strongest combined signal, with high TSS enrichment, high pooled FRiP, and a large number of conservative peaks. Mouse adrenal also showed relatively strong TSS enrichment and acceptable overall quality. In contrast, the uterus datasets, especially mouse uterus, showed weaker signal and lower confidence for downstream peak-based interpretation. These results justified prioritizing adrenal gland for promoter/enhancer partitioning and cross-species comparison.

### Functional Enrichment of OCR Sets

The `rGREAT` analyses linked OCRs to Gene Ontology categories, but the enriched terms were dominated by broad processes such as gene expression, biosynthesis, metabolism, and general regulation rather than sharply tissue-specific functions. This result suggests that the OCR sets captured authentic regulatory structure, but that whole-set enrichment using broad peak collections was not sufficient to produce a highly specific adrenal- or uterus-centered functional narrative. This outcome is broadly consistent with large-scale accessibility resources showing that whole-tissue OCR collections often contain mixtures of constitutive and cell-type-restricted regulatory elements.

### Promoter and Enhancer Composition of Adrenal OCRs

Promoter/enhancer partitioning revealed a difference in OCR composition between species. Human adrenal contained 206,765 total peaks, including 72,666 promoter peaks and 134,099 enhancer peaks. Mouse adrenal contained 48,263 total peaks, including 27,105 promoter peaks and 21,158 enhancer peaks. Thus, the human adrenal OCR set was relatively enhancer-enriched, whereas the mouse adrenal OCR set contained a comparatively larger promoter-proximal fraction. This pattern may reflect biological differences, but it may also be influenced by technical factors such as sequencing depth, annotation quality, or differences in peak calling sensitivity.

### Cross-Species Conservation of Adrenal Enhancers and Promoters

The Mouse-to-Human comparison recovered only a minority of mouse adrenal peaks as conserved under the current overlap-based definition. For enhancer-like OCRs, the mouse adrenal set contained 21,158 total peaks, of which 6,442 were lifted to the human genome in the available HALPER output. Among all mouse adrenal enhancer peaks, 2,516 were classified as conserved and 18,642 as mouse-specific after mapping back to native mouse coordinates. On the human side, 4,360 adrenal enhancer peaks had overlap from lifted mouse peaks, whereas 129,739 were classified as candidate human-specific.

For promoter-like OCRs, the mouse adrenal set contained 27,105 peaks, of which 3,533 lifted to the human genome. Among these, 2,727 were classified as conserved and 24,378 as mouse-specific in native mouse coordinates. On the human side, 8,323 promoter peaks had mouse ortholog support and 64,343 were classified as candidate human-specific. These values were taken directly from `results/Enhancer_and_Promoters/conserved_specific_summary.txt`.

Taken together, the results indicate limited overlap under the current one-directional conservation framework, with especially high fractions of candidate species-specific peaks on the human side. The enhancer results are broadly consistent with prior comparative studies suggesting that enhancer activity is more evolutionarily labile than promoter activity. However, the relatively modest overlap even for promoters suggests that technical and analytical factors contribute alongside true biological divergence.

### Motif Enrichment of Adrenal Enhancers and Promoters

TODO (waiting for results)

## Discussion

TODO

## References

(Need more for background and methods)

Carelli FN, Liechti A, Halbert J, Warnefors M, Kaessmann H. Repurposing of promoters and enhancers during mammalian evolution. *Nature Communications*. 2018;9:4066. https://doi.org/10.1038/s41467-018-06544-z

ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. *Nature*. 2012;489:57-74. https://doi.org/10.1038/nature11247

Hickey G, Paten B, Earl D, et al. HAL: a hierarchical format for storing and analyzing multiple genome alignments. *Bioinformatics*. 2013;29:1341-1342. https://doi.org/10.1093/bioinformatics/btt128

Gu Z, Hubschmann D. rGREAT: an R/Bioconductor package for functional enrichment on genomic regions. *Bioinformatics*. 2023;39(1):btac745. https://doi.org/10.1093/bioinformatics/btac745

Zhang X, Kaplow IM, Wirthlin M, Park TY, Pfenning AR. HALPER facilitates the identification of regulatory element orthologs across species. *Bioinformatics*. 2020;36:4339-4340. https://doi.org/10.1093/bioinformatics/btaa493
