library(rGREAT)
library(GenomicRanges)
library(IRanges)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(script_arg) > 0) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("scripts/2c_rGREAT.R", winslash = "/", mustWork = FALSE)
}
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

peak_inputs <- data.frame(
  peak_file = c(
    file.path(repo_root, "data", "idr_Optimal_Peaks", "Human_AdrenalGland_idr.optimal_peak.narrowPeak"),
    file.path(repo_root, "data", "idr_Optimal_Peaks", "Mouse_AdrenalGland_idr.optimal_peak.narrowPeak"),
    file.path(repo_root, "results", "Specific_and_Conserved", "human_coord_hg38.Human_AdrenalGland.conserved_overlap_mouse_MouseToHuman.bed"),
    file.path(repo_root, "results", "Specific_and_Conserved", "human_coord_hg38.Human_AdrenalGland.specific_no_overlap_mouse_MouseToHuman.bed"),
    file.path(repo_root, "results", "Specific_and_Conserved", "mouse_coord_mm10.Mouse_AdrenalGland.specific_no_MouseToHuman_overlap_human_adrenal.bed")
  ),
  genome = c("hg38", "mm10", "hg38", "hg38", "mm10"),
  stringsAsFactors = FALSE
)

output_dir <- file.path(repo_root, "results", "rGREAT")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_narrowpeak_as_granges <- function(path) {
  peak_df <- utils::read.delim(
    path,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  # Convert to 1-based start for GRanges
  GenomicRanges::GRanges(
    seqnames = peak_df[[1]],
    ranges = IRanges::IRanges(start = peak_df[[2]] + 1L, end = peak_df[[3]])
  )
}

run_single_ontology <- function(gr, ontology, genome) {
  rGREAT::great(
    gr = gr,
    gene_sets = ontology,
    tss_source = paste0("GREAT:", genome)
  )
}

# Ontology: BP, CC, MF
ontology_plan <- list(
  BP = list(label = "GO:BP", resolver = function(gr, genome) {
    list(object = run_single_ontology(gr, "GO:BP", genome), ontology_used = "GO:BP")
  }),
  CC = list(label = "GO:CC", resolver = function(gr, genome) {
    list(object = run_single_ontology(gr, "GO:CC", genome), ontology_used = "GO:CC")
  }),
  MF = list(label = "GO:MF", resolver = function(gr, genome) {
    list(object = run_single_ontology(gr, "GO:MF", genome), ontology_used = "GO:MF")
  })
)

for (i in seq_len(nrow(peak_inputs))) {
  peak_file <- peak_inputs$peak_file[[i]]
  genome <- peak_inputs$genome[[i]]
  file_name <- basename(peak_file)
  sample_name <- sub("\\.(bed|narrowPeak)$", "", file_name)
  sample_output_dir <- file.path(output_dir, sample_name)

  dir.create(sample_output_dir, recursive = TRUE, showWarnings = FALSE)

  message("\n[Sample] ", sample_name)
  message("  File: ", peak_file)
  message("  Genome: ", genome)

  gr <- read_narrowpeak_as_granges(peak_file)
  message("  Regions loaded: ", length(gr))

  for (ontology_name in names(ontology_plan)) {
    message("  Running ", ontology_name, " enrichment")
    resolved <- ontology_plan[[ontology_name]]$resolver(gr, genome)
    enrich_tbl <- rGREAT::getEnrichmentTable(resolved$object)

    enrich_tbl$sample <- sample_name
    enrich_tbl$genome <- genome
    enrich_tbl$ontology <- ontology_name
    enrich_tbl$ontology_input <- resolved$ontology_used
    
    # Write results
    out_tsv <- file.path(sample_output_dir, paste0("rGREAT_", ontology_name, ".tsv"))
    utils::write.table(
      enrich_tbl,
      file = out_tsv,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}

message("\nDone.")
message("Results written to: ", output_dir)
