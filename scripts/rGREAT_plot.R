library(ggplot2)

input_tasks <- data.frame(
  input_dir = c(
    "results/rGREAT/Human_AdrenalGland_idr.optimal_peak",
    "results/rGREAT/Human_Uterus_idr.optimal_peak",
    "results/rGREAT/Mouse_AdrenalGland_idr.optimal_peak",
    "results/rGREAT/Mouse_Uterus_idr.optimal_peak"
  ),
  task_name = c(
    "Human_AdrenalGland",
    "Human_Uterus",
    "Mouse_AdrenalGland",
    "Mouse_Uterus"
  ),
  stringsAsFactors = FALSE
)

top_n <- 20

tsv_files <- unlist(
  lapply(
    input_tasks$input_dir,
    function(input_dir) {
      list.files(path = input_dir, full.names = TRUE, recursive = TRUE)
    }
  ),
  use.names = FALSE
)

make_rgreat_plot <- function(tsv_file, task_name, top_n = 20) {
  enrich_tbl <- utils::read.delim(
    tsv_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c("description", "fold_enrichment", "observed_region_hits", "p_adjust")
  missing_cols <- setdiff(required_cols, colnames(enrich_tbl))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s) in ", tsv_file, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  enrich_tbl <- enrich_tbl[is.finite(enrich_tbl$p_adjust) & enrich_tbl$p_adjust > 0, , drop = FALSE]

  if (nrow(enrich_tbl) == 0) {
    message("Skipping empty/invalid file: ", tsv_file)
    return(invisible(NULL))
  }

  plot_tbl <- enrich_tbl[order(enrich_tbl$p_adjust, decreasing = FALSE), , drop = FALSE]
  plot_tbl <- utils::head(plot_tbl, top_n)
  plot_tbl$log10_p_adjust <- -log10(plot_tbl$p_adjust)
  plot_tbl$description <- factor(
    plot_tbl$description,
    levels = rev(plot_tbl$description[order(plot_tbl$fold_enrichment, decreasing = FALSE)])
  )

  ontology_name <- sub("^rGREAT_(.+)\\.tsv$", "\\1", basename(tsv_file))

  p <- ggplot(
    plot_tbl,
    aes(
      x = fold_enrichment,
      y = description,
      size = observed_region_hits,
      color = log10_p_adjust
    )
  ) +
    geom_point(alpha = 0.8) +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = paste0(task_name, " - ", ontology_name),
      x = "Fold enrichment",
      y = "GO term",
      size = "Region hits",
      color = "-log10(p_adjust)"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )

  out_file <- file.path(
    dirname(tsv_file),
    paste0(task_name, "_", ontology_name, "_dotplot.png")
  )

  ggplot2::ggsave(out_file, p, width = 12, height = 8, dpi = 300)
  message("Saved plot: ", out_file)
}

for (i in seq_len(nrow(input_tasks))) {
  input_dir <- input_tasks$input_dir[[i]]
  task_name <- input_tasks$task_name[[i]]
  current_files <- list.files(path = input_dir, full.names = TRUE, recursive = TRUE, pattern = "\\.tsv$")

  for (tsv_file in current_files) {
    make_rgreat_plot(tsv_file, task_name = task_name, top_n = top_n)
  }
}

message("Done. Processed ", length(tsv_files), " TSV file(s).")
