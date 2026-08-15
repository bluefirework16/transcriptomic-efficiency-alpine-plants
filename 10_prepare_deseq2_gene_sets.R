source(file.path("R", "_helpers.R"))
ensure_dirs()

library(DESeq2)
library(readr)
library(dplyr)
library(purrr)
library(tibble)

COUNT_PATHS <- c(
  Aj = Sys.getenv("AJ_RAW_COUNTS", unset = file.path("data", "expression", "Aj_raw_counts.tsv.gz")),
  Ce = Sys.getenv("CE_RAW_COUNTS", unset = file.path("data", "expression", "Ce_raw_counts.tsv.gz"))
)

DEG_PADJ <- 0.05
DEG_ABS_LOG2FC <- 1

is_deg <- function(padj, log2fc) {
  !is.na(padj) & padj < DEG_PADJ & !is.na(log2fc) & abs(log2fc) > DEG_ABS_LOG2FC
}

run_gene_sets <- function(species) {
  message("Preparing DESeq2 gene sets for ", species)

  counts <- read_count_matrix(COUNT_PATHS[[species]])
  meta <- species_metadata(species)
  counts <- match_matrix_to_metadata(counts, meta)
  rownames(meta) <- meta$Sample

  # ---------------------------------------------------------------------------
  # A. Pairwise warming and watering responses using treatment combinations
  # ---------------------------------------------------------------------------
  dds_group <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ Group_Combo
  )

  keep <- rowSums(counts(dds_group) >= 10) >= 4
  dds_group <- dds_group[keep, ]
  dds_group <- DESeq(
    dds_group,
    test = "Wald",
    fitType = "parametric",
    sfType = "ratio",
    minReplicatesForReplace = 7
  )

  get_group_contrast <- function(treat, control, response_class) {
    res <- results(
      dds_group,
      contrast = c("Group_Combo", treat, control),
      alpha = DEG_PADJ,
      independentFiltering = TRUE,
      pAdjustMethod = "BH"
    )

    as.data.frame(res) |>
      tibble::rownames_to_column("Gene_ID") |>
      dplyr::mutate(
        Species = species,
        Response_class = response_class,
        Contrast = paste0(treat, "_vs_", control),
        Significant = is_deg(padj, log2FoldChange)
      ) |>
      dplyr::relocate(Species, Response_class, Contrast, Gene_ID)
  }

  warming_grid <- expand.grid(
    Warming = c("T1", "T2", "T3"),
    Watering = c("Ambient", "Added"),
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      Treat = paste(Warming, Watering, sep = "_"),
      Control = paste("T0", Watering, sep = "_")
    )

  warming_results <- purrr::map2_dfr(
    warming_grid$Treat,
    warming_grid$Control,
    ~ get_group_contrast(.x, .y, "warming")
  )

  watering_grid <- data.frame(
    Warming = c("T0", "T1", "T2", "T3"),
    Treat = paste(c("T0", "T1", "T2", "T3"), "Added", sep = "_"),
    Control = paste(c("T0", "T1", "T2", "T3"), "Ambient", sep = "_")
  )

  watering_results <- purrr::map2_dfr(
    watering_grid$Treat,
    watering_grid$Control,
    ~ get_group_contrast(.x, .y, "watering")
  )

  # ---------------------------------------------------------------------------
  # B. Warming × Watering interaction-responsive genes
  # Warming and Watering are categorical factors. The comprehensive DEG set uses
  # the union of significant warming, watering, and interaction-responsive genes.
  # ---------------------------------------------------------------------------
  meta_int <- meta
  meta_int$Warming <- factor(meta_int$Warming, levels = c("T0", "T1", "T2", "T3"))
  meta_int$Watering <- factor(meta_int$Watering, levels = c("Ambient", "Added"))

  dds_int <- DESeqDataSetFromMatrix(
    countData = counts[rownames(dds_group), , drop = FALSE],
    colData = meta_int,
    design = ~ Warming * Watering
  )
  dds_int <- DESeq(
    dds_int,
    test = "Wald",
    fitType = "parametric",
    sfType = "ratio",
    minReplicatesForReplace = 7
  )

  interaction_names <- grep(
    "Warming.*Watering|Watering.*Warming",
    resultsNames(dds_int),
    value = TRUE
  )

  if (length(interaction_names) == 0) {
    stop("No Warming × Watering interaction coefficients were found in DESeq2 resultsNames().")
  }

  interaction_results <- purrr::map_dfr(interaction_names, function(coef_name) {
    res <- results(
      dds_int,
      name = coef_name,
      alpha = DEG_PADJ,
      independentFiltering = TRUE,
      pAdjustMethod = "BH"
    )

    as.data.frame(res) |>
      tibble::rownames_to_column("Gene_ID") |>
      dplyr::mutate(
        Species = species,
        Response_class = "interaction",
        Contrast = coef_name,
        Significant = is_deg(padj, log2FoldChange)
      ) |>
      dplyr::relocate(Species, Response_class, Contrast, Gene_ID)
  })

  # ---------------------------------------------------------------------------
  # C. Gene sets used downstream
  # ---------------------------------------------------------------------------
  all_results <- dplyr::bind_rows(
    warming_results,
    watering_results,
    interaction_results
  )

  readr::write_csv(
    all_results,
    file.path("results", "expression", paste0(species, "_DESeq2_all_tests.csv"))
  )

  warming_ids <- warming_results |>
    dplyr::filter(Significant) |>
    dplyr::pull(Gene_ID) |>
    unique()

  comprehensive_ids <- all_results |>
    dplyr::filter(Significant) |>
    dplyr::pull(Gene_ID) |>
    unique()

  writeLines(
    warming_ids,
    file.path("results", "expression", paste0(species, "_warming_DEG_ids.txt"))
  )
  writeLines(
    comprehensive_ids,
    file.path("results", "expression", paste0(species, "_comprehensive_DEG_ids.txt"))
  )

  data.frame(
    Species = species,
    Warming_responsive_genes = length(warming_ids),
    Comprehensive_DEG_set = length(comprehensive_ids)
  )
}

summary_table <- dplyr::bind_rows(lapply(c("Aj", "Ce"), run_gene_sets))
readr::write_csv(
  summary_table,
  file.path("results", "expression", "DEG_set_summary.csv")
)
