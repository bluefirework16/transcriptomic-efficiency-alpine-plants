source(file.path("R", "_helpers.R"))
ensure_dirs()

library(WGCNA)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

options(stringsAsFactors = FALSE)
WGCNA::allowWGCNAThreads()

EXPR_PATHS <- c(
  Aj = Sys.getenv("AJ_FPKM", unset = file.path("data", "expression", "Aj_FPKM.tsv.gz")),
  Ce = Sys.getenv("CE_FPKM", unset = file.path("data", "expression", "Ce_FPKM.tsv.gz"))
)

trait_map <- c(
  Biomass = "BMS_var",
  Photo = "Photo",
  Cond = "Cond",
  Chl = "Chl",
  MDA = "MDA",
  POD = "POD",
  APX = "APX",
  GR = "GR"
)

run_wgcna <- function(species) {
  message("Running WGCNA for ", species)

  gene_file <- file.path(
    "results", "expression", paste0(species, "_comprehensive_DEG_ids.txt")
  )
  assert_file(
    gene_file,
    "Comprehensive DEG list; run R/10_prepare_deseq2_gene_sets.R first"
  )
  comprehensive_genes <- unique(readLines(gene_file))

  expr <- read_expression_matrix(EXPR_PATHS[[species]])
  meta <- species_metadata(species)
  expr <- match_matrix_to_metadata(expr, meta)

  comprehensive_genes <- intersect(comprehensive_genes, rownames(expr))
  if (length(comprehensive_genes) < 30) {
    stop("Too few comprehensive DEGs are available for WGCNA in ", species)
  }

  # Rows = samples, columns = genes. The deposited normalized expression values
  # are used directly, consistent with the manuscript Methods.
  datExpr <- as.data.frame(t(expr[comprehensive_genes, , drop = FALSE]))

  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }

  # Supplementary Fig. S4: scale-free topology diagnostic.
  powers <- c(1:10, seq(12, 30, by = 2))
  sft <- WGCNA::pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = "unsigned",
    corFnc = "cor",
    corOptions = "use = 'p'",
    verbose = 0
  )

  pdf(
    file.path("figures", paste0("FigS4_", species, "_soft_threshold.pdf")),
    width = 7.2, height = 3.5
  )
  par(mfrow = c(1, 2))
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft threshold (power)",
    ylab = "Scale-free topology model fit, signed R²",
    type = "n",
    main = paste(species, "scale independence")
  )
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = 0.75
  )
  abline(h = 0.90, lty = 2)
  plot(
    sft$fitIndices[, 1], sft$fitIndices[, 5],
    xlab = "Soft threshold (power)",
    ylab = "Mean connectivity",
    type = "n",
    main = paste(species, "mean connectivity")
  )
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.75)
  dev.off()

  # Main network used in the manuscript.
  net <- WGCNA::blockwiseModules(
    datExpr,
    power = 6,
    maxBlockSize = 5000,
    TOMType = "unsigned",
    networkType = "unsigned",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.20,
    numericLabels = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    corType = "pearson",
    verbose = 2
  )

  module_colors <- as.character(net$colors)
  names(module_colors) <- colnames(datExpr)
  MEs <- WGCNA::orderMEs(net$MEs)

  assignments <- data.frame(
    Species = species,
    Gene_ID = names(module_colors),
    Module = unname(module_colors)
  )
  readr::write_csv(
    assignments,
    file.path("results", paste0("Fig3_", species, "_module_assignments.csv"))
  )

  module_counts <- assignments |>
    dplyr::count(Species, Module, name = "Gene_count")
  readr::write_csv(
    module_counts,
    file.path("results", paste0("Fig3_", species, "_module_counts.csv"))
  )

  traits <- meta |>
    dplyr::select(Sample, dplyr::all_of(unname(trait_map)))
  rownames(traits) <- traits$Sample
  traits$Sample <- NULL
  traits <- traits[rownames(MEs), , drop = FALSE]
  colnames(traits) <- names(trait_map)

  # Pearson correlation is used throughout the module–trait analysis.
  cor_rows <- list()
  k <- 1L
  for (module in colnames(MEs)) {
    for (trait in colnames(traits)) {
      keep <- stats::complete.cases(MEs[, module], traits[, trait])
      ct <- stats::cor.test(
        MEs[keep, module], traits[keep, trait],
        method = "pearson"
      )
      cor_rows[[k]] <- data.frame(
        Species = species,
        Module = module,
        Trait = trait,
        Pearson_r = unname(ct$estimate),
        P_value = ct$p.value,
        N = sum(keep)
      )
      k <- k + 1L
    }
  }

  cor_tab <- dplyr::bind_rows(cor_rows)
  readr::write_csv(
    cor_tab,
    file.path("results", paste0("Fig3_", species, "_module_trait.csv"))
  )

  pdf(
    file.path("figures", paste0("FigS4_", species, "_WGCNA_dendrogram.pdf")),
    width = 8, height = 5
  )
  WGCNA::plotDendroAndColors(
    net$dendrograms[[1]],
    module_colors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  dev.off()

  cor_tab
}

computed <- dplyr::bind_rows(lapply(c("Aj", "Ce"), run_wgcna))
readr::write_csv(computed, file.path("results", "Fig3_module_trait_recomputed.csv"))

# Redraw the final heatmap directly from Additional file 2.
source_cor <- read_source_table(source_sheet("module_trait")) |>
  dplyr::mutate(
    Module = factor(Module, levels = rev(unique(Module))),
    Trait = factor(Trait, levels = c("Biomass", "Photo", "Cond", "Chl", "MDA", "POD", "APX", "GR")),
    sig = dplyr::case_when(
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p <- ggplot(source_cor, aes(Trait, Module, fill = Pearson_r)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = paste0(sprintf("%.2f", Pearson_r), sig)), size = 2.6) +
  facet_wrap(~ Species, scales = "free_y", ncol = 1) +
  scale_fill_gradient2(
    low = "grey20", mid = "white", high = "grey75",
    midpoint = 0, limits = c(-1, 1)
  ) +
  labs(x = NULL, y = NULL, fill = "Pearson r") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path("figures", "Fig3_module_trait_heatmap_from_source.pdf"),
  p, width = 6.5, height = 7.5
)
