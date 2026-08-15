source(file.path("R", "_helpers.R"))
ensure_dirs()

library(DESeq2)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

COUNT_PATHS <- c(
  Aj = Sys.getenv("AJ_RAW_COUNTS", unset = file.path("data", "expression", "Aj_raw_counts.tsv.gz")),
  Ce = Sys.getenv("CE_RAW_COUNTS", unset = file.path("data", "expression", "Ce_raw_counts.tsv.gz"))
)

fig4 <- extract_functional_source_tables()
go <- fig4$go |>
  dplyr::mutate(
    GeneRatio = as.numeric(GeneRatio),
    Q_value = as.numeric(Q_value),
    Count = as.numeric(Count),
    logQ = -log10(Q_value)
  )
genesets <- fig4$genes

p_go <- ggplot(go, aes(GeneRatio, reorder(GO_term, logQ))) +
  geom_point(aes(size = Count, fill = logQ), shape = 21) +
  facet_grid(Module ~ Species, scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "grey90", high = "grey20") +
  labs(x = "Gene ratio", y = NULL, fill = expression(-log[10](Q)), size = "Gene count") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 7))

ggsave(file.path("figures", "Fig4ac_GO_enrichment.pdf"), p_go, width = 7.6, height = 6.6)

compute_subset_pc1 <- function(species, subset_label, trait) {
  meta <- species_metadata(species, complete_traits = trait)
  counts <- read_count_matrix(COUNT_PATHS[[species]])
  counts <- match_matrix_to_metadata(counts, meta)
  rownames(meta) <- meta$Sample

  ids <- genesets |>
    dplyr::filter(Species == species, Functional_subset == subset_label) |>
    dplyr::pull(Unigene_ID) |>
    unique()
  ids <- intersect(ids, rownames(counts))
  if (length(ids) < 2) stop("Too few Fig. 4 functional-subset genes were found for ", species)

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Group_Combo)
  vst_mat <- assay(vst(dds, blind = FALSE))
  expr_sub <- t(vst_mat[ids, , drop = FALSE])

  pca <- prcomp(expr_sub, center = TRUE, scale. = TRUE)
  pc1 <- pca$x[, 1]
  mean_expr <- rowMeans(expr_sub)
  if (cor(pc1, mean_expr) < 0) pc1 <- -pc1

  out <- data.frame(
    Species = species,
    Sample = rownames(expr_sub),
    Functional_subset = subset_label,
    PC1 = pc1,
    Trait = trait,
    Trait_value = meta[rownames(expr_sub), trait],
    Warming = meta[rownames(expr_sub), "Warming"],
    Watering = meta[rownames(expr_sub), "Watering"]
  )

  ct <- cor.test(out$PC1, out$Trait_value, method = "pearson")
  attr(out, "r") <- unname(ct$estimate)
  attr(out, "p") <- ct$p.value
  attr(out, "n_genes") <- length(ids)
  out
}

aj <- compute_subset_pc1("Aj", "ribosome", "BMS_var")
ce <- compute_subset_pc1("Ce", "response to heat + protein folding", "Chl")
source_reg <- dplyr::bind_rows(aj, ce)
readr::write_csv(source_reg, file.path("results", "Fig4bd_functional_subset_PC1.csv"))

plot_regression <- function(d, y_label, tag) {
  r <- attr(d, "r")
  p <- attr(d, "p")
  n_genes <- attr(d, "n_genes")
  lab <- sprintf("Pearson r = %.2f; P = %.3g; genes = %d", r, p, n_genes)

  ggplot(d, aes(PC1, Trait_value)) +
    geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 0.6) +
    geom_point(aes(shape = Watering), size = 2.0) +
    annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.02, vjust = 1.2, size = 3) +
    labs(x = "Functional-subset eigengene (PC1)", y = y_label, tag = tag) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank())
}

p_aj <- plot_regression(aj, "Biomass accumulation", "b")
p_ce <- plot_regression(ce, "Chlorophyll content", "d")
p_reg <- p_aj | p_ce

ggsave(file.path("figures", "Fig4bd_functional_subset_regressions.pdf"), p_reg, width = 7.4, height = 3.4)
