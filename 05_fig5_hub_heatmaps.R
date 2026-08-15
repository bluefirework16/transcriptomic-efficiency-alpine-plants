source(file.path("R", "_helpers.R"))
ensure_dirs()

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

hub <- read_source_table(source_sheet("hub_expression"))
expr_cols <- c("Ambient_T0", "Ambient_T1", "Ambient_T2", "Ambient_T3", "Added_T0", "Added_T1", "Added_T2", "Added_T3")

mat <- as.matrix(hub[, expr_cols])
storage.mode(mat) <- "numeric"
z <- zscore_rows(mat)

zdf <- cbind(hub[, c("Species", "Unigene_ID", "Gene_symbol", "Annotation")], as.data.frame(z)) |>
  tidyr::pivot_longer(dplyr::all_of(expr_cols), names_to = "Treatment", values_to = "Z_score") |>
  dplyr::mutate(Treatment = factor(Treatment, levels = expr_cols))

readr::write_csv(zdf, file.path("results", "Fig5_hub_expression_zscores.csv"))

plot_species <- function(sp) {
  d <- zdf |> dplyr::filter(Species == sp)
  d$Gene_label <- factor(d$Gene_symbol, levels = rev(unique(d$Gene_symbol)))
  ggplot(d, aes(Treatment, Gene_label, fill = Z_score)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "grey20", mid = "white", high = "grey80", midpoint = 0) +
    labs(title = sp, x = NULL, y = NULL, fill = "Z score") +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
}

p <- plot_species("Aj") | plot_species("Ce")
ggsave(file.path("figures", "Fig5_hub_gene_heatmaps.pdf"), p, width = 10, height = 5.5)

message(
  "Note: the arrows and pathway topology in Fig. 5 are a literature-supported working model, ",
  "not an inferred causal network. This script reproduces the expression heatmaps only."
)
