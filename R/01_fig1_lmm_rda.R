source(file.path("R", "_helpers.R"))
ensure_dirs()

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(vegan)

# Type III tests require sum-to-zero contrasts for categorical predictors.
options(contrasts = c("contr.sum", "contr.poly"))

env <- prepare_environment_source()
plant <- prepare_plant_source()

# =============================================================================
# 1. Environmental mixed-effects models
# Warming and Watering are categorical variables.
# Model: Y ~ Warming * Watering + (1 | Block)
# =============================================================================

env_vars <- c("ST", "AT", "VWC")

env_lmm <- purrr::map_dfr(env_vars, function(v) {
  d <- env |>
    dplyr::select(Block, Warming, Watering, all_of(v)) |>
    tidyr::drop_na()

  fit <- lmerTest::lmer(
    stats::as.formula(paste0(v, " ~ Warming * Watering + (1 | Block)")),
    data = d,
    REML = FALSE
  )

  clean_lmer_anova(
    lmerTest::anova(fit, type = 3),
    trait = v,
    analysis = "Environmental LMM"
  )
})
readr::write_csv(env_lmm, file.path("results", "LMM_environment_recomputed.csv"))

# =============================================================================
# 2. Within-species physiological mixed-effects models
# Warming and Watering are categorical variables.
# Model: Y ~ Warming * Watering + (1 | Block)
# =============================================================================

traits <- c(
  "H_var", "BMS_var", "Photo", "Cond", "Ci", "Tr", "Chl",
  "MDA", "POD", "APX", "CAT", "GR", "DR"
)

within_lmm <- purrr::map_dfr(c("Aj", "Ce"), function(sp) {
  purrr::map_dfr(traits, function(trait) {
    d <- plant |>
      dplyr::filter(Species == sp) |>
      dplyr::select(Block, Warming, Watering, all_of(trait)) |>
      tidyr::drop_na()

    fit <- lmerTest::lmer(
      stats::as.formula(paste0(trait, " ~ Warming * Watering + (1 | Block)")),
      data = d,
      REML = FALSE
    )

    clean_lmer_anova(
      lmerTest::anova(fit, type = 3),
      trait = trait,
      species = sp,
      analysis = "Within-species LMM"
    )
  })
})
readr::write_csv(within_lmm, file.path("results", "LMM_within_species_recomputed.csv"))

# T1/T2/T3 versus the corresponding T0 control within each watering regime.
within_emmeans <- purrr::map_dfr(c("Aj", "Ce"), function(sp) {
  purrr::map_dfr(traits, function(trait) {
    d <- plant |>
      dplyr::filter(Species == sp) |>
      dplyr::select(Block, Warming, Watering, all_of(trait)) |>
      tidyr::drop_na()

    fit <- lmerTest::lmer(
      stats::as.formula(paste0(trait, " ~ Warming * Watering + (1 | Block)")),
      data = d,
      REML = FALSE
    )

    emm <- emmeans::emmeans(fit, ~ Warming | Watering)
    emmeans::contrast(emm, method = "trt.vs.ctrl", ref = "T0", adjust = "fdr") |>
      summary(infer = TRUE) |>
      as.data.frame() |>
      dplyr::mutate(Species = sp, Trait = trait, .before = 1)
  })
})
readr::write_csv(within_emmeans, file.path("results", "LMM_within_species_T_vs_T0_emmeans.csv"))

# =============================================================================
# 3. Combined-species mixed-effects models
# Model: Y ~ Species * Warming * Watering + (1 | Block) + (1 | Sample)
# The Sample random intercept accounts for Aj and Ce measured in the same plot.
# =============================================================================

combined_lmm <- purrr::map_dfr(traits, function(trait) {
  d <- plant |>
    dplyr::select(Sample, Block, Species, Warming, Watering, all_of(trait)) |>
    tidyr::drop_na()

  fit <- lmerTest::lmer(
    stats::as.formula(
      paste0(trait, " ~ Species * Warming * Watering + (1 | Block) + (1 | Sample)")
    ),
    data = d,
    REML = FALSE
  )

  out <- clean_lmer_anova(
    lmerTest::anova(fit, type = 3),
    trait = trait,
    analysis = "Combined-species LMM"
  )

  out$Trait <- dplyr::recode(
    out$Trait,
    H_var = "HT_var",
    BMS_var = "Biomass_Var"
  )
  out
})
readr::write_csv(combined_lmm, file.path("results", "LMM_combined_species_recomputed.csv"))

# Supplementary Fig. S2: combined-species Type III ANOVA heatmap.
term_order <- c(
  "Species", "Warming", "Watering", "Warming:Watering",
  "Species:Warming", "Species:Watering", "Species:Warming:Watering"
)

trait_info <- tibble::tibble(
  Trait = c(
    "HT_var", "Biomass_Var", "Photo", "Cond", "Ci", "Tr", "Chl", "MDA",
    "POD", "APX", "CAT", "GR", "DR"
  ),
  Trait_label = c(
    "Height variation", "Biomass accumulation", "Photosynthetic rate",
    "Stomatal conductance", "Intercellular CO2", "Transpiration rate",
    "Chlorophyll", "MDA", "POD", "APX", "CAT", "GR", "DR"
  ),
  Trait_group = c(
    "Growth", "Growth", rep("Gas exchange", 4),
    rep("Pigment and oxidative status", 2), rep("Antioxidant enzymes", 5)
  )
)

heat_df <- combined_lmm |>
  dplyr::filter(Term %in% term_order) |>
  dplyr::left_join(trait_info, by = "Trait") |>
  dplyr::mutate(
    neg_log10_p = pmin(-log10(P_value), 8),
    sig = dplyr::case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Term = factor(Term, levels = term_order),
    Trait_label = factor(Trait_label, levels = rev(trait_info$Trait_label)),
    Trait_group = factor(
      Trait_group,
      levels = c("Growth", "Gas exchange", "Pigment and oxidative status", "Antioxidant enzymes")
    )
  )

p_heat <- ggplot(heat_df, aes(Term, Trait_label, fill = neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig), fontface = "bold", size = 3.4) +
  facet_grid(Trait_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient(
    low = "grey95", high = "grey25", limits = c(0, 8),
    name = expression(-log[10](italic(P)))
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_rect(fill = "grey90")
  )

ggsave(
  file.path("figures", "FigS2_combined_species_LMM_heatmap.pdf"),
  p_heat, width = 8.8, height = 6.8
)

# =============================================================================
# 4. Fig. 1: microclimate, physiological responses, and RDA
# =============================================================================

env_long <- env |>
  dplyr::select(Sample, Block, Warming, Watering, ST, AT, VWC) |>
  tidyr::pivot_longer(c(ST, AT, VWC), names_to = "Variable", values_to = "Value")

p_env <- ggplot(env_long, aes(Warming, Value, fill = Watering)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.62, outlier.shape = NA) +
  facet_wrap(~ Variable, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL, tag = "a") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), legend.position = "top")

plot_trait <- function(trait, label, tag) {
  ggplot(plant, aes(Warming, .data[[trait]], fill = Watering)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.62, outlier.shape = NA) +
    facet_wrap(~ Species, nrow = 1, scales = "free_y") +
    labs(x = "Warming", y = label, tag = tag) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), legend.position = "none")
}

p_bms <- plot_trait("BMS_var", "Biomass accumulation", "b")
p_photo <- plot_trait("Photo", "Net photosynthetic rate", "c")
p_apx <- plot_trait("APX", "APX activity", "d")
p_mda <- plot_trait("MDA", "MDA content", "e")

rda_traits <- c("BMS_var", "Photo", "Cond", "Chl", "MDA", "POD", "APX", "GR")
rda_dat <- plant |>
  dplyr::left_join(env |> dplyr::select(Sample, ST, AT, VWC), by = "Sample") |>
  dplyr::select(Species, Watering, all_of(rda_traits), ST, AT, VWC) |>
  tidyr::drop_na()

Y <- rda_dat |> dplyr::select(all_of(rda_traits))
X <- rda_dat |> dplyr::select(ST, AT, VWC)
rda_model <- vegan::rda(Y ~ ST + AT + VWC, data = X, scale = TRUE)

site_scores <- as.data.frame(vegan::scores(rda_model, display = "sites", choices = 1:2)) |>
  dplyr::mutate(Species = rda_dat$Species, Watering = rda_dat$Watering)
env_scores <- as.data.frame(vegan::scores(rda_model, display = "bp", choices = 1:2))
trait_scores <- as.data.frame(vegan::scores(rda_model, display = "species", choices = 1:2))

limit_val <- max(abs(c(site_scores$RDA1, site_scores$RDA2)))
arrow_max <- max(abs(c(env_scores$RDA1, env_scores$RDA2, trait_scores$RDA1, trait_scores$RDA2)))
scale_factor <- ifelse(arrow_max == 0, 1, 0.8 * limit_val / arrow_max)
env_scores <- env_scores * scale_factor
trait_scores <- trait_scores * scale_factor
env_scores$Label <- rownames(env_scores)
trait_scores$Label <- rownames(trait_scores)

rda_eig <- rda_model$CCA$eig
var_percent <- 100 * rda_eig / sum(rda_model$CA$eig + rda_eig)

p_rda <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey75") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey75") +
  stat_ellipse(
    data = site_scores,
    aes(RDA1, RDA2, color = Species, group = Species),
    level = 0.95, linewidth = 0.5
  ) +
  geom_point(
    data = site_scores,
    aes(RDA1, RDA2, color = Species, shape = Watering),
    size = 1.8
  ) +
  geom_segment(
    data = env_scores,
    aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
    arrow = arrow(length = grid::unit(0.18, "cm")), linewidth = 0.5, color = "grey40"
  ) +
  geom_text_repel(data = env_scores, aes(RDA1, RDA2, label = Label), size = 2.8, color = "grey30") +
  geom_segment(
    data = trait_scores,
    aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
    arrow = arrow(length = grid::unit(0.15, "cm")), linewidth = 0.45
  ) +
  geom_text_repel(data = trait_scores, aes(RDA1, RDA2, label = Label), size = 2.6) +
  labs(
    x = sprintf("RDA1 (%.1f%%)", var_percent[1]),
    y = sprintf("RDA2 (%.1f%%)", var_percent[2]),
    tag = "f"
  ) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank())

p_fig1 <- p_env / ((p_bms | p_photo) / (p_apx | p_mda)) / p_rda +
  patchwork::plot_layout(heights = c(0.8, 1.6, 1.0))

ggsave(file.path("figures", "Fig1_reproduced.pdf"), p_fig1, width = 8.4, height = 10.0)
readr::write_csv(site_scores, file.path("results", "Fig1_RDA_site_scores.csv"))
readr::write_csv(env_scores, file.path("results", "Fig1_RDA_environment_scores.csv"))
readr::write_csv(trait_scores, file.path("results", "Fig1_RDA_trait_scores.csv"))

# Supplementary Fig. S1: additional physiological and biochemical traits.
s1_traits <- c(
  Cond = "Stomatal conductance",
  GR = "GR activity",
  POD = "POD activity",
  Chl = "Chlorophyll content"
)

s1_plots <- purrr::imap(s1_traits, function(label, trait) {
  ggplot(plant, aes(Warming, .data[[trait]], fill = Watering)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.62, outlier.shape = NA) +
    facet_wrap(~ Species, nrow = 1, scales = "free_y") +
    labs(x = "Warming", y = label) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), legend.position = "top")
})

p_s1 <- (s1_plots[[1]] | s1_plots[[2]]) /
  (s1_plots[[3]] | s1_plots[[4]])

ggsave(
  file.path("figures", "FigS1_additional_physiology.pdf"),
  p_s1, width = 8.0, height = 6.0
)
