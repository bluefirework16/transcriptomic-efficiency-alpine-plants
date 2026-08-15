source(file.path("R", "_helpers.R"))
ensure_dirs()

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

clim <- read_source_table(source_sheet("climate"))

fit_temp <- lm(Mean_annual_air_temperature ~ Year, data = clim)
fit_prec <- lm(Growing_season_precipitation_JJA ~ Year, data = clim)

model_summary <- dplyr::bind_rows(
  data.frame(
    response = "Mean annual air temperature",
    slope_per_year = coef(fit_temp)[["Year"]],
    slope_per_decade = 10 * coef(fit_temp)[["Year"]],
    r_squared = summary(fit_temp)$r.squared,
    p_value = coef(summary(fit_temp))["Year", "Pr(>|t|)"]
  ),
  data.frame(
    response = "Growing-season precipitation (JJA)",
    slope_per_year = coef(fit_prec)[["Year"]],
    slope_per_decade = 10 * coef(fit_prec)[["Year"]],
    r_squared = summary(fit_prec)$r.squared,
    p_value = coef(summary(fit_prec))["Year", "Pr(>|t|)"]
  )
)
readr::write_csv(model_summary, file.path("results", "Fig6_climate_regressions.csv"))

p_temp <- ggplot(clim, aes(Year, Mean_annual_air_temperature)) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.4) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  labs(x = "Year", y = "Mean annual air temperature (°C)", tag = "b") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())

p_prec <- ggplot(clim, aes(Year, Growing_season_precipitation_JJA)) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.4) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  labs(x = "Year", y = "Growing-season precipitation (mm)", tag = "c") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())

p <- p_temp | p_prec
ggsave(file.path("figures", "Fig6bc_climate_trends.pdf"), p, width = 7.2, height = 3.2)
ggsave(file.path("figures", "Fig6bc_climate_trends.png"), p, width = 7.2, height = 3.2, dpi = 300)
