suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

root <- if (dir.exists("S19S_KGML_reproducible_workflow")) {
  "S19S_KGML_reproducible_workflow"
} else {
  "."
}

derived <- file.path(root, "derived")
figures <- file.path(root, "figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

kgml <- read_csv(
  file.path(derived, "section34_vh0_025_050_100_importance_family.csv"),
  show_col_types = FALSE
) %>%
  filter(Vh0 == 0.5) %>%
  transmute(
    Predictor_set,
    Model = "KGML",
    Importance_method = "Permutation importance",
    Feature_family,
    family_share
  )

pure <- read_csv(
  file.path(derived, "section34_pure_ml_importance_family.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    Predictor_set,
    Model = recode(Model,
                   `Pure RF` = "Pure RF",
                   `Pure XGBoost` = "Pure XGBoost"),
    Importance_method = recode(Importance_type,
                               Permutation = "Permutation importance",
                               Gain = "Gain importance"),
    Feature_family,
    family_share
  )

plot_data <- bind_rows(kgml, pure) %>%
  mutate(
    Predictor_set = factor(Predictor_set,
                           levels = c("Bulk", "FTICR", "Bulk+FTICR"),
                           labels = c("Bulk", "FTICR-MS", "Bulk + FTICR-MS")),
    Model = factor(Model, levels = c("KGML", "Pure RF", "Pure XGBoost")),
    Feature_family = factor(
      Feature_family,
      levels = c("NPOC", "Thermodynamic", "Sediment physical structure",
                 "Biomass", "OM molecular composition")
    )
  )

stopifnot(
  nrow(plot_data) > 0,
  all(abs(plot_data %>%
            group_by(Predictor_set, Model) %>%
            summarise(total = sum(family_share), .groups = "drop") %>%
            pull(total) - 1) < 1e-6)
)

family_colors <- c(
  "NPOC" = "#2d004b",
  "Thermodynamic" = "#542788",
  "Sediment physical structure" = "#8073ac",
  "Biomass" = "#b2abd2",
  "OM molecular composition" = "#e08214"
)

p <- ggplot(plot_data, aes(x = Model, y = family_share, fill = Feature_family)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  facet_wrap(~ Predictor_set, nrow = 1) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = family_colors, drop = FALSE) +
  labs(
    x = NULL,
    y = "Relative feature-family importance",
    fill = "Predictor family"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10.5),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    panel.spacing.x = unit(0.8, "lines"),
    plot.margin = margin(8, 8, 6, 8)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

write_csv(
  plot_data %>% mutate(across(where(is.factor), as.character)),
  file.path(derived, "supplement_feature_family_importance_plot_data.csv")
)

ggsave(
  file.path(figures, "supplement_feature_family_importance_comparison.pdf"),
  p, width = 8.2, height = 4.6, device = cairo_pdf
)
ggsave(
  file.path(figures, "supplement_feature_family_importance_comparison.png"),
  p, width = 8.2, height = 4.6, dpi = 400, bg = "white"
)

message("Wrote supplementary feature-family importance figure and plot data.")
