#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE)
analysis_dir <- dirname(dirname(script_path))
table_dir <- file.path(analysis_dir, "results", "tables")
figure_dir <- file.path(analysis_dir, "results", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Source Sans 3"

developmental <- read.delim(
  file.path(table_dir, "parikh_pgtd_developmental_expression.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
prespore <- read.delim(
  file.path(table_dir, "parikh_pgtd_prespore_enrichment.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

species_order <- c("D. discoideum", "D. purpureum")
developmental$species <- factor(developmental$species, levels = species_order)
prespore$species <- factor(prespore$species, levels = species_order)

theme_manuscript <- function() {
  ggplot2::theme_classic(base_size = 9, base_family = font_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = font_family, colour = "black"),
      axis.title = ggplot2::element_text(family = font_family, size = 9),
      axis.text = ggplot2::element_text(family = font_family, size = 8),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.30),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.30),
      strip.background = ggplot2::element_blank(),
      strip.text = ggtext::element_markdown(
        family = font_family,
        size = 9,
        face = "plain",
        colour = "black",
        margin = ggplot2::margin(b = 4)
      ),
      legend.position = "none",
      plot.title = ggplot2::element_blank()
    )
}

# Panel B: developmental trajectories and species means.
developmental_mean <- aggregate(
  log2_scaled_abundance ~ species + time_hour,
  data = developmental,
  FUN = mean
)
panel_b_labels <- c(
  "D. discoideum" = "<i>D. discoideum</i> (AX4)",
  "D. purpureum" = "<i>D. purpureum</i> (DpAX1)"
)

p_dev <- ggplot2::ggplot(
  developmental,
  ggplot2::aes(
    x = time_hour,
    y = log2_scaled_abundance,
    group = biological_replicate
  )
) +
  ggplot2::geom_line(colour = "grey70", linewidth = 0.35) +
  ggplot2::geom_point(colour = "grey70", size = 1.15) +
  ggplot2::geom_line(
    data = developmental_mean,
    ggplot2::aes(x = time_hour, y = log2_scaled_abundance),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.80
  ) +
  ggplot2::geom_point(
    data = developmental_mean,
    ggplot2::aes(x = time_hour, y = log2_scaled_abundance),
    inherit.aes = FALSE,
    colour = "black",
    size = 1.45
  ) +
  ggplot2::facet_wrap(
    ~species,
    nrow = 1,
    labeller = ggplot2::as_labeller(panel_b_labels)
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 24, 4),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.04, 0.07))
  ) +
  ggplot2::labs(
    x = "Hours after starvation",
    y = "log₂(scaled abundance + 1)"
  ) +
  theme_manuscript() +
  ggplot2::theme(panel.spacing.x = grid::unit(0.65, "cm"))

# Panel C: each line joins the two cell types from one biological replicate.
prespore_long <- rbind(
  data.frame(
    species = prespore$species,
    biological_replicate = prespore$biological_replicate,
    cell_type = "Prestalk",
    scaled_abundance = prespore$prestalk_scaled_abundance
  ),
  data.frame(
    species = prespore$species,
    biological_replicate = prespore$biological_replicate,
    cell_type = "Prespore",
    scaled_abundance = prespore$prespore_scaled_abundance
  )
)
prespore_long$species <- factor(prespore_long$species, levels = species_order)
prespore_long$cell_type <- factor(
  prespore_long$cell_type,
  levels = c("Prestalk", "Prespore")
)
prespore_long$log2_scaled_abundance <- log2(prespore_long$scaled_abundance + 1)

format_fold <- function(value) {
  sprintf("%.2g", value)
}
format_p_value <- function(value) {
  formatC(value, format = "g", digits = 3)
}

global_celltype_y_range <- range(prespore_long$log2_scaled_abundance)
annotation_y <- global_celltype_y_range[2L] +
  0.10 * diff(global_celltype_y_range)

prespore_annotations <- do.call(rbind, lapply(species_order, function(species_name) {
  pair_rows <- prespore[as.character(prespore$species) == species_name, , drop = FALSE]
  fold_range <- range(pair_rows$prespore_to_prestalk_ratio)
  published_p <- unique(pair_rows$whole_transcript_ebayes_p_value)
  data.frame(
    species = species_name,
    x = 1.5,
    y = annotation_y,
    label = sprintf(
      "atop(\"%s–%s-fold\", italic(P) == %s)",
      format_fold(fold_range[1L]),
      format_fold(fold_range[2L]),
      format_p_value(published_p)
    ),
    stringsAsFactors = FALSE
  )
}))
prespore_annotations$species <- factor(
  prespore_annotations$species,
  levels = species_order
)

panel_c_labels <- c(
  "D. discoideum" = "<i>D. discoideum</i> (NC4)",
  "D. purpureum" = "<i>D. purpureum</i> (DpAX1)"
)

p_prespore <- ggplot2::ggplot(
  prespore_long,
  ggplot2::aes(
    x = cell_type,
    y = log2_scaled_abundance,
    group = biological_replicate
  )
) +
  ggplot2::geom_line(colour = "black", linewidth = 0.45) +
  ggplot2::geom_point(colour = "black", size = 1.70) +
  ggplot2::geom_text(
    data = prespore_annotations,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    parse = TRUE,
    vjust = 1,
    family = font_family,
    size = 2.65,
    lineheight = 1.0,
    colour = "black"
  ) +
  ggplot2::facet_wrap(
    ~species,
    nrow = 1,
    labeller = ggplot2::as_labeller(panel_c_labels)
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.05, 0.12))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "log₂(scaled abundance + 1)"
  ) +
  theme_manuscript() +
  ggplot2::theme(
    panel.spacing.x = grid::unit(1.10, "cm"),
    axis.text.x = ggplot2::element_text(family = font_family, size = 8)
  )

# Export at final physical dimensions.
save_panel <- function(plot, name, height) {
  ggplot2::ggsave(
    file.path(figure_dir, paste0(name, ".svg")), plot,
    width = 7.5, height = height, units = "in",
    device = svglite::svglite, bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(name, ".png")), plot,
    width = 7.5, height = height, units = "in",
    dpi = 600, device = "png", bg = "white"
  )
}
save_panel(p_dev, "panelB", 2.1)
save_panel(p_prespore, "panelC", 1.6)

cat(
  "\nPlotting: scripts/02_make_expression_panels.R\n",
  paste(capture.output(sessionInfo()), collapse = "\n"), "\n",
  file = file.path(analysis_dir, "results", "session_info.txt"),
  append = TRUE
)

message("Wrote panelB and panelC as svg and png files to ", figure_dir)
