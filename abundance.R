# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            abundance.R                            ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : fungi-ITS-TEF1                                   ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-29                                       ║
# ║ Last Modified  : 2025-10-29                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

# Libraries

library(magrittr, include.only = "%<>%")
library(qiime2R)
library(readr)
library(tidyverse)
library(EnhancedVolcano)
library(cowplot)
library(patchwork)

## Colors and shapes

source("/home/sergio/projects/fungi-ITS-TEF1/colors.R")

## Functions

source("/home/sergio/projects/diversity-cereal/parse_ancombc.R")

volcanoFromAncombc <- function(qza_path,
                               log2fc_col,
                               pval_col,
                               up_color,
                               down_color,
                               up_shape,
                               down_shape,
                               up_legend,
                               down_legend,
                               ...,
                               lab_col = NA,
                               log2fc_cutoff = 2,
                               pval_cutoff = 0.05,
                               ns_color = "grey45",
                               ns_shape = 4,
                               ns_legend = "NS",
                               taxonomy_df = taxonomy)
{
  ancombc <- import_ancombc(qza_path)
  ancombc %<>% left_join(taxonomy_df)
  ancombc %<>% mutate(
    labnames = coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom),
    labnames = if_else(
      labnames %in% c(Species, Genus),
      paste0("italic('", labnames, "')"),
      paste0("'", labnames, "'")),
    labnames = gsub("_", " ", labnames))
  
  keyvals_col <- ifelse(
    ancombc[[log2fc_col]] < -log2fc_cutoff & ancombc[[pval_col]] < pval_cutoff, down_color,
    ifelse(ancombc[[log2fc_col]] > log2fc_cutoff & ancombc[[pval_col]] < pval_cutoff, up_color,
           ns_color))
  keyvals_col[is.na(keyvals_col)] <- ns_color
  names(keyvals_col)[keyvals_col == up_color] <- up_legend
  names(keyvals_col)[keyvals_col == ns_color] <- ns_legend
  names(keyvals_col)[keyvals_col == down_color] <- down_legend
  
  keyvals_shape <- keyvals_col
  keyvals_shape[keyvals_shape == up_color] <- up_shape
  keyvals_shape[keyvals_shape == ns_color] <- ns_shape
  keyvals_shape[keyvals_shape == down_color] <- down_shape
  keyvals_shape %<>% as.integer()
  names(keyvals_shape) <- names(keyvals_col)
  
  lab_arg <- if (is.na(lab_col)) NA else ancombc[[lab_col]]
  
  v_plot <- ancombc %>%
    EnhancedVolcano(lab = lab_arg,
                    x = {{log2fc_col}},
                    y = {{pval_col}},
                    pCutoff = pval_cutoff,
                    FCcutoff = log2fc_cutoff,
                    colCustom = keyvals_col,
                    shapeCustom = keyvals_shape,
                    ...) +
    guides(color = guide_legend("Combined Legend",
                                override.aes = list(alpha=1)),
           shape = guide_legend("Combined Legend")) +
    theme_classic() +
    theme(legend.title=element_blank(),
          legend.position="top")
  
  return(v_plot)
}

barplotFromAncombc <- function(qza_path,
                               log2fc_col,
                               pval_col,
                               pcutoff = 0.05,
                               logfccutoff = 2,
                               taxonomy_df = taxonomy,
                               ylim = NULL,
                               ylab = NULL,
                               down_color,
                               up_color,
                               down_legend,
                               up_legend) 
{
  
  ancombc_data <- import_ancombc(qza_path) %>%
    dplyr::left_join(taxonomy_df) %>% 
    dplyr::mutate(
      tax_label = dplyr::coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom),
      tax_label_clean = gsub("_", " ", tax_label),
      id_prefix = substr(id, 1, 6),
      feature_label_sort = paste0(tax_label_clean, " (", id_prefix, ")")
    )
  
  plot_data <- ancombc_data %>%
    dplyr::filter(.data[[pval_col]] < pcutoff)
  
  if (nrow(plot_data) == 0) {
    message("No features found with ", pval_col, " < ", pcutoff)
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = "No significant features found") +
        ggplot2::theme_classic()
    )
  }
  
  y_axis_label <- if (is.null(ylab)) {
    paste("Log2 Fold Change (", log2fc_col, ")")
  } else {
    ylab
  }
  
  b_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = reorder(feature_label_sort, .data[[log2fc_col]]),
      y = .data[[log2fc_col]]
    )
  ) +

    ggplot2::geom_col(ggplot2::aes(fill = Order), show.legend = TRUE) +
    ggplot2::scale_fill_manual(values = barplot_order_ITS_colors) +
    ggplot2::geom_hline(yintercept = logfccutoff, linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = -logfccutoff, linetype = "dashed", color = "grey50") +
    
    ggplot2::labs(
      y = y_axis_label,
      title = NULL
    ) +
    
    ggplot2::theme_classic() +
    
    ggplot2::theme(
      legend.position = "top",
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.x = element_blank()
    )
  
  if (nrow(plot_data) > 0) {
    n_neg <- sum(plot_data[[log2fc_col]] < 0)
    n_pos <- sum(plot_data[[log2fc_col]] > 0)
    
    if (n_neg > 0) {
      b_plot <- b_plot +
        ggplot2::annotate(
          "text",
          x = (1 + n_neg) / 2,
          y = 0.2,
          label = down_legend,
          color = down_color,
          size = 4,
          hjust = 0.5,
          vjust = -0.5
        )
    }
    
    if (n_pos > 0) {
      b_plot <- b_plot +
        ggplot2::annotate(
          "text",
          x = ((n_neg + 1) + (n_neg + n_pos)) / 2,
          y = -0.2,
          label = up_legend,
          color = up_color,
          size = 4,
          hjust = 0.5,
          vjust = 1.5
        )
    }
  }
  
  if (!is.null(ylim)) {
    b_plot <- b_plot + ggplot2::coord_cartesian(ylim = ylim)
  }
  
  return(b_plot)
}

## Import QIIME 2 files

project_name <- "grano_ITS"
out <- "grano-ITS-TEF1"
bar_tag <- "HordeumVulgare"
oat_tag <- "AvenaSativa"
whe_tag <- "Triticum"

readRenviron("/home/sergio/Renvs/.RenvBrigit")
brigit_IP <- Sys.getenv("IP_ADDRESS")
cluster_path <- paste0("/run/user/1001/gvfs/sftp:host=",
                       brigit_IP,
                       ",user=salias/mnt/lustre")
project_dir <- file.path(cluster_path,
                         "scratch/salias/projects",
                         project_name)
outdir <- file.path("/home/sergio/scratch",
                    out,
                    "abundance")

cereal_wheat_file_path <- file.path(project_dir,
                                    paste0("qiime2/abundance/Cereal_",
                                           whe_tag,
                                           "/filtered_ancombc.qza"))

cereal_barley_file_path <- file.path(project_dir,
                                     paste0("qiime2/abundance/Cereal_",
                                            bar_tag,
                                            "/filtered_ancombc.qza"))

taxonomy_file_path <- file.path(project_dir,
                                "qiime2/taxonomy/taxonomy.qza")

taxonomy <- read_qza(taxonomy_file_path)$data

taxonomy %<>% parse_taxonomy() %>% rownames_to_column("id")


## Volcano plots

### Barley vs wheat

v_bar_vs_whe <- volcanoFromAncombc(qza_path = cereal_wheat_file_path,
                                    lab_col = "labnames",
                                    log2fc_col = paste0("Cereal", bar_tag, "_lfc"),
                                    pval_col = paste0("Cereal", bar_tag, "_q_val"),
                                    up_color = cereal_colors[["Barley"]],
                                    down_color = cereal_colors[["Wheat"]],
                                    up_shape = cereal_shapes[["Barley"]],
                                    down_shape = cereal_shapes[["Wheat"]],
                                    up_legend = "DA (Barley)",
                                    down_legend = "DA (Wheat)",
                                    colAlpha = 1,
                                    cutoffLineCol = "grey70",
                                    ylab = bquote(~Log[10]~ "Q-value"),
                                    title = NULL,
                                    subtitle = NULL,
                                    caption = NULL,
                                    xlim = c(-8, 8),
                                    ylim = c(0, 60),
                                    drawConnectors = TRUE,
                                    typeConnectors = "open",
                                    widthConnectors = 0,
                                    colConnectors = "grey40",
                                    pointSize = 4.0,
                                    labSize = 3,
                                    boxedLabels = FALSE,
                                    parseLabels = TRUE)

pdf(file.path(outdir, "volcano_bar_vs_whe.pdf"))

v_bar_vs_whe

dev.off()

### Oat vs wheat

v_oat_vs_whe <- volcanoFromAncombc(qza_path = cereal_wheat_file_path,
                                   lab_col = "labnames",
                                   log2fc_col = paste0("Cereal", oat_tag, "_lfc"),
                                   pval_col = paste0("Cereal", oat_tag, "_q_val"),
                                   up_color = cereal_colors[["Oat"]],
                                   down_color = cereal_colors[["Wheat"]],
                                   up_shape = cereal_shapes[["Oat"]],
                                   down_shape = cereal_shapes[["Wheat"]],
                                   up_legend = "DA (Oat)",
                                   down_legend = "DA (Wheat)",
                                   colAlpha = 1,
                                   cutoffLineCol = "grey70",
                                   ylab = bquote(~Log[10]~ "Q-value"),
                                   title = NULL,
                                   subtitle = NULL,
                                   caption = NULL,
                                   xlim = c(-8, 8),
                                   ylim = c(0, 60),
                                   drawConnectors = TRUE,
                                   typeConnectors = "open",
                                   widthConnectors = 0,
                                   colConnectors = "grey40",
                                   pointSize = 4.0,
                                   labSize = 3,
                                   boxedLabels = FALSE,
                                   parseLabels = TRUE)

pdf(file.path(outdir, "volcano_oat_vs_whe.pdf"))

v_oat_vs_whe

dev.off()

### Oat vs barley

v_oat_vs_bar <- volcanoFromAncombc(qza_path = cereal_barley_file_path,
                                   lab_col = "labnames",
                                   log2fc_col = paste0("Cereal", oat_tag, "_lfc"),
                                   pval_col = paste0("Cereal", oat_tag, "_q_val"),
                                   up_color = cereal_colors[["Oat"]],
                                   down_color = cereal_colors[["Barley"]],
                                   up_shape = cereal_shapes[["Oat"]],
                                   down_shape = cereal_shapes[["Barley"]],
                                   up_legend = "DA (Oat)",
                                   down_legend = "DA (Barley)",
                                   colAlpha = 1,
                                   cutoffLineCol = "grey70",
                                   ylab = bquote(~Log[10]~ "Q-value"),
                                   title = NULL,
                                   subtitle = NULL,
                                   caption = NULL,
                                   xlim = c(-8, 8),
                                   ylim = c(0, 60),
                                   drawConnectors = TRUE,
                                   typeConnectors = "open",
                                   widthConnectors = 0,
                                   pointSize = 4.0,
                                   labSize = 3,
                                   boxedLabels = FALSE,
                                   parseLabels = TRUE)

pdf(file.path(outdir, "volcano_oat_vs_bar.pdf"))

v_oat_vs_bar

dev.off()

## Volcano legends

legend_data_v <- data.frame(
  Crop = factor(names(da_legend_colors), levels = names(da_legend_colors))
)

legend_plot_v <- ggplot(legend_data_v, aes(x = 1, y = 1, color = Crop, shape = Crop)) +
  geom_point(size = 4) + 
  scale_color_manual(
    name = NULL,
    values = da_legend_colors
  ) +
  scale_shape_manual(
    name = NULL,
    values = da_legend_shapes
  ) +
  theme_void() +
  theme(legend.direction = "horizontal")

v_legend <- plot_grid(cowplot::get_legend(legend_plot_v))

## Barplots

b_bar_vs_whe <- barplotFromAncombc(
  qza_path = cereal_wheat_file_path,
  log2fc_col = paste0("Cereal", bar_tag, "_lfc"),
  pval_col = paste0("Cereal", bar_tag, "_q_val"),
  taxonomy_df = taxonomy,
  pcutoff = 0.05,
  logfccutoff = 2,
  ylim = c(-5.5, 5.5), 
  ylab = bquote(~Log[2]~ "fold change"),
  down_color = cereal_colors[["Wheat"]],
  up_color = cereal_colors[["Barley"]],
  down_legend = "Wheat",
  up_legend = "Barley")


b_bar_vs_whe

b_oat_vs_whe <- barplotFromAncombc(
  qza_path = cereal_wheat_file_path,
  log2fc_col = paste0("Cereal", oat_tag, "_lfc"),
  pval_col = paste0("Cereal", oat_tag, "_q_val"),
  taxonomy_df = taxonomy,
  pcutoff = 0.05,
  logfccutoff = 2,
  ylim = c(-5.5, 5.5), 
  ylab = bquote(~Log[2]~ "fold change"),
  down_color = cereal_colors[["Wheat"]],
  up_color = cereal_colors[["Oat"]],
  down_legend = "Wheat",
  up_legend = "Oat")


b_oat_vs_whe

b_oat_vs_bar <- barplotFromAncombc(
  qza_path = cereal_barley_file_path,
  log2fc_col = paste0("Cereal", oat_tag, "_lfc"),
  pval_col = paste0("Cereal", oat_tag, "_q_val"),
  taxonomy_df = taxonomy,
  pcutoff = 0.05,
  logfccutoff = 2,
  ylim = c(-5.5, 5.5), 
  ylab = bquote(~Log[2]~ "fold change"),
  down_color = cereal_colors[["Barley"]],
  up_color = cereal_colors[["Oat"]],
  down_legend = "Barley",
  up_legend = "Oat")


b_oat_vs_bar

### Barplot legend

legend_data_b <- data.frame(
  Order = names(barplot_order_ITS_colors)
)

legend_plot_b <- ggplot(legend_data_b, aes(x = 1, fill = Order)) +
  geom_bar(alpha = 1) + 
  scale_fill_manual(
    name = "Order",
    values = barplot_order_ITS_colors,
    breaks = names(barplot_order_ITS_colors)
  ) +
  theme_void() +
  guides(fill = guide_legend(nrow = 2,
                             title.position = "left"))

b_legend <- plot_grid(cowplot::get_legend(legend_plot_b))

## Grouped plots

v_vol <- (v_bar_vs_whe + theme(legend.position="none") +
          plot_spacer() +
          v_oat_vs_whe + theme(legend.position="none") +
          plot_spacer() +
          v_oat_vs_bar + theme(legend.position="none") +
          theme(plot.tag.position = "topleft")) +
          plot_layout(axis_titles = "collect",
                      axes = "collect",
                      widths = c(9, 1, 9, 1, 9))

p_v <- v_legend / v_vol +
  plot_layout(heights = c(1, 15))

b_bars <- (b_bar_vs_whe + theme(legend.position="none") +
             plot_spacer() +
        b_oat_vs_whe + theme(legend.position="none") +
          plot_spacer() +
        b_oat_vs_bar + theme(legend.position="none")) +
        plot_layout(axis_titles = "collect",
                    axes = "collect",
                    widths = c(9, 1, 9, 1, 9))

p_b <- b_bars / b_legend / plot_spacer() + 
  plot_layout(heights = c(7, 1, 0.5))

v_b <- ggarrange(p_v, p_b, 
                 ncol = 1, 
                 nrow = 2, 
                 heights = c(4, 3))

## Save plots

pdf(file.path(outdir, "da.pdf"),
    width = 12,
    height = 8)

v_b

dev.off()

png(file.path(outdir, "da.png"),
    width = 12,
    height = 8,
    units = "in",
    res = 300)

v_b

dev.off()
