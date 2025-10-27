# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            rpca.R                                 ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : tenebrio-16S                                     ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-15                                       ║
# ║ Last Modified  : 2025-10-27                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(rlang)
library(magrittr, include.only = "%<>%")
library(tidyverse)
library(qiime2R)
library(patchwork)
library(ggrepel)
library(tibble)
library(ggnewscale)
library(viridisLite)

## Import QIIME 2 files

project_name <- "grano_ITS" # grano_ITS or grano_TEF1
amplicon <- "ITS" # ITS or TEF1
local_metadata <- project_name
out <- "grano-ITS-TEF1"

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
                    "beta")

jaccard_file_path <- file.path(project_dir,
                               "qiime2/diversity/jaccard_pcoa_results.qza")
bray_curtis_file_path <- file.path(project_dir,
                                   "qiime2/diversity/bray_curtis_pcoa_results.qza")
aitchison_file_path <- file.path(project_dir,
                                 "qiime2/diversity/aitchison_pcoa_results.qza")
gemelli_file_path <- file.path(project_dir,
                                 "qiime2/diversity/gemelli_rpca_results.qza")
taxonomy_file_path <- file.path(project_dir,
                                    "qiime2/taxonomy/taxonomy.qza")

jaccard <- read_qza(jaccard_file_path)
bray_curtis <- read_qza(bray_curtis_file_path)
aitchison <- read_qza(aitchison_file_path)
gemelli <- read_qza(gemelli_file_path)
taxonomy <- read_qza(taxonomy_file_path)$data %>% parse_taxonomy()
taxonomy %<>%
  tibble::rownames_to_column(var = "FeatureID") %>%
  mutate(
    feature_short = substr(FeatureID, 1, 6),
    taxa_level = case_when(
      !is.na(Species) ~ "Species",
      !is.na(Genus)   ~ "Genus",
      TRUE            ~ "Other"
    ),
    taxa_name = coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom, "Unknown"),
    taxa_name = gsub("_", " ", taxa_name),
    arrownames = if_else(
      taxa_level %in% c("Species", "Genus"),
      paste0("italic('", taxa_name, "') ~ '(", feature_short, ")'"),
      paste0("'", taxa_name, "' ~ '(", feature_short, ")'")
    )
  ) %>%
  select(-feature_short, -taxa_level, -taxa_name) %>% 
  tibble::column_to_rownames(var = "FeatureID")

metadata <- read.csv(file.path("/home/sergio/scratch",
                               local_metadata,
                               "metadata.tsv"),
                     sep = "\t")

colnames(metadata)[1] <- "SampleID"

metadata %<>% mutate(
  Cereal = case_when(
    Cereal == "AvenaSativa" ~ "Oat",
    Cereal == "HordeumVulgare" ~ "Barley",
    TRUE ~ "Wheat"
  )
)

## Get variance explained by each PCo 

jaccard_pco1 <- round(jaccard[["data"]]$ProportionExplained$PC1 * 100, 2)
jaccard_pco2 <- round(jaccard[["data"]]$ProportionExplained$PC2 * 100, 2)
jaccard_pco3 <- round(jaccard[["data"]]$ProportionExplained$PC3 * 100, 2)

bray_curtis_pco1 <- round(bray_curtis[["data"]]$ProportionExplained$PC1 * 100, 2)
bray_curtis_pco2 <- round(bray_curtis[["data"]]$ProportionExplained$PC2 * 100, 2)
bray_curtis_pco3 <- round(bray_curtis[["data"]]$ProportionExplained$PC3 * 100, 2)

aitchison_pco1 <- round(aitchison[["data"]]$ProportionExplained$PC1 * 100, 2)
aitchison_pco2 <- round(aitchison[["data"]]$ProportionExplained$PC2 * 100, 2)
aitchison_pco3 <- round(aitchison[["data"]]$ProportionExplained$PC3 * 100, 2)

gemelli_pco1 <- round(gemelli[["data"]]$ProportionExplained$PC1 * 100, 2)
gemelli_pco2 <- round(gemelli[["data"]]$ProportionExplained$PC2 * 100, 2)
gemelli_pco3 <- round(gemelli[["data"]]$ProportionExplained$PC3 * 100, 2)

## Colors and shapes

source("/home/sergio/projects/fungi-ITS-TEF1/colors.R")

## Functions

plot_pcoa_rpca <- function(pcoa_vectors = bray_curtis$data$Vectors, 
                           metadata_df = metadata, 
                           pc_x, 
                           pc_y, 
                           variance_x, 
                           variance_y,
                           grouping_var = "Cereal",
                           color_values = cereal_colors,
                           shape_values = cereal_shapes,
                           biplot = FALSE,
                           species_vectors = NULL,
                           prefix = "PCo-",
                           top_n_features = NULL,
                           taxonomy_df = NULL,
                           cleaner_plot = FALSE) {

  pc_x_sym <- sym(pc_x)
  pc_y_sym <- sym(pc_y)
  grouping_sym <- sym(grouping_var)
  
  cent_x_name <- paste0(pc_x, "_cent")
  cent_y_name <- paste0(pc_y, "_cent")
  cent_x_sym <- sym(cent_x_name)
  cent_y_sym <- sym(cent_y_name)
  
  plot_data <- pcoa_vectors %>%
    select(SampleID, !!pc_x_sym, !!pc_y_sym) %>%
    left_join(metadata_df, by = "SampleID")
  
  centroids <- plot_data %>%
    group_by(!!grouping_sym) %>%
    summarise(
      !!cent_x_name := mean(!!pc_x_sym),
      !!cent_y_name := mean(!!pc_y_sym),
      .groups = 'drop'
    )
  
  hull_data <- plot_data %>%
    group_by(!!grouping_sym) %>%
    slice(chull(!!pc_x_sym, !!pc_y_sym))
  
  plot_data_with_centroids <- plot_data %>%
    left_join(centroids, by = grouping_var)
  
  xlab_text <- paste0(prefix, sub("PC", "", pc_x), " (", variance_x, "%)")
  ylab_text <- paste0(prefix, sub("PC", "", pc_y), " (", variance_y, "%)")
  
  point_alpha_val <- if (cleaner_plot) 0.3 else 1.0
  
  p <- ggplot(plot_data_with_centroids,
              aes(x = !!pc_x_sym, y = !!pc_y_sym, color = !!grouping_sym, shape = !!grouping_sym))
  
  if (cleaner_plot) {
    p <- p + geom_polygon(data = hull_data, 
                          aes(fill = !!grouping_sym), 
                          alpha = 0.03, 
                          show.legend = FALSE,
                          linetype = "blank") 
  } else {
    p <- p + geom_polygon(data = hull_data, 
                          aes(fill = !!grouping_sym), 
                          alpha = 0.1, 
                          show.legend = FALSE,
                          linetype = "dotted", 
                          linewidth = 0.5) +
      geom_segment(aes(xend = !!cent_x_sym, yend = !!cent_y_sym), 
                   alpha = 0.2, 
                   show.legend = FALSE) +
      geom_point(data = centroids, 
                 aes(x = !!cent_x_sym, y = !!cent_y_sym),
                 size = 4, 
                 shape = 1,
                 show.legend = FALSE)
  }

  p <- p + 
    geom_point(alpha = point_alpha_val, size = 2) + 
    theme_classic() +
    scale_color_manual(values = color_values, name = NULL) +
    scale_shape_manual(values = shape_values, name = NULL) +
    scale_fill_manual(values = color_values, name = grouping_var) + 
    xlab(xlab_text) +
    ylab(ylab_text)
  
  if (cleaner_plot) {
    p <- p + guides(
      color = "none",
      shape = "none",
      fill = "none"
    )
  } else {
    p <- p + guides(
      color = guide_legend(override.aes = list(linetype = 0, alpha = 1, size = 4), order = 1),
      shape = guide_legend(order = 1) 
    )
  }

  if (biplot) {
    if (is.null(species_vectors)) {
      stop("Error: 'species_vectors' must be provided when biplot = TRUE.")
    }

    arrow_data_full <- species_vectors %>%
      select(FeatureID, !!pc_x_sym, !!pc_y_sym)
    
    max_sample_abs <- max(abs(c(plot_data_with_centroids[[pc_x]], 
                                plot_data_with_centroids[[pc_y]])), 
                          na.rm = TRUE)
    max_feature_abs_full <- max(abs(c(arrow_data_full[[pc_x]], 
                                      arrow_data_full[[pc_y]])), 
                                na.rm = TRUE)
    
    if (max_feature_abs_full > 0) {
      scaling_factor <- max_sample_abs / max_feature_abs_full
    } else {
      scaling_factor <- 1 
    }
    if (!is.null(top_n_features)) {
      if (!is.numeric(top_n_features) || top_n_features <= 0) {
        stop("Error: 'top_n_features' must be a positive number.")
      }
      
      arrow_data <- arrow_data_full %>%
        mutate(magnitude_sq = (!!pc_x_sym)^2 + (!!pc_y_sym)^2) %>%
        arrange(desc(magnitude_sq)) %>%
        slice_head(n = as.integer(top_n_features))
      
    } else {
      arrow_data <- arrow_data_full
    }

    if (!is.null(taxonomy_df)) {
      
      if (!"arrownames" %in% colnames(taxonomy_df)) {
        stop("Error: 'taxonomy_df' must contain a column named 'arrownames'.")
      }
      
      tax_data_to_join <- taxonomy_df %>%
        tibble::rownames_to_column(var = "FeatureID") %>%
        select(FeatureID, arrownames)
      
      arrow_data <- arrow_data %>%
        left_join(tax_data_to_join, by = "FeatureID") %>%
        rename(label_col = arrownames) 
      
    } else {
      arrow_data <- arrow_data %>%
        mutate(label_col = FeatureID)
    }

    arrow_data <- arrow_data %>%
      mutate(
        !!pc_x_sym := !!pc_x_sym * scaling_factor,
        !!pc_y_sym := !!pc_y_sym * scaling_factor
      )
    p <- p + 
      ggnewscale::new_scale_color() +
      geom_segment(data = arrow_data,
                   aes(x = 0, y = 0, xend = !!pc_x_sym, yend = !!pc_y_sym, color = label_col),
                   inherit.aes = FALSE, 
                   linewidth = 0.7, 
                   arrow = arrow(length = unit(0.2, "cm"))) +
      scale_color_viridis_d(
        name = NULL, 
        option = "turbo",  
        labels = function(l) parse(text = l) 
      ) +
      guides(color = guide_legend(override.aes = list(alpha = 1), order = 2))
  }
  return(p)
}

## PCoA plots

### Bray-Curtis (hull and centroids)



p_b <- plot_pcoa_hulls(pc_x = "PC1",
                       pc_y = "PC2",
                       variance_x = bray_curtis_pco1,
                       variance_y = bray_curtis_pco2)

p_b13 <- plot_pcoa_hulls(pc_x = "PC1",
                         pc_y = "PC3",
                         variance_x = bray_curtis_pco1,
                         variance_y = bray_curtis_pco3)

p_b23 <- plot_pcoa_hulls(pc_x = "PC2",
                         pc_y = "PC3",
                         variance_x = bray_curtis_pco2,
                         variance_y = bray_curtis_pco3)

p_b34 <- plot_pcoa_hulls(pc_x = "PC3",
                         pc_y = "PC4",
                         variance_x = bray_curtis_pco3,
                         variance_y = bray_curtis_pco4)

### Combined Bray-Curtis

p_b + (p_b13 / p_b23) + 
  plot_layout(widths = c(5, 3),
              guides = "collect")

### Aitchison (confidence ellipses)

p_a <- aitchison$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Cereal`, shape = `Cereal`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = cereal_colors, name = "Cereal") +
  scale_shape_manual(values = cereal_shapes, name = "Cereal") +
  xlab(paste0("PCo-1 | ", aitchison_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", aitchison_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Cereal`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1)))

pdf(file.path(outdir, "pcoa_aitchison.pdf"),
    height = 5,
    width = 7)

p_a

dev.off()

### rPCA (robust Aitchison) biplots

rpca_regular <- plot_pcoa_rpca(pcoa_vectors = gemelli$data$Vectors,
                               prefix = "PC",
                               pc_x = "PC1",
                               pc_y = "PC2",
                               variance_x = gemelli_pco1,
                               variance_y = gemelli_pco2)

rpca_biplot <- plot_pcoa_rpca(pcoa_vectors = gemelli$data$Vectors,
                              prefix = "PC",
                              pc_x = "PC1",
                              pc_y = "PC2",
                              variance_x = gemelli_pco1,
                              variance_y = gemelli_pco2,
                              biplot = TRUE,
                              species_vectors = gemelli$data$Species,
                              taxonomy_df = taxonomy,
                              top_n_features = 10,
                              cleaner_plot = TRUE)



rpca_regular + rpca_biplot +
  plot_layout(guides = "collect",
              axes = "collect")

## Save RDS

saveRDS(rpca_regular, file = file.path(outdir, "rPCA_regular.RDS"))
saveRDS(rpca_biplot, file = file.path(outdir, "rPCA_biplot.RDS"))
