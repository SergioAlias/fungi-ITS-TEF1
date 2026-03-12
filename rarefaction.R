# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            rarefaction.R                          ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : fungi-ITS-TEF1                                   ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2026-03-12                                       ║
# ║ Last Modified  : 2026-03-12                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(tidyverse)
library(magrittr, include.only = "%<>%")
library(patchwork)

## Colors and shapes

source("/home/sergio/projects/fungi-ITS-TEF1/colors.R")

## Functions

plot_rarefaction <- function(csv_file, metadata, y_title, threshold = NULL) {
  data <- read_csv(csv_file, show_col_types = FALSE)
  
  long_data <- data %>%
    pivot_longer(
      cols = -`sample-id`, 
      names_to = c("depth", "iteration"),
      names_pattern = "depth-(.*)_iter-(.*)",
      values_to = "observed_features",
      values_drop_na = TRUE
    ) %>%
    mutate(
      depth = as.numeric(depth),
      iteration = as.numeric(iteration)
    )
  
  summary_data <- long_data %>%
    group_by(`sample-id`, depth) %>%
    summarise(
      mean_features = mean(observed_features),
      .groups = 'drop'
    )
  
  summary_data <- summary_data %>%
    left_join(metadata, by = c("sample-id" = "SampleID"))
  
  p <- ggplot(summary_data, aes(x = depth, y = mean_features, group = `sample-id`, color = Cereal)) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    scale_color_manual(values = cereal_colors) + 
    labs(
      x = "Sequencing depth",
      y = y_title,            
      color = "Cereal"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.5, "cm"),
      plot.margin = margin(t = 20, r = 5, b = 5, l = 5) 
    )
  
  if (!is.null(threshold)) {
    p <- p + geom_vline(
      xintercept = threshold, 
      linetype = "dotted", 
      color = "black", 
      linewidth = 0.8,
      alpha = 0.8
    ) +
      annotate(
        "text", 
        x = threshold, 
        y = Inf,                         
        label = as.character(threshold), 
        angle = 0,                        
        vjust = -1,                      
        hjust = 0.5,
        size = 3.5                        
      ) +
      coord_cartesian(clip = "off")       
  }
  return(p)
}

clean_metadata <- function(path) {
  df <- read.csv(path, sep = "\t")
  df %<>%
    rename(SampleID = 1) %>% 
    mutate(
      Cereal = case_when(
        Cereal == "AvenaSativa" ~ "Oat",
        Cereal == "HordeumVulgare" ~ "Barley",
        TRUE ~ "Wheat"
      )
    )
  return(df)
}

## Import files

out <- "grano-ITS-TEF1"
outdir <- file.path("/home/sergio/scratch",
                    out,
                    "raref")
metadata_ITS <- clean_metadata("/home/sergio/scratch/grano_ITS/metadata.tsv")
metadata_TEF1 <- clean_metadata("/home/sergio/scratch/grano_TEF1/metadata.tsv")


## Rarefaction plots

of_its <- plot_rarefaction("/home/sergio/projects/fungi-ITS-TEF1/local/raref/observed_features_ITS.csv",
                       metadata_ITS,
                       y_title = "Observed features",
                       threshold = 7000)

s_its <- plot_rarefaction("/home/sergio/projects/fungi-ITS-TEF1/local/raref/shannon_ITS.csv",
                          metadata_ITS,
                          y_title = "Shannon",
                          threshold = 7000)

of_tef1 <- plot_rarefaction("/home/sergio/projects/fungi-ITS-TEF1/local/raref/observed_features_TEF1.csv",
                            metadata_TEF1,
                            y_title = "Observed features",
                            threshold = 500)

s_tef1 <- plot_rarefaction("/home/sergio/projects/fungi-ITS-TEF1/local/raref/shannon_TEF1.csv",
                           metadata_TEF1,
                           y_title = "Shannon",
                           threshold = 500)


title_its <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "ITS2", fontface = "bold", size = 5) + 
  theme_void()

title_tef1 <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "TEF1", fontface = "bold", size = 5) + 
  theme_void()

p_raref <- title_its / 
  (of_its + s_its) / 
  plot_spacer() / 
  title_tef1 / 
  (of_tef1 + s_tef1) + 
  plot_layout(
    guides = "collect",
    heights = c(0.1, 1, 0.05, 0.1, 1) 
  ) & 
  theme(
    legend.position = "right"
  )


## Save plots

pdf(file.path(outdir, "rarefaction_curves.pdf"),
    width = 9,
    height = 7)

p_raref

dev.off()

png(file.path(outdir, "rarefaction_curves.png"),
    width = 9,
    height = 7,
    units = "in",
    res = 300)

p_raref

dev.off()
