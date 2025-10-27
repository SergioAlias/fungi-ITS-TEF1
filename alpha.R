# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            alpha.R                                ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : fungi-ITS-TEF1                                   ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-09-23                                       ║
# ║ Last Modified  : 2025-10-27                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(tidyverse)
library(qiime2R)
library(ggpubr)
library(dplyr)
library(patchwork)
library(rstatix)

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
                    "alpha")

metadata <- read.csv(file.path("/home/sergio/scratch",
                               local_metadata,
                               "metadata.tsv"),
                     sep = "\t")

colnames(metadata)[1] <- "SampleID"

shannon_file_path <- file.path(project_dir,
                               "qiime2/diversity/shannon_vector.qza")
simpson_file_path <- file.path(project_dir,
                               "qiime2/diversity/simpson_vector.qza")
chao1_file_path <- file.path(project_dir,
                             "qiime2/diversity/chao1_vector.qza")

shannon <- read_qza(shannon_file_path)
shannon <- shannon$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(shannon)

simpson <- read_qza(simpson_file_path)
simpson <- simpson$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(simpson)

chao1 <- read_qza(chao1_file_path)
chao1 <- chao1$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(chao1)

metadata %<>% mutate(
  Cereal = case_when(
    Cereal == "AvenaSativa" ~ "Oat",
    Cereal == "HordeumVulgare" ~ "Barley",
    TRUE ~ "Wheat"
  )
)
metadata <- metadata[order(metadata$Cereal, decreasing = TRUE), ]

## Colors and shapes

source("/home/sergio/projects/fungi-ITS-TEF1/colors.R")

## Functions

plot_alpha <- function(alpha_metric,
                       ylab,
                       letters_df) {
  
  p <- metadata %>%
    ggboxplot(
      "Cereal", alpha_metric,
      color = "Cereal",
      fill = "Cereal",
      alpha = 0.1,
      palette = cereal_colors,
      add = "jitter",
      add.params = list(alpha = 0.6),
      shape = "Cereal"
    ) +
    scale_shape_manual(values = cereal_shapes, name = "Cereal") +
    ylab(ylab) +
    geom_text(
      data = letters_df,
      aes(x = Cereal, y = Position, label = Letter),
      inherit.aes = FALSE,
      color = "grey40",
    ) +
    coord_flip() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 11)
    )
  return(p)
}

### Shannon

shannon_letters <- data.frame(
  Cereal = c("Oat", "Barley", "Wheat"),
  Letter = c("a", "ab", "b"),
  Position = rep(0.4, 3)
)

pdf(file.path(outdir, "shannon_cereal.pdf"))

shannon_c <- plot_alpha(alpha_metric = "shannon_entropy",
                        ylab = "Shannon",
                        letters_df = shannon_letters)

shannon_c

dev.off()


### Simpson

simpson_letters <- data.frame(
  Cereal = c("Oat", "Barley", "Wheat"),
  Letter = c("a", "a", "b"),
  Position = rep(0.1, 3)
)

pdf(file.path(outdir, "simpson_cereal.pdf"))

simpson_c <- plot_alpha(alpha_metric = "simpson",
                        ylab = "Inverse Simpson",
                        letters_df = simpson_letters)

simpson_c

dev.off()

### Chao1

chao1_letters <- data.frame(
  Cereal = c("Oat", "Barley", "Wheat"),
  Letter = c("b", "a", "ab"),
  Position = rep(6, 3)
)

pdf(file.path(outdir, "chao1_cereal.pdf"))

chao1_c <- plot_alpha(alpha_metric = "chao1",
                        ylab = "Chao1",
                        letters_df = chao1_letters)

chao1_c

dev.off()

### Grouped plots

pdf(file.path(outdir, "patched_alhpa.pdf"),
    width = 15,
    height = 6.5)

(chao1_c + theme(legend.position="none") +
    shannon_c + theme(legend.position="none") +
    simpson_c + theme(legend.position="none") +
    theme(plot.tag.position = "topleft")) +
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = 'A')

dev.off()

### Save RDS

saveRDS(chao1_c, file = file.path(outdir, "chao1.RDS"))
saveRDS(shannon_c, file = file.path(outdir, "shannon.RDS"))
saveRDS(simpson_c, file = file.path(outdir, "simpson.RDS"))
