# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            alpha.R                                ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : fungi-ITS-TEF1                                   ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-09-23                                       ║
# ║ Last Modified  : 2025-09-23                                       ║
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
    Cereal == "AvenaSativa" ~ "Avena sativa",
    Cereal == "HordeumVulgare" ~ "Hordeum vulgare",
    TRUE ~ "Triticum sp."
  )
)

## Colors and shapes

source("/home/sergio/projects/fungi-ITS-TEF1/colors.R")


## Comparisons

comparisons_cereal <- combn(unique(metadata$Cereal), 2, simplify = FALSE)
comparisons_location <- combn(unique(metadata$Location), 2, simplify = FALSE)


## Alpha boxplots

if (amplicon == "ITS"){
  shannon_pos_stat <- 3.6
  shannon_pos_stat_paired <- c(3.1, 3.35, 3.2)
  simpson_pos_stat <- 0.92
  simpson_pos_stat_paired <- c(0.78, 0.85, 0.81)
  chao1_pos_stat <- 72
  chao1_pos_stat_paired <- c(61, 66.6, 63.1)
}
if (amplicon == "TEF1"){
  shannon_pos_stat <- 11.2
  shannon_pos_stat_paired <- c(10.55, 10.7, 10.625)
  simpson_pos_stat <- 0.99885
  simpson_pos_stat_paired <- c(0.99852, 0.99866, 0.99859)
  chao1_pos_stat <- 5500
  chao1_pos_stat_paired <- c(3100, 3500, 3300)
}

### Shannon

pdf(file.path(outdir, "shannon_cereal.pdf"))

shannon_c <- metadata %>%
  ggboxplot("Cereal", "shannon_entropy",
            color = "Cereal",
            fill = "Cereal",
            alpha = 0.1,
            palette = cereal_colors,
            add = "jitter",
            shape = "Cereal") +
  scale_shape_manual(values = cereal_shapes, name = "Cereal") +
  ylab("Shannon") +
  stat_compare_means(label.y = shannon_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons_cereal,
                     label.y = shannon_pos_stat_paired) +
  theme(axis.text.x = element_text(face = "italic"))

shannon_c

dev.off()


### Simpson

pdf(file.path(outdir, "simpson_cereal.pdf"))

simpson_c <- metadata %>%
  ggboxplot("Cereal", "simpson",
            color = "Cereal",
            fill = "Cereal",
            alpha = 0.1,
            palette = cereal_colors,
            add = "jitter",
            shape = "Cereal") +
  scale_shape_manual(values = cereal_shapes, name = "Cereal") +
  ylab("Inverse Simpson") +
  stat_compare_means(label.y = simpson_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons_cereal,
                     label.y = simpson_pos_stat_paired) +
  theme(axis.text.x = element_text(face = "italic"))

simpson_c

dev.off()

### Chao1

pdf(file.path(outdir, "chao1_cereal.pdf"))

chao1_c <- metadata %>%
  ggboxplot("Cereal", "chao1",
            color = "Cereal",
            fill = "Cereal",
            alpha = 0.1,
            palette = cereal_colors,
            add = "jitter",
            shape = "Cereal") +
  scale_shape_manual(values = cereal_shapes, name = "Cereal") +
  ylab("Chao1") +
  stat_compare_means(label.y = chao1_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons_cereal,
                     label.y = chao1_pos_stat_paired) +
  theme(axis.text.x = element_text(face = "italic"))

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
