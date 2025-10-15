# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            pcoa.R                                 ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : tenebrio-16S                                     ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-15                                       ║
# ║ Last Modified  : 2025-10-15                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(rlang)
library(magrittr, include.only = "%<>%")
library(tidyverse)
library(qiime2R)
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
                    "beta")

jaccard_file_path <- file.path(project_dir,
                               "qiime2/diversity/jaccard_pcoa_results.qza")
bray_curtis_file_path <- file.path(project_dir,
                                   "qiime2/diversity/bray_curtis_pcoa_results.qza")
aitchison_file_path <- file.path(project_dir,
                                 "qiime2/diversity/aitchison_pcoa_results.qza")

jaccard <- read_qza(jaccard_file_path)
bray_curtis <- read_qza(bray_curtis_file_path)
aitchison <- read_qza(aitchison_file_path)

metadata <- read.csv(file.path("/home/sergio/scratch",
                               local_metadata,
                               "metadata.tsv"),
                     sep = "\t")

colnames(metadata)[1] <- "SampleID"

metadata %<>% mutate(
  Cereal = case_when(
    Cereal == "AvenaSativa" ~ "Avena sativa",
    Cereal == "HordeumVulgare" ~ "Hordeum vulgare",
    TRUE ~ "Triticum sp."
  )
)

## Get variance explained by each PCo 

jaccard_pco1 <- round(jaccard[["data"]]$ProportionExplained$PC1 * 100, 2)
jaccard_pco2 <- round(jaccard[["data"]]$ProportionExplained$PC2 * 100, 2)

bray_curtis_pco1 <- round(bray_curtis[["data"]]$ProportionExplained$PC1 * 100, 2)
bray_curtis_pco2 <- round(bray_curtis[["data"]]$ProportionExplained$PC2 * 100, 2)

aitchison_pco1 <- round(aitchison[["data"]]$ProportionExplained$PC1 * 100, 2)
aitchison_pco2 <- round(aitchison[["data"]]$ProportionExplained$PC2 * 100, 2)


## Colors and shapes

source("/home/sergio/projects/fungi-ITS-TEF1/colors.R")

## PCoA plots

### Jaccard

p_j <- jaccard$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
    ggplot(aes(x = PC1, y = PC2, color = `Cereal`, shape = `Cereal`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = cereal_colors, name = "Cereal") +
  scale_shape_manual(values = cereal_shapes, name = "Cereal") +
  xlab(paste0("PCo-1 | ", jaccard_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", jaccard_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Cereal`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1))) +
  theme(legend.text = element_text(face = "italic"))

pdf(file.path(outdir, "pcoa_jaccard.pdf"),
    height = 5,
    width = 7)

p_j

dev.off()

### Bray-Curtis

p_b <- bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Cereal`, shape = `Cereal`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = cereal_colors, name = "Cereal") +
  scale_shape_manual(values = cereal_shapes, name = "Cereal") +
  xlab(paste0("PCo-1 | ", bray_curtis_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", bray_curtis_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Cereal`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1))) +
  theme(legend.text = element_text(face = "italic"))

pdf(file.path(outdir, "pcoa_bray_curtis.pdf"),
    height = 5,
    width = 7)

p_b

dev.off()

### Aitchison

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
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1))) +
  theme(legend.text = element_text(face = "italic"))

pdf(file.path(outdir, "pcoa_aitchison.pdf"),
    height = 5,
    width = 7)

p_a

dev.off()
