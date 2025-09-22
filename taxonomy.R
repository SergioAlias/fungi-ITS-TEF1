# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            taxonomy.R                             ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : fungi-ITS-TEF1                                   ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-09-19                                       ║
# ║ Last Modified  : 2025-09-22                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(file2meco)
library(microeco)
library(ggplot2)
library(ggnested)
library(ggh4x)
library(ggradar)
library(ggtern)


## Import QIIME 2 files

project_ITS <- "grano_ITS"
project_TEF1 <- "grano_TEF1"
color_palette <- "ITS" # 16S or ITS
out <- "grano-ITS-TEF1"

readRenviron("/home/sergio/Renvs/.RenvBrigit")
brigit_IP <- Sys.getenv("IP_ADDRESS")
cluster_path <- paste0("/run/user/1001/gvfs/sftp:host=",
                       brigit_IP,
                       ",user=salias/mnt/lustre")
project_dir_ITS <- file.path(cluster_path,
                             "scratch/salias/projects",
                              project_ITS)
project_dir_TEF1 <- file.path(cluster_path,
                             "scratch/salias/projects",
                             project_TEF1)
outdir <- file.path("/home/sergio/scratch",
                    out,
                    "taxonomy")


dada2_ITS_file_path <- file.path(project_dir_ITS,
                             "qiime2/feature_tables/filtered_table.qza")
dada2_TEF1_file_path <- file.path(project_dir_TEF1,
                                 "qiime2/feature_tables/filtered_table.qza")
metadata_ITS_file_path <- file.path("/home/sergio/scratch",
                                    project_ITS,
                                    "metadata.tsv")
metadata_TEF1_file_path <- file.path("/home/sergio/scratch",
                                    project_TEF1,
                                    "metadata.tsv")
taxonomy_ITS_file_path <- file.path(project_dir_ITS,
                                "qiime2/taxonomy/taxonomy.qza")
taxonomy_TEF1_file_path <- file.path(project_dir_TEF1,
                                    "qiime2/taxonomy/taxonomy.qza")

meco_ITS <- qiime2meco(dada2_ITS_file_path,
                       sample_table = metadata_ITS_file_path,
                       taxonomy_table = taxonomy_ITS_file_path)

meco_TEF1 <- qiime2meco(dada2_TEF1_file_path,
                        sample_table = metadata_TEF1_file_path,
                        taxonomy_table = taxonomy_TEF1_file_path)

## Colors and shapes

source("/home/sergio/projects/diversity-cereal/colors.R")

## Chamge amplicon to plot

meco <- meco_TEF1

## Relabel UNITE prefixes for cleaner plotting

meco$tax_table$Phylum <- gsub("p__", "", meco$tax_table$Phylum)
meco$tax_table$Family <- gsub("f__", "", meco$tax_table$Family)

# Create trans_abund objects and plot stuff

## Nested barplot (Phylum / Class)

t_stacked_phylum <- trans_abund$new(dataset = meco,
                                    taxrank = "Class",
                                    ntaxa = 20,
                                    delete_taxonomy_prefix = TRUE,
                                    high_level = "Phylum",
                                    prefix = "c__")

pdf(file.path(outdir, "barplot_class.pdf"),
    width = 9)

t_stacked_phylum$plot_bar(ggnested = TRUE,
                          high_level_add_other = TRUE,
                          xtext_keep = FALSE,
                          # xtext_angle = 90,
                          # xtext_size = 6,
                          facet = c("Cereal"),
                          others_color = "grey90") + 
  theme(ggh4x.facet.nestline = element_line(colour = "grey95"))

dev.off()

## Nested barplot (Family / Genus)

t_stacked_family <- trans_abund$new(dataset = meco,
                                    taxrank = "Genus",
                                    ntaxa = 15,
                                    delete_taxonomy_prefix = TRUE,
                                    high_level = "Family",
                                    prefix = "g__")

pdf(file.path(outdir, paste0(color_palette, "_barplot_genus.pdf")),
    width = 7.5,
    height = 5.5)

t_stacked_family$plot_bar(ggnested = TRUE,
                          high_level_add_other = TRUE,
                          xtext_keep = FALSE,
                          color_values = get(paste0("barplot_",
                                                    color_palette,
                                                    "_colors")),
                          # xtext_angle = 90,
                          # xtext_size = 6,
                          facet = c("Cereal")) + 
  theme(ggh4x.facet.nestline = element_line(colour = "grey95"))

dev.off()

## Boxplot

t_family <- trans_abund$new(dataset = meco,
                            taxrank = "Family",
                            ntaxa = 20)

pdf(file.path(outdir, "boxplot_family.pdf"))

t_family$plot_box(group = "Cereal",
                  xtext_size = 10,
                  xtext_angle = 30,
                  ytitle_size = 15) + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white", size = 0.5),
    legend.position = c(0.9, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(1, "cm")
  )

dev.off()

## Radar plot 

t_radar <- trans_abund$new(dataset = meco,
                           taxrank = "Class",
                           ntaxa = 8,
                           groupmean = "Cereal")

pdf(file.path(outdir, "radar_plot.pdf"))

t_radar$plot_radar(values.radar = c("0%", "25%", "50%"), grid.min = 0, grid.mid = 0.25, grid.max = 0.5)

dev.off()

## Ternary plot 

t_tern <- trans_abund$new(dataset = meco,
                          taxrank = "Genus",
                          ntaxa = 8,
                          groupmean = "Cereal")

pdf(file.path(outdir, "tern_plot.pdf"))

t_tern$plot_tern()

dev.off()