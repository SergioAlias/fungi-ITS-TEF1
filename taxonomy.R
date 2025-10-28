# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            taxonomy.R                             ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : fungi-ITS-TEF1                                   ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-09-19                                       ║
# ║ Last Modified  : 2025-10-28                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggstats)
library(reshape2)
library(qiime2R)
library(ggpubr)
library(ggtext)


## Colors and shapes

source("/home/sergio/projects/fungi-ITS-TEF1/colors.R")


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

## Custom barplot

df <- read_qza(dada2_ITS_file_path)$data
df_long <- data.frame(
  row = rep(rownames(df), times = ncol(df)),
  column = rep(colnames(df), each = nrow(df)),
  value = as.vector(df)
)

tax <- read_qza(taxonomy_ITS_file_path)$data %>% parse_taxonomy()

tax$row <- rownames(tax)
rownames(tax) <- NULL

df_long %<>%
  left_join(tax %>% select(row, Genus), by = "row")


metadata <- read.csv(metadata_ITS_file_path, sep = "\t", header = TRUE)
metadata$column <- metadata$ID
metadata$ID <- NULL

df_long %<>%
  left_join(metadata %>% select(column, Cereal), by = "column")

df_long %<>%
  mutate(Genus = replace_na(Genus, "Unassigned"))

top_14_genera <- df_long %>%
  filter(Genus != "Unassigned" & !grepl("Incertae_sedis", Genus)) %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(value, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 14) %>%
  pull(Genus)

df_long %<>%
  mutate(
    Genus_plot = case_when(
      Genus %in% top_14_genera ~ Genus,
      Genus == "Unassigned" ~ "Unassigned",
      TRUE ~ "Other"
    )
  )

plot_levels <- c(top_14_genera, "Other", "Unassigned")
df_long$Genus_plot <- factor(df_long$Genus_plot, levels = plot_levels)

plot_colors <- setNames(
  c(barplot_ITS_colors[1:15], "grey80"),
  plot_levels
)

legend_labels <- sapply(plot_levels, function(label) {
  if (label %in% c("Unassigned", "Other")) {
    return(label) # Sin cursiva
  } else {
    return(bquote(italic(.(label))))
  }
})

df_long %<>% mutate(Cereal = str_replace(Cereal, "Triticum", "Wheat"),
                    Cereal = str_replace(Cereal, "AvenaSativa", "Oat"),
                    Cereal = str_replace(Cereal, "HordeumVulgare", "Barley"))


custom_barplot <- ggplot(df_long, aes(fill = Genus_plot, x = value, y = column)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(
    name = "Genus",
    values = plot_colors,
    breaks = plot_levels,
    labels = legend_labels
  ) +
  labs(x = NULL, y = NULL) +
  theme(
    text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.position = NULL,
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 12, color = "black", margin = margin(r = -2)),
    axis.ticks.x = element_blank(),
    strip.text = element_text(
      size = 15,
      color = "black",
      hjust = 0.5,
      vjust = 0.5
    ),
    strip.placement = "outside",
    strip.background = element_rect(fill = NA, color = NA)
  ) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  scale_y_discrete(expand = c(0, 0.7)) +
  facet_grid(Cereal ~ ., switch = "y", scale = "free_y") +
  guides(fill = guide_legend(nrow = 2, ncol = 8))

pdf(file.path(outdir, "custom_barplot.pdf"),
    width = 9)

custom_barplot

dev.off()

## Custom toxigenic genus barplot

target_genera <- c("Aspergillus", "Fusarium", "Penicillium")

rel_abundance_df <- df_long %>%
  filter(Genus %in% target_genera) %>%
  group_by(Cereal, column, Genus) %>%
  summarise(genus_reads = sum(value, na.rm = TRUE), .groups = 'drop') %>%
  left_join(
    df_long %>% group_by(column) %>% summarise(total_reads = sum(value, na.rm = TRUE)),
    by = "column"
  ) %>%
  filter(total_reads > 0) %>%
  mutate(rel_abundance = genus_reads / total_reads)

custom_boxplot <- ggplot(rel_abundance_df, aes(x = rel_abundance, y = Genus, fill = Genus, color = Genus)) + # <-- Añadir color aquí
  geom_boxplot(outlier.alpha = 0.5, na.rm = TRUE, alpha = 0.2) +
  scale_fill_manual(values = c("Aspergillus" = "black", plot_colors), guide = "none") +
  scale_color_manual(values = c("Aspergillus" = "black", plot_colors), guide = "none") +
  scale_y_discrete(limits = rev(target_genera),
                   labels = parse(text = paste0("italic('", rev(target_genera), "')")),
                   position = "right") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = NULL, y = NULL) +
  facet_grid(Cereal ~ .) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    panel.spacing = unit(1, "lines")
  )


final_plot <- ggarrange(custom_barplot, custom_boxplot,
                        common.legend = TRUE,
                        legend = "top",
                        widths = c(6, 2))

final_plot <- annotate_figure(final_plot,
                              bottom = text_grob("Relative abundance", 
                                                 color = "black", 
                                                 size = 16))
pdf(file.path(outdir, "custom_combined.pdf"),
    width = 16,
    height = 10)

final_plot

dev.off()

png(file.path(outdir, "custom_combined.png"),
    width = 14,
    height = 10,
    units = "in",
    res = 300)

final_plot

dev.off()

## Bubble plot

tef1_wide <- read_qza("/home/sergio/projects/fungi-ITS-TEF1/local/TEF1_filtered_table.qza")$data
tef1_wide %<>% t() %>% as.data.frame()
tef1_taxa <- read_qza("/home/sergio/projects/fungi-ITS-TEF1/local/TEF1_taxonomy.qza")$data %>% parse_taxonomy()
tef1_taxa %<>% .[.$Genus == "Fusarium" & !is.na(.$Genus), ]
tef1_taxa[is.na(tef1_taxa$Species), "Species"] <- rep("sp", sum(is.na(tef1_taxa$Species)))
tef1_taxa %<>% mutate(Species = gsub("'", "", Species))
tef1_taxa %<>%mutate(Species = ifelse(!grepl("^F", Species),
                          paste(Genus, Species, sep = "_"),
                          Species))
tef1_wide %<>% .[, colnames(.) %in% rownames(tef1_taxa)] %>% {`colnames<-`(., tef1_taxa[colnames(.), "Species"])}
tef1_wide <- tef1_wide[, colnames(tef1_wide) != "Fusarium_dimerum"] # Fusarium dimerum is now Bisifusarium dimerum -> https://www.fusarium.org/page/SpeciesListBisifusarium
colnames(tef1_wide) <- gsub("\\.[0-9]+$", "", colnames(tef1_wide))
tef1_wide <- sapply(unique(colnames(tef1_wide)), function(x) rowSums(tef1_wide[names(tef1_wide) == x])) %>% as.data.frame()
tef1_meta <- read.csv(metadata_TEF1_file_path, header = TRUE, sep = "\t")
tef1_wide$ID <- rownames(tef1_wide)
tef1_wide %<>% merge(., tef1_meta[, c("ID", "Cereal")], by = "ID")
tef1_wide <- aggregate(. ~ Cereal, data = tef1_wide[,-1], FUN = sum)
rownames(tef1_wide) <- tef1_wide$Cereal
tef1_wide$Cereal <- NULL
tef1_wide %<>% cbind(rownames(.), .) %>% {`rownames<-`(., NULL)}
tef1_wide[, 2:ncol(tef1_wide)] <- t(apply(tef1_wide[, 2:ncol(tef1_wide)], 1, function(x) {
  row_sum <- sum(x)
  if(row_sum == 0) {
    return(rep(0, length(x)))
  } else {
    return((x / row_sum) * 100)
  }
}))
tef1_long <- melt(tef1_wide, id.vars = c("rownames(.)"))
colnames(tef1_long) <- c("Sample", "variable", "value")
tef1_long  %<>%
  mutate(Amplicon = case_when(
    TRUE ~ "TEF1"))


its_wide <- read_qza("/home/sergio/projects/fungi-ITS-TEF1/local/ITS_filtered_table.qza")$data
its_wide %<>% t() %>% as.data.frame()
its_taxa <- read_qza("/home/sergio/projects/fungi-ITS-TEF1/local/ITS_taxonomy.qza")$data %>% parse_taxonomy()
its_taxa %<>% .[.$Genus == "Fusarium" & !is.na(.$Genus), ]
its_taxa[is.na(its_taxa$Species), "Species"] <- rep("Fusarium_sp", sum(is.na(its_taxa$Species)))
its_wide %<>% .[, colnames(.) %in% rownames(its_taxa)] %>% {`colnames<-`(., its_taxa[colnames(.), "Species"])}
names(its_wide)[names(its_wide) == "Fusarium_lunatum"] <- "Fusarium_sp"
its_wide <- sapply(unique(colnames(its_wide)), function(x) rowSums(its_wide[names(its_wide) == x])) %>% as.data.frame()
its_meta <- read.csv(metadata_ITS_file_path, header = TRUE, sep = "\t")
its_wide$ID <- rownames(its_wide)
its_wide %<>% merge(., its_meta[, c("ID", "Cereal")], by = "ID")
its_wide <- aggregate(. ~ Cereal, data = its_wide[,-1], FUN = sum)
rownames(its_wide) <- its_wide$Cereal
its_wide$Cereal <- NULL
its_wide %<>% cbind(rownames(.), .) %>% {`rownames<-`(., NULL)}
its_wide[, 2:ncol(its_wide)] <- t(apply(its_wide[, 2:ncol(its_wide)], 1, function(x) {
  row_sum <- sum(x)
  if(row_sum == 0) {
    return(rep(0, length(x)))
  } else {
    return((x / row_sum) * 100)
  }
}))
its_long <- melt(its_wide, id.vars = c("rownames(.)"))
colnames(its_long) <- c("Sample", "variable", "value")
its_long  %<>%
  mutate(Amplicon = case_when(
                    TRUE ~ "ITS2"))

comb_long <- rbind(its_long, tef1_long)

comb_long  %<>%
  mutate(sc = case_when(
  variable == "Fusarium_sp" ~ " ",
  variable == "Fusarium_tricinctum" ~ "FTSC",
  variable == "Fusarium_oxysporum" ~ "FOSC",
  variable == "Fusarium_sporotrichioides" ~ "FSAMSC",
  variable == "Fusarium_redolens" ~ "FRSC",
  variable == "Fusarium_poae" ~ "FSAMSC",
  variable == "Fusarium_subglutinans" ~ "FFSC",
  variable == "Fusarium_scirpi" ~ "FIESC",
  variable == "Fusarium_coffeatum" ~ "FIESC",
  variable == "Fusarium_hostae" ~ "FRSC",
  variable == "Fusarium_iranicum" ~ "FTSC",
  variable == "Fusarium_clavum" ~ "FIESC",
  variable == "Fusarium_avenaceum" ~ "FTSC",
  variable == "Fusarium_equiseti" ~ "FIESC",
  variable == "Fusarium_pseudograminearum" ~ "FSAMSC",
  variable == "Fusarium_acuminatum" ~ "FTSC",
  variable == "Fusarium_langsethiae" ~ "FSAMSC",
  variable == "Fusarium_torolosum" ~ "FTSC",
  variable == "Fusarium_flagelliforme" ~ "FIESC",
  variable == "FSAMSC26" ~ "FSAMSC",
  variable == "Fusarium_proliferatum" ~ "FFSC",
  variable == "Fusarium_graminearum" ~ "FSAMSC",
  variable == "Fusarium_gracilipes" ~ "FIESC",
  variable == "FFSCundesc" ~ "FFSC",
  variable == "Fusarium_flocciferum" ~ "FTSC",
  variable == "Fusarium_gamsii" ~ "FTSC",
  variable == "FSAMSC25" ~ "FSAMSC",
  variable == "Fusarium_lactis" ~ "FFSC",
  variable == "Fusarium_culmorum" ~ "FSAMSC",
  variable == "FTSCundesc" ~ "FTSC",
  variable == "Fusarium_cerealis" ~ "FSAMSC",
  variable == "Fusarium_nodosum" ~ "FSAMSC",
  variable == "Fusarium_nelsonii" ~ "FCSC",
  TRUE ~ "YOU SHOULD NOT BE SEEING THIS"
))

comb_long %<>% mutate(variable = str_replace(variable, "^Fusarium_", "F. "),
                      Sample = str_replace(Sample, "Triticum", "Wheat"),
                      Sample = str_replace(Sample, "AvenaSativa", "Oat"),
                      Sample = str_replace(Sample, "HordeumVulgare", "Barley"))


comb_bubplot <- ggplot(comb_long, aes(x = Amplicon, y = variable)) +
  geom_point(aes(size = value, fill = Amplicon, color = Amplicon), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1, 10, 50, 75)) +
  scale_fill_manual(values = amplicon_colors) +
  scale_color_manual(values = amplicon_colors) +
  labs(x = "", y = "", size = "Relative abundance (%)", fill = "") +
    guides(fill = guide_legend(override.aes = list(size = 6,
                                                   shape = 21,
                                                   color = amplicon_colors,
                                                   fill  = amplicon_colors)),
           color = "none",
           size = guide_legend()) +
  facet_grid(sc~Sample, scales = "free_y", space="free_y", switch = "y") +
  theme(text = element_text(size = 15),
        legend.key = element_blank(),
        legend.title = element_text(size = 12),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, linewidth = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "italic"),
        strip.placement = "outside",
        strip.text = element_text(
          size = 11,                   
          color = "black",                 
          hjust = 0.5,                     
          vjust = 0.5                      
        ),
        strip.text.y.left = element_text(
          angle = 0,
          hjust = 1),
        strip.background = element_rect(
          fill = NA,
          color = NA
        )) +
  scale_y_discrete(limits = rev, position = "right")

pdf(file.path(outdir, "bubble_plot.pdf"),
    width = 9,
    height = 12)

comb_bubplot

dev.off()

png(file.path(outdir, "bubble_plot.png"),
    width = 9,
    height = 12,
    units = "in",
    res = 300)

comb_bubplot

dev.off()
