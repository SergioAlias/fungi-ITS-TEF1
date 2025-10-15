# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            taxonomy.R                             ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : fungi-ITS-TEF1                                   ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-09-19                                       ║
# ║ Last Modified  : 2025-10-15                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(dplyr)
library(stringr)
library(file2meco)
library(microeco)
library(ggplot2)
library(ggstats)
library(ggnested)
library(ggh4x)
library(reshape2)
library(qiime2R)

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

meco_ITS <- qiime2meco(dada2_ITS_file_path,
                       sample_table = metadata_ITS_file_path,
                       taxonomy_table = taxonomy_ITS_file_path)

meco_TEF1 <- qiime2meco(dada2_TEF1_file_path,
                        sample_table = metadata_TEF1_file_path,
                        taxonomy_table = taxonomy_TEF1_file_path)

## Change amplicon to plot

meco <- meco_ITS # meco_TEF1

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
  group_by(Genus) %>%
  mutate(total_abundance = sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Genus = reorder(Genus, -total_abundance))

top <- c(levels(df_long$Genus)[!grepl("Incertae_sedis", levels(df_long$Genus))][1:15], NA)

legend_labels <- sapply(top, function(label) {
  if (is.na(label)) {
    "Unassigned"
  } else {
    bquote(italic(.(label)))
  }
})

df_long %<>% mutate(Cereal = str_replace(Cereal, "Triticum", "Triticum sp."),
                    Cereal = str_replace(Cereal, "AvenaSativa", "Avena sativa"),
                    Cereal = str_replace(Cereal, "HordeumVulgare", "Hordeum vulgare"))

custom_barplot <- ggplot(df_long, aes(fill=Genus, y=value, x=column)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(
    values = setNames(c(barplot_ITS_colors, "grey"), top),
    breaks = top,
    labels = legend_labels
  ) +
  labs(x = "", y = "Relative abundance (%)") +
  theme(text = element_text(size = 15),
        legend.title = element_text(size = 15),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 14, margin = margin(r =-2)),
        axis.ticks.y = element_blank(), 
        strip.text = element_text(
          size = 15,                       
          face = "italic",                   
          color = "black",                 
          hjust = 0.5,                     
          vjust = 0.5                      
        ),
        strip.background = element_rect(
          fill = NA,
          color = NA
        )) + 
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) x * 100) +
  scale_x_discrete(expand = c(0, 0.7)) +
  facet_grid(~Cereal, scale = "free_x", space = "free_x")

pdf(file.path(outdir, "custom_barplot.pdf"),
    width = 9)

custom_barplot

dev.off()


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
    width = 14,
    height = 5.5)

t_stacked_family$plot_bar(ggnested = TRUE,
                          high_level_add_other = TRUE,
                          xtext_keep = FALSE,
                          color_values = get(paste0("barplot_",
                                                    color_palette,
                                                    "_colors")),
                          # xtext_angle = 90,
                          # xtext_size = 6,
                          facet = c("Location")) + 
  theme(ggh4x.facet.nestline = element_line(colour = "grey95"))

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
                      Sample = str_replace(Sample, "Triticum", "Triticum sp."),
                      Sample = str_replace(Sample, "AvenaSativa", "Avena sativa"),
                      Sample = str_replace(Sample, "HordeumVulgare", "Hordeum vulgare"))


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
          face = "bold.italic",                   
          color = "black",                 
          hjust = 0.5,                     
          vjust = 0.5                      
        ),
        strip.text.y.left = element_text(
          angle = 0,
          hjust = 1,
          face = "bold"),
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
