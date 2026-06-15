# ╔═══════════════════════════════════════════════════════════════════╗
# ║                 cross_amplicon_correlations.R                     ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : fungi-ITS-TEF1                                   ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2026-06-11                                       ║
# ║ Last Modified  : 2026-06-15                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

if (!requireNamespace("CompoCor", quietly = TRUE)) {
  library(devtools)
  install_github("IbTJensen/CompoCor")
}

if (!requireNamespace("qiime2R", quietly = TRUE)) {
  library(devtools)
  install_github("jbisanz/qiime2R")
}

library(magrittr, include.only = "%<>%")
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(CompoCor)
library(pheatmap)
library(gridExtra)
library(ggplotify)
library(patchwork)
library(ggpubr)

## Functions

fromQ2toPhyloTable <- function(project) {
  table <- file.path(cluster_path,
                     base_path,
                     project,
                     table_path)
  table %<>% read_qza() %>%
    {.$data} %>%
    as.matrix()
  colnames(table) <- sapply(strsplit(colnames(table), "_"), function(x) {
    paste(head(x, -1), collapse = "_")
  })
  table %<>% otu_table(taxa_are_rows = TRUE)
  return(table)
}

doPheatmap <- function(mat,
                       gaps_row = NULL,
                       gaps_col = NULL,
                       breaksList = seq(-1, 1, by = 0.02),
                       legendbreaks = c(-0.95, -0.5, 0, 0.5, 0.95),
                       legendlabels = c(-1, -0.5, 0, 0.5, 1),
                       angle = 90,
                       fontsize = 6,
                       bordercolor = "grey60",
                       labelsrow = NULL,
                       labelscol = NULL,
                       drawlegend = TRUE) {
  return(pheatmap(mat,
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  legend = drawlegend,
                  gaps_row = gaps_row,
                  gaps_col = gaps_col,
                  color = colorRampPalette(c("orange",
                                             "white",
                                             "purple"))(100),
                  scale = "none",
                  show_rownames = TRUE,
                  show_colnames = TRUE,
                  angle_col = angle,
                  fontsize = fontsize,
                  fontsize_row = fontsize,
                  fontsize_col = fontsize,
                  legend_breaks = legendbreaks,
                  legend_labels = legendlabels,
                  breaks = breaksList,
                  border_color = bordercolor,
                  labels_row = labelsrow,
                  labels_col = labelscol)
  )
}


## Import QIIME 2 files

set.seed(1234)
n_cores <- 10

fungi_project <- "grano_ITS"
fusarium_project <- "grano_TEF1"
base_path <- "scratch/salias/projects"
table_path <- "qiime2/feature_tables/filtered_table.qza"

readRenviron("/home/sergio/Renvs/.RenvBrigit")
brigit_IP <- Sys.getenv("IP_ADDRESS")
cluster_path <- paste0("/run/user/1001/gvfs/sftp:host=",
                       brigit_IP,
                       ",user=salias/mnt/lustre")

metadata <- sample_data(read.table(file.path(cluster_path,
                                             "home/salias/projects/cross_amplicon/metadata.tsv"),
                                   row.names = 1,
                                   sep = "\t",
                                   quote = "",
                                   header = TRUE)
)

fungi_table <- fromQ2toPhyloTable(fungi_project)
sample_names(fungi_table) <- sample_names(metadata)

fusarium_table <- fromQ2toPhyloTable(fusarium_project)
sample_names(fusarium_table) <- sample_names(metadata)

fungi_table %<>% as.data.frame() %>% t()
fusarium_table %<>% as.data.frame() %>% t()

fungi_table_bar <- fungi_table[grepl("_BAR", rownames(fungi_table)), ]
fusarium_table_bar <- fusarium_table[grepl("_BAR", rownames(fusarium_table)), ]

fungi_table_oat <- fungi_table[grepl("_OAT", rownames(fungi_table)), ]
fusarium_table_oat <- fusarium_table[grepl("_OAT", rownames(fusarium_table)), ]

fungi_table_whe <- fungi_table[grepl("_WHE", rownames(fungi_table)), ]
fusarium_table_whe <- fusarium_table[grepl("_WHE", rownames(fusarium_table)), ]

list2env(lapply(list(fungi_table_bar = fungi_table_bar, 
                     fusarium_table_bar = fusarium_table_bar, 
                     fungi_table_oat = fungi_table_oat, 
                     fusarium_table_oat = fusarium_table_oat,
                     fungi_table_whe = fungi_table_whe, 
                     fusarium_table_whe = fusarium_table_whe),
                function(mat) mat[, colSums(mat) != 0]),
         envir = .GlobalEnv)


sparxcc_bar <- SparXCC_base(fungi_table_bar,
                            fusarium_table_bar,
                            pseudo_count = 1,
                            cores = n_cores)

sparxcc_oat <- SparXCC_base(fungi_table_oat,
                            fusarium_table_oat,
                            pseudo_count = 1,
                            cores = n_cores)

sparxcc_whe <- SparXCC_base(fungi_table_whe,
                            fusarium_table_whe,
                            pseudo_count = 1,
                            cores = n_cores)


cor_bar <- sparxcc_bar[["cor"]]
cor_oat <- sparxcc_oat[["cor"]]
cor_whe <- sparxcc_whe[["cor"]]

fungi_to_keep <- c("Cladosporium",
                   "Alternaria",
                   "Vishniacozyma",
                   "Stemphylium",
                   "Nigrospora",
                   "Aureobasidium",
                   "Filobasidium",
                   "Penicillium",
                   "Cystofilobasidium",
                   "Sarocladium",
                   "Acremonium",
                   "Dioszegia")

# (I know the following code is a mess)

taxa_fungi <- parse_taxonomy(read_qza(file.path(cluster_path,
                                                base_path,
                                                fungi_project,
                                                "qiime2/taxonomy/taxonomy.qza"))[["data"]])

taxa_fusarium <- parse_taxonomy(read_qza(file.path(cluster_path,
                                                   base_path,
                                                   fusarium_project,
                                                   "qiime2/taxonomy/taxonomy.qza"))[["data"]])

rownames(cor_bar) <- gsub("_", " ",
                          paste0(taxa_fungi$Species[match(rownames(cor_bar),
                                                          rownames(taxa_fungi))], 
                                 " (",
                                 substr(rownames(cor_bar), 1, 7),
                                 ")"))

colnames(cor_bar) <- gsub("_", " ",
                          paste0(taxa_fusarium$Genus[match(colnames(cor_bar),
                                                           rownames(taxa_fusarium))],
                                 " - ",
                                 taxa_fusarium$Species[match(colnames(cor_bar),
                                                             rownames(taxa_fusarium))], 
                                 " (",
                                 substr(colnames(cor_bar), 1, 7),
                                 ")"))

colnames(cor_bar) <- gsub("'", "", colnames(cor_bar))

rownames(cor_oat) <- gsub("_", " ",
                          paste0(taxa_fungi$Species[match(rownames(cor_oat),
                                                          rownames(taxa_fungi))], 
                                 " (",
                                 substr(rownames(cor_oat), 1, 7),
                                 ")"))

colnames(cor_oat) <- gsub("_", " ",
                          paste0(taxa_fusarium$Genus[match(colnames(cor_oat),
                                                           rownames(taxa_fusarium))],
                                 " - ",
                                 taxa_fusarium$Species[match(colnames(cor_oat),
                                                             rownames(taxa_fusarium))], 
                                 " (",
                                 substr(colnames(cor_oat), 1, 7),
                                 ")"))

colnames(cor_oat) <- gsub("'", "", colnames(cor_oat))

rownames(cor_whe) <- gsub("_", " ",
                          paste0(taxa_fungi$Species[match(rownames(cor_whe),
                                                          rownames(taxa_fungi))], 
                                 " (",
                                 substr(rownames(cor_whe), 1, 7),
                                 ")"))

colnames(cor_whe) <- gsub("_", " ",
                          paste0(taxa_fusarium$Genus[match(colnames(cor_whe),
                                                           rownames(taxa_fusarium))],
                                 " - ",
                                 taxa_fusarium$Species[match(colnames(cor_whe),
                                                             rownames(taxa_fusarium))], 
                                 " (",
                                 substr(colnames(cor_whe), 1, 7),
                                 ")"))

colnames(cor_whe) <- gsub("'", "", colnames(cor_whe))

cor_bar <- cor_bar[grepl(paste(fungi_to_keep, collapse = "|"), rownames(cor_bar)), ]
cor_oat <- cor_oat[grepl(paste(fungi_to_keep, collapse = "|"), rownames(cor_oat)), ]
cor_whe <- cor_whe[grepl(paste(fungi_to_keep, collapse = "|"), rownames(cor_whe)), ]


cor_bar <- cor_bar[order(rownames(cor_bar)), order(colnames(cor_bar))]
cor_oat <- cor_oat[order(rownames(cor_oat)), order(colnames(cor_oat))]
cor_whe <- cor_whe[order(rownames(cor_whe)), order(colnames(cor_whe))]

fusariumRenamer <- function(x) {
  if (grepl("- NA", x)) {
    return(gsub("- NA", "sp", x))
  } else if (grepl("SC", x)) {
    return(strsplit(x, " - ")[[1]][2])
  } else {
    return(gsub(" - ", " ", x))
  }
}

colnames(cor_bar) <- sapply(colnames(cor_bar), fusariumRenamer)
colnames(cor_oat) <- sapply(colnames(cor_oat), fusariumRenamer)
colnames(cor_whe) <- sapply(colnames(cor_whe), fusariumRenamer)

italic_rownames_bar <- as.expression(lapply(
  rownames(cor_bar),
  function(x) bquote(italic(.(x)))))

italic_colnames_bar <- as.expression(lapply(
  colnames(cor_bar),
  function(x) bquote(italic(.(x)))))

italic_rownames_oat <- as.expression(lapply(
  rownames(cor_oat),
  function(x) bquote(italic(.(x)))))

italic_colnames_oat <- as.expression(lapply(
  colnames(cor_oat),
  function(x) bquote(italic(.(x)))))

italic_rownames_whe <- as.expression(lapply(
  rownames(cor_whe),
  function(x) bquote(italic(.(x)))))

italic_colnames_whe <- as.expression(lapply(
  colnames(cor_whe),
  function(x) bquote(italic(.(x)))))

bar_plot <- doPheatmap(mat = cor_bar,
                       angle = 45,
                       fontsize = 8,
                       bordercolor = "black",
                       labelsrow = italic_rownames_bar,
                       labelscol = italic_colnames_bar,
                       gaps_col = c(2, 6, 7, 11, 12, 15, 16, 17, 18, 19, 20, 22, 24, 26, 27, 31, 32, 37, 39, 41, 43, 44),
                       gaps_row = c(4, 9, 12, 22, 26, 33, 39, 40, 44, 52, 57))
oat_plot <- doPheatmap(mat = cor_oat,
                       angle = 45,
                       fontsize = 8,
                       bordercolor = "black",
                       labelsrow = italic_rownames_oat,
                       labelscol = italic_colnames_oat,
                       gaps_col = c(2, 3, 10, 11, 12, 17, 18, 19, 20, 22, 24, 25, 27, 29, 33, 34, 35, 46, 47, 50, 51, 52, 53, 54),
                       gaps_row = c(1, 6, 9, 23, 28, 33, 37, 38, 40, 45, 48))
whe_plot <- doPheatmap(mat = cor_whe,
                       angle = 45,
                       fontsize = 8,
                       bordercolor = "black",
                       labelsrow = italic_rownames_whe,
                       labelscol = italic_colnames_whe,
                       gaps_col = c(2, 3, 9, 10, 14, 16, 17, 18, 21, 22, 23, 25, 28, 29, 31, 34, 35, 38, 40, 41),
                       gaps_row = c(1, 5, 7, 20, 24, 32, 37, 41, 42, 48, 53))

png("/home/salias/scratch/cross_amplicon/asv_heatmap.png",
    width = 18,
    height = 22,
    units = "in",
    res = 300)

(plot_spacer() | wrap_plots(list(as.ggplot(bar_plot) + ggtitle("Barley"),  
                                 as.ggplot(oat_plot) + ggtitle("Oat"),  
                                 as.ggplot(whe_plot) + ggtitle("Wheat")), 
                            ncol = 1)) + 
  plot_layout(widths = c(1, 20))

dev.off()

aggregate_matrix <- function(matrix_data, row_groups) {
  col_names <- colnames(matrix_data)
  col_groups <- unique(sapply(strsplit(col_names, " (", fixed = TRUE), `[`, 1))
  row_indices <- lapply(row_groups, function(pattern) grep(paste0("^", pattern), rownames(matrix_data)))
  col_indices <- lapply(col_groups, function(pattern) grep(paste0("^", pattern), colnames(matrix_data)))
  aggregated_matrix <- matrix(NA, nrow = length(row_groups), ncol = length(col_groups),
                              dimnames = list(row_groups, col_groups))
  for (i in seq_along(row_groups)) {
    for (j in seq_along(col_groups)) {
      selected_values <- matrix_data[row_indices[[i]], col_indices[[j]], drop = FALSE]
      aggregated_matrix[i, j] <- mean(selected_values, na.rm = TRUE)
    }
  }
  return(aggregated_matrix)
}


aggr_cor_bar <- aggregate_matrix(cor_bar, fungi_to_keep)
aggr_cor_oat <- aggregate_matrix(cor_oat, fungi_to_keep)
aggr_cor_whe <- aggregate_matrix(cor_whe, fungi_to_keep)


italic_rownames_aggr_bar <- as.expression(lapply(
  rownames(aggr_cor_bar),
  function(x) bquote(italic(.(x)))))

italic_colnames_aggr_bar <- as.expression(lapply(
  colnames(aggr_cor_bar),
  function(x) bquote(italic(.(x)))))

italic_rownames_aggr_oat <- as.expression(lapply(
  rownames(aggr_cor_oat),
  function(x) bquote(italic(.(x)))))

italic_colnames_aggr_oat <- as.expression(lapply(
  colnames(aggr_cor_oat),
  function(x) bquote(italic(.(x)))))

italic_rownames_aggr_whe <- as.expression(lapply(
  rownames(aggr_cor_whe),
  function(x) bquote(italic(.(x)))))

italic_colnames_aggr_whe <- as.expression(lapply(
  colnames(aggr_cor_whe),
  function(x) bquote(italic(.(x)))))

aggr_bar_plot <- doPheatmap(mat = aggr_cor_bar,
                            angle = 45,
                            fontsize = 8,
                            bordercolor = "black",
                            labelsrow = italic_rownames_aggr_bar,
                            labelscol = italic_colnames_aggr_bar)

aggr_oat_plot <- doPheatmap(mat = aggr_cor_oat,
                            angle = 45,
                            fontsize = 8,
                            bordercolor = "black",
                            labelsrow = italic_rownames_aggr_oat,
                            labelscol = italic_colnames_aggr_oat)

aggr_whe_plot <- doPheatmap(mat = aggr_cor_whe,
                            angle = 45,
                            fontsize = 8,
                            bordercolor = "black",
                            labelsrow = italic_rownames_aggr_whe,
                            labelscol = italic_colnames_aggr_whe)


png("/home/salias/scratch/cross_amplicon/aggregated_heatmap.png",
    width = 10,
    height = 10,
    units = "in",
    res = 300)

(plot_spacer() | wrap_plots(list(as.ggplot(aggr_bar_plot) + ggtitle("Barley"),  
                                 as.ggplot(aggr_oat_plot) + ggtitle("Oat"),  
                                 as.ggplot(aggr_whe_plot) + ggtitle("Wheat")), 
                            ncol = 1)) + 
  plot_layout(widths = c(1, 20))

dev.off()
