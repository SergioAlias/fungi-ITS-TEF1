# ╔═══════════════════════════════════════════════════════════════════╗
# ║                      diversity_patchwork.R                        ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : fungi-ITS-TEF1                                   ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-24                                       ║
# ║ Last Modified  : 2026-01-21                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(patchwork)

## Import plots

out <- "grano-ITS-TEF1"
alpha_outdir <- file.path("/home/sergio/scratch",
                          out,
                          "alpha")
beta_outdir <- file.path("/home/sergio/scratch",
                         out,
                         "beta")
outdir <- file.path("/home/sergio/scratch",
                    out,
                    "diversity")

chao1_c <- readRDS(file.path(alpha_outdir, "chao1.RDS"))
shannon_c <- readRDS(file.path(alpha_outdir, "shannon.RDS"))
simpson_c <- readRDS(file.path(alpha_outdir, "simpson.RDS"))
rpca_regular <- readRDS(file.path(beta_outdir, "rPCA_regular.RDS"))
rpca_biplot <- readRDS(file.path(beta_outdir, "rPCA_biplot.RDS"))

## Compose plot

alpha <- (chao1_c + theme(legend.position="none") +
          shannon_c + theme(legend.position="none") +
          simpson_c + theme(legend.position="none")) +
          plot_layout(ncol = 1)

adonis_r2 <- 0.2484     # hardcoded Adonis results for beta div plot
adonis_r2_adj <- 0.1463 # same

adonis_text <- paste0("Adonis: R² = ",
                      round(adonis_r2, 3), 
                      ", Adjusted R² = ",
                      round(adonis_r2_adj, 3))

rpca_regular <- rpca_regular + 
  labs(subtitle = adonis_text) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 11, face = "plain"))

beta <- rpca_regular / rpca_biplot +
  plot_layout(guides = "collect",
              axes = "collect")

div <- (alpha | plot_spacer() | beta) +
  plot_layout(widths = c(2, 0.1, 3))

## Save plot

pdf(file.path(outdir, "patched_div.pdf"),
    width = 10,
    height = 6)

div

dev.off()

png(file.path(outdir, "patched_div.png"),
    width = 10,
    height = 6,
    units = "in",
    res = 300)

div

dev.off()
