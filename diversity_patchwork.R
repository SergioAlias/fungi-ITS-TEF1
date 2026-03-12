# ╔═══════════════════════════════════════════════════════════════════╗
# ║                      diversity_patchwork.R                        ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : fungi-ITS-TEF1                                   ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-24                                       ║
# ║ Last Modified  : 2026-03-12                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(patchwork)

## Import plots

out <- "grano-ITS-TEF1"
amplicon <- "TEF1"
alpha_outdir <- file.path("/home/sergio/scratch",
                          out,
                          "alpha")
beta_outdir <- file.path("/home/sergio/scratch",
                         out,
                         "beta")
outdir <- file.path("/home/sergio/scratch",
                    out,
                    "diversity")

chao1_c <- readRDS(file.path(alpha_outdir, paste0("chao1_", amplicon, ".RDS")))
shannon_c <- readRDS(file.path(alpha_outdir, paste0("shannon_", amplicon, ".RDS")))
simpson_c <- readRDS(file.path(alpha_outdir, paste0("simpson_", amplicon, ".RDS")))
rpca_regular <- readRDS(file.path(beta_outdir, paste0("rPCA_regular_", amplicon, ".RDS")))
rpca_biplot <- readRDS(file.path(beta_outdir, paste0("rPCA_biplot_", amplicon, ".RDS")))

## Compose plot

alpha <- (chao1_c + theme(legend.position="none") +
          shannon_c + theme(legend.position="none") +
          simpson_c + theme(legend.position="none")) +
          plot_layout(ncol = 1)

if (amplicon == "ITS") {
  adonis_r2 <- 0.2484
  adonis_r2_adj <- 0.1463
} else if (amplicon == "TEF1") {
  adonis_r2 <- 0.0425
  adonis_r2_adj <- 0.057994
}

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

beta <- rpca_regular + 
  inset_element(
    rpca_biplot + 
      theme_void() + 
      theme(plot.background = element_rect(color = "black", fill = "white", linewidth = 0.5)), 
    left = 0.7, bottom = 0.015, right = 1, top = 0.255
  ) +
  plot_layout(guides = "collect",
              axes = "collect")

div <- (alpha | plot_spacer() | beta) +
  plot_layout(widths = c(2, 0.1, 5))

## Save plot

pdf(file.path(outdir, paste0("patched_div_", amplicon, ".pdf")),
    width = 11,
    height = 6)

div

dev.off()

png(file.path(outdir, paste0("patched_div_", amplicon, ".png")),
    width = 11,
    height = 6,
    units = "in",
    res = 300)

div

dev.off()
