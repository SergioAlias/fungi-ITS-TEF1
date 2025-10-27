# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            map.R                                  ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : fungi-ITS-TEF1                                   ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-10-27                                       ║
# ║ Last Modified  : 2025-10-27                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/fungi-ITS-TEF1    ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(sf)
library(mapSpain)
library(magrittr, include.only = "%<>%")
library(tidyverse)
library(ggrepel)

## Setup

out <- "grano-ITS-TEF1"
outdir <- file.path("/home/sergio/scratch",
                    out,
                    "maps")
source("/home/sergio/projects/fungi-ITS-TEF1/colors.R")

## Load map resources

esp <- esp_get_prov(moveCAN = FALSE)
esp %<>%
  mutate(n_samples = case_when(
    iso2.prov.code == "ES-VI" ~ "1", # Álava
    iso2.prov.code == "ES-BU" ~ "2", # Burgos
    iso2.prov.code == "ES-CU" ~ "2", # Cuenca
    iso2.prov.code == "ES-GU" ~ "1", # Guadalajara
    iso2.prov.code == "ES-M" ~ "2", # Madrid
    iso2.prov.code == "ES-SA" ~ "5", # Salamanca
    iso2.prov.code == "ES-SG" ~ "3", # Segovia
    iso2.prov.code == "ES-SO" ~ "2", # Soria
    iso2.prov.code == "ES-TE" ~ "1", # Teruel
    iso2.prov.code == "ES-VA" ~ "10", # Valladolid
    iso2.prov.code == "ES-ZA" ~ "1", # Zamora
    .default = "0"
  ))

locations <- data.frame(
  lng = c(-2.67, -3.70, -2.13, -3.16, -3.70, -5.66, -4.12, -2.46, -1.10, -4.72, -5.75),
  lat = c(42.85, 42.34, 40.07, 40.63, 40.41, 40.97, 40.94, 41.76, 40.34, 41.65, 41.50),
  ID = c("Álava (n = 1)\n· 1 wheat sample",
         "Burgos (n = 2)\n· 1 barley sample\n· 1 wheat sample",
         "Cuenca (n = 2)\n· 2 barley samples",
         "Guadalajara (n = 1)\n· 1 barley sample",
         "Madrid (n = 2)\n· 1 barley sample\n· 1 wheat sample",
         "Salamanca (n = 5)\n· 2 barley samples\n· 1 oat sample\n· 2 wheat samples",
         "Segovia (n = 3)\n· 1 barley sample\n· 1 oat sample\n· 1 wheat sample",
         "Soria (n = 2)\n· 1 barley sample\n· 1 wheat sample",
         "Teruel (n = 1)\n· 1 wheat sample",
         "Valladolid (n = 10)\n· 2 barley samples\n· 4 oat samples\n· 4 wheat samples",
         "Zamora (n = 1)\n· 1 oat sample"))

locations %<>% st_as_sf(coords = c("lng", "lat"), remove = FALSE, 
                        crs = 4326, agr = "constant")

# A vector of 11 x-nudges
nudge_x_values <- c(
  0.1,  # Alava
  0,    # Burgos
  0,    # Cuenca
  0.2,  # Guadalajara (push right)
  0,    # Madrid
  -0.2, # Salamanca (push left)
  0.1,  # Segovia (push right)
  0,    # Soria
  0,    # Teruel
  0,    # Valladolid
  -0.2  # Zamora (push left)
)

# A vector of 11 y-nudges
nudge_y_values <- c(
  0.1,  # Alava (push up)
  -0.1, # Burgos (push down)
  -0.1, # Cuenca (push down)
  0,    # Guadalajara
  -0.2, # Madrid (push down)
  0,    # Salamanca
  0.2,  # Segovia (push up)
  0,    # Soria
  0,    # Teruel
  0.1,  # Valladolid (push up)
  0     # Zamora
)

## Plot map

pdf(file.path(outdir, "sample_map.pdf"),
    width = 9)

### PNG version

# png(file.path(outdir, "sample_map.png"),
#     width = 9,
#     height = 7,
#     units = "in",
#     res = 300)

ggplot(esp) +
  geom_sf(aes(fill = n_samples)) +
  scale_fill_discrete(type = sample_map_colors,
                      guide = "none") +
  geom_sf(data = locations) +
  geom_label_repel(data = locations,
                   aes(x = lng, y = lat,
                       label = ID),
                   parse = FALSE,
                   hjust = 0,
                   max.overlaps = Inf,
                   force = 2,
                   box.padding = unit(0.5, "lines"),
                   min.segment.length = unit(0, "lines")) +
                   # nudge_x = nudge_x_values,
                   # nudge_y = nudge_y_values) +
  xlab(NULL) + ylab(NULL) +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.2),
        panel.background = element_rect(fill = "aliceblue")) +
  coord_sf(xlim = c(-11, 6), ylim = c(35.5, 44.5), expand = TRUE) # + theme_void()

dev.off()

