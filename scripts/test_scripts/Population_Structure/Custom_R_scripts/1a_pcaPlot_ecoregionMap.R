# ==========================================================
# SET WORKING DIRECTORY
# ==========================================================
setwd("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics")

# ==========================================================
# LOAD PACKAGES
# ==========================================================
library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(rnaturalearth)
library(ggnewscale)

# ==========================================================
# LOAD ECOREGIONS (FULL PATH)
# ==========================================================
eco <- st_read("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared-data/shapefiles/us_eco_l3_state_boundaries/us_eco_l3_state_boundaries.shp",
               quiet = TRUE)

eco_west <- eco %>%
  filter(STATE_NAME %in% c("California", "Oregon", "Washington"))

# ==========================================================
# LOAD ALL METADATA
# ==========================================================
all_points_df <- bind_rows(
  read_csv("data/Chamaea/Chamaea_meta.csv", show_col_types = FALSE) %>% mutate(species = "Chamaea"),
  read_csv("data/Picoides/Picoides_meta.csv", show_col_types = FALSE) %>% mutate(species = "Picoides"),
  read_csv("data/Poecile/Poecile_meta.csv", show_col_types = FALSE) %>% mutate(species = "Poecile"),
  read_csv("data/Toxostoma/Toxostoma_meta.csv", show_col_types = FALSE) %>% mutate(species = "Toxostoma")
) %>%
  filter(!is.na(lat), !is.na(long))

# ==========================================================
# CONVERT TO SF
# ==========================================================
points_sf <- st_as_sf(all_points_df, coords = c("long", "lat"), crs = 4326)
points_sf <- st_transform(points_sf, st_crs(eco_west))

# ==========================================================
# ECOREGIONS USED
# ==========================================================
eco_levels <- sort(unique(all_points_df$ecoregion))

eco_sampled <- eco_west %>%
  filter(US_L3NAME %in% eco_levels)

# ==========================================================
# COLOR PALETTE
# ==========================================================
base_cols <- c(
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Paired")
)

eco_colors <- setNames(base_cols[1:length(eco_levels)], eco_levels)

# ==========================================================
# STATES
# ==========================================================
states <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(name %in% c("California", "Oregon", "Washington")) %>%
  st_transform(st_crs(eco_west))

# ==========================================================
# MAP PLOT (LEGEND FIXED)
# ==========================================================
map_plot <- ggplot() +
  geom_sf(data = states, fill = "gray90", color = "black", size = 0.4) +
  
  geom_sf(
    data = eco_sampled,
    aes(fill = US_L3NAME),
    color = "black",
    size = 0.2
  ) +
  
  scale_fill_manual(
    values = eco_colors,
    guide = guide_legend(order = 2, ncol = 3, byrow = TRUE)
  ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_sf(
    data = points_sf %>%
      mutate(point_type = ifelse(species %in% c("Chamaea", "Toxostoma"),
                                 "Lowland species", "Highland species")),
    aes(
      shape = ifelse(point_type == "Lowland species", "square", species),
      fill = point_type
    ),
    color = "black",
    size = 3.5,
    stroke = 0.3,
    alpha = 0.9,
    show.legend = c(shape = FALSE)
  ) +
  
  scale_shape_manual(values = c("square" = 22, "Picoides" = 21, "Poecile" = 21, "Chamaea" = 21, "Toxostoma" = 21)) +
  
  scale_fill_manual(
    values = c(
      "Lowland species" = "gray60",
      "Highland species" = "white"
    ),
    name = "Species group",
    guide = guide_legend(
      order = 1,
      override.aes = list(
        shape = c(22, 21),
        fill = c("gray60", "white"),
        color = "black",
        stroke = 0.4
      )
    )
  ) +
  
  coord_sf() +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal"
  )
map_plot


# SAVE MAP
ggsave("Figs/Fig_1_map.pdf",
       map_plot,
       width = 10,
       height = 9,
       dpi = 600)



###done with map









































# ==========================================================
# SET WORKING DIRECTORY
# ==========================================================
setwd("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics")

# ==========================================================
# LOAD PACKAGES
# ==========================================================
library(sf)
library(dplyr)
library(readr)

# ==========================================================
# LOAD ECOREGIONS
# ==========================================================
eco <- st_read("/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared-data/shapefiles/us_eco_l3_state_boundaries/us_eco_l3_state_boundaries.shp",
               quiet = TRUE)

eco_west <- eco %>%
  dplyr::filter(STATE_NAME %in% c("California", "Oregon", "Washington"))

# ==========================================================
# FUNCTION: ADD ONLY ONE COLUMN (ecoregion)
# ==========================================================
add_ecoregion <- function(meta_file, output_file) {
  
  # ---- LOAD META ----
  meta <- readr::read_csv(meta_file, show_col_types = FALSE)
  
  # ---- CONVERT TO SF ----
  pts <- sf::st_as_sf(
    meta,
    coords = c("long", "lat"),
    crs = 4326,
    remove = FALSE
  )
  
  pts <- sf::st_transform(pts, sf::st_crs(eco_west))
  
  # ---- SPATIAL JOIN ----
  pts_join <- sf::st_join(pts, eco_west)
  
  # ---- ADD ONLY ONE COLUMN ----
  out <- pts_join %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(ecoregion = US_L3NAME) %>%
    dplyr::select(dplyr::all_of(names(meta)), ecoregion)
  
  # ---- WRITE NEW FILE ----
  readr::write_csv(out, output_file)
  
  cat("Wrote:", output_file, "\n")
}

# ==========================================================
# RUN FOR EACH SPECIES
# ==========================================================
add_ecoregion(
  "data/Chamaea/Chamaea_meta.csv",
  "data/Chamaea/Chamaea_meta.csv"
)

add_ecoregion(
  "data/Picoides/Picoides_meta.csv",
  "data/Picoides/Picoides_meta.csv"
)

add_ecoregion(
  "data/Poecile/Poecile_meta.csv",
  "data/Poecile/Poecile_meta.csv"
)

add_ecoregion(
  "data/Toxostoma/Toxostoma_meta.csv",
  "data/Toxostoma/Toxostoma_meta.csv"
)





























# ==========================================================
# PCA FUNCTION (FIXED)
# ==========================================================
plot_pca <- function(meta_file, eigenval_file, eigenvec_file, species_name) {
  
  eigvals <- scan(eigenval_file, quiet = TRUE)
  pc_percent <- 100 * eigvals / sum(eigvals)
  
  pca <- read.table(eigenvec_file, header = TRUE)
  
  # 🔥 FIX: USE FULL ID (NO CLEANING)
  pca$sampleID <- pca$IID
  
  meta <- read_csv(meta_file, show_col_types = FALSE)
  
  df <- left_join(pca, meta, by = "sampleID")
  
  p <- ggplot(df, aes(PC1, PC2)) +
    geom_point(
      aes(fill = ecoregion),
      shape = 21,
      color = "black",
      size = 3.5,
      stroke = 0.4,
      alpha = 0.95
    ) +
    scale_fill_manual(values = eco_colors, na.translate = FALSE) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      title = paste("PCA -", species_name),
      x = paste0("PC1 (", round(pc_percent[1], 1), "%)"),
      y = paste0("PC2 (", round(pc_percent[2], 1), "%)")
    )
  
  # 🔥 NEW: DISPLAY IN R PLOT PANEL
  print(p)
  
  ggsave(
    paste0("Figs/Fig_2_PCA_", species_name, ".pdf"),
    p,
    width = 4,
    height = 4,
    dpi = 600
  )
}

# ==========================================================
# RUN PCA
# ==========================================================
plot_pca("data/Chamaea/Chamaea_meta.csv",
         "out/Structure/1_pca/chamaea_pca.eigenval",
         "out/Structure/1_pca/chamaea_pca.eigenvec",
         "Chamaea")

plot_pca("data/Picoides/Picoides_meta.csv",
         "out/Structure/1_pca/picoides_pca.eigenval",
         "out/Structure/1_pca/picoides_pca.eigenvec",
         "Picoides")

plot_pca("data/Poecile/Poecile_meta.csv",
         "out/Structure/1_pca/poecile_pca.eigenval",
         "out/Structure/1_pca/poecile_pca.eigenvec",
         "Poecile")

plot_pca("data/Toxostoma/Toxostoma_meta.csv",
         "out/Structure/1_pca/toxostoma_pca.eigenval",
         "out/Structure/1_pca/toxostoma_pca.eigenvec",
         "Toxostoma")

dev.off()








# ==========================================================
# PCA FUNCTION (FIXED)
# ==========================================================
plot_pca <- function(meta_file, eigenval_file, eigenvec_file, species_name) {
  
  eigvals <- scan(eigenval_file, quiet = TRUE)
  pc_percent <- 100 * eigvals / sum(eigvals)
  
  pca <- read.table(eigenvec_file, header = TRUE)
  
  # 🔥 FIX: USE FULL ID (NO CLEANING)
  pca$sampleID <- pca$IID
  
  meta <- read_csv(meta_file, show_col_types = FALSE)
  
  df <- left_join(pca, meta, by = "sampleID")
  
  p <- ggplot(df, aes(PC1, PC2)) +
    geom_point(
      aes(fill = ecoregion),
      shape = 21,
      color = "black",
      size = 3.5,
      stroke = 0.4,
      alpha = 0.95
    ) +
    scale_fill_manual(values = eco_colors, na.translate = FALSE) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank()   # 🔥 removes ALL gridlines
    ) +
    labs(
      x = paste0("PC1 (", round(pc_percent[1], 1), "%)"),
      y = paste0("PC2 (", round(pc_percent[2], 1), "%)")
    )
  
  # 🔥 DISPLAY IN R PLOT PANEL
  print(p)
  
  ggsave(
    paste0("Figs/Fig_2_PCA_", species_name, ".pdf"),
    p,
    width = 3,
    height = 2,
    dpi = 600
  )
}

# ==========================================================
# RUN PCA
# ==========================================================
plot_pca("data/Chamaea/Chamaea_meta.csv",
         "out/Structure/1_pca/chamaea_pca.eigenval",
         "out/Structure/1_pca/chamaea_pca.eigenvec",
         "Chamaea")

plot_pca("data/Picoides/Picoides_meta.csv",
         "out/Structure/1_pca/picoides_pca.eigenval",
         "out/Structure/1_pca/picoides_pca.eigenvec",
         "Picoides")

plot_pca("data/Poecile/Poecile_meta.csv",
         "out/Structure/1_pca/poecile_pca.eigenval",
         "out/Structure/1_pca/poecile_pca.eigenvec",
         "Poecile")

plot_pca("data/Toxostoma/Toxostoma_meta.csv",
         "out/Structure/1_pca/toxostoma_pca.eigenval",
         "out/Structure/1_pca/toxostoma_pca.eigenvec",
         "Toxostoma")

dev.off()







