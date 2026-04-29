# ==========================================================
# LOAD LIBRARIES
# ==========================================================
library(tidyverse)
library(stringr)
library(sf)

# ==========================================================
# FILE PATHS (.het)
# ==========================================================
het_files <- c(
  Chamaea   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/41-Chamaea.diversity.het",
  Picoides  = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/picoides.diversity.het",
  Poecile   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/poecile.diversity.autosomes.het",
  Toxostoma = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/Toxostoma_filteredQC_2.het"
)

# ==========================================================
# FILE PATHS (METADATA WITH ECOREGIONS)
# ==========================================================
meta_files <- c(
  Chamaea   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/data/Chamaea/Cfa-metadata-with-ecoregion.csv",
  Picoides  = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/data/Leuconotopicus/picoides_Meta_with_ecoregion.csv",
  Poecile   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/data/MountainChickadee/poecile_metadata_with_ecoregion.csv",
  Toxostoma = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/data/Toxostoma/toxostoma_meta_with_ecoregion.csv"
)

# ==========================================================
# ECOREGION SHAPEFILE (FOR ORDERING)
# ==========================================================
ecoregions_path <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared-data/shapefiles/us_eco_l3_state_boundaries/us_eco_l3_state_boundaries.shp"

eco <- st_read(ecoregions_path, quiet = TRUE)

eco_west <- eco %>%
  filter(STATE_NAME %in% c("California","Oregon","Washington"))

# ==========================================================
# NORTH → SOUTH ORDER
# ==========================================================
eco_order <- eco_west %>%
  group_by(US_L3NAME) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_centroid() %>%
  mutate(lat = st_coordinates(.)[,2]) %>%
  arrange(desc(lat)) %>%
  pull(US_L3NAME)

# ==========================================================
# FUNCTION: READ HET + JOIN ECOREGION
# ==========================================================
read_het_with_eco <- function(het_path, meta_path, species){
  
  het <- read.table(het_path, header = TRUE, check.names = FALSE) %>%
    mutate(
      H_obs = 1 - `O(HOM)` / `N(NM)`,
      sampleID = str_replace(IID, "_S.*$", ""),
      species = species
    ) %>%
    dplyr::select(sampleID, H_obs, species)
  
  meta <- read_csv(meta_path, show_col_types = FALSE) %>%
    dplyr::select(sampleID, US_L3NAME)
  
  left_join(het, meta, by = "sampleID")
}

# ==========================================================
# BUILD DATASET
# ==========================================================
het_all <- purrr::imap_dfr(het_files, function(path, sp){
  read_het_with_eco(path, meta_files[[sp]], sp)
})

het_all <- het_all %>%
  filter(!is.na(US_L3NAME))

# ==========================================================
# SUMMARIZE
# ==========================================================
het_summary <- het_all %>%
  group_by(species, US_L3NAME) %>%
  summarise(
    mean_H = mean(H_obs),
    sd_H   = sd(H_obs),
    n      = n(),
    se_H   = sd_H / sqrt(n),
    .groups = "drop"
  )

# ==========================================================
# APPLY ORDER
# ==========================================================
het_summary$US_L3NAME <- factor(
  het_summary$US_L3NAME,
  levels = eco_order
)

# ==========================================================
# DEFINE POINT GROUPS
# ==========================================================
het_summary <- het_summary %>%
  mutate(
    point_group = case_when(
      species %in% c("Chamaea","Toxostoma") ~ "Lowland",
      TRUE ~ "Highland"
    )
  )

# ==========================================================
# FINAL PLOT
# ==========================================================
ggplot(het_summary, aes(x = US_L3NAME, y = mean_H, group = species)) +
  
  # lines
  geom_line(aes(color = species), size = 1) +
  
  # points (custom shapes + fill)
  geom_point(
    aes(shape = point_group, fill = point_group),
    size = 3,
    color = "black",
    stroke = 0.6
  ) +
  
  # error bars
  geom_errorbar(
    aes(
      ymin = mean_H - se_H,
      ymax = mean_H + se_H,
      color = species
    ),
    width = 0.2,
    size = 0.6
  ) +
  
  # shape mapping
  scale_shape_manual(values = c(
    "Lowland" = 22,   # square
    "Highland" = 21   # circle
  )) +
  
  # fill mapping
  scale_fill_manual(values = c(
    "Lowland" = "gray60",
    "Highland" = "white"
  )) +
  
  theme_bw() +
  
  labs(
    x = "Ecoregion (North → South)",
    y = "Mean genome-wide heterozygosity",
    color = "Species"
  ) +
  
  theme(
    panel.grid.major = element_line(color = "gray85", size = 0.4),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title  = element_text(size = 13),
    
    legend.position = "top"
  )

