# ==========================================================
# LOAD PACKAGES
# ==========================================================
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(readr)
library(sf)
library(RColorBrewer)

# ==========================================================
# FILE PATHS
# ==========================================================
tree_file <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/iqtree_run_bootstrap_4.varsites.phy.treefile'

picoides_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/data/Picoides/Picoides_meta.csv"

chamaea_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/data/Chamaea/Chamaea_meta.csv"
poecile_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/data/Poecile/Poecile_meta.csv"
toxostoma_file <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/data/Toxostoma/Toxostoma_meta.csv"

ecoregions_path <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared-data/shapefiles/us_eco_l3_state_boundaries/us_eco_l3_state_boundaries.shp"

# ==========================================================
# LOAD ECOREGIONS (UNCHANGED)
# ==========================================================
eco <- sf::st_read(ecoregions_path, quiet = TRUE)

eco_west <- eco %>%
  dplyr::filter(STATE_NAME %in% c("California", "Oregon", "Washington"))

# ==========================================================
# BUILD ECOREGION COLOR ORDER (UNCHANGED)
# ==========================================================
get_ecos <- function(file){
  readr::read_csv(file, show_col_types = FALSE) %>%
    dplyr::pull(ecoregion)
}

eco_levels <- sort(unique(c(
  get_ecos(chamaea_file),
  get_ecos(picoides_file),
  get_ecos(poecile_file),
  get_ecos(toxostoma_file)
)))

base_cols <- c(
  RColorBrewer::brewer.pal(8, "Set2"),
  RColorBrewer::brewer.pal(8, "Dark2"),
  RColorBrewer::brewer.pal(8, "Paired")
)

eco_colors <- stats::setNames(
  base_cols[1:length(eco_levels)],
  eco_levels
)

# ==========================================================
# LOAD TREE
# ==========================================================
tree <- ape::read.tree(tree_file)

# ==========================================================
# SAFE ROOTING (OPTIONAL)
# ==========================================================
outgroup_id <- "Paradoxornis"

if(outgroup_id %in% tree$tip.label){
  tree <- ape::root(tree, outgroup = outgroup_id, resolve.root = TRUE)
  tree <- ape::drop.tip(tree, outgroup_id)
}

# ==========================================================
# LOAD META (PICOIDES)
# ==========================================================
meta_raw <- readr::read_csv(chamaea_file, show_col_types = FALSE) %>%
  dplyr::transmute(
    label = sampleID,
    US_L3NAME = ecoregion
  )

# ==========================================================
# MATCH META TO TREE (NO DROPPING)
# ==========================================================
meta <- meta_raw[match(tree$tip.label, meta_raw$label), ]

# ==========================================================
# SANITY CHECKS
# ==========================================================
if(any(is.na(meta$label))){
  missing <- tree$tip.label[is.na(meta$label)]
  stop("Missing samples in metadata:\n", paste(missing, collapse = "\n"))
}

if(any(is.na(meta$US_L3NAME))){
  stop("Some samples have missing ecoregion values")
}

# ==========================================================
# FACTOR LEVELS
# ==========================================================
tree_eco_levels <- eco_levels[eco_levels %in% unique(meta$US_L3NAME)]
meta$US_L3NAME <- factor(meta$US_L3NAME, levels = tree_eco_levels)

# ==========================================================
# BUILD TREE
# ==========================================================
p <- ggtree::ggtree(tree, layout = "rectangular") %<+% meta +
  ggtree::geom_tippoint(aes(color = US_L3NAME), size = 3) +
  ggtree::geom_tiplab(size = 0.001, offset = 0) +
  
  # 🔥 ASTERISK FOR LOW SUPPORT (<95)
  ggtree::geom_text2(
    aes(label = ifelse(as.numeric(label) < 95, "*", "")),
    size = 10,
    hjust = -0.3
  ) +
  
  ggplot2::scale_color_manual(
    values = eco_colors[tree_eco_levels],
    drop = FALSE
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "top",
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 9),
    plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
  ) +
  ggplot2::labs(color = "Ecoregion")

# tighten right edge
p <- p + ggplot2::coord_cartesian(
  xlim = c(0, max(p$data$x, na.rm = TRUE) + 0.02)
)

# ==========================================================
# PLOT
# ==========================================================
p

# ==========================================================
# SAVE
# ==========================================================
ggplot2::ggsave(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/Figs/Fig_2_tree_chamaea.pdf",
  plot = p,
  width = 4,
  height = 4
)
