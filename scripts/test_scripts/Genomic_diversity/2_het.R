# ===============================
# Load libraries
# ===============================
library(tidyverse)
library(stringr)

# ===============================
# File paths
# ===============================
files <- c(
  Chamaea   = '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/41-Chamaea.diversity.het',
  Picoides  = '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/picoides.diversity.het',
  Poecile   = '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/poecile.diversity.autosomes.het',
  Toxostoma = '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/Toxostoma_filteredQC_2.het'
)

# ===============================
# Function to read PLINK .het files
# ===============================
read_het <- function(path, species){
  read.table(path, header = TRUE, check.names = FALSE) %>%
    mutate(
      H_obs = 1 - `O(HOM)` / `N(NM)`,
      sampleID = str_replace(IID, "_S.*$", ""),
      species = species
    ) %>%
    select(sampleID, H_obs, species)
}

# ===============================
# Read and combine all species
# ===============================
het_all <- imap_dfr(files, read_het)

# ===============================
# Add clade color grouping
# ===============================
het_all <- het_all %>%
  mutate(
    color_group = case_when(
      species %in% c("Chamaea", "Toxostoma") ~ "Wrentit / Thrasher clade",
      species %in% c("Picoides", "Poecile")  ~ "Woodpecker / Chickadee clade"
    )
  )

# Order genera for plotting
het_all$species <- factor(
  het_all$species,
  levels = c("Chamaea", "Toxostoma", "Picoides", "Poecile")
)

# ===============================
# Summary table (mean ± SD)
# ===============================
het_summary <- het_all %>%
  group_by(species, color_group) %>%
  summarise(
    mean_H = mean(H_obs),
    sd_H   = sd(H_obs),
    n      = n(),
    .groups = "drop"
  )

print(het_summary)

# ===============================
# Boxplot (genus on x-axis)
# ===============================
ggplot(het_all, aes(x = species, y = H_obs, fill = color_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c(
    "Wrentit / Thrasher clade" = "#D55E00",
    "Woodpecker / Chickadee clade" = "#0072B2"
  )) +
  theme_bw() +
  labs(
    x = "Genus",
    y = "Genome-wide heterozygosity",
    fill = "Evolutionary clade"
  ) +
  theme(
    axis.text = element_text(size = ),
    axis.title = element_text(size = 13),
    legend.position = "top"
  )

