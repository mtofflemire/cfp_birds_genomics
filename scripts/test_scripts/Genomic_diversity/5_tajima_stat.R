# ===============================
# Load libraries
# ===============================
library(tidyverse)

# ===============================
# File paths to Tajima's D files
# ===============================
files <- c(
  Chamaea   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/tajima/chamaea_tajima_100kb.Tajima.D",
  Picoides  = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/tajima/picoides_tajima_100kb.Tajima.D",
  Poecile   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/tajima/poecile_tajima_100kb.Tajima.D",
  Toxostoma = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/tajima/toxostoma_tajima_100kb.Tajima.D"
)

# ===============================
# Function to read Tajima's D
# ===============================
read_tajima <- function(path, species){
  read.table(path, header = TRUE) %>%
    mutate(species = species) %>%
    select(species, TajimaD)
}

# ===============================
# Combine all species
# ===============================
taj_all <- purrr::imap_dfr(files, read_tajima)

# ===============================
# Add clade grouping (same as Hobs)
# ===============================
taj_all <- taj_all %>%
  mutate(
    color_group = case_when(
      species %in% c("Chamaea", "Toxostoma") ~ "Wrentit / Thrasher clade",
      species %in% c("Picoides", "Poecile")  ~ "Woodpecker / Chickadee clade"
    )
  )

# Order species for plotting
taj_all$species <- factor(
  taj_all$species,
  levels = c("Chamaea", "Toxostoma", "Picoides", "Poecile")
)

# ===============================
# Summary table (mean ± SD)
# ===============================
taj_summary <- taj_all %>%
  group_by(species, color_group) %>%
  summarise(
    mean_D = mean(TajimaD, na.rm = TRUE),
    sd_D   = sd(TajimaD, na.rm = TRUE),
    n      = n(),
    .groups = "drop"
  )

print(taj_summary)

# ===============================
# Boxplot styled like heterozygosity
# ===============================
ggplot(taj_all, aes(x = species, y = TajimaD, fill = color_group)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = c(
    "Wrentit / Thrasher clade" = "#D55E00",
    "Woodpecker / Chickadee clade" = "#0072B2"
  )) +
  theme_bw() +
  labs(
    y = "Tajima's D (10 kb windows)",
    fill = "Evolutionary clade"
  ) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "top"
  )
