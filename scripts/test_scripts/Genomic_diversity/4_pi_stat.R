library(tidyverse)

# ===============================
# File paths to per-site π
# ===============================
files <- c(
  Chamaea   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/pi/chamaea.sites.sites.pi",
  Picoides  = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/pi/picoides.sites.sites.pi",
  Poecile   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/pi/poecile.sites.sites.pi",
  Toxostoma = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/pi/toxostoma.sites.sites.pi"
)

# ===============================
# Read and summarize π
# ===============================
read_pi <- function(path, species){
  dat <- read.table(path, header = TRUE)
  
  tibble(
    species = species,
    mean_pi = mean(dat$PI, na.rm = TRUE),
    sd_pi   = sd(dat$PI, na.rm = TRUE)
  )
}

pi_summary <- purrr::imap_dfr(files, read_pi)

# ===============================
# Add clade grouping (same as Hobs)
# ===============================
pi_summary <- pi_summary %>%
  mutate(
    color_group = case_when(
      species %in% c("Chamaea", "Toxostoma") ~ "Wrentit / Thrasher clade",
      species %in% c("Picoides", "Poecile")  ~ "Woodpecker / Chickadee clade"
    )
  )

# Order species
pi_summary$species <- factor(
  pi_summary$species,
  levels = c("Chamaea", "Toxostoma", "Picoides", "Poecile")
)

print(pi_summary)

# ===============================
# Plot styled like heterozygosity
# ===============================
ggplot(pi_summary, aes(x = species, y = mean_pi, fill = color_group)) +
  geom_point(
    shape = 21,        # allows fill + outline
    size = 4,
    color = "black",   # outline color
    stroke = 0.8       # outline thickness
  ) +
  geom_errorbar(aes(
    ymin = mean_pi - sd_pi,
    ymax = mean_pi + sd_pi
  ), width = 0.15, linewidth = 0.3) +
  scale_fill_manual(values = c(
    "Wrentit / Thrasher clade" = "#D55E00",
    "Woodpecker / Chickadee clade" = "#0072B2"
  )) +
  theme_bw() +
  labs(
    x = "Genus",
    y = expression(nucleotide~diversity~(Per-site~pi)),
    fill = "Evolutionary clade"
  ) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "top"
  )
