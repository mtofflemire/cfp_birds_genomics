# ==========================================================
# LOAD LIBRARIES
# ==========================================================
library(tidyverse)

# ==========================================================
# FILE PATHS
# ==========================================================
files <- c(
  Chamaea   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/pi/chamaea.sites.sites.pi",
  Picoides  = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/pi/picoides.sites.sites.pi",
  Poecile   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/pi/poecile.sites.sites.pi",
  Toxostoma = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/pi/toxostoma.sites.sites.pi"
)

# ==========================================================
# SPECIES COLORS (EXACT GGPlot DEFAULTS YOU USED BEFORE)
# ==========================================================
species_colors <- c(
  "Chamaea"   = "#F8766D",
  "Toxostoma" = "#7CAE00",
  "Picoides"  = "#00BFC4",
  "Poecile"   = "#C77CFF"
)

# ==========================================================
# FUNCTION TO READ π FILES
# ==========================================================
read_pi <- function(path, species){
  
  dat <- read.table(path, header = TRUE)
  
  tibble(
    species = species,
    mean_pi = mean(dat$PI, na.rm = TRUE),
    sd_pi   = sd(dat$PI, na.rm = TRUE)
  )
}

# ==========================================================
# COMBINE ALL SPECIES
# ==========================================================
pi_summary <- purrr::imap_dfr(files, read_pi)

# ==========================================================
# ORDER SPECIES
# ==========================================================
pi_summary$species <- factor(
  pi_summary$species,
  levels = c("Chamaea", "Toxostoma", "Picoides", "Poecile")
)

# ==========================================================
# FINAL PLOT
# ==========================================================
ggplot(pi_summary, aes(x = species, y = mean_pi)) +
  
  geom_point(
    aes(fill = species, color = species),
    shape = 21,
    size = 4,
    stroke = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = mean_pi - sd_pi,
      ymax = mean_pi + sd_pi,
      color = species
    ),
    width = 0.15,
    linewidth = 0.4
  ) +
  
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  
  theme_bw() +
  
  labs(
    x = "Genus",
    y = expression(nucleotide~diversity~(Per-site~pi))
  ) +
  
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "none"
  )












# ==========================================================
# LOAD LIBRARIES
# ==========================================================
library(tidyverse)
library(stringr)

# ==========================================================
# FILE PATHS
# ==========================================================
files <- c(
  Chamaea   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/41-Chamaea.diversity.het",
  Picoides  = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/picoides.diversity.het",
  Poecile   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/poecile.diversity.autosomes.het",
  Toxostoma = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/Toxostoma_filteredQC_2.het"
)

# ==========================================================
# SPECIES COLORS
# ==========================================================
species_colors <- c(
  "Chamaea"   = "#F8766D",
  "Toxostoma" = "#7CAE00",
  "Picoides"  = "#00BFC4",
  "Poecile"   = "#C77CFF"
)

# ==========================================================
# READ HET FILES
# ==========================================================
read_het <- function(path, species){
  read.table(path, header = TRUE, check.names = FALSE) %>%
    mutate(
      H_obs = 1 - `O(HOM)` / `N(NM)`,
      sampleID = str_replace(IID, "_S.*$", ""),
      species = species
    ) %>%
    dplyr::select(sampleID, H_obs, species)
}

# ==========================================================
# COMBINE DATA
# ==========================================================
het_all <- purrr::imap_dfr(files, read_het)

# ==========================================================
# ORDER SPECIES
# ==========================================================
het_all$species <- factor(
  het_all$species,
  levels = c("Chamaea", "Toxostoma", "Picoides", "Poecile")
)

# ==========================================================
# PLOT
# ==========================================================
ggplot(het_all, aes(x = species, y = H_obs, fill = species)) +
  
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  
  geom_jitter(width = 0.12, size = 2, alpha = 0.6) +
  
  scale_fill_manual(values = species_colors) +
  
  theme_bw() +
  
  labs(
    x = "Genus",
    y = "Genome-wide heterozygosity"
  ) +
  
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "none"
  )













# ==========================================================
# LOAD LIBRARIES
# ==========================================================
library(tidyverse)

# ==========================================================
# FILE PATHS
# ==========================================================
files <- c(
  Chamaea   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/tajima/chamaea_tajima_100kb.Tajima.D",
  Picoides  = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/tajima/picoides_tajima_100kb.Tajima.D",
  Poecile   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/tajima/poecile_tajima_100kb.Tajima.D",
  Toxostoma = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/tajima/toxostoma_tajima_100kb.Tajima.D"
)

# ==========================================================
# SPECIES COLORS
# ==========================================================
species_colors <- c(
  "Chamaea"   = "#F8766D",
  "Toxostoma" = "#7CAE00",
  "Picoides"  = "#00BFC4",
  "Poecile"   = "#C77CFF"
)

# ==========================================================
# READ TAJIMA FILES
# ==========================================================
read_tajima <- function(path, species){
  read.table(path, header = TRUE) %>%
    mutate(species = species) %>%
    dplyr::select(species, TajimaD)
}

# ==========================================================
# COMBINE DATA
# ==========================================================
taj_all <- purrr::imap_dfr(files, read_tajima)

# ==========================================================
# ORDER SPECIES
# ==========================================================
taj_all$species <- factor(
  taj_all$species,
  levels = c("Chamaea", "Toxostoma", "Picoides", "Poecile")
)

# ==========================================================
# PLOT
# ==========================================================
ggplot(taj_all, aes(x = species, y = TajimaD, fill = species)) +
  
  geom_boxplot(alpha = 0.8) +
  
  scale_fill_manual(values = species_colors) +
  
  theme_bw() +
  
  labs(
    x = "Genus",
    y = "Tajima's D (100 kb windows)"
  ) +
  
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "none"
  )
























# ==========================================================
# LOAD LIBRARIES
# ==========================================================
library(tidyverse)
library(stringr)
library(patchwork)

# ==========================================================
# SPECIES COLORS
# ==========================================================
species_colors <- c(
  "Chamaea"   = "#F8766D",
  "Toxostoma" = "#C77CFF",
  "Picoides"  = "#7CAE00",
  "Poecile"   = "#00BFC4"
)

# ==========================================================
# -------- π DATA --------
# ==========================================================
pi_files <- c(
  Chamaea   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/pi/chamaea.sites.sites.pi",
  Picoides  = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/pi/picoides.sites.sites.pi",
  Poecile   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/pi/poecile.sites.sites.pi",
  Toxostoma = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/pi/toxostoma.sites.sites.pi"
)

read_pi <- function(path, species){
  dat <- read.table(path, header = TRUE)
  tibble(
    species = species,
    mean_pi = mean(dat$PI, na.rm = TRUE),
    sd_pi   = sd(dat$PI, na.rm = TRUE)
  )
}

pi_summary <- purrr::imap_dfr(pi_files, read_pi)

pi_summary$species <- factor(pi_summary$species,
                             levels = c("Chamaea","Toxostoma","Picoides","Poecile"))

p_pi <- ggplot(pi_summary, aes(x = species, y = mean_pi)) +
  geom_point(aes(fill = species, color = species),
             shape = 21, size = 4, stroke = 0.8) +
  geom_errorbar(aes(ymin = mean_pi - sd_pi,
                    ymax = mean_pi + sd_pi,
                    color = species),
                width = 0.15, linewidth = 0.4) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  theme_bw() +
  labs(y = expression(nucleotide~diversity~(pi)), x = NULL) +
  theme(legend.position = "none")

# ==========================================================
# -------- HET DATA --------
# ==========================================================
het_files <- c(
  Chamaea   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/41-Chamaea.diversity.het",
  Picoides  = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/picoides.diversity.het",
  Poecile   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/poecile.diversity.autosomes.het",
  Toxostoma = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/het/Toxostoma_filteredQC_2.het"
)

read_het <- function(path, species){
  read.table(path, header = TRUE, check.names = FALSE) %>%
    mutate(
      H_obs = 1 - `O(HOM)` / `N(NM)`,
      sampleID = str_replace(IID, "_S.*$", ""),
      species = species
    ) %>%
    dplyr::select(sampleID, H_obs, species)   # 🔥 FIXED
}

het_all <- purrr::imap_dfr(het_files, read_het)

het_all$species <- factor(het_all$species,
                          levels = c("Chamaea","Toxostoma","Picoides","Poecile"))

p_het <- ggplot(het_all, aes(x = species, y = H_obs, fill = species)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = species_colors) +
  theme_bw() +
  labs(y = "Heterozygosity", x = NULL) +
  theme(legend.position = "none")

# ==========================================================
# -------- TAJIMA DATA --------
# ==========================================================
taj_files <- c(
  Chamaea   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/tajima/chamaea_tajima_100kb.Tajima.D",
  Picoides  = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/tajima/picoides_tajima_100kb.Tajima.D",
  Poecile   = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/tajima/poecile_tajima_100kb.Tajima.D",
  Toxostoma = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics/out/tajima/toxostoma_tajima_100kb.Tajima.D"
)

read_tajima <- function(path, species){
  read.table(path, header = TRUE) %>%
    mutate(species = species) %>%
    dplyr::select(species, TajimaD)   # 🔥 FIXED
}

taj_all <- purrr::imap_dfr(taj_files, read_tajima)

taj_all$species <- factor(taj_all$species,
                          levels = c("Chamaea","Toxostoma","Picoides","Poecile"))

p_taj <- ggplot(taj_all, aes(x = species, y = TajimaD, fill = species)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(values = species_colors) +
  theme_bw() +
  labs(y = "Tajima's D", x = NULL) +
  theme(legend.position = "none")

# ==========================================================
# -------- COMBINE --------
# ==========================================================
final_plot <-  p_het | p_pi | p_taj
 
final_plot

