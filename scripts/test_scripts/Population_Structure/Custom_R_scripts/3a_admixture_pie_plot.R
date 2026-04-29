# ==========================================================
# LOAD PACKAGES
# ==========================================================
library(conStruct)
library(readr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

sf::sf_use_s2(FALSE)

# ==========================================================
# BASE DIRECTORY + OUTPUT
# ==========================================================
base_dir <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/cfp_comparative_genomics"
out_dir  <- file.path(base_dir, "Figs")

# ==========================================================
# LOAD RANGE SHAPEFILE
# ==========================================================
range_shp <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared-data/shapefiles/SppDataRequest/SppDataRequest.shp"
ranges <- st_read(range_shp, quiet = TRUE)

# ==========================================================
# COUNTRY BORDERS
# ==========================================================
countries <- ne_countries(scale = "medium", returnclass = "sf")
countries <- countries[countries$admin %in% c("United States of America", "Canada", "Mexico"), ]
countries <- st_transform(countries, 4326)
countries <- st_make_valid(countries)

# ==========================================================
# 🔥 STATES FROM NATURAL EARTH (THIS FIXES EVERYTHING)
# ==========================================================
states <- ne_states(country = "united states of america", returnclass = "sf")
states <- st_transform(states, 4326)

west_states <- states[states$name %in% c("California","Oregon","Washington"), ]
west_states <- st_make_valid(west_states)

# dissolve into one clean mask
west_mask <- st_union(west_states)

# ==========================================================
# SPECIES SETUP
# ==========================================================
species_list <- list(
  chamaea = list(
    q = "out/Structure/admixture/chamaea.2.Q",
    fam = "out/Structure/admixture/chamaea.fam",
    meta = "data/Chamaea/Chamaea_meta.csv",
    sci = "Chamaea fasciata"
  ),
  picoides = list(
    q = "out/Structure/admixture/picoides.2.Q",
    fam = "out/Structure/admixture/picoides.fam",
    meta = "data/Picoides/Picoides_meta.csv",
    sci = "Leuconotopicus albolarvatus"
  ),
  poecile = list(
    q = "out/Structure/admixture/poecile.2.Q",
    fam = "out/Structure/admixture/poecile.fam",
    meta = "data/Poecile/Poecile_meta.csv",
    sci = "Poecile gambeli"
  ),
  toxostoma = list(
    q = "out/Structure/admixture/toxostoma.2.Q",
    fam = "out/Structure/admixture/toxostoma.fam",
    meta = "data/Toxostoma/Toxostoma_meta.csv",
    sci = "Toxostoma redivivum"
  )
)

# ==========================================================
# COLORS
# ==========================================================
colors <- c("#1b9e77", "#d95f02")
colors_with_alpha <- adjustcolor(colors, alpha.f = 0.9)

# ==========================================================
# LOOP
# ==========================================================
for (sp in names(species_list)) {
  
  cat("Saving:", sp, "\n")
  
  q_file   <- file.path(base_dir, species_list[[sp]]$q)
  fam_file <- file.path(base_dir, species_list[[sp]]$fam)
  meta_file<- file.path(base_dir, species_list[[sp]]$meta)
  sci_name <- species_list[[sp]]$sci
  
  # ----------------------------------------------------------
  # LOAD ADMIXTURE
  # ----------------------------------------------------------
  admix <- as.matrix(read.table(q_file))
  
  fam <- read.table(fam_file)
  sampleID <- fam$V2
  
  meta <- read_csv(meta_file, show_col_types = FALSE)
  meta <- meta[match(sampleID, meta$sampleID), ]
  
  coords <- as.matrix(meta[, c("long", "lat")])
  
  keep <- complete.cases(coords)
  coords <- coords[keep, ]
  admix  <- admix[keep, ]
  
  # ----------------------------------------------------------
  # RANGE
  # ----------------------------------------------------------
  range_sp <- ranges[grepl(sci_name, ranges$SCI_NAME, ignore.case = TRUE), ]
  range_sp <- st_transform(range_sp, 4326)
  range_sp <- st_make_valid(range_sp)
  
  # 🔥 CLEAN INTERSECTION
  range_sp <- st_intersection(range_sp, west_mask)
  
  # ----------------------------------------------------------
  # COUNTRIES MASKED
  # ----------------------------------------------------------
  countries_mask <- st_intersection(countries, west_mask)
  
  # ----------------------------------------------------------
  # PDF
  # ----------------------------------------------------------
  pdf(file.path(out_dir, paste0("Fig_2_Admixture_", sp, "_K2_map.pdf")),
      width = 6, height = 8)
  
  par(mar = c(0, 0, 0, 0))
  
  # ----------------------------------------------------------
  # BASE WINDOW (SAFE NOW)
  # ----------------------------------------------------------
  bbox <- st_bbox(west_mask)
  
  plot(
    NA,
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    xlab = "",
    ylab = "",
    axes = FALSE,
    asp = 1
  )
  
  # ----------------------------------------------------------
  # STATES
  # ----------------------------------------------------------
  plot(
    st_geometry(west_states),
    col = NA,
    border = "black",
    lwd = 3,
    add = TRUE
  )
  
  # ----------------------------------------------------------
  # COUNTRIES
  # ----------------------------------------------------------
  plot(
    st_geometry(countries_mask),
    col = NA,
    border = "black",
    lwd = 2,
    add = TRUE
  )
  
  # ----------------------------------------------------------
  # RANGE
  # ----------------------------------------------------------
  plot(
    st_geometry(range_sp),
    col = adjustcolor("grey70", alpha.f = 0.5),
    border = "grey40",
    lwd = 1.5,
    add = TRUE
  )
  
  # ----------------------------------------------------------
  # ADMIXTURE
  # ----------------------------------------------------------
  make.admix.pie.plot(
    admix,
    coords,
    colors_with_alpha,
    radii = 3,
    add = TRUE
  )
  
  dev.off()
}

cat("Done — WA/OR/CA mask fully working.\n")

