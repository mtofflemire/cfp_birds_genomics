# ==========================================================
# 1) GLOBAL CLIMATE PCA (ALL FOUR SPECIES)
# ==========================================================

library(terra)
library(sf)
library(dplyr)

project_dir <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs"

# ---- Read coordinates ----
cham <- read.table(file.path(project_dir, "data/chamaea/coords.txt"), header=FALSE)
pico <- read.table(file.path(project_dir, "data/picoides/coords.txt"), header=FALSE)
toxo <- read.table(file.path(project_dir, "data/toxostoma/coords.txt"), header=FALSE)
poec <- read.table(file.path(project_dir, "data/poecile/coords.txt"), header=FALSE)

colnames(cham)[1:2] <- c("lat","lon")
colnames(pico)[1:2] <- c("lat","lon")
colnames(toxo)[1:2] <- c("lat","lon")
colnames(poec)[1:2] <- c("lat","lon")

cham$species <- "chamaea"
pico$species <- "picoides"
toxo$species <- "toxostoma"
poec$species <- "poecile"

coords_all <- rbind(cham, pico, toxo, poec)

coords_all_sf   <- st_as_sf(coords_all, coords=c("lon","lat"), crs=4326)
coords_all_vect <- vect(coords_all_sf)

# ---- Load climate ----
clim_files <- list.files(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared-data/climate/wc2",
  pattern="bio_.*\\.tif$",
  full.names=TRUE
)

climate <- rast(clim_files)

env_vals_all <- terra::extract(climate, coords_all_vect, ID=FALSE)
keep_rows <- complete.cases(env_vals_all)

env_vals_all <- env_vals_all[keep_rows, ]
coords_all   <- coords_all[keep_rows, ]

# ---- Global PCA ----
env_scaled_all <- scale(env_vals_all)

pca_env_global <- prcomp(env_scaled_all,
                         center=FALSE,
                         scale.=FALSE)

pc_scores_all <- as.data.frame(pca_env_global$x[,1:3])
colnames(pc_scores_all) <- c("env_PC1","env_PC2","env_PC3")
pc_scores_all$species <- coords_all$species

# ---- Save for reuse ----
saveRDS(pca_env_global, file.path(project_dir, "global_climate_pca.rds"))
saveRDS(pc_scores_all,  file.path(project_dir, "global_climate_pc_scores.rds"))













# ==========================================================
# 2) SPECIES-SPECIFIC RDA (PICOIDES)
# ==========================================================

library(vcfR)
library(vegan)
library(dplyr)

# ---- Load global PCA objects ----
pca_env_global <- readRDS(file.path(project_dir, "global_climate_pca.rds"))
pc_scores_all  <- readRDS(file.path(project_dir, "global_climate_pc_scores.rds"))

# ---- Extract Picoides environmental PCs ----
env_pc_pico <- pc_scores_all %>%
  filter(species == "toxostoma") %>%
  select(env_PC1, env_PC2, env_PC3)

env_pc_pico <- as.data.frame(env_pc_pico)

# ---- Build genotype matrix ----
setwd('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs')

vcf <- read.vcfR('data/toxostoma/toxostoma.pruned.selective.snps.vcf.gz')
gt  <- extract.gt(vcf, element="GT", as.numeric=TRUE)

gen.imp <- apply(gt, 2, function(x)
  replace(x, is.na(x),
          as.numeric(names(which.max(table(x)))))
)


gen.mat <- t(gen.imp)

# ---- Safety check ----
if(nrow(gen.mat) != nrow(env_pc_pico)){
  stop("Mismatch between genotype rows and env PC rows.")
}

# ---- Run RDA ----
picoides.rda <- rda(gen.mat ~ env_PC1 + env_PC2 + env_PC3,
                    data = env_pc_pico,
                    scale = TRUE)

summary(picoides.rda)
anova.cca(picoides.rda)
anova.cca(picoides.rda, by="axis")
anova.cca(picoides.rda, by="terms")







# ==========================================================
# 3) SNP LOADINGS + OUTLIERS
# ==========================================================

load.rda <- scores(picoides.rda,
                   choices=1:3,
                   display="species",
                   scaling=3)

outliers <- function(x, z){
  lims <- mean(x) + c(-1,1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}

cand1 <- outliers(load.rda[,1], 3)
cand2 <- outliers(load.rda[,2], 3)
cand3 <- outliers(load.rda[,3], 3)

cand1_df <- data.frame(axis=1, snp=names(cand1), loading=cand1)
cand2_df <- data.frame(axis=2, snp=names(cand2), loading=cand2)
cand3_df <- data.frame(axis=3, snp=names(cand3), loading=cand3)

cand <- rbind(cand1_df, cand2_df, cand3_df)
cand <- cand[!duplicated(cand$snp),]

nrow(cand)







foo <- matrix(nrow=nrow(cand),
              ncol=ncol(env_pc_pico))
colnames(foo) <- colnames(env_pc_pico)

for(i in 1:nrow(cand)){
  snp.gen <- gen.mat[, cand$snp[i]]
  foo[i,] <- apply(env_pc_pico, 2,
                   function(x) cor(x, snp.gen))
}

cand <- cbind(cand, foo)

cand$predictor   <- NA
cand$correlation <- NA

for(i in 1:nrow(cand)){
  vals <- abs(cand[i,4:(3+ncol(env_pc_pico))])
  cand$predictor[i]   <- names(which.max(vals))
  cand$correlation[i] <- max(vals)
}

table(cand$predictor)




getwd()

dir.create("out/rda", recursive = TRUE, showWarnings = FALSE)

write.csv(
  cand,
  file = "out/rda/toxostoma_RDA_outlier_snps_full.csv",
  row.names = FALSE
)








# Split on last underscore
outlier_split <- do.call(
  rbind,
  strsplit(cand$snp, "_(?=[^_]+$)", perl = TRUE)
)

outlier_df <- data.frame(
  CHROM = outlier_split[,1],
  POS   = outlier_split[,2],
  stringsAsFactors = FALSE
)

# Sort by chromosome + position
outlier_df <- outlier_df[
  order(outlier_df$CHROM,
        as.numeric(outlier_df$POS)),
]

write.table(
  outlier_df,
  file = "out/rda/toxostoma_RDA_outliers.pos.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

getwd()






all_snps     <- colnames(gen.mat)
neutral_snps <- setdiff(all_snps, cand$snp)

neutral_split <- do.call(
  rbind,
  strsplit(neutral_snps, "_(?=[^_]+$)", perl = TRUE)
)

neutral_df <- data.frame(
  CHROM = neutral_split[,1],
  POS   = neutral_split[,2],
  stringsAsFactors = FALSE
)

neutral_df <- neutral_df[
  order(neutral_df$CHROM,
        as.numeric(neutral_df$POS)),
]

write.table(
  neutral_df,
  file = "out/rda/toxostoma_RDA_neutral.pos.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)






























# ==========================================================
# LOAD LIBRARIES
# ==========================================================

library(vegan)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)

# ==========================================================
# 1) PARTIAL RDA (PC1 | PC2 + PC3 + geography)
# ==========================================================

picoides.rda <- rda(
  gen.mat ~ env_PC1 +
    Condition(env_PC2 + env_PC3 + lon + lat),
  data  = env_geo_toxo,
  scale = TRUE
)

# ==========================================================
# 2) EXTRACT RAW RDA1 SITE SCORES
# ==========================================================

site_scores <- scores(picoides.rda,
                      display = "sites",
                      choices = 1,
                      scaling = 3)

coords_toxo_sf$rda1_partial <- site_scores[,1]

# ==========================================================
# 3) PROJECT GLOBAL ENVIRONMENTAL PC1 RASTER
# ==========================================================

names(climate) <- rownames(pca_env_global$rotation)

pc1_raster <- predict(climate,
                      pca_env_global,
                      index = 1,
                      na.rm = TRUE)

# ==========================================================
# 4) CROP TO CALIFORNIA
# ==========================================================

usa_states <- ne_states(country = "United States of America",
                        returnclass = "sf")

california <- usa_states[usa_states$name == "California", ]

california_vect <- vect(california)

pc1_ca <- crop(pc1_raster, california_vect)
pc1_ca <- mask(pc1_ca, california_vect)

pc1_df <- as.data.frame(pc1_ca, xy = TRUE, na.rm = TRUE)
colnames(pc1_df) <- c("x", "y", "PC1")

# ==========================================================
# 5) FINAL MAP
# ==========================================================

ggplot() +
  
  geom_raster(data = pc1_df,
              aes(x = x, y = y, fill = PC1),
              alpha = 1) +
  
  scale_fill_gradientn(
    colors = terrain.colors(100),
    name = "Env PC1"
  ) +
  
  geom_sf(data = california,
          fill = NA,
          color = "black",
          linewidth = 0.6) +
  
  geom_sf(data = coords_toxo_sf,
          aes(color = rda1_partial),
          size = 3.5) +
  
  scale_color_viridis_c(
    option = "plasma",
    name = "Genetic (RDA1)",
    limits = quantile(coords_toxo_sf$rda1_partial,
                      probs = c(0.05, 0.95),
                      na.rm = TRUE),
    oob = scales::squish
  ) +
  
  coord_sf(expand = FALSE) +
  
  theme_bw(base_size = 10) +
  
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      linewidth = 1
    ),
    
    legend.position = c(1, 1),
    legend.justification = c("right", "top"),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.4
    )
  )




















# ==========================================================
# PARTIAL RDA SETUP FOR PC1, PC2, PC3
# ==========================================================

library(vegan)

# ----------------------------------------------------------
# Make sure this object exists:
# env_geo_toxo must contain:
# env_PC1, env_PC2, env_PC3, lon, lat
# Rows must match gen.mat
# ----------------------------------------------------------

stopifnot(nrow(gen.mat) == nrow(env_geo_toxo))

# ==========================================================
# 1) RUN PARTIAL RDAs
# ==========================================================

# ---- PC1 only ----
rda_pc1 <- rda(
  gen.mat ~ env_PC1 +
    Condition(env_PC2 + env_PC3 + lon + lat),
  data  = env_geo_toxo,
  scale = TRUE
)

# ---- PC2 only ----
rda_pc2 <- rda(
  gen.mat ~ env_PC2 +
    Condition(env_PC1 + env_PC3 + lon + lat),
  data  = env_geo_toxo,
  scale = TRUE
)

# ---- PC3 only ----
rda_pc3 <- rda(
  gen.mat ~ env_PC3 +
    Condition(env_PC1 + env_PC2 + lon + lat),
  data  = env_geo_toxo,
  scale = TRUE
)

# ==========================================================
# 2) TEST SIGNIFICANCE
# ==========================================================

cat("\nPC1 model:\n")
print(anova.cca(rda_pc1))

cat("\nPC2 model:\n")
print(anova.cca(rda_pc2))

cat("\nPC3 model:\n")
print(anova.cca(rda_pc3))

# ==========================================================
# 3) EXTRACT RAW SITE SCORES
# ==========================================================

coords_toxo_sf$rda1_partial <-
  scores(rda_pc1, display = "sites", choices = 1, scaling = 3)[,1]

coords_toxo_sf$rda2_partial <-
  scores(rda_pc2, display = "sites", choices = 1, scaling = 3)[,1]

coords_toxo_sf$rda3_partial <-
  scores(rda_pc3, display = "sites", choices = 1, scaling = 3)[,1]

# ==========================================================
# 4) QUICK DIAGNOSTIC CHECKS
# ==========================================================

summary(coords_toxo_sf$rda1_partial)
summary(coords_toxo_sf$rda2_partial)
summary(coords_toxo_sf$rda3_partial)

sd(coords_toxo_sf$rda1_partial)
sd(coords_toxo_sf$rda2_partial)
sd(coords_toxo_sf$rda3_partial)






# Make sure raster names match PCA loadings
names(climate) <- rownames(pca_env_global$rotation)

# Project PCs
pc1_raster <- predict(climate, pca_env_global, index = 1, na.rm = TRUE)
pc2_raster <- predict(climate, pca_env_global, index = 2, na.rm = TRUE)
pc3_raster <- predict(climate, pca_env_global, index = 3, na.rm = TRUE)







library(rnaturalearth)
library(rnaturalearthdata)

usa_states <- ne_states(country = "United States of America",
                        returnclass = "sf")

california <- usa_states[usa_states$name == "California", ]
california_vect <- vect(california)

# Function to crop + mask
crop_to_ca <- function(r){
  r2 <- crop(r, california_vect)
  r2 <- mask(r2, california_vect)
  df <- as.data.frame(r2, xy = TRUE, na.rm = TRUE)
  colnames(df) <- c("x","y","PC")
  return(df)
}

pc1_df <- crop_to_ca(pc1_raster)
pc2_df <- crop_to_ca(pc2_raster)
pc3_df <- crop_to_ca(pc3_raster)









plot_partial_map <- function(pc_df, genetic_vector, title_text){
  
  ggplot() +
    
    # Environmental surface
    geom_raster(data = pc_df,
                aes(x = x, y = y, fill = PC),
                alpha = 1) +
    
    scale_fill_gradientn(
      colors = terrain.colors(100),
      name = "Env Axis"
    ) +
    
    # California border
    geom_sf(data = california,
            fill = NA,
            color = "black",
            linewidth = 0.6) +
    
    # Genetic points
    geom_sf(data = coords_toxo_sf,
            aes(color = genetic_vector),
            size = 3.5) +
    
    scale_color_viridis_c(
      option = "plasma",
      name = "Genetic Axis",
      limits = quantile(genetic_vector,
                        probs = c(0.05, 0.95),
                        na.rm = TRUE),
      oob = scales::squish
    ) +
    
    coord_sf(expand = FALSE) +
    
    theme_bw(base_size = 8) +
    
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill  = NA,
        linewidth = 1.5
      ),
      legend.position = c(1,1),
      legend.justification = c("right","top"),
      legend.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.4
      )
    ) +
    
    labs(title = title_text)
}






p1 <- plot_partial_map(pc1_df,
                       coords_toxo_sf$rda1_partial,
                       "PC1 (Env | geo removed)")

p2 <- plot_partial_map(pc2_df,
                       coords_toxo_sf$rda2_partial,
                       "PC2 (Env | geo removed)")

p3 <- plot_partial_map(pc3_df,
                       coords_toxo_sf$rda3_partial,
                       "PC3 (Env | geo removed)")



library(patchwork)

p1 / p2 / p3




rda_ibd <- rda(
  gen.mat ~ lon + lat +
    Condition(env_PC1 + env_PC2 + env_PC3),
  data  = env_geo_toxo,
  scale = TRUE
)


anova.cca(rda_ibd)





coords_toxo_sf$rda_ibd <-
  scores(rda_ibd,
         display = "sites",
         choices = 1,
         scaling = 3)[,1]





p_ibd <- ggplot() +
  
  geom_sf(data = california,
          fill = "gray95",
          color = "black",
          linewidth = 0.6) +
  
  geom_sf(data = coords_toxo_sf,
          aes(color = rda_ibd),
          size = 3.5) +
  
  scale_color_viridis_c(
    option = "plasma",
    name = "IBD Axis",
    limits = quantile(coords_toxo_sf$rda_ibd,
                      probs = c(0.05, 0.95),
                      na.rm = TRUE),
    oob = scales::squish
  ) +
  
  coord_sf(expand = FALSE) +
  
  theme_bw(base_size = 8) +
  
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      linewidth = 1.5
    ),
    legend.position = c(1,1),
    legend.justification = c("right","top"),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.4
    )
  ) +
  
  labs(title = "Isolation by Distance (Env Controlled)")


p_ibd


(p1 / p2 / p3 / p_ibd)


































# ==========================================================
# TOXOSTOMA — CLEAN PARTIAL RDA + MAP WORKFLOW
# ==========================================================

library(vcfR)
library(vegan)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(scales)

# ----------------------------------------------------------
# 1) FILTER GLOBAL PCA SCORES FOR TOXOSTOMA
# ----------------------------------------------------------

env_pc_toxo <- pc_scores_all %>%
  dplyr::filter(species == "picoides") %>%
  dplyr::select(env_PC1, env_PC2, env_PC3)

# ----------------------------------------------------------
# 2) FILTER COORDINATES FOR TOXOSTOMA
# ----------------------------------------------------------

coords_toxo_sf <- coords_all_sf[
  coords_all$species == "picoides", ]

# ----------------------------------------------------------
# 3) LOAD OUTLIER VCF (TOXOSTOMA)
# ----------------------------------------------------------

vcf <- read.vcfR(
  "data/toxostoma/picoides.pruned.selective.snps.vcf.gz"
)

gt  <- extract.gt(vcf, element="GT", as.numeric=TRUE)

gen.imp <- apply(gt, 2, function(x)
  replace(x, is.na(x),
          as.numeric(names(which.max(table(x)))))
)

gen.mat <- t(gen.imp)

stopifnot(nrow(gen.mat) == nrow(env_pc_toxo))

# ----------------------------------------------------------
# 4) BUILD ENV + GEO DATAFRAME
# ----------------------------------------------------------

env_geo_toxo <- data.frame(
  env_PC1 = env_pc_toxo$env_PC1,
  env_PC2 = env_pc_toxo$env_PC2,
  env_PC3 = env_pc_toxo$env_PC3,
  lon     = st_coordinates(coords_toxo_sf)[,1],
  lat     = st_coordinates(coords_toxo_sf)[,2]
)

# ----------------------------------------------------------
# 5) PARTIAL RDAs
# ----------------------------------------------------------

rda_pc1 <- rda(
  gen.mat ~ env_PC1 +
    Condition(env_PC2 + env_PC3 + lon + lat),
  data  = env_geo_toxo,
  scale = TRUE
)

rda_pc2 <- rda(
  gen.mat ~ env_PC2 +
    Condition(env_PC1 + env_PC3 + lon + lat),
  data  = env_geo_toxo,
  scale = TRUE
)

rda_pc3 <- rda(
  gen.mat ~ env_PC3 +
    Condition(env_PC1 + env_PC2 + lon + lat),
  data  = env_geo_toxo,
  scale = TRUE
)

rda_ibd <- rda(
  gen.mat ~ lon + lat +
    Condition(env_PC1 + env_PC2 + env_PC3),
  data  = env_geo_toxo,
  scale = TRUE
)

# ----------------------------------------------------------
# 6) EXTRACT SITE SCORES
# ----------------------------------------------------------

coords_toxo_sf$rda1_partial <-
  scores(rda_pc1, display="sites", choices=1, scaling=3)[,1]

coords_toxo_sf$rda2_partial <-
  scores(rda_pc2, display="sites", choices=1, scaling=3)[,1]

coords_toxo_sf$rda3_partial <-
  scores(rda_pc3, display="sites", choices=1, scaling=3)[,1]

coords_toxo_sf$rda_ibd <-
  scores(rda_ibd, display="sites", choices=1, scaling=3)[,1]

# ----------------------------------------------------------
# 7) PROJECT GLOBAL PCA RASTER (PC1 ONLY SHOWN HERE)
# ----------------------------------------------------------

names(climate) <- rownames(pca_env_global$rotation)

pc1_raster <- predict(climate, pca_env_global, index=1)

# ----------------------------------------------------------
# 8) CROP TO CALIFORNIA
# ----------------------------------------------------------

usa_states <- ne_states(country="United States of America",
                        returnclass="sf")

california <- usa_states[usa_states$name=="California", ]
california_vect <- vect(california)

pc1_ca <- crop(pc1_raster, california_vect)
pc1_ca <- mask(pc1_ca, california_vect)

pc1_df <- as.data.frame(pc1_ca, xy=TRUE, na.rm=TRUE)
colnames(pc1_df) <- c("x","y","PC")

library(ggplot2)
library(patchwork)
library(scales)

# ----------------------------------------------------------
# COMMON MAP THEME
# ----------------------------------------------------------

map_theme <- theme_bw(base_size = 8) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      linewidth = 1.5
    ),
    legend.position = c(1,1),
    legend.justification = c("right","top"),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.4
    )
  )

# ----------------------------------------------------------
# IBD PANEL (top)
# ----------------------------------------------------------

p_ibd <- ggplot() +
  geom_sf(data = california,
          fill = "gray95",
          color = "black",
          linewidth = 0.6) +
  geom_sf(data = coords_toxo_sf,
          aes(color = rda_ibd),
          size = 4) +
  scale_color_viridis_c(
    option = "plasma",
    name   = "IBD",
    limits = quantile(coords_toxo_sf$rda_ibd,
                      probs = c(0.05, 0.95),
                      na.rm = TRUE),
    oob = squish
  ) +
  coord_sf(expand = FALSE) +
  map_theme

# ----------------------------------------------------------
# PC1 PANEL
# ----------------------------------------------------------

p1 <- ggplot() +
  geom_raster(data = pc1_df,
              aes(x = x, y = y, fill = PC)) +
  scale_fill_gradientn(
    colors = terrain.colors(100),
    name   = "Env PC1"
  ) +
  geom_sf(data = california,
          fill = NA,
          color = "black",
          linewidth = 0.6) +
  geom_sf(data = coords_toxo_sf,
          aes(color = rda1_partial),
          size = 3.5) +
  scale_color_viridis_c(
    option = "plasma",
    name   = "Genetic",
    limits = quantile(coords_toxo_sf$rda1_partial,
                      probs = c(0.05, 0.95),
                      na.rm = TRUE),
    oob = squish
  ) +
  coord_sf(expand = FALSE) +
  map_theme

# ----------------------------------------------------------
# PC2 PANEL
# ----------------------------------------------------------

p2 <- ggplot() +
  geom_raster(data = pc2_df,
              aes(x = x, y = y, fill = PC)) +
  scale_fill_gradientn(
    colors = terrain.colors(100),
    name   = "Env PC2"
  ) +
  geom_sf(data = california,
          fill = NA,
          color = "black",
          linewidth = 0.6) +
  geom_sf(data = coords_toxo_sf,
          aes(color = rda2_partial),
          size = 3.5) +
  scale_color_viridis_c(
    option = "plasma",
    name   = "Genetic",
    limits = quantile(coords_toxo_sf$rda2_partial,
                      probs = c(0.05, 0.95),
                      na.rm = TRUE),
    oob = squish
  ) +
  coord_sf(expand = FALSE) +
  map_theme

# ----------------------------------------------------------
# PC3 PANEL
# ----------------------------------------------------------

p3 <- ggplot() +
  geom_raster(data = pc3_df,
              aes(x = x, y = y, fill = PC)) +
  scale_fill_gradientn(
    colors = terrain.colors(100),
    name   = "Env PC3"
  ) +
  geom_sf(data = california,
          fill = NA,
          color = "black",
          linewidth = 0.6) +
  geom_sf(data = coords_toxo_sf,
          aes(color = rda3_partial),
          size = 3.5) +
  scale_color_viridis_c(
    option = "plasma",
    name   = "Genetic",
    limits = quantile(coords_toxo_sf$rda3_partial,
                      probs = c(0.05, 0.95),
                      na.rm = TRUE),
    oob = squish
  ) +
  coord_sf(expand = FALSE) +
  map_theme

# ----------------------------------------------------------
# STACK VERTICALLY (IBD on top)
# ----------------------------------------------------------

final_plot <- p_ibd / p1 / p2 / p3

final_plot

# ----------------------------------------------------------
# SAVE OUTPUT
# ----------------------------------------------------------

ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/out/rda_outliers/toxostoma_partial_RDA_maps.png",
  plot     = final_plot,
  width    = 6,
  height   = 16,
  dpi      = 600,
  bg       = "white"
)




































library(vcfR)
library(vegan)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(scales)

# ----------------------------------------------------------
# 1) FILTER GLOBAL PCA SCORES FOR TOXOSTOMA
# ----------------------------------------------------------

env_pc_toxo <- pc_scores_all %>%
  dplyr::filter(species == "poecile") %>%
  dplyr::select(env_PC1, env_PC2, env_PC3)

# ----------------------------------------------------------
# 2) FILTER COORDINATES FOR TOXOSTOMA
# ----------------------------------------------------------

coords_toxo_sf <- coords_all_sf[
  coords_all$species == "poecile", ]

# ----------------------------------------------------------
# 3) LOAD OUTLIER VCF (TOXOSTOMA)
# ----------------------------------------------------------

vcf <- read.vcfR(
  "data/poecile/poecile.pruned.selective.snps.vcf.gz"
)

gt  <- extract.gt(vcf, element="GT", as.numeric=TRUE)

gen.imp <- apply(gt, 2, function(x)
  replace(x, is.na(x),
          as.numeric(names(which.max(table(x)))))
)

gen.mat <- t(gen.imp)

stopifnot(nrow(gen.mat) == nrow(env_pc_toxo))

# ----------------------------------------------------------
# 4) BUILD ENV + GEO DATAFRAME
# ----------------------------------------------------------

env_geo_toxo <- data.frame(
  env_PC1 = env_pc_toxo$env_PC1,
  env_PC2 = env_pc_toxo$env_PC2,
  env_PC3 = env_pc_toxo$env_PC3,
  lon     = st_coordinates(coords_toxo_sf)[,1],
  lat     = st_coordinates(coords_toxo_sf)[,2]
)

# ----------------------------------------------------------
# 5) PARTIAL RDAs
# ----------------------------------------------------------

rda_pc1 <- rda(
  gen.mat ~ env_PC1 +
    Condition(env_PC2 + env_PC3 + lon + lat),
  data  = env_geo_toxo,
  scale = TRUE
)

rda_pc2 <- rda(
  gen.mat ~ env_PC2 +
    Condition(env_PC1 + env_PC3 + lon + lat),
  data  = env_geo_toxo,
  scale = TRUE
)

rda_pc3 <- rda(
  gen.mat ~ env_PC3 +
    Condition(env_PC1 + env_PC2 + lon + lat),
  data  = env_geo_toxo,
  scale = TRUE
)

rda_ibd <- rda(
  gen.mat ~ lon + lat +
    Condition(env_PC1 + env_PC2 + env_PC3),
  data  = env_geo_toxo,
  scale = TRUE
)

# ----------------------------------------------------------
# 6) EXTRACT SITE SCORES
# ----------------------------------------------------------

coords_toxo_sf$rda1_partial <-
  scores(rda_pc1, display="sites", choices=1, scaling=3)[,1]

coords_toxo_sf$rda2_partial <-
  scores(rda_pc2, display="sites", choices=1, scaling=3)[,1]

coords_toxo_sf$rda3_partial <-
  scores(rda_pc3, display="sites", choices=1, scaling=3)[,1]

coords_toxo_sf$rda_ibd <-
  scores(rda_ibd, display="sites", choices=1, scaling=3)[,1]

# ----------------------------------------------------------
# 7) PROJECT GLOBAL PCA RASTERS
# ----------------------------------------------------------

names(climate) <- rownames(pca_env_global$rotation)

pc1_raster <- predict(climate, pca_env_global, index=1)
pc2_raster <- predict(climate, pca_env_global, index=2)
pc3_raster <- predict(climate, pca_env_global, index=3)

# ----------------------------------------------------------
# 8) WEST COAST EXTENT (CA + OR + WA)
# ----------------------------------------------------------

usa_states <- ne_states(
  country = "United States of America",
  returnclass = "sf"
)

west_coast <- usa_states[
  usa_states$name %in% c("California","Oregon","Washington"),
]

west_coast_vect <- vect(west_coast)

# Crop & mask each PC raster
pc1_west <- mask(crop(pc1_raster, west_coast_vect), west_coast_vect)
pc2_west <- mask(crop(pc2_raster, west_coast_vect), west_coast_vect)
pc3_west <- mask(crop(pc3_raster, west_coast_vect), west_coast_vect)

pc1_df <- as.data.frame(pc1_west, xy=TRUE, na.rm=TRUE)
pc2_df <- as.data.frame(pc2_west, xy=TRUE, na.rm=TRUE)
pc3_df <- as.data.frame(pc3_west, xy=TRUE, na.rm=TRUE)

colnames(pc1_df) <- c("x","y","PC")
colnames(pc2_df) <- c("x","y","PC")
colnames(pc3_df) <- c("x","y","PC")

# ----------------------------------------------------------
# COMMON MAP THEME
# ----------------------------------------------------------

map_theme <- theme_bw(base_size = 5) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      linewidth = 1.2
    ),
    legend.position = c(1,1),
    legend.justification = c("right","top"),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.4
    )
  )

# ----------------------------------------------------------
# IBD PANEL
# ----------------------------------------------------------

p_ibd <- ggplot() +
  geom_sf(data = west_coast,
          fill = "gray95",
          color = "black",
          linewidth = 0.6) +
  geom_sf(data = coords_toxo_sf,
          aes(color = rda_ibd),
          size = 3.8) +
  scale_color_viridis_c(
    option = "plasma",
    name   = "IBD",
    limits = quantile(coords_toxo_sf$rda_ibd,
                      probs = c(0.05, 0.95),
                      na.rm = TRUE),
    oob = squish
  ) +
  coord_sf(expand = FALSE) +
  map_theme

# ----------------------------------------------------------
# PC1 PANEL
# ----------------------------------------------------------

p1 <- ggplot() +
  geom_raster(data = pc1_df,
              aes(x = x, y = y, fill = PC)) +
  scale_fill_gradientn(
    colors = terrain.colors(100),
    name   = "Env PC1"
  ) +
  geom_sf(data = west_coast,
          fill = NA,
          color = "black",
          linewidth = 0.6) +
  geom_sf(data = coords_toxo_sf,
          aes(color = rda1_partial),
          size = 3.2) +
  scale_color_viridis_c(
    option = "plasma",
    name   = "Genetic",
    limits = quantile(coords_toxo_sf$rda1_partial,
                      probs = c(0.05, 0.95),
                      na.rm = TRUE),
    oob = squish
  ) +
  coord_sf(expand = FALSE) +
  map_theme

# ----------------------------------------------------------
# PC2 PANEL
# ----------------------------------------------------------

p2 <- ggplot() +
  geom_raster(data = pc2_df,
              aes(x = x, y = y, fill = PC)) +
  scale_fill_gradientn(
    colors = terrain.colors(100),
    name   = "Env PC2"
  ) +
  geom_sf(data = west_coast,
          fill = NA,
          color = "black",
          linewidth = 0.6) +
  geom_sf(data = coords_toxo_sf,
          aes(color = rda2_partial),
          size = 3.2) +
  scale_color_viridis_c(
    option = "plasma",
    name   = "Genetic",
    limits = quantile(coords_toxo_sf$rda2_partial,
                      probs = c(0.05, 0.95),
                      na.rm = TRUE),
    oob = squish
  ) +
  coord_sf(expand = FALSE) +
  map_theme

# ----------------------------------------------------------
# PC3 PANEL
# ----------------------------------------------------------

p3 <- ggplot() +
  geom_raster(data = pc3_df,
              aes(x = x, y = y, fill = PC)) +
  scale_fill_gradientn(
    colors = terrain.colors(100),
    name   = "Env PC3"
  ) +
  geom_sf(data = west_coast,
          fill = NA,
          color = "black",
          linewidth = 0.6) +
  geom_sf(data = coords_toxo_sf,
          aes(color = rda3_partial),
          size = 3.2) +
  scale_color_viridis_c(
    option = "plasma",
    name   = "Genetic",
    limits = quantile(coords_toxo_sf$rda3_partial,
                      probs = c(0.05, 0.95),
                      na.rm = TRUE),
    oob = squish
  ) +
  coord_sf(expand = FALSE) +
  map_theme

# ----------------------------------------------------------
# STACK VERTICALLY
# ----------------------------------------------------------

final_plot <- p_ibd / p1 / p2 / p3

final_plot

# ----------------------------------------------------------
# SAVE OUTPUT
# ----------------------------------------------------------

ggsave(
  filename = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/out/rda_outliers/poecile_partial_RDA_maps_westcoast.png",
  plot     = final_plot,
  width    = 6,
  height   = 16,
  dpi      = 600,
  bg       = "white"
)
