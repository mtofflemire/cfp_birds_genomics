library(terra)
library(sf)
library(dplyr)

# ==========================================================
# 0) DEFINE PROJECT ROOT (YOUR ACTUAL PATH)
# ==========================================================

project_dir <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs"


# ==========================================================
# 1) READ ALL FOUR SPECIES COORDINATES
# ==========================================================

cham <- read.table(file.path(project_dir, "data/chamaea/coords.txt"), header = FALSE)
pico <- read.table(file.path(project_dir, "data/picoides/coords.txt"), header = FALSE)
toxo <- read.table(file.path(project_dir, "data/toxostoma/coords.txt"), header = FALSE)
poec <- read.table(file.path(project_dir, "data/poecile/coords.txt"), header = FALSE)

colnames(cham)[1:2] <- c("lat","lon")
colnames(pico)[1:2] <- c("lat","lon")
colnames(toxo)[1:2] <- c("lat","lon")
colnames(poec)[1:2] <- c("lat","lon")


pico
cham$species <- "chamaea"
pico$species <- "picoides"
toxo$species <- "toxostoma"
poec$species <- "poecile"

coords_all <- rbind(cham, pico, toxo, poec)

coords_all_sf <- st_as_sf(coords_all, coords = c("lon","lat"), crs = 4326)
coords_all_vect <- vect(coords_all_sf)





library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

library(ggplot2)
library(rnaturalearth)

usa <- ne_states(country = "United States of America",
                 returnclass = "sf")

west_coast <- usa[usa$name %in% c("California",
                                  "Washington",
                                  "Oregon"), ]

ggplot() +
  geom_sf(data = west_coast,
          fill = "gray95",
          color = "black") +
  geom_sf(data = coords_all_sf,
          aes(color = species),
          size = 2,
          alpha = 0.8) +
  coord_sf(xlim = c(-125, -116),
           ylim = c(32, 49),
           expand = FALSE) +
  theme_minimal(base_size = 12) +
  labs(title = "Sampling Localities – CFP Species",
       color = "Species")



# ==========================================================
# 2) LOAD WORLDCLIM BIO LAYERS (YOUR SHARED-DATA PATH)
# ==========================================================

climate_files <- list.files(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared-data/climate/wc2",
  pattern = "bio_.*\\.tif$",
  full.names = TRUE
)

climate <- rast(climate_files)


# ==========================================================
# 3) EXTRACT CLIMATE FOR ALL LOCALITIES
# ==========================================================

env_vals_all <- terra::extract(climate, coords_all_vect, ID = FALSE)

keep_rows <- complete.cases(env_vals_all)
env_vals_all <- env_vals_all[keep_rows, ]
coords_all   <- coords_all[keep_rows, ]


# ==========================================================
# 4) GLOBAL PCA
# ==========================================================

env_scaled_all <- scale(env_vals_all)

pca_env_global <- prcomp(env_scaled_all,
                         center = FALSE,
                         scale. = FALSE)

var_exp <- 100 * pca_env_global$sdev^2 /
  sum(pca_env_global$sdev^2)

print(round(var_exp[1:5],2))


# ==========================================================
# 5) EXTRACT PC SCORES
# ==========================================================

pc_scores_all <- as.data.frame(pca_env_global$x[,1:3])
colnames(pc_scores_all) <- c("PC1","PC2","PC3")

pc_scores_all$species <- coords_all$species


# ==========================================================
# 6) SPLIT BY SPECIES
# ==========================================================

pc_cham <- pc_scores_all %>% filter(species == "chamaea") %>% select(PC1,PC2,PC3)
pc_pico <- pc_scores_all %>% filter(species == "picoides") %>% select(PC1,PC2,PC3)
pc_toxo <- pc_scores_all %>% filter(species == "toxostoma") %>% select(PC1,PC2,PC3)
pc_poec <- pc_scores_all %>% filter(species == "poecile") %>% select(PC1,PC2,PC3)


# ==========================================================
# 7) ENVIRONMENTAL DISTANCE MATRICES
# ==========================================================
round(pca_env_global$rotation[,1:3], 3)

library(ggplot2)
library(ggrepel)

# Extract loadings
loadings_12 <- as.data.frame(pca_env_global$rotation[, c(1,2)])

# Scale by eigenvalues (important for proper correlation circle)
loadings_12$PC1 <- loadings_12$PC1 * pca_env_global$sdev[1]
loadings_12$PC2 <- loadings_12$PC2 * pca_env_global$sdev[2]

loadings_12$BIO <- rownames(loadings_12)

# Contribution magnitude
loadings_12$contrib <- sqrt(loadings_12$PC1^2 + loadings_12$PC2^2)

# Unit circle
circle <- data.frame(
  x = cos(seq(0, 2*pi, length.out = 200)),
  y = sin(seq(0, 2*pi, length.out = 200))
)

p12 <- ggplot() +
  geom_path(data = circle, aes(x, y),
            color = "grey70", linewidth = 1) +
  geom_hline(yintercept = 0, linewidth = 1) +
  geom_vline(xintercept = 0, linewidth = 1) +
  geom_segment(
    data = loadings_12,
    aes(x = 0, y = 0, xend = PC1, yend = PC2, color = contrib),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 1
  ) +
  geom_text_repel(
    data = loadings_12,
    aes(PC1, PC2, label = BIO, color = contrib),
    size = 4,
    fontface = "bold"
  ) +
  scale_color_gradientn(
    colors = c("#2CB1A1", "#F1C40F", "red2"),
    name = "Contribution"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Global Climate PCA",
    x = xlab,
    y = ylab
  )



ggsave(
  filename = "global_climate_PCA_PC1_PC2.png",
  plot = last_plot(),
  width = 7,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)








loadings_13 <- as.data.frame(pca_env_global$rotation[, c(1,3)])

loadings_13$PC1 <- loadings_13$PC1 * pca_env_global$sdev[1]
loadings_13$PC3 <- loadings_13$PC3 * pca_env_global$sdev[3]

loadings_13$BIO <- rownames(loadings_13)

loadings_13$contrib <- sqrt(loadings_13$PC1^2 + loadings_13$PC3^2)

p13 <- ggplot() +
  geom_path(data = circle, aes(x, y),
            color = "grey70", linewidth = 1) +
  geom_hline(yintercept = 0, linewidth = 1) +
  geom_vline(xintercept = 0, linewidth = 1) +
  
  geom_segment(
    data = loadings_13,
    aes(x = 0, y = 0, xend = PC1, yend = PC3, color = contrib),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 1
  ) +
  
  geom_text_repel(
    data = loadings_13,
    aes(PC1, PC3, label = BIO, color = contrib),
    size = 4,
    fontface = "bold"
  ) +
  
  scale_color_gradientn(
    colors = c("#2CB1A1", "#F1C40F", "red2"),
    name = "Contribution"
  ) +
  
  coord_fixed() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Global Climate PCA",
    x = paste0("PC1 (", round(var_exp[1],1), "%)"),
    y = paste0("PC3 (", round(var_exp[3],1), "%)")
  )


ggsave("global_climate_PCA_PC1_PC3.png",
       plot = p13,
       width = 7,
       height = 7,
       units = "in",
       dpi = 600,
       bg = "white")








# ==========================================================
# 7) ENVIRONMENTAL DISTANCE MATRICES (SEPARATE AXES)
# ==========================================================

# ----- CHAMAEA -----
pc1_cham <- as.matrix(dist(pc_cham$PC1))
pc2_cham <- as.matrix(dist(pc_cham$PC2))
pc3_cham <- as.matrix(dist(pc_cham$PC3))

# ----- PICOIDES -----
pc1_pico <- as.matrix(dist(pc_pico$PC1))
pc2_pico <- as.matrix(dist(pc_pico$PC2))
pc3_pico <- as.matrix(dist(pc_pico$PC3))

# ----- TOXOSTOMA -----
pc1_toxo <- as.matrix(dist(pc_toxo$PC1))
pc2_toxo <- as.matrix(dist(pc_toxo$PC2))
pc3_toxo <- as.matrix(dist(pc_toxo$PC3))

# ----- POECILE -----
pc1_poec <- as.matrix(dist(pc_poec$PC1))
pc2_poec <- as.matrix(dist(pc_poec$PC2))
pc3_poec <- as.matrix(dist(pc_poec$PC3))














coords_cham_sf <- coords_all_sf[coords_all$species == "chamaea", ]
coords_pico_sf <- coords_all_sf[coords_all$species == "picoides", ]
coords_toxo_sf <- coords_all_sf[coords_all$species == "toxostoma", ]
coords_poec_sf <- coords_all_sf[coords_all$species == "poecile", ]

geo_cham <- as.matrix(dist(st_coordinates(coords_cham_sf)))
geo_pico <- as.matrix(dist(st_coordinates(coords_pico_sf)))
geo_toxo <- as.matrix(dist(st_coordinates(coords_toxo_sf)))
geo_poec <- as.matrix(dist(st_coordinates(coords_poec_sf)))


geo_cham



gen_cham <- as.matrix(read.table(
  file.path(project_dir, "data/chamaea/mmrr/chamaea.gendist.neutral.snp.txt")
))

gen_pico <- as.matrix(read.table(
  file.path(project_dir, "data/picoides/mmrr/picoides.gendist.selective.snps.txt")
))


gen_toxo <- as.matrix(read.table(
  file.path(project_dir, "data/toxostoma/mmrr/toxostoma.gendist.selective.snp.txt")
))

gen_poec <- as.matrix(read.table(
  file.path(project_dir, "data/poecile/mmrr/poecile.gendist.selective.snp.txt")
))




dim(gen_cham); dim(geo_cham); dim(pc1_cham)
dim(gen_pico); dim(geo_pico); dim(pc1_pico)
dim(gen_toxo); dim(geo_toxo); dim(pc1_toxo)
dim(gen_poec); dim(geo_poec); dim(pc1_poec)









unfold <- function(X){
  x <- vector()
  for(i in 2:nrow(X)) x <- c(x, X[i,1:i-1])
  scale(x, center=TRUE, scale=TRUE)
}

MMRR <- function(Y, X, nperm=999){
  
  nrowsY <- nrow(Y)
  y <- unfold(Y)
  
  Xmats <- sapply(X, unfold)
  fit <- lm(y ~ Xmats)
  summ <- summary(fit)
  
  coeffs <- fit$coefficients
  r.squared <- summ$r.squared
  tstat <- summ$coefficients[,"t value"]
  
  tprob <- rep(1,length(tstat))
  
  for(i in 1:nperm){
    rand <- sample(1:nrowsY)
    Yperm <- Y[rand,rand]
    yperm <- unfold(Yperm)
    fitp <- lm(yperm ~ Xmats)
    summp <- summary(fitp)
    tprob <- tprob +
      as.numeric(abs(summp$coefficients[,"t value"]) >= abs(tstat))
  }
  
  tp <- tprob/(nperm+1)
  
  names(coeffs) <- c("Intercept", names(X))
  names(tp) <- paste0(c("Intercept",names(X))," (p)")
  
  list(r.squared=r.squared,
       coefficients=coeffs,
       tpvalue=tp)
}






run_mmrr <- function(gen, geo, pc1, pc2, pc3){
  X <- list(
    geodist = geo,
    PC1 = pc1,
    PC2 = pc2,
    PC3 = pc3
  )
  set.seed(12)
  MMRR(gen, X, nperm=999)
}

mmrr_cham <- run_mmrr(gen_cham, ibr_cham,
                      pc1_cham, pc2_cham, pc3_cham)

mmrr_pico <- run_mmrr(gen_pico, geo_pico,
                      pc1_pico, pc2_pico, pc3_pico)

mmrr_toxo <- run_mmrr(gen_toxo, geo_toxo,
                      pc1_toxo, pc2_toxo, pc3_toxo)

mmrr_poec <- run_mmrr(gen_poec, geo_poec,
                      pc1_poec, pc2_poec, pc3_poec)










extract_results <- function(obj, sp){
  data.frame(
    Species = sp,
    R2 = obj$r.squared,
    beta_geo = obj$coefficients["geodist"],
    p_geo = obj$tpvalue["geodist (p)"],
    beta_PC1 = obj$coefficients["PC1"],
    p_PC1 = obj$tpvalue["PC1 (p)"],
    beta_PC2 = obj$coefficients["PC2"],
    p_PC2 = obj$tpvalue["PC2 (p)"],
    beta_PC3 = obj$coefficients["PC3"],
    p_PC3 = obj$tpvalue["PC3 (p)"]
  )
}

results_table <- rbind(
  extract_results(mmrr_cham, "chamaea"),
  extract_results(mmrr_pico, "picoides"),
  extract_results(mmrr_toxo, "toxostoma"),
  extract_results(mmrr_poec, "poecile")
)

results_table



# Round values nicely
results_clean <- results_table

results_clean$R2        <- round(results_clean$R2, 3)
results_clean$beta_geo  <- round(results_clean$beta_geo, 3)
results_clean$beta_PC1  <- round(results_clean$beta_PC1, 3)
results_clean$beta_PC2  <- round(results_clean$beta_PC2, 3)
results_clean$beta_PC3  <- round(results_clean$beta_PC3, 3)

# Print as tab-separated
write.table(
  results_clean,
  file = "",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)










plot_partial <- function(genMat, ibr_cham,
                         pc1_mat, pc2_mat, pc3_mat,
                         mmrr_obj,
                         effect = c("geo","PC1","PC2","PC3")) {
  
  effect <- match.arg(effect)
  
  gen_s <- unfold(genMat)
  geo_s <- unfold(ibr_cham)
  pc1_s <- unfold(pc1_mat)
  pc2_s <- unfold(pc2_mat)
  pc3_s <- unfold(pc3_mat)
  
  if(effect == "geo") {
    gen_resid <- resid(lm(gen_s ~ pc1_s + pc2_s + pc3_s))
    x_resid   <- resid(lm(geo_s ~ pc1_s + pc2_s + pc3_s))
    pval <- mmrr_obj$tpvalue["geodist (p)"]
    xlab <- "IBD (Geo | PCs removed)"
  }
  
  if(effect == "PC1") {
    gen_resid <- resid(lm(gen_s ~ geo_s + pc2_s + pc3_s))
    x_resid   <- resid(lm(pc1_s ~ geo_s + pc2_s + pc3_s))
    pval <- mmrr_obj$tpvalue["PC1 (p)"]
    xlab <- "IBE PC1"
  }
  
  if(effect == "PC2") {
    gen_resid <- resid(lm(gen_s ~ geo_s + pc1_s + pc3_s))
    x_resid   <- resid(lm(pc2_s ~ geo_s + pc1_s + pc3_s))
    pval <- mmrr_obj$tpvalue["PC2 (p)"]
    xlab <- "IBE PC2"
  }
  
  if(effect == "PC3") {
    gen_resid <- resid(lm(gen_s ~ geo_s + pc1_s + pc2_s))
    x_resid   <- resid(lm(pc3_s ~ geo_s + pc1_s + pc2_s))
    pval <- mmrr_obj$tpvalue["PC3 (p)"]
    xlab <- "IBE PC3"
  }
  
  df <- data.frame(Genetic = gen_resid,
                   Predictor = x_resid)
  
  fit <- lm(Genetic ~ Predictor, data = df)
  
  beta <- round(coef(fit)[2], 3)
  r2   <- round(summary(fit)$r.squared, 3)
  pval <- signif(pval, 3)
  
  ggplot(df, aes(Predictor, Genetic)) +
    geom_point(shape=21, fill="gray70", color="gray70",
               size=2, stroke=0) +
    geom_smooth(method="lm", color="black") +
    annotate("text",
             x = Inf, y = Inf,
             hjust = 1.1, vjust = 1.5,
             size = 4.5,
             label = paste0(
               "beta == ", beta,
               "*','~~R^2 == ", r2,
               "*','~~p == ", pval
             ),
             parse = TRUE) +
    theme_bw(base_size = 12) +
    labs(
      x = xlab,
      y = "Genetic Distance (others removed)"
    )
}





# ======================================
# PARTIAL REGRESSION PLOTS — PICOIDES
# ======================================
p_geo  <- plot_partial(gen_cham, ibr_cham,
                       pc1_cham, pc2_cham, pc3_cham,
                       mmrr_cham,
                       effect = "geo")
p_geo
p_pc1  <- plot_partial(gen_cham, geo_cham,
                       pc1_cham, pc2_cham, pc3_cham,
                       mmrr_cham,
                       effect = "PC1")

p_pc2  <- plot_partial(gen_cham, geo_cham,
                       pc1_cham, pc2_cham, pc3_cham,
                       mmrr_cham,
                       effect = "PC2")

p_pc3  <- plot_partial(gen_cham, geo_cham,
                       pc1_cham, pc2_cham, pc3_cham,
                       mmrr_cham,
                       effect = "PC3")






p_geo  <- plot_partial(gen_pico, geo_pico,
                       pc1_pico, pc2_pico, pc3_pico,
                       mmrr_pico,
                       effect = "geo")

p_pc1  <- plot_partial(gen_pico, geo_pico,
                       pc1_pico, pc2_pico, pc3_pico,
                       mmrr_pico,
                       effect = "PC1")

p_pc2  <- plot_partial(gen_pico, geo_pico,
                       pc1_pico, pc2_pico, pc3_pico,
                       mmrr_pico,
                       effect = "PC2")

p_pc3  <- plot_partial(gen_pico, geo_pico,
                       pc1_pico, pc2_pico, pc3_pico,
                       mmrr_pico,
                       effect = "PC3")




# ======================================
# PARTIAL REGRESSION PLOTS — POECILE
# ======================================

p_geo  <- plot_partial(gen_toxo, geo_toxo,
                       pc1_toxo, pc2_toxo, pc3_toxo,
                       mmrr_toxo,
                       effect = "geo")

p_pc1  <- plot_partial(gen_poec, geo_poec,
                       pc1_poec, pc2_poec, pc3_poec,
                       mmrr_poec,
                       effect = "PC1")

p_pc2  <- plot_partial(gen_poec, geo_poec,
                       pc1_poec, pc2_poec, pc3_poec,
                       mmrr_poec,
                       effect = "PC2")

p_pc3  <- plot_partial(gen_poec, geo_poec,
                       pc1_poec, pc2_poec, pc3_poec,
                       mmrr_poec,
                       effect = "PC3")







# ======================================
# PARTIAL REGRESSION PLOTS — TOXOSTOMA
# ======================================

p_geo  <- plot_partial(gen_toxo, geo_toxo,
                       pc1_toxo, pc2_toxo, pc3_toxo,
                       mmrr_toxo,
                       effect = "geo")

p_pc1  <- plot_partial(gen_toxo, geo_toxo,
                       pc1_toxo, pc2_toxo, pc3_toxo,
                       mmrr_toxo,
                       effect = "PC1")

p_pc2  <- plot_partial(gen_toxo, geo_toxo,
                       pc1_toxo, pc2_toxo, pc3_toxo,
                       mmrr_toxo,
                       effect = "PC2")

p_pc3  <- plot_partial(gen_toxo, geo_toxo,
                       pc1_toxo, pc2_toxo, pc3_toxo,
                       mmrr_toxo,
                       effect = "PC3")
library(patchwork)

combined_poec <- (p_geo | p_pc1) /
  (p_pc2 | p_pc3)

combined_poec


p <- combined_poec +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(
      family = "Helvetica",
      size   = 18
    )
  )

p


getwd()

setwd('/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/out/mmrr_regression_plots')

ggsave("chamaea_neutral_snps_partial_MMRR.png",
       plot = p,
       width = 9.5,
       height = 10,
       units = "in",
       dpi = 600,
       bg = "white")





cor(unfold(pc1_pico), unfold(geo_pico))






ibr_cham_log <- log(ibr_cham + 1)
diag(ibr_cham_log) <- 0
ibr_cham_log <- (ibr_cham_log + t(ibr_cham_log)) / 2
round(pca_env_global$rotation[,1], 3)



cor(unfold(gen_cham), unfold(ibr_cham))


cor(unfold(gen_cham), unfold(ibr_cham_log))

