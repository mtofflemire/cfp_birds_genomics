library(terra)
library(sf)
library(ggplot2)
library(ggrepel)


coords <- read.table(
  '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/chamaea/coords.txt',
  header = FALSE
)

coords_sf <- st_as_sf(coords, coords = c(2,1), crs = 4326)
coords_vect <- vect(coords_sf)



climate_files <- list.files(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Shared-data/climate/wc2",
  pattern = "bio_.*\\.tif$",
  full.names = TRUE
)

climate <- rast(climate_files)









env_vals <- terra::extract(climate, coords_vect, ID = FALSE)
head(env_vals)



env_scaled <- scale(env_vals)

pca_env <- prcomp(env_scaled, center = TRUE, scale. = TRUE)

pc_scores <- as.data.frame(pca_env$x[, 1:3])
colnames(pc_scores) <- c("PC1", "PC2", "PC3")

loadings <- as.data.frame(pca_env$rotation[, 1:3])
loadings$BIO <- rownames(loadings)
var_exp <- summary(pca_env)$importance[2, 1:3] * 100


library(ggplot2)
library(ggrepel)

loadings$BIO_num <- sub(".*bio_", "BIO", loadings$BIO)
# % variance
var_exp <- 100 * pca_env$sdev^2 / sum(pca_env$sdev^2)

# Loadings for PC1/PC2 on correlation circle
loadings <- as.data.frame(pca_env$rotation[, 1:2])

loadings$PC1 <- loadings$PC1 * pca_env$sdev[1]
loadings$PC2 <- loadings$PC2 * pca_env$sdev[2]

loadings$BIO <- paste0("BIO", 1:19)

# magnitude of contribution to PC1+PC2
loadings$contrib <- sqrt(loadings$PC1^2 + loadings$PC2^2)
circle <- data.frame(
  x = cos(seq(0, 2*pi, length.out = 200)),
  y = sin(seq(0, 2*pi, length.out = 200))
)

ggplot() +
  geom_path(data = circle, aes(x, y), color = "grey70", linewidth = 1) +
  geom_hline(yintercept = 0, linewidth = 1.2) +
  geom_vline(xintercept = 0, linewidth = 1.2) +
  
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = PC1, yend = PC2, color = contrib),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 1
  ) +
  
  geom_text_repel(
    data = loadings,
    aes(PC1, PC2, label = BIO, color = contrib),
    size = 4,
    fontface = "bold"
  ) +
  
  scale_color_gradientn(
    colors = c("#2CB1A1", "#F1C40F", "red2"),
    name = "Contribution"
  ) +
  
  coord_fixed() +
  theme_minimal(base_size = 10) +
  labs(
    title = "PCA of WorldClim BIO variables",
    x = paste0("PC1 (", round(var_exp[1],1), "%)"),
    y = paste0("PC2 (", round(var_exp[2],1), "%)")
  )





# Loadings for PC1 & PC3
load_13 <- as.data.frame(pca_env$rotation[, c(1,3)])

load_13$PC1 <- load_13$PC1 * pca_env$sdev[1]
load_13$PC3 <- load_13$PC3 * pca_env$sdev[3]

load_13$BIO <- paste0("BIO", 1:19)

load_13$contrib <- sqrt(load_13$PC1^2 + load_13$PC3^2)

circle <- data.frame(
  x = cos(seq(0, 2*pi, length.out = 200)),
  y = sin(seq(0, 2*pi, length.out = 200))
)


ggplot() +
  geom_path(data = circle, aes(x, y), color = "grey70", linewidth = 1) +
  geom_hline(yintercept = 0, linewidth = 1.2) +
  geom_vline(xintercept = 0, linewidth = 1.2) +
  
  geom_segment(
    data = load_13,
    aes(x = 0, y = 0, xend = PC1, yend = PC3, color = contrib),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 1
  ) +
  
  geom_text_repel(
    data = load_13,
    aes(PC1, PC3, label = BIO, color = contrib),
    size = 4,
    fontface = "bold"
  ) +
  
  scale_color_gradientn(
    colors = c("#2CB1A1", "#F1C40F", "red2"),
    name = "Contribution"
  ) +
  
  coord_fixed() +
  theme_minimal(base_size = 10) +
  labs(
    title = "PCA of WorldClim BIO variables",
    x = paste0("PC1 (", round(var_exp[1],1), "%)"),
    y = paste0("PC3 (", round(var_exp[3],1), "%)")
  )












pc1_mat <- as.matrix(dist(pc_scores$PC1))
pc2_mat <- as.matrix(dist(pc_scores$PC2))
pc3_mat <- as.matrix(dist(pc_scores$PC3))



setdiff(rownames(genMat), rownames(geoMat))
all(rownames(genMat) == rownames(geoMat))


geoMat <- as.matrix(dist(st_coordinates(coords_sf)))
genMat <- as.matrix(
  read.table(
    '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/chamaea/mmrr/chamaea.gendist.txt'
  )
)

dim(genMat)
dim(geoMat)
dim(pc1_mat)
dim(pc2_mat)
dim(pc3_mat)






MMRR<-function(Y,X,nperm=999){
  #compute regression coefficients and test statistics
  nrowsY<-nrow(Y)
  y<-unfold(Y)
  if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
  Xmats<-sapply(X,unfold)
  fit<-lm(y~Xmats)
  coeffs<-fit$coefficients
  summ<-summary(fit)
  r.squared<-summ$r.squared
  tstat<-summ$coefficients[,"t value"]
  Fstat<-summ$fstatistic[1]
  tprob<-rep(1,length(tstat))
  Fprob<-1
  
  #perform permutations
  for(i in 1:nperm){
    rand<-sample(1:nrowsY)
    Yperm<-Y[rand,rand]
    yperm<-unfold(Yperm)
    fit<-lm(yperm~Xmats)
    summ<-summary(fit)
    Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
    tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
  }
  
  #return values
  tp<-tprob/(nperm+1)
  Fp<-Fprob/(nperm+1)
  names(r.squared)<-"r.squared"
  names(coeffs)<-c("Intercept",names(X))
  names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
  names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
  names(Fstat)<-"F-statistic"
  names(Fp)<-"F p-value"
  return(list(r.squared=r.squared,
              coefficients=coeffs,
              tstatistic=tstat,
              tpvalue=tp,
              Fstatistic=Fstat,
              Fpvalue=Fp))
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR

unfold<-function(X){
  x<-vector()
  for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
  x<-scale(x, center=TRUE, scale=TRUE)  # Comment this line out if you wish to perform the analysis without standardizing the distance matrices! 
  return(x)
}



X_list <- list(
  geodist = geoMat,
  PC1 = pc1_mat,
  PC2 = pc2_mat,
  PC3 = pc3_mat
)




set.seed(12)
mmrr_poecile <- MMRR(genMat, X_list, nperm = 999)

mmrr_poecile




print_mmrr_row_for_excel <- function(res, species_name){
  
  row <- c(
    Species = species_name,
    R2      = round(res$r.squared, 3),
    
    b_Geo = round(res$coefficients["geodist"], 3),
    p_Geo = res$tpvalue["geodist(p)"],
    
    b_PC1 = round(res$coefficients["PC1"], 3),
    p_PC1 = res$tpvalue["PC1(p)"],
    
    b_PC2 = round(res$coefficients["PC2"], 3),
    p_PC2 = res$tpvalue["PC2(p)"],
    
    b_PC3 = round(res$coefficients["PC3"], 3),
    p_PC3 = res$tpvalue["PC3(p)"]
  )
  
  cat(paste(row, collapse = "\t"), "\n")
}


cat("Species\tR2\tb_Geo\tp_Geo\tb_PC1\tp_PC1\tb_PC2\tp_PC2\tb_PC3\tp_PC3\n")
print_mmrr_row_for_excel(mmrr_poecile,  "Wrentit")











gen_s <- unfold(genMat)
geo_s <- unfold(geoMat)
pc1_s <- unfold(pc1_mat)
pc2_s <- unfold(pc2_mat)
pc3_s <- unfold(pc3_mat)   # add this if you haven’t yet











df_ibd <- data.frame(Genetic = gen_s, Geo = geo_s)

fit <- lm(Genetic ~ Geo, data = df_ibd)

beta <- round(coef(fit)[2], 3)
r2   <- round(summary(fit)$r.squared, 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)

label_txt <- paste0("β = ", beta,
                    "\nR² = ", r2,
                    "\np = ", pval)

ggplot(df_ibd, aes(Geo, Genetic)) +
  geom_point(shape=21, fill="gray65", color="gray65",
             alpha=1, size=1, stroke=1) +
  geom_smooth(method="lm", color="black") +
  
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1.5,
           size = 5,
           label = paste0(
             "beta == ", beta,
             "*','~~R^2 == ", r2,
             "*','~~p == ", pval
           ),
           parse = TRUE) +
  
  theme_bw() +
  labs(x="IBD", y="Genetic Distance") +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 12)
  )









df_pc1 <- data.frame(Genetic = gen_s, PC1 = pc1_s)

fit <- lm(Genetic ~ PC1, data = df_pc1)

beta <- round(coef(fit)[2], 3)
r2   <- round(summary(fit)$r.squared, 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)

ggplot(df_pc1, aes(PC1, Genetic)) +
  geom_point(
    shape = 21, fill = "gray65", color = "gray65",
    alpha = 1, size = 1, stroke = 1
  ) +
  geom_smooth(method = "lm", color = "black") +
  
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1,
           size = 4.25,
           label = paste0(
             "beta == ", beta,
             "*','~~R^2 == ", r2,
             "*','~~p == ", pval
           ),
           parse = TRUE) +
  
  theme_bw() +
  labs(
    x = "IBE (env_PC1)",
    y = "Genetic Distance"
  ) +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 12)
  )





df_pc2 <- data.frame(Genetic = gen_s, PC2 = pc2_s)

fit <- lm(Genetic ~ PC2, data = df_pc2)

beta <- round(coef(fit)[2], 3)
r2   <- round(summary(fit)$r.squared, 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)

ggplot(df_pc2, aes(PC2, Genetic)) +
  geom_point(
    shape = 21, fill = "gray65", color = "gray65",
    alpha = 1, size = 1, stroke = 1
  ) +
  geom_smooth(method = "lm", color = "black") +
  
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1,
           size = 4.25,
           label = paste0(
             "beta == ", beta,
             "*','~~R^2 == ", r2,
             "*','~~p == ", pval
           ),
           parse = TRUE) +
  
  theme_bw() +
  labs(
    x = "IBE (env_PC2)",
    y = "Genetic Distance"
  ) +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 12)
  )







df_pc3 <- data.frame(Genetic = gen_s, PC3 = pc3_s)

fit <- lm(Genetic ~ PC3, data = df_pc3)

beta <- round(coef(fit)[2], 3)
r2   <- round(summary(fit)$r.squared, 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)

ggplot(df_pc3, aes(PC3, Genetic)) +
  geom_point(
    shape = 21, fill = "gray65", color = "gray65",
    alpha = 1, size = 1, stroke = 1
  ) +
  geom_smooth(method = "lm", color = "black") +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1,
           size = 4.25,
           label = paste0(
             "beta == ", beta,
             "*','~~R^2 == ", r2,
             "*','~~p == ", pval
           ),
           parse = TRUE) +
  theme_bw() +
  labs(
    x = "IBE (env_PC3)",
    y = "Genetic Distance"
  ) +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 12)
  )









plot_geo_vs_env(pc1_s, "PC1")
plot_geo_vs_env(pc2_s, "PC2")
plot_geo_vs_env(pc3_s, "PC3")












pc1_load <- pca_env$rotation[, "PC1"]
pc2_load <- pca_env$rotation[, "PC2"]
pc3_load <- pca_env$rotation[, "PC3"]




pc1_raster <- sum(climate * pc1_load)
pc2_raster <- sum(climate * pc2_load)
pc3_raster <- sum(climate * pc3_load)



library(sf)
library(terra)
library(tigris)


options(tigris_use_cache = TRUE)

states <- states(cb = TRUE, year = 2023)

ca_or <- states[states$STUSPS %in% c("CA", "OR", "WA"), ]


crs(climate)   # check raster CRS

ca_or <- st_transform(ca_or, crs(climate))



ca_or <- st_union(ca_or)
ca_or <- st_as_sf(ca_or)


plot(vect(ca_or))
plot(climate[[1]], add = TRUE)


pc1_mask <- mask(crop(pc1_raster, vect(ca_or)), vect(ca_or))
pc2_mask <- mask(crop(pc2_raster, vect(ca_or)), vect(ca_or))
pc3_mask <- mask(crop(pc3_raster, vect(ca_or)), vect(ca_or))



pc1_df <- as.data.frame(pc1_mask, xy = TRUE)
names(pc1_df)[3] <- "PC1"

pc2_df <- as.data.frame(pc2_mask, xy = TRUE)
names(pc2_df)[3] <- "PC2"

pc3_df <- as.data.frame(pc3_mask, xy = TRUE)
names(pc3_df)[3] <- "PC3"


ggplot() +
  geom_raster(data = pc1_df, aes(x, y, fill = PC1)) +
  geom_sf(data = ca_or, fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_gradientn(
    colours = terrain.colors(100),
    name = "PC1"
  ) +
  theme_void() +
  labs(title = "PC1 Climate Gradient")


ggplot() +
  geom_raster(data = pc2_df, aes(x, y, fill = PC2)) +
  geom_sf(data = ca_or, fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_gradientn(
    colours = terrain.colors(100),
    name = "PC2"
  ) +
  theme_void() +
  labs(title = "PC2 Climate Gradient")


ggplot() +
  geom_raster(data = pc3_df, aes(x, y, fill = PC3)) +
  geom_sf(data = ca_or, fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_gradientn(
    colours = terrain.colors(100),
    name = "PC3"
  ) +
  theme_void() +
  labs(title = "PC3 Climate Gradient")





scale01 <- function(r) {
  mm <- minmax(r)
  (r - mm[1]) / (mm[2] - mm[1])
}


pc1_s <- scale01(pc1_mask)
pc2_s <- scale01(pc2_mask)
pc3_s <- scale01(pc3_mask)



rgb_stack <- c(pc1_s, pc2_s, pc3_s)
names(rgb_stack) <- c("R", "G", "B")



plotRGB(rgb_stack, r=1, g=2, b=3, stretch="lin")
plot(vect(ca_or), add=TRUE)



# Make sure points are in same CRS as raster (WGS84)
coords_sf <- st_transform(coords_sf, crs(rgb_stack))

# Extract lon/lat
pts <- st_coordinates(coords_sf)

plotRGB(rgb_stack, r=1, g=2, b=3, stretch="lin")
plot(vect(ca_or), add=TRUE, lwd=1.5)

points(
  pts[,1], pts[,2],
  pch = 21,          # circle with fill
  bg  = "white",     # fill color
  col = "black",     # outline
  cex = 1.4,
  lwd = 1.2
)











library(sf)

# Helper function
read_species_coords <- function(path, species_name) {
  df <- read.table(path, header = FALSE)
  sf_obj <- st_as_sf(df, coords = c(2,1), crs = 4326)
  sf_obj$species <- species_name
  return(sf_obj)
}

# Read each species
chamaea_sf <- read_species_coords(
  '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/chamaea/coords.txt',
  "Wrentit"
)

picoides_sf <- read_species_coords(
  '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/picoides/coords.txt',
  "Picoides"
)

poecile_sf <- read_species_coords(
  '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/poecile/coords.txt',
  "Poecile"
)

toxostoma_sf <- read_species_coords(
  '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/California_birds_wgs/data/toxostoma/coords.txt',
  "Toxostoma"
)

# Combine
all_species_sf <- rbind(
  chamaea_sf,
  picoides_sf,
  poecile_sf,
  toxostoma_sf
)



all_species_sf <- st_transform(all_species_sf, crs(rgb_stack))

pts_all <- st_coordinates(all_species_sf)
species_factor <- factor(all_species_sf$species)

default_cols <- scales::hue_pal()(length(levels(species_factor)))
names(default_cols) <- levels(species_factor)



plotRGB(rgb_stack, r=1, g=2, b=3, stretch="lin")
plot(vect(ca_or), add=TRUE, lwd=1.5)

# Shapes per species
shapes <- c(21, 22, 24, 23)
names(shapes) <- levels(species_factor)

# Plot each species
for(sp in levels(species_factor)) {
  
  idx <- which(all_species_sf$species == sp)
  
  points(
    pts_all[idx,1],
    pts_all[idx,2],
    pch = shapes[sp],
    bg  = default_cols[sp],
    col = "black",
    cex = 1.8,
    lwd = 1
  )
}





# Open clean plotting device


par(mar = c(2,2,2,2))

plot.new()

legend(
  "center",
  legend = levels(species_factor),
  pch = shapes,
  pt.bg = default_cols,
  col = "black",
  pt.cex = 2,
  cex = 1.3,
  bty = "n",
  title = "Species"
)




















library(rnaturalearth)
library(sf)
library(terra)

na <- ne_countries(continent = "North America", returnclass = "sf")
na <- st_union(na)
na <- st_transform(na, crs(climate))
na_vect <- vect(na)





scale01 <- function(r){
  mm <- minmax(r)
  (r - mm[1]) / (mm[2] - mm[1])
}

pc1_s <- scale01(pc1_raster)
pc2_s <- scale01(pc2_raster)
pc3_s <- scale01(pc3_raster)

rgb_stack <- c(pc1_s, pc2_s, pc3_s)
names(rgb_stack) <- c("R","G","B")




rgb_crop <- crop(rgb_stack, na_vect)
rgb_na   <- mask(rgb_crop, na_vect)




aea <- "EPSG:5070"  # North America Albers Equal Area

rgb_aea <- project(rgb_na, aea)
na_aea  <- project(na_vect, aea)


plotRGB(rgb_na, r=1, g=2, b=3, stretch="lin")
plot(na_vect, add=TRUE, lwd=1.5)























#testing




# --- PC3 from PCA ---
pc3_load   <- pca_env$rotation[, "PC3"]
pc3_raster <- sum(climate * pc3_load)


library(tigris)
library(sf)
library(terra)

states <- states(cb = TRUE, year = 2023)
ca_or  <- states[states$STUSPS %in% c("CA","OR"), ]

ca_or  <- st_transform(ca_or, crs(climate))
ca_or_vect <- vect(st_union(ca_or))

pc3_mask <- mask(crop(pc3_raster, ca_or_vect), ca_or_vect)



pc3_df <- as.data.frame(pc3_mask, xy = TRUE)
colnames(pc3_df)[3] <- "PC3"



colnames(coords) <- c("lon","lat")

coords_vect <- vect(coords,
                    geom = c("lon","lat"),
                    crs  = "EPSG:4326")

coords_vect <- project(coords_vect, crs(pc3_mask))

pts <- as.data.frame(crds(coords_vect))
colnames(pts) <- c("x","y")


plot(pc3_mask)
points(pts$x, pts$y, col="red", pch=20)



library(ggplot2)

ggplot() +
  geom_raster(data = pc3_df, aes(x, y, fill = PC3)) +
  geom_point(data = pts, aes(x, y),
             shape = 21, fill = "white",
             color = "black", size = 2.5, stroke = 1) +
  geom_sf(data = ca_or, fill = NA, color = "black", linewidth = 0.6) +
  scale_fill_gradientn(
    colours = terrain.colors(100),
    name = "PC3\nThermal amplitude"
  ) +
  theme_void() +
  labs(title = "PC3 Climate Gradient (Thermal Amplitude / Continentality)")





# extract PC3 climate score at each bird location
pc3_vals <- terra::extract(pc3_mask, coords_vect, ID = FALSE)

coords$PC3 <- pc3_vals[,1]




library(ggplot2)

ggplot(coords, aes(x = reorder(row.names(coords), PC3), y = PC3)) +
  geom_point(size = 3) +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Individuals (ordered by PC3 climate)",
    y = "PC3 climate score",
    title = "Birds ordered along PC3 thermal/continental gradient"
  )





ggplot() +
  geom_raster(data = pc3_df, aes(x, y, fill = PC3)) +
  geom_point(
    data = coords,
    aes(x = lon, y = lat, color = PC3),
    size = 3
  ) +
  scale_color_viridis_c(option = "plasma") +
  scale_fill_gradientn(colours = terrain.colors(100)) +
  theme_void()


pca_env$rotation[,2]



library(sf)
library(ggplot2)

# coords file you just made (in VCF order)
coords <- read.table(
  "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile_coords_ordered.txt",
  col.names = c("lat","lon")
)

# PC scores from your PCA
pc_df <- as.data.frame(pca_env$x[,1:3])
pc_df$PC2 <- pc_df[,2]

# combine
map_df <- cbind(coords, pc_df)


pts_sf <- st_as_sf(map_df,
                   coords = c("lon","lat"),
                   crs = 4326)


library(rnaturalearth)
states <- ne_states(country = "United States of America", returnclass = "sf")


ggplot() +
  geom_sf(data = states, fill = "gray95", color = "gray70") +
  geom_sf(data = pts_sf,
          aes(color = PC2),
          size = 4) +
  scale_color_viridis_c(option = "plasma") +
  theme_bw() +
  labs(color = "PC2 (Seasonality axis)")


ggplot() +
  geom_sf(data = west_states,
          fill = "gray95",
          color = "gray60",
          linewidth = 0.4) +
  
  geom_sf(data = pts_sf,
          aes(fill = PC2),
          shape = 21,
          color = "black",
          size = 4,
          stroke = 0.6) +
  
  scale_fill_viridis_c(option = "plasma") +
  
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"])
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  
  labs(
    fill = "PC2\n(Seasonality axis)"
  )



library(rnaturalearth)
library(dplyr)

states <- ne_states(
  country = "United States of America",
  returnclass = "sf"
)

west_states <- states %>%
  filter(name %in% c("California", "Oregon", "Washington"))


# Tight crop around your samples
bbox <- st_bbox(pts_sf)

ggplot() +
  geom_sf(data = west_states,
          fill = "gray95",
          color = "gray60",
          linewidth = 0.4) +
  
  geom_sf(data = pts_sf,
          aes(fill = PC2),
          shape = 21,
          color = "black",
          size = 4,
          stroke = 0.6) +
  
  scale_fill_viridis_c(option = "plasma") +
  
  coord_sf(
    xlim = c(bbox["xmin"] - 1, bbox["xmax"] + 1),
    ylim = c(bbox["ymin"] - 1, bbox["ymax"] + 1)
  ) +
  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  
  labs(fill = "PC2\n(Seasonality axis)")







df_pc1_geo <- data.frame(Geo = geo_s, Env = pc1_s)

fit <- lm(Env ~ Geo, data = df_pc1_geo)

beta <- round(coef(fit)[2], 3)
r2   <- round(summary(fit)$r.squared, 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)

ggplot(df_pc1_geo, aes(Geo, Env)) +
  geom_point(shape=21, fill="gray65", color="gray65",
             alpha=1, size=1, stroke=1) +
  geom_smooth(method="lm", color="black") +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1,
           size = 4.25,
           label = paste0(
             "beta == ", beta,
             "*','~~R^2 == ", r2,
             "*','~~p == ", pval
           ),
           parse = TRUE) +
  theme_bw() +
  labs(
    x = "Geographic Distance",
    y = "Environmental Distance (PC1)"
  ) +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 12)
  )






df_pc2_geo <- data.frame(Geo = geo_s, Env = pc2_s)

fit <- lm(Env ~ Geo, data = df_pc2_geo)

beta <- round(coef(fit)[2], 3)
r2   <- round(summary(fit)$r.squared, 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)

ggplot(df_pc2_geo, aes(Geo, Env)) +
  geom_point(shape=21, fill="gray65", color="gray65",
             alpha=1, size=1, stroke=1) +
  geom_smooth(method="lm", color="black") +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1,
           size = 4.25,
           label = paste0(
             "beta == ", beta,
             "*','~~R^2 == ", r2,
             "*','~~p == ", pval
           ),
           parse = TRUE) +
  theme_bw() +
  labs(
    x = "Geographic Distance",
    y = "Environmental Distance (PC2)"
  ) +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 12)
  )




df_pc3_geo <- data.frame(Geo = geo_s, Env = pc3_s)

fit <- lm(Env ~ Geo, data = df_pc3_geo)

beta <- round(coef(fit)[2], 3)
r2   <- round(summary(fit)$r.squared, 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)

ggplot(df_pc3_geo, aes(Geo, Env)) +
  geom_point(shape=21, fill="gray65", color="gray65",
             alpha=1, size=1, stroke=1) +
  geom_smooth(method="lm", color="black") +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1,
           size = 4.25,
           label = paste0(
             "beta == ", beta,
             "*','~~R^2 == ", r2,
             "*','~~p == ", pval
           ),
           parse = TRUE) +
  theme_bw() +
  labs(
    x = "Geographic Distance",
    y = "Environmental Distance (PC3)"
  ) +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 12)
  )













gen_s <- unfold(genMat)
geo_s <- unfold(geoMat)
pc1_s <- unfold(pc1_mat)
pc2_s <- unfold(pc2_mat)
pc3_s <- unfold(pc3_mat)









gen_resid <- resid(lm(gen_s ~ pc1_s + pc2_s + pc3_s))
geo_resid <- resid(lm(geo_s ~ pc1_s + pc2_s + pc3_s))

df_ibd <- data.frame(Genetic = gen_resid, Geo = geo_resid)

fit <- lm(Genetic ~ Geo, data = df_ibd)

beta <- round(coef(fit)[2], 3)
r2   <- round(summary(fit)$r.squared, 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)

ggplot(df_ibd, aes(Geo, Genetic)) +
  geom_point(shape=21, fill="gray75", color="gray75",
             alpha=1, size=1, stroke=1) +
  geom_smooth(method="lm", color="black") +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1.5,
           size = 5,
           label = paste0(
             "beta == ", beta,
             "*','~~R^2 == ", r2,
             "*','~~p == ", pval
           ),
           parse = TRUE) +
  theme_bw() +
  labs(x="IBD (Geo | PCs removed)", y="Genetic Distance (PCs removed)", title="Poecile")









gen_resid <- resid(lm(gen_s ~ geo_s + pc2_s + pc3_s))
pc1_resid <- resid(lm(pc1_s ~ geo_s + pc2_s + pc3_s))

df_pc1 <- data.frame(Genetic = gen_resid, PC1 = pc1_resid)

fit <- lm(Genetic ~ PC1, data = df_pc1)

beta <- round(coef(fit)[2], 3)
r2   <- round(summary(fit)$r.squared, 3)
pval <- signif(mmrr_poecile$tpvalue["PC1(p)"], 3)


ggplot(df_pc1, aes(PC1, Genetic)) +
  geom_point(shape=21, fill="gray75", color="gray75",
             alpha=1, size=1, stroke=1) +
  geom_smooth(method="lm", color="black") +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1.5,
           size = 5,
           label = paste0(
             "beta == ", beta,
             "*','~~R^2 == ", r2,
             "*','~~p == ", pval
           ),
           parse = TRUE) +
  theme_bw() +
  labs(
    x="IBE PC1 (Geo, PC2, PC3 removed)",
    y="Genetic Distance (Geo, PC2, PC3 removed)",
    title="Poecile"
  )



gen_s <- unfold(genMat)
geo_s <- unfold(geoMat)
pc1_s <- unfold(pc1_mat)
pc2_s <- unfold(pc2_mat)
pc3_s <- unfold(pc3_mat)



gen_resid <- resid(lm(gen_s ~ geo_s + pc1_s + pc3_s))
pc2_resid <- resid(lm(pc2_s ~ geo_s + pc1_s + pc3_s))

df_pc2 <- data.frame(Genetic = gen_resid, PC2 = pc2_resid)

fit <- lm(Genetic ~ PC2, data = df_pc2)

beta <- round(coef(fit)[2], 3)
r2   <- round(summary(fit)$r.squared, 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)

ggplot(df_pc2, aes(PC2, Genetic)) +
  geom_point(shape=21, fill="gray75", color="gray75",
             alpha=1, size=1, stroke=1) +
  geom_smooth(method="lm", color="black") +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1.5,
           size = 5,
           label = paste0(
             "beta == ", beta,
             "*','~~R^2 == ", r2,
             "*','~~p == ", pval
           ),
           parse = TRUE) +
  theme_bw() +
  labs(
    x="IBE PC2 (Geo, PC1, PC3 removed)",
    y="Genetic Distance (Geo, PC1, PC3 removed)",
    title="Poecile"
  )






gen_resid <- resid(lm(gen_s ~ geo_s + pc1_s + pc2_s))
pc3_resid <- resid(lm(pc3_s ~ geo_s + pc1_s + pc2_s))

df_pc3 <- data.frame(Genetic = gen_resid, PC3 = pc3_resid)

fit <- lm(Genetic ~ PC3, data = df_pc3)

beta <- round(coef(fit)[2], 3)
r2   <- round(summary(fit)$r.squared, 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)

ggplot(df_pc3, aes(PC3, Genetic)) +
  geom_point(shape=21, fill="gray75", color="gray75",
             alpha=1, size=1, stroke=1) +
  geom_smooth(method="lm", color="black") +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1.5,
           size = 5,
           label = paste0(
             "beta == ", beta,
             "*','~~R^2 == ", r2,
             "*','~~p == ", pval
           ),
           parse = TRUE) +
  theme_bw() +
  labs(
    x="IBE PC3 (Geo, PC1, PC2 removed)",
    y="Genetic Distance (Geo, PC1, PC2 removed)",
    title="Poecile"
  )

