# paths
eigenval_path <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/1_pca/toxostoma_pca.eigenval'
eigenvec_path <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/1_pca/toxostoma_pca.eigenvec'

# eigenvalues → % variance
eigvals <- scan(eigenval_path, quiet = TRUE)
pc_percent <- 100 * eigvals / sum(eigvals)

# eigenvectors → coordinates
pca <- read.table(eigenvec_path, header = TRUE, stringsAsFactors = FALSE)



library(ggplot2)

ggplot(pca, aes(PC1, PC2)) +
  geom_point(color = "gray40", size = 3) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", round(pc_percent[1], 2), "%)"),
    y = paste0("PC2 (", round(pc_percent[2], 2), "%)")
  )




library(ggplot2)
library(ggrepel)

ggplot(pca, aes(PC1, PC2)) +
  geom_point(color = "gray40", size = 3) +
  geom_text_repel(
    aes(label = IID),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", round(pc_percent[1], 2), "%)"),
    y = paste0("PC2 (", round(pc_percent[2], 2), "%)")
  )



geom_text_repel(
  aes(label = IID),
  size = 2.5,
  color = "gray20",
  max.overlaps = Inf
)
