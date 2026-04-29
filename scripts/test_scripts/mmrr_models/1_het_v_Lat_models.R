library(tidyverse)

het_file  <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/41-Chamaea.diversity.het'
meta_file <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/chamaea/Cfa-metadata-nsamp159copy.csv'

# Read PLINK .het
het <- read.table(het_file, header = TRUE, check.names = FALSE)

het <- het %>%
  mutate(
    H_obs = (`N(NM)` - `O(HOM)`) / `N(NM)`
  ) %>%
  rename(sampleID = IID)

# Read metadata
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Join
dat <- left_join(het, meta, by = "sampleID")

# Keep rows with latitude
dat <- dat %>%
  filter(!is.na(lat))

# Regression
model <- lm(H_obs ~ lat, data = dat)
summary(model)

# Plot
ggplot(dat, aes(x = lat, y = H_obs)) +
  geom_point(size = 3, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Latitude",
    y = "Observed heterozygosity",
    title = "Heterozygosity vs Latitude — Chamaea"
  )


ggplot(dat, aes(x = lat, y = H_obs)) +
  geom_point(size = 3, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey70") +
  theme_bw() +
  labs(
    x = "Latitude",
    y = "Observed heterozygosity",
    title = "Heterozygosity vs Latitude — Chamaea"
  ) +
  annotate("text",
           x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.2,
           size = 4,
           parse = TRUE,
           label = "beta == -0.00424*','~~R^2 == 0.47*','~~p < 2 %*% 10^-16")



















library(tidyverse)
library(stringr)

het_file  <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/Toxostoma_filteredQC_2.het'
meta_file <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/toxostoma/toxostoma-meta.csv'

# Read PLINK .het
het <- read.table(het_file, header = TRUE, check.names = FALSE)

het <- het %>%
  mutate(
    H_obs   = (`N(NM)` - `O(HOM)`) / `N(NM)`,
    sampleID = str_remove(IID, "_S\\d+_L\\d+$")
  )

# Read metadata
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Join
dat <- left_join(het, meta, by = "sampleID") %>%
  filter(!is.na(decimallatitude))

# Regression
model <- lm(H_obs ~ decimallatitude, data = dat)
summary(model)

# Plot
ggplot(dat, aes(x = decimallatitude, y = H_obs)) +
  geom_point(size = 3, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Latitude",
    y = "Observed heterozygosity",
    title = "Heterozygosity vs Latitude — Toxostoma"
  )


ggplot(dat, aes(x = decimallatitude, y = H_obs)) +
  geom_point(size = 3, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey70") +
  theme_bw() +
  labs(
    x = "Latitude",
    y = "Observed heterozygosity",
    title = "Heterozygosity vs Latitude — Toxostoma"
  ) +
  annotate("text",
           x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.2,
           size = 4,
           parse = TRUE,
           label = "beta == -0.00320*','~~R^2 == 0.37*','~~p == 2.36 %*% 10^-4")

















library(tidyverse)

het_file  <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/picoides.diversity.het'
meta_file <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/picoides/picoides_Meta.csv'

# Read PLINK .het
het <- read.table(het_file, header = TRUE, check.names = FALSE)

het <- het %>%
  mutate(
    H_obs = (`N(NM)` - `O(HOM)`) / `N(NM)`
  ) %>%
  rename(sampleID2 = IID)

# Read metadata
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Join
dat <- left_join(het, meta, by = "sampleID2") %>%
  filter(!is.na(decimallatitude))

# Regression
model <- lm(H_obs ~ decimallatitude, data = dat)
summary(model)

# Plot
ggplot(dat, aes(x = decimallatitude, y = H_obs)) +
  geom_point(size = 3, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Latitude",
    y = "Observed heterozygosity",
    title = "Heterozygosity vs Latitude — Picoides"
  )



ggplot(dat, aes(x = decimallatitude, y = H_obs)) +
  geom_point(size = 3, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey70") +
  theme_bw() +
  labs(
    x = "Latitude",
    y = "Observed heterozygosity",
    title = "Heterozygosity vs Latitude — Picoides"
  ) +
  annotate("text",
           x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.2,
           size = 4,
           parse = TRUE,
           label = "beta == 0.00146*','~~R^2 == 0.11*','~~p == 0.060")







library(tidyverse)
library(stringr)

het_file  <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/poecile.diversity.autosomes.het'
meta_file <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile_metadata.csv'

# Read PLINK .het
het <- read.table(het_file, header = TRUE, check.names = FALSE)

het <- het %>%
  mutate(
    H_obs = (`N(NM)` - `O(HOM)`) / `N(NM)`,
    SampleID = str_remove(IID, "_S\\d+_L\\d+$")
  )

# Read metadata
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Join
dat <- left_join(het, meta, by = "SampleID") %>%
  filter(!is.na(decimallatitude))

# Regression
model <- lm(H_obs ~ decimallatitude, data = dat)
summary(model)

# Plot
ggplot(dat, aes(x = decimallatitude, y = H_obs)) +
  geom_point(size = 3, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Latitude",
    y = "Observed heterozygosity",
    title = "Heterozygosity vs Latitude — Poecile"
  )











library(tidyverse)
library(stringr)

het_file  <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/out/het/poecile.diversity.autosomes.het'
meta_file <- '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Cfp_birds_wgs/data/poecile/poecile_metadata.csv'

# Read PLINK .het
het <- read.table(het_file, header = TRUE, check.names = FALSE)

het <- het %>%
  mutate(
    H_obs   = (`N(NM)` - `O(HOM)`) / `N(NM)`,
    sampleID = str_remove(IID, "_S\\d+_L\\d+$")
  )

# Read metadata
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Join
dat <- left_join(het, meta, by = "sampleID") %>%
  filter(!is.na(decimallatitude))

# Regression
model <- lm(H_obs ~ decimallatitude, data = dat)
summary(model)







ggplot(dat, aes(x = decimallatitude, y = H_obs)) +
  geom_point(size = 3, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey70") +
  theme_bw() +
  labs(
    x = "Latitude",
    y = "Observed heterozygosity",
    title = "Heterozygosity vs Latitude — Poecile"
  ) +
  annotate("text",
           x = -Inf, y = Inf,
           hjust = -0.1, vjust = 1.2,
           size = 4,
           parse = TRUE,
           label = "beta == 0.00380*','~~R^2 == 0.48*','~~p == 5.94 %*% 10^-5")
