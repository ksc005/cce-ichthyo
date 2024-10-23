# Script Header -----------------------------------------------------------


# Set Up ------------------------------------------------------------------
setwd("/Volumes/petrik-lab/kchen/cce-ichthyo")

# load libraries
library(plyr)
library(tidyverse)
library(reshape2)
library(readxl)
library(fastDummies)


# Read in eco traits table -------------------------------------------------
eco_traits <- read_excel("nhl-eco-traits.xlsx", sheet = 1)

eco_traits <- eco_traits %>%
  dplyr::rename(distribution = cross_shore_distribution) %>%
  mutate(season_4_NHL = as.factor(season_4_NHL),
         calculated_yearstart_NHL = as.factor(calculated_yearstart_NHL),
         group.3333_NHL = as.factor(group.3333_NHL),
         adult_habitat = as.factor(adult_habitat),
         distribution = as.factor(distribution),
         biogeographic_affinity = as.factor(biogeographic_affinity),
         fishing_status = as.factor(fishing_status),
         month_maximum = match(eco_traits$month_maximum_NHL, month.name)) %>%
  na.omit()

order_count <- eco_traits %>%
  dplyr::count(order) %>%
  dplyr::filter(n >= 3) #if there's less than three spp in the order, we'll have the order as 'Misc'

eco_traits$order <- as.factor(ifelse(eco_traits$order %in% order_count$order,
                                     as.character(eco_traits$order),
                                     "Misc"))


# PCA ---------------------------------------------------------------------
# create the dummy columns, and remove columns not needed for PCA
pca_data <- eco_traits %>%
  select(scientific_name, month_maximum, order, adult_habitat, distribution, biogeographic_affinity) %>%
  arrange(scientific_name) %>%
  dummy_cols(select_columns = c("adult_habitat", "distribution",
                                "biogeographic_affinity", "order"),
             remove_selected_columns =  T) %>%
  mutate(order_Misc = NULL,
         `biogeographic_affinity_Warm-water` = NULL) %>%
  dplyr::rename(Demersal = adult_habitat_Demersal,
                Epipelagic = adult_habitat_Epipelagic,
                Mesopelagic = adult_habitat_Mesopelagic,
                Coastal = distribution_Coastal,
                `Coastal-Oceanic` = `distribution_Coastal-Oceanic`,
                Oceanic = distribution_Oceanic,
                `Cool-water` = `biogeographic_affinity_Cool-water`,
                `Wide Distribution` = `biogeographic_affinity_Wide distribution`,
                Myctophiformes = order_Myctophiformes,
                `Perciformes/Cottoidei` = `order_Perciformes/Cottoidei`,
                Pleuronectiformes = order_Pleuronectiformes)

# perform PCA
pca.df <- pca_data[,-1]

scaled_pca.df <- as.data.frame(scale(pca.df)) #standardize variables

fit <- princomp(scaled_pca.df, cor=TRUE)

get_eig(fit)
summary(fit)

scores = as.data.frame(fit$scores)
scores$scientific_name = pca_data$scientific_name
names(scores) <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9",
                   "PC10", "PC11","PC12","SciName")

set.seed(70524)
fit_kmeans <- kmeans(scores, 3, nstart = 25) # the cluster number will be contained in the vector drop.fit_kmeans$cluster

# save files needed for plotting
eco_traits <- eco_traits %>%
  arrange(scientific_name)

markers <- cbind(scores, eco_traits$trend_NHL_ddec) %>%
  dplyr::rename(trend = `eco_traits$trend_NHL_ddec`) %>%
  mutate(direction = as.factor(ifelse(trend > 0, "positive", "negative")))

write.csv(markers, "pca-markers.nhl.csv", row.names = F)

correlations = as.data.frame(cor(scaled_pca.df, fit$scores))
write.csv(correlations, "pca-plots.nhl.csv", row.names = T)


# Species clustering closely together -------------------------------------
# attach the clusters to the species names; i.e., which species are in which clusters?
spp_cta <- read.csv("nhl-ct-summary.csv")
spp_cta <- spp_cta %>%
  group_by(Taxa) %>%
  mutate(ID = cur_group_id())

clusters <- data.frame(index = 1:length(unique(spp_cta$Taxa)),
                       cluster = fit_kmeans$cluster)

spp_cta$cluster <- clusters$cluster[match(spp_cta$ID, clusters$index)]

# cluster 1: positive on PC1, negative on PC2. epipelagic CO
cluster1 <- unique(spp_cta$Taxa[spp_cta$cluster == 1])

cluster1.table <- spp_cta %>%
  filter(cluster == 1)
cluster1.group <- cluster1.table %>%
  group_by(Taxa) %>%
  dplyr::count(group)
table(cluster1.group$group) # 0 spp shifting earlier, 0 later, 1 no change

cluster1_scores <- scores %>%
  select(PC1, PC2, SciName) %>%
  filter(SciName %in% cluster1)

cluster1_ct <- spp_cta %>%
  filter(Taxa %in% cluster1)

cluster1_regression <- lm(ct_anomaly_days ~ as.numeric(as.character(Period_3y)),
                               data = cluster1_ct)
print(summary(cluster1_regression))

# cluster 2: primarily positive on PC1,PC2. mesopelagic oceanic, myctophiformes
cluster2 <- unique(spp_cta$Taxa[spp_cta$cluster == 2])

cluster2.table <- spp_cta %>%
  filter(cluster == 2)
cluster2.group <- cluster2.table %>%
  group_by(Taxa) %>%
  dplyr::count(group)
table(cluster2.group$group) # 1 spp shifting earlier, 0 shifting later, 3 no change

cluster2_scores <- scores %>%
  select(PC1, PC2, SciName) %>%
  filter(SciName %in% cluster2)

cluster2_ct <- spp_cta %>%
  filter(Taxa %in% cluster2)

cluster2_regression <- lm(ct_anomaly_days ~ as.numeric(as.character(Period_3y)),
                               data = cluster2_ct)
print(summary(cluster2_regression))

# looks like the first 3y period has abnormally low CTa, let's repeat the regression and plot without it
cluster2_ct.excluded <- cluster2_ct %>%
  filter(Period_3y != 1)

cluster2_regression.excluded <- lm(ct_anomaly_days ~ as.numeric(as.factor(Period_3y)),
                                        data = cluster2_ct.excluded)
print(summary(cluster2_regression.excluded))

# third cluster: primarily negative on PC1. pleuronectiformes, p.cottoidei, coastal, demersal, coolwater
cluster3 <- unique(spp_cta$Taxa[spp_cta$cluster == 3])

cluster3.table <- spp_cta %>%
  filter(cluster == 3)
cluster3.group <- cluster3.table %>%
  group_by(Taxa) %>%
  dplyr::count(group)
table(cluster3.group$group) # 6 spp shifting earlier, 2 later, 12 no change

cluster3_scores <- scores %>%
  select(PC1, PC2, SciName) %>%
  filter(SciName %in% cluster3)

cluster3_ct <- spp_cta %>%
  filter(Taxa %in% cluster3)

cluster3_regression <- lm(ct_anomaly_days ~ as.numeric(as.factor(Period_3y)),
                               data = cluster3_ct)
print(summary(cluster3_regression))

# save a csv of cluster, decade, mean cta (days), se cta (days)
cluster1.n <- length(unique(cluster1.table$Taxa))
cluster1_cta <- cluster1.table %>%
  ddply(.(Period_3y), summarise,
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(cluster1.n),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Period_3y) %>%
  mutate(cluster = "Epipelagic") %>%
  select(cluster, Period_3y, mean_cta_days, se_cta_days) %>%
  distinct()

cluster2.n <- length(unique(cluster2.table$Taxa))
cluster2_cta <- cluster2.table %>%
  ddply(.(Period_3y), summarise,
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(cluster2.n),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Period_3y) %>%
  mutate(cluster = "Mesopelagic") %>%
  select(cluster, Period_3y, mean_cta_days, se_cta_days) %>%
  distinct()

cluster3.n <- length(unique(cluster3.table$Taxa))
cluster3_cta <- cluster3.table %>%
  ddply(.(Period_3y), summarise,
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(cluster3.n),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Period_3y) %>%
  mutate(cluster = "Demersal") %>%
  select(cluster, Period_3y, mean_cta_days, se_cta_days) %>%
  distinct()

cluster.plots <- rbind(cluster1_cta, cluster2_cta, cluster3_cta)
write.csv(cluster.plots, "nhl-cluster-plots.csv", row.names = F)