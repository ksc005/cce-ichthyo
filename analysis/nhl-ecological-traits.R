# Script Header -----------------------------------------------------------
# Topic: Newport Hydrographic Line Ichthyoplankton Ecological Traits Analysis
# Author: Kathryn Chen
# Contact: ksc005@ucsd.edu
# Date: November 2024

# Objective: This script takes a table of Newport Hydrographic Line species ecological traits, performs principal components analysis (PCA) to identify major modes of variability in those traits, assesses trends in central tendency for species clustering closely together in the PCA, and outputs three files used to create manuscript figures.

# Set Up ------------------------------------------------------------------
# location declaration
here::i_am("analysis/nhl-enso.R")

# load libraries
library(plyr)
library(tidyverse)
library(reshape2)
library(readxl)
library(fastDummies)
library(here)

# load data
eco_traits <- read_excel(here("data", "nhl-eco-traits.xlsx"), sheet = 1)
spp_cta <- read.csv(here("data", "nhl-ct-summary.csv"))

# Ecological Traits Table ----------------------------------------------------
# Process the ecological traits table, including classifying taxonomic orders with < 3 species represented as miscellaneous.
eco_traits <- eco_traits %>%
  rename(distribution = cross_shore_distribution) %>%
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
  count(order) %>%
  filter(n >= 3)

eco_traits$order <- as.factor(ifelse(eco_traits$order %in% order_count$order,
                                     as.character(eco_traits$order),
                                     "Misc"))


# Principal Components Analysis -----------------------------------------------
# Create dummy columns, remove unneeded columns, standardize the data, and perform PCA on the species' ecological traits.
pca_data <- eco_traits %>%
  select(scientific_name, month_maximum, order, adult_habitat,
         distribution, biogeographic_affinity) %>%
  arrange(scientific_name) %>%
  dummy_cols(select_columns = c("adult_habitat", "distribution",
                                "biogeographic_affinity", "order"),
             remove_selected_columns =  T) %>%
  mutate(order_Misc = NULL,
         `biogeographic_affinity_Warm-water` = NULL) %>%
  rename(Demersal = adult_habitat_Demersal,
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

pca.df <- pca_data[,-1]
scaled_pca.df <- as.data.frame(scale(pca.df))

fit <- princomp(scaled_pca.df, cor = TRUE)

summary(fit)
scores = as.data.frame(fit$scores)

# Interpretation of visualizations identifies 3 clusters of species; assign cluster numbers.
set.seed(70524)
fit_kmeans <- kmeans(scores, 3, nstart = 25)



# Save files needed for visualizing the PCA
eco_traits <- eco_traits %>%
  arrange(scientific_name)

markers <- cbind(scores, eco_traits$trend_NHL_d3y) %>%
  rename(trend = `eco_traits$trend_NHL_d3y`) %>%
  mutate(direction = as.factor(ifelse(trend > 0, "positive", "negative")))

write.csv(markers, here("data", "nhl-pca-markers.csv"), row.names = F)

correlations = as.data.frame(cor(scaled_pca.df, fit$scores))

write.csv(correlations, here("data", "nhl-pca-plots.csv"), row.names = T)

scores$scientific_name = pca_data$scientific_name
names(scores) <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9",
                   "PC10", "PC11","PC12","SciName")


# Species Clustering Closely Together ---------------------------------------
# Group the species clustering closely together and regress their CTa against time.
spp_cta <- spp_cta %>%
  group_by(Taxa) %>%
  mutate(ID = cur_group_id())

clusters <- data.frame(index = 1:length(unique(spp_cta$Taxa)),
                       cluster = fit_kmeans$cluster)

spp_cta$cluster <- clusters$cluster[match(spp_cta$ID, clusters$index)]

# Cluster 1: positive on PC1, negative on PC2. Epipelagic, Coastal-oceanic
cluster1 <- unique(spp_cta$Taxa[spp_cta$cluster == 1])

cluster1.table <- spp_cta %>%
  filter(cluster == 1)
cluster1.group <- cluster1.table %>%
  group_by(Taxa) %>%
  count(group)
table(cluster1.group$group) # 0 spp shifting earlier, 0 later, 1 no change

cluster1_scores <- scores %>%
  select(PC1, PC2, SciName) %>%
  filter(SciName %in% cluster1)

cluster1_ct <- spp_cta %>%
  filter(Taxa %in% cluster1)

cluster1_regression <- lm(ct_anomaly_days ~ as.numeric(as.character(Period_3y)),
                               data = cluster1_ct)
print(summary(cluster1_regression))

# Cluster 2: primarily positive on PC1,PC2. Mesopelagic, Oceanic, Myctophiformes
cluster2 <- unique(spp_cta$Taxa[spp_cta$cluster == 2])

cluster2.table <- spp_cta %>%
  filter(cluster == 2)
cluster2.group <- cluster2.table %>%
  group_by(Taxa) %>%
  count(group)
table(cluster2.group$group) # 0 spp shifting earlier, 0 shifting later, 4 no change

cluster2_scores <- scores %>%
  select(PC1, PC2, SciName) %>%
  filter(SciName %in% cluster2)

cluster2_ct <- spp_cta %>%
  filter(Taxa %in% cluster2)

cluster2_regression <- lm(ct_anomaly_days ~ as.numeric(as.character(Period_3y)),
                               data = cluster2_ct)
print(summary(cluster2_regression))

# The first 3-y period has abnormally low CTa for this group. Repeat the regression without it.
cluster2_ct.excluded <- cluster2_ct %>%
  filter(Period_3y != 1)

cluster2_regression.excluded <- lm(ct_anomaly_days ~ as.numeric(as.factor(Period_3y)),
                                        data = cluster2_ct.excluded)
print(summary(cluster2_regression.excluded))

# Cluster 3: primarily negative on PC1. Pleuronectiformes, P.Cottoidei, Coastal, Demersal, Cool-water
cluster3 <- unique(spp_cta$Taxa[spp_cta$cluster == 3])

cluster3.table <- spp_cta %>%
  filter(cluster == 3)
cluster3.group <- cluster3.table %>%
  group_by(Taxa) %>%
  count(group)
table(cluster3.group$group) # 5 spp shifting earlier, 1 later, 14 no change

cluster3_scores <- scores %>%
  select(PC1, PC2, SciName) %>%
  filter(SciName %in% cluster3)

cluster3_ct <- spp_cta %>%
  filter(Taxa %in% cluster3)

cluster3_regression <- lm(ct_anomaly_days ~ as.numeric(as.factor(Period_3y)),
                               data = cluster3_ct)
print(summary(cluster3_regression))

# Save a csv of cluster, decade, mean cta (days), se cta (days) for a figure.
cluster1_cta <- cluster1.table %>%
  ddply(.(Period_3y), summarise,
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(length(unique(cluster1.table$Taxa))),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Period_3y) %>%
  mutate(cluster = "Epipelagic") %>%
  select(cluster, Period_3y, mean_cta_days, se_cta_days) %>%
  distinct()

cluster2_cta <- cluster2.table %>%
  ddply(.(Period_3y), summarise,
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(length(unique(cluster2.table$Taxa))),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Period_3y) %>%
  mutate(cluster = "Mesopelagic") %>%
  select(cluster, Period_3y, mean_cta_days, se_cta_days) %>%
  distinct()

cluster3_cta <- cluster3.table %>%
  ddply(.(Period_3y), summarise,
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(length(unique(cluster3.table$Taxa))),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Period_3y) %>%
  mutate(cluster = "Demersal") %>%
  select(cluster, Period_3y, mean_cta_days, se_cta_days) %>%
  distinct()

cluster.plots <- rbind(cluster1_cta, cluster2_cta, cluster3_cta)
write.csv(cluster.plots, here("data", "nhl-cluster-plots.csv"), row.names = F)

# More Data for Figures -------------------------------------------------------
new.design <- eco_traits %>%
  select(order, adult_habitat, distribution, biogeographic_affinity, fishing_status, cor_NHL, trend_NHL_d3y)

# get trend in d/y
new.design$trend_NHL_d3y <- new.design$trend_NHL_d3y / 3

# 3 minimum species per category, so we are dropping: Clupeiformes, Misc, Perciformes/Zoarcoidei, Warm-water

new.design.min3 <- new.design %>%
  mutate(order = as.character(order),
         order = ifelse(order %in%
                          c("Clupeiformes", "Misc", "Perciformes/Zoarcoidei"),
                        NA,
                        order),
         order = as.factor(order)) %>%
  mutate(biogeographic_affinity = as.character(biogeographic_affinity),
         biogeographic_affinity = ifelse(biogeographic_affinity == "Warm-water",
                                         NA,
                                         biogeographic_affinity),
         biogeographic_affinity = as.factor(biogeographic_affinity))

write.csv(new.design.min3, "design.spp.nhl.csv")