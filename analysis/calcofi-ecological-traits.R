# Script Header -----------------------------------------------------------
# Topic: CalCOFI Ichthyoplankton Ecological Traits Analysis
# Author: Kathryn Chen
# Contact: ksc005@ucsd.edu
# Date: November 2024

# Objective: This script takes a table of CalCOFI species ecological traits, performs principal components analysis (PCA) to identify major modes of variability in those traits, assesses trends in central tendency for species clustering closely together in the PCA, and outputs three files used to create manuscript figures.

# Set Up ------------------------------------------------------------------
# location declaration
here::i_am("analysis/calcofi-ecological-traits.R")

# load libraries
library(plyr)
library(tidyverse)
library(reshape2)
library(readxl)
library(fastDummies)
library(here)

# load data
eco_traits <- read.csv(here("data", "calcofi-eco-traits.csv"),
                       stringsAsFactors = T)
analysis_summary <- read.csv(here("data", "calcofi-ct-summary.csv"))

# Processing Trait Data ----------------------------------------------------
eco_traits <- eco_traits[,c(1:12)] %>%
  dplyr::rename(distribution = cross_shore_distribution)

order_count <- eco_traits %>%
  count(order) %>%
  filter(n >= 3)

eco_traits$order <- as.factor(ifelse(eco_traits$order %in% order_count$order,
                                     as.character(eco_traits$order),
                                     "Misc"))

pca_data <- eco_traits %>%
  dummy_cols(select_columns = c("adult_habitat", "distribution",
                                "biogeographic_affinity",
                                "fishing_status", "order"),
             remove_selected_columns =  T) %>%
  mutate(month_maximum = match(month_maximum, month.name),
         fishing_status_Unfished = NULL,
         r = NULL,
         mean_ct = NULL,
         sd_ct = NULL,
         trend_in_phenology = NULL,
         order_Misc = NULL,
         order_Ophidiiformes = NULL,
         order_Beloniformes = NULL,
         `order_Carangaria/misc` = NULL,
         order_Carangiformes = NULL,
         order_Centrarchiformes = NULL,
         order_Clupeiformes = NULL,
         order_Gadiformes = NULL,
         `order_Ovalentaria/misc` = NULL,
         `order_Perciformes/Cottoidei` = NULL) %>%
  dplyr::rename(Demersal = adult_habitat_Demersal,
                Epipelagic = adult_habitat_Epipelagic,
                Mesopelagic = adult_habitat_Mesopelagic,
                Coastal = distribution_Coastal,
                `Coastal-Oceanic` = `distribution_Coastal-Oceanic`,
                Oceanic = distribution_Oceanic,
                `Cool-water` = `biogeographic_affinity_Cool-water`,
                `Warm-water` = `biogeographic_affinity_Warm-water`,
                `Wide Distribution` = `biogeographic_affinity_Wide distribution`,
                Fished = fishing_status_Fished,
                Argentiniformes = order_Argentiniformes,
                Myctophiformes = order_Myctophiformes,
                `Perciformes/Scorpaenoidei` = `order_Perciformes/Scorpaenoidei`,
                Pleuronectiformes = order_Pleuronectiformes,
                Scombriformes = order_Scombriformes,
                Stomiiformes = order_Stomiiformes)

write.csv(pca_data, "pca-data.calcofi.csv", row.names = F)

# Perform PCA -----------------------------------------------------------------
drop.pca.df <- pca_data[,-c(1,3,13)]

drop.scaled_pca.df <- as.data.frame(scale(drop.pca.df))

drop.fit <- princomp(drop.scaled_pca.df, cor=TRUE)

summary(drop.fit) # print variance accounted for. the first two PCs account for 49.49% of the variance in this dataset

# create data frame with scores and loadings
drop.scores = as.data.frame(drop.fit$scores)
drop.loadings = as.data.frame(drop.fit$loadings[,1:2])

# Interpretation of visualizations identifies 3 clusters of species; assign cluster numbers.
set.seed(70524)
drop.fit_kmeans <- kmeans(drop.scores, 3, nstart = 25) # the cluster number will be contained in the vector fit_kmeans$cluster

# save needed items for figures
drop.markers <- cbind(drop.scores, eco_traits$trend_in_phenology) %>%
  dplyr::rename(trend = `eco_traits$trend_in_phenology`) %>%
  mutate(direction = as.factor(ifelse(trend > 0, "positive", "negative")))

write.csv(drop.markers, here("data", "calcofi-pca-markers.csv"), row.names = F)

drop.correlations3 = as.data.frame(cor(drop.scaled_pca.df, drop.fit$scores))
write.csv(drop.correlations3, here("data", "calcofi-pca-plots.csv"), row.names = T)

# summary scores
drop.scores$scientific_name = pca_data$scientific_name
names(drop.scores) <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9",
                        "PC10", "PC11","PC12","PC13","PC14", "PC15", "PC16","SciName")


# Species Clustering Closely Together ---------------------------------------
# Group the species clustering closely together and regress their CTa against time.
spp_cta <- analysis_summary

kept_for_analysis <- unique(as.character(eco_traits$scientific_name))

spp_cta <- spp_cta %>%
  mutate(scientific_name = as.character(scientific_name)) %>%
  filter(scientific_name %in% kept_for_analysis)

spp_cta$ID <- as.numeric(as.factor(spp_cta$scientific_name))

spp_cta <- spp_cta %>%
  mutate(Decade = as.numeric(as.factor(Decade),
                                 levels = c(1955, 1965, 1975, 1985,
                                            1995, 2005, 2015)))

# attach the clusters to the species names; i.e., which species are in which clusters?
drop.clusters <- data.frame(index = 1:57,
                            cluster = drop.fit_kmeans$cluster)

spp_cta$drop.cluster <- drop.clusters$cluster[match(spp_cta$ID, drop.clusters$index)]

# one cluster: primarily negative loadings on PC2. Scombriformes; Coastal-oceanic, epipelagic, wide distribution, warm-water. 9 spp
drop.cluster1 <- unique(spp_cta$scientific_name[spp_cta$drop.cluster == 1])

drop.cluster1.table <- spp_cta %>%
  filter(drop.cluster == 1)
drop.cluster1.group <- drop.cluster1.table %>%
  group_by(scientific_name) %>%
  count(group)
table(drop.cluster1.group$group) # 3 spp shifting earlier, 2 later, 4 no change

drop.cluster1_scores <- drop.scores %>%
  dplyr::select(PC1, PC2, SciName) %>%
  dplyr::filter(SciName %in% drop.cluster1)

drop.cluster1_ct <- spp_cta %>%
  filter(scientific_name %in% drop.cluster1)

drop.cluster1_regression <- lm(ct_anomaly_days ~ Decade,
                               data = drop.cluster1_ct)
print(summary(drop.cluster1_regression))
#y = 3.3339 - 0.8301x, F = 0.3616, df = 60, p = 0.5499

drop.set1 <- eco_traits %>%
  filter(scientific_name %in% drop.cluster1)


# cluster 2: primarily positive loadings on both PC1 and PC2. Pleuronectiformes, Scorpaeniformes; coastal, demersal, cool-water. 25 spp
drop.cluster2 <- unique(spp_cta$scientific_name[spp_cta$drop.cluster == 2])

drop.cluster2.table <- spp_cta %>%
  filter(drop.cluster == 2)
drop.cluster2.group <- drop.cluster2.table %>%
  group_by(scientific_name) %>%
  count(group)
table(drop.cluster2.group$group) # 10 spp shifting earlier, 5 later, 10 no change

drop.cluster2_scores <- drop.scores %>%
  select(PC1, PC2, SciName) %>%
  filter(SciName %in% drop.cluster2)

drop.cluster2_ct <- spp_cta %>%
  filter(scientific_name %in% drop.cluster2)

drop.cluster2_regression <- lm(ct_anomaly_days ~ as.numeric(as.factor(Decade)),
                               data = drop.cluster2_ct)
print(summary(drop.cluster2_regression))
#y = 6.904 - 1.695x, F = 3.979, df = 164, p = 0.04774

drop.set2 <- eco_traits %>%
  filter(scientific_name %in% drop.cluster2)


# cluster 3: primarily negative loadings on PC1 and positive on PC2. Myctophiformes, Stomiiformes, Argentiniformes; Mesopelagic, oceanic distribution. 23 spp
drop.cluster3 <- unique(spp_cta$scientific_name[spp_cta$drop.cluster == 3])

drop.cluster3.table <- spp_cta %>%
  filter(drop.cluster == 3)
drop.cluster3.group <- drop.cluster3.table %>%
  group_by(scientific_name) %>%
  count(group)
table(drop.cluster3.group$group) # 14 spp shifting earlier, 2 shifting later, 7 no change

drop.cluster3_scores <- drop.scores %>%
  select(PC1, PC2, SciName) %>%
  filter(SciName %in% drop.cluster3)

drop.cluster3_ct <- spp_cta %>%
  filter(scientific_name %in% drop.cluster3)

drop.cluster3_regression <- lm(ct_anomaly_days ~ as.numeric(as.factor(Decade)),
                               data = drop.cluster3_ct)
print(summary(drop.cluster3_regression))
#y = 10.1094 - 2.4590x, F = 9.394, df = 151, p = 0.002579

# look at what species are in this group
drop.set3 <- eco_traits %>%
  filter(scientific_name %in% drop.cluster3)

# save a csv of cluster, decade, mean cta (days), se cta (days)
drop.cluster1.n <- length(unique(drop.cluster1.table$scientific_name))
drop.cluster1_cta <- drop.cluster1.table %>%
  ddply(.(Decade), summarise,
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(drop.cluster1.n),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Decade) %>%
  mutate(cluster = "Epipelagic") %>%
  select(cluster, Decade, mean_cta_days, se_cta_days) %>%
  distinct()

drop.cluster2.n <- length(unique(drop.cluster2.table$scientific_name))
drop.cluster2_cta <- drop.cluster2.table %>%
  ddply(.(Decade), summarise,
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(drop.cluster2.n),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Decade) %>%
  mutate(cluster = "Demersal") %>%
  select(cluster, Decade, mean_cta_days, se_cta_days) %>%
  distinct()

drop.cluster3.n <- length(unique(drop.cluster3.table$scientific_name))
drop.cluster3_cta <- drop.cluster3.table %>%
  ddply(.(Decade), summarise,
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(drop.cluster3.n),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Decade) %>%
  mutate(cluster = "Mesopelagic") %>%
  select(cluster, Decade, mean_cta_days, se_cta_days) %>%
  distinct()

cluster.plots <- rbind(drop.cluster1_cta, drop.cluster2_cta, drop.cluster3_cta)
write.csv(cluster.plots, here("data", "calcofi-cluster-plots.csv"), row.names = F)

# More Data for Figures ------------------------------------------------------
summary(eco_traits)
new.design <- eco_traits %>%
  select(order, adult_habitat, distribution, biogeographic_affinity,
         fishing_status, r, trend_in_phenology)
summary(new.design)

# get trend in d/y
new.design$trend_in_phenology <- new.design$trend_in_phenology / 10

# 3 minimum species per category, so we are dropping Misc

new.design.min3 <- new.design %>%
  mutate(order = as.character(order),
         order = ifelse(order == "Misc",
                        NA,
                        order),
         order = as.factor(order))

write.csv(new.design.min3, "design.spp.calcofi.csv")