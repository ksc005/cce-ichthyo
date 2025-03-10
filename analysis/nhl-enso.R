# Script Header -----------------------------------------------------------
# Topic: Newport Hydrographic Line Ichthyoplankton-ENSO Analysis
# Author: Kathryn Chen
# Contact: ksc005@ucsd.edu
# Date: November 2024

# Objective: This script takes a data file of Newport Hydrographic Line ichthyoplankton monthly abundances from 1996-2023, bins the data by phase of the El Niño-Southern Oscillation and by species' phenology change group, and runs statistical tests comparing group central tendency anomaly (CTa) against ENSO phase.

# Set Up ------------------------------------------------------------------
# location declaration
here::i_am("analysis/nhl-enso.R")

# load libraries
library(plyr)
library(tidyverse)
library(reshape2)
library(readxl)
library(FSA)
library(here)

# load data
spp_data <- read.csv(here("data", "fall_NHL_shifted.csv"))
ct_analysis_summary <- read.csv(here("data", "nhl-ct-summary.csv"))
enso <- read_excel(here("data", "nhl-enso.xlsx"))

# ENSO Analysis -----------------------------------------------------------
# For each species, bin abundance data by ENSO phase and calculate CT and CTa.
enso <- enso %>%
  melt(id.vars = "Year",
       variable.name = "Month",
       value.name = "enso") %>%
  mutate(Month = as.numeric(Month))

ichthyo_enso <- left_join(spp_data, enso, by = c("Year" = "Year",
                                                "Month" = "Month"))

ichthyo_enso_ct <- ichthyo_enso %>%
  mutate(numerator = New_Month * Density,
         enso = as.factor(enso)) %>%
  ddply(.(enso, Taxa), reframe,
        numerator = sum(numerator),
        denominator = sum(Density),
        ct = numerator / denominator) %>%
  na.omit() %>%
  ddply(.(Taxa), reframe,
        enso = enso,
        ct = ct,
        mean_ct_months = mean(ct),
        sd_ct_months = sd(ct),
        ct_anomaly_months = ct - mean_ct_months,
        ct_anomaly_days = ct_anomaly_months * 30.44)

ichthyo_enso_ct$group <- ct_analysis_summary$group[match(ichthyo_enso_ct$Taxa, ct_analysis_summary$Taxa)]


# Statistical Comparisons -------------------------------------------------
# For each phenology change group, test if there's a significant difference in CTa according to ENSO phase. Check ANOVA assumptions: normality (Shapiro-Wilk), homoscedasticity (Bartlett's); if violated perform Kruskal-Wallis instead. Post-hoc: Tukey's HSD for ANOVA, Dunn for Kruskal-Wallis.

# earlier phenology change group
enso_earlier_taxa <- ichthyo_enso_ct %>%
  filter(group == "earlier")

shapiro.test(enso_earlier_taxa$ct_anomaly_days) # p = 0.5647, violated

kruskal.test(ct_anomaly_days ~ enso,
             data = enso_earlier_taxa) # p = 0.5326, no significant difference


# no long-term linear change group
enso_no_change_taxa <- ichthyo_enso_ct %>%
  filter(group == "no_change")

shapiro.test(enso_no_change_taxa$ct_anomaly_days) # p = 0.008119, normal

bartlett.test(ct_anomaly_days ~ enso, data = enso_no_change_taxa) # p = 0.1898, violated

kruskal.test(ct_anomaly_days ~ enso,
             data = enso_no_change_taxa) # p = 0.004192, significant difference

dunnTest(ct_anomaly_days ~ enso,
         data = enso_no_change_taxa,
         method = "holm") # significant difference for El Niño vs Neutral


# later phenology change group
enso_later_taxa <- ichthyo_enso_ct %>%
  filter(group == "later")

shapiro.test(enso_later_taxa$ct_anomaly_days) # p = 0.2102, violated

kruskal.test(ct_anomaly_days ~ enso,
             data = enso_later_taxa) # p = 0.3679, no significant difference


# Data for Figures ------------------------------------------------------------
climate.nhl <- ichthyo_enso_ct %>%
  mutate(mode = "enso") %>%
  rename(phase = enso)

write.csv(climate.nhl, here("data", "climate.nhl.csv"))