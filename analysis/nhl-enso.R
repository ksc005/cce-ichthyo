# Script Header -----------------------------------------------------------


# Set Up ------------------------------------------------------------------
setwd("/Volumes/petrik-lab/kchen/cce-ichthyo")

# load libraries
library(plyr)
library(tidyverse)
library(reshape2)
library(readxl)
library(FSA)

# load data
ichthyo3yAvgs <- read.csv("nhl-ichthyo-3y.csv")
enso <- read_excel("nhl-enso.xlsx")


# ENSO Analysis -----------------------------------------------------------
enso <- enso %>%
  melt(id.vars = "Year",
       variable.name = "Month",
       value.name = "enso") %>%
  mutate(Month = as.numeric(Month))

# label ichthyo data with enso phase
ichthyo_enso <- left_join(ichthyo3yAvgs, enso, by = c("Old_Year" = "Year",
                                                "Old_Month" = "Month"))

# calculate central tendency for each species for each ENSO phase
ichthyo_enso_ct <- ichthyo_enso %>%
  mutate(numerator = New_Month * Density_Avg_3y,
         enso = as.factor(enso)) %>%
  ddply(.(enso, Taxa), reframe,
        numerator = sum(numerator),
        denominator = sum(Density_Avg_3y),
        ct = numerator / denominator) %>%
  na.omit() %>%
  ddply(.(Taxa), reframe,
        enso = enso,
        ct = ct,
        mean_ct_months = mean(ct),
        sd_ct_months = sd(ct),
        ct_anomaly_months = ct - mean_ct_months,
        ct_anomaly_days = ct_anomaly_months * 30.44)

# attach phenology change groups
ct_analysis_summary <- read.csv("nhl-ct-summary.csv")

ichthyo_enso_ct$group <- ct_analysis_summary$group[match(ichthyo_enso_ct$Taxa, ct_analysis_summary$Taxa)]

# statistical tests for ENSO comparisons -----------------------------------------
# Check ANOVA assumptions: normality (Shapiro-Wilk), homoscedasticity (Bartlett's). If violated perform Kruskal-Wallis instead. Post-hoc: Tukey's HSD for ANOVA, Dunn for Kruskal-Wallis.

# earlier phenology change group
enso_earlier_taxa <- ichthyo_enso_ct %>%
  filter(group == "earlier")

shapiro.test(enso_earlier_taxa$ct_anomaly_days) # p = 0.06089, violated

kruskal.test(ct_anomaly_days ~ enso,
             data = enso_earlier_taxa) # p = 0.0316, significant difference

dunnTest(ct_anomaly_days ~ enso,
         data = enso_earlier_taxa,
         method = "holm") # significant difference for El Nino vs La Nina


# no long-term linear change group
enso_no_change_taxa <- ichthyo_enso_ct %>%
  filter(group == "no_change")

shapiro.test(enso_no_change_taxa$ct_anomaly_days) # p = 0.004386, normal

bartlett.test(ct_anomaly_days ~ enso, data = enso_no_change_taxa) # p = 0.2899, violated

kruskal.test(ct_anomaly_days ~ enso,
             data = enso_no_change_taxa) # p = 0.02478, significant difference

dunnTest(ct_anomaly_days ~ enso,
         data = enso_no_change_taxa,
         method = "holm") # significant difference for El Nino vs Neutral


# later phenology change group
enso_later_taxa <- ichthyo_enso_ct %>%
  filter(group == "later")

shapiro.test(enso_later_taxa$ct_anomaly_days) # p = 0.3736, violated

kruskal.test(ct_anomaly_days ~ enso,
             data = enso_later_taxa) # p = 0.2765, no significant difference
