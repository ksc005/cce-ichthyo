# Script Header -----------------------------------------------------------
# Topic: CalCOFI Ichthyoplankton Central Tendency Analysis
# Author: Kathryn Chen
# Contact: ksc005@ucsd.edu
# Date: November 2024

# Objective: This script takes a data file of CalCOFI ichthyoplankton monthly abundances from 1951-2022 and outputs a data file with species' central tendencies (CT), central tendency anomalies (CTa), Pearson's correlations, phenology change groups, and estimated trends in phenology. It also includes linear regression of CTa against time for phenology change groups and the entire ichthyoplankton assemblage.

# Set Up ------------------------------------------------------------------
# location declaration
here::i_am("analysis/calcofi-ichthyo-ct.R")

# load libraries
library(plyr)
library(tidyverse)
library(lmtest)
library(here)

# load data
ichthyoDecadalAvgs <- read.csv(here("data", "calcofi-ichthyo-decadal.csv"))

# Calculate Species' CT, CTa ----------------------------------------------
# For each species, take their decadal average monthly abundances and calculate CT and CTa.
speciesCTa <- ichthyoDecadalAvgs %>%
  group_by(scientific_name, Decade) %>%
  mutate(New_Month = as.numeric(New_Month),
         Numerator = New_Month * Avg_Abundance) %>%
  ddply(.(Decade, scientific_name), summarise,
        Numerator = sum(Numerator),
        Denominator = sum(Avg_Abundance),
        ct = Numerator / Denominator) %>%
  na.omit() %>%
  ddply(.(scientific_name), reframe,
        Decade = Decade,
        ct = ct,
        mean_ct_months = mean(ct),
        sd_ct_months = sd(ct),
        ct_anomaly_months = ct - mean_ct_months,
        ct_anomaly_days = ct_anomaly_months * 30.44)

# Assemblage CTa x Time ---------------------------------------------------
# For the entire assemblage, regress CTa against time, check for autocorrelation, and print the results.
decades_assemblage_regr <- lm(ct_anomaly_days ~ as.numeric(as.factor(Decade)),
                              data = speciesCTa)
dwtest(decades_assemblage_regr) # no autocorrelation
print(summary(decades_assemblage_regr))

# Classify Phenology Change Groups ----------------------------------------
# For each species, run a Pearson's correlation between their CTa and time. Classify species with Pearson's r ≥ 0.33 as later/delaying, r ≥ 0.33 as earlier/advancing, and those intermediate as having no long-term linear change.
correlations <- speciesCTa %>%
  group_by(scientific_name) %>%
  mutate(cor = cor(ct_anomaly_days,
                   as.numeric(as.character(Decade))),
         group = if(any(cor >= 0.3333)) {
           "later"
         } else{
           if(any(cor <= -0.3333)) {
             "earlier"
           } else{
             "no_change"
           }
         })

earlierTaxa <- correlations %>%
  filter(group == "earlier")
countEarlierTaxa <- length(unique(earlierTaxa$scientific_name))

laterTaxa <- correlations %>%
  filter(group == "later")
countLaterTaxa <- length(unique(laterTaxa$scientific_name))

noChangeTaxa <- correlations %>%
  filter(group == "no_change")
countNoChangeTaxa <- length(unique(noChangeTaxa$scientific_name))

# Group CTa x Time --------------------------------------------------------
# For each phenology change group, regress their CTa against time, check for autocorrelation, and print the results.
earlier_regr_3y <- lm(ct_anomaly_days ~ as.numeric(as.factor(Decade)),
                      data = earlierTaxa)
dwtest(earlier_regr_3y)
print(summary(earlier_regr_3y))

later_regr_3y <- lm(ct_anomaly_days ~ as.numeric(as.factor(Decade)),
                    data = laterTaxa)
dwtest(later_regr_3y)
print(summary(later_regr_3y))

no_change_regr_3y <- lm(ct_anomaly_days ~ as.numeric(as.factor(Decade)),
                        data = noChangeTaxa)
dwtest(no_change_regr_3y)
print(summary(no_change_regr_3y))

# Species' CTa x Time -----------------------------------------------------
# For each species, regress their CTa against time and add the slope of the regression to a new column 'trend_ddec', estimated phenological trend in days per decade. Save a data file for use in other analyses.
for(i in unique(correlations$scientific_name)) {
  df <- correlations %>% filter(scientific_name == i)

  spp_regr <- lm(ct_anomaly_days ~ as.numeric(Decade),
                 data = df)

  correlations$trend_ddec[correlations$scientific_name == i] <- spp_regr$coefficients[2]
}

write.csv(correlations, here("data", "calcofi-ct-summary.csv"), row.names = F)