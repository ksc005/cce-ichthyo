# Script Header -----------------------------------------------------------
# Topic: Newport Hydrographic Line Ichthyoplankton Central Tendency Analysis
# Author: Kathryn Chen
# Contact: ksc005@ucsd.edu
# Date: November 2024

# Objective: This script takes a data file of Newport Hydrographic Line ichthyoplankton monthly abundances from 1996-2023 and outputs a data file with species' central tendencies (CT), central tendency anomalies (CTa), Pearson's correlations, phenology change groups, and estimated trends in phenology. It also includes linear regression of CTa against time for phenology change groups and the entire ichthyoplankton assemblage.

# Set Up ------------------------------------------------------------------
# location declaration
here::i_am("analysis/nhl-ichthyo-ct.R")

# load libraries
library(plyr)
library(tidyverse)
library(lmtest)
library(here)

# load data
ichthyo3yAvgs <- read.csv(here("data", "nhl-ichthyo-3y.csv"))

# Calculate Species' CT, CTa ----------------------------------------------
# For each species, take their 3-year average monthly abundances and calculate CT and CTa.
speciesCTa <- ichthyo3yAvgs %>%
  group_by(Taxa, Period_3y) %>%
  mutate(New_Month = as.numeric(New_Month),
         Numerator = New_Month * Density_Avg_3y) %>%
  ddply(.(Period_3y, Taxa), summarise,
        Numerator = sum(Numerator),
        Denominator = sum(Density_Avg_3y),
        ct = Numerator / Denominator) %>%
  na.omit() %>%
  ddply(.(Taxa), reframe,
        Period_3y = Period_3y,
        ct = ct,
        mean_ct_months = mean(ct),
        sd_ct_months = sd(ct),
        ct_anomaly_months = ct - mean_ct_months,
        ct_anomaly_days = ct_anomaly_months * 30.44)

# Assemblage CTa x Time ---------------------------------------------------
# For the entire assemblage, regress CTa against time, check for autocorrelation, and print the results.
assemblage_regr_3y <- lm(ct_anomaly_days ~ as.numeric(Period_3y),
                         data = speciesCTa)
dwtest(assemblage_regr_3y) # no autocorrelation
print(summary(assemblage_regr_3y))

# Classify Phenology Change Groups ----------------------------------------
# For each species, run a Pearson's correlation between their CTa and time. Classify species with Pearson's r ≥ 0.33 as later/delaying, r ≥ 0.33 as earlier/advancing, and those intermediate as having no long-term linear change.
correlations <- speciesCTa %>%
  group_by(Taxa) %>%
  mutate(cor = cor(ct_anomaly_days,
                   as.numeric(as.character(Period_3y))),
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
countEarlierTaxa <- length(unique(earlierTaxa$Taxa))

laterTaxa <- correlations %>%
  filter(group == "later")
countLaterTaxa <- length(unique(laterTaxa$Taxa))

noChangeTaxa <- correlations %>%
  filter(group == "no_change")
countNoChangeTaxa <- length(unique(noChangeTaxa$Taxa))

# Group CTa x Time --------------------------------------------------------
# For each phenology change group, regress their CTa against time, check for autocorrelation, and print the results.
earlier_regr_3y <- lm(ct_anomaly_days ~ as.numeric(Period_3y),
                      data = earlierTaxa)
dwtest(earlier_regr_3y)
print(summary(earlier_regr_3y))

later_regr_3y <- lm(ct_anomaly_days ~ as.numeric(Period_3y),
                    data = laterTaxa)
dwtest(later_regr_3y)
print(summary(later_regr_3y))

no_change_regr_3y <- lm(ct_anomaly_days ~ as.numeric(Period_3y),
                        data = noChangeTaxa)
dwtest(no_change_regr_3y)
print(summary(no_change_regr_3y))

# Species' CTa x Time -----------------------------------------------------
# For each species, regress their CTa against time and add the slope of the regression to a new column 'trend_d3y', estimated phenological trend in days per 3 year period. Save a data file for use in other analyses.
for(i in unique(correlations$Taxa)) {
  df <- correlations %>% filter(Taxa == i)

  spp_regr <- lm(ct_anomaly_days ~ Period_3y,
                 data = df)

  correlations$trend_d3y[correlations$Taxa == i] <- spp_regr$coefficients[2]
}

write.csv(correlations, here("data", "nhl-ct-summary.csv"), row.names = F)