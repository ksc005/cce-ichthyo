# Script Header -----------------------------------------------------------


# Set Up ------------------------------------------------------------------
# load libraries
library(plyr)
library(tidyverse)
library(lmtest)
library(here)

# load data
ichthyo3yAvgs <- read.csv(here("nhl-ichthyo-3y.csv"))

# Calculate Species' CT, CTa ----------------------------------------------
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
assemblage_regr_3y <- lm(ct_anomaly_days ~ as.numeric(Period_3y),
                         data = speciesCTa)
dwtest(assemblage_regr_3y) # no autocorrelation
print(summary(assemblage_regr_3y))

# Classify Phenology Change Groups ----------------------------------------
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
earlier_regr_3y <- lm(ct_anomaly_days ~ as.numeric(Period_3y),
                      data = earlierTaxa)
dwtest(earlier_regr_3y) # no autocorrelation
print(summary(earlier_regr_3y))

later_regr_3y <- lm(ct_anomaly_days ~ as.numeric(Period_3y),
                    data = laterTaxa)
dwtest(later_regr_3y) # no autocorrelation
print(summary(later_regr_3y))

no_change_regr_3y <- lm(ct_anomaly_days ~ as.numeric(Period_3y),
                        data = noChangeTaxa)
dwtest(no_change_regr_3y) # no autocorrelation
print(summary(no_change_regr_3y))

# Species' CTa x Time -----------------------------------------------------
for(i in unique(correlations$Taxa)) {
  df <- correlations %>% filter(Taxa == i)

  spp_regr <- lm(ct_anomaly_days ~ Period_3y,
                 data = df)

  correlations$trend_ddec[correlations$Taxa == i] <- spp_regr$coefficients[2]

}

write.csv(correlations, "nhl-ct-summary.csv", row.names = F)
