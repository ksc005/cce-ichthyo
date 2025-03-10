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
ichthyo3yAvgs <- ichthyo3yAvgs %>%
  select(New_Month, Taxa, Period_3y, Density_Avg_3y) %>%
  unique()

ichthyo3yAvgs <- ichthyo3yAvgs %>%
  filter(Period_3y < 8) # restrict analysis to only 1997-2017. comment out if doing full time series

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

  correlations$trend_dy <- correlations$trend_d3y / 3
}

# write.csv(correlations, here("data", "nhl-ct-summary.csv"), row.names = F)

#write.csv(correlations, here("data", "comp-trends-nhl.csv"), row.names = F)


# Group Mean, SE, Range of Trend --------------------------------------------
e <- correlations %>%
  filter(group == "earlier") %>%
  select(Taxa, cor, group, trend_d3y) %>%
  unique()

mean_se(e$trend_d3y)[1] / 3
(mean_se(e$trend_d3y)[3] - mean_se(e$trend_d3y)[1]) / 3
max(e$trend_d3y) / 3
min(e$trend_d3y) / 3

l <- correlations %>%
  filter(group == "later") %>%
  select(Taxa, cor, group, trend_d3y) %>%
  unique()

mean_se(l$trend_d3y)[1] / 3
(mean_se(l$trend_d3y)[3] - mean_se(l$trend_d3y)[1]) / 3
max(l$trend_d3y) / 3
min(l$trend_d3y) / 3

n <- correlations %>%
  filter(group == "no_change") %>%
  select(Taxa, cor, group, trend_d3y) %>%
  unique()

mean_se(n$trend_d3y)[1] / 3
(mean_se(n$trend_d3y)[3] - mean_se(n$trend_d3y)[1]) / 3
max(n$trend_d3y) / 3
min(n$trend_d3y) / 3

# Data for Figures ------------------------------------------------------------
earlierPanel <- earlierTaxa %>%
  ddply(.(Period_3y), summarise,
        mean_ct = mean(ct),
        mean_cta_months = mean(ct_anomaly_months),
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(25),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Period_3y) %>%
  distinct() %>%
  mutate(group = "earlier")

noChangePanel <- noChangeTaxa %>%
  ddply(.(Period_3y), summarise,
        mean_ct = mean(ct),
        mean_cta_months = mean(ct_anomaly_months),
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(25),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Period_3y) %>%
  distinct() %>%
  mutate(group = "no_change")

laterPanel <- laterTaxa %>%
  ddply(.(Period_3y), summarise,
        mean_ct = mean(ct),
        mean_cta_months = mean(ct_anomaly_months),
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(25),
        se_cta_days = se_cta_months * 30.44) %>%
  arrange(Period_3y) %>%
  distinct() %>%
  mutate(group = "later")

assemblagePanel <- speciesCTa %>%
  ddply(.(Period_3y), dplyr::summarise,
        Period_3y = Period_3y,
        mean_ct = mean(ct),
        mean_cta_months = mean(ct_anomaly_months),
        mean_cta_days = mean(ct_anomaly_days),
        sd_cta_months = sd(ct_anomaly_months),
        se_cta_months = sd_cta_months/sqrt(25),
        se_cta_days = se_cta_months * 30.44) %>%
  unique() %>%
  mutate(group = "assemblage")

# nhl.fig1 <- rbind(assemblagePanel,
#                       earlierPanel,
#                       laterPanel,
#                       noChangePanel)
#
# write.csv(nhl.fig1,
#           here("data", "fig1.nhl.csv"))

nhl.2017 <- rbind(assemblagePanel,
                  earlierPanel,
                  laterPanel,
                  noChangePanel)

write.csv(nhl.2017,
          here("data", "nhl-3y.csv"))

## 1.17.25 ASIDE: plotting NHL 1997-2017 only ----------------------------
my.ggplot.theme <- list(
  theme_bw(),
  theme(axis.text = element_text(size = 15, family = "Helvetica"),
        axis.title = element_text(size = 16, family = "Helvetica"),
        axis.ticks = element_line(linewidth = 1.2),

        legend.position = "none",

        panel.border = element_rect(linewidth = 1.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        strip.background = element_rect(linewidth = 1.5),
        strip.text = element_text(size = 16, family = "Helvetica"),

        plot.title = element_text(size = 18, family = "Helvetica"))
)


p1 <- assemblagePanel %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-50, 50),
                     breaks = c(-30, 0, 30)) +
  scale_x_continuous(limits = c(.4, 9.9),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "NH Line Assemblage",
       x = "3-y Period",
       y = "CTa (days)")

p2 <- earlierPanel %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-90, 90),
                     breaks = c(-60, 0, 60)) +
  scale_x_continuous(limits = c(.4, 9.9),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "Shifting Earlier",
       x = "3-y Period",
       y = "CTa (days)")

p3 <- noChangePanel %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-90, 90),
                     breaks = c(-60, 0, 60)) +
  scale_x_continuous(limits = c(.4, 9.9),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "No Linear Change",
       x = "3-y Period",
       y = "CTa (days)")

p4 <- laterPanel %>%
  ggplot(aes(x = as.numeric(as.character(Period_3y)),
             y = mean_cta_days)) +
  geom_point(size = 2.5) +
  geom_path(linewidth = 1.2) +
  geom_linerange(aes(ymin = mean_cta_days - se_cta_days,
                     ymax = mean_cta_days + se_cta_days),
                 linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  my.ggplot.theme +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.6)) +
  scale_y_continuous(limits = c(-90, 90),
                     breaks = c(-60, 0, 60)) +
  scale_x_continuous(limits = c(.4, 9.9),
                     breaks = seq(from = .5, to = 9.5, by = 1),
                     labels = seq(from = 1997, to = 2024, by = 3)) +
  labs(title = "Shifting Later",
       x = "3-y Period",
       y = "CTa (days)")

plot <- ggarrange(p1, p2, p3, p4,
                      labels = "AUTO",
                      font.label = list(size = 20),
                      hjust = -1,
                      vjust = 1.25,
                      ncol = 2, nrow = 2, align = "hv")

# ggsave(filename = "nhl_1997_2017.png",
#        plot = plot,
#        path = "~/Desktop", height = 13, width = 12, units = "in")
