# Script Header -----------------------------------------------------------
# Topic: CalCOFI Ichthyoplankton-ENSO Analysis
# Author: Kathryn Chen
# Contact: ksc005@ucsd.edu
# Date: November 2024

# Objective: This script takes a data file of CalCOFI ichthyoplankton monthly abundances from 1951-2022, bins the data by phase of the El Niño-Southern Oscillation and by species' phenology change group, and runs statistical tests comparing group central tendency anomaly (CTa) against ENSO phase.

# Set Up ------------------------------------------------------------------
# location declaration
here::i_am("analysis/calcofi-enso.R")

# load libraries
library(plyr)
library(tidyverse)
library(reshape2)
library(readxl)
library(FSA)
library(here)

# load data
spp_data <- read.csv(here("data", "CalCOFI_by_month_and_spp.csv"))
ct_analysis_summary <- read.csv(here("data", "calcofi-ct-summary.csv"))
enso <- read.csv(here("data", "calcofi-enso.csv"), stringsAsFactors = T)

# Processing Ichthyo Data ----------------------------------------------------
# remove species with insufficient data
dropped_spp <- c("Nansenia crassa", "Ophiodon elongatus",
                 "Psettichthys melanostictus", "Seriola lalandi",
                 "Hemilepidotus spinosus", "Embassichthys bathybius",
                 "Symbolophorus evermanni", "Loweina rara",
                 "Diogenichthys laternatus")

spp_data <- spp_data %>%
  filter(!scientific_name %in% dropped_spp)

# some species also had insufficient data in certain decades, drop those:
spp_data$decade <- as.factor(
  ifelse(spp_data$year >=1950 & spp_data$year <=1959, 1955,
         ifelse(spp_data$year >=1960 & spp_data$year <=1969, 1965,
                ifelse(spp_data$year >=1970 & spp_data$year <=1979, 1975,
                       ifelse(spp_data$year >=1980 & spp_data$year <=1989, 1985,
                              ifelse(spp_data$year >=1990 & spp_data$year <=1999, 1995,
                                     ifelse(spp_data$year >=2000 & spp_data$year <=2009, 2005, 2015))))))
)

spp_data <- spp_data[!(spp_data$scientific_name %in% c("Electrona risso", "Scorpaenichthys marmoratus") & spp_data$decade == 1955),]

spp_data <- spp_data[!(spp_data$scientific_name %in% c("Pseudobathylagus milleri","Loweina rara") & spp_data$decade == 1965),]

spp_data <- spp_data[!(spp_data$scientific_name == "Symphurus atricaudus" & spp_data$decade == 1975),]

spp_data <- spp_data[!(spp_data$scientific_name %in% c("Aristostomias scintillans", "Notolychnus valdiviae", "Hippoglossina stomata") & spp_data$decade %in% c(1955, 1975)),]

spp_data <- spp_data[!(spp_data$scientific_name == "Diogenichthys laternatus" & spp_data$decade %in% c(1955, 1965, 1975, 1985)),]

unique(spp_data$scientific_name) # confirm = 57 spp

# make sure to year-shift for our fall-winter spawning species
fall_winter <- c("Argentina sialis", "Aristostomias scintillans",
                 "Bathylagus pacificus", "Diogenichthys atlanticus",
                 "Diogenichthys laternatus", "Electrona risso",
                 "Embassichthys bathybius", "Hemilepidotus spinosus",
                 "Hippoglossina stomata", "Leuroglossus stilbius",
                 "Lipolagus ochotensis", "Loweina rara",
                 "Merluccius productus", "Ophiodon elongatus",
                 "Oxylebius pictus", "Pleuronichthys decurrens",
                 "Protomyctophum crockeri", "Pseudobathylagus milleri",
                 "Scorpaenichthys marmoratus", "Sebastes goodei",
                 "Sebastes jordani", "Sebastes levis",
                 "Sebastes paucispinis", "Stomias atriventer",
                 "Tetragonurus cuvieri")

shift_july <- spp_data %>%
  filter(scientific_name %in% fall_winter) %>%
  mutate(New_Year =
           ifelse(month >= 7, year, year - 1),
         New_Month =
           ifelse(month >= 7, month - 6,
                  ifelse(month < 7, month + 6, NA)),
         Year_Start = "July")
unique(shift_july$scientific_name) # confirm = 20 spp

no_shift_jan <- spp_data %>%
  filter(!scientific_name %in% fall_winter) %>%
  mutate(New_Year = year,
         New_Month = month,
         Year_Start = "January")
unique(no_shift_jan$scientific_name) # confirm = 37 spp

year_shifted <- rbind(no_shift_jan, shift_july) %>%
  arrange(scientific_name)
unique(year_shifted$scientific_name) # confirm total = 57 spp

year_shifted$decade <- factor(
  ifelse(year_shifted$year >=1950 & year_shifted$year <=1959, 1955,
         ifelse(year_shifted$year >=1960 & year_shifted$year <=1969, 1965,
                ifelse(year_shifted$year >=1970 & year_shifted$year <=1979, 1975,
                       ifelse(year_shifted$year >=1980 & year_shifted$year <=1989, 1985,
                              ifelse(year_shifted$year >=1990 & year_shifted$year <=1999, 1995,
                                     ifelse(year_shifted$year >=2000 & year_shifted$year <=2009, 2005, 2015))))))
)

year_shifted$New_Decade <- factor(
  ifelse(year_shifted$New_Year >=1950 & year_shifted$New_Year <=1959, 1955,
         ifelse(year_shifted$New_Year >=1960 & year_shifted$New_Year <=1969, 1965,
                ifelse(year_shifted$New_Year >=1970 & year_shifted$New_Year <=1979, 1975,
                       ifelse(year_shifted$New_Year >=1980 & year_shifted$New_Year <=1989, 1985,
                              ifelse(year_shifted$New_Year >=1990 & year_shifted$New_Year <=1999, 1995,
                                     ifelse(year_shifted$New_Year >=2000 & year_shifted$New_Year <=2009, 2005, 2015))))))
)

# ENSO Analysis -----------------------------------------------------------
# For each species, bin abundance data by ENSO phase and calculate CT and CTa.
enso <- enso[,c(1:13)] %>%
  melt(id.vars = "Year",
       variable.name = "Month",
       value.name = "enso") %>%
  mutate(Month = as.numeric(Month))

ichthyo_enso <- left_join(year_shifted, enso, by = c("year" = "Year",
                                                      "month" = "Month"))

ichthyo_enso_ct <- ichthyo_enso %>%
  mutate(numerator = New_Month * abundance,
         enso = as.factor(enso)) %>%
  ddply(.(enso, scientific_name), reframe,
        numerator = sum(numerator),
        denominator = sum(abundance),
        ct = numerator / denominator) %>%
  na.omit() %>%
  ddply(.(scientific_name), reframe,
        enso = enso,
        ct = ct,
        mean_ct_months = mean(ct),
        sd_ct_months = sd(ct),
        ct_anomaly_months = ct - mean_ct_months,
        ct_anomaly_days = ct_anomaly_months * 30.44)

ichthyo_enso_ct$group <- ct_analysis_summary$group[match(ichthyo_enso_ct$scientific_name,
                                                         ct_analysis_summary$scientific_name)]


# Statistical Comparisons -------------------------------------------------
# For each phenology change group, test if there's a significant difference in CTa according to ENSO phase. Check ANOVA assumptions: normality (Shapiro-Wilk), homoscedasticity (Bartlett's); if violated perform Kruskal-Wallis instead. Post-hoc: Tukey's HSD for ANOVA, Dunn for Kruskal-Wallis.

# earlier phenology change group
enso_earlier_taxa <- ichthyo_enso_ct %>%
  filter(group == "earlier")

shapiro.test(enso_earlier_taxa$ct_anomaly_days) # p = 0.1312, violated

kruskal.test(ct_anomaly_days ~ enso,
             data = enso_earlier_taxa) # p = 0.1615, insignificant


# no long-term linear change group
enso_no_change_taxa <- ichthyo_enso_ct %>%
  filter(group == "no_change")

shapiro.test(enso_no_change_taxa$ct_anomaly_days) # p = 0.0004265, normal

bartlett.test(ct_anomaly_days ~ enso, data = enso_no_change_taxa) # p = 0.04739, homoscedastic

no_change_enso_anova <- aov(ct_anomaly_days ~ enso,
                            data = enso_no_change_taxa)
summary(no_change_enso_anova) # p = 0.389, no significant difference


# later phenology change group
enso_later_taxa <- ichthyo_enso_ct %>%
  filter(group == "later")

shapiro.test(enso_later_taxa$ct_anomaly_days) # p = 0.01502, normal

bartlett.test(ct_anomaly_days ~ enso, data = enso_later_taxa) # p = 0.07913, violated

kruskal.test(ct_anomaly_days ~ enso,
             data = enso_later_taxa) # p = 0.7703, no significant difference

# Data for Figures ------------------------------------------------------------
## pull out the data needed to produce these figures
climate.calcofi <- ichthyo_enso_ct %>%
  mutate(mode = "enso") %>%
  rename(phase = enso)

write.csv(climate.calcofi,
          here("data", "climate.calcofi.csv"))