# Script Header -----------------------------------------------------------
# Topic: Newport Hydrographic Line Local Oceanic Variables Central Tendency Analysis
# Author: Kathryn Chen
# Contact: ksc005@ucsd.edu
# Date: November 2024

# Objective: This script takes (raw and 3-y averaged) data files of Newport Hydrographic Line upper ocean (100m) temperature, northern and southern copepod biomass, and positive values of the Bakun coastal upwelling index at 45˚N, 125˚W. For each variable, it calculates central tendency (CT) and central tendency anomaly (CTa) with associated standard errors and regresses CTa against time. CT and CTa are calculated three times, where month 1 is set as: January, July, and October (to match ichthyoplankton species whose CT analyses were shifted to begin in these months), and the script outputs data files for each (3 starting months x 3 variables = 9 output data files).

# Set Up ------------------------------------------------------------------
# location declaration
here::i_am("analysis/nhl-enso.R")

# load libraries
library(plyr)
library(tidyverse)
library(reshape2)
library(here)

# Temperature Analysis -------------------------------------------------------
# load data
temp <- read.csv(here("data", "nhl-temp.csv"))
temp_3y <- read.csv(here("data", "nhl-temp-3y.csv"))

## CT/CTa when month 1 = January ----------------------------------------------
# subtract the rounded minimum from all values:
temp_3y <- temp_3y %>%
  mutate(adjusted_t = avg_t - floor(min(temp_3y$avg_t)))

# calculate ct and cta:
temp_ct <- temp_3y %>%
  group_by(Period_3y) %>%
  summarise(ct = (sum(Month * adjusted_t) / sum(adjusted_t))) %>%
  ungroup() %>%
  mutate(mean_ct_months = mean(ct),
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44)

# ct_function <- function(temp_3y, i){
#   x <- temp_3y[i,] %>%
#     group_by(Period_3y, Month) %>%
#     reframe(adjusted_t = avg_t - floor(min(temp_3y$avg_t))) %>%
#     group_by(Period_3y) %>%
#     reframe(ct = (sum(Month * adjusted_t) / sum(adjusted_t)))
#   return(x$ct)
# }
#
# set.seed(123)
# temp_ct_boot = by(temp_3y, temp_3y$Period_3y,
#                   function(i)boot(i, ct_function, R = 1000))
# temp_ct_boot

# from the print output of the bootstrap, we get the standard errors:
temp_ct_se <- c(0.1440763, 0.1230981, 0.06477802,
               0.1026722, 0.1071475, 0.1401642,
               0.09183508, 0.08943311, 0.09143974)

temp_ct <- temp_ct %>%
  mutate(se_ct_months = temp_ct_se,
         lower_ct = ct - 1.96 * se_ct_months,
         upper_ct = ct + 1.96 * se_ct_months)

# cta_function <- function(temp_3y, i){
#   x <- temp_3y[i,] %>%
#     group_by(Period_3y, Month) %>%
#     reframe(adjusted_t = avg_t - floor(min(temp_3y$avg_t))) %>%
#     group_by(Period_3y) %>%
#     reframe(ct = (sum(Month * adjusted_t) / sum(adjusted_t)))
#   y <- x %>%
#     mutate(mean_ct_months = mean(ct),
#            ct_anomaly_months = ct - mean_ct_months,
#            ct_anomaly_days = ct_anomaly_months * 30.44)
#   return(y$ct_anomaly_days)
# }
#
# set.seed(124)
# temp_cta_boot <- boot(temp_3y, cta_function, R = 1000)
# temp_cta_boot

# from the print output, we get the standard errors of CTa:
temp_cta_se <- c(3.644214, 1.381582, 1.373281,
                1.493738, 1.462603, 2.504765,
                2.615278, 1.680094, 1.655903)

# add the CTa SEs to the dataframe:
temp_ct$se_cta_days <- temp_cta_se

# match labels for plots:
period_years <- tibble(year = seq(from = 1997, to = 2021, by = 3),
                       period = 1:9)
temp_ct$Year <- period_years$year[match(temp_ct$Period_3y,
                                        period_years$period)]

write.csv(temp_ct, here("data", "nhl-temp-ct.csv"), row.names = F)

## Regression CTa x Time ----------------------------------------------------
temp_regression <- lm(ct_anomaly_months ~ Year,
                     data = temp_ct)
print(summary(temp_regression))

## CT/CTa when month 1 = October or July -------------------------------------
# shift T to start in Oct
temp_shift_oct <- temp %>%
  mutate(New_Year =
           ifelse(Month >= 10, Year, Year - 1),
         New_Month =
           ifelse(Month >= 10, Month - 9,
                  ifelse(Month < 10, Month + 3, NA)),
         Year_Start = "October") # assign new month and year values

temp_shift_oct2 <- temp_shift_oct %>%
  mutate(Period_3y = period_years$period[match(temp_shift_oct$New_Year,
                                             period_years$year)]) # match shifted years to periods

temp_shift_oct_3y <- temp_shift_oct2 %>%
  group_by(Period_3y, Month) %>%
  summarise(avg_t = mean(Upper100mTemp)) # average over periods

temp_shift_oct_3y <- temp_shift_oct_3y %>%
  mutate(adjusted_t = avg_t - floor(min(temp_shift_oct_3y$avg_t))) # subtract rounded minimum value

temp_shift_oct_3y <- temp_shift_oct_3y %>%
  group_by(Period_3y) %>%
  summarise(ct = (sum(Month * adjusted_t) / sum(adjusted_t))) %>%
  ungroup() %>%
  mutate(mean_ct_months = mean(ct),
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44) # get ct and cta

write.csv(temp_shift_oct_3y, here("data", "nhl-temp-oct-ct.csv"), row.names = F)

# shift T to start in July
temp_shift_july <- temp %>%
  mutate(New_Year =
           ifelse(Month >= 7, Year, Year - 1),
         New_Month =
           ifelse(Month >= 7, Month - 6,
                  ifelse(Month < 7, Month + 6, NA)),
         Year_Start = "July")

temp_shift_july2 <- temp_shift_july %>%
  mutate(Period_3y = period_years$period[match(temp_shift_july$New_Year,
                                             period_years$year)]) # match shifted years to periods

temp_shift_july_3y <- temp_shift_july2 %>%
  group_by(Period_3y, Month) %>%
  summarise(avg_t = mean(Upper100mTemp)) # average over periods

temp_shift_july_3y <- temp_shift_july_3y %>%
  mutate(adjusted_t = avg_t - floor(min(temp_shift_july_3y$avg_t))) # subtract rounded minimum value

temp_shift_july_3y <- temp_shift_july_3y %>%
  group_by(Period_3y) %>%
  summarise(ct = (sum(Month * adjusted_t) / sum(adjusted_t))) %>%
  ungroup() %>%
  mutate(unshift_ct = ct - 6,
         mean_ct_months = mean(ct),
         unshift_mean_ct = mean_ct_months - 6,
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44) #g et ct and cta

write.csv(temp_shift_july_3y, here("data", "nhl-temp-july-ct.csv"), row.names = F)


# Copepod Analysis --------------------------------------------------------
# load data
cope <- read.csv(here("data", "nhl-copepods.csv"))
cope_3y <- read.csv(here("data", "nhl-copepods-3y.csv"))

## CT/CTa when month 1 = January ----------------------------------------------
cope_ct <- cope_3y %>%
  group_by(Period_3y, Assemblage) %>%
  summarise(ct = (sum(Month * avg_Cope) / sum(avg_Cope))) %>%
  ungroup() %>%
  group_by(Assemblage) %>%
  mutate(mean_ct_months = mean(ct),
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44)

# ct_function <- function(cope_3y, i){
#   x <- cope_3y[i,] %>%
#     group_by(Period_3y, Assemblage, Month) %>%
#     reframe(avg_Cope = mean(Log10Biomass)) %>%
#     group_by(Period_3y, Assemblage) %>%
#     reframe(ct = (sum(Month * avg_Cope) / sum(avg_Cope)))
#   return(x$ct)
# }
#
# set.seed(124)
# cope_ct_boot = by(cope_3y, cope_3y[,c(3,5)],
#                   function(i)boot(i, ct_function, R = 1000))
# cope_ct_boot

# get the standard errors:
n_cope_ct_se <- data.frame(se = c(0.2688833, 0.1294764, 0.1607315,
                                  0.1679882, 0.176898, 0.1397594,
                                  0.2507668, 0.1343653, 0.1597159),
                           Period_3y = 1:9)
s_cope_ct_se <- data.frame(se = c(0.3686187, 0.3655488, 0.1956784,
                                  0.409006, 0.2682012, 0.2803371,
                                  0.1469943, 0.3432209, 0.3328616),
                           Period_3y = 1:9)

cope_ct <- cope_ct %>%
  mutate(se_ct_months = ifelse(Assemblage == "Northern",
                               n_cope_ct_se$se[match(Period_3y, n_cope_ct_se$Period_3y)],
                               s_cope_ct_se$se[match(Period_3y, s_cope_ct_se$Period_3y)]),
         lower_ct = ct - 1.96 * se_ct_months,
         upper_ct = ct + 1.96 * se_ct_months)

# cta_function <- function(cope_3y, i){
#   x <- cope_3y[i,] %>%
#     group_by(Period_3y, Assemblage, Month) %>%
#     reframe(avg_Cope = mean(Log10Biomass)) %>%
#     group_by(Period_3y, Assemblage) %>%
#     reframe(ct = (sum(Month * avg_Cope) / sum(avg_Cope)))
#   y <- x %>%
#     group_by(Assemblage) %>%
#     mutate(mean_ct_months = mean(ct),
#            ct_anomaly_months = ct - mean_ct_months,
#            ct_anomaly_days = ct_anomaly_months * 30.44)
#   return(y$ct_anomaly_days)
# }
#
# set.seed(126)
# cope_cta_boot <- boot(cope_3y, cta_function, R = 1000)
# cope_cta_boot

n_cope_cta_se <- data.frame(se = c(7.369330, 4.268114, 4.704061,
                                   4.873223, 4.874256, 4.278765,
                                   6.830186, 4.146559, 4.720907),
                            Period_3y = 1:9)
s_cope_cta_se <- data.frame(se = c(10.537443, 10.502526, 6.185150,
                                   11.570690, 8.126599, 8.756003,
                                   4.969566, 9.323560, 10.237621),
                            Period_3y = 1:9)

cope_ct <- cope_ct %>%
  mutate(se_cta_days = ifelse(Assemblage == "Northern",
                              n_cope_cta_se$se[match(Period_3y, n_cope_cta_se$Period_3y)],
                              s_cope_cta_se$se[match(Period_3y, s_cope_cta_se$Period_3y)]))

period_years <- tibble(year = seq(from = 1997, to = 2021, by = 3),
                       period = 1:9)
cope_ct$LabelYear <- period_years$year[match(cope_ct$Period_3y,
                                                period_years$period)]

write.csv(cope_ct, here("data", "nhl-cope-ct.csv"))

## Regression CTa x Time ----------------------------------------------------
n_cope_regression <- lm(ct_anomaly_days ~ Period_3y,
                        data = cope_ct[cope_ct$Assemblage == "Northern",])
print(summary(n_cope_regression))

s_cope_regression <- lm(ct_anomaly_days ~ Period_3y,
                        data = cope_ct[cope_ct$Assemblage == "Southern",])
print(summary(s_cope_regression))


## CT/CTa when month 1 = October or July --------------------------------------
# shift n/s copepods to start in Oct
ns_cope2 <- cope %>%
  dplyr::select(Year, Month, SouthernLog10Biomass..mg.m3., NorthernLog10Biomass..mg.m3.) %>%
  dplyr::rename(Southern = SouthernLog10Biomass..mg.m3.,
                Northern = NorthernLog10Biomass..mg.m3.) %>%
  melt(id.vars = c("Year", "Month"),
       value.name = "Log10Biomass",
       variable.name = "Assemblage")

ns_cope_shift_oct <- ns_cope2 %>%
  mutate(New_Year =
           ifelse(Month >= 10, Year, Year - 1),
         New_Month =
           ifelse(Month >= 10, Month - 9,
                  ifelse(Month < 10, Month + 3, NA)),
         Year_Start = "October") # assign new month and year values

ns_cope_shift_oct2 <- ns_cope_shift_oct %>%
  mutate(Period_3y = period_years$period[match(ns_cope_shift_oct$New_Year,
                                             period_years$year)]) # match shifted years to periods

ns_cope_shift_oct_3y <- ns_cope_shift_oct2 %>%
  group_by(Period_3y, Assemblage, Month) %>%
  summarise(avg_Cope = mean(Log10Biomass)) # average over periods

ns_cope_shift_oct_3y <- ns_cope_shift_oct_3y %>%
  group_by(Period_3y, Assemblage) %>%
  summarise(ct = (sum(Month * avg_Cope) / sum(avg_Cope))) %>%
  ungroup() %>%
  drop_na() %>%
  group_by(Assemblage) %>%
  mutate(mean_ct_months = mean(ct),
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44) # get ct and cta

write.csv(ns_cope_shift_oct_3y, here("data", "nhl-cope-oct-ct.csv"), row.names = F)

# shift n/s cope to start in July
ns_cope_shift_july <- ns_cope2 %>%
  mutate(New_Year =
           ifelse(Month >= 7, Year, Year - 1),
         New_Month =
           ifelse(Month >= 7, Month - 6,
                  ifelse(Month < 7, Month + 6, NA)),
         Year_Start = "July")

ns_cope_shift_july2 <- ns_cope_shift_july %>%
  mutate(Period_3y = period_years$period[match(ns_cope_shift_july$New_Year,
                                             period_years$year)]) # match shifted years to periods

ns_cope_shift_july_3y <- ns_cope_shift_july2 %>%
  group_by(Period_3y, Assemblage, Month) %>%
  summarise(avg_Cope = mean(Log10Biomass)) # average over periods

ns_cope_shift_july_3y <- ns_cope_shift_july_3y %>%
  group_by(Period_3y, Assemblage) %>%
  summarise(ct = (sum(Month * avg_Cope) / sum(avg_Cope))) %>%
  ungroup() %>%
  drop_na() %>%
  group_by(Assemblage) %>%
  mutate(mean_ct_months = mean(ct),
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44) # get ct and cta

write.csv(ns_cope_shift_july_3y, here("data", "nhl-cope-july-ct.csv"), row.names = F)


# BCUI Analysis -----------------------------------------------------------
bcui <- read.csv(here("data", "nhl-bcui.csv"))
bcui_3y <- read.csv(here("data", "nhl-bcui-3y.csv"))

## CT/CTa when month 1 = January ----------------------------------------------
onlypos_bcui_ct <- bcui_3y %>%
  group_by(Period_3y) %>%
  summarise(ct = (sum(MONTH * Avg_BCUI) / sum(Avg_BCUI))) %>%
  ungroup() %>%
  mutate(ct_anomaly_months = ct - mean(ct),
         ct_anomaly_days = ct_anomaly_months * 30.44) # get ct and cta

# ct_function <- function(bcui_3y, i){
#   x <- bcui_3y[i,] %>%
#     group_by(Period_3y) %>%
#     reframe(ct = (sum(MONTH * avg_BCUI) / sum(avg_BCUI)))
#   return(x$ct)
# }
#
# set.seed(125)
# onlypos_bcui_ct_boot = by(bcui_3y, bcui_3y$Period_3y,
#                   function(i)boot(i, ct_function, R = 1000))
# onlypos_bcui_ct_boot

# from the print output of the bootstrap, we get the standard errors:
onlypos_bcui_ct_se <- c(0.2022297, 0.1805425, 0.1925684,
                        0.2557991, 0.180777, 0.2128862,
                        0.1810176, 0.2340665, 0.2151879)

onlypos_bcui_ct <- onlypos_bcui_ct %>%
  mutate(se_ct_months = onlypos_bcui_ct_se,
         lower_ct = ct - 1.96 * se_ct_months,
         upper_ct = ct + 1.96 * se_ct_months)

# cta_function <- function(bcui_3y, i){
#   x <- bcui_3y[i,] %>%
#     group_by(Period_3y) %>%
#     reframe(ct = (sum(MONTH * avg_BCUI) / sum(avg_BCUI)))
#   y <- x %>%
#     mutate(mean_ct_months = mean(ct),
#            ct_anomaly_months = ct - mean_ct_months,
#            ct_anomaly_days = ct_anomaly_months * 30.44)
#   return(y$ct_anomaly_days)
# }
#
# set.seed(126)
# onlypos_bcui_cta_boot <- boot(onlypos_bcui_3y, cta_function, R = 1000)
# onlypos_bcui_cta_boot

# from the print output, we get the standard errors of CTa:
onlypos_bcui_cta_se <- c(6.260224, 5.613621, 5.888801,
                         7.870936, 5.605435, 6.622123,
                         5.532436, 6.999252, 6.490995)

# add the CTa SEs to the dataframe:
onlypos_bcui_ct$se_cta_days <- onlypos_bcui_cta_se

period_years <- tibble(year = seq(from = 1997, to = 2021, by = 3),
                       period = 1:9)
onlypos_bcui_ct$Year <- period_years$year[match(onlypos_bcui_ct$Period_3y,
                                                period_years$period)]

write.csv(onlypos_bcui_ct, here("data", "nhl-bcui-ct.csv"))

## Regression CTa x Time ---------------------------------------------------
onlypos_bcui_regression <- lm(ct_anomaly_days ~ Period_3y,
                              data = onlypos_bcui_ct)
print(summary(onlypos_bcui_regression))

## CT/CTa when month 1 = October or July ------------------------------------
# shift BCUI to start in October
onlypos_bcui_shift_oct <- bcui %>%
  mutate(New_Year =
           ifelse(MONTH >= 10, YEAR, YEAR - 1),
         New_Month =
           ifelse(MONTH >= 10, MONTH - 9,
                  ifelse(MONTH < 10, MONTH + 3, NA)),
         Year_Start = "October") # assign new month and year values

onlypos_bcui_shift_oct_3y <- onlypos_bcui_shift_oct %>%
  mutate(Period_3y = period_years$period[match(onlypos_bcui_shift_oct$New_Year,
                                             period_years$year)]) # attach the periods to the years

onlypos_bcui_shift_oct_3y <- onlypos_bcui_shift_oct_3y %>%
  group_by(Period_3y, MONTH) %>%
  summarise(Avg_BCUI = mean(pos_bcui)) # average x month x period

onlypos_bcui_shift_oct_3y <- onlypos_bcui_shift_oct_3y %>%
  group_by(Period_3y) %>%
  summarise(ct = (sum(MONTH * Avg_BCUI) / sum(Avg_BCUI))) %>%
  ungroup() %>%
  mutate(ct_anomaly_months = ct - mean(ct),
         ct_anomaly_days = ct_anomaly_months * 30.44) # get ct and cta

write.csv(onlypos_bcui_shift_oct_3y, here("data", "nhl-bcui-oct-ct.csv"),
          row.names = F)

# shift BCUI to start in July
onlypos_bcui_shift_july <- bcui %>%
  mutate(New_Year =
           ifelse(MONTH >= 7, YEAR, YEAR - 1),
         New_Month =
           ifelse(MONTH >= 7, MONTH - 6,
                  ifelse(MONTH < 7, MONTH + 6, NA)),
         Year_Start = "July") #assign new month and year values

onlypos_bcui_shift_july_3y <- onlypos_bcui_shift_july %>%
  mutate(Period_3y = period_years$period[match(onlypos_bcui_shift_july$New_Year,
                                             period_years$year)]) # attach the periods to the years

onlypos_bcui_shift_july_3y <- onlypos_bcui_shift_july_3y %>%
  group_by(Period_3y, MONTH) %>%
  summarise(Avg_BCUI = mean(pos_bcui)) # average x month x period

onlypos_bcui_shift_july_3y <- onlypos_bcui_shift_july_3y %>%
  group_by(Period_3y) %>%
  summarise(ct = (sum(MONTH * Avg_BCUI) / sum(Avg_BCUI))) %>%
  ungroup() %>%
  mutate(ct_anomaly_months = ct - mean(ct),
         ct_anomaly_days = ct_anomaly_months * 30.44) # get ct and cta

write.csv(onlypos_bcui_shift_july_3y, here("data", "nhl-bcui-july-ct.csv"),
          row.names = F)