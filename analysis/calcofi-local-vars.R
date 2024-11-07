# Script Header -----------------------------------------------------------
# Topic: CalCOFI Local Oceanic Variables Central Tendency Analysis
# Author: Kathryn Chen
# Contact: ksc005@ucsd.edu
# Date: November 2024

# Objective: This script takes data files of CalCOFI bottle data, zooplankton displacement volume, and the Bakun coastal upwelling index at 45˚N, 125˚W. For each variable, it calculates central tendency (CT) and central tendency anomaly (CTa) with associated standard errors and regresses CTa against time. CT and CTa are calculated two times, where month 1 is set as January or July (to match ichthyoplankton species whose CT analyses were shifted to begin in these months), and the script outputs data files for each (2 starting months x 3 variables = 6 output data files).

# Set Up ------------------------------------------------------------------
# location declaration
here::i_am("analysis/calcofi-local-vars.R")

# load libraries
library(plyr)
library(tidyverse)
library(reshape2)
library(here)

# Temperature Analysis -------------------------------------------------------
# load CalCOFI temperature dataset
bottle_data <- read.csv(here("data", "calcofi-bottle.csv"),
                        check.names = F)

# extract the year, month, and line from the ID variable; filter for depth <= 100m since we're focusing on epipelagic T, also filter for the 66 standard stations and years >= 1951
temp <- bottle_data %>%
  dplyr::select(Depth_ID, T_degC, R_Depth) %>%
  filter(R_Depth <= 100) %>%
  mutate(ID = substr(Depth_ID, 1, 7),
         date = gsub("-","", ID),
         year = as.numeric(substr(date, 1, 4)),
         month = as.numeric(substr(date, 5, 6)),
         line = as.numeric(substr(Depth_ID, 23, 25))) %>%
  filter(line >= 76.7 & line <= 93.3) %>%
  dplyr::select(year, month, line, T_degC, R_Depth) %>%
  filter(year >= 1951) %>%
  na.omit()

# we want to average decadally to compare with ichthyoplankton decadal central tendencies. add a column for decade:
temp$decade <-
  as.factor(
    ifelse(temp$year <=1959, 1955,
           ifelse(temp$year >=1960 & temp$year <=1969, 1965,
                  ifelse(temp$year >=1970 & temp$year <=1979, 1975,
                         ifelse(temp$year >=1980 & temp$year <=1989, 1985,
                                ifelse(temp$year >=1990 & temp$year <=1999, 1995,
                                       ifelse(temp$year >=2000 & temp$year <=2009, 2005, 2015))))))
  )

## CT/CTa when month 1 = January ----------------------------------------------
# next, we will start calculating central tendency (CT) and central tendency anomalies (CTa) for SST. we start with taking decadal averages:
temp_averages <- temp %>%
  group_by(decade, month) %>%
  summarise(avg_t = mean(T_degC)) %>%
  ungroup()

# from Asch (2015): "Unlike upwelling and larval fish abundance, measurements of SST and zooplankton volume never approached zero during any month of the year. This pattern weakened the weighting of these variables by month when calculating the CT, resulting in a CT skewed toward the middle of the year. To increase the influence of monthly weights, the minimum monthly mean value (e.g., 13 °C for SST and 38 cm3/1,000 m3 of seawater strained for zooplankton volume) was subtracted from the decadally averaged time series before calculating the CT"

min(temp_averages$avg_t) #our min value is 12.53957˚C
temp_averages[which(temp_averages$avg_t == min(temp_averages$avg_t)),] #this minimum value is in April of decade 1955

# subtract the rounded minimum from all values:
temp_averages2 <- temp_averages %>%
  mutate(adjusted_t = avg_t - 12)

# now we can calculate central tendency:
temp_ct <- temp_averages2 %>%
  group_by(decade) %>%
  summarise(ct = (sum(month * adjusted_t) / sum(adjusted_t))) %>%
  ungroup() %>%
  mutate(mean_ct_months = mean(ct),
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44)

# # calculate standard errors and 95% confidence intervals for central tendency:
# ct_function <- function(temp, i){
#   x <- temp[i,] %>%
#     group_by(decade, month) %>%
#     reframe(adjusted_t = mean(T_degC) - 12) %>%
#     group_by(decade) %>%
#     reframe(ct = (sum(month * adjusted_t) / sum(adjusted_t)))
#   return(x$ct)
# }
#
# set.seed(123)
# temp_ct_boot = by(temp, temp$decade,
#                   function(i)boot(i, ct_function, R = 1000))

# from the print output of the bootstrap, we get the standard errors. save these as a new vector, then we can calculate the 95% CI as ± 1.96*SE:
temp_ct_se <- c(0.03414565, 0.06182696, 0.05594333, 0.02504391,
               0.0252116, 0.02647175, 0.02877941)

temp_ct <- temp_ct %>%
  mutate(se_ct_months = temp_ct_se,
         lower_ct = ct - 1.96 * se_ct_months,
         upper_ct = ct + 1.96 * se_ct_months)

# # we also want to calculate SEs for CTa:
# cta_function <- function(temp, i){
#   x <- temp[i,] %>%
#     group_by(decade, month) %>%
#     reframe(avg_t = mean(T_degC) - 12) %>%
#     group_by(decade) %>%
#     reframe(ct = (sum(month * avg_t) / sum(avg_t)))
#   y <- x %>%
#     mutate(mean_ct_months = mean(ct),
#            ct_anomaly_months = ct - mean_ct_months,
#            ct_anomaly_days = ct_anomaly_months * 30.44)
#   return(y$ct_anomaly_days)
# }
#
# set.seed(124)
# temp_cta_boot <- boot(temp, cta_function, R = 1000)

# from the print output, we get the standard errors of CTa:
temp_cta_se <- c(0.9752994, 1.6038802, 1.5713348, 0.7639873, 0.7941501, 0.8492110,
                0.8808168)

# add the CTa SEs to the dataframe:
temp_ct$se_cta_days <- temp_cta_se

## pull out the dataset needed for this figure
write.csv(temp_ct, here("data", "calcofi-temp-ct.csv"))

## Regression CTa x Time ----------------------------------------------------
temp_regression <- lm(ct_anomaly_months ~ as.numeric(as.character(decade)),
                     data = temp_ct)

print(summary(temp_regression))
# Asch (2015) had: Y = 33.7200 − 0.0170X, F = 11.0, df = 4, P < 0.05
# we get: y = 24.410776 - 0.012298x, F = 6.849, df = 5, p = 0.04726, this is comparable

## CT/CTa when month 1 = July -----------------------------------------------
## shift SST to start in July
temp_shift_july <- temp %>%
  mutate(New_Year =
           ifelse(month >= 7, year, year - 1),
         New_Month =
           ifelse(month >= 7, month - 6,
                  ifelse(month < 7, month + 6, NA)),
         Year_Start = "July")

temp_shift_july2 <- temp_shift_july %>%
  mutate(New_Decade =
           as.factor(
             ifelse(New_Year <=1959, 1955,
                    ifelse(New_Year >=1960 & New_Year <=1969, 1965,
                           ifelse(New_Year >=1970 & New_Year <=1979, 1975,
                                  ifelse(New_Year >=1980 & New_Year <=1989, 1985,
                                         ifelse(New_Year >=1990 & New_Year <=1999, 1995,
                                                ifelse(New_Year >=2000 & New_Year <=2009, 2005, 2015))))))
           )
  ) # match shifted years to decades

temp_shift_july3 <- temp_shift_july2 %>%
  group_by(New_Decade, month) %>%
  summarise(avg_t = mean(T_degC)) #average over periods

min(temp_shift_july3$avg_t) #our min value is 12.55345˚C

temp_shift_july3 <- temp_shift_july3 %>%
  mutate(adjusted_t = avg_t - 12) #subtract rounded minimum value

temp_shift_july3 <- temp_shift_july3 %>%
  group_by(New_Decade) %>%
  summarise(ct = (sum(month * adjusted_t) / sum(adjusted_t))) %>%
  ungroup() %>%
  mutate(mean_ct_months = mean(ct),
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44) #get ct and cta

write.csv(temp_shift_july3, here("data", "calcofi-temp-july-ct.csv"), row.names = F)

# Zooplankton Analysis -------------------------------------------------------
# load CalCOFI zooplankton dataset
zdv_data <- read.csv(here("data", "calcofi-zooplankton.csv"))

# we are interested in the column "Sml_PVolC3" which is the small (organisms <5mL) zooplankton displacement volume (cm^3) per 1000 m^3 water strained.
# we also need the year and month, which is contained in variable "Cruise", and we must filter for the 66 standard stations.
zdv_data2 <- zdv_data %>%
  dplyr::select("Cruise", "Sml_PVolC3", "St_Line") %>%
  filter(St_Line >= 76.7 & St_Line <= 93.3) %>%
  mutate(year = as.numeric(substr(Cruise, 1, 4)),
         month = as.numeric(substr(Cruise, 5, 6))) %>%
  filter(year >= 1951) %>%
  dplyr::select("year", "month", "St_Line", "Sml_PVolC3") %>%
  na.omit()

# now we can add a column for denoting decade, since we will need to take decadal averages
zdv_data2$decade <- as.factor(
  ifelse(zdv_data2$year <=1959, 1955,
         ifelse(zdv_data2$year >=1960 & zdv_data2$year <=1969, 1965,
                ifelse(zdv_data2$year >=1970 & zdv_data2$year <=1979, 1975,
                       ifelse(zdv_data2$year >=1980 & zdv_data2$year <=1989, 1985,
                              ifelse(zdv_data2$year >=1990 & zdv_data2$year <=1999, 1995,
                                     ifelse(zdv_data2$year >=2000 & zdv_data2$year <=2009, 2005, 2015))))))
)

# from Asch (2015): "zooplankton volume from 1969 to 1977 was multiplied by a correction factor of 1.366 to account for changes in net type and tow depth"
zdv_data3 <- zdv_data2 %>%
  mutate(Sml_PVolC3 = ifelse(year >= 1969 & year <= 1977,
                             Sml_PVolC3 * 1.366,
                             Sml_PVolC3))


## CT/CTa when month 1 = January ----------------------------------------------
# next, in following with the methods of Asch (2015) we will subtract the rounded minimum monthly value of ZDV from all the data points. the reason is that ZDV does not ever approach 0, which would result in CT values skewed towards the middle of the year. to correct for this:
zdv_averages <- zdv_data3 %>%
  group_by(decade, month) %>%
  summarise(avg_zdv = mean(Sml_PVolC3)) %>%
  ungroup()

min(zdv_averages$avg_zdv) #our min value is 28.10484
zdv_averages[which(zdv_averages$avg_zdv == min(zdv_averages$avg_zdv)),] #this minimum value is in September of decade 1995

zdv_averages2 <- zdv_averages %>%
  mutate(adjusted_zdv = avg_zdv - 28)


# now, we can calculate ZDV decadal central tendencies (CT) and central tendency anomalies (CTa)
zdv_ct <- zdv_averages2 %>%
  group_by(decade) %>%
  summarise(ct = (sum(month * adjusted_zdv) / sum(adjusted_zdv))) %>%
  ungroup() %>%
  mutate(mean_ct_months = mean(ct),
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44)

## we want to calculate standard error and 95% confidence intervals for ZDV CT. this section includes a bootstrap function. it is commented out to reduce sourcing time.
# ct_function <- function(zdv_data3, i){
#   x <- zdv_data3[i,] %>%
#     group_by(decade, month) %>%
#     reframe(adjusted_zdv = mean(Sml_PVolC3) - 28) %>%
#     group_by(decade) %>%
#     reframe(ct = (sum(month * adjusted_zdv) / sum(adjusted_zdv)))
#   return(x$ct)
# }
#
# set.seed(126)
# zdv_ct_boot = by(zdv_data3, zdv_data3$decade,
#                   function(i)boot(i, ct_function, R = 1000))

# from the print output of the bootstrap object, we get the standard errors, which we can save as a new vector. then we can calculate the 95% CI using ± 1.96*SE rule:
zdv_ct_se <- c(0.0812132, 0.09585983, 0.1488053, 0.166341, 0.2223582, 0.1404659,
               0.2240938)

zdv_ct <- zdv_ct %>%
  mutate(se_ct_months = zdv_ct_se,
         lower_ct = ct - 1.96 * se_ct_months,
         upper_ct = ct + 1.96 * se_ct_months)


## we also want to calculate SEs for CTa. again, the bootstrap here is commented out to reduce sourcing time.
# cta_function <- function(zdv_data3, i){
#   x <- zdv_data3[i,] %>%
#     group_by(decade, month) %>%
#     reframe(adjusted_zdv = mean(Sml_PVolC3) - 28) %>%
#     group_by(decade) %>%
#     reframe(ct = (sum(month * adjusted_zdv) / sum(adjusted_zdv)))
#   y <- x %>%
#     mutate(mean_ct_months = mean(ct),
#            ct_anomaly_months = ct - mean_ct_months,
#            ct_anomaly_days = ct_anomaly_months * 30.44)
#   return(y$ct_anomaly_days)
# }
#
# set.seed(127)
# zdv_cta_boot <- boot(zdv_data3, cta_function, R = 1000)

# from the print output, we get the standard errors:
zdv_cta_se <- c(2.769301, 3.189663, 4.270011, 4.580393, 6.074944,
                4.086008, 6.131126)

zdv_ct$cta_se <- zdv_cta_se

# save these as csv
write.csv(zdv_ct, file = here("data", "calcofi-zdv-ct.csv"), row.names = F)

## Regression CTa x Time ----------------------------------------------------
# here we want to regress ZDV CTa against year. in Asch (2015), the individual linear regression is not significant.
zdv_regression <- lm(ct_anomaly_days ~ as.numeric(as.character(decade)),
                     data = zdv_ct)
print(summary(zdv_regression)) #not significant, F = 0.6324, df = 5, p = 0.4625

## CT/CTa when month 1 = July -----------------------------------------------
zdv_shift_july <- zdv_data3 %>%
  mutate(New_Year =
           ifelse(month >= 7, year, year - 1),
         New_Month =
           ifelse(month >= 7, month - 6,
                  ifelse(month < 7, month + 6, NA)),
         Year_Start = "July")

zdv_shift_july2 <- zdv_shift_july %>%
  mutate(New_Decade =
           as.factor(
             ifelse(New_Year <=1959, 1955,
                    ifelse(New_Year >=1960 & New_Year <=1969, 1965,
                           ifelse(New_Year >=1970 & New_Year <=1979, 1975,
                                  ifelse(New_Year >=1980 & New_Year <=1989, 1985,
                                         ifelse(New_Year >=1990 & New_Year <=1999, 1995,
                                                ifelse(New_Year >=2000 & New_Year <=2009, 2005, 2015))))))
           )
  ) # match shifted years to decades

zdv_shift_july_3y <- zdv_shift_july2 %>%
  group_by(New_Decade, month) %>%
  summarise(avg_zdv = mean(Sml_PVolC3)) #average over periods

min(zdv_shift_july_3y$avg_zdv) #our min value is 28.10484˚C

zdv_shift_july_3y <- zdv_shift_july_3y %>%
  mutate(adjusted_zdv = avg_zdv - 28) #subtract rounded minimum value

zdv_shift_july_3y <- zdv_shift_july_3y %>%
  group_by(New_Decade) %>%
  summarise(ct = (sum(month * adjusted_zdv) / sum(adjusted_zdv))) %>%
  ungroup() %>%
  mutate(mean_ct_months = mean(ct),
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44) #get ct and cta

write.csv(zdv_shift_july_3y, here("data", "calcofi-zdv-july-ct.csv"), row.names = F)

# Upwelling Analysis -------------------------------------------------------
# read in Bakun upwelling index data. this is the BCUI at 33˚N and 119˚W provided by NOAA's Environmental Research Division (https://oceanview.pfeg.noaa.gov/products/upwelling/intro)
bcui <- read.csv(here("data", "calcofi-bcui.csv"))

# organize the data and filter for our desired time period >= 1951
bcui2 <- bcui %>%
  melt(id.vars = c("Position", "Year"),
       variable = "Month",
       value.name = "bcui") %>%
  mutate(Month = as.numeric(Month)) %>%
  filter(Year >= 1951) %>%
  na.omit()

# attach decade
bcui2$decade <-
  as.factor(
    ifelse(bcui2$Year <=1959, 1955,
           ifelse(bcui2$Year >=1960 & bcui2$Year <=1969, 1965,
                  ifelse(bcui2$Year >=1970 & bcui2$Year <=1979, 1975,
                         ifelse(bcui2$Year >=1980 & bcui2$Year <=1989, 1985,
                                ifelse(bcui2$Year >=1990 & bcui2$Year <=1999, 1995,
                                       ifelse(bcui2$Year >=2000 & bcui2$Year <=2009, 2005, 2015))))))
  )

# starting with the year x month x BCUI dataframe, convert any negative BCUI values to 0. this is to keep same methods as the Newport Line analysis (BCUI at NHL had a lot of negative values, which we chose to coerce to 0)
onlypos_bcui <- bcui2 %>%
  mutate(bcui = ifelse(bcui < 0, 0, bcui))

## CT/CTa when month 1 = January ----------------------------------------------
# next, we can calculate the decadal central tendency (CT) and central tendency anomalies (CTa)
bcui_averages <- onlypos_bcui %>%
  group_by(decade, Month) %>%
  summarise(avg_bcui = mean(bcui)) %>%
  ungroup()

bcui_ct <- bcui_averages %>%
  group_by(decade) %>%
  summarise(ct = (sum(Month * avg_bcui) / sum(avg_bcui))) %>%
  ungroup() %>%
  mutate(mean_ct_months = mean(ct),
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44)

## we will also calculate the standard errors and 95% CI for the CTs using bootstrap analysis. this code has been commented out for faster sourcing.
# ct_function <- function(onlypos_bcui, i){
#   x <- onlypos_bcui[i,] %>%
#     group_by(decade, Month) %>%
#     reframe(avg_bcui = mean(bcui)) %>%
#     group_by(decade) %>%
#     reframe(ct = (sum(Month * avg_bcui) / sum(avg_bcui)))
#   return(x$ct)
# }
#
# set.seed(123)
# bcui_ct_boot = by(onlypos_bcui, onlypos_bcui$decade,
#               function(i)boot(i, ct_function, R = 1000))

# from the print output of the bootstrap, we get the standard errors. save these as a new vector, then we can calculate the 95% CI as ± 1.96*SE:
bcui_ct_se <- c(0.07324911, 0.06557856, 0.05307275,
                0.05358981, 0.0609192, 0.08329756,
                0.04544613)

bcui_ct <- bcui_ct %>%
  mutate(se_ct_months = bcui_ct_se,
         lower_ct = ct - 1.96 * se_ct_months,
         upper_ct = ct + 1.96 * se_ct_months)

## we also want to calculate the SEs for CTa. the bootstrap here has also been commented out for faster sourcing.
# cta_function <- function(onlypos_bcui, i){
#   x <- onlypos_bcui[i,] %>%
#     group_by(decade, Month) %>%
#     reframe(avg_bcui = mean(bcui)) %>%
#     group_by(decade) %>%
#     reframe(ct = (sum(Month * avg_bcui) / sum(avg_bcui)))
#   y <- x %>%
#     mutate(mean_ct_months = mean(ct),
#            ct_anomaly_months = ct - mean_ct_months,
#            ct_anomaly_days = ct_anomaly_months * 30.44)
#   return(y$ct_anomaly_days)
# }
#
# set.seed(124)
# bcui_cta_boot <- boot(onlypos_bcui, cta_function, R = 1000)

#from the print output, we get the standard errors:
bcui_cta_se <- c(2.090939, 1.874919, 1.597823,
                 1.555605, 1.727100, 2.270023,
                 1.384211)

bcui_ct$cta_se <- bcui_cta_se

# save calculations as a csv file
write.csv(bcui_ct, here("data", "calcofi-bcui-ct.csv"), row.names = F)

## Regression CTa x Time ----------------------------------------------------
bcui_regression <- lm(ct_anomaly_months ~ as.numeric(as.character(decade)),
                      data = bcui_ct)
print(summary(bcui_regression))
#as expected it is not significant here, we get: F = 0.1112, df = 5, p = 0.7523

## CT/CTa when month 1 = July -----------------------------------------------
bcui_shift_july <- onlypos_bcui %>%
  mutate(New_Year =
           ifelse(Month >= 7, Year, Year - 1),
         New_Month =
           ifelse(Month >= 7, Month - 6,
                  ifelse(Month < 7, Month + 6, NA)),
         Year_Start = "July")

bcui_shift_july2 <- bcui_shift_july %>%
  mutate(New_Decade =
           as.factor(
             ifelse(New_Year <=1959, 1955,
                    ifelse(New_Year >=1960 & New_Year <=1969, 1965,
                           ifelse(New_Year >=1970 & New_Year <=1979, 1975,
                                  ifelse(New_Year >=1980 & New_Year <=1989, 1985,
                                         ifelse(New_Year >=1990 & New_Year <=1999, 1995,
                                                ifelse(New_Year >=2000 & New_Year <=2009, 2005, 2015))))))
           )
  ) # match shifted years to decades

bcui_shift_july_3y <- bcui_shift_july2 %>%
  group_by(New_Decade, Month) %>%
  summarise(avg_bcui = mean(bcui)) #average over periods

bcui_shift_july_3y <- bcui_shift_july_3y %>%
  group_by(New_Decade) %>%
  summarise(ct = (sum(Month * avg_bcui) / sum(avg_bcui))) %>%
  ungroup() %>%
  mutate(mean_ct_months = mean(ct),
         sd_ct_months = sd(ct),
         ct_anomaly_months = ct - mean_ct_months,
         ct_anomaly_days = ct_anomaly_months * 30.44) #get ct and cta

write.csv(bcui_shift_july_3y, here("data", "calcofi-bcui-july-ct.csv"), row.names = F)
