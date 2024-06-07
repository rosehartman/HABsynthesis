#get teh dayflow data

library(tidyverse)
library(lubridate)
library(smonitr)
#I need to read in the Dayflow data from the CNRA portal
# https://data.cnra.ca.gov/dataset/dayflow
#Still needs a little fiddling, but much better!.

Dayflow = get_odp_data(pkg_id = "dayflow", fnames = "Dayflow Results")



DF1997_2020 =  Dayflow$`Dayflow Results 1997 - 2023` %>%
  mutate( Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  select(Date, OUT, EXPORTS, SJR, GCD, SAC, CVP, SWP, X2, TOT)

#now I can put them all together!
DF = DF1997_2020
save(DF, file = "Dayflow1997_2023.RData")

#what's inflow in the summer?

summer = filter(DF, month(Date) %in% c(6,7,8,9))

mean(summer$TOT)*.02832


#Mean annual flow

meanflow = mutate(DF, Year = year(Date)) %>%
  group_by(Year) %>%
  summarise(across(c(OUT:TOT), mean, .names = "{.col}_mean"))

#mean spring flow
meanspringflow = mutate(DF, Year = year(Date), Month = month(Date)) %>%
  filter(Month %in% c(2,3,4,5)) %>%
  group_by(Year) %>%
  summarise(across(c(OUT:TOT), mean, .names = "{.col}_spr_mean"))

write.csv(meanflow, "outputs/meanflow")

write.csv(meanspringflow, "outputs/meanspringflow")
