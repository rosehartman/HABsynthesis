#get teh dayflow data

library(tidyverse)
library(lubridate)
library(smonitr)
#I need to read in the Dayflow data from the CNRA portal
# https://data.cnra.ca.gov/dataset/dayflow
#Still needs a little fiddling, but much better!.

Dayflow = get_odp_data(pkg_id = "dayflow", fnames = "Dayflow Results")



DF1997_2020 =  Dayflow$`Dayflow Results 1997 - 2020` %>%
  mutate( Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  select(Date, OUT, EXPORTS, SJR, GCD, SAC, CVP, SWP, X2)

DF2021 =  Dayflow$`Dayflow Results 2021` %>%
  select(Date, OUT, EXPORTS, SJR, GCD, SAC, CVP, SWP, X2)

DF2022 =  Dayflow$`Dayflow Results 2022` %>%
  select(Date, OUT, EXPORTS, SJR, GCD, SAC, CVP, SWP, X2)

#now I can put them all together!
DF = bind_rows(DF1997_2020, DF2021, DF2022)
save(DF, file = "Dayflow1997_2021.RData")
