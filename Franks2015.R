#quick check of visual index in 2015

library(tidyverse)

load("data/HABrestime.RData")

hab2015 = filter(HABrestime, Year == 2015)

ggplot(hab2015, aes(x = as.factor(Month), y = Microcystis))+
  geom_boxplot()+
  facet_wrap(~Region)

hab2021 = filter(HABrestime, Year == 2021)

ggplot(hab2021, aes(x = as.factor(Month), y = Microcystis))+
  geom_boxplot()+
  facet_wrap(~Region)


#also check on Peggy's data
micro2015 = filter(micro, EMP_Site == "D19") %>%
  mutate(Year = year(Collection_Date))
ggplot(micro2015, aes(x = Collection_Date, y = log(dwr.MIC_totbvL+1)))+ geom_point()+
  facet_wrap(~Year, scales = "free_x")

ggplot(micro2015, aes(x = Collection_Date, y = log(ucd.qpcr.total.MIC+1)))+ geom_point()+
  facet_wrap(~Year, scales = "free_x")

