#OK, i can't find a good reference for the Sacramento having higher velocities than the san joaquin

library(dataRetrieval)
library(tidyverse)
library(lubridate)


SacFlow = readNWISdata(site = c("11455420", "11455485", "11455495"),service = "iv",
                       startDate = "2014-05-01T00:00Z", endDate = "2022-05-01T12:00Z",
                              parameterCd = c("72255", "00060", "72254"))

SJRFlow =  readNWISdata(site = c("11304810", "11337190", "11313452", "11313315"), 
                                        parameterCd = c("72255", "00060", "72254"), service = "iv",
                        startDate = "2017-05-01T00:00Z", endDate = "2022-05-01T12:00Z")

SacFlow = mutate(SacFlow, River = "Sacramento")
SJRFlow = mutate(SJRFlow, River = "San Joaquin")

allflow = bind_rows(SacFlow, SJRFlow) %>%
  mutate(Year = year(dateTime), DOY = yday (dateTime))

ggplot(allflow, aes(x = DOY, y = X_72254_00000)) +
  geom_point()+
  facet_wrap(~site_no)


ggplot(allflow, aes(x = DOY, y = X_72255_00000)) +
  geom_point()+
  facet_wrap(River~site_no)+
  ylab("Water Velocity")
