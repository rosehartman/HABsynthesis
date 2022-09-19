#integrate Ellen's data 

library(tidyverse)
library(readxl)
library(lubridate)

load("data/HABs.RData")

#read in Prop1 data
Prop1 = read_excel("data/Prop 1 field data_8.22.22.xlsx")
str(Prop1)
unique(Prop1$Sample_ID)
stations = read_excel("data/Prop 1 field data_8.22.22.xlsx", sheet = "Stations") %>%
  rename(Station = `Station Name`)

#rename things and add station coordinates
Prop1 = rename(Prop1, Station = Sample_ID, Microcystis = `Visual Index`, Turbidity = turbidity, Salinity = salinity, Temperature = Temp,
               Date = Collection_Date) %>%
  mutate(Source = "Prop1", Year = year(Date), Month = month(Date)) %>%
  left_join(stations) %>%
  select(Source, Station, Month, Year, Microcystis, Turbidity, Temperature, Salinity, Tide, Latitude, Longitude)

HABsP = bind_rows(HABs, Prop1)

#OK, now stockton
Stock = read_excel("data/stocktonIVIdata.xlsx")
str(Stock)
#Huh. NO water quality with it. That makes it less than useful. I wonder if there are other places to pull WQ

Stock = mutate(Stock, Station = paste("Stockton", Site, sep = ""), Year = year(Date), Month = month(Date), Source = "StocktonWF") %>%
  rename(Microcystis = `Visual index`, Latitude = latitude, Longitude = longitude)

HABsPS = bind_rows(HABsP, Stock) %>%
  select(Source, Station, Month, Year, Microcystis, Turbidity, Temperature, Salinity, Tide, Latitude, Longitude, Secchi)

HABs = HABsPS
save(HABs, file = "data/HABs.RData")
