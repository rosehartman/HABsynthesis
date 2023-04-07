#integrate Ellen's data 

library(tidyverse)
library(readxl)
library(lubridate)
# 
# load("data/HABsw2022.RData")
# HABs = HABs2022
# #read in Prop1 data
#  Prop1 = read_excel("data/Prop 1 field data_8.22.22.xlsx")
#  str(Prop1)
#  unique(Prop1$Sample_ID)
# stations = read_excel("data/Prop 1 field data_8.22.22.xlsx", sheet = "Stations") %>%
#   rename(Station = `Station Name`)
# 
# #rename things and add station coordinates
# Prop1 = rename(Prop1, Station = Sample_ID, Microcystis = `Visual Index`, Turbidity = turbidity, Salinity = salinity, Temperature = Temp,
#                Date = Collection_Date) %>%
#   mutate(Source = "Prop1", Year = year(Date), Month = month(Date)) %>%
#   left_join(stations) %>%
#   select(Source, Station, Month, Year, Microcystis, Turbidity, Temperature, Salinity, Tide, Latitude, Longitude)
# 
# HABsP = bind_rows(HABs, Prop1)
# 
# #OK, now stockton
# Stock = read_excel("data/stocktonIVIdata.xlsx")
# str(Stock)
# #Huh. NO water quality with it. That makes it less than useful. I wonder if there are other places to pull WQ
# # 
# # Stock = mutate(Stock, Station = paste("Stockton", Site, sep = ""), Year = year(Date), Month = month(Date), Source = "StocktonWF") %>%
# #   rename(Microcystis = `Visual index`, Latitude = latitude, Longitude = longitude)
# # 
# # HABsPS = bind_rows(HABsP, Stock) %>%
# #   select(Source, Station, Month, Year, Microcystis, Turbidity, Temperature, Salinity, Tide, Latitude, Longitude, Secchi)
# 
# #new stockton data with all the water quality
# 
# Stock = read_excel("data/2020,2021,2022 vi data.xlsx")
# Stocksites = read_excel("data/stocktonIVIdata.xlsx") %>%
#   dplyr::select(Site, latitude, longitude) %>%
#   distinct()
# 
# str(Stock)
# Stock = mutate(Stock, Station = paste("Stockton", `Site #`, sep = ""), Year = year(Date), Month = month(Date), Source = "StocktonWF") %>%
#   rename(Site = `Site #`) %>%
#   left_join(Stocksites) %>%
#      rename(Microcystis = `Visual index`, Latitude = latitude, Longitude = longitude, Salinity = `calculated salinity`)
# HABs = filter(HABs, Source != "StocktonWF")
#  HABsPS = bind_rows(HABs, Stock) %>%
#    select(Source, Station, Date, Month, Year, Microcystis, Temperature, Salinity, Latitude, Longitude, Secchi)
# 
# HABs = HABsPS
# save(HABs, file = "data/HABs.RData")
# 
# ##################################################################
# unique(HABs$Source)
# 

#######################################################################
#Strat over

library(discretewq)

Alldata = wq(Source = c("STN", "NCRO", "FMWT", "EMP"), Start_year = 2007, End_year = 2022)


Stock = read_excel("data/2020,2021,2022 vi data.xlsx")
Stocksites = read_excel("data/stocktonIVIdata.xlsx") %>%
  dplyr::select(Site, latitude, longitude) %>%
  distinct()

str(Stock)
Stock = mutate(Stock, Station = paste("Stockton", `Site #`, sep = ""), Year = year(Date), Month = month(Date), Source = "StocktonWF") %>%
  rename(Site = `Site #`) %>%
  left_join(Stocksites) %>%
  rename(Microcystis = `Visual index`, Latitude = latitude, Longitude = longitude, Salinity = `calculated salinity`,
         DissNitrateNitrite = NO3, DissAmmonia = NH3, DissOrthophos = PO4, TotPhos = TP, TKN = TN) %>%
  mutate(DissAmmonia = as.numeric(DissAmmonia))

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
  select(Source, Station, Month, Year, Microcystis, Turbidity, Temperature, Salinity, Tide, Latitude, Longitude, Date)

HABsall = bind_rows(Alldata, Stock, Prop1)

######################################################

res1 = read_csv("data/EX_2020_locRT.csv") %>%
  pivot_longer(cols = -SimPeriod, names_to = "Location", values_to = "ResTime") %>%
  mutate(Date = as.Date(SimPeriod, format = "%d%b%Y"), Region = str_extract(Location, "[A-Z]+" ))

Yrs = read_csv("data/yearassignments.csv") %>%
  select(Year, Yr_type, Index)

#In other circumstances I'd merge these by water years, but the october of a year 
#after a dry summer is probalby more like September of the same year than it is like February of the
#next year, so we'll go by calendar year instead.


res1 = mutate(res1, Month = month(Date), Year = year(Date)) %>%
  filter(Region != "SDWS") %>%
  left_join(Yrs)

resave = group_by(res1, SimPeriod, Region, Date, Yr_type, Month, Year, Index) %>%
  summarize(ResTime = mean(ResTime)) %>%
  mutate(DOY = yday(Date))


HABregions15 = st_read("HABregions15.shp") %>%
  st_transform(crs = 4326) %>%
  st_make_valid()
#summarize residence time by water year type
resave2 = group_by(resave, Region, Yr_type, Month) %>%
  summarize(ResTime = mean(ResTime))  %>%
  rename(RegionDSM = Region)

#switch the names so they  match the other region names
reglook = read_csv("regionlookup.csv")

resave2 = left_join(resave2, reglook)

#Make the data frame a sf object and join it to the regions
HABsReg = st_as_sf(filter(HABsall, !is.na(Latitude)), coords = c("Longitude", "Latitude"), crs = 4326, remove = F) %>%
  st_join(HABregions15) %>%
  mutate(Year = case_when(Source == "FMWTx" ~ 2021,
                          TRUE ~ Year)) %>%
  left_join(Yrs) %>%
  st_drop_geometry()

#join restime and HABs by region and water year type
HABrestime = left_join(HABsReg, resave2) %>%
 # filter(!is.na(Microcystis), Source != "DOP", !is.na(Region), Month %in% c(6,7,8,9,10)) %>%
  mutate(Mic = factor(round(Microcystis), levels = c(1,2,3,4,5), labels = c("Absent", "Low", "Med", "High", "V.High")),
         Yr_type = factor(Yr_type, levels = c("Critical", "Dry", "Below Normal", "Above Normal", "Wet"))) 

ggplot(HABrestime, aes(x = Date, y = TKN))+ geom_point()+
  facet_wrap(~Source)
ggplot(HABrestime, aes(x = Date, y = Microcystis))+ geom_point()+
  facet_wrap(~Source)

ggplot(HABrestime, aes(x = Date, y = ResTime))+ geom_point()+
  facet_grid(Region~Source)


ggplot(HABrestime, aes(x = Date, y = Salinity))+ geom_point()+
  facet_grid(Region~Source)

save(HABrestime, file = "HABrestime.RData")
write.csv(HABrestime, "AlltheData.csv", row.names = F)
