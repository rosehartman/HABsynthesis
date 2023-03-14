# Maps maps maps

library(tidyverse)
library(lubridate)
library(sf)
library(ggmap)
library(readxl)
library(deltamapr)

load("stations.RData")

dcd_nodes = st_read("dsm2_8_2_1_shapefiles/dsm2_v8_2_1_historical_and_dcd_nodes.shp")

smcd_nodes = st_read("dsm2_8_2_1_shapefiles/dsm2_v8_2_1_historical_and_smcd_nodes.shp")

centerline = st_read("dsm2_8_2_1_shapefiles/dsm2_v8_2_1_historical_centerline_chan_norest.shp")

ggplot() + geom_sf(data = WW_Delta) + geom_sf(data = dcd_nodes, color = "blue")+
  geom_sf(data = smcd_nodes, color = "green") + geom_sf(data = centerline)

#import group designations
DSM2_chan <- read_excel("data/DSM2_chan.xlsx")

centerline = left_join(centerline, DSM2_chan, by = c("id" = "channel number"))


ggplot() + geom_sf(data = WW_Delta) + 
  #geom_sf(data = dcd_nodes, color = "blue")+
  #geom_sf(data = smcd_nodes, color = "green") + 
  geom_sf(data = centerline, aes(color = tag), size = 2)+
   scale_color_viridis_d(option = "turbo", na.value = "grey")+
  geom_sf(data = stashap, aes(size = Microcystis))+
  coord_sf(xlim = c(-122.2, -121.2), ylim = c(37.6, 38.6))
 

#just plot the Microcystis average
ggplot() + geom_sf(data = WW_Delta) + 
  scale_color_viridis_c(option = "A", na.value = "grey")+
  geom_sf(data = stashap, aes(color = Microcystis), size = 5)+
  coord_sf(xlim = c(-122.2, -121.2), ylim = c(37.6, 38.6))


#add the regions
regions = st_read("data/shpExport/shpExport.shp")
regions = st_transform(regions, crs = st_crs(centerline)) %>%
  mutate(Region = c("WSuisun", "ESuisun", "Montezuma",  "LowerSac", "Lower SJ", "Mid SJ", 
                    "Stockton", "Franks", "Mildred", "Clifton", "SDelta", "East", "UpperSac", "Cache", "SDWSC"))

ggplot()+ geom_sf(data = WW_Delta)+ geom_sf(data = regions, aes(fill = Region), alpha = 0.3)+
 # geom_sf_label(data = regions, aes(label = Region))+
  geom_sf(data = stashap, aes(color = Microcystis))+
  scale_color_viridis_c(option = "B")+
  coord_sf(xlim = c(-122.2, -121.2), ylim = c(37.6, 38.6))+
  theme_bw()

Center2 = st_join(regions, centerline) 
stashap = st_transform(stashap, crs = st_crs(centerline))

Stations = st_join(stashap, regions)

st_write(regions, "HABregions15.shp")
write.csv(Center2, "DSM2Channels.csv", row.names = FALSE)
write.csv(Stations, "StationswRegions.csv")
