#try a CART analysis to see which regions have similar water quality
#But maybe just limited to summer because the coverage is better?
#Actually, maybe I should do just EMP...
#UGH I dunno.

library(tidyverse)
library(lubridate)
library(brms)
library(mgcv)
library(sf)
library(stars)
library(deltamapr)
library(sp)
library(tree)
#load integrated HABs dataset

load("data/HABs.RData") 

#filter to just summer
HABs = HABs %>%
  mutate(Source = case_when(Source == "DWR_EMP" ~"EMP",
                            TRUE ~ Source)) %>%
  filter(Month %in% c(7,8,9,10),
         Source != "DOP",
         !Station %in%  c("EZ2","EZ2-SJR","EZ6","EZ6-SJR" ))

#Average by station and month, i guess

HABave = group_by(HABs, Station) %>%
  summarize(Microcystis = mean(Microcystis, na.rm = T), Temperature = mean(Temperature, na.rm = T),
   Secchi = mean(Secchi, na.rm = T), Salinity = mean(Salinity, na.rm = T), N = n())


stas = HABs %>%
  group_by(Station, Latitude, Longitude) %>%
  summarize(N = n()) %>%
  filter(N >1) %>%
  mutate(N = NULL)

stashap = left_join(HABave, stas) %>%
  filter(!is.na(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

save(stashap, file = "stations.RData")

row.names(HABave) = HABave$Station
#tempmaxmin$Station...368 = NULL

######################################################
#Map with other regions
HABregions15 = st_read("HABregions15.shp") %>%
  st_transform(crs = 4326) %>%
  st_make_valid()
#read in shapefile of the delta
delta = WW_Delta
load("stations.RData")

ggplot()+  geom_sf(data = delta)+
  geom_sf(data = HABregions15, aes(fill = Region), alpha = 0.2)+
  geom_sf(data = Stations, aes(color = Microcystis))+
  theme_bw() +
  scale_color_viridis_c(option = "B")+
  coord_sf(xlim = c(-122.2, -121.2), ylim = c(37.7, 38.6))

Stations2 = st_join(stashap, st_transform(HABregions15, crs = st_crs(stashap)))

###########################################################
#calculate distance and cluster
tempdist2 = dist(HABave, "euclidean")
tempfit2 = hclust(tempdist2, method = "ward.D")
plot(tempfit2, main = "Clusters based WQ")

#Now cut the tree

cutday = as.data.frame(cutree(tempfit2, k = c(2,4,6,8,12)))
cutday$Station = row.names(cutday)


####################################

#read in shapefile of the delta
delta = WW_Delta

#add lat/longs for the stations
stas = HABs %>%  group_by(Station, Latitude, Longitude) %>%
  summarize(N = n()) %>%
  filter(N >1)

cutday = left_join(cutday, stas, all.x =T)
names(cutday) = c("grps2", "grps4", "grps6", "grps8", "grps12","Station", "Latitude", "Longitude", "N")

#turn it into a spatial object
stashap = filter(cutday, !is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
                   crs = 4326)
save(stashap, file ="stations.rdata")


mypal2 = c(brewer.pal(8, "Dark2"), brewer.pal(5, "Set3"))

#Make plots of the stations divided into four groups, six groups
#or twelve groups based on teh cluster analysis. 
gp4 = ggplot() +
  geom_sf(data = delta, alpha = 0.5)+
  geom_sf(data = stashap, mapping =aes(color = as.factor(grps4)))+
  theme_bw() + guides(color = "none") + 
  coord_sf(xlim = c(-122.4, -121.2), ylim = c(37.6, 38.6))+
  ggtitle("4 groups")

gp4

library(RColorBrewer)
mypal = c(brewer.pal(9, "Blues"), brewer.pal(9, "BuGn"), 
          brewer.pal(9, "Greys"), brewer.pal(9, "PuBu"), brewer.pal(9, "YlGn"))
gp6 = ggplot() +
  scale_fill_brewer(palette = "Set3")+
  geom_sf(data = delta, alpha = 0.5)+
  geom_sf(data = stashap, mapping =aes(color = as.factor(grps6)))+
  
  # geom_sf_label(aes(label = Station), stashap) +
  scale_color_brewer(palette = "Set1")+
  theme_bw() + guides(color = FALSE, fill = FALSE) + 
  coord_sf(xlim = c(-122.4, -121.2), ylim = c(37.6, 38.6))+
  ggtitle("6 groups")

gp6


gp8 = ggplot() +
  scale_fill_brewer(palette = "Set3")+
  geom_sf(data = delta, alpha = 0.5)+
  geom_sf(data = stashap, mapping =aes(color = as.factor(grps8)))+
  
  # geom_sf_label(aes(label = Station), stashap) +
  scale_color_brewer(palette = "Set1")+
  theme_bw() + guides(color = FALSE, fill = FALSE) + 
  coord_sf(xlim = c(-122.4, -121.2), ylim = c(37.6, 38.6))+
  ggtitle("8 groups")

gp8




gp6a = ggplot() +
  geom_sf(data = delta, alpha = 0.5)+
  geom_sf(data = stashap, mapping =aes(color = as.factor(grps6)))+
  theme_bw() + guides(color = FALSE, fill = FALSE) + 
  #coord_sf(xlim = c(-122.2, -121), ylim = c(37.6, 38.8))+
  ggtitle("6 groups")

gp6a


gp12 = ggplot() +
  geom_sf(data = delta, alpha = 0.5)+
  geom_sf(data = stashap, mapping =aes(color = as.factor(grps12)))+
  theme_bw() + scale_color_manual(values = mypal2)+ 
  coord_sf(xlim = c(-122.4, -121.2), ylim = c(37.6, 38.6))+
  ggtitle("12 groups")
gp12

library(gridExtra)

foo = grid.arrange(gp4, gp6, gp12, nrow = 1, ncol = 3)


st_write(stashap, "HABstationgroups.shp")

#now do it just based on microcystis
##################################################################


#Average by station and month, i guess

HABlevels = filter(HABs, !is.na(Microcystis), Month %in% c(7,8,9), Year >2015) %>%
  group_by(Year, Station) %>%
  summarize(Microcystis = mean(Microcystis, na.rm = T)) %>%
  pivot_wider(id_cols = Station, names_from = Year, values_from = Microcystis) %>%
  filter(!is.na(`2021`))


row.names(HABlevels) = HABlevels$Station
#tempmaxmin$Station...368 = NULL


#calculate distance and cluster
dist2 = dist(HABlevels, "euclidean")
fit2 = hclust(dist2, method = "ward.D")
plot(fit2, main = "Clusters based on microcystis")

#Now cut the tree

cutday = as.data.frame(cutree(fit2, k = c(2,4,6,8,12)))
cutday$Station = row.names(cutday)


####################################

#read in shapefile of the delta
delta = WW_Delta

#add lat/longs for the stations
stas = HABs %>%
  group_by(Station, Latitude, Longitude) %>%
  summarize(N = n(), Microcystis = mean(Microcystis, na.rm = T)) %>%
  filter(N >1)

cutday = left_join(cutday, stas, all.x =T)
names(cutday) = c("grps2", "grps4", "grps6", "grps8", "grps12","Station", "Latitude", "Longitude", "N")

#Add the average Microcysis score and stuff
cutday2 = left_join(cutday, HABave)


#turn it into a spatial object
stashap = filter(cutday2, !is.na(Longitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326)
save(stashap, file ="stations.rdata")


mypal2 = c(brewer.pal(8, "Dark2"), brewer.pal(5, "Set3"))

#Make plots of the stations divided into four groups, six groups
#or twelve groups based on teh cluster analysis. 
gp4 = ggplot() +
  geom_sf(data = delta, alpha = 0.5)+
  geom_sf(data = stashap, mapping =aes(color = as.factor(grps4)))+
  theme_bw() + guides(color = "none") + 
  coord_sf(xlim = c(-122.4, -121.2), ylim = c(37.6, 38.6))+
  ggtitle("4 groups")

gp4

library(RColorBrewer)
mypal = c(brewer.pal(9, "Blues"), brewer.pal(9, "BuGn"), 
          brewer.pal(9, "Greys"), brewer.pal(9, "PuBu"), brewer.pal(9, "YlGn"))
gp6 = ggplot() +
  scale_fill_brewer(palette = "Set3")+
  geom_sf(data = delta, alpha = 0.5)+
  geom_sf(data = stashap, mapping =aes(color = as.factor(grps6)))+
  
  # geom_sf_label(aes(label = Station), stashap) +
  scale_color_brewer(palette = "Set1")+
  theme_bw() + guides(color = FALSE, fill = FALSE) + 
  coord_sf(xlim = c(-122.4, -121.2), ylim = c(37.6, 38.6))+
  ggtitle("6 groups")

gp6


gp8 = ggplot() +
  scale_fill_brewer(palette = "Set3")+
  geom_sf(data = delta, alpha = 0.5)+
  geom_sf(data = stashap, mapping =aes(color = as.factor(grps8)), size = 4)+
  
  # geom_sf_label(aes(label = Station), stashap) +
  scale_color_brewer(palette = "Set1")+
  theme_bw() + guides(color = FALSE, fill = FALSE) + 
  coord_sf(xlim = c(-122.4, -121.2), ylim = c(37.6, 38.6))+
  ggtitle("8 groups")

gp8




gp6a = ggplot() +
  geom_sf(data = delta, alpha = 0.5)+
  geom_sf(data = stashap, mapping =aes(color = as.factor(grps6)))+
  theme_bw() + guides(color = FALSE, fill = FALSE) + 
  #coord_sf(xlim = c(-122.2, -121), ylim = c(37.6, 38.8))+
  ggtitle("6 groups")

gp6a


gp12 = ggplot() +
  geom_sf(data = delta, alpha = 0.5)+
  geom_sf(data = stashap, mapping =aes(color = as.factor(grps12)), size = 4)+
  theme_bw() + scale_color_manual(values = mypal2)+ 
  coord_sf(xlim = c(-122.4, -121.2), ylim = c(37.6, 38.6))+
  ggtitle("12 groups")
gp12


