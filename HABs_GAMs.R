#Look at the HAB data and experiment with models that include all regions and seasons

library(tidyverse)
library(lubridate)
library(brms)
library(mgcv)
library(sf)
library(stars)
library(deltamapr)
library(sp)
#load integrated HABs dataset

load("data/HABs.RData") 

#Create a data frame with all the stations
# Habsstas = select(HABs, Station, Latitude, Longitude) %>%
#   distinct() %>%
#   filter(!is.na(Latitude), Longitude > -122.2)
# 
# #create a convex hull around all the stations
# ch = chull(Habsstas$Longitude, Habsstas$Latitude)
# hull = Habsstas[ch,]
# coordinates(Habsstas) = ~ Longitude + Latitude
# proj4string(Habsstas) <- "+proj=longlat +datum=NAD83"
# 
# 
# delta = st_transform(WW_Delta,crs=4326)
# stas = st_as_sf(Habsstas) %>%
#   st_transform(Habsstas, crs=4326) 
# 
# 
# 
# #map of the delta with just areas close to stations
# 
# library(sfheaders)
# hullp = sf_polygon(
#    obj = hull
#    , x = "Longitude"
#    , y = "Latitude")
# st_crs(hullp) = 4326
# 
# #buffer the hull by 0.02 degrees from the nearest station
# chbuff = st_as_sf(hullp) %>%
# st_buffer(dist = 0.02)
# 
# # buffer the waterways a bit and find the intersection
# deltabuff = st_buffer(delta, dist = 300) %>%
# st_intersection(chbuff)
# # 
# deltabuff = st_make_valid(deltabuff)
# 
# ggplot() + geom_sf(data = deltabuff)#+ geom_sf(data = WW_Delta)
# 
# save(deltabuff, file = "deltabuff.RData")

load("deltabuff.RData")




ggplot(HABs, aes(x = Date, y = Microcystis, color = Source)) + geom_point()

ggplot(HABs, aes(x = Source, y = Microcystis, color = Source)) + geom_boxplot()+
  facet_wrap(~Month)
#DOP is consistently higher than the other surveys, which could be problemataic. Let's
#leave that one out.

HABs2 = filter(HABs, Source != "DOP") %>%
  mutate(Year = case_when(Source == "FMWTx" ~ 2021,
                          TRUE ~ Year),
    Source = case_when(Source == "FMWTx" ~ "FMWT",
                       Source == "DWR_NCRO" ~ "NCRO",
                       Source == "DWR_EMP" ~ "EMP",
                            TRUE ~ Source),
    HAB = case_when(Microcystis == 1 ~ 1,
                    Microcystis %in% c(2,3)~ 2,
                    Microcystis %in% c(4,5) ~ 3)) %>%
    filter(!is.na(Year), !is.na(Microcystis), Temperature >5) 

ggplot(HABs2, aes(x = Source, y = Microcystis, color = Source)) + geom_boxplot()+
  facet_wrap(~Year)

#Let's try an ordered categorical GAM
#I need to include random effects at some point, but I dont think it will work here.

GAM1 = gam(Microcystis ~ s(Month, bs = "cc", k = 10)+ s(Temperature, k = 5) + 
             te(Latitude, Longitude), data = HABs2, family = ocat(R = 5))
summary(GAM1)
plot(GAM1)
gam.check(GAM1)

#Matrix of values over which to predict
Temp = 5:35
Month = 1:12
Latlong = distinct(select(HABs2, Latitude, Longitude)) %>%
  mutate(Temperature = median(HABs2$Temperature, na.rm = T), Month = 7)
Latlong2 = distinct(select(HABs2, Latitude, Longitude)) %>%
  mutate(Temperature = median(HABs2$Temperature, na.rm = T), Month = 5)
Latlong3 = distinct(select(HABs2, Latitude, Longitude)) %>%
  mutate(Temperature = median(HABs2$Temperature, na.rm = T), Month = 9)

Newdat = bind_rows(Latlong, Latlong2, Latlong3)

preds = predict(GAM1, newdata = Newdat, type = "response")
preds = as.data.frame(preds)

toplot = bind_cols(Newdat, preds)
ggplot(toplot, aes(x = Longitude, y = Latitude, color = V1)) + geom_point()
ggplot(toplot, aes(x = Longitude, y = Latitude, color = V2)) + geom_point()
ggplot(toplot, aes(x = Longitude, y = Latitude, color = V3)) + geom_point()
ggplot(toplot, aes(x = Longitude, y = Latitude, color = V4)) + geom_point()
ggplot(toplot, aes(x = Longitude, y = Latitude, color = V5)) + geom_point()

vis.gam(GAM1, type = "response", view = c("Temperature", "Month"))
vis.gam(GAM1, type = "response", view = c("Longitude", "Latitude"), plot.type = "contour")
#Bleh. 

#now let's do it on three categories instead of 5
GAM2 = gam(HAB ~ s(Month, bs = "cc", k = 10)+ s(Temperature, k = 5) + 
             te(Latitude, Longitude, k = 10), data = HABs2, family = ocat(R = 3))
summary(GAM2)
plot(GAM2)
vis.gam(GAM2)
gam.check(GAM2)

preds2 = predict(GAM2, newdata = Newdat, type = "response")
preds2 = as.data.frame(preds2)


colors = colorRampPalette(colors = c("red", "yellow", "green", "orange"), bias = 0.5)

toplot2 = bind_cols(Newdat, preds2) %>%
  pivot_longer(cols = c(V1, V2, V3), names_to = "Microcystis", values_to = "Probability") %>%
  mutate(Microcystis = factor(Microcystis, levels = c("V1", "V2", "V3"), labels = c("Absent", "Low", "High"))) %>%
  filter(Longitude > -122.5)
ggplot(toplot2, aes(x = Longitude, y = Latitude, color = Probability)) + geom_point()+
  facet_grid(Month~Microcystis, scales = "free")+
 scale_color_viridis_c()


##########################################
#redo with rasterization of the whole Delta

#####################################################################
#predicions
rasterbuff = st_rasterize(deltabuff)
Points = as.data.frame(rasterbuff) %>%
  filter(!is.na(id)) %>%
  mutate(Location=1:nrow(.))%>%
  dplyr::select(Longitude=x, Latitude=y, Location) 

ggplot() + geom_stars(data = rasterbuff)

# Create full dataset for predictions
newdata<-expand.grid(Location=Points$Location,
                     Month=c(7,8,9), Temperature = c(20,25,30))%>% # Create all combinations of predictor variables
  left_join(Points)%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=FALSE)


#run predictions on the model
preds2 = predict(GAM2, newdata = newdata, type = "response")
preds2 = as.data.frame(preds2) 

#organize predictions and input dataset
toplot = bind_cols(newdata, preds2) %>%
  pivot_longer(cols = c(V1, V2, V3), names_to = "Microcystis", values_to = "Probability") %>%
  mutate(Microcystis = factor(Microcystis, levels = c("V1", "V2", "V3"), labels = c("Absent", "Low", "High"))) %>%
  filter(Longitude > -122.5) 


ggplot(filter(toplot, Month == 7), aes(x = Longitude, y = Latitude, color = Probability)) + geom_point()+
  facet_grid(Temperature~Microcystis, scales = "free")+
  scale_color_viridis_c() + ggtitle("July")

ggplot(filter(toplot, Month == 8), aes(x = Longitude, y = Latitude, color = Probability)) + geom_point()+
  facet_grid(Temperature~Microcystis, scales = "free")+
  scale_color_viridis_c()+ ggtitle("August")

ggplot(filter(toplot, Month == 9), aes(x = Longitude, y = Latitude, color = Probability)) + geom_point()+
  facet_grid(Temperature~Microcystis, scales = "free")+
  scale_color_viridis_c()+ ggtitle("September")

