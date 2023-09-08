#Map for X2 paper

library(sf)
library(tidyverse)
library(DroughtData)
library(readxl)
library(ggspatial)

library(cowplot)
library(grid)

library(ggsn)
library(deltamapr)

#add points for Sacramento, Stockton, Martinez, and the CVP and SWP pumps

Points =read_excel("data/points.xlsx")
Points = st_as_sf(Points, coords = c("Longitude", "Latitude"), crs = 4326)
load("data/stations.RData")
stations = mutate(stations, Years = factor(years, levels = c(2, 4, 6)))
stations$Latitude= st_coordinates(stations)[,2]
stations$Longitude= st_coordinates(stations)[,1]
#X2 locations
P_X2_x = dplyr::filter(P_X2, RKI/10 == round(RKI/10), RKI <101, RKI >30)

#segments for X2  locations
P_X2_x$Latitude= st_coordinates(P_X2_x)[,2]
P_X2_x$Longitude= st_coordinates(P_X2_x)[,1]
P_X2_x = mutate(P_X2_x, yend = Latitude+.12, Y = Latitude -.12)

X2map = ggplot()+
  geom_sf(data = WW_Delta, fill = "lightblue", color = "lightgrey", alpha = .5)+
  geom_sf(data = filter(WW_Delta, HNAME == "SACRAMENTO RIVER"), fill = "seagreen", color = "seagreen", alpha = .5)+
  geom_sf(data = filter(WW_Delta, HNAME == "SAN JOAQUIN RIVER"), fill = "cadetblue", color = "cadetblue", alpha = .5)+
  
  geom_sf(data = filter(Points, Type != "flow", Type != "Island"), size = 3, color = "grey50")+
      geom_sf_text(data =  filter(Points, Type == "City"), aes(label = Label), 
                   nudge_x = 0.02, nudge_y = -0.02, fontface = "bold", size = 3, color = "grey50")+
  geom_sf_text(data =  filter(Points, Type == "POI"), aes(label = Label), 
               nudge_x = 0.01, nudge_y = -0.01, size = 3,  color = "grey50")+
  geom_spatial_segment(
    mapping = aes(x = Longitude, xend = Longitude,
                 y = Y, yend = yend),
    data =P_X2_x,
    crs = 4326)+
  geom_point(data = stations2, aes(x = Longitude, y = Latitude, fill =as.factor(years)),
             size =2, shape = 21)+
  scale_fill_brewer(palette = "YlOrRd", name = "Sampling Stations\nYears Sampled")+
  geom_sf_text(data = P_X2_x, aes(label = RKI),fontface = "bold" )+
  theme_bw()+
  scalebar( y.min = 37.8, y.max = 38.3, x.min = -122.2, x.max = -121.2, 
            transform = TRUE, dist = 10, st.size = 4, st.dist = 0.03,
            dist_unit = "km", model = "WGS84", location = "bottomleft") +
  north(y.min = 37.8, y.max = 38.3, x.min = -122.2, x.max = -121.2,  symbol = 2) +
  annotate("text", x = -122.0, y = 38.2, label = "X2 Distances", size =5)+
  theme_bw()+ylab("")+xlab("")+
  
  annotate("segment", x = -121.6, xend = -121.65, y = 38.3, yend = 38.2,
            linewidth = 2, arrow = arrow()) +
  annotate("text", x = -121.61, y = 38.3, label = "Sacramento\nRiver Inflow", size =3, hjust =1) +
  annotate("segment", x = -121.35, xend = -121.40, y = 37.85, yend = 37.9,
           linewidth = 1, arrow = arrow(length = unit(.15, "inches"))) +
  annotate("text", x = -121.35, y = 37.82, label = "San Joaquin\nRiver Inflow", size =3) +
  annotate("segment", x = -122.0, xend = -122.15, y = 37.9, yend = 37.9,
           linewidth = 2, arrow = arrow()) +
  annotate("text", x = -122, y = 37.88, label = "Delta Outflow", size =4) +
  annotate("segment", x = -121.65, xend = -121.68, y = 37.80, yend = 37.75,
           linewidth = 1, arrow = arrow(length = unit(.15, "inches"))) +
  annotate("text", x = -121.6, y = 37.83, label = "Project\nExports", size =3, hjust =1) +
 coord_sf(xlim = c(-122.2, -121.25), ylim = c(37.75, 38.3))


X2map
ggsave("plots/x2map.tiff", device = "tiff", width = 7, height = 8)
ggsave("plots/x2map.pdf", device = "pdf", width = 7, height = 7)

####################################################################
