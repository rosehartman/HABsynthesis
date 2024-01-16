library(plyr)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(sf)
library(ggmap)
library(RColorBrewer)
library(scales)
library(readxl)


#-------------Nutrient basic graphs-------------------------
#quick check of hte nurient data
load("data/Regions.RData")
load("data/NewHABregions.RData")
hab_nutr_chla_mvi <- read_csv("data/hab_nutr_chla_mvi.csv")

#Use the data set that Dave made for us, change some names
#Remove Nitrate outliers and start with just the data from April to September
nc = hab_nutr_chla_mvi %>%
  rename(Ammonium_mgL = DissAmmonia,
         Nitrate_mgL = DissNitrateNitrite,
         Orthophosphate_mgL = DissOrthophos,
         Chla_ugL = Chlorophyll) 
nc <- nc[-c(289,354,671,795),]

NutsRL <- nc %>%
  mutate(., Month = month(Date), 
                Year = year(Date)) %>% 
  dplyr::filter(Month %in% c(4, 5, 6, 7, 8, 9, 10, 11, 12)) 

#make a color pallet
pal = c(brewer.pal(8, "Set2"), brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), brewer.pal(12, "Paired"), "black", "grey")

#Add regional assignments
NutsRL<- st_as_sf(NutsRL, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_join(Newregions) %>%
  st_drop_geometry() %>%
  mutate(Station = case_when(Source == "USGS_CAWSC" ~str_sub(Station, start = 6),
                             TRUE ~ Station))

save(NutsRL, file = "data/Nutrients.RData")

NutsRL <- NutsRL %>%
  mutate(WYType = case_when(
    Year == 2014 ~ "Critical",
    Year == 2015 ~ "Critical",
    Year == 2021 ~ "Critical",
    Year == 2016 ~ "Below Normal",
    Year == 2018 ~ "Below Normal",
    Year == 2020 ~ "Dry",
    Year == 2017 ~ "Wet",
    Year == 2019 ~ "Wet"))

NutsRL <- NutsRL %>%
  mutate(Season = case_when(
    Month == 4 ~ "Spring",
    Month == 5 ~ "Spring",
    Month == 6 ~ "Spring",
    Month == 7 ~ "Summer",
    Month == 8 ~ "Summer",
    Month == 9 ~ "Summer"))x

NutsRL_newregion <- NutsRL %>% pivot_wider(names_from = Region.y, values_from = Datetime) %>% rename(., "Franks/OMR" = "OMR/Franks") %>% pivot_longer(20:25, names_to = "Region.y", values_to = "Datetime", values_drop_na = TRUE) %>% relocate(Region.y, .after = Year) %>% relocate(Datetime, .after = Date)

#NutsRL_newregion %>% write_excel_csv(file = "HAB_nutr_raw_2014-2021.csv")

Nuts_df <- NutsRL %>% pivot_longer(cols = c(7,9,11,13), names_to = "par", values_to = "val") 
Nuts_sum <- Nuts_df %>% 
  group_by(Year, Season, Region.y, par) %>% 
  summarize(min = quantile(val, probs = .25, na.rm = TRUE), 
            max = quantile(val, probs = .75, na.rm = TRUE),
            mean = mean(val, na.rm = TRUE), 
            WYType = max(WYType)) 

Nuts_sum <- Nuts_sum %>% pivot_wider(names_from = Region.y, values_from = min) %>% rename(., "Franks/OMR" = "OMR/Franks") %>% pivot_longer(7:12, names_to = "Region.y", values_to = "min", values_drop_na = TRUE) %>% relocate(Region.y, .after = Season) %>% relocate(min, .after = par)

Nuts_sum$WYType <- factor(Nuts_sum$WYType, levels = c('Critical','Dry','Below Normal','Wet'))

#Nuts_sum %>% write_excel_csv(file = "HAB_nutr_sum_2014-2021.csv")


Nuts_sum %>% filter(par %in% "Nitrate_mgL") %>% 
  ggplot(., aes(x= Year, y = mean, fill = WYType)) + 
  geom_col(colour = "black") + 
  geom_errorbar(aes(ymin = min, ymax = max), width = .25) +
  facet_grid(Region.y~Season, scales = "free_y") + 
  labs(y = "Nitrate + Nitrite (mg/L)") + 
  scale_x_continuous(breaks = c(2014,2015,2016,2017,2018,2019,2020,2021)) + 
  theme_bw() +
  theme(legend.position = "top") +
  scale_fill_manual(values= c("firebrick", "darkorange2", "lightgoldenrod", "lightskyblue"), 
                    name= "Water Year")
#ggsave(filename = "Nitrate_mgL.tiff", device = "tiff", width = 12, height = 10)


Nuts_sum %>% filter(par %in% "Ammonium_mgL") %>% 
  ggplot(., aes(x= Year, y = mean, fill = WYType)) + 
  geom_col(colour = "black") + 
  geom_errorbar(aes(ymin = min, ymax = max), width = .25) +
  facet_grid(Region.y~Season, scales = "free_y") + 
  labs(y = "Ammonium (mg/L)") + 
  scale_x_continuous(breaks = c(2014,2015,2016,2017,2018,2019,2020,2021)) + 
  theme_bw() +
  theme(legend.position = "top") +
  scale_fill_manual(values= c("firebrick", "darkorange2", "lightgoldenrod", "lightskyblue"), 
                    name= "Water Year")
#ggsave(filename = "Ammonium_mgL.tiff", device = "tiff", width = 12, height = 10)

Nuts_sum %>% filter(par %in% "Orthophosphate_mgL") %>% 
  ggplot(., aes(x= Year, y = mean, fill = WYType)) + 
  geom_col(colour = "black") + 
  geom_errorbar(aes(ymin = min, ymax = max), width = .25) +
  facet_grid(Region.y~Season, scales = "free_y") + 
  labs(y = "Orthophosphate (mg/L)") + 
  scale_x_continuous(breaks = c(2014,2015,2016,2017,2018,2019,2020,2021)) + 
  theme_bw() +
  theme(legend.position = "top") +
  scale_fill_manual(values= c("firebrick", "darkorange2", "lightgoldenrod", "lightskyblue"), 
                    name= "Water Year")
#ggsave(filename = "Orthophosphate_mgL.tiff", device = "tiff", width = 12, height = 10)


#---------------------------------N:P ratio----------------
#N:P ratio

#Use the cleaned data set from line 21
#add regional assignments
ncr <- st_as_sf(nc, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_join(Newregions) %>%
  st_drop_geometry() %>%
  mutate(Station = case_when(Source == "USGS_CAWSC" ~str_sub(Station, start = 6),
                             TRUE ~ Station))

# Calculate molar nitrogen and phosphorus. 
ncln = mutate(ncr, NH4_umolL= Ammonium_mgL*71.43,
              NO3_umolL= Nitrate_mgL*71.43,
              PO4_umolL= Orthophosphate_mgL*32.26,
              DIN_umolL= NH4_umolL + NO3_umolL,
              NPratio = DIN_umolL / PO4_umolL, 
              Year = year(Date),
              Month = month(Date)) %>%
  dplyr::filter(Month %in% c(4, 5, 6, 7, 8, 9)) %>%
  select(1:5,15:17,22:24)

ncln <- ncln %>%
  mutate(WYType = case_when(
    Year == 2014 ~ "Critical",
    Year == 2015 ~ "Critical",
    Year == 2021 ~ "Critical",
    Year == 2016 ~ "Below Normal",
    Year == 2018 ~ "Below Normal",
    Year == 2020 ~ "Dry",
    Year == 2017 ~ "Wet",
    Year == 2019 ~ "Wet"))

ncln <- ncln %>%
  mutate(Season = case_when(
    Month == 4 ~ "Spring",
    Month == 5 ~ "Spring",
    Month == 6 ~ "Spring",
    Month == 7 ~ "Summer",
    Month == 8 ~ "Summer",
    Month == 9 ~ "Summer"))

Nuts_np <- ncln %>% pivot_longer(cols = 9, names_to = "par", values_to = "val") %>%
  group_by(Year, Season, Region.y, par) %>% 
  summarize(min = min(val, na.rm = TRUE),
            max = max(val, na.rm = TRUE),
            mean = mean(val, na.rm = TRUE), 
            minln = log1p(min),
            maxln = log1p(max),
            meanln = log1p(mean),
            WYType = max(WYType)) 

Nuts_np <- Nuts_np %>% pivot_wider(names_from = Region.y, values_from = min) %>% rename(., "Franks/OMR" = "OMR/Franks") %>% pivot_longer(10:15, names_to = "Region.y", values_to = "min", values_drop_na = TRUE) %>% relocate(Region.y, .after = Season) %>% relocate(min, .after = par)

Nuts_np$WYType <- factor(Nuts_np$WYType, levels = c('Critical','Dry','Below Normal','Wet'))

#Nuts_np %>% write_excel_csv(file = "HAB_NPratio_sum_2014-2021.csv")

#Average spring/summer NP ratios by Year and Region 
Nuts_np %>% ggplot(aes(x = Year, y = mean, fill = WYType)) +
  geom_col(colour = "black") + 
  geom_errorbar(aes(ymin = min, ymax = max), width = .25) +
  facet_grid(vars(Region.y), vars(Season), scales = "free_y") + 
  scale_x_continuous(breaks = c(2014,2015,2016,2017,2018,2019,2020,2021)) + 
  labs(y = "N:P Ratio") +
  geom_hline(yintercept = 16, linetype = 2, color = "red") + 
  theme_bw() +
  theme(legend.position = "top") +
  scale_fill_manual(values= c("firebrick", "darkorange2", "lightgoldenrod", "lightskyblue"), 
                    name= "Water Year")
#ggsave(filename = "NPratio_seasonal.tiff", device = "tiff", width = 12, height = 10)

#Average spring/summer Log NP ratios by Year and Region 
Nuts_np %>% ggplot(aes(x = Year, y = meanln, fill = WYType)) +
  geom_col(colour = "black") + 
  geom_errorbar(aes(ymin = minln, ymax = maxln), width = .25) +
  facet_grid(vars(Region.y), vars(Season), scales = "free_y") + 
  scale_x_continuous(breaks = c(2014,2015,2016,2017,2018,2019,2020,2021)) + 
  labs(y = "N:P Ratio") +
  geom_hline(yintercept = 2.833213, linetype = 2, color = "red") + 
  theme_bw() +
  theme(legend.position = "top") +
  scale_y_continuous(labels = function(x)(round(exp(x)-1)), breaks = breaks_pretty()) +
  scale_fill_manual(values= c("firebrick", "darkorange2", "lightgoldenrod", "lightskyblue"), 
                    name= "Water Year")
#ggsave(filename = "NPratioln_seasonal.tiff", device = "tiff", width = 12, height = 10)

Nuts_sum_NP <- full_join(Nuts_sum,Nuts_np)
#Nuts_sum_NP %>% write_excel_csv(file = "HAB_nutrients_2014-2021.csv")

