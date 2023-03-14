
#Let's try and recreate Peggy's Microcystis and X2 analysis

#first I"ll load a few useful libraries
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(MuMIn)


#now let's grab the microcystis data from online:
microcystis = read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.1076.1&entityid=86da4465dde063c89afb0dc256bfa619")

#first I'll explore the dataset a bit
str(microcystis)

#plot microcystis biovolume versus time
ggplot(microcystis, aes(x = Collection_Date, y = dwr.MIC_totbvL))+ geom_point()

#check for normality
hist(microcystis$dwr.MIC_totbvL)
#nope!

#try log-transforming
hist(log(microcystis$dwr.MIC_totbvL))
#oh, that looks quite lovely!http://127.0.0.1:16619/graphics/plot_zoom_png?width=1920&height=1137

ggplot(microcystis, aes(x = field.Water.temp, y = log(dwr.MIC_totbvL+1), 
                        color = as.factor(Survey_Year)))+
  geom_point()

summary(microcystis$dwr.MIC_totbvL)

#THe dayflow dataset is a bit difficult to deal with when you first download it,
#so I have a version I got ealier and cleaned up
load("Dayflow1997_2021.RData")
str(DF)

#bind dayflow data to microcystis data
micro = left_join(microcystis, DF, by = c("Collection_Date" = "Date"))

#plot microcystis versus X2
ggplot(micro, aes(x = X2, y = log(dwr.MIC_totbvL+1), 
                        color = as.factor(Survey_Year)))+
  geom_point()

ggplot(micro, aes(x = X2, y = log(dwr.MIC_totbvL+1) ))+
  geom_smooth(method = "lm")+
  geom_point(aes(color = as.factor(DWR_Site)))

read_excel()

#now set up a linear regression that is actually good
micro = mutate(micro, logBV = log(dwr.MIC_totbvL+1), Year = as.factor(Survey_Year)) %>%
  filter(!is.na(field.Water.temp), !is.na(X2), !is.na(OUT), !is.na(bryte.NH4.mgL), !is.na(bryte.NO3.mgL))
mglobal = lmer(logBV ~ field.Water.temp + X2 + OUT + bryte.NH4.mgL + bryte.NO3.mgL + field.Salinity+(1|Year),
             data = micro, na.action = "na.fail")

#Warning messages:Some predictor variables are on very different scales: consider rescaling

#so we might want to rescale variables. But for now, let's go through allthe models

dredge(mglobal)
#so the best model includes ammonium, nitrate, water temp, and X2

mbest = lmer(logBV ~ field.Water.temp + X2 + bryte.NH4.mgL + bryte.NO3.mgL + (1|Year),
             data = micro, na.action = "na.fail")

#check the diagnostic plots
plot(mbest)
#hm. All those zeros make it ugly. I'm not sur ehow much of a problem that is.
plot(hist(residuals(mbest)))
#residuals are nicely normal

#look at the model outputs
summary(mbest)

#so we still have a very strong effect of X2, and it works better than outflow, which is interesting

library(effects)
plot(allEffects(mbest))


#I wonder if water year index would be useful to add...

Yrs = read_csv("data/yearassignments.csv")
micro2 = left_join(micro, Yrs, by = c("Survey_Year" = "Year"))


mglobal2 = lmer(logBV ~ field.Water.temp + X2 + OUT + bryte.NH4.mgL + bryte.NO3.mgL + field.Salinity+Index + (1|Year),
               data = micro2, na.action = "na.fail")
dredge(mglobal2)

mbest2 = lmer(logBV ~ field.Water.temp +  bryte.NH4.mgL + Index + bryte.NO3.mgL + (1|Year),
             data = micro2, na.action = "na.fail")
mbest2.1 = lmer(logBV ~ field.Water.temp + X2 + bryte.NH4.mgL +  bryte.NO3.mgL + (1|Year),
              data = micro2, na.action = "na.fail")
AICc(mbest2)
AICc(mbest2.1)
#Huh, X2 has more explanatory power than water index


######################################################################
#Is X2 a good measure of residence time?


res1 = read_csv("data/EX_2020_locRT.csv") %>%
  pivot_longer(cols = -SimPeriod, names_to = "Location", values_to = "ResTime") %>%
  mutate(Date = as.Date(SimPeriod, format = "%d%b%Y"), Region = str_extract(Location, "[A-Z]+" ),
         Month = month(Date), Year = year(Date))

resave = group_by(res1, Region, Location, Month, Year) %>%
  summarize(Restime = mean(ResTime, na.rm =T)) %>%
  group_by( Region, Month, Year) %>%
  summarize(Restime = mean(Restime, na.rm =T)) %>%
  filter(Region != "SDWS")



#historic X2 from Hutton et al
library(readxl)
library(RColorBrewer)
X2s = read_excel("supplemental_data_wr.1943-5452.0000617_hutton3.xlsx", skip =1)
X2s = rename(X2s, X2 = `Sacramento River X2 Position (km from GGB)`) %>%
  mutate(Date = date(Date), Month = month(Date), Year = year(Date))

restime = full_join(resave, X2s, by = c("Year", "Month"))

mypal = c(brewer.pal(8, "Dark2"), brewer.pal(10, "Set2"))
ggplot(restime, aes(x = X2, y = Restime, color = Region))+ geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values = mypal)+
  facet_wrap(~Region)
