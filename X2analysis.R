
#Let's try and recreate Peggy's Microcystis and X2 analysis

#first I"ll load a few useful libraries
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(effects)
library(MuMIn)
library(DHARMa)
library(car)
library(AICcmodavg)
library(readxl)
library(RColorBrewer)

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
micro = left_join(microcystis, DF, by = c("Collection_Date" = "Date")) %>%
  mutate(CVPSWP = CVP + SWP)

#plot microcystis versus X2
ggplot(micro, aes(x = X2, y = log(dwr.MIC_totbvL+1), 
                        color = as.factor(Survey_Year)))+
  geom_point()

ggplot(micro, aes(x = X2, y = log(dwr.MIC_totbvL+1) ))+
  geom_smooth(method = "lm")+
  geom_point(aes(color = as.factor(DWR_Site)))


#now set up a linear regression that is actually good
micro = mutate(micro, logBV = log(dwr.MIC_totbvL/100000+1), Year = as.factor(Survey_Year), DOY = yday(Collection_Date)) %>%
  filter(!is.na(field.Water.temp), !is.na(X2), !is.na(OUT), !is.na(bryte.NH4.mgL), !is.na(bryte.NO3.mgL))

mglobal = lmer(logBV ~ field.Water.temp + X2 + OUT + DOY+bryte.NH4.mgL +CVPSWP+ bryte.NO3.mgL + field.Salinity+(1|Year),
             data = micro, na.action = "na.fail")

#Warning messages:Some predictor variables are on very different scales: consider rescaling

micro = mutate(micro, Temp = scale(field.Water.temp), X2s = scale(X2), Outs = scale(log(OUT)),
               CVPSWPs = scale(CVPSWP),
               NH4s = scale(bryte.NH4.mgL),
               Salinitys = scale(field.Salinity), NO3s = scale(bryte.NO3.mgL),
               DWR_Site = str_remove(DWR_Site, "rep"),
               Station = paste(DWR_Site, EMP_Site),
               Station = case_when(Station %in% c("OR D28A", "OR NA") ~ "OR",
                                   TRUE ~ Station))

mglobal = lmer(logBV ~ Temp + X2s + Outs+ DOY+NH4s +CVPSWPs+ NO3s + Salinitys+ (1|Year) + (1|Station),
               data = micro, na.action = "na.fail")


vif(mglobal)
#DOY and temperature has highest vif()

mglobal2 = lmer(logBV ~ Temp + X2s + Outs+CVPSWPs+ NH4s + NO3s + Salinitys+(1|Year),
               data = micro, na.action = "na.fail")

vif(mglobal2)
#much better
#so we might want to rescale variables. But for now, let's go through allthe models

dredge(mglobal2)
#so the best model includes  nitrate, water temp, and X2, and outflow, and exports
 
#but outflow and X2 are probably tooo highly correlated
ggplot(micro, aes(x = X2s, y = Outs)) + geom_point() + geom_smooth(method = "lm")
ggplot(micro, aes(x = X2s, y = CVPSWPs)) + geom_point() + geom_smooth(method = "lm")

mbest = lmer(logBV ~ Temp+X2s+NO3s+ CVPSWPs + (1|Year),
             data = micro, na.action = "na.fail")

#check the diagnostic plots
plot(mbest)
#hm. All those zeros make it ugly. I'm not sur ehow much of a problem that is.
plot(hist(residuals(mbest)))
#residuals are nicely normal
acf(residuals(mbest))
#hm. Temporal autocorrelation? I'll include DOY and station as random effects too.
#but DOY is continuous. 

ggplot(micro, aes(x = DOY, y = logBV))+geom_point()+ geom_smooth()


ggplot(micro, aes(x = EMP_Site, y = NO3s))+ geom_point()
#there appears to be one nitrate value thowing everythign off.
#let's try just empterature and x2

mbest = lmer(logBV ~ Temp+X2s+ CVPSWPs +(1|Year)+ (1|DOY) + (1|Station),
             data = micro, na.action = "na.fail")
plot(allEffects(mbest))
plot(mbest)
mbestr = simulateResiduals(mbest)
plot(mbestr)
acf(residuals(mbest))
summary(mbest)
r.squaredGLMM(mbest)
foo = summary(mbest)
write.csv(foo$coefficients, "X2model.csv")
write.csv(foo$varcor, "X2modelrandom.csv")

#comparea a model with just X2 and temp versus just exports and temp
mx2 = lmer(logBV ~ Temp+X2s+(1|DOY)+ (1|Year)+ (1|Station),
           data = micro, na.action = "na.fail")
mexp = lmer(logBV ~ Temp+CVPSWPs+(1|Year)+ (1|DOY) + (1|Station),
           data = micro, na.action = "na.fail")

AIC(mx2)
AIC(mexp)
AIC(mbest)
BIC(mx2)
BIC(mexp)
BIC(mbest)


#jereme says you need to use LL instead of REML when doing model selsection. LEt's see if that matters. 


mbest = lmer(logBV ~ Temp+X2s+ CVPSWPs +(1|Year)+ (1|DOY) + (1|Station),
             data = micro, na.action = "na.fail", REML = FALSE)
mx2 = lmer(logBV ~ Temp+X2s+(1|DOY)+ (1|Year)+ (1|Station),
           data = micro, na.action = "na.fail", REML = FALSE)
mexp = lmer(logBV ~ Temp+CVPSWPs+(1|Year)+ (1|DOY) + (1|Station),
            data = micro, na.action = "na.fail", REML = FALSE)

aictab(list(mx2, mexp, mbest))
bictab(list(mx2, mexp, mbest))

#yeah, let's go with that.
#################################################################################
#try again and look at more different flow parameters

#calculate spring outflow
DFspring = DF %>%
  mutate(Month = month(Date), Year = year(Date)) %>%
  filter(Month %in% c(2,3,4,5)) %>%
  group_by(Year) %>%
  summarize(OUTm = mean(OUT), INm = mean(SJR + SAC), SJRm = mean(SJR), SACm = mean(SAC)) %>%
  ungroup() %>%
  mutate(OUTms = scale(OUTm), INms = scale(INm), 
         SJRms = scale(SJRm), SACms = scale(SACm), Year = as.factor(Year))

micro = mutate(micro, SACs = scale(SAC), SJRs = scale(SJR)) %>%
  left_join(DFspring)
mglobal3 = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+ Year+(1|Station),
                data = micro, na.action = "na.fail")

dredge(mglobal3)
vif(mglobal3)


mglobal4 = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+INms + Year+(1|Station),
                data = micro, na.action = "na.fail", REML = FALSE)

mglobal4a = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+SACms+Year+(1|Station),
                data = micro, na.action = "na.fail", REML = FALSE)

mglobal4b = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+ SJRms+Year+(1|Station),
                 data = micro, na.action = "na.fail", REML = FALSE)

mods1 = dredge(mglobal4)
mods2 = dredge(mglobal4a)
mods3 = dredge(mglobal4b)


mX2= lmer(logBV ~ Temp +X2+(1|Year)+ (1|Station),
                data = micro, na.action = "na.fail")

mX2a= lmer(logBV ~ Temp +X2+OUTms+  (1|Station),
          data = micro, na.action = "na.fail")

AICc(mX2)
AICc(mX2a)

#SO, spring outflow and inflow are about the same. 
#they are really good predictors if you don't include year as a random effect.
#san joaquin is a better predictor than Sacramento flow. Exports are also good.
#X2 is better than all the combined flow parameters


mglobal5 = lmer(logBV ~ Temp +X2+CVPSWPs+ SACs+SJRs+INms+(1|Station),
                data = micro, na.action = "na.fail")

dredge(mglobal5)
vif(mglobal5)
#yuck, can't have all those in there

#let's see if there are any combos that work

m1 = lmer(logBV ~ Temp +X2s+ SACs+SJRs+Year+(1|Station),
          data = micro, na.action = "na.fail")
m1a = lmer(logBV ~ Temp+CVPSWPs +SACs+SJRs+INms+Year+(1|Station),
          data = micro, na.action = "na.fail")
m2 = lmer(logBV ~ Temp +X2s+CVPSWPs+ SJRs+INms+(1|Station),
          data = micro, na.action = "na.fail")
m3 = lmer(logBV ~ Temp +X2s+CVPSWPs+ SACs+INms+(1|Station),
          data = micro, na.action = "na.fail")
m4 = lmer(logBV ~ Temp + SACs+SJRs+INms+(1|Station),
          data = micro, na.action = "na.fail")
m5 = lmer(logBV ~ Temp +X2s+SJRs+INms+(1|Station),
          data = micro, na.action = "na.fail")
m6 = lmer(logBV ~ Temp +CVPSWPs+ SACs+INms+(1|Station),
          data = micro, na.action = "na.fail")
m7 = lmer(logBV ~ Temp +CVPSWPs+ SJRs+INms+(1|Station),
          data = micro, na.action = "na.fail")
m8 = lmer(logBV ~ Temp +X2s+SJRs+INms+(1|Station),
          data = micro, na.action = "na.fail")
m9 = lmer(logBV ~ Temp +X2s+SACs+INms+(1|Station),
          data = micro, na.action = "na.fail")
m10 = lmer(logBV ~ Temp +X2s+CVPSWPs+ SACs+SJRs+(1|Station),
           data = micro, na.action = "na.fail")
m11 = lmer(logBV ~ Temp+CVPSWPs +X2s+ SACs+SJRs+INms+(1|Station),
          data = micro, na.action = "na.fail")
m12 = lmer(logBV ~ Temp+CVPSWPs +X2s+INms+(1|Station),
           data = micro, na.action = "na.fail")
m13 = lmer(logBV ~ Temp+X2s+ SACs+INms+(1|Station),
           data = micro, na.action = "na.fail")
m14 = lmer(logBV ~ Temp+X2s+ SJRs+INms+(1|Station),
           data = micro, na.action = "na.fail")

vif(m1) #no
vif(m1a) #yes
vif(m2) #no
vif(m3) #no
vif(m4) #yes - no x2
vif(m5) #no
vif(m6) #yes - no x2
vif(m7)#yes - no x2
vif(m8) # no
vif(m9) # no
vif(m10) #no
vif(m11) #no
vif(m12) #yes
vif(m13) #maybe
vif(m14) #no

#so, here are all the moddels that are legit
m1a = lmer(logBV ~ Temp+CVPSWPs +SACs+SJRs+INms+Year+(1|Station),
           data = micro, na.action = "na.fail", REML = FALSE)
m12 = lmer(logBV ~ Temp+CVPSWPs +X2s+INms+Year+(1|Station),
           data = micro, na.action = "na.fail", REML = FALSE)
m12a = lmer(logBV ~ Temp+CVPSWPs +X2s+SACms+Year+(1|Station),
           data = micro, na.action = "na.fail", REML = FALSE)
m12b = lmer(logBV ~ Temp+CVPSWPs +X2s+SJRms+Year+(1|Station),
            data = micro, na.action = "na.fail", REML = FALSE)

mglobal3 = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+ (1|Station),
                data = micro, na.action = "na.fail", REML = FALSE)
mglobal4 = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+INms+(1|Station),
                data = micro, na.action = "na.fail", REML = FALSE)


mods1 = dredge(mglobal4, rank = "BIC")
mods2 = dredge(mglobal4a, rank = "BIC")
mods3 = dredge(mglobal4b, rank = "BIC")
mods4 = dredge(m12, rank = "BIC")
mods5 = dredge(m12a, rank = "BIC")
mods6 = dredge(m12b, rank = "BIC")

allmods = bind_rows(mods1, mods2, mods3, mods4, mods5, mods6) 

write.csv(as.data.frame(allmods), file = "outputs/AllMicroModels.csv", row.names = FALSE)


#this is what happens with year in there

bestwyear = lmer(logBV ~ Temp + CVPSWPs+ SACs+ SJRs + Year + (1|Station), data = micro)
vif(bestwyear)
summary(bestwyear)
plot(allEffects(bestwyear))
#ok, there should not be an increase in microcystis with increased sac flow. That has to be an issue with colinearity
plot(simulateResiduals(bestwyear))



bestwyeart = lmer(logBV ~ Temp + CVPSWPs+ X2s + Year + (1|Station), data = micro)
vif(bestwyeart)
summary(bestwyeart)
plot(allEffects(bestwyeart))


#this is the best model
bestmod = lmer(logBV ~ Temp+X2s+SACms+(1|Station),
               data = micro, na.action = "na.fail")
summary(bestmod)
vif(bestmod)
plot(allEffects(bestmod))
plot(simulateResiduals(bestmod))


#second best model was with total inflow instead of sac inflow, and i think i'll just skip that one
bestmod2 = lmer(logBV ~ Temp+CVPSWPs+X2s+SACms+(1|Station),
               data = micro, na.action = "na.fail")
summary(bestmod2)
vif(bestmod2)
plot(allEffects(bestmod2))
plot(simulateResiduals(bestmod2))
#OK, now exports aren't even signficant. That's surprising

######################################################################################

#what does a 3km change in x2 get you?
X2scaled = (rep(c(60,63,66,69,72,75, 78, 81, 84, 87, 90), 6)-mean(micro$X2))/sd(micro$X2)
Tempscaled = data.frame(Temp = (c(15,17,19,21,23,25,27)-mean(micro$field.Water.temp))/sd(micro$field.Water.temp),
                        TempNotScaled = c(15,17,19,21,23,25,27))

CVPSWPscaled = (rep(seq(500, 10000, length.out=11), each =66)-mean(micro$CVPSWP))/sd(micro$CVPSWP)
SACscaled = (rep(seq(9000, 70000, length.out=11), each =66)-mean(micro$SACm))/sd(micro$SACm)

newdata = data.frame(DOY = 254,
                     Year = rep(c(rep(2014, 11), rep(2015, 11), rep(2016, 11), 
                              rep(2017, 11),rep(2018, 11),rep(2019, 11)),11),
                     Station = "CV D4",
                     X2s = rep(X2scaled,11),
                     X2 = rep(rep(c(60,63,66,69,72,75, 78, 81, 84, 87, 90), 6),11),
                     CVPSWPs =CVPSWPscaled,
                     CVPSWP = rep(seq(500, 10000, length.out=11), each =66),
                     SACms =SACscaled,
                     SACm = rep(seq(9000, 70000, length.out=11), each =66))
                     
test = merge(newdata, Tempscaled)

test$predictions = predict(bestmod, newdata = test)


newdata = mutate(test, BV = exp(predictions), lagBV = lag(BV), percent = (BV+lagBV)/lagBV)


############################################################################
#what are the range in X2 you get every year?
DFyear = DF %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  filter(Month %in%c(6:12), Year %in%c(2014:2019))%>%
  group_by(Year) %>%
  summarize(MinX2 = min(X2, na.rm =T), MaxX2 = max(X2, na.rm =T), MeanX2 = mean(X2, na.rm =T),
            MinExp = min((CVP+SWP), na.rm =T), MaxExp = max((CVP+SWP), na.rm =T))


yrtypes = data.frame(Year = c(2014:2019), yeartype = c("Critical", "Critical", "Below Normal", "Wet", "Dry", "Wet"))
#add water year type to these.
#biovolume is 'micrometer cubed per leter.
# micrometerCubedPerLiter

p = ggplot(newdata, aes(x = X2, y = BV)) + 
  geom_smooth(aes(color = as.factor(TempNotScaled)), se = FALSE)+
  scale_color_brewer(name = "Temperature", palette = "OrRd")+
  #geom_text(data = yrtypes, aes(x = 70, y = 230000, label = yeartype))+
  ylab("Biovolume (um3/L)")+
  xlab("X2 (km)")+
  facet_wrap(~SACm)+
  theme_bw()

p
# 
# library(grid)
# 
# g <- ggplot_gtable(ggplot_build(p))
# stripr <- which(grepl('strip-t', g$layout$name))
# fills <- c("skyblue","tan1","skyblue",  "indianred", "indianred","gold")
# k <- 1
# 
# for (i in stripr) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid.draw(g)



p2 = ggplot(newdata, aes(x = SACms, y = BV)) + 
  geom_smooth(aes(color = as.factor(TempNotScaled)), se = FALSE)+
  scale_color_brewer(name = "Temperature", palette = "OrRd")+
  ylab("Biovolume (um3/L)")+
  facet_wrap(~X2)+
  xlab("Average Spring (Feb-May)\n Sacramento River Inflow (cfs)")+
  #geom_text(data = yrtypes, aes(x = 6000, y = 180000, label = yeartype))+
  theme_bw()
  #geom_rect(data = DFyear, aes(ymin = 0, ymax = 200000, xmin = MinExp, xmax =  MaxExp), 
  #          inherit.aes = F, alpha = 0.3)
p2

ggplot(newdata, aes(x = SACm, y = log(BV))) + 
  geom_smooth(aes(color = as.factor(TempNotScaled)), se = FALSE)+
  scale_color_brewer(name = "Temperature", palette = "OrRd")+
  ylab("Log Biovolume (um3/L)")+
  facet_wrap(~X2)+
  xlab("Average Spring (Feb-May)\n Sacramento River Inflow (cfs)")+
  #geom_text(data = yrtypes, aes(x = 6000, y = 180000, label = yeartype))+
  theme_bw()


p3 = ggplot(newdata, aes(x = TempNotScaled, y = BV)) + 
  geom_smooth(aes(color = as.factor(SACm)), se = FALSE)+
  scale_color_viridis_d(name = "Outflow", option = "C")+
  ylab("Biovolume (um3/L)")+
  facet_wrap(~X2)+
  xlab("Water Temperature")+
  #geom_text(data = yrtypes, aes(x = 6000, y = 180000, label = yeartype))+
  theme_bw()
#geom_rect(data = DFyear, aes(ymin = 0, ymax = 200000, xmin = MinExp, xmax =  MaxExp), 
#          inherit.aes = F, alpha = 0.3)
p3

# g2 <- ggplot_gtable(ggplot_build(p2))
# stripr <- which(grepl('strip-t', g2$layout$name))
# fills <- c("skyblue","tan1","skyblue",  "indianred", "indianred","gold")
# k <- 1
# 
# for (i in stripr) {
#   j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid.draw(g2)


ggplot(filter(ungroup(newdata), CVPSWP == median(newdata$CVPSWP)), 
       aes(x = X2s, y = predictions, color = as.factor(Year))) + geom_point() + geom_line()+
  facet_wrap(~TempNotScaled)
ggplot(newdata, aes(x = X2s, y = percent, color = as.factor(Year))) + geom_point() + geom_line()+
  facet_wrap(~TempNotScaled)
#so we still have a very strong effect of X2, and it works better than outflow, which is interesting
library(effects)
plot(allEffects(mbest))

unique(micro$EMP_Site)
library(visreg)
#make all the partial residual plots
p2 = visreg(mbest, gg = T)


ptemp = p2[[1]]
pX2 = p2[[2]]
pexp = p2[[3]]
preds = ptemp$data
preds = bind_cols(micro, preds)

ggplot(preds, aes(x = field.Water.temp, y = y))+
  geom_point(aes(color = Year, shape = Year))+
  scale_color_viridis_d(option = "B")+
  geom_smooth(method = "lm")+
  theme_bw()+
  ylab("log-transforme Microcystis biovolume\npartial residuals")+
  xlab("Water temperature")

ggsave("plots/MicTempParResid.tiff", device = "tiff", width =6, height =5.5)

preds2 = pX2$data
preds2 = bind_cols(micro, preds2)

ggplot(preds2, aes(x = X2, y = y))+
  geom_point(aes(color = Year,  shape = Year))+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_color_viridis_d(option = "B")+
  ylab("log-transforme Microcystis biovolume\npartial residuals")+
  xlab("X2 (km)")

ggsave("plots/MicX2ParResid.tiff", device = "tiff", width =6, height =5.5)

preds3 = pexp$data
preds3 = bind_cols(micro, preds3)

ggplot(preds3, aes(x = CVPSWP, y = y))+
  geom_point(aes(color = Year,  shape = Year))+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_color_viridis_d(option = "B")+
  ylab("log-transforme Microcystis biovolume\npartial residuals")+
  xlab("CVP+SWP exports (cfs)")

ggsave("plots/MicExpParResid.tiff", device = "tiff", width =6, height =5.5)
#R-squared
library(MuMIn)
r.squaredGLMM(mbest)

predsall = bind_rows(mutate(preds, parameter = "Temperature (C)", Value = field.Water.temp),
                     mutate(preds2, parameter = "X2 (km)", Value = X2),
                     mutate(preds3, parameter = "Project Exports (cfs)", Value = CVPSWP))

ggplot(predsall, aes(x =Value, y = y))+
  geom_point(aes(color = Year,  shape = Year))+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_color_viridis_d(option = "B")+
  facet_wrap(~parameter, scales = "free_x")+
  ylab("log-transforme Microcystis biovolume\npartial residuals")

ggsave("plots/PartialResid.tiff", device = "tiff", width =8, height =5)

############################################################################

#ok, now just do stations in the South Delta/San JOaquin

microSJR = filter(micro, Station %in% c("NA D26", "SJ NA", "NA MD10A", "MI NA", "RR P8", "RR NA", "VC NA", "OR"))


mbestSJR = lmer(logBV ~ Temp+X2s+(1|Year)+ (1|DOY) + (1|Station),
             data = microSJR, na.action = "na.fail")
plot(allEffects(mbestSJR))
plot(mbestSJR)
mbestr = simulateResiduals(mbestSJR)
acf(residuals(mbestSJR))
summary(mbestSJR)
r.squaredGLMM(mbestSJR)
fooSJR = summary(mbestSJR)
write.csv(fooSJR$coefficients, "outputs/X2modelSJR.csv")
write.csv(fooSJR$varcor, "outputs/X2modelrandomSJR.csv")

############################################################################33
#look by region
library(deltamapr)
library(sf)

regions = R_EDSM_Regions_1617P1%>%
  st_transform(crs = 4326)

microsf = st_as_sf(micro, coords = c("Lon", "Lat"), crs = 4326) %>%
  st_join(regions)

ggplot()+
  geom_sf(data = WW_Delta)+
  geom_sf(data = microsf, aes(color = Region))

microsf = mutate(microsf, Region = case_when(Region == "West" ~ "North",
                                             TRUE ~ Region))


stations = select(microsf, Station, Region, Survey_Year) %>%
  distinct() %>%
  group_by(Station, Region) %>%
  summarize(years = length(unique(Survey_Year))) %>%
  mutate(Station = str_remove(Station, "NA "),
         Station = str_remove(Station, " NA")) %>%
  filter(Station != "RR")

stations2 = select(micro, Station,Survey_Year, Lat, Lon) %>%
  distinct() %>%
  group_by(Station) %>%
  summarize(years = length(unique(Survey_Year)), Latitude = first(Lat), Longitude = first(Lon)) %>%
  mutate(Station = str_remove(Station, "NA "),
         Station = str_remove(Station, " NA")) %>%
  filter(Station != "RR")

save(microsf, stations, stations2, file ="data/stations.RData")

ggplot()+
  geom_sf(data = WW_Delta)+
  geom_sf(data = stations, aes(color = Region))+
  geom_sf_label(data = stations, aes(label = Station))+
  coord_sf(xlim = c(-122.1, -121.2), ylim = c(37.8, 38.2))

mbest2 = lmer(logBV ~ Temp+X2s*Region+(1|Year)+ (1|DOY),
             data = microsf, na.action = "na.fail")


summary(mbest2)
plot(allEffects(mbest2))

X2scaled = (rep(c(60,63,66,69,72,75, 78, 81, 84, 87, 90), 12)-mean(micro$X2))/sd(micro$X2)
newdata = data.frame(Temp = mean(micro$Temp, na.rm = T), 
                     bryte.NH4.mgL = mean(micro$bryte.NH4.mgL, na.rm =T),
                     bryte.NO3.mgL = mean(micro$bryte.NO3.mgL, na.rm =T),
                     Region = c(rep("North", 66), rep("South", 66)),
                     DOY = 192,
                     Year = c(rep(2014, 11), rep(2015, 11), rep(2016, 11), 
                              rep(2017, 11),rep(2018, 11),rep(2019, 11), rep(2014, 11), rep(2015, 11), rep(2016, 11), 
                              rep(2017, 11),rep(2018, 11),rep(2019, 11)),
                     X2s = X2scaled,
                     X2 = rep(c(60,63,66,69,72,75, 78, 81, 84, 87, 90), 12))

newdata$predictions = predict(mbest2, newdata = newdata)

newdata = mutate(newdata, BV = exp(predictions), lagBV = lag(BV), percent = (BV+lagBV)/lagBV)

ggplot(newdata, aes(x = X2, y = BV, color = as.factor(Year))) + geom_point() + geom_line()+
  facet_wrap(~Region)
ggplot(newdata, aes(x = X2, y = predictions, color = as.factor(Year))) + geom_point() + geom_line()+
  facet_wrap(~Region)
ggplot(newdata, aes(x = X2, y = percent, color = as.factor(Year))) + geom_point() + geom_line()+
  facet_wrap(~Region)


ggplot(micro, aes(x = X2, y = field.Water.temp))+ geom_point()+ geom_smooth()


#################################################################################################

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

X2s = read_excel("data/supplemental_data_wr.1943-5452.0000617_hutton3.xlsx", skip =1)
X2s = rename(X2s, X2 = `Sacramento River X2 Position (km from GGB)`) %>%
  mutate(Date = date(Date), Month = month(Date), Year = year(Date)) %>%
  select(Date, Month, Year, X2) 

#more recent X2 from Dayflow
X2s2 = DF %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  group_by(Year, Month) %>%
  summarize(X2 = mean(X2, na.rm =T)) %>%
  filter(Year >2012)

X2sall = bind_rows(X2s, X2s2)

ggplot(X2sall, aes(x = Date, y = X2))+geom_line()
x2017 = filter(X2sall, Year == 2017)

ggplot(dplyr::filter(DF, year(Date) == 2017), aes(x = Date, y = X2))+geom_line()

restime = full_join(resave, X2sall, by = c("Year", "Month")) %>%
  filter(!is.na(Region))

mypal = c(brewer.pal(8, "Dark2"), brewer.pal(10, "Set2"))
ggplot(restime, aes(x = X2, y = Restime, color = Month))+ geom_point()+
  geom_smooth(method = "lm")+
  scale_color_viridis_c()+
  facet_wrap(~Region, scales = "free_y")
#wow, that's not good at all
#but i need to see whether Joey used historical operations or current operations. THis might not actually be appropriate.


mypal = c(brewer.pal(8, "Dark2"), brewer.pal(10, "Set2"))
ggplot(filter(restime, Year>1990), aes(x = X2, y = Restime, color = Month))+ geom_point()+
  geom_smooth(method = "lm")+
  scale_color_viridis_c()+
  facet_wrap(~Region, scales = "free_y")
##############################################################################################
#What if we use the residence time data set I made from Hammock et al for the drought paper?

load("data/ResidenceTime.RData")
RTs = left_join(DFRTall, X2sall, by = c("Year", "Month2"="Month"))

ggplot(RTs, aes(x = X2, y = log(SACRT), color = Month2))+ geom_point()+
  scale_color_viridis_c(name = "Month\nof year")+
  theme_bw()+
  geom_smooth(method = "lm")+ ylab("log-transformed Sacramento \nResidence Time (Days)")
#OK, that's pretty good

ggsave("plots/SacRTvX2.tiff", device = "tiff", width = 5, height = 4)

ggplot(RTs, aes(x = X2, y = log(SJRT)))+ geom_point(aes(color = Month2))+
  scale_color_viridis_c(name = "Month\nof year")+
  geom_smooth(method = "lm")+ ylab("log-transformed San Joaquin\n Residence Time (Days)")+
  theme_bw()
#Much less good. 

ggsave("plots/SJRTvX2.tiff", device = "tiff", width = 5, height = 4)

#now let's actually run a GAM on it

library(mgcv)

gam1 = gam(SJRT ~ s(X2), data = RTs)
summary(gam1)
gam.check(gam1)

plot(gam1)


gam2 = gam(SACRT ~ s(X2), data = RTs)
summary(gam2)
gam.check(gam2)

plot(gam2)

#What if I just did a linear model? Maybe simpler?
library(lmer)
library(rsq)
lmsj = lmer(log(SJRT) ~X2 + (1|Year)+ (1|Month), data = RTs)
summary(lmsj)
acf(residuals(lmsj))
r.squaredGLMM(lmsj)
plot(lmsj)
plot(allEffects(lmsj))

lmsac = lmer(log(SACRT) ~X2  + (1|Year)+ (1|Month), data = RTs)
summary(lmsac)
r.squaredGLMM(lmsac)

plot(lmsac)

sj = as.data.frame(summary(lmsj)$coefficients) %>%
  mutate(River = "San Joaquin")
sac = as.data.frame(summary(lmsac)$coefficients) %>%
  mutate(River = "Sacramento")

lms = bind_rows(sj, sac)
write.csv(lms, "outputs/residenceTimeModels.csv")
########################################################################
#check residence time versus some other factors
ggplot(RTs, aes(x = pump, y = SJRT))+ geom_point(aes(color = Month2))+
  scale_color_viridis_c()+
  geom_smooth()+ ylab("San Joaquin Residence Time (Days)")+
  theme_bw()
#I mean, pumping's part of the residence time model, so this isn't surprising, butit's kinda neat.

#what about outflow?
RTs2 = left_join(RTs, select(DF,Date, OUT)) 

ggplot(RTs2, aes(x = log(OUT), y = SJRT))+ geom_point(aes(color = Month2))+
  scale_color_viridis_c()+
  geom_smooth()+ ylab("San Joaquin Residence Time (Days)")+
  theme_bw()

ggplot(RTs2, aes(x = log(OUT), y = SACRT))+ geom_point(aes(color = Month2))+
  scale_color_viridis_c()+
  geom_smooth()+ ylab("Sacramento Residence Time (Days)")+
  theme_bw()

ggplot(filter(DF, OUT>0), aes(x = log(OUT), y = CVP+SWP))+ geom_point(aes(color = year(Date)))+
  scale_color_viridis_c()+
  geom_smooth()+ ylab("CVP+SWP pumping")+
  theme_bw()


#######################################################################################
#visual index versus X2

HABsReg = st_drop_geometry(HABsReg)

HAB2 = HABsReg %>%
  st_drop_geometry() %>%
  filter(!is.na(Microcystis)) %>%
  group_by(Year, Month, Region) %>%
  summarize(Mic = mean(Microcystis, na.rm =T)) %>%
  filter(!is.na(Region))

#X2 from dayflow
DFHAB = left_join(HABsReg, DF)%>%
  filter(X2>60)

HAB2 = left_join(HAB2, RTs, by = c("Year", "Month"="Month2")) 

ggplot(DFHAB, aes(x = X2, y = Microcystis, color = Region))+ geom_point()+ geom_smooth(method = "lm")


#just monthly data
ggplot(HAB2, aes(x = X2, y = Mic))+ geom_point(aes(color = Region))+ geom_smooth(method = "lm")

#Model with DOY and stuff
DFHAB = mutate(DFHAB, DOY = yday(Date), MicF = case_when(Microcystis ==1 ~ "Absent",
                                                         Microcystis %in% 2:3~ "Mid",
                                                         Microcystis %in% 4:5 ~ "High"),
               MicF = factor(MicF, levels = c("Absent", "Mid", "High"), ordered = T),
               Temp = scale(Temperature), X2a = scale(X2), DOYa = scale(DOY)) %>%
  filter(!is.na(Temp), !is.na(X2a))

m4 = clmm(MicF ~ Temp+ X2a+DOYa + (1|Station)+ (1|Year), data = DFHAB)
summary(m4)

#not as exciting as I'd think
#Also, i need to work on how to plot this.

DFHAB$fits = fitted(m4)

ggplot(DFHAB, aes( X2a, fits, color = MicF))+ geom_point()+ geom_smooth(method = "lm")+
  ylab("Probability")+ xlab("X2")

#Hm.

#is X2 a better indicator than local residence time?

DFHAB2 = left_join(DFHAB, resave2) %>%
  mutate(resscale = scale(ResTime)) %>%
  filter(!is.na(Region), !is.na(ResTime))

m5 = clmm(MicF ~ Temp+ resscale+X2+DOYa + (1|Station)+ (1|Year), data = DFHAB2)
summary(m5)


DFHAB2$fits = fitted(m5)

ggplot(DFHAB2, aes( X2, fits, color = MicF))+ geom_point()+ geom_smooth(method = "lm")+
  ylab("Probability")+ xlab("X2")

ggplot(DFHAB2, aes(ResTime, fits, color = MicF))+ geom_point()+ geom_smooth(method = "lm")+
  ylab("Probability")+ xlab("ResidenceTime")


#try an interaction with region
m6 = clmm(MicF ~ Temp+ X2a*Region +DOYa + (1|Station)+ (1|Year), data = DFHAB2)
summary(m6)


DFHAB2$fits = fitted(m6)

ggplot(DFHAB2, aes( X2, fits, color = MicF))+ geom_point()+ geom_smooth(method = "lm")+
  ylab("Probability")+ xlab("X2")+
  facet_wrap(~Region)

ggplot(DFHAB2, aes(ResTime, fits, color = MicF))+ geom_point()+ geom_smooth(method = "lm")+
  ylab("Probability")+ xlab("ResidenceTime")

#bigger regions?

DFHAB2 = mutate(DFHAB2, Region2 = case_when(RegionDSM %in% c("USAC","LSAC", "CACHE", 'ESUISUN') ~ "Sacramento",
                                              RegionDSM %in% c("CLIFTON", "LSJR", "MSJR", "FRANKS", "STOCKTON", "EAST", "MILDRED")~ "San Joaquin",
                                              RegionDSM %in% c("MONTZ", "WSUISUN") ~ "Suisun"))



#try an interaction with region
m7 = clmm(MicF ~ Temp+ X2a*Region2 +DOYa + (1|Station)+ (1|Year), data = DFHAB2)
summary(m7)


DFHAB2$fits = fitted(m7)

ggplot(DFHAB2, aes( X2, fits, color = MicF))+ geom_point()+ geom_smooth(method = "lm")+
  ylab("Probability")+ xlab("X2")+
  facet_wrap(~Region2)

ggplot(DFHAB2, aes(ResTime, fits, color = MicF))+ geom_point()+ geom_smooth(method = "lm")+
  ylab("Probability")+ xlab("ResidenceTime")

##########################################################################
#####################
#segmented regresstion - is there a threshold at 85km?
library(segmented)

mbest = lmer(logBV ~ Temp+X2s+ CVPSWPs +(1|Year)+ (1|DOY) + (1|Station),
             data = micro, na.action = "na.fail")
summary(mbest)

mbestlme = lme(logBV ~ Temp+CVPSWPs+X2s, random = list(Year =~1|Year),
             data = micro)

summary(mbestlme)
seg_fit <- segmented(mbestlme, 
                     seg.Z=~X2s, 
                     data=micro) 
summary(seg_fit)
#Ugh, this is awful.

#let's try a janky way

micro = mutate(micro, X2bin = case_when(X2 >=85 ~ "High",
                                        X2 < 85 ~ "Low"))

mbestbreak = lmer(logBV ~ Temp+X2s+ CVPSWPs+ X2bin +(1|Year)+ (1|DOY) + (1|Station),
             data = micro, na.action = "na.fail")
summary(mbestbreak)

BIC(mbest)
BIC(mbestbreak)

mbest = lmer(logBV ~ Temp+X2s+ CVPSWPs +(1|Year)+ (1|DOY) + (1|Station),
             data = micro, na.action = "na.fail", REML = FALSE)
mbestbreak = lmer(logBV ~ Temp+X2s*X2bin+ CVPSWPs +(1|Year)+ (1|DOY) + (1|Station),
                  data = micro, na.action = "na.fail", REML = FALSE)
BIC(mbest)
BIC(mbestbreak)

ggplot(micro, aes(x = X2, y = logBV))+ geom_point()+geom_smooth()
