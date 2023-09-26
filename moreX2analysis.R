#Microcystis analysis
#Rosemary Hartman
#9/26/2023


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
#oh, that looks quite lovely!

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


#now set up a linear regression 
micro = mutate(micro, logBV = log(dwr.MIC_totbvL/100000+1), Year = as.factor(Survey_Year), DOY = yday(Collection_Date)) %>%
  filter(!is.na(field.Water.temp), !is.na(X2), !is.na(OUT)) %>%
  ungroup()

mglobal = lmer(logBV ~ field.Water.temp + X2 + OUT + DOY+CVPSWP+ field.Salinity+(1|Year),
               data = micro, na.action = "na.fail")

#Warning messages:Some predictor variables are on very different scales: consider rescaling

#rescale variables
micro = mutate(micro, Temp = scale(field.Water.temp), X2s = scale(X2), Outs = scale(log(OUT)),
               CVPSWPs = scale(CVPSWP),
               Salinitys = scale(field.Salinity), 
               DWR_Site = str_remove(DWR_Site, "rep"),
               Station = paste(DWR_Site, EMP_Site),
               Station = case_when(Station %in% c("OR D28A", "OR NA") ~ "OR",
                                   TRUE ~ Station))

#first global model
mglobal = lmer(logBV ~ Temp + X2s + Outs+ DOY+CVPSWPs+ Salinitys+ (1|Year) + (1|Station),
               data = micro, na.action = "na.fail")


vif(mglobal)
#DOY and temperature has highest vif()

mglobal2 = lmer(logBV ~ Temp + X2s + Outs+CVPSWPs+ + Salinitys+(1|Year),
                data = micro, na.action = "na.fail", REML = FALSE)

vif(mglobal2)
#much better

#Evaluate all combinations of these variables

dredge(mglobal2, rank = "BIC")

#the best model has Exports, Outflow, Temperature, and X2


#but outflow and X2 are probably tooo highly correlated. I wonder why that didn't show up in the vif test?
ggplot(micro, aes(x = X2s, y = Outs)) + geom_point() + geom_smooth(method = "lm")
ggplot(micro, aes(x = X2s, y = CVPSWPs)) + geom_point() + geom_smooth(method = "lm")


#try again and look at more different flow parameters

#calculate spring outflow
DFspring = DF %>%
  mutate(Month = month(Date), Year = year(Date)) %>%
  filter(Month %in% c(2,3,4,5)) %>%
  group_by(Year) %>%
  summarize(OUTm = mean(OUT), INm = mean(SJR + SAC), SJRm = mean(SJR), SACm = mean(SAC)) %>%
  ungroup()%>%
  mutate(OUTms = scale(OUTm), INms = scale(INm), 
         SJRms = scale(SJRm), SACms = scale(SACm), Year = as.factor(Year))

micro = mutate(micro, SACs = scale(SAC), SJRs = scale(SJR)) %>%
  left_join(DFspring) %>%
  select(logBV, Temp, CVPSWPs, SACs, SJRs, OUTms, Outs, SJRms, SACms, INms, Year, DOY, Collection_Date, X2s, CVPSWPs, Station)

mglobal3 = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+SACms+X2s+ (1|Station),
                data = micro, na.action = "na.fail")

vif(mglobal3)
#Sacramento river flow and X2 can't be used together.

mglobal3a = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+SACms+ (1|Station),
                 data = micro, na.action = "na.fail", REML = FALSE)
dredge(mglobal3a, rank = "BIC")

mglobal3b = lmer(logBV ~ Temp +CVPSWPs+ X2s+SACms+ (1|Station),
                 data = micro, na.action = "na.fail", REML = FALSE)

dredge(mglobal3b, rank = "BIC")

#now without spring outflow and with an effect of year


mglobal3c = lmer(logBV ~ Temp +CVPSWPs+ X2s+Year+ (1|Station),
                 data = micro, na.action = "na.fail", REML = FALSE)

dredge(mglobal3c, rank = "BIC")


mglobal3d = lmer(logBV ~ Temp +CVPSWPs+SACs + SJRs +Year+ (1|Station),
                 data = micro, na.action = "na.fail", REML = FALSE)

dredge(mglobal3d, rank = "BIC")
#but why does this say that sacramento flow is positively correlated with microcystis? That make sno sense. 

#what else is highly correlated?
ggplot(micro, aes(x = SACs, y = Outs)) + geom_point() + geom_smooth(method = "lm") #highly correlated
ggplot(micro, aes(x = SACs, y = CVPSWPs)) + geom_point() + geom_smooth(method = "lm") #not too bad
ggplot(micro, aes(x = SACs, y = X2s)) + geom_point() + geom_smooth(method = "lm")#highly correlated
ggplot(micro, aes(x = CVPSWPs, y =  Outs)) + geom_point() + geom_smooth(method = "lm")#not too bad
ggplot(micro, aes(x = SACms, y = X2s)) + geom_point() + geom_smooth(method = "lm")#not terrible, not great
ggplot(micro, aes(x = SACms, y = Outs)) + geom_point() + geom_smooth(method = "lm")#not terrible, not great
ggplot(micro, aes(x = SACms, y = CVPSWPs)) + geom_point() + geom_smooth(method = "lm")#not too bad
ggplot(micro, aes(x = X2s, y = CVPSWPs)) + geom_point() + geom_smooth(method = "lm")#not too bad
ggplot(micro, aes(x = SJRs, y = CVPSWPs)) + geom_point() + geom_smooth(method = "lm")#not too bad, not great tho

#I'm starting to feel like I should only do models with one flow variable 

#let's calculate the correlation coefficient for all the variables, maybe that is better than vif
is.correlated <- function(i, j, data, conf.level = 0.95, cutoff = 0.5, ...) {
  	if(j >= i) return(NA)
   	ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
   	ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}

 # Need vectorized function to use with 'outer'
vCorrelated <- Vectorize(is.correlated, c("i", "j"))

smat <- outer(1:11, 1:11, vCorrelated, data = as.data.frame(micro[,c(2:10, 12,14)]))

nm <- colnames(micro[c(2:10, 12,14)])

dimnames(smat) <- list(nm, nm)
smat

#So, DOY and temperature can't be used together. 

mglobal4 = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+Outs+ INms+X2s+SJRms+SACms+DOY+X2s+ (1|Station),
                data = micro, na.action = "na.fail", REML = FALSE)

modsSub <- dredge(mglobal4, subset = smat)

#Now try adding year 

mglobal4a = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+Outs+ INms+X2s+SJRms+SACms+DOY+X2s+(1|Year)+ (1|Station),
                data = micro, na.action = "na.fail", REML = FALSE)

modsSub <- dredge(mglobal4a, subset = smat)

#I think i'm going crazy

#########################################################################
#relax my correlation coefficient cuttoff


smat <- outer(1:11, 1:11, vCorrelated, data = as.data.frame(micro[,c(2:10, 12,14)]), cutoff = 0.6)

nm <- colnames(micro[c(2:10, 12,14)])

dimnames(smat) <- list(nm, nm)

modsSub <- dredge(mglobal4, subset = smat)
modsSub <- dredge(mglobal4a, subset = smat)

#every time I run this I get something different

