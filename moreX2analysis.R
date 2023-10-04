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
  select(logBV, Temp, CVPSWPs, SACs, SJRs, OUTms, OUT,Outs, SJRms, SACms, INms, Year, 
         DOY, Collection_Date, X2s, CVPSWPs, Station, SAC, SJR)


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

smat <- outer(1:11, 1:11, vCorrelated, data = as.data.frame(micro[,c(2:6,8:11, 13,15)]))

nm <- colnames(micro[c(2:6,8:11, 13,15)])


dimnames(smat) <- list(nm, nm)
smat


mglobal4a = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+Outs+ INms+X2s+SJRms+SACms+DOY+X2s+(1|Year)+ (1|Station),
                data = micro, na.action = "na.fail", REML = FALSE)

modsSub <- dredge(mglobal4a, subset = smat)

#I think i'm going crazy

#########################################################################
#relax my correlation coefficient cuttoff


smat <- outer(1:11, 1:11, vCorrelated, data = as.data.frame(micro[,c(2:6,8:11, 13,15)]), cutoff = 0.7)

nm <- colnames(micro[c(2:6,8:11, 13,15)])

dimnames(smat) <- list(nm, nm)

modsSuba <- dredge(mglobal4, subset = smat, rank = "BIC")
modsSubb <- dredge(mglobal4a, subset = smat, rank = "BIC")

#every time I run this I get something different
######################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#this is the best version so far
#Month as a random effect?
#Definitely year and station as random effect. 

micro = mutate(micro, Month = month(Collection_Date)) %>%
  filter(Month !=6) 

smatx = smat[c(1:9,11), c(1:9,11)]
mglobal4x1 = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+Outs+ INms+X2s+SJRms+SACms+DOY+X2s+(1|Year)+ (1|Month)+(1|Station),
                 data = micro, na.action = "na.fail", REML = FALSE)


modsSubx1 <- dredge(mglobal4x1, subset = smat, rank = "BIC")
#And this is yet another answer
bx1a =get.models(modsSubx1, subset = 1)[[1]]
bx1a2 =get.models(modsSubx1, subset = 2)[[1]]

bx1a = update(bx1a, REML = TRUE)
summary(bx1a)
vif(bx1a)
plot(simulateResiduals(bx1a))
plot(allEffects(bx1a))
plot(allEffects(bx1a2))
visreg(bx1a)
r.squaredGLMM(bx1a)


foo1 = predict(bx1a)
micro$predictions = foo1

#look at this model in more detail. 
ggplot(micro, aes(x = logBV, y = predictions, color = Year))+ geom_point()+
  geom_abline(a = 1, b =0)

ggplot(micro, aes(x = Outs, y = predictions, color = Year))+ geom_point()

ggplot(micro, aes(x = INms, y = predictions, color = Year))+ geom_point()

ggplot(micro, aes(x = Temp, y = predictions, color = Year))+ geom_point()

ggplot(micro, aes(x = CVPSWPs, y = predictions, color = Year))+ geom_point()


#how critical is month as a random effect? what's my sample distribution like/

ggplot(micro, aes(x = Station, fill = as.factor(Month)))+
  facet_wrap(~Year)+
  geom_histogram(stat = "count")

ggplot(micro, aes(x = Month, fill = as.factor(Month)))+
  facet_wrap(~Year)+
  geom_histogram(stat = "count")
# 
# #how about month as a fixed effect?
# mglobal4x2 = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+Outs+ INms+X2s+SJRms+SACms+X2s+Month+(1|Year)+ (1|Station),
#                   data = micro, na.action = "na.fail", REML = FALSE)
# 
# smat2 = cbind(smat,rep(NA, 11))
# smat3x = rbind(smat2,rep(TRUE, 12))
# dimnames(smat3x)[[1]]=c(dimnames(smat)[[1]], "Month")
# dimnames(smat3x)[[2]]=c(dimnames(smat)[[2]], "Month")
# modsSubx2 <- dredge(mglobal4x2, subset = smat3x)
# bxa = get.models(modsSubx2, subset = 1)[[1]]
# bxb = get.models(modsSubx2, subset = 2)[[1]]
# plot(allEffects(bxa))
# plot(allEffects(bxb))
# vif(bxa)
# vif(bxb)
# 
# #but temperature and month have that colinerity problem
# 
# #yup, so lost.
# 
# #I have much more consistant sampling for 2016-2019. Let's try just those
# microsub =  filter(micro, Year %in% c("2016", "2017","2018", "2019"))
# mglobalsub = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+Outs+ INms+X2s+SJRms+SACms+X2s+Month+(1|Year)+ (1|Station),
#                   data =microsub, na.action = "na.fail", REML = FALSE)
# modssub2016 = dredge(mglobalsub, rank = "BIC")
# bx2 = get.models(modssub2016, subset = 1)[[1]]
# #boundary (singular) fit: see help('isSingular')
# plot(allEffects(bx2))
# 
# modssub2016b = dredge(mglobalsub, rank = "BIC", subset = smat3x)
# bx3 = get.models(modssub2016b, subset = 1)[[1]]
# plot(allEffects(bx3))
# #that makes so much sense.
# #I could also log-transform Sac flow. and maybe SJR and Outflow
# 
# # I could do all the things> this is super confusing.
# 
# #################################################################################
# #Now again, with log-transformed flow variables
# 
# microlog = mutate(micro, logSAC = log(SAC), logOUT = log(OUT), logSJR = log(SJR))
# 
# smatlog <- outer(1:10, 1:10, vCorrelated, data = as.data.frame(microlog[,c(2,3,9,10,11, 15,19:22)]), cutoff = 0.7)
# nm <- colnames(microlog[,c(2,3,9,10,11, 15,19:22)])
# 
# dimnames(smatlog) <- list(nm, nm)
# 
# mglobal4x3 = lmer(logBV ~ Temp +CVPSWPs+ logSAC+logSJR+logOUT+ INms+X2s+SJRms+SACms+X2s+(1|Year)+ (1|Month)+(1|Station),
#                   data = microlog, na.action = "na.fail", REML = FALSE)
# 
# 
# modsSubx3 <- dredge(mglobal4x3, subset = smatlog, rank = "BIC")
# #And this is yet another answer
# bx1 =get.models(modsSubx3, subset = 1)[[1]]
# bx1 = update(bx1, REML = TRUE)
# summary(bx1)
# vif(bx1)
# plot(simulateResiduals(bx1))
# plot(allEffects(bx1))
# visreg(bx1)
# #basically identical to the version without log-transformation. 
# 
# ###########################################################################################################
# #what about month as a factor?
# microlog = mutate(microlog, Month1 = as.factor(Month))
# mglobal4x4 = lmer(logBV ~ Temp +CVPSWPs+ logSAC+logSJR+logOUT+ INms+X2s+SJRms+SACms+X2s+Month1+(1|Year)+ (1|Station),
#                   data = microlog, na.action = "na.fail", REML = FALSE)
# 
# smatlog2 = cbind(smatlog,rep(NA, 10))
# 
# smatlog3x = rbind(smatlog2,rep(TRUE, 11))
# dimnames(smatlog3x)[[1]]=c(dimnames(smatlog)[[1]], "Month1")
# dimnames(smatlog3x)[[2]]=c(dimnames(smatlog)[[2]], "Month1")
# 
# modsSubx4 <- dredge(mglobal4x4, subset = smatlog3x, rank = "BIC")
# #And this is yet another answer
# bx1b =get.models(modsSubx4, subset = 1)[[1]]
# bx1b = update(bx1b, REML = TRUE)
# summary(bx1b)
# vif(bx1b)
# plot(simulateResiduals(bx1b))
# plot(allEffects(bx1b))
# visreg(bx1)
# #basically identical to the version without log-transformation. 
