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
library(visreg)
library(sf)
library(corrplot)
library(ggcorrplot)

#now let's grab the microcystis data from online:
#microcystis = read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.1076.1&entityid=86da4465dde063c89afb0dc256bfa619")
#save(microcystis, file = "Data/microcystis.RData")
load("Data/microcystis.RData")
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


#THe dayflow dataset is a bit difficult to deal with when you first download it,
#so I have a version I got ealier and cleaned up
load("Dayflow1997_2021.RData")
str(DF)

#bind dayflow data to microcystis data
microX = left_join(microcystis, DF, by = c("Collection_Date" = "Date")) %>%
  mutate(CVPSWP = CVP + SWP)

#now add regions to the microcystis dataset
micro2<- microX %>%
  mutate(Year = Survey_Year, Month = month(Collection_Date)) 


#plot microcystis versus X2
ggplot(micro2, aes(x = X2, y = log(dwr.MIC_totbvL+1), 
                  color = as.factor(Survey_Year)))+
  geom_point()+
  ylab("Log-transformed biovolume")

ggplot(micro2, aes(x = X2, y = log(dwr.MIC_totbvL+1) ))+
  geom_smooth(method = "lm")+
  geom_point(aes(color = as.factor(DWR_Site)))

#now set up a linear regression 
micro = mutate(micro2, logBV = log(dwr.MIC_totbvL/100000+1), 
               Year = as.factor(Survey_Year), DOY = yday(Collection_Date)) %>%
  filter(!is.na(field.Water.temp)) %>%
  ungroup()
#Warning messages:Some predictor variables are on very different scales: consider rescaling

#rescale variables
micro = mutate(micro, Temp = scale(field.Water.temp), X2s = scale(X2), Outs = scale(log(OUT)),
               CVPSWPs = scale(CVPSWP), WESTs = scale(abs(WEST)),
               Salinitys = scale(field.Salinity), 
               DWR_Site = str_remove(DWR_Site, "rep"),
               Station = paste(DWR_Site, EMP_Site),
               Station = case_when(Station %in% c("OR D28A", "OR NA") ~ "OR",
                                   TRUE ~ Station))
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
         DOY, Collection_Date, X2s, CVPSWPs, Station, SAC, SJR, field.Water.temp, CVPSWP, SACm)


#what else is highly correlated?
ggplot(micro, aes(x = SACs, y = Outs)) + geom_point() + geom_smooth(method = "lm") #highly correlated
ggplot(micro, aes(x = SACs, y = CVPSWPs)) + geom_point() + geom_smooth(method = "lm") #not too bad
ggplot(micro, aes(x = SACs, y = X2s)) + geom_point() + geom_smooth(method = "lm")#highly correlated
ggplot(micro, aes(x = CVPSWPs, y =  Outs)) + geom_point() + geom_smooth(method = "lm")#not too bad
ggplot(micro, aes(x = SACms, y = X2s)) + geom_point() + geom_smooth(method = "lm")#not terrible, not great
ggplot(micro, aes(y = SACms, x = SJRs)) + geom_point(aes(color = as.factor(Year))) + geom_smooth(method = "lm")

ggplot(micro, aes(x = SACms, y = Outs)) + geom_point() + geom_smooth(method = "lm")#not terrible, not great
ggplot(micro, aes(x = SACms, y = CVPSWPs)) + geom_point() + geom_smooth(method = "lm")#not too bad
ggplot(micro, aes(x = X2s, y = CVPSWPs)) + geom_point() + geom_smooth(method = "lm")#not too bad
ggplot(micro, aes(x = SJRs, y = CVPSWPs)) + geom_point() + geom_smooth(method = "lm")#not too bad, not great tho
ggplot(micro, aes(x = abs(WEST), y = CVPSWPs)) + geom_point() + geom_smooth(method = "lm")#not too bad, not great tho
ggplot(micro, aes(x = abs(WEST), y = logBV)) + geom_point() + geom_smooth(method = "lm")#not too bad, not great tho

#I'm starting to feel like I should only do models with one flow variable 

#let's calculate the correlation coefficient for all the variables, maybe that is better than vif
is.correlated <- function(i, j, data, conf.level = 0.95, cutoff = 0.5, ...) {
  	if(j >= i) return(NA)
   	ct <- cor.test(data[, i], data[, j], conf.level = conf.level,...)
   	ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}

 # Need vectorized function to use with 'outer'
vCorrelated <- Vectorize(is.correlated, c("i", "j"))

#########################################################################
#relax my correlation coefficient cuttoff

micro = mutate(micro, Month = month(Collection_Date))%>%
  filter(Month !=6) 
microtest = as.data.frame(micro[,c(2:6,8:10, 15)])
smat <- outer(1:9, 1:9, vCorrelated, data = microtest, cutoff = 0.71)
nm <- colnames(microtest)
dimnames(smat) <- list(nm, nm)

#correlation coefficient matrix

Corrs = cor(microtest)
corp = corrplot.mixed(Corrs, upper = "ellipse")

#make a matrix of whethe ror not things are included
smatpvals = apply(smat,c(1,2), FUN = function(x){ if(is.na(x)) 1 else if(x) 0 else 2})

ggcorrplot(corp$corr, method = "square", type = "lower", lab = TRUE, p.mat = smatpvals)

corp2 = corp$corrPos
corpPs = mutate(as.data.frame(smatpvals), xName = row.names(smatpvals)) %>%
  pivot_longer(cols = c(Temp:X2s), names_to = "yName", values_to = "included") 
corp2 = left_join(corp2, corpPs) %>%
  filter(included !=1) %>%
  mutate(xName = factor(xName, levels = c("CVPSWPs", "SACs", "SJRs", "OUTms", "Outs", "SJRms", "SACms", "X2s"),
                        labels = c("CVP+SWP", "Sac", "SJR", "Spring Outflow", "Outflow", "Spring SJR", "Spring Sac", "X2")))

ggplot(corp2, aes(x = x, y = y))+ geom_tile(aes(fill = corr), color = "black")+
  geom_tile(data = filter(corp2, included ==2), fill = "grey", color = "black")+
  scale_alpha(range = c(1,0))+
  geom_text(aes(label = round(corr, 2)))+
  scale_fill_gradient2(name = "Correlation\nCoefficient")+
  scale_x_continuous(breaks = c(2:9), labels = c(levels(corp2$xName)))+
  scale_y_continuous(breaks = c(2:9), labels = rev(c("Temperature", "CVP+SWP", "Sac", "SJR", 
                                                     "Spring Outflow", "Outflow", "Spring SJR", "Spring Sac")))+
  ylab(NULL)+xlab(NULL)+ theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust =1))

ggsave("plots/corrplot.tiff", device = "tiff", width =6, height =5)
######################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#this is the best version so far
#Month as a random effect?
#Definitely year and station as random effect. 

# 
# #smatx = smat[c(1:9,11), c(1:9,11)]
# #I tried adding nitrogen for Shawn et al, but I don't like it.
# 
# mglobal4x1 = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+ Outs+ X2s+SJRms+SACms+(1|Year)+ (1|Month)+(1|Station),
#                    data = micro, na.action = "na.fail", REML = FALSE)
# 
# write.csv(smat, "Correlations.csv")
# modsSubx1 <- dredge(mglobal4x1, subset = smat, rank = "BIC")
# write.csv(as.data.frame(modsSubx1), "outputs/MicroModelsOct2023.csv", row.names = F)
# 
# #And this is yet another answer
# bx1a =get.models(modsSubx1, subset = 1)[[1]]
# bx1a = update(bx1a, REML = T)
# summary(bx1a)
# plot(allEffects(bx1a))
# 
# #The absolute value of QWEEST is now in the second-best model instead of spring sac flow. but it's only boarder-line significant
# #More close to significant than spring flow was though. 
# #honestly I don't love it because it's hard to explain. 
# bx1a3 =get.models(modsSubx1, subset = 2)[[1]]
# summary(bx1a3)
# plot(allEffects(bx1a3))
# 
# #Ellen wants the one with Sacramento spring flow. 
# bx1a2 =get.models(modsSubx1, subset = 3)[[1]]
# summary(bx1a2)
# 
# 
# bx1a = update(bx1a2, REML = TRUE)
# summary(bx1a2)
# vif(bx1a2)
# plot(simulateResiduals(bx1a2))
# plot(allEffects(bx1a2))
# plot(allEffects(bx1a))
# visreg(bx1a2)
# 
# r.squaredGLMM(bx1a)
# 
# 
# foo1 = predict(bx1a)
# micro$predictions = foo1
# 
# #look at this model in more detail. 
# ggplot(micro, aes(x = logBV, y = predictions, color = Year))+ geom_point()+
#   geom_abline(a = 1, b =0)
# 
# ggplot(micro, aes(x = Outs, y = predictions, color = Year))+ geom_point()
# 
# ggplot(micro, aes(x = INms, y = predictions, color = Year))+ geom_point()
# 
# ggplot(micro, aes(x = Temp, y = predictions, color = Year))+ geom_point()
# 
# ggplot(micro, aes(x = CVPSWPs, y = predictions, color = Year))+ geom_point()
# 
# 
# #how critical is month as a random effect? what's my sample distribution like/
# microcystis = mutate(microcystis, Month = as.factor(month(Collection_Date)),
#                      Year = year(Collection_Date),
#                      DWR_Site = str_remove(DWR_Site, "rep"),
#                      Station = paste(DWR_Site, EMP_Site),
#                      Station = case_when(Station %in% c("OR D28A", "OR NA") ~ "OR",
#                                          TRUE ~ Station)) %>%
#   filter(!is.na(Collection_Date))
# ggplot(microcystis, aes(x = Station, fill = Month))+
#   facet_wrap(~Year)+
#   geom_histogram(stat = "count")
# 
# ggplot(micro, aes(x = Month, fill = as.factor(Month)))+
#   facet_wrap(~Year)+
#   geom_histogram(stat = "count")
# 
## 
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
# 
# 
# 
# ######################################################################################
# 
# #better plots of the effects
# Tempscaled = data.frame(Temp = (c(15,17,19,21,23,25,27)-mean(micro$field.Water.temp))/sd(micro$field.Water.temp),
#                         TempNotScaled = c(15,17,19,21,23,25,27))
# 
# CS = data.frame(CVPSWPs= (seq(500, 12000, length.out=11)-mean(micro$CVPSWP))/sd(micro$CVPSWP),
#                 CVPSWP = seq(500, 12000, length.out=11))
# 
# SJ =data.frame(SJRs =  (c(150,200,250,300,400,500,1000,2000,3000,5000)-mean(micro$SJR, na.rm =T))/sd(micro$SJR, na.rm =T),
#                SJR = c(150,200,250,300,400,500,1000,2000,3000,5000))
# 
# # SAC =data.frame(SACms =  (c(10000,15000, 20000, 30000, 40000, 65000)-mean(micro$SACm, na.rm =T))/sd(micro$SACm, na.rm =T),
# #                SACm = c(10000,15000, 20000, 30000, 40000, 65000))
# 
# 
# newdata = data.frame(Month = 8,
#                      Year = 2014:2019,
#                      Station = "FT D19")
# 
# test = merge(newdata, Tempscaled)
# test2 = merge(test, CS) %>%
#   merge(SJ)
# test2$predictions = predict(bx1a, newdata = test2)
# 
# 
# newdata = mutate(test2, BV = exp(predictions), lagBV = lag(BV), percent = (BV+lagBV)/lagBV)
# 
# 
# ############################################################################
# #what are the range in X2 you get every year?
DFyear = DF %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  filter(Month ==8, Year %in%c(2014:2019))%>%
  group_by(Year) %>%
  summarize(MinX2 = min(X2, na.rm =T), MaxX2 = max(X2, na.rm =T), MeanX2 = mean(X2, na.rm =T),
            MinExp = min((CVP+SWP), na.rm =T), MaxExp = max((CVP+SWP), na.rm =T),
            MinSJR = min(SJR, na.rm =T), MaxSJR = max(SJR, na.rm =T))
# 
# newdata2 = left_join(newdata, DFyear) %>%
#   mutate(REALexp = case_when(CVPSWP>= MinExp & CVPSWP <= MaxExp ~TRUE,
#                           TRUE ~ FALSE),
#          REALSJR = case_when(SJR>= MinSJR & SJR <= MaxSJR ~TRUE,
#                              TRUE ~ FALSE))
# 
# newtest = filter(newdata2, REALexp, REALSJR)
# 
# yrtypes = data.frame(Year = c(2014:2019), yeartype = c("Critical", "Critical", "Below Normal", "Wet", "Dry", "Wet"))
# #add water year type to these.
# #biovolume is 'micrometer cubed per leter.
# # micrometerCubedPerLiter
# 
# p = ggplot(newdata, aes(x = CVPSWP, y = BV)) + 
#   geom_line(aes(color = as.factor(TempNotScaled)))+
#   scale_color_brewer(name = "Temperature", palette = "OrRd")+
#   #geom_text(data = yrtypes, aes(x = 70, y = 230000, label = yeartype))+
#   ylab("Biovolume (um3/L)")+
#   xlab("Exports (cfs)")+
#   facet_grid(SJR~Year)+
#   theme_bw()
# 
# p
# 
library(grid)
YearLablePlot = function(plot){
  g <- ggplot_gtable(ggplot_build(plot))
  stripr <- which(grepl('strip-t', g$layout$name))
  fills <- c("skyblue","tan1","skyblue","indianred", "indianred","gold")
  k <- 1

  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  grid.draw(g)
}
# 
# ################Plot just values that actually happened
# pt = ggplot(newtest, aes(x = CVPSWP, y = BV)) + 
#   geom_line(aes(color = as.factor(TempNotScaled)))+
#   scale_color_brewer(name = "Temperature", palette = "OrRd")+
#   #geom_text(data = yrtypes, aes(x = 70, y = 230000, label = yeartype))+
#   ylab("Biovolume (um3/L)")+
#   xlab("Exports (cfs)")+
#   facet_grid(SJR~Year)+
#   theme_bw()
# 
# YearLablePlot(pt)
# 
# #can I shade the values that actually happened instead?
# 
# p = ggplot(newdata2, aes(x = CVPSWP, y = BV)) + 
#   
#   geom_rect(data = filter(newdata2, REALSJR), aes(xmin = MinExp, xmax = MaxExp, ymin = 0, ymax = max(BV)),
#             alpha = 0.2, fill = "grey90")+
#   scale_color_brewer(name = "Temperature", palette = "OrRd")+
#   #geom_text(data = yrtypes, aes(x = 70, y = 230000, label = yeartype))+
#   geom_line(aes(color = as.factor(TempNotScaled)))+
#   ylab("Biovolume (um3/L)")+
#   xlab("Exports (cfs)")+
#   facet_grid(SJR~Year)+
#   theme_bw()
# 
# YearLablePlot(p)
# 
# 
# 
# #now plot different ways
# p2 = ggplot(filter(newdata, SJR == 500), aes(x = CVPSWP, y = BV)) + 
#   geom_line(aes(color = as.factor(TempNotScaled)))+
#   scale_color_brewer(name = "Temperature", palette = "OrRd")+
#   scale_x_continuous(breaks =c(0,4000, 8000,  12000), labels = c(0, 4,8,12))+
#   ylab("Biovolume (um3/L)")+
#   facet_grid(~Year)+
#   xlab("CVP/SWP Exports (thousand cfs)")+
#   #geom_text(data = yrtypes, aes(x = 6000, y = 180000, label = yeartype))+
#   theme_bw()+
# geom_rect(data = DFyear, aes(ymin = 0, ymax = 450000, xmin = MinExp, xmax =  MaxExp), 
#           inherit.aes = F, alpha = 0.3)+
#   theme(legend.position = "bottom")
# YearLablePlot(p2)
# ggsave("Plots/ExporstbyBVbyYear.tiff", device = "tiff", width = 7, height =5)
# 
# summary(micro$CVPSWP)
# p3 = ggplot(filter(newdata, CVPSWP == 5100), aes(x = SJR, y = BV)) + 
#   geom_line(aes(color = as.factor(TempNotScaled)))+
#   scale_color_brewer(name = "Temperature", palette = "OrRd")+
#   scale_x_continuous(breaks =c(0,1000, 2000, 3000, 4000, 5000), labels = c(0, 1,2,3,4,5))+
#   ylab("Biovolume (um3/L)")+
#   facet_wrap(~Year, nrow =1)+
#   xlab("San Joaquin Flow (thousand cfs)")+
#   #geom_text(data = yrtypes, aes(x = 6000, y = 180000, label = yeartype))+
#   theme_bw()+
# geom_rect(data = DFyear, aes(ymin = 0, ymax = 350000, xmin = MinSJR, xmax =  MaxSJR), 
#           inherit.aes = F, alpha = 0.3)+
#   theme(legend.position = "bottom")
# p3
# YearLablePlot(p3)
# ggsave("Plots/SJRBVbyYear.tiff", device = "tiff", width = 7, height =5)
# 
# 
# summary(micro$CVPSWP)
# p3 = ggplot(filter(newdata, CVPSWP == 5100), aes(x = TempNotScaled, y = BV)) + 
#   geom_smooth(aes(color = as.factor(SJR)), se = FALSE)+
#   scale_color_viridis_d(name = "San Joaquin Flow", option = "C")+
#   ylab("Biovolume (um3/L)")+
#   facet_wrap(~Year)+
#   xlab("Water Temperature")+
#   #geom_text(data = yrtypes, aes(x = 6000, y = 180000, label = yeartype))+
#   theme_bw()
# #geom_rect(data = DFyear, aes(ymin = 0, ymax = 200000, xmin = MinExp, xmax =  MaxExp), 
# #          inherit.aes = F, alpha = 0.3)
# p3
# 
# ##############################################
# #maybe hold X2 or SJR constant and plot the rest
# 
# 
# p3 = ggplot(filter(newdata, SJR == 500, SACm == 30000), aes(x = TempNotScaled, y = BV)) + 
#   geom_smooth(aes(color = as.factor(CVPSWP)), se = FALSE)+
#   scale_color_viridis_d(name = "Exports", option = "C")+
#   ylab("Biovolume (um3/L)")+
#   facet_wrap(~Year)+
#   xlab("Water Temperature")+
#   #geom_text(data = yrtypes, aes(x = 6000, y = 180000, label = yeartype))+
#   theme_bw()
# #geom_rect(data = DFyear, aes(ymin = 0, ymax = 200000, xmin = MinExp, xmax =  MaxExp), 
# #          inherit.aes = F, alpha = 0.3)
# p3
# 
# 
# 
# 
# #####################################################
# #partial residual plots
# 
# library(visreg)
# 
# visreg(bx1a)
# #make all the partial residual plots
# p2 = visreg(bx1a, gg = T)
# 
# pSJR = p2[[2]]
# ptemp = p2[[3]]
# pexp = p2[[1]]
# preds = ptemp$data
# preds = bind_cols(micro, preds)
# 
# 
# preds2 = pSJR$data
# preds2 = bind_cols(micro, preds2)
# 
# preds3 = pexp$data
# preds3 = bind_cols(micro, preds3)
# 
# r.squaredGLMM(bx1a)
# 
# predsall = bind_rows(mutate(preds, parameter = "Temperature", Value = field.Water.temp),
#                      mutate(preds2, parameter = "San Joaquin Flow", Value = SJR/1000),
#                      mutate(preds3, parameter = "Project Exports", Value = CVPSWP/1000))
# 
# ylabM =  expression(paste("log-transformed ",italic(Microcystis), "biovolume rediduals"))
# ggplot(predsall, aes(x =Value, y = y))+
#   geom_point(aes(color = Year,  shape = Year))+
#   geom_smooth(method = "lm")+
#   theme_bw()+
#   scale_color_viridis_d(option = "B")+
#   facet_wrap(~parameter, scales = "free_x")+
#   ylab(ylabM)+
#   xlab("Thousand cfs                                Thousand cfs                               degrees C")+
#   theme(legend.position = "bottom")
# 
# ggsave("plots/PartialResid.tiff", device = "tiff", width =8, height =6)
# 
# ###########################################
# #export the model
# 
# fooSJR = summary(bx1a2)
# write.csv(fooSJR$coefficients, "outputs/ExpSJRSACm.csv")
# write.csv(fooSJR$varcor, "outputs/ExpSJRSACm_random.csv")
# 
# 
# ###################
# 
# #is biovolume greater int he South Delta?
# 
# ggplot(micro, aes(x = Station, y = logBV)) + geom_boxplot()
# 
# ggplot(microcystis, aes(x = Lon, y = Lat, size = 
# ggplot(microcystis, aes(x = Lon, y = Lat, size = log(dwr.MIC_totbvL+1))) + geom_point(position = "jitter")
# )) + geom_point(position = "jitter")
# 
# microcystis = mutate(microcystis, Region = case_when(Station %in% c("NA D26", "SJ NA", "NA MD10A", "MI NA", "RR P8", "RR NA", "VC NA", "OR")~ "South",
#                                                      TRUE ~ "North"))
# 
# ggplot(microcystis, aes(x = Region, y = log(dwr.MIC_totbvL+1))) + geom_boxplot()+
#   facet_wrap(~Year)
# 
# ggplot(microcystis, aes(x = Region, y = log(dwr.MIC_totbvL+1))) + geom_boxplot()
# 
# 
# ggplot(microcystis, aes(x = Station, fill = Region, y = log(dwr.MIC_totbvL+1))) + geom_boxplot()
# 
# 
# 
# lm1 = lmer(log(dwr.MIC_totbvL+1) ~ Region  + (1|Year), data = microcystis)
# summary(lm1)


############################################
#take year out as a random effect


mglobal4xyr = lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+ Outs+ X2s+SJRms+SACms+(1|Month)+(1|Station) ),
                  data = micro, na.action = "na.fail", REML = FALSE)
getAllTerms(mglobal4xyr)


modsSubx1yr <- dredge(mglobal4xyr, m.max =3,  rank = "BIC")

modsSubx1yr <- dredge(mglobal4xyr, subset = smat,  rank = "BIC")
modsSubx1yrall <- dredge(mglobal4xyr,  rank = "BIC")
write.csv(as.data.frame(modsSubx1yr), "outputs/MicroModelsJan2024.csv", row.names = F)
write.csv(as.data.frame(modsSubx1yrall), "outputs/MicroModelsJan2024all.csv", row.names = F)

#And this is yet another answer
bx1ayr =get.models(modsSubx1yr, subset = 1)[[1]]
bx1ayr = update(bx1ayr, REML = T)
summary(bx1ayr)
plot(allEffects(bx1ayr))

bx1ayr2 =get.models(modsSubx1yr, subset = 2)[[1]]
bx1ayr2 = update(bx1ayr2, REML = T)
summary(bx1ayr2)
plot(allEffects(bx1ayr2))


mod = lmer(logBV ~ Temp +SJRs+SACms+(1|Month)+(1|Station),
                   data = micro, na.action = "na.fail", REML = FALSE)
summary(mod)

#better plots of the effects
Tempscaled = data.frame(Temp = (c(15,17,19,21,23,25,27)-mean(micro$field.Water.temp))/sd(micro$field.Water.temp),
                        TempNotScaled = c(15,17,19,21,23,25,27))

#CS = data.frame(CVPSWPs= (seq(500, 12000, length.out=11)-mean(micro$CVPSWP))/sd(micro$CVPSWP),
#                CVPSWP = seq(500, 12000, length.out=11))

SJ =data.frame(SJRs =  (c(150,200,250,300,400,500,1000,2000,3000,5000)-mean(micro$SJR, na.rm =T))/sd(micro$SJR, na.rm =T),
               SJR = c(150,200,250,300,400,500,1000,2000,3000,5000))

Month = data.frame(Month = 7:12)

 SAC = group_by(micro, Year) %>%
   summarize(SACms = mean(SACms),
             SACm = round(mean(SACm)))
   
unique(micro$Year)
unique(micro$SACms)
newdata = data.frame(Month =7:12,
                     Station = "OR")

test = merge(newdata, Tempscaled)
test2 = merge(test, SAC) %>%
  merge(SJ) 
test2$predictions = predict(bx1ayr, newdata = test2)

ylabM2 =  expression(paste(italic(Microcystis), " biovolume, thousand um3/L"))
newdata = mutate(test2, BV = exp(predictions), lagBV = lag(BV), percent = (BV+lagBV)/lagBV)

ggplot(filter(newdata, Year == 2016), aes(x = SJR, y = BV, color = as.factor(TempNotScaled))) + geom_line()+
  facet_wrap(~Month)

###########$$$$$$$$$$$$$$$$$ FIgure in paper!
DFyear = DF %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  filter(Month ==8, Year %in%c(2014:2019))%>%
  group_by(Year) %>%
  summarize(MinSJR = min(SJR, na.rm =T), MaxSJR = max(SJR, na.rm =T)) %>%
  mutate(Year = as.factor(Year))%>%
  mutate(faceter = case_when(Year == 2014 ~ "2014, Spring Flow @ 10,217 cfs",
                             Year == 2015 ~ "2015, Spring Flow @ 9,706 cfs",
                             Year == 2016 ~ "2016, Spring Flow @ 27,088 cfs",
                             Year == 2017 ~ "2017, Spring Flow @ 65,731 cfs",
                             Year == 2018 ~ "2018, Spring Flow @ 21,698 cfs",
                             Year == 2019 ~ "2019, Spring Flow @ 21,698 cfs",))



realvalues = filter(micro, Month == 8, Station == "FT D19") %>%
  mutate( TempNotScaled =  (2 * ceiling(field.Water.temp/2) -1))%>%
  mutate(BV = exp(logBV), SACm = round(SACm))%>%
  mutate(faceter = case_when(Year == 2014 ~ "2014, Spring Flow @ 10,217 cfs",
                             Year == 2015 ~ "2015, Spring Flow @ 9,706 cfs",
                             Year == 2016 ~ "2016, Spring Flow @ 27,088 cfs",
                             Year == 2017 ~ "2017, Spring Flow @ 65,731 cfs",
                             Year == 2018 ~ "2018, Spring Flow @ 21,698 cfs",
                             Year == 2019 ~ "2019, Spring Flow @ 21,698 cfs",))

newdataPlot = filter(newdata, Month ==8) %>%
  mutate(faceter = case_when(Year == 2014 ~ "2014, Spring Flow @ 10,217 cfs",
                              Year == 2015 ~ "2015, Spring Flow @ 9,706 cfs",
                              Year == 2016 ~ "2016, Spring Flow @ 27,088 cfs",
                              Year == 2017 ~ "2017, Spring Flow @ 65,731 cfs",
                             Year == 2018 ~ "2018, Spring Flow @ 21,698 cfs",
                             Year == 2019 ~ "2019, Spring Flow @ 21,698 cfs",))

px = ggplot(filter(newdataPlot, Month == 8), aes(x = SJR/1000, y = BV/1000, color = as.factor(TempNotScaled))) + 
  geom_line()+
  facet_wrap(~faceter)+
  scale_color_brewer(name = "Temperature", palette = "OrRd")+
  geom_rect(data = DFyear, aes(ymin = 0, ymax = 100, xmin = MinSJR/1000, xmax =  MaxSJR/1000), 
                      inherit.aes = F, alpha = 0.3)+
  theme_bw()+ylab(ylabM2)+xlab("San Joaquin Flow, thousand cfs")+
  geom_point(data = realvalues, aes(x = SJR/1000, y = BV/1000, color = as.factor(TempNotScaled)))+
  theme(legend.position = "bottom")

tiff('plots/SJRBVbyYear.tiff', units="in", width=8, height=8, res=300, compression = 'lzw')
YearLablePlot(px)

dev.off()

#other esperimental plots
ggplot(filter(newdata, Month == 8), aes(x = SACm, y = BV, color = as.factor(TempNotScaled))) + geom_line()+
  facet_wrap(~SJR)

ggplot(filter(newdata, TempNotScaled == 21), aes(x = SJR, y = BV, color = as.factor(Month))) + geom_line()+
  facet_wrap(~Year)



ggplot(microcystis, aes(x = Collection_Date, y = field.Water.temp))+ geom_point(aes(color = as.factor(Collection_Month)))


microcystis = mutate(microcystis, doy = yday(Collection_Date))

ggplot(microcystis, aes(x = doy, y = field.Water.temp, color = as.factor(Survey_Year)))+ geom_point()+
  geom_smooth()

#http://127.0.0.1:14563/graphics/plot_zoom_png?width=1200&height=900
############################################################################

ggplot(micro, aes(x = Station, y = logBV))+ geom_boxplot()

#make a new version of everythign


newdata2 = left_join(newdata, DFyear) %>%
  filter(Month ==8)

yrtypes = data.frame(Year = c(2014:2019), yeartype = c("Critical", "Critical", "Below Normal", "Wet", "Dry", "Wet"))
#add water year type to these.

ggplot(newdata2, aes(x = SJR, y = BV, color = as.factor(TempNotScaled))) + geom_line()+
  geom_rect(data = DFyear, aes(ymin = 0, ymax = 30000, xmin = MinSJR, xmax =  MaxSJR), 
            inherit.aes = F, alpha = 0.3)+
  facet_wrap(~Year+SACm)

plot(allEffects(bx1ayr))

library(visreg)

visreg(bx1ayr)
#make all the partial residual plots
p2 = visreg(bx1ayr, gg = T)

pSJR = p2[[2]]
ptemp = p2[[3]]
psacm = p2[[1]]
preds = ptemp$data
preds = bind_cols(micro, preds)


preds2 = pSJR$data
preds2 = bind_cols(micro, preds2)

preds3 = psacm$data
preds3 = bind_cols(micro, preds3)

r.squaredGLMM(bx1ayr)

predsall = bind_rows(mutate(preds, parameter = "Temperature", Value = field.Water.temp),
                     mutate(preds2, parameter = "San Joaquin Flow", Value = SJR/1000),
                     mutate(preds3, parameter = "Spring Sacramento Flow", Value = SACm/1000))

ylabM =  expression(paste("log-transformed ",italic(Microcystis), "biovolume rediduals"))

###########$$$$$$$$$$$$$$$$$ FIgure in paper!
ggplot(predsall, aes(x =Value, y = y))+
  geom_point(aes(color = Year,  shape = Year))+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_color_viridis_d(option = "B")+
  facet_wrap(~parameter, scales = "free_x")+
  ylab(ylabM)+
  xlab("Thousand cfs                                Thousand cfs                               degrees C")+
  theme(legend.position = "bottom")

ggsave("plots/PartialResid_noyear.tiff", device = "tiff", width =8, height =6)

###########################################
#export the model

fooSJR = summary(bx1ayr)
write.csv(fooSJR$coefficients, "outputs/SJRSACm_year.csv")
write.csv(fooSJR$varcor, "outputs/SJRSACm_random_year.csv")


###################

#is biovolume greater int he South Delta?
