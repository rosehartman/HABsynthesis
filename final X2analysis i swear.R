#let's start over and walk through this one step at a time

library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(MuMIn)
library(readxl)
library(RColorBrewer)
library(visreg)
library(sf)
library(corrplot)
library(ggcorrplot)
library(effects)
library(grid)

#load dataset
microcystis = read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.1076.1&entityid=86da4465dde063c89afb0dc256bfa619")


#THe dayflow dataset is a bit difficult to deal with when you first download it,
#so I have a version I got ealier and cleaned up
load("Dayflow1997_2021.RData")
str(DF)


yrs = read_csv("data/yearassignments.csv") %>%
  filter(Year %in% c(2014:2019)) %>%
  mutate(Year = as.factor(Year)) %>%
  select(Year, Yr_type, Index)

#calculate spring outflow
DFspring = DF %>%
  mutate(Month = month(Date), Year = year(Date)) %>%
  filter(Month %in% c(2,3,4,5)) %>%
  group_by(Year) %>%
  summarize(OUTm = mean(OUT),  SJRm = mean(SJR), SACm = mean(SAC)) %>%
  ungroup()%>%
  mutate(OUTms = scale(OUTm), 
         SJRms = scale(SJRm), SACms = scale(SACm), Year = as.factor(Year))

#bind dayflow data to microcystis data
microX = left_join(microcystis, DF, by = c("Collection_Date" = "Date")) %>%
  ungroup() %>%
  mutate(CVPSWP = CVP + SWP, logBV = log(dwr.MIC_totbvL/100000+1), 
         Year = as.factor(year(Collection_Date)), 
         Month = month(Collection_Date),
         Temp = scale(field.Water.temp), 
         X2s = scale(X2), 
         SACs = scale(SAC),
         SJRs = scale(SJR),
         Outs = scale(OUT),
         CVPSWPs = scale(CVPSWP), 
         DWR_Site = str_remove(DWR_Site, "rep"),
         Station = paste(DWR_Site, EMP_Site),
         Station = case_when(Station %in% c("OR D28A", "OR NA") ~ "OR",
                             TRUE ~ Station)) %>%
  filter(!is.na(field.Water.temp), Month !=6)

#add spring flow and slect important variables
micro = left_join(microX, DFspring) %>%
  select(logBV, Temp, CVPSWPs, SACs, SJRs, OUTms, OUT,Outs, SJRms, SACms,  Year, 
         Month, Collection_Date, X2s, CVPSWPs, Station, SAC, SJR, field.Water.temp, CVPSWP, SACm) %>%
  left_join(yrs)

#calculate correlations

#let's calculate the correlation coefficient for all the variables, maybe that is better than vif
is.correlated <- function(i, j, data, conf.level = 0.95, cutoff = 0.5, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level,...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}

# Need vectorized function to use with 'outer'
vCorrelated <- Vectorize(is.correlated, c("i", "j"))

#####################################################################

#funciton to color the plots based on water year
YearLablePlot = function(plot){
  g <- ggplot_gtable(ggplot_build(plot))
  stripr <- which(grepl('strip-t', g$layout$name))
  fills <- c("skyblue","tan1","steelblue2","red2", "indianred","gold")
  k <- 1
  
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  grid.draw(g)
}

#########################################################################
#create the correlation matrix

microtest = as.data.frame(micro[,c(2:5,8, 14, 22)])
smat <- outer(1:7, 1:7, vCorrelated, data = microtest, cutoff = 0.71)
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
  mutate(xName = factor(xName, levels = c("CVPSWPs", "SACs", "SJRs", "Outs",  "X2s", "Index"),
                        labels = c("CVP+SWP", "Sac", "SJR",  "Outflow",  "X2", "Index")))

#Plot the correlation matrix
ggplot(corp2, aes(x = x, y = y))+ geom_tile(aes(fill = corr), color = "black")+
  geom_tile(data = filter(corp2, included ==2), fill = "grey", color = "black")+
  scale_alpha(range = c(1,0))+
  geom_text(aes(label = round(corr, 2)))+
  scale_fill_gradient2(name = "Correlation\nCoefficient")+
  scale_x_continuous(breaks = c(2:7), labels = c(levels(corp2$xName)))+
  scale_y_continuous(breaks = c(2:7), labels = rev(c("Temperature", "CVP+SWP", "Sac", "SJR", 
                                                      "Outflow", "X2")))+
  ylab(NULL)+xlab(NULL)+ theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust =1))

ggsave("plots/corrmatrix.tiff", device = "tiff", width =6, height =5)

###########################################################################
#let's do this - with year as random effect

# 
# mglobal= lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+ Outs+ X2s+SJRms+SACms+(1|Month)+(1|Station) +(1|Year),
# data = micro, na.action = "na.fail", REML = FALSE)
# 
# 
# mods <- dredge(mglobal, subset = smat,  rank = "BIC")
# bestmod =get.models(mods, subset = 1)[[1]]
# bestmod = update(bestmod, REML = T)
# summary(bestmod)
# plot(allEffects(bestmod))
# 
# 
# #no without year as a random effect
# 
# mglobal2= lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+OUTms+ Outs+ X2s+SJRms+SACms+(1|Month)+(1|Station),
#               data = micro, na.action = "na.fail", REML = FALSE)
# 
# 
# mods2 <- dredge(mglobal2, subset = smat,  rank = "BIC")
# #Yeah, this is saying the best model is one without any anual factors
# #But I don't know if I can have a model without some factor that keeps year in there and still be OK. 
# 
# #what if i force it to have spring outflow in there?
# 
# mods2.1 = dredge(mglobal2, fixed = "OUTms", subset = smat,  rank = "BIC")
# mods2.2 = dredge(mglobal2, fixed = "SACms", subset = smat,  rank = "BIC")
# mods2.3 = dredge(mglobal2, fixed = "SJRms", subset = smat,  rank = "BIC")
#mods2.4 = bind_rows(mods2.1, mods2.2, mods2.3)

#now have water year index in there in all models



mglobal3= lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+ Outs+ X2s+ Index+(1|Month)+(1|Station),
               data = micro, na.action = "na.fail", REML = FALSE)


mods3 <- dredge(mglobal3, subset = smat, fixed = "Index",  rank = "BIC")
write.csv(as.data.frame(mods3), file = "outputs/MicoModelswIndex.csv")

mod3best = get.models(mods3, subset = 1)[[1]]
mod3best = update(mod3best, REML = TRUE)
summary(mod3best)
plot(allEffects(mod3best))
foo = summary(mod3best)
write.csv(foo$coefficients, "outputs/mod3best.csv")

#compare water year index to spring sac flow
wyi= lmer(logBV ~ Temp +SJRs+ Index+(1|Month)+(1|Station),
               data = micro, na.action = "na.fail", REML = FALSE)
springsac= lmer(logBV ~  Temp +SJRs+ SACms+(1|Month)+(1|Station),
               data = micro, na.action = "na.fail", REML = FALSE)
BIC(wyi)
BIC(springsac)
# 
# mglobal3.1= lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+ Outs+ X2s+ SACms+(1|Month)+(1|Station),
#                data = micro, na.action = "na.fail", REML = FALSE)
# 
# 
# mods3.1 <- dredge(mglobal3.1, subset = smat, fixed = "SACms",  rank = "BIC")
# mod3best.1 = get.models(mods3.1, subset = 1)[[1]]
# summary(mod3best.1)
# plot(allEffects(mod3best.1))

############################################
#partial resisduals plot of model with the water year index in there.
p3 = visreg(mod3best, gg = T)

ptemp = p3[[2]]$data
pindex = p3[[3]]$data
psjr = p3[[1]]$data

preds = bind_cols(micro, ptemp)
preds2 = bind_cols(micro, psjr)
preds3 = bind_cols(micro, pindex)

r.squaredGLMM(mod3best)

predsall = bind_rows(mutate(preds, parameter = "Temperature", Value = field.Water.temp),
                     mutate(preds2, parameter = "San Joaquin Flow", Value = SJR/1000),
                     mutate(preds3, parameter = "Sac Valley Index", Value = Index))

ylabM =  expression(paste("log-transformed ",italic(Microcystis), "biovolume rediduals"))

pyear = predictorEffect("Index", mod3best, focal.levels = c(4.07, 4, 6.71, 14.1, 7.14, 10.3), residuals=TRUE)
pyear2 = data.frame(Year = yrs$Year, Fit = pyear$fit, SE = pyear$se)
pyearpoints = preds3 %>%
  left_join(pyear2) %>%
  mutate(Parameter = "Year")

ggplot(pyearpoints )+ geom_point(aes(x = Year, y = y, color = Yr_type), position = "jitter")+
  geom_point(aes(x = Year, y = Fit))+
  facet_wrap(~Parameter)+
  scale_color_manual(values = c("darkred","darkorange", "skyblue"))+
  geom_errorbar(aes(x = Year, ymin = Fit-SE, ymax = Fit+SE))+
  ylab(ylabM)+
  theme_bw()
ggsave("plots/PartialResid_year.tiff", device = "tiff", width =5, height =5)


###########

ggplot(predsall, aes(x =Value, y = y))+
  geom_point(aes(color = Year,  shape = Year))+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_color_viridis_d(option = "B")+
  facet_wrap(~parameter, scales = "free_x")+
  ylab(ylabM)+
  xlab("                                Thousand cfs                               degrees C")+
  theme(legend.position = "bottom")
ggsave("plots/PartialResid_exports.tiff", device = "tiff", width =8, height =6)

#############################################
#THIS IS FIGURE figure 5
ggplot(predsall, aes(x =Value, y = y))+
  geom_point(aes(color = Year,  shape =  Yr_type))+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_color_manual(values = c("darkred","indianred","tan1", "skyblue", "gold", "steelblue"))+
  facet_wrap(~parameter, scales = "free_x")+
  ylab(ylabM)+
  xlab("                                Thousand cfs                               degrees C")+
  theme(legend.position = "bottom")
ggsave("plots/PartialResid_exports_wy.tiff", device = "tiff", width =8, height =6)

##############################################################
#what about year as a fixed effect?

# 
# mglobal3= lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+ Outs+ X2s+Year +(1|Month)+(1|Station) ,
#               data = micro, na.action = "na.fail", REML = FALSE)
# 
# mods3 <- dredge(mglobal3,  rank = "BIC")
# #This basically has the same result as having year as a random effect. 
# 
# bestmod2 = lmer(logBV ~ Temp +CVPSWPs+SJRs+Year+(1|Month)+(1|Station) ,
#                 data = micro, na.action = "na.fail", REML = TRUE)
# summary(bestmod2)
# plot(allEffects(bestmod2))


##############################################################################
#effectsplot

#better plots of the effects
Tempscaled = data.frame(Temp = (c(15,16, 17,18,19,21,23,25,27)-mean(micro$field.Water.temp))/sd(micro$field.Water.temp),
                        TempNotScaled = c(15,16, 17,18,19,21,23,25,27))

Index = select(yrs, Year, Yr_type, Index) %>%
  mutate(Year = as.factor(Year))

SJ =data.frame(SJRs =  (c(150,200,250,300,400,500,1000,2000,3000,5000)-mean(micro$SJR, na.rm =T))/sd(micro$SJR, na.rm =T),
               SJR = c(150,200,250,300,400,500,1000,2000,3000,5000))


newdata = data.frame(Month =7:12,
                     Station = "FT D19")

test = merge(newdata, Tempscaled)
test2 = merge(test, Index) %>%
  merge(SJ) 
  
test2$predictions = predict(mod3best, newdata = test2)


ylabM2 =  expression(paste(italic(Microcystis), " biovolume, thousand um3/L"))
newdata = mutate(test2, BV = exp(predictions), lagBV = lag(BV), percent = (BV+lagBV)/lagBV)

ggplot(filter(newdata, Year == 2016), aes(x = SJR, y = BV, color = as.factor(TempNotScaled))) + geom_line()+
  facet_wrap(~Month)

###########
#Plot of predictions
DFyear = DF %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  filter(Month ==8, Year %in%c(2014:2019))%>%
  group_by(Year) %>%
  summarize(MinSJR = min(SJR, na.rm =T), MaxSJR = max(SJR, na.rm =T)) %>%
  mutate(Year = as.factor(Year))%>%
  left_join(yrs) %>%
  mutate(faceter = paste(Year, "Index:", Index))


realvalues = filter(micro, Month == 8, Station == "FT D19") %>%
  mutate( TempNotScaled =  (2 * ceiling(field.Water.temp/2) -1))%>%
  mutate(BV = exp(logBV))%>%
  mutate(faceter = paste(Year, "Index:", Index))

newdataPlot = filter(newdata, Month ==8) %>%
  mutate(faceter = paste(Year, "Index:", Index))

###########################################
#THIS IS FIGURE 6
px = ggplot(newdataPlot, aes(x = SJR/1000, y = BV/1000, color = as.factor(TempNotScaled))) + 
  geom_line()+
  facet_wrap(~faceter)+
  scale_color_brewer(name = "Temperature", palette = "OrRd")+
  scale_fill_manual(name = "Temperature", values = c("#D7301F", "#B30000" ), guide = NULL)+
  geom_rect(data = DFyear, aes(ymin = 0, ymax = 300, xmin = MinSJR/1000, xmax =  MaxSJR/1000), 
            inherit.aes = F, alpha = 0.3)+
  theme_bw()+ylab(ylabM2)+xlab("San Joaquin Flow, thousand cfs")+
  geom_point(data = realvalues, aes(x = SJR/1000, y = BV/1000, fill = as.factor(TempNotScaled)),
             shape = 21, size =3, color = "black")+
  theme(legend.position = "bottom")
px

tiff('plots/SJRBVbyYear.tiff', units="in", width=8, height=7, res=300, compression = 'lzw')
YearLablePlot(px)

dev.off()



ggplot(microX, aes(x = SJR*0.028316831998814504, y = logBV)) + geom_point()
