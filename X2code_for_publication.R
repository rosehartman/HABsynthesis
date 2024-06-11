#clean data and code for publication
# Rosemary Hartman, Department of Water Resources
# 2024 Feb 21

#load necessary packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(MuMIn)
library(readxl)
library(effects)
library(corrplot)
library(ggcorrplot)
library(visreg)
library(grid)
library(gridExtra)
library(patchwork)

##################################################################################

# Data Access and integration

#Data on Microcystis biovolume came from Kurobe et al. 
# Kurobe, T., P.W. Lehman, S. Lesmeister, and S.J. Teh. 2022. 
# Cyanobacteria abundance, cyanotoxin concentration, and 
# water quality data for the upper San Francisco Estuary, 
# California, USA: 2014-2019 ver 1. Environmental Data Initiative. 
# https://doi.org/10.6073/pasta/b4e180cb925266171285375a3e202635 (Accessed 2023-09-06). 

#read in the dataset from the EDI publication
microcystis = read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.1076.1&entityid=86da4465dde063c89afb0dc256bfa619")


#Dayflow
#Data on flow parameters came from the Department of Water Resources Dayflow model.
# https://data.cnra.ca.gov/dataset/dayflow

#Import the dayflow dataset and subset it to times and parameters of interest


Dayflow2021 = read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/83122ce7-e7f5-4ad1-b5e9-6a7032cb1117/download/dayflowcalculations2021.csv")

#the older data uses a different formula for some of the parameters, so they have different names. read the documentation for details on the differences.
Dayflow29_39 = read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/ab12e85f-82f4-4723-9973-deeed41b2057/download/dayflow-results-1929-1939.csv") %>%
  rename(YOLO = YOLO1, OUT = OUT1, RIO = RIO1, DIVER = DIVER1, EFFDIV = EFFDIV1, WEST = WEST1, Mo = Month) %>%
  mutate(Date = mdy(Date), DIVER = as.numeric(DIVER), EFFDIV = as.numeric(EFFDIV))

Dayflow40_49 = read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/bf58c67c-63b4-47d4-9a25-2b95e5479a0c/download/dayflow-results-1940-1949.csv")%>%
  rename(YOLO = YOLO2, OUT = OUT2, RIO = RIO2, DIVER = DIVER2, EFFDIV = EFFDIV2, WEST = WEST2, TOT = TOT2)%>%
  mutate(Date = mdy(Date))

Dayflow50_55 = read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/9225dbe7-54a6-4466-b360-e66f51407683/download/dayflow-results-1950-1955.csv")%>%
  rename(YOLO = YOLO2, OUT = OUT2, RIO = RIO2, DIVER = DIVER2, EFFDIV = EFFDIV2, WEST = WEST2, TOT = TOT2)%>%
  mutate(Date = mdy(Date))

Dayflow56_69 = read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/3109f3ef-b77b-4288-9ece-3483899d10da/download/dayflow-results-1956-1969.csv")%>%
  mutate(Date = mdy(Date))

Dayflow70_83 = read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/a0a46a1d-bec5-4db9-b331-655e306860ba/download/dayflow-results-1970-1983.csv")%>%
  mutate(Date = mdy(Date))

Dayflow84_96 = read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/cb04e626-9729-4105-af81-f6e5a37f116a/download/dayflow-results-1984-1996.csv")%>%
  mutate(Date = mdy(Date))

Dayflow97_20 = read_csv("https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/21c377fe-53b8-4bd6-9e1f-2025221be095/download/dayflow-results-1997-2020.csv")%>%
  mutate(Date = mdy(Date))


DF = bind_rows(Dayflow29_39, Dayflow40_49, Dayflow50_55, Dayflow56_69, Dayflow70_83, Dayflow84_96, Dayflow97_20,
                    Dayflow2021)



#X2

#historic X2 from Hutton et al. 2016
# https://ascelibrary.org/doi/10.1061/%28ASCE%29WR.1943-5452.0000617
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

#Residence Time

#Residence time came from Hartman at al. 2024
# Hartman, R., E. Stumpner, C. Burdi, D. Bosworth, A. Maguire, and IEP 
# Drought Synthesis Team. 2024. Dry me a river: Ecological effects 
# of drought in the upper San Francisco Estuary. San Francisco 
# Estuary and Watershed Science.  in press. 

#Read in the data. Month2 = numeric month, Year = calendar year, WY = water year, 
#SACRT = modeled scaramento river residence time. SJRT = modeled San Joaquin river residence time.

DFRTall = read_csv("data/ResidenceTime.csv")

RTs = left_join(DFRTall, X2sall, by = c("Year", "Month2"="Month"))


#Year types and water year indices are from DWR"s hisorical water year data
# https://cdec.water.ca.gov/reportapp/javareports?name=WSIHIST 

yrs = read_csv("data/yearassignments.csv") %>%
  filter(Year %in% c(2014:2019)) %>%
  mutate(Year = as.factor(Year)) %>%
  select(Year, Yr_type, Index)


####################################################################################
#X2=residence time analysis

#This is Figure 3


SacRTplot = ggplot(RTs, aes(x = X2, y = log(SACRT), color = Month2))+ geom_point()+
  scale_color_viridis_c(name = "Month\nof year", guide = NULL)+
  theme_bw()+
  annotate("text", x = 40, y = 5, label = "A", fontface = "bold", size =6)+
  geom_smooth(method = "lm")+ ylab("log-transformed Sacramento \nResidence Time (Days)")
#OK, that's pretty good

SJRRTplot = ggplot(RTs, aes(x = X2, y = log(SJRT)))+ geom_point(aes(color = Month2))+
  scale_color_viridis_c(name = "Month\nof year")+
  geom_smooth(method = "lm")+ ylab("log-transformed San Joaquin\n Residence Time (Days)")+
  theme_bw()+
  annotate("text", x = 40, y = 5.5, label = "B", fontface = "bold", size =6)
#Much less good. 

SacRTplot + SJRRTplot+ plot_layout(widths = c(3.8,4))
  

ggsave("plots/bothrestime.tiff", device = "tiff", width = 9, height = 5)

#This is Figure 4

ggsave("plots/SJRTvX2.tiff", device = "tiff", width = 5, height = 4)

#here is the model of San Joaquin residence time versus X2

lmsj = lmer(log(SJRT) ~X2 + (1|Year)+ (1|Month), data = RTs)
summary(lmsj)

#the R-squared tells you how much of the variance is explained by X2
r.squaredGLMM(lmsj)


#here is the model of Sacramento residence time versus X2
lmsac = lmer(log(SACRT) ~X2  + (1|Year)+ (1|Month), data = RTs)
summary(lmsac)

#the R-squared tells you how much of the variance is explained by X2
r.squaredGLMM(lmsac)

#output of model results
sj = as.data.frame(summary(lmsj)$coefficients) %>%
  mutate(River = "San Joaquin")
sac = as.data.frame(summary(lmsac)$coefficients) %>%
  mutate(River = "Sacramento")

lms = bind_rows(sj, sac)

#This is Table 1
write.csv(lms, "outputs/residenceTimeModels.csv")

##################################################################################
#Microcystis versus temperature and flow parameters

# Step 1. Combine the data

#bind dayflow data to microcystis data
#convert CFS to CMS
micro = left_join(microcystis, DF, by = c("Collection_Date" = "Date")) %>%
  ungroup() %>%
  mutate(CVPSWP = CVP + SWP, logBV = log(dwr.MIC_totbvL/100000+1), 
         Year = as.factor(year(Collection_Date)), 
         Month = month(Collection_Date),
         Temp = scale(field.Water.temp), 
         X2s = scale(X2), 
         SACcms = SAC/35.315,
         SACs = scale(SACcms),
         SJRcms = SJR/35.315,
         SJRs = scale(SJRcms),
         OUTcms = OUT/35.315,
         Outs = scale(OUTcms),
         CVPSWPcms = CVPSWP/35.315,
         CVPSWPs = scale(CVPSWPcms), 
         DWR_Site = str_remove(DWR_Site, "rep"),
         Station = paste(DWR_Site, EMP_Site),
         Station = case_when(Station %in% c("OR D28A", "OR NA") ~ "OR",
                             TRUE ~ Station)) %>%
  filter(!is.na(field.Water.temp), Month !=6) %>%
  select(logBV, Temp, CVPSWPs, SACs, SJRs,  OUTcms, Outs,  Year, 
         Month, Collection_Date, X2s, CVPSWPs, Station, SACcms, SJRcms, 
         field.Water.temp, CVPSWPcms) %>%
  left_join(yrs)


# Step 2. Calculate correlations

#This function takes a data frame and calculates the correlatoin coefficient for all of them
#and produces as matrix to input to the model selection analysis
is.correlated <- function(i, j, data, conf.level = 0.95, cutoff = 0.5, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level,...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}

# Need vectorized function to use with 'outer'
vCorrelated <- Vectorize(is.correlated, c("i", "j"))


microtest = as.data.frame(micro[,c(2:5,7, 11, 18)])
smat <- outer(1:7, 1:7, vCorrelated, data = microtest, cutoff = 0.71)
nm <- colnames(microtest)
dimnames(smat) <- list(nm, nm)

#correlation coefficient matrix

Corrs = cor(microtest)
corp = corrplot.mixed(Corrs, upper = "ellipse")

#make a matrix of whethe ror not things are included
smatpvals = apply(smat,c(1,2), FUN = function(x){ if(is.na(x)) 1 else if(x) 0 else 2})

#Reorganize the data and fix the labels.
corp2 = corp$corrPos
corpPs = mutate(as.data.frame(smatpvals), xName = row.names(smatpvals)) %>%
  pivot_longer(cols = c(Temp:X2s), names_to = "yName", values_to = "included") 
corp2 = left_join(corp2, corpPs) %>%
  filter(included !=1) %>%
  mutate(xName = factor(xName, levels = c("CVPSWPs", "SACs", "SJRs", "Outs",  "X2s", "Index"),
                        labels = c("CVP+SWP", "Sac", "SJR",  "Outflow",  "X2", "Index")))

#Plot the correlation matrix.
#This is Figure 2
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


# Step 3. Set up global model and run all the possible models that don't have highly correlated variables.

#Global model. all variables have been scaled and centered
# Random effects of month and station included
# Temp = Water Temperature 
# CVPSWPs = Central valley project and state water project exports
# SACs = Sacramento Valley flow 
# SJRs = San Joaquin River Flow
# Outs = Net Delta Outflow
# Index = Sacramento Valley Index

mglobal3= lmer(logBV ~ Temp +CVPSWPs+ SACs+SJRs+ Outs+ X2s+ Index+(1|Month)+(1|Station),
               data = micro, na.action = "na.fail", REML = FALSE)


#Run through all possible models. The 'smat' is the matrix of which variables should be inclued in teh same model.
#Include the water year index in all models.
smat[7,] = TRUE
mods3 <- dredge(mglobal3, subset = smat, fixed = "Index",  rank = "BIC")

#This is table 2
write.csv(as.data.frame(mods3), file = "outputs/MicoModelswIndex.csv")

#Extract the top models
mod3best = get.models(mods3, subset = 1)[[1]]

#Update the model to use restricted maximum likelihood estimator
mod3best = update(mod3best, REML = TRUE)

#check out the results of the model
summary(mod3best)
plot(allEffects(mod3best))

#Print the coeffeicients. This is Table 3
foo = summary(mod3best)
write.csv(foo$coefficients, "outputs/mod3best.csv")



# Step 4. Plot the results of the final model

#Here are the parital residuals plots for each of teh predictor variables. 


#Use the 'visreg' package to calculate the partial residuals
p3 = visreg(mod3best, gg = T)

#seperate the data for each of the predictors
ptemp = p3[[2]]$data
pindex = p3[[3]]$data
psjr = p3[[1]]$data

#bind the datasets togeter. 
preds = bind_cols(micro, ptemp)
preds2 = bind_cols(micro, psjr)
preds3 = bind_cols(micro, pindex)

predsall = bind_rows(mutate(preds, parameter = "Temperature", Value = field.Water.temp),
                     mutate(preds2, parameter = "San Joaquin Flow", Value = SJRcms/100),
                     mutate(preds3, parameter = "Sac Valley Index", Value = Index))


#This is Figure 5.
ylabM =  expression(paste("log-transformed ",italic(Microcystis), "biovolume rediduals"))

ggplot(predsall, aes(x =Value, y = y))+
  geom_point(aes(color = Year,  shape =  Yr_type))+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_shape(name = "Year Type")+
  scale_color_manual(values = c("darkred","indianred","tan1", "skyblue", "gold", "steelblue"))+
  facet_wrap(~parameter, scales = "free_x")+
  ylab(ylabM)+
  xlab("                                      m3/sec                                   degrees C")+
  theme(legend.position = "bottom")
ggsave("plots/PartialResid_exports_wy.tiff", device = "tiff", width =8, height =6)


#Now we will plot the range of microcystis we would predict for various levels of temperature and flow


#function to color the plots based on water year
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

#Set up a dataset with a range of temperatures and flows. We need the scaled and centered
#version to use with 'predict' as well as the non-scaled version to plot.

#start with temperature
Tempscaled = data.frame(Temp = (c(15,16, 17,18,19,21,23,25,27)-mean(micro$field.Water.temp))/sd(micro$field.Water.temp),
                        TempNotScaled = c(15,16, 17,18,19,21,23,25,27))
#one index for each year
Index = select(yrs, Year, Yr_type, Index) %>%
  mutate(Year = as.factor(Year))

#Range of San Joaquin flows
SJ =data.frame(SJRs =  (c(4.3,5.7, 7.1,8.5,11.5,14,28.5,57,86,143)-mean(micro$SJRcms, na.rm =T))/sd(micro$SJRcms, na.rm =T),
               SJRcms = c(4.3,5.7, 7.1,8.5,11.5,14,28.5,57,86,143))

#All the months. And choose Franks Tract as an example station.
newdata = data.frame(Month =7:12,
                     Station = "FT D19")

#put all the data together
test = merge(newdata, Tempscaled)
test2 = merge(test, Index) %>%
  merge(SJ) 

#add a column with predictions from the model
test2$predictions = predict(mod3best, newdata = test2)

#Set things up to plot
ylabM2 =  expression(paste(italic(Microcystis), " biovolume, thousand um3/L"))
newdata = mutate(test2, BV = exp(predictions), lagBV = lag(BV), percent = (BV+lagBV)/lagBV)

#calculate minimum and maximum flow parameters seen in each year
DFyear = DF %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  filter(Month ==8, Year %in%c(2014:2019))%>%
  group_by(Year) %>%
  summarize(MinSJR = min(SJR, na.rm =T)/35, MaxSJR = max(SJR, na.rm =T)/35) %>%
  mutate(Year = as.factor(Year))%>%
  left_join(yrs) %>%
  mutate(faceter = paste(Year, "Index:", Index))

#subset the actual microcystis values so we can put them on the plot
realvalues = filter(micro, Month == 8, Station == "FT D19") %>%
  mutate( TempNotScaled =  (2 * ceiling(field.Water.temp/2) -1))%>%
  mutate(BV = exp(logBV))%>%
  mutate(faceter = paste(Year, "Index:", Index))

#We'll just plot august, and put a label on each year/index
newdataPlot = filter(newdata, Month ==8) %>%
  mutate(faceter = paste(Year, "Index:", Index))

#THIS IS FIGURE 6
px = ggplot(newdataPlot, aes(x = SJRcms, y = BV/1000, color = as.factor(TempNotScaled))) + 
  geom_line()+
  facet_wrap(~faceter)+
  scale_color_brewer(name = "Temperature", palette = "OrRd")+
  scale_fill_manual(name = "Temperature", values = c("#D7301F", "#B30000" ), guide = NULL)+
  geom_rect(data = DFyear, aes(ymin = 0, ymax = 300, xmin = MinSJR, xmax =  MaxSJR), 
            inherit.aes = F, alpha = 0.3)+
  theme_bw()+ylab(ylabM2)+xlab("San Joaquin Flow, m3/sec")+
  geom_point(data = realvalues, aes(x = SJRcms, y = BV/1000, fill = as.factor(TempNotScaled)),
             shape = 21, size =3, color = "black")+
  theme(legend.position = "bottom")
px

tiff('plots/SJRBVbyYear.tiff', units="in", width=8, height=7, res=300, compression = 'lzw')
YearLablePlot(px)

dev.off()


