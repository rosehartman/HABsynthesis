#get teh dayflow data

library(tidyverse)
library(lubridate)
library(smonitr)
#I need to read in the Dayflow data from the CNRA portal
# https://data.cnra.ca.gov/dataset/dayflow
#Still needs a little fiddling, but much better!.

Dayflow = get_odp_data(pkg_id = "dayflow", fnames = "Dayflow Results")



DF1997_2020 =  Dayflow$`Dayflow Results 1997 - 2023` %>%
  mutate( Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  select(Date, OUT, EXPORTS, SJR, GCD, SAC, CVP, SWP, X2, TOT)

#now I can put them all together!
DF = DF1997_2020
save(DF, file = "Dayflow1997_2023.RData")

#what's inflow in the summer?

summer = filter(DF, month(Date) %in% c(6,7,8,9))

mean(summer$TOT)*.02832


#Mean annual flow

meanflow = mutate(DF, Year = year(Date)) %>%
  group_by(Year) %>%
  summarise(across(c(OUT:TOT), mean, .names = "{.col}_mean"))

#mean spring flow
meanspringflow = mutate(DF, Year = year(Date), Month = month(Date)) %>%
  filter(Month %in% c(2,3,4,5)) %>%
  group_by(Year) %>%
  summarise(across(c(OUT:TOT), mean, .names = "{.col}_spr_mean"))

write.csv(meanflow, "outputs/meanflow")

write.csv(meanspringflow, "outputs/meanspringflow")

#SJR flow, exports, delta outflow, and sac flow

DFlong = DF %>%
  mutate(CVPSWP = CVP+SWP) %>%
  pivot_longer(cols = c(-Date), names_to = "Parameter", values_to = "Flow") %>%
  filter(!(Flow <0 & Parameter == "OUT"), Parameter %in% c("OUT", "CVPSWP", "SJR", "SAC", "X2"), year(Date) %in% c(2014:2019)) %>%
  mutate(Param = factor(Parameter, levels = c("CVPSWP", "SAC", "SJR", "OUT", "X2"),
         labels = c("CVP+SWP (m3/sec)", "Sac (m3/sec)", "SJR (m3/sec)","Outflow (m3/sec)",  "X2 (km)"))) %>%
  mutate(CMS = case_when(Parameter == "X2" ~ Flow,
                          TRUE ~ Flow*0.028316831998814504))


test = mutate(micro, DOY = yday(Collection_Date)) %>%
  group_by(Year) %>%
  summarize(startdate = min(Collection_Date), enddate = max(Collection_Date))

minmax = DFlong %>%
  filter()%>%
  group_by(Param) %>%
  summarize(Ymin = min(Flow), Ymax = max(Flow)) %>%
  merge(test)

fg = ggplot(DFlong,
       aes(x = Date, y = Flow)) +
  geom_line()+
  geom_rect(data = minmax, aes(xmin = startdate, xmax = enddate, ymin = Ymin, ymax = Ymax), alpha = 0.5, inherit.aes = FALSE)+
  facet_wrap(~Param, scales = "free", nrow = 5, strip.position = "left")+
  theme_bw()+
  theme(strip.placement = "outside")+
  ylab(NULL)+xlab(NULL)


##add water year index

yrs = read_csv("data/yearassignments.csv") %>%
  filter(Year %in% c(2014:2019)) %>%
  select(Year, Yr_type, Index) %>%
  mutate(Param = "Sac Valley Index",
         Yr_type = factor(Yr_type, levels = c("Critical", "Below Normal", "Wet"),
                          labels = c("Critical", "Below\nNormal","Wet"))) 

ytg = ggplot(yrs, aes(x = Year, y = Index))+
  geom_col(aes(fill = Yr_type), color = "black")+
  facet_wrap(~Param, scales = "free",strip.position = "left")+
  theme_bw()+
  scale_fill_manual(values =c( "tomato","orange", "lightblue"), guide = NULL)+
  geom_text(aes(label = Yr_type), nudge_y = -2.1, size =3.5)+
  theme(strip.placement = "outside")+
  ylab(NULL)+xlab(NULL)

library(patchwork)

fg/ytg+
  plot_layout(heights = unit(c(8, 2.2), c('in', 'null')))

ggsave("plots/flows.tiff", device = "tiff", width =8, height =10)


#summary statistics on summer flows

summer = filter(DFlong, month(Date) %in% c(6,7,8,9)) %>%
  group_by(Param)%>%
           summarize(min = min(Flow), max = max(Flow), mean = mean(Flow), median = median(Flow))
