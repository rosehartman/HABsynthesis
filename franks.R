#franks tract

library(tidyverse)
library(cder)
library(lubridate)

FRK = cdec_query("FRK", sensors = c(28, 100, 25, 61, 62, 27), duration = "E",
                 start.date = ymd("2022-06-01"), end.date = ymd("2022-07-31"))
FRK = filter(FRK, Value < 100 & SensorNumber == 28 | SensorNumber != 28) %>%
  filter(Value < 100 & SensorNumber == 27 | SensorNumber != 27)

ggplot(FRK, aes(x = DateTime, y = Value))+
  geom_point()+facet_wrap(~SensorType, scales = "free_y")
