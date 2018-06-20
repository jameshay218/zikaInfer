## Using Cuevas et al. 2016 data
dat <- read.csv("~/Documents/Zika/Data/colombia/raw_births.csv")
daily_births <- rep(dat$births/dat$buckets,dat$buckets,each=TRUE)

births <- round(colSums(matrix(daily_births[3:681],nrow=7)))
write.csv(births,"weekly_births.csv")


## I used the Cuevas data to get monthly births for the period covered. However, I didn't have monthly births for the second half of November 2016 onwards. I therefore used the average number of weekly births from 2016 (647521/365)*7 for the remaining weeks

## OR using total data
## 2015 births
births_2015 <- 659255
births_2016 <- 647521

daily_2015 <- births_2015/365
daily_2016 <- births_2016/365

## EW 1 in 2016 starts on 02/01/2016
days_surveyed <- as.numeric(as.Date("2016-12-31",origin="2013-01-01")) - as.numeric(as.Date("2016-01-02",origin="2013-01-01"))
days_start <- as.numeric(as.Date("2016-01-02",origin="2013-01-01")) - as.numeric(as.Date("2016-01-01",origin="2013-01-01"))
days_end <- as.numeric(as.Date("2016-12-31",origin="2013-01-01")) - as.numeric(as.Date("2016-01-01",origin="2013-01-01"))
weekly_2016 <- round(rep(daily_2016*7, 52))
