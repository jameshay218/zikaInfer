## ----setup, echo=FALSE---------------------------------------------------
library(kableExtra)
library(knitr)
library(xtable)
options(knitr.table.format = "markdown")
vital_statistics <- read.csv("vital_statistics.csv")

## ---- echo=FALSE---------------------------------------------------------
data_sources <- read.csv("data_sources.csv", stringsAsFactors=FALSE)
data_sources <- data_sources[,colnames(data_sources) != "FullSource"]
colnames(data_sources) <- c("Country","Location","Incidence type","Resolution","Digitised?","Source")
kable(data_sources, caption="Summary of datasets included in the analysis")

## ---- echo=FALSE---------------------------------------------------------
parameters <- read.csv("model_parameters2.csv", stringsAsFactors=FALSE, check.names=FALSE)
kable(parameters, caption="Summary of model parameters, sources and assumed parameter ranges")

## ---- echo=FALSE---------------------------------------------------------
colnames(vital_statistics) <- c("State/Country","Population size (2015)a,b","Life expectancy (years) (2014)c,d")
vital_statistics <- vital_statistics[!(vital_statistics$`State/Country` %in% c("Colombia 2014","Colombia 2016","Ceará")),]
row.names(vital_statistics) <- NULL
kable(vital_statistics, caption="Summary of vital statistics and sources")

