## ----setup, echo=FALSE---------------------------------------------------
library(kableExtra)
library(knitr)
options(knitr.table.format = "html")
vital_statistics <- read.csv("vital_statistics.csv")

## ---- echo=FALSE---------------------------------------------------------
data_sources <- read.csv("data_sources.csv", stringsAsFactors=FALSE)
data_sources <- data_sources[,colnames(data_sources) != "FullSource"]
colnames(data_sources) <- c("Country","Location","Incidence type","Resolution","Digitised?","Source")
kable(data_sources,"html", caption="Table 1. Summary of datasets included in the analysis") %>%
kable_styling(font_size=10)

## ---- echo=FALSE---------------------------------------------------------
parameters <- read.csv("model_parameters.csv", stringsAsFactors=FALSE)
kable(parameters,"html", caption="Table 2. Summary of model parameters") %>%
kable_styling(font_size=10)

## ---- echo=FALSE---------------------------------------------------------
colnames(vital_statistics) <- c("State/Country","Population size (2015)a,b","Life expectancy (years) (2014)c,d")
vital_statistics <- vital_statistics[!(vital_statistics$`State/Country` %in% c("Colombia 2014","Colombia 2016","Ceará")),]
kable(vital_statistics,"html", caption="Table 3. Vital statistics")

