## ----setup, echo=FALSE, message=FALSE, warning=FALSE---------------------
library(zikaProj)
library(lazymcmc)
library(kableExtra)
library(knitr)
library(zoo)
options(knitr.table.format = "html")
knitr::opts_chunk$set(fig.width=6, fig.height=4) 

## ----message=FALSE, warning=FALSE----------------------------------------
data(exampleParTab)

## Extract parameter values to solve model. Note that 
## these parameters MUST be named
pars <- exampleParTab$values
names(pars) <- exampleParTab$names

## ----message=FALSE, warning=FALSE----------------------------------------
## Generate starting population sizes based on human population size
## and mosquito density
y0s <- generate_y0s(pars["N_H"], pars["density"],iniI=10)
## Times to solve model over. Note that the unit of time is in days
ts <- seq(0,3003,by=1)
y <- solveSEIRModel_rlsoda(ts, y0s, pars, TRUE)
plot(y[,"I_H"],type='l',col="red",
     xlab="Time (days)", ylab="ZIKV incidence in humans")

## ----message=FALSE, warning=FALSE----------------------------------------
## Calculate R0
print(r0.calc(pars))

## Find the mosquito density required for a desired R0
print(density.calc(pars, R0=3))

print(calculate_AR(r0=3))

## ------------------------------------------------------------------------
risk <- generate_micro_curve(pars)
plot(risk,type='l',col="blue",
     xlab="Gestational time at infection (days)",ylab="Prob of congenital abnormality")

## ----message=FALSE, warning=FALSE----------------------------------------
probM <- generate_probM(y[,"I_M"],pars["N_H"],risk,pars["b"],pars["p_MH"],pars["baselineProb"],1)
plot(probM,type='l',col="green",
     xlab="Time (days)",ylab="Proportion of microcephaly affected births")

## ----message=FALSE, warning=FALSE----------------------------------------
library(ggplot2)

## Generate simulated data with known parameters
simDat <- generate_multiple_data(ts, parTab=exampleParTab, weeks=FALSE, 
                                  dataRangeMicro=c(600,2000),dataRangeInc=c(600,1500), 
                                  noise=FALSE,peakTimeRange=60)
microDat <- simDat[[1]]
incDat <- simDat[[2]]
peakTimes <- simDat[[3]]

## Save simulated data
write.table(microDat,"sim_microDat.csv",sep=",",row.names=FALSE)
write.table(incDat,"sim_incDat.csv",sep=",",row.names=FALSE)

## Check the generated data
ggplot(incDat) + 
  geom_line(aes(x=startDay,y=inc/N_H), col="red") +
  geom_line(data=microDat,aes(x=startDay,y=microCeph*0.001/births),col="blue") + 
  facet_wrap(~local) +
  ylab("Per capita ZIKV incidence (red)") +
  xlab("Time (days)") +
  scale_y_continuous(sec.axis=sec_axis(~.*1000,name="Per birth microcephaly\n incidence (blue)")) +
  theme_bw()

