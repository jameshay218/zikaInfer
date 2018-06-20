setwd("~/net/home/zika")
source("~/Documents/Zika/zikaInfer/scripts/job_submissions/model_fitting_script.R")
source("~/Documents/Zika/zikaInfer/scripts/job_submissions/cluster_setup.R")

#setwd("~/Documents/Zika/zikaInfer/scripts/job_submissions/runs_18112017")
setwd("~/net/home/zika")

mcmcPars <- c("adaptive_period"=300000,"iterations"=700000,"opt_freq"=1000,
              "thin"=50,"save_block"=100,"popt"=0.44)
mcmcPars2 <- c("adaptive_period"=750000,"iterations"=2000000,"opt_freq"=1000,
               "thin"=50,"save_block"=100,"popt"=0.234)

## Single state jobs
source("bahia.R") # Yes
source("nejm.R") # Yes
source("colombia.R") # Yes - note jobs from here are a list
source("rio.R") # Yes
source("pernambuco.R") # Yes

## Multi state jobs
source("multi_3.R") # Yes - list of jobs
source("reports_3.R") # Yes - list of jobs
source("reports_2.R") # Yes - list of jobs


mcmcPars1 <- c("iterations"=1000000,"adaptive_period"=500000,
              "popt"=0.44,"thin"=10,"save_block"=1000,
              "opt_freq"=1000)
mcmcPars2 <- c("iterations"=2000000,"adaptive_period"=500000,
              "popt"=0.234,"thin"=10,"save_block"=1000,
              "opt_freq"=1000)

## Forecasts
source("bahia_forecasts.R") # Yes
source("nejm_forecasts.R") # Yes
source("rio_forecasts.R") # Yes

source("bahia_forecasts_noreportingchange.R") # Yes
source("nejm_forecasts_noreportingchange.R") # Yes

all_jobs <- list(jobs_bahia, jobs_nejm, jobs_colombia,
                 jobs_rio,jobs_pernambuco,
                 jobs_multi3, jobs_reports3, jobs_reports2,
                 jobs_bahia_forecast, jobs_nejm_forecasts,
                 jobs_rio_forecasts,  jobs_nejm_forecastsB,
                 jobs_bahia_forecastsB)

statuses <- NULL
for(i in 1:length(all_jobs)){
    x <- all_jobs[[i]]
    if(is.list(x)){
        for(j in 1:length(x)){
            y <- x[[j]]
            statuses <- c(statuses, y$status())
        }
    } else {
        statuses <- c(statuses, x$status())
    }
}

which(statuses == "ERROR")
