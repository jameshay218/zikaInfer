library(data.table)
library(lazymcmc)
library(zikaProj)
library(coda)

chainwd <- "/media/james/JH USB/outputs/northeast_forecast"
parTab <- read_inipars(chainwd)
chain <- lazymcmc::load_mcmc_chains(chainwd, parTab, TRUE, thin=1,burnin=50000,
                                    TRUE, FALSE, FALSE)

ess <- ess_diagnostics(chain[[1]],200)
gelman <- gelman_diagnostics(chain[[1]], 1.1)
