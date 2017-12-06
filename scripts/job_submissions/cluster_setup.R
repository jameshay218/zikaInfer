## Install the necessary packages from Rich's github
#devtools::install_github(c("gaborcsardi/progress","dide-tools/context","richfitz/queuer","dide-tools/didewin"))
## Submit credentials to didewin cluster
options(didehpc.credentials = "~/.smbcredentials",
        #didehpc.cluster = "fi--didemrchnb")
        didehpc.cluster = "fi--dideclusthn")

src <- provisionr::package_sources(local = c("~/Documents/Zika/zikaInfer", "~/Documents/R_Packages/lazymcmc","~/src/rlsoda"))
sources <- c("scripts/model_fitting_script.R","scripts/model_forecasting.R")
## Setup contexts
context::context_log_start()
root <- "contexts_test"         
ctx <- context::context_save(packages=c("zikaInfer","lazymcmc","plyr"), path=root, sources=sources,package_sources=src)
## Submit setup to cluster
obj1 <- didehpc::queue_didehpc(ctx)
#t <- obj1$enqueue(sessionInfo())
