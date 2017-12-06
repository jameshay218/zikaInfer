## CRS fitting
create_f <- function(parTab,data,PRIOR_FUNC,...){

  startDays <- data$startDay
  endDays <- data$endDay
  buckets <- data$buckets
  births <- data$births
  microCeph <- data$CRS_confirmed
  
  inc <- data$inc
  nh <- data$N_H
  inc_buckets <- data$buckets
  inc_start <- data$startDay
  inc_end <- data$endDay
  
  names_pars <- parTab$names
  f <- function(values){
    names(values) <- names_pars
    lik <- posterior_CRS(values, startDays, endDays,
                                    buckets, microCeph, births,
                                    inc, nh, inc_buckets,
                                    inc_start, inc_end)
    if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(values)
    return(lik)
  }
}
devtools::load_all("~/Documents/Zika/zikaInfer")

mcmcPars <- c("iterations"=100000,"adaptive_period"=50000,
               "popt"=0.44,"thin"=10,"save_block"=1000,
               "opt_freq"=1000)
data <- read.csv("~/Documents/Zika/Data/vietnam_data.csv")
parTab <- read.csv("~/net/home/zika/inputs/parTab_CRS.csv",stringsAsFactors = FALSE)
f <- create_f(parTab, data, NULL)
f(parTab$values)

omg <- antibodyKinetics::run_MCMC(parTab, data=data, mcmcPars, "test",create_f, NULL,NULL,0.2)
chain <- read.csv(omg$file)
#plot(coda::as.mcmc(chain))
tmp <- extra_microceph_calculations(chain,1,0.001,1)
#plot_random_microceph_curves(chain,100)


dat <- make_micro_dat(chain,1000)
dat$Model <- "Rubella"
chain <- cbind(chain, tmp)
wow <- reshape2::melt(chain)
wow <- wow[wow$variable %in% c("tr1","tr2","tr3"),]
wow$variable <- as.character(wow$variable)
wow <- plyr::ddply(wow,~variable, function(x) quantile(x$value,c(0.025,0.5,0.975)))
wow$Model <- "Rubella"

colnames(wow) <- c("variable","low","mid","high","Model")
