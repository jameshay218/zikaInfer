

#' Plot single best trajectory
#'
#' For a given run name and state, plots the best zika and microcephaly incidence plot
#' @param local the name of the state to be fitted
#' @param chain the MCMC chain to use
#' @param realDat the data frame of real data that was fitted to
#' @param parTab the parameter table as returned by \code{\link{setupParTable}}
#' @param ts vector of time parameters
#' @param runs the number of runs to use for prediction intervals
#' @param incDat incidence data set, if to be included in the plot
#' @param ylabel boolean whether or not to add y axis label
#' @param xlabel boolean whether or not to add x axis label
#' @param mcmcPars a named vector of the burnin, adaptive period and thinning
#' @param ylimM ylimit for the microcephaly plot
#' @param ylimI ylimit for the incidence plot
#' @param startDay if realDat is not provided, need to provide the first day on which we want to see predicted microcephaly numbers
#' @param months if realDat not provided, need to provide the number of months of forecasted data that we wish to see
#' @param weeks if no incidence data provided, number of weeks over which we should plot the data
#' @param xlim OPTIONAL x coordinate limits
#' @return a ggplot object with the incidence plots
#' @export
plot_best_trajectory_single <- function(local, chain=NULL, realDat=NULL, parTab=NULL, ts=NULL, runs=100,incDat=NULL, ylabel=TRUE,xlabel=TRUE, mcmcPars=c("burnin"=50000,"adaptive_period"=100000,"thin"=50),ylimM=NULL, ylimI=NULL, startDay=NULL,months=NULL,weeks=NULL, xlim=NULL){
    allDat <- plot_setup_data(chain,realDat, incDat,parTab, ts,local,runs, startDay, months,weeks)
    bestMicro <- allDat$bestMicro
    bestInc <- allDat$bestInc
    incBounds <- allDat$incBounds
    microBounds <- allDat$microBounds
    dat <- allDat$data
    incDat <- allDat$incDat
    if(is.null(xlim)) xlim <- c(min(dat[,"startDay"]),max(dat[,"endDay"]))
    
    quantiles <- unique(microBounds[,"quantile"])
    botM <- microBounds[microBounds[,"quantile"]==quantiles[1],c("time","micro")]
    topM <- microBounds[microBounds[,"quantile"]==quantiles[2],c("time","micro")]

    botI <- incBounds[incBounds[,"quantile"]==quantiles[1],c("time","inc")]
    topI <-incBounds[incBounds[,"quantile"]==quantiles[2],c("time","inc")]
    polygonM <- create_polygons(botM, topM)
    polygonI <- create_polygons(botI, topI)

    tmp_p_I <- polygonI
    tmp_p_I <- tmp_p_I[tmp_p_I$x >= xlim[1] & tmp_p_I$x <= xlim[2],]
    polygonI <- tmp_p_I

    polygonM <- polygonM[polygonM$x >= xlim[1] & polygonM <= xlim[2],]

    xlabs <- generate_x_labels(xlim[1],xlim[2])

    myPlot <- microceph_plot(dat,microBounds,bestMicro,polygonM,local,xlim,ylimM,xlabs)
    incPlot <- inc_plot(incBounds,bestInc,polygonI,ylimI,xlim,incDat)

    if(!ylabel){
        myPlot <- myPlot + ylab("") + theme(plot.margin=unit(c(0,0.4,-0.2,-0.3),"cm"))
        incPlot <- incPlot + ylab("") + theme(plot.margin=unit(c(0,0.4,-0.2,-0.3),"cm"))
    }
    if(!xlabel) myPlot <- myPlot + xlab("")
    
    if(is.null(ylimM)) ylimM <- max(1.2*microBounds$micro)
    if(is.null(ylimI)) ylimI <- max(1.2*incBounds$inc)
    ylabInc <- NULL
    if(ylabel) ylabInc <- "Per capita incidence"
    myPlot <- add_inc_plot(myPlot, ylimM, incBounds,bestInc,polygonI,ylimI,ylabInc,incDat)
    return(myPlot)
}

microceph_plot <- function(dat, microBounds, bestMicro, polygonM, local, xlim, ylim, xlabs){
    xlabels <- xlabs$labels
    xlabBreaks <- xlabs$positions
    
    myPlot <- ggplot() + 
        geom_line(data=microBounds,aes_string(y="micro",x="time",group="quantile"),lwd=0.5,linetype=2,col="blue",alpha=0.5) +
        geom_line(data=bestMicro,aes_string(y="number",x="day"),col="blue",lwd=0.5) +
        geom_polygon(data=polygonM,aes_string(x="x",y="y"),alpha=0.2,fill="blue")+
        scale_x_continuous(labels=xlabels,breaks=xlabBreaks,limits=xlim)+
        theme_bw() +
        ylab("Microcephaly cases (blue)") +
        xlab("Date") +
        ggtitle(convert_name_to_state(local)) + 
        theme(
            panel.grid.minor=element_blank(),
            plot.title=element_text(size=12,hjust=0.5),
            axis.text.x=element_text(size=10,hjust=1,angle=45),
            axis.text.y=element_text(size=10),
            axis.title.x=element_text(size=10),
            axis.title.y=element_text(size=10),
            plot.margin=unit(c(0.1,0.8,0.1,0.1),"cm")
        )
    if(!is.null(dat$microCeph)) myPlot <- myPlot + geom_point(data=dat,aes_string(y="microCeph",x="meanDay"),col="black",size=1)
    if(!is.null(ylim)) myPlot <- myPlot + scale_y_continuous(limits=c(0,ylim))
    return(myPlot)
}

add_inc_plot <- function(p1, ylimMicro,incBounds, bestInc, polygonI, ylimI,ylab=NULL,incDat=NULL){
    modscale <- ylimI/ylimMicro

    tmpInc <- bestInc
    tmpInc$inc <- tmpInc$inc/modscale

    tmpBounds <- incBounds
    tmpBounds$inc <- tmpBounds$inc/modscale

    tmpPoly <- polygonI
    tmpPoly$y <- tmpPoly$y/modscale

    p2 <- p1 +  geom_line(data=tmpBounds,aes_string(y="inc",x="time",group="quantile"),linetype=2,col="red",size=0.5,alpha=0.5)+  
        geom_line(data=tmpInc,aes_string(y="inc",x="time"),col="red",lwd=0.5)+
        geom_polygon(data=tmpPoly,aes_string(x="x",y="y"),alpha=0.2,fill="red")
    if(!is.null(ylab)){
        p2 <- p2 +        
            scale_y_continuous(limits=c(0,ylimMicro),sec.axis=sec_axis(~.*modscale,name=ylab))
    } else {
        p2 <- p2 + scale_y_continuous(limits=c(0,ylimMicro),sec.axis=sec_axis(~.*modscale))

    }
    if(!is.null(incDat$inc)){
        incDat$meanDay <- rowMeans(incDat[,c("startDay","endDay")])
        incDat$perCapInc <- incDat[,"inc"]/incDat[,"N_H"]
        incDat$perCapInc <- incDat$perCapInc/modscale
        p2 <- p2 + geom_point(data=incDat,aes_string(x="meanDay",y="perCapInc"),col="black",size=1, shape=3)
    }
    return(p2)
}

inc_plot <- function(incBounds, bestInc, polygonI, ylimI,xlim, incDat=NULL){
    incPlot <- ggplot() +
        geom_line(data=incBounds,aes_string(y="inc",x="time",group="quantile"),linetype=2,col="red",size=0.5,alpha=0.5) +
        geom_line(data=bestInc,aes_string(y="inc",x="time"),col="red",lwd=0.5)+
        geom_polygon(data=polygonI,aes_string(x="x",y="y"),alpha=0.2,fill="red")+
        scale_x_continuous(limits=xlim)
    
    if(!is.null(ylimI)) incPlot <- incPlot + scale_y_continuous(limits=c(0,ylimI))

    incPlot <- incPlot + 
        ylab("\nPer capita incidence (red)")+
        xlab("")+
        theme(
            panel.background=element_rect(fill=NA),
            panel.grid=element_blank(),
            axis.line.y = element_line(colour="black"),
            axis.text.y=element_text(size=8,colour="black"),
            axis.title.y=element_text(size=10,angle=-90),
            axis.text.x=element_blank(),
            axis.title.x=element_text(size=10),
            plot.margin=unit(c(0.1,0.8,0.1,0.1),"cm")
            
        )
    
    if(!is.null(incDat$inc)){
        incDat$meanDay <-rowMeans(incDat[,c("startDay","endDay")])
        incDat$perCapInc <- incDat[,"inc"]/incDat[,"N_H"]
        incPlot <- incPlot + geom_point(data=incDat,aes_string(x="meanDay",y="perCapInc"),col="black",size=1, shape=3)
    }
    return(incPlot)
}


overlapPlots <- function(p1,p2, centre_plot=TRUE){
    
                                        # extract gtable
    g1 <- ggplot_gtable(ggplot_build(p1))
    g2 <- ggplot_gtable(ggplot_build(p2))
    
                                        # overlap the panel of 2nd plot on that of 1st plot
    
    name=se=t=r=NULL
    pp <- c(subset(g1$layout, name == "panel", se = t:r))
    g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
                         pp$l, pp$b, pp$l)
    if(centre_plot) return(g)
    
                                        # axis tweaks
    ia <- which(g2$layout$name == "axis-l")
    if(!is.null(ia) & !is.null(ia)){
        ga <- g2$grobs[[ia]]
        ax <- ga$children[[2]]
        ax$widths <- rev(ax$widths)
        ax$grobs <- rev(ax$grobs)
        ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.1, "cm")
    }

    g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
    g <- gtable_add_grob(g, ax, pp$t, length(g$widths)-1, pp$b)
    g <- gtable_add_grob(g, g2$grobs[[7]], pp$t, length(g$widths), pp$b)
    
                                        # draw it
    return(g)
}

#' Plot all best trajectoroes
#'
#' For a given run name and state, plots the best zika and microcephaly incidence for all states
#' @param chain the MCMC chain of estimated parameters
#' @param realDat the data frame of data used for fitting
#' @param parTab the parameter table used for fitting
#' @param ts vector of times to solve model over
#' @param runs the number of runs to use for prediction intervals
#' @param incDat optional data frame of actual incidence data
#' @param mcmcPars a named vector of the burnin, adaptive period and thinning
#' @param ylimM ylimit for the microcephaly plot
#' @param ylimI ylimit for the incidence plot
#' @param startDay if realDat is not provided, need to provide the first day on which we want to see predicted microcephaly numbers
#' @param months if realDat not provided, need to provide the number of months of forecasted data that we wish to see
#' @param weeks if no incidence data provided, number of weeks over which we should plot the data
#' @return a ggplot object with the incidence plots
#' @export
plot_best_trajectory_multi <- function(chain, realDat, parTab, ts, runs=100, incDat=NULL, mcmcPars=c("burnin"=50000,"adaptive"=100000,"thin"=50), ylimM=NULL,ylimI=NULL,startDay=NULL,months=NULL,weeks=NULL,ncol=4,xlim=NULL){
    ps <- NULL
    states <- unique(parTab$local)
    states <- states[states != "all"]

    for(i in 1:length(states)){
        ps[[i]] <- plot_best_trajectory_single(states[i], chain, realDat, parTab, ts, runs, incDat=incDat, ylabel=FALSE, xlabel=FALSE, mcmcPars,ylimM,ylimI,startDay,months,weeks,xlim)
    }

    if(length(states) <= 1) return(ps[[1]])

    ncols <- ceiling(length(states)/ncol)
    order <- get_correct_order(FALSE,TRUE)
    order1 <- as.numeric(sapply(order, function(x) which(x == states)))
    order1 <- order1[!is.na(order1)]
    ps <- ps[order1]
    allPlot <- do.call("plot_grid",c(ps,ncol=ncols))
    return(allPlot)
}




create_polygons <- function(lower,upper){
    bounds <- NULL
    bounds <- rbind(lower[rev(rownames(lower)),],upper)
    colnames(bounds) <- c("x","y")
    return(bounds)    
}

    
#' Adaptive Metropolis-within-Gibbs Random Walk Algorithm.
#'
#' The Adaptive Metropolis-within-Gibbs algorithm. Given a starting point and the necessary MCMC parameters as set out below, performs a random-walk of the posterior space to produce an MCMC chain that can be used to generate MCMC density and iteration plots. The algorithm undergoes an adaptive period, where it changes the step size of the random walk for each parameter to approach the desired acceptance rate, popt. After this, a burn in period is established, and the algorithm then uses \code{\link{proposalfunction}} to explore the parameter space, recording the value and posterior value at each step. The MCMC chain is saved in blocks as a .csv file at the location given by filename.
#' @param data the microcephaly data against which the likelihood is calculated
#' @param ts vector of times to solve the ODE model over
#' @param parTab a table of parameter data used for information such as bounds and prior function pointers.
#' @param mcmcPars a named vector with parameters for the MCMC procedure. Iterations, popt, opt_freq, thin, burnin, adaptive_period and save_block.
#' @param filename the full filepath at which the MCMC chain should be saved. "_chain.csv" will be appended to the end of this, so filename should have no file extensions
#' @param mvrPars a list of parameters if using a multivariate proposal. Must contain an initial covariance matrix, weighting for adapting cov matrix, and an initial scaling parameter (0-1)
#' @param incDat optional data frame of incidence data if including incidence data in the likelihood function
#' @param peakTimes optional parameter - data frame of peak times for Zika incidence for each state
#' @param allPriors user function of prior for model parameters. Should take values, names and local from parTab
#' @param version usually just leave this as "normal". I've added a "forecast" version which expects parameters to do with the lack of second wave.
#' @return a list with: 1) full file path at which the MCMC chain is saved as a .csv file; 2) the last used covarianec matrix; 3) the last used scale size
#' @seealso \code{\link{posterior_complex_buckets}}, \code{\link{proposalfunction}}
run_metropolis_MCMC <- function(data=NULL,
                                ts,
                                parTab,
                                mcmcPars,
                                filename,
                                mvrPars=NULL,
                                incDat=NULL,
                                peakTimes=NULL,
                                allPriors=NULL,
                                truePars=NULL,
                                version="normal"
                                ){
                                        # MCMC par setup ---------------------------------------------------------- 
    ## Allowable error in scale tuning
    TUNING_ERROR <- 0.1
    OPT_TUNING  <- 0.2
    
    ## Extract MCMC parameters
    iterations <- mcmcPars["iterations"]
    popt <- mcmcPars["popt"]
    opt_freq<- mcmcPars["opt_freq"]
    thin <- mcmcPars["thin"]
    adaptive_period<- mcmcPars["adaptive_period"]
    save_block <- mcmcPars["save_block"]

    
                                        # Parameter par setup -------------------------------------
    ## Time vector for ODE model
    ## Extract parameter table into individual vectors for speed
    tmp_non_fixed_params <- non_fixed_params <- which(parTab$fixed==0)
    par_names <- parTab$names
    current_params <- start_values <- parTab$values
    par_labels <- parTab$local
    lower_bounds <- parTab$lower_bounds
    upper_bounds <- parTab$upper_bounds
    steps <- parTab$steps
    fixed <- parTab$fixed
    names(current_params) <- names(start_values) <- par_names
    unique_states <- unique(par_labels)
    unique_states <- unique_states[unique_states != "all"]
    all_states <- unique(parTab[parTab$fixed==0,"local"])
    ## Get some parameters useful for indexing
    non_fixed_params_length <- length(non_fixed_params)
    all_param_length <- nrow(parTab)

  
    
    ## Arrays to store acceptance rates
    ## If univariate proposals
    if(is.null(mvrPars)){
        tempaccepted <- tempiter <- integer(all_param_length)
        reset <- integer(all_param_length)
        reset[] <- 0
    } else { # If multivariate proposals
        tempaccepted <- tempiter <- 0
        covMat <- mvrPars[[1]][non_fixed_params,non_fixed_params]
        scale <- mvrPars[[3]]
    }
    
                                        # Data extraction ---------------------------------------------------------
    ## Extract data into vectors
    startDays <- NULL
    endDays <- NULL
    buckets <- NULL
    microCeph <- NULL
    births <- NULL
    data_locals <- NULL

    ## Microceph data
    if(!is.null(data)){
        startDays <- data[,"startDay"]
        endDays <- data[,"endDay"]
        buckets <- data[,"buckets"]
        microCeph <- data[,"microCeph"]
        births <- data[,"births"]
        data_locals <- data[,"local"]
    }

    ## Peak times
    peak_startDays <- NULL
    peak_endDays <- NULL
    peak_locals <- NULL
    if(!is.null(peakTimes)){
        peak_locals <- peakTimes$local
        peak_startDays <- peakTimes$start
        peak_endDays <- peakTimes$end
    }

    ## Incidence data
    inc_startDays <- NULL
    inc_endDays <- NULL
    inc_locals <- NULL
    inc_buckets <- NULL
    inc_ZIKV <- NULL
    inc_NH <- NULL
    if(!is.null(incDat)){
        inc_startDays <- incDat[,"startDay"]
        inc_endDays <- incDat[,"endDay"]
        inc_locals <- incDat[,"local"]
        inc_buckets <- incDat[,"buckets"]
        inc_ZIKV <- incDat[,"inc"]
        inc_NH <- incDat[,"N_H"]
    }
                                        # Posterior setup ---------------------------------------------------------
    ## Create posterior function with closures for neatness
    posterior_new <- create_posterior(ts, current_params, par_names, par_labels, 
                                      startDays, endDays, buckets, microCeph, births, 
                                      data_locals, inc_startDays,inc_endDays,inc_locals,
                                      inc_buckets,inc_ZIKV,inc_NH, peak_startDays, 
                                      peak_endDays,peak_locals, unique_states, allPriors,
                                      version
                                      )

    posterior_simp <- protect(posterior_new)
    
                                        # Chain setups ------------------------------------------------------------
    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste(filename,"_chain.csv",sep="")
    chain_colnames <- c("sampno",parTab$names,"lnlike")
    
    ## Create empty chain to store every iteration for the adaptive period
    opt_chain <- matrix(nrow=adaptive_period,ncol=non_fixed_params_length)
    chain_index <- 1
    
    ## Create empty chain to store "save_block" iterations at a time
    save_chain <- empty_save_chain <- matrix(nrow=save_block,ncol=all_param_length+2)
                                        # Initial conditions ------------------------------------------------------
    ## Initial likelihood
    probab <- posterior_simp(current_params)

    log_probab <- 0
    ## If available, find the true parameter posterior for comparison
    true_probab <- NULL
    if(!is.null(truePars)){
        true_probab <- posterior_simp(truePars)
        message(cat("True parameter posterior: ",true_probab,sep="\t"))
    }
    ## Set up initial csv file
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
  
    tmp_table[1,] <- c(1,as.numeric(current_params),probab)
    colnames(tmp_table) <- chain_colnames
    
    ## Write starting conditions to file
    write.table(tmp_table,file=mcmc_chain_file,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)

    ## Initial indexing parameters
    no_recorded <- 1
    sampno <- 2
    par_i <- 1
                                        # Main MCMC algorithm -----------------------------------------------------
    ## Go through chain
    for (i in 1:(iterations+adaptive_period)){
        ## If using univariate proposals
        if(is.null(mvrPars)) {
            ## For each parameter (Gibbs)
            j <- non_fixed_params[par_i]
            par_i <- par_i + 1
            if(par_i > length(non_fixed_params)) par_i <- 1
            proposal <- proposalfunction(current_params, lower_bounds, upper_bounds, steps,j)
            tempiter[j] <- tempiter[j] + 1
            ## If using multivariate proposals
        } else {
            proposal <- mvr_proposal(current_params, non_fixed_params, scale*covMat)
            tempiter <- tempiter + 1
        }
        names(proposal) <- names(current_params)
        ## Propose new parameters and calculate posterior
        ## Check that all proposed parameters are in allowable range
        if(!any(
                proposal[non_fixed_params] < lower_bounds[non_fixed_params] |
                proposal[non_fixed_params] > upper_bounds[non_fixed_params]
            )
           ){
            ## Calculate new likelihood and find difference to old likelihood
            new_probab <- posterior_simp(proposal)
            log_prob <- min(new_probab-probab,0)
            if(!is.finite(log_prob)){
                #message("Not finite")
                #message(cat(proposal[c("density","t0")]," "))
            }
            ## Accept with probability 1 if better, or proportional to
            ## difference if not
            if(is.finite(log_prob) && log(runif(1)) < log_prob){
                current_params <- proposal
                probab <- new_probab
                
                ## Store acceptances
                if(is.null(mvrPars)){
                    tempaccepted[j] <- tempaccepted[j] + 1
                } else {
                    tempaccepted <- tempaccepted + 1
                }
            }
        }
        
        ## If current iteration matches with recording frequency, store in the chain. If we are at the limit of the save block,
        ## save this block of chain to file and reset chain
        if(i %% thin ==0){
            save_chain[no_recorded,1] <- sampno
            save_chain[no_recorded,2:(ncol(save_chain)-1)] <- current_params
            save_chain[no_recorded,ncol(save_chain)] <- probab
            no_recorded <- no_recorded + 1
        }
        
        ## If within adaptive period, need to do some adapting!

        if(i <= adaptive_period){
            ## Current acceptance rate
            pcur <- tempaccepted/tempiter
            ## Save each step
            opt_chain[chain_index,] <- current_params[non_fixed_params]
            
            ## If in an adaptive step
            if(chain_index %% opt_freq == 0){
                ## If using univariate proposals
                if(is.null(mvrPars)){
                    ## For each non fixed parameter, scale the step size
                    for(x in non_fixed_params) steps[x] <- scaletuning(steps[x],popt,pcur[x])
                    message(cat("Optimisation iteration: ", i,sep="\t"))
                    message(cat("Pcur: ", pcur[non_fixed_params],sep="\t"))
                    message(cat("Step sizes: ", steps[non_fixed_params],sep="\t"))
                    tempaccepted <- tempiter <- reset
                } else {       ## If using multivariate proposals
                    if(chain_index > OPT_TUNING*adaptive_period & chain_index < (0.9*adaptive_period)){
                        #covMat <- scale*cov(opt_chain[1:chain_index,])
                        covMat <- cov(opt_chain[1:chain_index,])
                        tempiter <- tempaccepted <- 0
                        message(cat("Optimisation iteration: ", i,sep="\t"))
                        ## Print acceptance rate
                        message(cat("Pcur: ", pcur,sep="\t"))
                        message(cat("Step size: ", scale,sep="\t"))
                    }
                    if(chain_index > (0.9)*adaptive_period){
                        scale <- scaletuning(scale, popt,pcur)
                        message(cat("Scale: ",scale,sep=""))
                    }
                }
            }
            chain_index <- chain_index + 1
        }
        if(i %% save_block == 0) message(cat("Current iteration: ", i, sep="\t"))
        if(no_recorded == save_block){
            write.table(save_chain[1:(no_recorded-1),],file=mcmc_chain_file,col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
            save_chain <- empty_save_chain
            no_recorded <- 1
        }
        sampno <- sampno + 1
    }
    
    ## If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
    ## that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
    ## rather than an array
    if(no_recorded > 2){
        write.table(save_chain[1:(no_recorded-1),],file=mcmc_chain_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
    }

    if(is.null(mvrPars)){
        covMat <- NULL
        scale <- NULL
    } else {
        steps <- NULL
    }
    return(list("file"=mcmc_chain_file,"covMat"=covMat,"scale"=scale, "steps"=steps,"truePosterior"=true_probab))
}






#' Scale step sizes
#'
#' Scales the given step size (between 0 and 1) based on the current acceptance rate to get closed to the desired acceptance rate
#' @param step the current step size
#' @param popt the desired acceptance rate
#' @param pcur the current acceptance rate
#' @return the scaled step size
scaletuning <- function(step, popt,pcur){
    if(pcur ==1) pcur <- 0.99
    if(pcur == 0) pcur <- 0.01
    step = (step*qnorm(popt/2))/qnorm(pcur/2)
    if(step > 1) step <- 1
    return(step)
}




#' Step size scaler
#'
#' Given current step scale, desired optimum step size etc, adapts the MCMC step size
#' @param step_scale the current step size
#' @param mc the number of MCMC iterations run
#' @param popt the desired acceptance rate
#' @param log_prob the current log probability
#' @param N_adapt the total number of adaptive steps to be taken
#' @export
rm_scale <- function(step_scale, mc, popt,log_prob, N_adapt)
{
    dd <- exp(log_prob)
    if( dd < -30 ){ dd <- 0 }
    dd <- min( dd, 1 )

    rm_temp <- ( dd - popt )/( (mc+1)/(0.01*N_adapt+1) )
    
    out <- step_scale*exp(rm_temp)
    
    out <- max( out, 0.02 )
    out <- min( out, 2)
    out
}

#' Multivariate proposal function
#'
#' Given the current parameters and a covariance matrix, returns a vector for a proposed jump from a multivariate normal distribution
#' @param values the vector of current parameter values
#' @param fixed set of flags corresponding to the parameter vector indicating which parameters are fixed
#' @param covMat the 2D covariance matrix for all of the parameters
#' @return a parameter vector of a proposed move. Note that these may fall outside the allowable ranges.
#' @export
mvr_proposal <- function(values, fixed, covMat){
    proposed <- values
    proposed[fixed] <- MASS::mvrnorm(n=1,mu=proposed[fixed],Sigma=(5.6644/length(fixed))*covMat)
    return(proposed)
}

#' MCMC proposal function
#'
#' Proposal function for MCMC random walk, taking random steps of a given size. Random walk may be on a linear or log scale
#' @param values a vector of the parameters to be explored
#' @param lower_bounds a vector of the low allowable bounds for the proposal
#' @param upper_bounds a vector of the upper allowable bounds for the proposal
#' @param steps a vector of step sizes for the proposal
#' @param index numeric value for the index of the parameter to be moved from the param table and vector
#' @return the parameter vector after step
#' @export
proposalfunction <- function(values, lower_bounds, upper_bounds,steps, index){
    mn <- lower_bounds[index]
    mx <- upper_bounds[index]

    rtn <- values
    
    x <- toUnitScale(values[index],mn,mx)

    ## 5th index is step size
    stp <- steps[index]

    rv <- runif(1)
    rv <- (rv-0.5)*stp
    x <- x + rv

    ## Bouncing boundary condition
    if (x < 0) x <- -x
    if (x > 1) x <- 2-x

    ## Cyclical boundary conditions
    ##if (x < 0) x <- 1 + x	
    ##if (x > 1) x <- x - 1
    
    if(x < 0 | x > 1) print("Stepped outside of unit scale. Something went wrong...")

    rtn[index] <- fromUnitScale(x,mn,mx)
    rtn
}
