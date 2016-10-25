# PDF - Rich's function to print to device without potential for bad errors
to.pdf <- function(expr, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}
# PNG - Rich's function to print to device without potential for bad errors
to.png <- function(expr, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  png(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}


#' Plot MCMC chains
#'
#' Given an MCMC chain and optionally the parameter table, plots the MCMC chains and densities using coda
#' @param chains list of MCMC chains
#' @param parTab optional parameter table
#' @param filename the name of the file to save to
#' @return nothing - saves plots to pdf
#' @export
plot_MCMC <- function(chains, parTab=NULL,filename="mcmcplots.pdf"){
    for(i in 1:length(chains)){
        if(!is.null(parTab)){
            indices <- which(parTab$fixed==0)
            chains[[i]] <- chains[[i]][,c(indices+1, ncol(chains[[i]]))]
        }
        chains[[i]] <- coda::as.mcmc(chains[[i]])
    }
    to.pdf(plot(coda::as.mcmc.list(chains)),filename)
}

#' Plot microceph curves
#'
#' Given an MCMC chain, plots the best fitting microcephaly risk curve and some faint, random draws from the posterior
#' @param chain The MCMC chain. Must have the gamma distribution mean, variance and scale
#' @param runs number of random draws to plot
#' @return a ggplot object with the microcephaly risk curve
#' @export
plot_random_microceph_curves <- function(chain, runs){
    ## Get best fitting parameters
    bestPars <- as.numeric(chain[which.max(chain[,"lnlike"]),])
   
    names(bestPars) <- colnames(chain)
    bestPars["tstep"] <- 1
    probs <- generate_micro_curve(bestPars)
    bestProbs <- data.frame(prob=probs,week=seq(0,279,by=1))
    myPlot <- ggplot() +
        geom_line(data=bestProbs, aes_string(x="week",y="prob"),col="blue",lwd=1) +
        ylab("Probability of microcephaly given infection")+
        xlab("Week of gestation at infection") +
        ggtitle("Microcephaly Risk Curve") +
        theme_bw()+
        theme(axis.text.x=element_text(size=6),
              panel.grid.minor=element_blank(),
              plot.title=element_text(size=10),
              axis.text.x=element_text(size=6),
              axis.text.y=element_text(size=6),
              axis.title.x=element_text(size=8),
              axis.title.y=element_text(size=8))
   
    ## Add some random lines
    samples <- sample(nrow(chain),runs)
    index <- 1
    allProbs <- NULL
    for(i in samples){
        tmpPars <- as.numeric(chain[i,])
        names(tmpPars) <- colnames(chain)
        tmpPars["tstep"] <- 1
        probs <- generate_micro_curve(tmpPars)
        allProbs[[index]] <- probs <- data.frame(prob=probs,week=seq(0,279,by=1))
        myPlot <- myPlot + geom_line(data=probs, aes_string(x="week",y="prob"),alpha=0.3,lwd=0.5,col="red")
        index <- index+1
    }
  return(myPlot)
}

#' Translates state code to name
#'
#' As described
#' @param state the state code to be changed
#' @return a string
#' @export
get_state_name <- function(state){
    country_names <- c("pernambuco"="Pernambuco","amapa"="Amapá","amazonas" ="Amazonas",
                     "distritofederal" = "Distrito Federal","bahia"="Bahia","saopaulo"="São Paulo",
                     "paraiba"="Paraíba","maranhao"="Maranhão","ceara"="Ceará","sergipe"="Sergipe",
                     "riodejaneiro"="Rio de Janeiro","piaui"="Piauí","riograndedonorte"="Rio Grande Norte",
                     "minasgerais"="Minas Gerais", "matogrosso"="Mato Grosso","alagoas"="Alagoas",
                     "para"="Pará","acre"="Acre","espiritosanto"="Espírito Santo","goias"="Goiás",
                     "tocantins"="Tocantins","matogrossodosul"="Mato Grosso do Sul",
                     "matogrossdosul"="Mato Grosso do Sul","parana"="Paraná","riograndedosul"="Rio Grande do Sul",
                     "rondonia"="Rondônia","roraima"="Roraima","santacatarina"="Santa Catarina")
    return(country_names[state])
}


create_polygons <- function(lower,upper){
    bounds <- NULL
    bounds <- rbind(lower[rev(rownames(lower)),],upper)
    colnames(bounds) <- c("x","y")
    return(bounds)    
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
plot_best_trajectory_multi <- function(chain, realDat, parTab, ts, runs=100, incDat=NULL, mcmcPars=c("burnin"=50000,"adaptive"=100000,"thin"=50), ylimM=NULL,ylimI=NULL,startDay=NULL,months=NULL,weeks=NULL){
    ps <- NULL
    states <- unique(parTab$local)
    states <- states[states != "all"]

    for(i in 1:length(states)){
        ps[[i]] <- plot_best_trajectory_single(states[i], chain, realDat, parTab, ts, runs, incDat=incDat, ylabel=FALSE, xlabel=FALSE, mcmcPars,ylimM,ylimI,startDay,months,weeks)
    }
    ncols <- ceiling(length(states)/4)
    allPlot <- do.call("grid.arrange",c(ps,ncol=ncols))
    return(allPlot)
}

generate_x_labels <- function(startDay, endDay, rep_=6){
    days <- getDaysPerMonth(12)
    buckets <- rep(days, 5)
    xpositions <- cumsum(buckets) - buckets
    indices <- which(xpositions >= startDay & xpositions <= endDay)
    years <- c("2013","2014","2015","2016","2017")
    months <- c("01/","02/","03/","04/","05/","06/","07/","08/","09/","10/","11/","12/")
    labels <- apply(expand.grid(months,years),1,paste,collapse="")    

    new_i <- seq(indices[1],indices[length(indices)],by=rep_)
       
    return(list("labels"=labels[new_i],"positions"=xpositions[new_i]))
}


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
#' @return a ggplot object with the incidence plots
#' @export
plot_best_trajectory_single <- function(local, chain=NULL, realDat=NULL, parTab=NULL, ts=NULL, runs=100,incDat=NULL, ylabel=TRUE,xlabel=TRUE, mcmcPars=c("burnin"=50000,"adaptive_period"=100000,"thin"=50),ylimM=NULL, ylimI=NULL, startDay=NULL,months=NULL,weeks=NULL){
    allDat <- plot_setup_data(chain,realDat, incDat,parTab, ts,local,runs, startDay, months,weeks)
    bestMicro <- allDat$bestMicro
    bestInc <- allDat$bestInc
    incBounds <- allDat$incBounds
    microBounds <- allDat$microBounds
    dat <- allDat$data
    incDat <- allDat$incDat

    xlim <- c(min(dat[,"startDay"]),max(dat[,"endDay"]))
    
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

    #xlim <- c(min(dat[,"startDay"]),max(dat[,"endDay"]))
    xlabs <- generate_x_labels(xlim[1],xlim[2])
    myPlot <- microceph_plot(dat,microBounds,bestMicro,polygonM,local,xlim,ylimM,xlabs)
    incPlot <- inc_plot(incBounds,bestInc,polygonI,ylimI,xlim,incDat)

    if(!ylabel){
        myPlot <- myPlot + ylab("") + theme(plot.margin=unit(c(0,0.4,-0.2,-0.3),"cm"))
        incPlot <- incPlot + ylab("") + theme(plot.margin=unit(c(0,0.4,-0.2,-0.3),"cm"))
    }
    if(!xlabel) myPlot <- myPlot + xlab("")
    
    combined <- overlapPlots(myPlot,incPlot,FALSE)
    
    return(combined)
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
        ggtitle(get_state_name(local)) + 
        theme(
            panel.grid.minor=element_blank(),
            plot.title=element_text(size=12),
            axis.text.x=element_text(size=8,hjust=1,angle=45),
            axis.text.y=element_text(size=8),
            axis.title.x=element_text(size=10),
            axis.title.y=element_text(size=10),
            plot.margin=unit(c(0.1,0.8,0.1,0.1),"cm")
        )
    if(!is.null(dat$microCeph)) myPlot <- myPlot + geom_point(data=dat,aes_string(y="microCeph",x="meanDay"),col="black",size=1)
    if(!is.null(ylim)) myPlot <- myPlot + scale_y_continuous(limits=c(0,ylim))
    return(myPlot)
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
            axis.text.x=element_text(size=8),
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
    ga <- g2$grobs[[ia]]
    ax <- ga$children[[2]]
    ax$widths <- rev(ax$widths)
    ax$grobs <- rev(ax$grobs)
    ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.1, "cm")
    
    ib <- which(g2$layout$name=="ylab")
    
    gb <- g2$grobs[[ib]]
    
    g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
    g <- gtable_add_grob(g, ax, pp$t, length(g$widths)-1, pp$b)
    g <- gtable_add_grob(g, g2$grobs[[7]], pp$t, length(g$widths), pp$b)
    
                                        # draw it
    return(g)
}


#' Get plot data
#'
#' Given an MCMC chain and other parameters, returns a list of data frames that can be used to generate an incidence plot
#' @param chain the MCMC chain
#' @param dat the data frame of real microcephaly data
#' @param incDat the data frame of real ZIKV incidence data
#' @param parTab the parameter table
#' @param ts vector of times
#' @param local string for state under consideration
#' @param runs the number of samples to take
#' @param startDay if realDat is not provided, need to provide the first day on which we want to see predicted microcephaly numbers
#' @param noMonths if realDat not provided, need to provide the number of months of forecasted data that we wish to see
#' @param noWeeks if no incidence data provided, number of weeks over which we should plot the data
#' @export
plot_setup_data <- function(chain, dat=NULL, incDat=NULL, parTab, ts, local, runs=NULL,startDay=NULL, noMonths=12,noWeeks=52){
    ## Get the subset of data for this particular state and find the middle day for each month
    ## If not provided, generate artificial sampling times
    tmpDat <- NULL
    if(!is.null(dat)){
        tmpDat <- dat[dat$local==local,]
        tmpDat$meanDay <- rowMeans(tmpDat[,c("startDay","endDay")])
    } else {
        months <- rep(getDaysPerMonth(),pmax(1,noMonths/12))
        startDays <- startDay + cumsum(months) - months[1]
        endDays <- startDay + cumsum(months)
        tmpDat <- data.frame(startDay=startDays,endDay=endDays)
        tmpDat$buckets <- months
        tmpDat$meanDay <- rowMeans(tmpDat[,c("startDay","endDay")])
    }

    # Get subset of incidence data for this state.
    tmpInc <- NULL
    if(!is.null(incDat)){
        tmpInc <- incDat[incDat$local == local,]
        tmpInc$meanDay <- rowMeans(tmpInc[,c("startDay","endDay")])
    } else {
        inc_weeks <- rep(7,noWeeks)

        startDays <- startDay + cumsum(inc_weeks) - inc_weeks[1]
        endDays <- startDay + cumsum(inc_weeks)
        tmpInc <- data.frame(startDay=startDays,endDay=endDays)
        tmpInc$buckets <- inc_weeks
        tmpInc$meanDay <- rowMeans(tmpInc[,c("startDay","endDay")])
    }        

    ## Get the set of best fitting parameters. These will be used to generate the best fit incidence and
    ## microcephaly line
    bestPars <- get_best_pars(chain)
    tmp <- plot_setup_data_auxiliary(bestPars,tmpDat,tmpInc,parTab,ts,local)
    bestInc <- tmp$inc
    bestMicro <- tmp$micro
    
    ## Get sample indices with which to plot prediction intervals
    samples <- sample(nrow(chain),runs)

    allInc <- NULL
    allMicro <- NULL
    ## For each sample, get the sample parameters from the chain and calculate the trajectory
    for(i in samples){
        pars <- get_index_pars(chain,i)
        tmp <- plot_setup_data_auxiliary(pars,tmpDat,tmpInc,parTab,ts,local)
        inc <- tmp$inc
        micro <- tmp$micro
        allInc <- rbind(allInc, inc)
        allMicro <- rbind(allMicro, micro)
    }

    incBounds <- as.data.frame(reshape2::melt(sapply(unique(allInc$time),function(x) quantile(allInc[allInc$time==x,"inc"],c(0.025,0.5,0.975)))[c(1,3),]))
    microBounds <- as.data.frame(reshape2::melt(sapply(unique(allMicro$day),function(x) quantile(allMicro[allMicro$day==x,"number"],c(0.025,0.5,0.975)))[c(1,3),]))
    colnames(microBounds) <- c("quantile","time","micro")
    colnames(incBounds) <- c("quantile","time","inc")
    
    microBounds[,"time"] <- tmpDat[,"meanDay"][microBounds[,"time"]]
    incBounds[,"time"] <- tmpInc[,"meanDay"][incBounds[,"time"]]
    
    return(list("bestInc"=bestInc,"incBounds"=incBounds,"bestMicro"=bestMicro,"microBounds"=microBounds, "data"=tmpDat,"incDat"=tmpInc))
}

plot_setup_data_auxiliary <- function(pars, dat,incDat, parTab, ts, local){
    ## As the model has many parameters with the same name, need to find the index in the MCMC colnames
    ## that corresponds to this state
    number <- which(unique(parTab$local)==local) - 2
  
    ## Format parameter vector correctly
    state_pars <- parTab[parTab$local==local,"names"]
    if(number >= 1) state_pars <- paste(state_pars,".",number,sep="")
    state_pars <- pars[state_pars]
    names(state_pars) <- parTab[parTab$local==local,"names"]
    all_pars <- pars[parTab[parTab$local=="all","names"]]
    pars <- c(state_pars, all_pars)

    ## Generate actual incidence data for these parameters
    y0s <- generate_y0s(pars["N_H"],pars["density"])
    y <- solveModelSimple_lsoda(ts, y0s, pars,TRUE)
    
    ## Generate predicted microcephaly incidence for these parameters - need to restrict to predictions within the data range
    probs <- generate_micro_curve(pars)
    probM <- generate_probM(y[,"I_M"], pars["N_H"], probs, pars["b"], pars["p_MH"], pars["baselineProb"], 1)
    probM <- probM[which(y[,"time"] >= min(dat[,"startDay"]) & y[,"time"] <= max(dat[,"endDay"]))]

    probM <- average_buckets(probM, dat[,"buckets"])

    ## Generate predicted microcephaly cases or per birth incidence depending on what was provided
    if(!is.null(dat$births)){
        predicted <- probM*dat[,"births"]*pars["propn"]
    } else {
        predicted <- probM*pars["propn"]
    }
    predicted <- data.frame(day=dat[,"meanDay"],number=predicted)
    
    ## Generate predicted incidence cases
    N_H <- average_buckets(rowSums(y[,5:8]), incDat$buckets)
    
    tmpY <- y[which(y[,"time"] >= min(incDat[,"startDay"]) & y[,"time"] <= max(incDat[,"endDay"])),]
    inc <- diff(tmpY[,"incidence"])
                                     
    ## At the moment this really needs to be in weeks
    inc <- sum_buckets(inc, incDat$buckets)

    perCapInc <- (1-(1-(inc/N_H))*(1-pars["baselineInc"]))*pars["incPropn"]
    y <- data.frame(time=incDat$meanDay,inc=perCapInc)
    
    return(list("micro"=predicted,"inc"=y))
}

#' Combine MCMC chains
#'
#' Given a list of MCMC chains, rbinds them all to give one big chain
#' @param chains the list of chains
#' @return a single data frame
#' @export
combine_chains <- function(chains){
    return(do.call("rbind",chains))
}
