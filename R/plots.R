
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
plot_setup_data <- function(chain, dat=NULL, incDat=NULL, parTab, ts, local, runs=NULL,
                            startDay=0, noMonths=12,noWeeks=52,perCap=FALSE,forecast=FALSE,
                            forecastPostChange=FALSE){
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
    if(!is.null(incDat[incDat$local == local,])){
        if(!(nrow(incDat[incDat$local == local,]) == 0)){
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
    ## If we want to assume that no behaviour changed then we need to manually set these parameters
    ## for forecasting method
    if(forecastPostChange){
      bestPars["incPropn2"] <- bestPars["incPropn"]
      bestPars["abortion_rate"] <- bestPars["avoided_births"] <- bestPars["birth_reduction"] <- 0
      bestPars["propn2"] <- bestPars["propn"]
    }
    f <- create_forecast_function(parTab,tmpDat,tmpInc,ts,!forecast)
    tmp <- f(bestPars,perCap)
    bestInc <- tmp$ZIKV
    bestMicro <- tmp$microCeph
    
    ## Get sample indices with which to plot prediction intervals
    samples <- sample(nrow(chain),runs)

    allInc <- NULL
    allMicro <- NULL
    ## For each sample, get the sample parameters from the chain and calculate the trajectory
    for(i in samples){
        pars <- get_index_pars(chain,i)
        ## If we want to assume that no behaviour changed then we need to manually set these parameters
        ## for forecasting method
        if(forecastPostChange){
          pars["incPropn2"] <- pars["incPropn"]
          pars["abortion_rate"] <- pars["avoided_births"] <- pars["birth_reduction"] <- 0
          pars["propn2"] <- pars["propn"]
        }
        tmp <- f(pars,perCap)
        inc <- tmp$ZIKV
        micro <- tmp$microCeph
        allInc <- rbind(allInc, inc)
        allMicro <- rbind(allMicro, micro)
    }
    
    #incBounds <- as.data.frame(reshape2::melt(sapply(unique(allInc$time),function(x) quantile(allInc[allInc$time==x,"inc"],c(0.025,0.5,0.975)))[c(1,3),]))
    #microBounds <- as.data.frame(reshape2::melt(sapply(unique(allMicro$day),function(x) quantile(allMicro[allMicro$day==x,"number"],c(0.025,0.5,0.975)))[c(1,3),]))
    incBounds <- as.data.frame(t(sapply(unique(allInc$time),function(x) quantile(allInc[allInc$time==x,"inc"],c(0.025,0.5,0.975)))[c(1,3),]))
    microBounds <- as.data.frame(t(sapply(unique(allMicro$time),function(x) quantile(allMicro[allMicro$time==x,"microCeph"],c(0.025,0.5,0.975)))[c(1,3),]))
    incBounds <- cbind(incBounds, bestInc$inc)
    microBounds <- cbind(microBounds, bestMicro$microCeph)
    colnames(microBounds) <- c("lower","upper","best")
    colnames(incBounds) <- c("lower","upper","best")
    
    incBounds$x <- bestInc$time
    microBounds$x <- bestMicro$time
    incBounds$local <- local
    microBounds$local <- local
    #colnames(microBounds) <- c("quantile","time","micro")
    #colnames(incBounds) <- c("quantile","time","inc")
    
                                        # microBounds[,"time"] <- tmpDat[,"meanDay"][microBounds[,"time"]]
                                        # incBounds[,"time"] <- tmpInc[,"meanDay"][incBounds[,"time"]]
    
    return(list("incBounds"=incBounds,"microBounds"=microBounds, "data"=tmpDat,"incDat"=tmpInc))
}


#' Microcephaly parameter density plot
#'
#' Given a melted MCMC chain with a column for variable name, plots the posterior density with shaded regions for 95% CI.
#' @param varName the name of the variable to plot
#' @param dat the melted MCMC chain
#' @param xlabtext optional x-axis title
#' @param weeks if set, this divides everything by this amount. ie. if we want value in weeks rather than days
#' @param xlim optional vector of lower and upper x axis limits
#' @param xbreaks optional vector of x axis break points
#' @return a ggplot object density plot
#' @export
density_plot_func <- function(varName, dat,xlabtext="",weeks=1,xlim=NULL,xbreaks=NULL){
  ## Mode plot - note that this isn't the same as the calculated
    ## "peak week"
    var <- varName
    x <- dat[,"value"]/weeks
    q5 <- quantile(x,.025)
    q95 <- quantile(x,.975)
    medx <- median(x)
    x.dens <- density(x)
    df.dens <- data.frame(x = x.dens$x, y = x.dens$y)
    
    p1 <- ggplot(data =dat) + 
        theme_bw() + 
        geom_density(aes(x=value/weeks, y = ..density..), color = 'black',alpha=0.3,fill="blue")+ 
        geom_area(data = subset(df.dens, x >= q5 & x <= q95), aes(x=x,y=y), fill = 'blue',alpha=0.5) +
        geom_vline(xintercept = medx)+
        theme(axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),
              axis.text.x=element_text(size=8),axis.title.x=element_text(size=10),panel.grid.minor=element_blank()) +
        xlab(xlabtext)
    
    if(!is.null(xlim)){
        p1 <- p1 + coord_cartesian(xlim=xlim)+
            scale_x_continuous(breaks=xbreaks)
    }
    p1
}


#' Posterior density plot
#'
#' Plots posterior density of a given variable, facetted by state
#' @param meltedChain the melted MCMC chain
#' @param variableName the variable to plot
#' @param xlimits optional vector of lower and upper x axis range
#' @param title the title to give to the plot
#' @param breaks option x axis breaks
#' @param optional x axis labels. Must match breaks length
#' @return a ggplot object of posterior density
#' @export
density_plot_state <- function(meltedChain, variableName, xlimits=NULL, title="",breaks=NULL,labels=NULL){
    tmpDat <- meltedChain[meltedChain$variable==variableName,]
    plot <- ggplot(tmpDat,aes(x=value))+
        stat_density(aes(ymax=..density..,ymin=-..density..,fill=state,color=state),trim=TRUE,geom="ribbon",position="identity")+
        facet_grid(state~.)+

    ylab("") +
    ggtitle(title)+
    xlab("Value")+
    theme_bw() +
    theme(
        plot.margin=unit(c(0.1,0.5,0.1,0),"cm"),
        legend.position="none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
    if(!is.null(xlimits)){
        if(!is.null(breaks) & !is.null(labels)) plot <- plot + scale_x_continuous(limits=xlimits,breaks=breaks,labels=labels)+
                                                    theme(plot.margin=unit(c(0.1,0.5,0.1,0),"cm"),
                                                          legend.position="none",
                                                          axis.text.y=element_blank(),
                                                          axis.ticks.y=element_blank(),
                                                          axis.text.x=element_text(angle=45,hjust=1))
        else plot <- plot +
                 scale_x_continuous(limits=xlimits)
    }
    return(plot)
}


#' Generate date labels
#'
#' Given a data frame of microcephaly data, generates a vector of x-labels that represent dates
#' @param birth_dat the data frame of microcephaly data
#' @return the vector of labels
#' @export
generate_plot_date_labels <- function(birth_dat){
    x <- (unique(birth_dat$startDay) + unique(birth_dat$endDay))/2 - 365
    months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    labels <- unlist(lapply(c("2014","2015","2016"), function(x) paste(months,x,sep=" ")))
    return(labels)
}


#' @export
generate_country_boundaries <- function(country="Brazil"){
    brazil1 <- raster::getData("GADM",country=country,level=1)
    map <- ggplot2::fortify(brazil1)
    
    map$id <- as.integer(map$id)
    dat <- data.frame(id=1:(length(brazil1@data$NAME_1)),state=brazil1@data$NAME_1)
    map_df <- dplyr::inner_join(map,dat,by="id")
    map_df$state <- get_states()[as.character(map_df$state)]
    return(map_df)
}

#' @export
generate_country_plot <- function(var,map_df,centers, title){
  brazil1 <- raster::getData("GADM",country="Brazil",level=1)
  map <- ggplot2::fortify(brazil1)
  
  map$id <- as.integer(map$id)
                                        # centers[,var] <- signif(centers[,var],3)
  plot <- ggplot() + 
    geom_map(data=map_df,map=map,aes_string(x="long",y="lat",map_id="id",group="group",fill=var),col="black")+
    ggtitle(title)+
    geom_text(data = centers, aes_string(label = var, x = "x", y = "y"), size = 3)+
    mytheme + theme(plot.title=element_text(size=14),legend.text=element_text(size=16))
}

#' @export
generate_centroids <- function(country="Brazil",map_df,plotLabels){
  brazil1 <- raster::getData("GADM",country=country,level=1)
  ggplot2::map <- fortify(brazil1)
  
  map$id <- as.integer(map$id)
  dat <- data.frame(id=1:(length(brazil1@data$NAME_1)),state=brazil1@data$NAME_1)
  dat$state <- get_states()[as.character(dat$state)]
  
  # Finding centroids of each state as location for labels
  centers <- data.frame(gCentroid(brazil1, byid = TRUE))
  centers$state <- dat$state
  centers <- join(centers,plotLabels,by="state")
  return(centers)
}





#' Plot microcephaly curves with bounds
#'
#' Given an MCMC chain, plots microcephaly risk curves with 95% prediction intervals
#' @param chain the MCMC chain
#' @return the ggplot object
#' @export
plot_micro_bounds <- function(chain, samp_no=10000,ylim=0.3,chain2=NULL, colombia_chain=NULL){
    samples <- sample(nrow(chain), samp_no)
    microCurves <- matrix(nrow=samp_no,ncol=40*7)
    trimesterCurves <- matrix(nrow=samp_no,ncol=40*7)
    tm1 <- tm2 <- tm3 <- numeric(samp_no)
    i <- 1
    for(samp in samples){
        pars <- chain[samp,]
         pars <- as.numeric(chain[samp,])
        names(pars) <- colnames(chain)
        probs <- generate_micro_curve(pars)/0.8
        tm1[i] <- mean(probs[1:14*7])
        tm2[i] <- mean(probs[(14*7 + 1):(28*7)])
        tm3[i] <- mean(probs[(28*7 + 1):length(probs)])
        #tmp <- average_buckets(probs, rep(7,40))
        microCurves[i,] <- probs
        i <- i + 1
    }

    if(!is.null(colombia_chain)){
        samples <- sample(nrow(colombia_chain),samp_no)
        microCurves1 <- matrix(nrow=samp_no,ncol=40*7)
        i <- 1
        for(samp in samples){
            pars <- as.numeric(colombia_chain[samp,])
            names(pars) <- colnames(colombia_chain)
            probs <- generate_micro_curve(pars)
            microCurves1[i,] <- probs
            i <- i + 1
        }
        colMeans <- colMeans(microCurves1)
        colLower <- apply(microCurves1,2,function(x) quantile(x, 0.025))
        colUpper <- apply(microCurves1,2,function(x) quantile(x,0.975))
        colDat <- data.frame(weeks=1:(40*7),colMeans,colUpper,colLower)
        colDat[,c("colMeans","colUpper","colLower")] <- colDat[,c("colMeans","colUpper","colLower")]/max(colUpper)
    }

    
    mean_tm1 <- mean(tm1)
    upper_tm1 <- quantile(tm1, 0.025)
    lower_tm1 <- quantile(tm1,0.975)

    mean_tm2 <- mean(tm2)
    upper_tm2 <- quantile(tm2, 0.025)
    lower_tm2 <- quantile(tm2,0.975)

    mean_tm3 <- mean(tm3)
    upper_tm3 <- quantile(tm3, 0.025)
    lower_tm3 <- quantile(tm3,0.975)
    
    tm_dat <- data.frame(x=c(0,14*7, 14*7, 28*7, 28*7, 40*7),
                         zero=rep(0,6),
                         means=c(mean_tm1, mean_tm1, mean_tm2, mean_tm2,  mean_tm3, mean_tm3),
                         lower=c(lower_tm1, lower_tm1, lower_tm2, lower_tm2, lower_tm3, lower_tm3),
                         upper=c(upper_tm1, upper_tm1, upper_tm2, upper_tm2, upper_tm3, upper_tm3))
    
    
    if(!is.null(chain2)){
        samples <- sample(nrow(chain2),samp_no)
        i <- 1
        for(samp in samples){
            pars <- as.numeric(chain2[samp,])
            names(pars) <- colnames(chain2)
            probs <- generate_micro_curve(pars)/0.8
            trimesterCurves[i,] <- probs
            i <- i + 1
        }
        trimMeans <- colMeans(trimesterCurves)
        trimUpper <- apply(trimesterCurves,2,function(x) quantile(x, 0.025))
        trimLower <- apply(trimesterCurves,2,function(x) quantile(x,0.975))
        trimDat <- data.frame(weeks=1:(40*7),trimMeans,trimUpper,trimLower)
        trimDat[,2:4] <- trimDat[,2:4]/max(trimUpper)
    }
        
    
    means <- colMeans(microCurves)
    lower <- apply(microCurves,2,function(x) quantile(x, 0.025))
    upper <- apply(microCurves,2,function(x) quantile(x,0.975))
    dat <- data.frame(weeks=1:(40*7),means,lower,upper)
    dat[,2:4] <- dat[,2:4]/max(upper)

    p1 <- ggplot(dat)
    
    if(!is.null(chain2)){
        barDat <- data.frame(x=c(7,21,35),ymax=c(upper_tm1,upper_tm2,upper_tm3),ymin=c(lower_tm1,lower_tm2,lower_tm3))
        ##       print(barDat)
        ##        barDat <- data.frame(x=c(7,21,35),ymax=unique(trimUpper),ymin=unique(trimLower))
        ##      print(barDat)
        p1 <- p1 +
            geom_errorbar(data=barDat,aes(x=c(7*7,21*7,35*7),ymax=ymax,ymin=ymin),width=40) +
            ##geom_ribbon(data=trimDat,aes(x=weeks,ymax=trimMeans,ymin=rep(0,nrow(trimDat))),fill="blue",alpha=0.1,col="black")
            geom_ribbon(data=tm_dat,aes(x=x,ymax=means,ymin=zero),fill="blue",alpha=0.5,col="black")
        
        
    }
    if(!is.null(colombia_chain)){
        p1 <- p1 +
            geom_ribbon(data=colDat,aes(x=weeks,ymax=colUpper,ymin=colLower),fill="green",alpha=0.3) +
            geom_line(data=colDat,aes(x=weeks,y=colMeans),col="black",lty="dashed")
        
    }
    p1 <- p1 + 
        geom_ribbon(data=dat,aes(x=weeks,ymax=upper,ymin=lower),fill="firebrick2",alpha=0.8) +
        geom_line(aes(x=weeks,y=means),col="black",lty="dashed") +
        ylab("Probability of microcephaly given ZIKV infection")+
        xlab("Week of gestation at time of infection") +
                                        #ggtitle("Microcephaly Risk Curve") +
        theme_bw() +
        theme(axis.text.x=element_text(size=10),
              panel.grid.minor=element_blank(),
              plot.title=element_text(size=10),
              axis.text.y=element_text(size=10),
              axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=10)) +
        scale_y_continuous(expand=c(0,0),limits=c(0,ylim)) +
        scale_x_continuous(expand=c(0,0),limits=c(0,40*7),breaks=seq(0,280,by=7*4), labels=seq(0,40,by=4))+
        geom_vline(xintercept=c(14*7, 28*7),col="grey",lty="dashed")
    

  
    return(p1)

}

    
# PDF - Rich's function to print to device without potential for bad errors
to.pdf <- function(expr, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}
#' PNG - Rich's function to print to device without potential for bad errors
#'
#' Prints to png, but turns dev off if fails
#' @param expr expression to give plot
#' @param filename filename to print to
#' @export
to.png <- function(expr, filename, ..., verbose=TRUE) {
    if ( verbose )
        cat(sprintf("Creating %s\n", filename))
    png(filename, ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
}

#' SVG - Rich's function to print to device without potential for bad errors
#'
#' Prints to SVG, but turns dev off if fails
#' @param expr expression to give plot
#' @param filename filename to print to
#' @export
to.svg <- function(expr, filename, ..., verbose=TRUE) {
    if ( verbose )
        cat(sprintf("Creating %s\n", filename))
    svg(filename, ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
}


#' Plot MCMC chains
#'
#' Given a list of MCMC chains and optionally the parameter table, plots the MCMC chains and densities using coda
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
    probs <- average_buckets(probs,rep(7,40))
    bestProbs <- data.frame(prob=probs,week=seq(0,279/7,by=1))
    myPlot <- ggplot() +
        geom_line(data=bestProbs, aes_string(x="week",y="prob"),col="blue",lwd=1) +
        ylab("Probability of microcephaly given infection")+
        xlab("Week of gestation at infection") +
        ggtitle("Microcephaly Risk Curve") +
        theme_bw()+
        theme(axis.text.x=element_text(size=10),
              panel.grid.minor=element_blank(),
              plot.title=element_text(size=10),
              axis.text.y=element_text(size=10),
              axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=10))
   
    ## Add some random lines
    samples <- sample(nrow(chain),runs)
    index <- 1
    allProbs <- NULL
    for(i in samples){
        tmpPars <- as.numeric(chain[i,])
        names(tmpPars) <- colnames(chain)
        tmpPars["tstep"] <- 1
        probs <- generate_micro_curve(tmpPars)
        probs <- average_buckets(probs,rep(7,40))
        allProbs[[index]] <- probs <- data.frame(prob=probs,week=seq(0,279/7,by=1))
        myPlot <- myPlot + geom_line(data=probs, aes_string(x="week",y="prob"),alpha=0.3,lwd=0.5,col="red")
        index <- index+1
    }
  return(myPlot)
}


#' Format MCMC chain
#'
#' Formats an MCMC chain to be plotted in Figure 1
#' @export
make_micro_dat <- function(chain, samp_no,scale=FALSE){
  samples <- sample(nrow(chain), samp_no)
  microCurves <- matrix(nrow=samp_no,ncol=40*7)
  trimesterCurves <- matrix(nrow=samp_no,ncol=40*7)
  tm1 <- tm2 <- tm3 <- numeric(samp_no)
  i <- 1
  for(samp in samples){
    pars <- chain[samp,]
    pars <- as.numeric(chain[samp,])
    names(pars) <- colnames(chain)
    probs <- generate_micro_curve(pars)
    probs <- probs
    tm1[i] <- mean(probs[1:14*7])
    tm2[i] <- mean(probs[(14*7 + 1):(28*7)])
    tm3[i] <- mean(probs[(28*7 + 1):length(probs)])
    microCurves[i,] <- probs
    i <- i + 1
  }
  #microCurves <- microCurves/max(microCurves)
   means <- colMeans(microCurves)
  lower <- apply(microCurves,2,function(x) quantile(x, 0.025))
  upper <- apply(microCurves,2,function(x) quantile(x,0.975))
  if(scale){
    means <- means/max(upper)
    lower <- lower/max(upper)
    upper <- upper/max(upper)
  }
  dat <- data.frame(weeks=1:(40*7),means,upper,lower)
  return(dat)
}
