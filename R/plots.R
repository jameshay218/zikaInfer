#' Plot MCMC chains
#'
#' Given an MCMC chain and optionally the parameter table, plots the MCMC chains and densities using coda
#' @param chains list of MCMC chains
#' @param parTab optional parameter table
#' @return nothing - saves plots to pdf
#' @export
#' @useDynLib zikaProj
plot_MCMC <- function(chains, parTab=NULL){
    for(i in 1:length(chains)){
        if(!is.null(parTab)){
            indices <- which(parTab$fixed==0)
            chains[[i]] <- chains[[i]][,c(indices+1, ncol(chains[[i]]))]
        }
        chains[[i]] <- as.mcmc(chains[[i]])
    }

    pdf("mcmcplots.pdf")
    plot(as.mcmc.list(chains))
    dev.off()
}

#' Plot microceph curves
#'
#' Given an MCMC chain, plots the best fitting microcephaly risk curve and some faint, random draws from the posterior
#' @param chain The MCMC chain. Must have the gamma distribution mean, variance and scale
#' @param runs number of random draws to plot
#' @return a ggplot object with the microcephaly risk curve
#' @export
#' @useDynLib zikaProj
plot_random_microceph_curves <- function(chain, runs){
    ## Get best fitting parameters
    bestPars <- as.numeric(chain[which.max(chain[,"lnlike"]),])
   
    names(bestPars) <- colnames(chain)
    bestPars["tstep"] <- 1

    probs <- generate_micro_curve(bestPars)
 
    bestProbs <- data.frame(prob=probs,week=seq(0,39,by=1))

    myPlot <- ggplot() +
        geom_line(dat=bestProbs, aes(x=week,y=prob),col="blue",lwd=1) +
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
        allProbs[[index]] <- probs <- data.frame(prob=probs,week=seq(0,39,by=1))
        myPlot <- myPlot + geom_line(dat=probs, aes(x=week,y=prob),alpha=0.3,lwd=0.5,col="red")
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
#' @useDynLib zikaProj
get_state_name <- function(state){
    country_names <- c(
        "pernambuco"="Pernambuco",
        "bahia"="Bahia",
        "saopaulo"="Sao Paulo",
        "paraiba"="Paraiba",
        "maranhao"="Maranhao",
        "ceara"="Ceara",
        "sergipe"="Sergipe",
        "riodejaneiro"="Rio de Janeiro",
        "piaui"="Piaui",
        "riograndedonorte"="Rio Grande Norte",
        "minasgerais"="Minas Gerais",
        "matogrosso"="Mato Grosso",
        "alagoas"="Algoas",
        "para"="Para",
        "acre"="Acre",
        "espiritosanto"="Espirito Santo",
        "goias"="Goias",
        "tocantins"="Tocantins"
    )
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
#' @param t_pars the time parameters used for solving the ODEs
#' @param runs the number of runs to use for prediction intervals
#' @param incDat optional data frame of actual incidence data
#' @param mcmcPars a named vector of the burnin, adaptive period and thinning
#' @param ylimM ylimit for the microcephaly plot
#' @param ylimI ylimit for the incidence plot
#
#' @return a ggplot object with the incidence plots
#' @export
#' @useDynLib zikaProj
plot_best_trajectory_multi <- function(chain, realDat, parTab, t_pars, runs=100, incDat=NULL, mcmcPars=c("burnin"=50000,"adaptive"=100000,"thin"=50), ylimM=200,ylimI=0.01){
    ps <- NULL
    
    states <- unique(parTab$local)
    states <- states[states != "all"]

    for(i in 1:length(states)){
        ps[[i]] <- plot_best_trajectory_single(states[i], chain, realDat, parTab, t_pars, runs, incDat=incDat, ylabel=FALSE, xlabel=FALSE, mcmcPars,ylimM,ylimI)
    }
    ncols <- ceiling(length(states)/4)
    allPlot <- do.call("grid.arrange",c(ps,ncol=ncols))
    return(allPlot)
}

generate_x_labels <- function(startDay, endDay){
    days <- getDaysPerMonth(4)
    buckets <- rep(days, 4)
    xpositions <- cumsum(buckets) - buckets
    indices <- which(xpositions >= startDay & xpositions <= endDay)
    
    labels <- c("01/2013","04/2013","07/2013","10/2013","01/2014","04/2014","07/2014","10/2014","01/2015","04/2015","07/2015","10/2015","01/2016","04/2016","07/2016","10/2016")
    return(list("labels"=labels[indices],"positions"=xpositions[indices]))
}


#' Plot single best trajectory
#'
#' For a given run name and state, plots the best zika and microcephaly incidence plot
#' @param runName the type of model fitting done
#' @param state the name of the state to be fitted
#' @param runs the number of runs to use for prediction intervals
#' @param ylabel boolean whether or not to add y axis label
#' @param xlabel boolean whether or not to add x axis label
#' @param mcmcPars a named vector of the burnin, adaptive period and thinning
#' @param ylimM ylimit for the microcephaly plot
#' @param ylimI ylimit for the incidence plot
#' @param incDat optional data frame of actual incidence data
#' @return a ggplot object with the incidence plots
#' @export
#' @useDynLib zikaProj
plot_best_trajectory_single <- function(local, chain=NULL, realDat=NULL, parTab=NULL, t_pars=NULL, runs=100,incDat=NULL, ylabel=TRUE,xlabel=TRUE, mcmcPars=c("burnin"=50000,"adaptive_period"=100000,"thin"=50),ylimM=200, ylimI=0.02){
    allDat <- plot_setup_data(chain,realDat, parTab, t_pars,local,runs)

    bestMicro <- allDat$bestMicro
    bestInc <- allDat$bestInc
    incBounds <- allDat$incBounds
    microBounds <- allDat$microBounds
    dat <- allDat$data

    xlim <- c(min(dat[,"startDay"]),max(dat[,"endDay"]))

    quantiles <- unique(microBounds[,"quantile"])
    botM <- microBounds[microBounds[,"quantile"]==quantiles[1],c("time","micro")]
    topM <- microBounds[microBounds[,"quantile"]==quantiles[2],c("time","micro")]

    botI <- incBounds[incBounds[,"quantile"]==quantiles[1],c("time","inc")]
    topI <-incBounds[incBounds[,"quantile"]==quantiles[2],c("time","inc")]
    polygonM <- create_polygons(botM, topM)
    polygonI <- create_polygons(botI, topI)

    xlabs <- generate_x_labels(xlim[1],xlim[2])
    myPlot <- microceph_plot(dat,microBounds,bestMicro,polygonM,local,xlim,ylimM,xlabs)
    incPlot <- inc_plot(incBounds,bestInc,polygonI,ylimI,incDat)
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
    
    myPlot <- ggplot() + geom_point(dat=dat,aes(y=microCeph,x=meanDay),col="black",size=2) +
        geom_line(dat=microBounds,aes(y=micro,x=time,group=quantile),lwd=0.5,linetype=2,col="blue",alpha=0.5) +
        geom_line(dat=bestMicro,aes(y=number,x=day),col="blue",lwd=0.5) +
        geom_polygon(data=polygonM,aes(x=x,y=y),alpha=0.2,fill="blue")+
        scale_y_continuous(limits=c(0,ylim))+
        scale_x_continuous(labels=xlabels,breaks=xlabBreaks,limits=xlim)+
        theme_bw() +
        ylab("Microcephaly cases (blue)") +
        xlab("Date") +
        ggtitle(get_state_name(local)) + 
        theme(
            panel.grid.minor=element_blank(),
            plot.title=element_text(size=10),
            axis.text.x=element_text(size=6,hjust=1,angle=45),
            axis.text.y=element_text(size=6),
            axis.title.x=element_text(size=8),
            axis.title.y=element_text(size=8),
            plot.margin=unit(c(0.1,0.8,0.1,0.1),"cm")
        )
    return(myPlot)
}

inc_plot <- function(incBounds, bestInc, polygonI, ylimI,incDat=NULL){
    incPlot <- ggplot() +
        geom_line(dat=incBounds,aes(y=inc,x=time,group=quantile),linetype=2,col="red",size=0.5,alpha=0.5) +
        geom_line(dat=bestInc,aes(y=I_H,x=times),col="red",lwd=0.5)+
        geom_polygon(data=polygonI,aes(x=x,y=y),alpha=0.2,fill="red")+
        scale_y_continuous(limits=c(0,ylimI))+
        ylab("\nPer capita incidence (red)")+
        xlab("")+
        theme(
            panel.background=element_rect(fill=NA),
            panel.grid=element_blank(),
            axis.line.y = element_line(colour="black"),
            axis.text.y=element_text(size=6,colour="black"),
            axis.title.y=element_text(size=8,angle=-90),
            axis.text.x=element_text(size=6),
            axis.title.x=element_text(size=8),
            plot.margin=unit(c(0.1,0.8,0.1,0.1),"cm")
            
        )
    if(!is.null(incDat)){
        incDat$meanDay <- rowMeans(cbind(incDat[,"startDay"],incDat[,"endDay"]))
        incDat$perCapInc <- incDat[,"inc"]/incDat[,"N_H"]
        incPlot <- incPlot + geom_point(dat=incDat,aes(x=meanDay,y=perCapInc),col="black",size=1, shape=3)
    }
    return(incPlot)
}


overlapPlots <- function(p1,p2, centre_plot=TRUE){
                                        # extract gtable
  g1 <- ggplot_gtable(ggplot_build(p1))
  g2 <- ggplot_gtable(ggplot_build(p2))
  
  # overlap the panel of 2nd plot on that of 1st plot
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
#' @param dat the data frame of real data
#' @param parTab the parameter table
#' @param t_pars time parameters
#' @param local string for state under consideration
#' @param number the number to be appended to parameter names in the MCMC chain
#' @export
#' @useDynLib zikaProj
plot_setup_data <- function(chain, dat, parTab, t_pars, local, runs=NULL){
    ## Get the subset of data for this particular state and find the middle day for each month
    tmpDat <- dat[dat$local==local,]
    tmpDat$meanDay <- rowMeans(cbind(tmpDat[,"startDay"],tmpDat[,"endDay"]))
    
    ## Get the set of best fitting parameters. These will be used to generate the best fit incidence and
    ## microcephaly line
    bestPars <- get_best_pars(chain)
    tmp <- plot_setup_data_auxiliary(bestPars,tmpDat,parTab,t_pars,local)
    bestInc <- tmp$inc
    bestMicro <- tmp$micro

    ## Get sample indices with which to plot prediction intervals
    samples <- sample(nrow(chain),runs)

    allInc <- NULL
    allMicro <- NULL
    ## For each sample, get the sample parameters from the chain and calculate the trajectory
    for(i in samples){
        pars <- get_index_pars(chain,i)
        tmp <- plot_setup_data_auxiliary(pars,tmpDat,parTab,t_pars,local)
        inc <- tmp$inc
        micro <- tmp$micro
        allInc <- rbind(allInc, inc)
        allMicro <- rbind(allMicro, micro)
    }

    incBounds <- as.data.frame(melt(sapply(unique(allInc$times),function(x) quantile(allInc[allInc$times==x,"I_H"],c(0.025,0.5,0.975)))[c(1,3),]))
    microBounds <- as.data.frame(melt(sapply(unique(allMicro$day),function(x) quantile(allMicro[allMicro$day==x,"number"],c(0.025,0.5,0.975)))[c(1,3),]))
    colnames(microBounds) <- c("quantile","time","micro")
    colnames(incBounds) <- c("quantile","time","inc")
    
    microBounds[,"time"] <- tmpDat[,"meanDay"][microBounds[,"time"]]

    return(list("bestInc"=bestInc,"incBounds"=incBounds,"bestMicro"=bestMicro,"microBounds"=microBounds, "data"=tmpDat))  
}

plot_setup_data_auxiliary <- function(pars, dat, parTab, t_pars, local){
    ## As the model has many parameters with the same name, need to find the index in the MCMC colnames
    ## that corresponds to this state
    dat$meanDay <- rowMeans(cbind(dat[,"startDay"],dat[,"endDay"]))

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
    y <- solveModelSimple(t_pars, y0s, pars)
   
    ## Generate predicted microcephaly incidence for these parameters - need to restrict to predictions within the data range
    probs <- generate_micro_curve(pars)
    probM <- generate_probM(y[,"I_M"],probs, pars["N_H"], pars["b"], pars["p_MH"], pars["baselineProb"], 1)
    probM <- probM[which(y[,"times"] >= min(dat[,"startDay"]) & y[,"times"] <= max(dat[,"endDay"]))]
    probM <- average_buckets(probM, dat[,"buckets"])

    predicted <- probM*dat[,"births"]*pars["propn"]
    predicted <- data.frame(day=dat[,"meanDay"],number=predicted)
    
    y[,"I_H"] <- pars["incPropn"]*y[,"I_H"]/rowSums(y[,c("I_H","E_H","S_H","R_H")])
    y[,"I_H"] <- 1- (1-(y[,"I_H"]))*(1-pars["baselineInc"])
    y <- as.data.frame(y)
    y$times <- as.integer(y$times)
    
    return(list("micro"=predicted,"inc"=y))
}

#' Combine MCMC chains
#'
#' Given a list of MCMC chains, rbinds them all to give one big chain
#' @param chains the list of chains
#' @return a single data frame
#' @export
#' @useDynLib zikaProj
combine_chains <- function(chains){
    return(do.call("rbind",chains))
}

#' Volcano plots
#'
#' Given a list of MCMC chains, produces a grid of volcano plots for the unfixed state-specific parameters
#' @param chains the list of chains
#' @param parTab the parameter table used for correct indexing
#' @return a ggplot table
#' @export
#' @useDynLib
volcano_plots <- function(chains, parTab){
    allChains <- combine_chains(chains)

    states <- unique(parTab$local)
    states <- states[states != "all"]
    tmp <- sort(states,index.return=TRUE)
    numbers <- tmp$ix - 1
    states <- tmp$x

    country_names <- c(
        "pernambuco"="Pernambuco",
        "bahia"="Bahia",
        "saopaulo"="Sao Paulo",
        "paraiba"="Paraiba",
        "maranhao"="Maranhao",
        "ceara"="Ceara",
        "sergipe"="Sergipe",
        "riodejaneiro"="Rio de Janeiro",
        "piaui"="Piaui",
        "riograndedonorte"="Rio Grande Norte",
        "minasgerais"="Minas Gerais",
        "matogrosso"="Mato Grosso",
        "alagoas"="Algoas",
        "para"="Para",
        "acre"="Acre",
        "espiritosanto"="Espirito Santo",
        "goias"="Goias",
        "tocantins"="Tocantins"
    )

    allR0 <- allAttack <- allEpi <- allDensity <- allPropn <- NULL
    
    for(i in 1:length(states)){
        local <- states[i]
                                        #name <- get_state_name(local)
        name <- as.character(states[i])
        unfixed_pars <- parTab[parTab$local==local,"names"]
        number <- numbers[i]
        if(number >= 1) unfixed_pars <- paste(unfixed_pars,".",number,sep="")
        unfixed_pars <- allChains[,unfixed_pars]
        names(unfixed_pars) <- parTab[parTab$local==local,"names"]
        
        fixed_pars <- allChains[,parTab[parTab$local=="all","names"]]
        pars <- cbind(fixed_pars, unfixed_pars, row.names=NULL)
        pars <- pars[sample(seq(1,nrow(pars),by=1),1000,replace=FALSE),]
        
        R0s <- data.frame("value"=r0.vector(pars),"dat"=name, row.names=NULL)
        attack <- data.frame("value"=sapply(R0s$value, function(x) nleqslv(0.8, zikaProj::simeq, R0=x)$x),"dat"=name, row.names=NULL)
        tmpEpi <- data.frame("value"=pars$epiStart,"dat"=name, row.names=NULL)
        tmpDensity <- data.frame("value"=pars$density,"dat"=name, row.names=NULL)
        tmpPropn <- data.frame("value"=pars$propn,"dat"=name, row.names=NULL)

        allR0 <- rbind(allR0, R0s,row.names=NULL)
        allAttack <- rbind(allAttack, attack,row.names=NULL)
        allEpi <- rbind(allEpi, tmpEpi,row.names=NULL)
        allDensity <- rbind(allDensity, tmpDensity,row.names=NULL)
        allPropn <- rbind(allPropn, tmpPropn,row.names=NULL)
    }
    
    mytheme <- theme(
        panel.grid.minor=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour="gray20",linetype="dashed",size=0.25),
        text=element_text(size=10,colour="gray20"),
        axis.line=element_line(colour="gray20"),
        axis.line.y=element_line(colour="gray20"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y=element_text(size=10,angle=0),
        strip.text = element_blank(),
        legend.position="none",
        panel.margin=unit(0,"lines"),
        plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"),
        panel.grid = element_line(colour="gray20"),
        panel.background=element_rect(colour="gray40"),
        strip.background=element_rect(colour="black")
    )

    R0plot <- ggplot(allR0,aes(x=value))+
        stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
                                        #coord_cartesian(xlim=c(0,10))+
        scale_x_continuous(limits=c(0,10),breaks=seq(0,10,by=1))+
        facet_grid(dat~.,labeller=as_labeller(country_names))+
                                        #coord_flip()+ 
        ylab("") +
        ggtitle("Basic Reproductive Number, R0")+
        xlab("Value")+
        mytheme +
        theme(
            plot.margin=unit(c(0.1,0.5,0.1,0),"cm")
        )
    
    densityPlot <- ggplot(allDensity,aes(x=value))+
        stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
        scale_x_continuous(limits=c(0,15),breaks=seq(0,15,by=3))+
        facet_grid(dat~.)+
        ylab("") +
        ggtitle("Mosquito Density per Person")+
        xlab("Value")+
        mytheme +
        theme(
            plot.margin=unit(c(0.1,0.5,0.1,0),"cm")
        )

    xlabels <- c("01/2014","04/2014","07/2014","10/2014","01/2015","04/2015","07/2015","10/2015","01/2016")
    days <- getDaysPerMonth()
    years <- cumsum(c(0,days,days,days[1:2]))
    years <- years[seq(1,length(years),by=3)]
    epiplot <- ggplot(allEpi,aes(x=value))+
        stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
        facet_grid(dat~.)+
        scale_x_continuous(breaks=years,labels=xlabels)+
        ylab("") +
        ggtitle("Epidemic Start Time")+
        xlab("Value")+
        mytheme

    propnplot <- ggplot(allPropn,aes(x=value))+
        stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
      #  scale_x_continuous(limits=c(0,1))+
        facet_grid(dat~.)+
        ylab("") +
        ggtitle("Proportion of affected births")+
        xlab("Value")+
        mytheme+
        theme(
            plot.margin=unit(c(0.1,0.5,0.1,0),"cm")
  )


    #return(list(R0plot,densityPlot,epiplot, propnplot))
    strips <- R0plot + theme(strip.text=element_text())
    tmp <- ggplot_gtable(ggplot_build(strips))
    tmp1 <- gtable_filter(tmp,"strip-right")
    
g3 <- grid.arrange(arrangeGrob(textGrob(label="State",gp=gpar(fontsize=16,colour="grey20")),
                               tmp1,rectGrob(gp=gpar(col="white")),heights=c(0.55,7,0.8)),
                   R0plot,densityPlot,arrangeGrob(textGrob(label="State",gp=gpar(fontsize=16,colour="grey20")),
                                                 tmp1,rectGrob(gp=gpar(col="white")),heights=c(0.55,7,0.8)),
                   epiplot,propnplot,ncol=3,widths=c(1.7,6,6))
    return(g3)
}


