#' Auxiliary setup
#'
#' Given the run name, sets up all of the parameters needed to produce all plots
#' @param runName the name of the runto be plotted. This should match the filename of the MCMC chains
#' @param runDir the full directory path where the MCMC chains are located
#' @param allDatFile local filename of actual data OPTIONAL
#' @param nchain number of chains run
#' @param mcmcPars named vector (burnin, adaptive, thin) to format the chain correctly
#' @return a list of all the parameters needed for plotting
#' @export
#' @useDynLib zikaProj
plot_setup <- function(runName, runDir="~/net/home/zika/outputs3", allDatFile = "~/net/home/zika/allDat20.04.16.csv",nchains=3, mcmcPars=c("burnin"=50000,"adaptive"=100000,"thin"=50)){
    curWD <- getwd()
    setwd(runDir)
    realDat <- read.csv(allDatFile)
    parTab <- setupParTable(version=3,realDat=realDat)
    pars <- setupListPars(duration = 1200,N_H = 9000000,N_M=27000000,version=3)

    ## Get filenames of chains
    filenames <- NULL
    for(i in 1:nchains) filenames[[i]] <- paste(sprintf("%s_%d", runName, i),"B_chain.csv",sep="_")

    priors <- NULL
    peakTimes <- NULL
    
    if(runName == "simulation_single_unconstrained"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco")
        realDat <- testDat
    } else if(runName == "simulation_single_priors"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco")
        peakTimes <- data.frame("start"=400,"end"=500,"local"=as.character(unique(states[states != "all"])))
    } else if(runName == "simulation_multi_unconstrained"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia","saopaulo")
    } else if(runName=="simulation_multi_prior"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia","saopaulo")
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
    } else if(runName=="real_A_unconstrained"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia")
    } else if(runName=="real_A_prior"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia")
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
    } else if(runName=="real_A_propn"){
        states <- c("all","pernambuco","bahia")
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
        parTab[parTab$names=="propn","fixed"] <- 0
    } else if(runName=="real_B_unconstrained"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia","saopaulo")
    } else if(runName=="real_B_prior"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia","saopaulo")
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
    } else if(runName=="real_B_propn"){
        states <- c("all","pernambuco","bahia","saopaulo")
        parTab[parTab$names=="propn","fixed"] <- 0
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
    } else if(runName == "real_3_all"){
        states <- c("all","pernambuco", "bahia","saopaulo", "paraiba", "maranhao", "ceara","sergipe","riodejaneiro","piaui","riograndedonorte","minasgerais", "matogrosso", "alagoas", "para", "acre", "goias", "espiritosanto","tocantins")
        parTab[parTab$names=="propn","fixed"] <- 0
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
    }
    
    parTab <- parTab[parTab$local %in% states,]
    testDat <- generate_multiple_data(pars[[1]],parTab,NULL)
    if(grepl("simulation",runName)) realDat <- testDat
    places <- states[states != "all"]
    
    if(!is.null(peakTimes)) peakTimes$local <- as.character(peakTimes$local)

    realDat <- realDat[realDat$local %in% places,colnames(realDat) %in% colnames(testDat)]
    
    parTab[parTab$names=="b","fixed"] <- 1
    parTab[parTab$names=="mean","fixed"] <- 0
    parTab[parTab$names=="var","fixed"] <- 0
    parTab[parTab$names=="scale","fixed"] <- 0

    allChains <- NULL
    for(i in 1:length(filenames)){
        tmpChain <- read.csv(filenames[[i]])
        tmpChain <- tmpChain[tmpChain[,1] > (mcmcPars["adaptive"]+mcmcPars["burnin"]),]
        tmpChain <- cbind(tmpChain, i)
        allChains[[i]] <- tmpChain
    }
    setwd(curWD)
    return(list("chains"=allChains,"data"=realDat,"parTab"=parTab,"states"=places,"pars"=pars))
}

#' Plot MCMC chains
#'
#' Given an MCMC chain and optionally the parameter table, plots the MCMC chains and densities using coda
#' @param chains list of MCMC chains
#' @param parTab optional parameter table
#' @return nothing - saves plots to pdf
#' @export
#' @useDynLib zikaProj
plot_MCMC <- function(chain, parTab=NULL){
    for(i in 1:length(chain)){
        if(!is.null(parTab)){
            indices <- which(parTab$fixed==0)
            chain[[i]] <- chain[[i]][,c(indices+1, ncol(chain[[i]])-1)]
        }
        chain[[i]] <- as.mcmc(chain[[i]])
    }

    pdf("mcmcplots.pdf")
    plot(as.mcmc.list(chain))
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
    tmpChain <- chain[,c("mean","var","scale","lnlike")]

    ## Get best fitting parameters
    bestPars <- as.numeric(tmpChain[which.max(tmpChain[,"lnlike"]),])
    bestMean <- bestPars[1]
    bestVar <- bestPars[2]
    bestRate <- bestMean/bestVar
    bestShape <- bestMean/bestRate
    probs <- dgamma(0:39,bestShape,bestRate)*bestPars[3]
    probs[probs > 1] <- 1
    probs <- data.frame(prob=probs,week=seq(0,39,by=1))
    myPlot <- ggplot() +
        geom_line(dat=probs, aes(x=week,y=prob),col="blue",lwd=1) +
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
    samples <- sample(nrow(tmpChain),runs)
    index <- 1
    for(i in samples){
        tmpPars <- as.numeric(tmpChain[i,])
        gammaMean <- tmpPars[1]
        gammaVar <- tmpPars[2]
        
        rate <- gammaMean/gammaVar
        shape <- gammaMean*rate
        probs <- dgamma(0:39,shape,rate)*tmpPars[3]
        probs[probs > 1] <- 1
        probs <- data.frame(prob=probs,week=seq(0,39,by=1))
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
#' @param runName the type of model fitting done
#' @param state the name of the state to be fitted
#' @param runs the number of runs to use for prediction intervals
#' @param ylabel boolean whether or not to add y axis label
#' @param xlabel boolean whether or not to add x axis label
#' @param runDir the full directory path where the MCMC chains are located
#' @param allDatFile local filename of actual data OPTIONAL
#' @param nchain number of chains run
#' @return a ggplot object with the incidence plots
#' @export
#' @useDynLib zikaProj
plot_best_trajectory_multi <- function(runName, chain=NULL, realDat=NULL, parTab=NULL, pars=NULL, runs=100, runDir="~/net/home/zika/outputs3", allDatFile = "allDat20.04.16.csv",nchains=3, mcmcPars=c("burnin"=50000,"adaptive"=100000,"thin"=50)){
    ps <- NULL
    if(is.null(chain) | is.null(realDat) | is.null(parTab) | is.null(pars)){
        allPars <- plot_setup(runName, runDir, allDatFile, nchains, mcmcPars)
        chain <- allPars$chains
        parTab <- allPars$parTab
        realDat <- allPars$dat
        pars <- allPars$pars
    }
    
    states <- unique(parTab$local)
    states <- states[states != "all"]

    for(i in 1:length(states)){
        ps[[i]] <- plot_best_trajectory_single(runName, states[i], chain, realDat, parTab, pars, runs, ylabel=FALSE, xlabel=FALSE, runDir, allDatFile, nchains, mcmcPars)
    }
    allPlot <- do.call("grid.arrange",c(ps,ncol=4))
    return(allPlot)
}


#' Plot single best trajectory
#'
#' For a given run name and state, plots the best zika and microcephaly incidence plot
#' @param runName the type of model fitting done
#' @param state the name of the state to be fitted
#' @param runs the number of runs to use for prediction intervals
#' @param ylabel boolean whether or not to add y axis label
#' @param xlabel boolean whether or not to add x axis label
#' @param runDir the full directory path where the MCMC chains are located
#' @param allDatFile local filename of actual data OPTIONAL
#' @param nchain number of chains run
#' @return a ggplot object with the incidence plots
#' @export
#' @useDynLib zikaProj
plot_best_trajectory_single <- function(runName, local, chain=NULL, realDat=NULL, parTab=NULL, pars=NULL, runs=100,ylabel=TRUE,xlabel=TRUE, runDir="~/net/home/zika/outputs3", allDatFile = "allDat20.04.16.csv",nchains=3, mcmcPars=c("burnin"=50000,"adaptive"=100000,"thin"=50)){
    if(is.null(chain) | is.null(realDat) | is.null(parTab) | is.null(pars)){
        allPars <- plot_setup(runName, runDir, allDatFile, nchains, mcmcPars)
        chain <- allPars$chains
        parTab <- allPars$parTab
        realDat <- allPars$dat
        pars <- allPars$pars
    }

    allDat <- plot_setup_data(chain,realDat, parTab, pars[[1]],local,runs)

    bestMicro <- allDat$bestMicro
    bestInc <- allDat$bestInc
    incBounds <- allDat$incBounds
    microBounds <- allDat$microBounds
    dat <- allDat$data

    xlabels <- c("01/2014","04/2014","07/2014","10/2014","01/2015","04/2015","07/2015","10/2015","01/2016")
    xlabBreaks <- dat[,"meanDay"][seq(1,length(dat[,"meanDay"]),by=3)]

    quantiles <- unique(microBounds[,"quantile"])
    botM <- microBounds[microBounds[,"quantile"]==quantiles[1],c("time","micro")]
    topM <- microBounds[microBounds[,"quantile"]==quantiles[2],c("time","micro")]

    botI <- incBounds[incBounds[,"quantile"]==quantiles[1],c("time","inc")]
    topI <-incBounds[incBounds[,"quantile"]==quantiles[2],c("time","inc")]
    
    polygonM <- create_polygons(botM, topM)
    polygonI <- create_polygons(botI, topI)
    
    myPlot <- ggplot() + geom_point(dat=dat,aes(y=microCeph,x=meanDay),col="black",size=2) +
      geom_line(dat=microBounds,aes(y=micro,x=time,group=quantile),lwd=0.5,linetype=2,col="blue",alpha=0.5) +
      geom_line(dat=bestMicro,aes(y=number,x=day),col="blue",lwd=0.5) +
        geom_polygon(data=polygonM,aes(x=x,y=y),alpha=0.2,fill="blue")+
        scale_y_continuous(limits=c(0,200))+
        scale_x_continuous(labels=xlabels,breaks=xlabBreaks)+
        theme_bw() +
        ylab("Microcephaly cases (blue)") +
        xlab("Date") +
        ggtitle(get_state_name(local)) + 
        theme(axis.text.x=element_text(hjust=1,angle=45,size=6),
              panel.grid.minor=element_blank(),
              plot.title=element_text(size=10),
              axis.text.x=element_text(size=6),
              axis.text.y=element_text(size=6),
              axis.title.x=element_text(size=8),
              axis.title.y=element_text(size=8),
              plot.margin=unit(c(0.1,0.8,0.1,0.1),"cm")
              )

     incPlot <- ggplot() +
         geom_line(dat=incBounds,aes(y=inc,x=time,group=quantile),linetype=2,col="red",size=0.5,alpha=0.5) +
         geom_line(dat=bestInc,aes(y=I_H,x=times),col="red",lwd=0.5)+
         geom_polygon(data=polygonI,aes(x=x,y=y),alpha=0.2,fill="red")+
         scale_y_continuous(limits=c(0,0.2))+
         ylab("\nPer capita incidence (red)")+
         xlab("")+
         theme(
             panel.background=element_rect(fill=NA),
             panel.grid=element_blank(),
             axis.line.y = element_line(colour="black"),
             axis.text.y=element_text(size=6,colour="black"),
             axis.title.y=element_text(size=8,angle=-90),
             axis.text.x=element_text(size=6),
             axis.title.x=element_text(size=8)
             )
    if(!ylabel){
        myPlot <- myPlot + ylab("") + theme(plot.margin=unit(c(0,0.4,-0.2,-0.3),"cm"))
        incPlot <- incPlot + ylab("") + theme(plot.margin=unit(c(0,0.4,-0.2,-0.3),"cm"))
    }
    if(!xlabel) myPlot <- myPlot + xlab("")
    
    combined <- overlapPlots(myPlot,incPlot,FALSE)
    
    return(combined)
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
  #ax2 <- gb$children[[2]]
  #ax2$widths <- rev(a2x$widths)
 #ax2$grobs <- rev(ax2$grobs)
  #ax2$grobs[[1]]$x <- ax2$grobs[[1]]$x - unit(1, "npc") + unit(0.1, "cm")
  
  g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
  g <- gtable_add_grob(g, ax, pp$t, length(g$widths)-1, pp$b)
  g <- gtable_add_grob(g, g2$grobs[[7]], pp$t, length(g$widths), pp$b)
  
  # draw it
  return(g)
}

#' Best pars
#'
#' Given an MCMC chain, returns the set of best fitting parameters
#' @param chain the MCMC chain
#' @return a name vector of the best parameters
#' @export
#' @useDynLib zikaProj
get_best_pars <- function(chain){
    tmpNames <- colnames(chain)[2:(ncol(chain)-1)]
    bestPars <- as.numeric(chain[which.max(chain[,"lnlike"]),2:(ncol(chain)-1)])
    names(bestPars) <- tmpNames
    return(bestPars)
}

#' Index pars
#'
#' Given an MCMC chain, returns the parameters at the specified index
#' @param chain the MCMC chain
#' @param index the index
#' @return a named vector of the best parameters
#' @export
#' @useDynLib zikaProj
get_index_pars <- function(chain, index){
    tmpNames <- colnames(chain)[2:(ncol(chain)-1)]
    pars <- as.numeric(chain[index,2:(ncol(chain)-1)])
    names(pars) <- tmpNames
    return(pars)
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
    unfixed_pars <- parTab[parTab$local==local,"names"]
    if(number >= 1) unfixed_pars <- paste(unfixed_pars,".",number,sep="")
    unfixed_pars <- pars[unfixed_pars]
    names(unfixed_pars) <- parTab[parTab$local==local,"names"]
    fixed_pars <- pars[parTab[parTab$local=="all","names"]]
    pars <- c(unfixed_pars, fixed_pars)

    y0s <- generate_y0s(pars["N_H"],pars["density"])

    y <- solveModelSimple(t_pars, y0s, pars)
    y <- y[y[,"times"] >= min(dat[,"startDay"]) & y[,"times"] <= max(dat[,"endDay"]),]

    probs <- generate_micro_curve(pars)

    probM <- generate_probM(y[,"I_M"],probs, pars["N_H"], pars["b"], pars["p_MH"], pars["baselineProb"], 1)
    probM <- average_buckets(probM, dat[,"buckets"])

    predicted <- probM*dat[,"births"]*pars["propn"]
    predicted <- data.frame(number=predicted, day=dat[,"meanDay"])
    
    y[,"I_H"] <- y[,"I_H"]/rowSums(y[,c("I_H","E_H","S_H","R_H")])
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
    numbers <- tmp$ix - 2
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

r0.vector <- function(params){
    NH <- params[,"N_H"]
    NM <- params[,"N_H"]*params[,"density"]
    muM <- 1/params[,"L_M"]
    sigmaM <- 1/params[,"D_EM"]

    muH <- 1/params[,"L_H"]
    gammaH <- 1/params[,"D_IH"]

    b <- params[,"b"]
    pHM <- params[,"p_HM"]
    pMH <- params[,"p_MH"]

    R0 <- (b^2*pHM*pMH*NM*sigmaM)/((sigmaM+muM)*muM*(gammaH+muH)*NH)
    return(R0)
}








#' SEIR dynamics plot
#'
#' Plots SEIR dynamics given a data frame of solved ODEs
#' @param y The data frame or matrix of times and population sizes
#' @param N_H human population size
#' @param N_M mosquito population size
#' @param file.name optional filename at which to save the plot. Must be a PNG. If NULL, does not open the png device.
#' @export
#' @useDynLib zikaProj
plot_dynamics <- function(y, N_H, N_M, file.name = NULL){
    y <- as.data.frame(y)
    n <- ncol(y)
    cols <- c("times","S_M","E_M","I_M","S_C","S_A","S_F","E_C","E_A","E_F","I_C","I_A","I_F","R_C","R_A","R_F")
    y <- y[,1:length(cols)]
    colnames(y) <- cols

    if(!is.null(file.name)){
        png(file.name)
    }
    par(mfrow=c(2,2))
    plot(y$S_M~y$times,col="green",ylim=c(0,N_M),type='l',main="Mosquito Dynamics",xlab="Time",ylab="Incidence")
    lines(y$E_M~y$times,col="red")
    lines(y$I_M~y$times,col="blue")

    plot(y$S_C~y$times,col="green",ylim=c(0,0.3*N_H),type='l',main="Children Dynamics",xlab="Time",ylab="Incidence")
    lines(y$E_C~y$times,col="red")
    lines(y$I_C~y$times,col="blue")
    lines(y$R_C~y$times,col="purple")

    plot(y$S_A~y$times,col="green",ylim=c(0,0.8*N_H),type='l',main="Adult Dynamics",xlab="Time",ylab="Incidence")
    lines(y$E_A~y$times,col="red")
    lines(y$I_A~y$times,col="blue")
    lines(y$R_A~y$times,col="purple")

    plot(y$S_F~y$times,col="green",ylim=c(0,0.004*N_H),type='l',main="First Trimester Dynamics",xlab="Time",ylab="Incidence")
    lines(y$E_F~y$times,col="red")
    lines(y$I_F~y$times,col="blue")
    lines(y$R_F~y$times,col="purple")
    if(!is.null(file.name)){
        dev.off()
    }
    
}

#' Head circumference heatmap
#'
#' Given a matrix or data frame of head sizes over time (rows represent sampling times), plots a heatmap showing distribution and mean head sizes over time.
#' @param dat matrix of head count data. Rows represent sampling times and columns represent individual measurements.
#' @return A ggplot object with the heatmap of head sizes over time. White line shows mean head size.
#' @export
plotDataHeatMap <- function(dat){
tmp <- createCounts(dat)
meanDat <- tmp[[2]]
tmp <- tmp[[1]]

plot <- ggplot(tmp) + geom_raster(aes(x=Day,y=Size,fill=Proportion),interpolate=FALSE) +
    geom_line(data=meanDat,aes(y=y,x=x),linetype=2,colour="white",size=1)+
    scale_fill_gradientn(colours=c("darkblue","red")) +
    scale_y_continuous(expand=c(0,0),breaks=seq(0,max(tmp$Size),by=1),limits=c(19,50),labels=seq(0,max(tmp$Size),by=1))+
    scale_x_continuous(expand=c(0,0),breaks=seq(0,max(tmp$Day),by=max(tmp$Day)/15),labels=round(seq(0,max(tmp$Day),by=max(tmp$Day)/15),digits=0))+
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        text=element_text(size=16,colour="gray20"),
        axis.line=element_line(colour="gray20"),
        axis.line.x = element_line(colour = "gray20"),
        axis.line.y=element_line(colour="gray20")
    )
return(plot)
}
