plot_all <- function(){
 # skip <- c("bahia_2","espiritosanto_2","goias_1","para_1","riodejaneiro_3","saopaulo_2")
  countries <- c("pernambuco",
                 "bahia",
                 "saopaulo",
                 "paraiba",
                 "maranhao",
                 "ceara",
                 "sergipe",
                 "riodejaneiro",
                 "piaui",
                 "riograndedonorte",
                 "minasgerais",
                 "matogrosso",
                 "alagoas",
                 "para",
                 "acre",
                 "goias",
                 "espiritosanto",
                 "tocantins")
  
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
  
  densityM <- 3
  allDat <- read.csv("~/Dropbox/Zika/Data/allDat07.04.16.csv")
  
  xlabels <- NULL
  xlab <- c("01","02","03","04","05","06","07","08","09","10","11","12")
  for(i in xlab){
    xlabels <- c(xlabels,paste(i,"/2014",sep=""))
  }
  for(i in xlab){
    xlabels <- c(xlabels,paste(i,"/2015",sep=""))
  }
  xlabels <- c(xlabels,"01/2016","02/2016")
  
  xlabels <- c("01/2014","04/2014","07/2014","10/2014","01/2015","04/2015","07/2015","10/2015","01/2016")
  
  #correctOrder <- sort(by(allDat[,"microCeph"],allDat[,"local"],sum,simplify=TRUE),decreasing=TRUE)
  correctOrder <- sort(unique(allDat[allDat$local %in% countries,"local"]))
  correctOrder <- as.factor(correctOrder)
  places <- correctOrder
  adaptive <- 30000
  p <- NULL
  setwd("~/tmpZika6")
  for(i in 1:17){
    print(places[i])
    tmpDat <- allDat[allDat$local == places[i],]
    start <- seq(0,(nrow(tmpDat)/12) - 1/12,by=1/12)
    end <- seq(1/12,nrow(tmpDat)/12,by=1/12)
    N_H <- tmpDat[1,"NH"]
    L_H <- tmpDat[1,"LH"]
 
    pars <- setupListPars(duration=end[length(end)], N_H = N_H, N_M = N_H*densityM)  
    pars[[1]][1] <- 4
    pars[[3]]["sampFreq"] <- 30
    pars[[3]]["probMicro"] <- 0.2
    pars[[3]]["baselineProb"] <- 0.01
    pars[[3]]["epiStart"] <- 1.2
    pars[[3]]["L_H"] <- L_H
    pars[[3]]["constSeed"] <- 0
    
    filename <- as.character(places[i])
    chains <- NULL
    index <- 1
    print("lol")
    for(j in 1:3){
      tmp <- paste(filename,"_",j,sep="")
      #if(tmp %in% skip) next
      chains[[index]] <- paste(filename,"_",j,"_chain.csv",sep="")
      index <- index + 1
    }
    start <- seq(0,(nrow(tmpDat)/12) - 1/12,by=1/12)
    end <- seq(1/12,nrow(tmpDat)/12,by=1/12)
    buckets <- cbind(start,end)
    tmpDat <- tmpDat[,c("microCeph","births")]
    print(places[i])
    print(country_names[filename])
    tmpPlot <- plot_model(chains,tmpDat,buckets,pars,adaptive,100,0.95,0,c(0,500),country_names[filename],xlabels,200,TRUE,TRUE,TRUE)
    p[[i]] <- overlapPlots(tmpPlot[[1]],tmpPlot[[2]],FALSE)
  }
  ps <- do.call("grid.arrange",c(p,ncol=5))
  return(ps)
}

generate_prediction_intervals <- function(chains, data, times, pars, burnin=NULL, runs=100, level=0.95, smoothing=0){
  pars[[1]]["step"] <- 0.001
    ## Get quantiles for prediction
    upper_level <- level+((1-level)/2)
    lower_level <- (1-level)/2

    
    if(!is.null(burnin)) burnin_start<- burnin
    else burnin_start <- 1

    ## Combine chains so can sample from all
    mcmc <- NULL
    for(i in 1:length(chains)){
        tmpMCMC <- fread(chains[[i]],data.table=FALSE)
        mcmc <- rbind(mcmc, tmpMCMC[burnin_start:nrow(tmpMCMC),])
    }
    colnames(mcmc) <- c("sampno",names(pars[[3]]),"r0","lnlike")
    ################################################################
    ## GET BEST FIT TRAJECTORY
    ################################################################
    bestPars <- as.numeric(mcmc[which.max(mcmc[,"lnlike"]),c(names(pars[[3]]),"lnlike")])
    names(bestPars) <- c(names(pars[[3]]),"lnlike")
    bestTrajectory <- solveModel(list(pars[[1]],pars[[2]],bestPars))
    
    N_H <- sum(bestTrajectory[1,4:ncol(bestTrajectory)])
    
    bestY <- cbind(bestTrajectory[,"times"],(bestTrajectory[,"I_F"]+bestTrajectory[,"I_C"]+bestTrajectory[,"I_A"]))
    bestY[,1] <- bestY[,1] - 0.5
    bestY <- bestY[bestY[,1] < max(buckets) & bestY[,1] >= min(buckets),]
    bestY[,2] <- bestY[,2]/N_H
    
    probMicro <- bestPars["probMicro"]
    baselineProb <- bestPars["baselineProb"]
    
    ## Calculate probability of microcephaly for these data
    bestAlpha <- calculate_alphas_prob_buckets(
        as.matrix(unname(bestTrajectory[,c("times","I_F","S_F","E_F","R_F")])),
        probMicro,
        baselineProb,
        times
    )

    ## Using probability of microcephaly and total number of births, get predicted number of microcephaly births
    bestMicro <- data[,2]*bestAlpha
    bestMicro <- cbind(times[,2],bestMicro)
################################################################
    
    ## Get sample indices
    samples <- sample(1:nrow(mcmc),runs,replace=T)

    ## Get test model output for times and structure
    y1 <- tmp <- solveModel(pars)

    temp_matrix <- matrix(nrow=nrow(buckets),ncol=runs)
    temp_matrix_incidence <- matrix(nrow=nrow(bestY),ncol=runs)
    
    ## For each sample
    for(i in 1:runs){
        ## Extract sample parameters from chain
        params <- as.numeric(mcmc[samples[i],names(pars[[3]])])
        names(params) <- names(pars[[3]])
        probMicro <- params["probMicro"]
        baselineProb <- params["baselineProb"]

        ## Solve model with these parameters
        y <- solveModel(list(pars[[1]],pars[[2]],params))
        
        ## Calculate probability of microcephaly for these data
        tmpAlpha <- calculate_alphas_prob_buckets(
            as.matrix(unname(y[,c("times","I_F","S_F","E_F","R_F")])),
            probMicro,
            baselineProb,
            times
        )

        ## Using probability of microcephaly and total number of births, get predicted number of microcephaly births
        numberMicro <- data[,2]*tmpAlpha
       
        temp_matrix[,i] <- numberMicro
        y[,"times"] <- y[,"times"] - 0.5
        
        y <- y[y[,"times"] < max(buckets) & y[,"times"] >= min(buckets),]
        if(nrow(y) > nrow(bestY)) y <- y[1:nrow(bestY),]
        else if(nrow(y) < nrow(bestY)){
                  wow <- nrow(y)
          wow1 <- nrow(bestY)
          while(wow < wow1){
            y <- rbind(y, c(bestY[wow,1],rep(0,ncol(bestY)-1)))
            wow <- wow + 1
          }
        
        }
        
        incidence <- y[,"I_F"] + y[,"I_C"] + y[,"I_A"]
        incidence <- incidence/N_H
        temp_matrix_incidence[,i] <- incidence
    }

    ## Create empty matrices for upper and lower bounds of these trajectories
    upper <- lower <- matrix(nrow=nrow(buckets),ncol=2)
    upper[,1] <- lower[,1] <- times[,2]
    for(i in 1:nrow(upper)){
        upper[i,2] <- quantile(temp_matrix[i,], upper_level)
        lower[i,2] <- quantile(temp_matrix[i,],lower_level)
    }
    upper_smooth <- as.data.frame(predict(smooth.spline(upper,spar=smoothing)))
    lower_smooth <- as.data.frame(predict(smooth.spline(lower,spar=smoothing)))
    
    upperI <- lowerI <- matrix(nrow=nrow(bestY),ncol=2)
    upperI[,1] <- lowerI[,1] <- bestY[,1]
    
    for(i in 1:nrow(upperI)){
      upperI[i,2] <- quantile(temp_matrix_incidence[i,],upper_level)
      lowerI[i,2] <- quantile(temp_matrix_incidence[i,],lower_level)
    }
    upper_smoothI <- as.data.frame(predict(smooth.spline(upperI,spar=smoothing)))
    lower_smoothI <- as.data.frame(predict(smooth.spline(lowerI,spar=smoothing)))
    
    colnames(upper_smoothI) <- colnames(lower_smoothI) <- colnames(bestY) <- colnames(bestMicro) <- colnames(upper_smooth) <- colnames(lower_smooth) <- c("variable","value")
    
    return(list(best=as.data.frame(bestMicro),upper=as.data.frame(upper_smooth), lower=as.data.frame(lower_smooth),bestI=as.data.frame(bestY),upperI=as.data.frame(upper_smoothI),lowerI=lower_smoothI))
}

create_polygons <- function(lower,upper){
    bounds <- NULL
    bounds <- rbind(lower[rev(rownames(lower)),],upper)
    colnames(bounds) <- c("x","y")
    return(bounds)    
}

plot_model <- function(chains, data, times, pars, burnin=NULL, runs=100, level=0.95,smoothing=0,bounds=c(0,1000),plot_title=NULL,xlabels=NULL, upperLim=NULL,
                       yflag1 =FALSE,yflag2=FALSE,xflag=FALSE){
    low <- bounds[1]
    high <- bounds[2]
    
    predictionIntervals <- generate_prediction_intervals(chains, data, times, pars, burnin, runs, level, smoothing)
    predictionIntervals$best[predictionIntervals$best[,2] < low,] <- low
    predictionIntervals$upper[predictionIntervals$upper[,2] < low,] <- low
    predictionIntervals$lower[predictionIntervals$lower[,2] < low,] <- low

    predictionIntervals$best[predictionIntervals$best[,2] > high,] <- high
    predictionIntervals$upper[predictionIntervals$upper[,2] > high,] <- high
    predictionIntervals$lower[predictionIntervals$lower[,2] > high,] <- high

    predictionIntervals$upper[,1] <- predictionIntervals$upper[,1] -1/24
    predictionIntervals$lower[,1] <- predictionIntervals$lower[,1] -1/24
    
    polygon <- create_polygons(predictionIntervals$lower, predictionIntervals$upper)

    polygonI <- create_polygons(predictionIntervals$lowerI, predictionIntervals$upperI)
    
    
    title <- ""
    if(!is.null(plot_title)){
        title <- plot_title
    }

    data <- as.data.frame(cbind(times[,2],data[,1]))
    colnames(data) <- c("variable","value")
    data$variable <- data$variable - 1/24
    predictionIntervals$best$variable <- predictionIntervals$best$variable - 1/24
    upperLimit <- max(max(data$value,predictionIntervals$upper$value)) * 1.2
    if(!is.null(upperLim)) upperLimit <- upperLim
       plot <- ggplot() +
        geom_point(data=data,aes(x=variable,y=value),size=2) +
        geom_line(data=predictionIntervals$best,aes(x=variable,y=value),colour="blue",size=0.8)+
        
        scale_y_continuous(limits=c(0,upperLimit))+
        ggtitle(title) +
              theme_bw()+
        theme(axis.text.x=element_text(hjust=1,angle=45),panel.grid.minor=element_blank(),
              text=element_text(size=14,colour="gray20"),
              plot.title=element_text(size=16),
              legend.text=element_text(size=12,colour="gray20"),
              panel.grid.minor = element_blank(),
                        axis.line=element_line(colour="gray20"),
              
              plot.margin=unit(c(0,0.5,-0.5,-0.5),"cm"),
              panel.border=element_rect(fill=NA,colour="gray20"))+
        geom_polygon(data=polygon,aes(x=x,y=y),alpha=0.2,fill="blue")
       if(xflag){
        if(!is.null(xlabels)) plot <- plot + scale_x_continuous(labels=xlabels,breaks=seq(0,25/12,by=3/12))
        plot <- plot + xlab("") + theme(axis.line.x = element_line(colour = "gray20"),axis.text.x=element_text(colour="gray20",size=14))
       } else {
         plot <- plot + theme(axis.line.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + xlab("")
       }
       if(yflag1) plot <- plot + ylab("Reported Cases of Microcephaly") + ylab("") + theme(axis.line.y=element_line(colour="gray20"),axis.text.y=element_text(colour="gray20",size=14))
        else plot <- plot + theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(), axis.text.y=element_blank()) + ylab("")
       
      
       incPlot <- ggplot() + geom_line(dat=predictionIntervals$bestI,aes(variable,value),colour="red") + xlab("") + ylab("") + 
     theme(
       panel.background=element_rect(fill=NA),
       panel.grid=element_blank(),axis.text.y=element_text(colour="gray20",size=14),
       axis.line.y=element_line(colour="gray20"),
       axis.title.y=element_text(colour="gray20",size=16,angle=-90)) + 
     scale_y_continuous(limits=c(0,0.25),labels=c(" 0.00"," 0.05"," 0.10", " 0.15"," 0.20"," 0.25"))+
     geom_line(dat=predictionIntervals$upperI,aes(variable,value),colour="red",linetype="dashed")+
     geom_line(dat=predictionIntervals$lowerI,aes(variable,value),colour="red",linetype="dashed")
     if(yflag2) incPlot <- incPlot + ylab("\nPer Capita Incidence") + ylab("")
     else incPlot <- incPlot + theme( axis.line.y=element_blank(),
                                      axis.line.x=element_blank(),
                                      axis.ticks.y=element_blank(),
                                      axis.ticks.x=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.text.y=element_blank())
       
    return(list(plot,incPlot))
}


plot_line_dat <- function(pars, tmpDat, title=""){
  xlabels <- NULL
  xlab <- c("01","02","03","04","05","06","07","08","09","10","11","12")
  for(i in xlab){
    xlabels <- c(xlabels,paste(i,"/2014",sep=""))
  }
  for(i in xlab){
    xlabels <- c(xlabels,paste(i,"/2015",sep=""))
  }
  xlabels <- c(xlabels,"01/2016","02/2016")
  start <- seq(0,(nrow(tmpDat)/12) - 1/12,by=1/12)
  end <- seq(1/12,nrow(tmpDat)/12,by=1/12)
  times <- cbind(start,end)
  
  probMicro <- pars[[3]]["probMicro"]
  baselineProb <- pars[[3]]["baselineProb"]
  
  ## Solve model with these parameters
  y <- solveModel(pars)
  ## Calculate probability of microcephaly for these data

  tmpAlpha <- calculate_alphas_prob_buckets(
    as.matrix(unname(y[,c("times","I_F","S_F","E_F","R_F")])),
    probMicro,
    baselineProb,
    times
  )
  
  numberMicro <- tmpDat[,2]*tmpAlpha
  modelDat <- data.frame("variable"=times[,2],"value"=numberMicro)
  modelDat$variable <- modelDat$variable - 1/24
  data <- as.data.frame(cbind(buckets[,2],tmpDat[,1]))
  colnames(data) <- c("variable","value")
  data$variable <- data$variable - 1/24
  upperLimit <- max(max(data$value,modelDat$value)) * 1.2
  plot <- ggplot() +
    geom_point(data=data,aes(x=variable,y=value),size=4) +
    geom_line(data=modelDat,aes(x=variable,y=value),colour="blue",size=0.8)+
    xlab("Date") +
    ylab("Reported Cases of Microcephaly") +
    scale_y_continuous(limits=c(0,upperLimit))+
    ggtitle(title) +
    theme(axis.text.x=element_text(angle=90,hjust=1),panel.grid.minor=element_blank(),
          text=element_text(size=16,colour="gray20"),
          plot.title=element_text(size=28),
          legend.text=element_text(size=14,colour="gray20"),
          panel.grid.minor = element_blank(),
          axis.line=element_line(colour="gray20"),
          axis.line.x = element_line(colour = "gray20"),
          axis.line.y=element_line(colour="gray20"),
          axis.text.x=element_text(colour="gray20"),
          axis.text.y=element_text(colour="gray20",size=14))
  return(plot)
  
}

get_best_pars <- function(chains, burnin){
  mcmc <- NULL
  for(i in 1:length(chains)){
    tmpMCMC <- fread(chains[[i]],data.table=FALSE)
    mcmc <- rbind(mcmc, tmpMCMC[burnin:nrow(tmpMCMC),])
  }
  colnames(mcmc) <- c("sampno",names(pars[[3]]),"r0","lnlike")
  ################################################################
  ## GET BEST FIT TRAJECTORY
  ################################################################
  bestPars <- as.numeric(mcmc[which.max(mcmc[,"lnlike"]),c(names(pars[[3]]),"lnlike")])
  names(bestPars) <- c(names(pars[[3]]),"lnlike")
  return(bestPars)
  
}

get_best_trajectory <- function(chains, burnin){
  bestPars <- get_best_pars(chains, burnin)
  bestTrajectory <- solveModel(list(pars[[1]],pars[[2]],bestPars))
  return(bestTrajectory)
}
get_best_micro <- function(chains, burnin, times, data){
  bestPars <- get_best_pars(chains, burnin)
  bestTrajectory <- solveModel(list(pars[[1]],pars[[2]],bestPars))
  probMicro <- bestPars["probMicro"]
  baselineProb <- bestPars["baselineProb"]
  
  ## Calculate probability of microcephaly for these data
  bestAlpha <- calculate_alphas_prob_buckets(
    as.matrix(unname(bestTrajectory[,c("times","I_F","S_F","E_F","R_F")])),
    probMicro,
    baselineProb,
    times
  )
  
  ## Using probability of microcephaly and total number of births, get predicted number of microcephaly births
  bestMicro <- data[,2]*bestAlpha
  bestMicro <- cbind(times[,2],bestMicro)
  return(bestMicro)
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

