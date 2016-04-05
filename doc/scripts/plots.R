plot_all <- function(){
  allDat <- read.csv("~/Dropbox/Zika/Data/allDat.csv")
  
  xlabels <- NULL
  xlab <- c("01","02","03","04","05","06","07","08","09","10","11","12")
  for(i in xlab){
    xlabels <- c(xlabels,paste(i,"/2014",sep=""))
  }
  for(i in xlab){
    xlabels <- c(xlabels,paste(i,"/2015",sep=""))
  }
  xlabels <- c(xlabels,"01/2016","02/2016")
  
  correctOrder <- sort(by(allDat[,"microCeph"],allDat[,"local"],sum,simplify=TRUE),decreasing=TRUE)
  correctOrder <- as.factor(names(correctOrder))
  places <- unique(allDat$local)
  places <- correctOrder
  adaptive <- 30000
  p <- NULL
  setwd("~/tmpZika5")
  for(i in 1:length(places)){
    print(places[i])
    tmpDat <- allDat[allDat$local == places[i],]
    filename <- as.character(places[i])
    chains <- NULL
    for(j in 1:2){
      chains[[j]] <- paste(filename,j,"_chain.csv",sep="")
    }
    start <- seq(tmpDat[1,"start"],(nrow(tmpDat)/12) - 1/12,by=1/12)
    end <- seq(1/12,nrow(tmpDat)/12,by=1/12)
    buckets <- cbind(start,end)
    tmpDat <- tmpDat[,c("microCeph","births")]
    tmpDat[,2] <- tmpDat[,2] - tmpDat[,1]
    p[[i]] <- plot_model(chains,tmpDat,buckets,pars,adaptive,100,0.95,0,c(0,500),filename,xlabels)
  }
  ps <- do.call("grid.arrange",c(p,ncol=2))
  return(ps)
}
generate_prediction_intervals <- function(chains, data, times, pars, burnin=NULL, runs=100, level=0.95, smoothing=0){

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
    tmp <- solveModel(pars)

    temp_matrix <- matrix(nrow=nrow(buckets),ncol=runs)

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
    

    colnames(bestMicro) <- colnames(upper_smooth) <- colnames(lower_smooth) <- c("variable","value")
    
    return(list(best=as.data.frame(bestMicro),upper=as.data.frame(upper_smooth), lower=as.data.frame(lower_smooth)))
}

create_polygons <- function(lower,upper){
    bounds <- NULL
    bounds <- rbind(lower[rev(rownames(lower)),],upper)
    colnames(bounds) <- c("x","y")
    return(bounds)    
}

plot_model <- function(chains, data, times, pars, burnin=NULL, runs=100, level=0.95,smoothing=0,bounds=c(0,1000),plot_title=NULL,xlabels=NULL){

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

    title <- ""
    if(!is.null(plot_title)){
        title <- plot_title
    }

    data <- as.data.frame(cbind(times[,2],data[,1]))
    colnames(data) <- c("variable","value")
    data$variable <- data$variable - 1/24
    predictionIntervals$best$variable <- predictionIntervals$best$variable - 1/24

       plot <- ggplot() +
        geom_point(data=data,aes(x=variable,y=value),size=4) +
        geom_line(data=predictionIntervals$best,aes(x=variable,y=value),colour="blue",size=0.8)+
        xlab("Date") +
        ylab("Reported Cases of Microcephaly") +
        scale_y_continuous(limits=c(0,200),breaks=seq(0,200,by=25))+
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
                         axis.text.y=element_text(colour="gray20",size=14))+
        geom_polygon(data=polygon,aes(x=x,y=y),alpha=0.2,fill="blue")
       if(!is.null(xlabels)) plot <- plot + scale_x_continuous(labels=xlabels,breaks=seq(0,25/12,by=1/12))
    return(plot)
}
