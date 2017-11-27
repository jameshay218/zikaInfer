generate_peak_time_table <- function(dat, incDat){
    all_states <- get_epidemic_locations(dat)$include
    peakTimes <- rep(858,length(all_states))
    
    peakTimes[which(all_states == "pernambuco")] <- 804
    peakTimes[which(all_states == "bahia")] <- 855
    peakTimes[which(all_states == "riograndedonorte")] <- 862
    
    peakWidths <- rep(120,length(all_states))
    peakWidths[which(all_states == "pernambuco")] <- 60
    peakWidths[which(all_states == "bahia")] <- 60
    peakWidths[which(all_states == "riograndedonorte")] <- 60

    ## Actual peak times
    data(locationInfo)
    actualPeaks <- data.frame(local=locationInfo[locationInfo$rawName %in%
                                                 c("bahia","pernambuco","riograndedonorte"),"fullName"],
                              peakTime=c(855,804,862),
                              stringsAsFactors = FALSE)
    
    peakTimes <- data.frame(local=locationInfo[locationInfo$rawName %in%
                                               c("bahia","pernambuco","riograndedonorte"),"fullName"],
                            start=peakTimes-peakWidths/2,
                            end=peakTimes+peakWidths/2,
                            stringsAsFactors = FALSE)
    
    peakTimes <- merge(actualPeaks,peakTimes,by="local",all=TRUE)
    
    if(!is.null(incDat)){
        fariaPeaks <- ddply(incDat, c("local"), .fun=function(x) x[which.max(x$inc),"meanDay"])
        colnames(fariaPeaks) = c("local","peakTimeFaria")
        fariaPeaks <- data.frame(fariaPeaks,startFaria=fariaPeaks[,2]-30,endFaria=fariaPeaks[,2]+30,stringsAsFactors=FALSE)
        peakTimes <- merge(peakTimes, fariaPeaks, by="local",all=TRUE)
    }
    return(peakTimes)
}

#' @export
main_model_fits <- function(chainWD = "~/Documents/Zika/28.02.2017_chains/multi_all1/model_1",
                            datFile = "~/Documents/Zika/Data/allDat28.02.17.csv",
                            incDatFile = "~/Documents/Zika/Data/inc_data_120317.csv",
                            runs=200,
                            datCutOff=1214,
                            incScale=2000
                            ){
    setwd(chainWD)
    
    chain <- zikaProj::load_mcmc_chains(
        location = chainWD,
        asList = FALSE,
        convertMCMC = FALSE,
        unfixed = FALSE,
        thin = 10,
        burnin = 750000)
    
    parTab <- read_inipars()
    ts <- seq(0,3003,by=1)
    
    dat <- read.csv(datFile,stringsAsFactors=FALSE)
    incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)
    tmp <- NULL
    tmpDat <- dat[dat$local %in% unique(parTab$local),]
    states <- unique(tmpDat$local)
    for(state in states){
        tmp[[state]] <- plot_setup_data(
            chain, tmpDat, incDat,
            parTab, ts, state,
            runs=runs,startDay=365,
            noWeeks=150,
            perCap=TRUE)
    }
    
    incDat$meanDay <- (incDat$startDay + incDat$endDay)/2
    incDat$local <- convert_name_to_state_factor(incDat$local)
    
    peakTimes <- generate_peak_time_table(dat,incDat)
    
    labels <- rep(getDaysPerMonth(3),4)
    labels <- c(0,cumsum(labels))
    labels_names <- as.yearmon(as.Date(labels,origin="2013-01-01"))
    
    microBounds <- NULL
    for(i in 1:length(tmp)) microBounds <- rbind(microBounds, tmp[[i]][["microBounds"]])
    incBounds <- NULL
    for(i in 1:length(tmp)) incBounds <- rbind(incBounds, tmp[[i]][["incBounds"]])
   
    dat$meanDay <- rowMeans(dat[,c("startDay","endDay")])
    dat <- dat[dat$startDay < datCutOff,]
    dat$local <- convert_name_to_state_factor(dat$local)
    microBounds$local <- convert_name_to_state_factor(microBounds$local)
    incBounds$local <- convert_name_to_state_factor(incBounds$local)

    inc_scales <- ddply(microBounds, "local", function(x) max(x$upper))

    ## Need to scale all of the incidence data for plotting purposes
    for(state in unique(incBounds$local)){
        scale <- as.numeric(inc_scales[inc_scales$local == state, 2])
        tmp_bounds <- incBounds[incBounds$local == state,]
        tmp_bounds[,c("lower","upper","best")] <- tmp_bounds[,c("lower","upper","best")] * incScale
        incBounds[incBounds$local == state,] <- tmp_bounds        
    }
    
    states <- unique(dat$local)
    #peakTimes <- peakTimes[peakTimes$local %in% states,]
    p <- ggplot() +
        geom_rect(data=peakTimes,
                  aes(xmin=start,xmax=end,ymin=0,ymax=Inf,group=local),
                  alpha=0.5,fill="red") +
        geom_rect(data=peakTimes,
                  aes(xmin=startFaria,xmax=endFaria,ymin=0,ymax=Inf,group=local),
                  alpha=0.5,fill="orange") +
        geom_vline(data=peakTimes,
                   aes(xintercept=peakTime,group=local),col="black",lty="dashed") +
        geom_ribbon(data=microBounds,aes(ymin=lower,ymax=upper,x=time),fill="blue",alpha=0.5) +
        geom_ribbon(data=incBounds,aes(ymin=lower,ymax=upper,x=time),fill="green",alpha=0.5) +
        geom_line(data=microBounds,aes(x=time,y=best),colour="blue") +
        geom_line(data=incBounds,aes(x=time,y=best),colour="green") +
        geom_point(data=dat, aes(x=meanDay, y = microCeph/births), size=0.6) +
        facet_wrap(~local,ncol=4) +
        scale_y_continuous(limits=c(0,0.02),breaks=seq(0,0.025,by=0.005),expand=c(0,0),
                           sec.axis=sec_axis(~.*(1/incScale)))+
        theme_bw() + 
        theme(axis.text.x=element_text(angle=90,vjust=0.5,size=8,family="Arial"),
              axis.text.y=element_text(size=8,family="Arial"),
              axis.title=element_text(size=10,family="Arial"),
              strip.text=element_text(size=8,family="Arial"),
              panel.grid.minor = element_blank()) +
        ylab("Per birth microcephaly incidence (blue)") +
        xlab("") +
        scale_x_continuous(breaks=labels,labels=labels_names)

     return(p)
}

#' @export
indiv_model_fit <- function(datFile = "~/Documents/Zika/Data/northeast_microceph.csv",
                            incFile = "~/Documents/Zika/Data/northeast_zikv.csv",
                            local = "bahia",
                            localName = "Northeast Brazil NEJM",
                            incScale=50,
                            runs=1000,
                            ylim=0.1,
                            xlim=NULL,
                            bot=FALSE,
                            standalone=FALSE,
                            parTab=NULL,
                            chain=NULL,
                            forecast=FALSE,
                            forecastPostChange=FALSE,
                            weeks=FALSE){
    ts <- seq(0,3003,by=1)

    microDat <- read.csv(datFile,stringsAsFactors = FALSE)
    microDat <- microDat[microDat$local == local,]

    incDat <- NULL
    if(!is.null(incFile)){
        incDat <- read.csv(incFile,stringsAsFactors=FALSE)    
        incDat <- incDat[incDat$local == local,]
        incDat$meanDay <- rowMeans(incDat[,c("startDay","endDay")])
        peakTime <- incDat[which.max(incDat$inc),"meanDay"]
    }
    
    if(forecast){
        f <- create_forecast_function(parTab, microDat, incDat, ts, FALSE)
    } else {
        f <- create_forecast_function(parTab,microDat,incDat,ts,TRUE)
    }
    
    samples <-  sample(nrow(chain), runs)
    microCurves <- NULL
    
    f1 <- create_forecast_function(parTab, microDat, incDat,ts,TRUE)
    for(i in 1:length(samples)){
        pars <- get_index_pars(chain,samples[i])
        ## If we want to assume that no behaviour changed then we need to manually set these parameters
        ## for forecasting method
        if(forecastPostChange){
          pars["incPropn2"] <- pars["incPropn"]
          pars["abortion_rate"] <- pars["avoided_births"] <- pars["birth_reduction"] <- 0
          pars["propn2"] <- pars["propn"]
        }
        ## If forecasting version, we had baseline prob on exponential scale
        if(forecast) pars["baselineProb"] <- exp(pars["baselineProb"])
        microCurves <- rbind(microCurves, f1(pars,TRUE)$microCeph)
    }
    predict_bounds <- as.data.frame(t(sapply(unique(microCurves$time),function(x) quantile(microCurves[microCurves$time==x,"microCeph"],c(0.025,0.5,0.975)))[c(1,3),]))
    bestPars <- get_best_pars(chain)
    ## If forecasting version, we had baseline prob on exponential scale
    if(forecast) bestPars["baselineProb"] <- exp(bestPars["baselineProb"])
    ## If we want to assume that no behaviour changed then we need to manually set these parameters
    ## for forecasting method
    if(forecastPostChange){
      bestPars["incPropn2"] <- bestPars["incPropn"]
      bestPars["abortion_rate"] <- bestPars["avoided_births"] <- bestPars["birth_reduction"] <- 0
      bestPars["propn2"] <- bestPars["propn"]
    }
    best_predict <- f1(bestPars,TRUE)$microCeph
    predict_bounds <- cbind(predict_bounds, best_predict[,2])
    colnames(predict_bounds) <- c("lower","upper","best")
    predict_bounds$time <- best_predict[,1]
    
    if(!(local %in% parTab$local)){
        microDat$local <- "bahia"
        if(!is.null(incDat)) incDat$local <- "bahia"
        local <- "bahia"
    }
    microDat <- microDat[,c("startDay","endDay","microCeph","buckets","births","local")]
    if(!is.null(incDat)) incDat <- incDat[,c("startDay","endDay","buckets","inc","N_H","local")]

    labels <- rep(getDaysPerMonth(3),4)
    labels <- c(0,cumsum(labels))
    labels_names <- as.yearmon(as.Date(labels,origin="2013-01-01"))
    tmp <- plot_setup_data(chain, microDat, incDat,parTab, ts, local,200,365, noMonths=36,noWeeks=150,perCap=TRUE, forecast=forecast, forecastPostChange=forecastPostChange)
    microBounds <- tmp[["microBounds"]]
    incBounds <- tmp[["incBounds"]]
    microDat$meanDay <- rowMeans(microDat[,c("startDay","endDay")])
    if(!is.null(incDat)){
        incDat$meanDay <- rowMeans(incDat[,c("startDay","endDay")])
        incDat$local <- localName
    }
    
    microBounds$local <- localName
    microDat$local <- localName    
    incBounds[,c("lower","upper","best")] <- incBounds[,c("lower","upper","best")]
    incBounds$local <- localName
    data(locationInfo)
    peakTimes <- locationInfo[locationInfo$rawName == local,c("peakTime","peakTimeRange")]
    peakTimes$local <- localName
    peakTimes$peakUpper <- peakTimes$peakTime + peakTimes$peakTimeRange/2
    peakTimes$peakLower <- peakTimes$peakTime - peakTimes$peakTimeRange/2
    
    p <- ggplot() +
        geom_ribbon(data=microBounds,aes(ymin=lower,ymax=upper,x=x,fill="Model microcephaly (two waves)"),alpha=0.5) +
        geom_ribbon(data=incBounds,aes(ymin=lower/incScale,ymax=upper/incScale,x=x,fill="Model ZIKV"),alpha=0.5) +
        geom_line(data=microBounds,aes(x=x,y=best,colour="Model microcephaly (two waves)")) +
        geom_line(data=incBounds,aes(x=x,y=best/incScale,colour="Model ZIKV"))
    if(!is.null(incFile)){
        p <- p + geom_line(data=incDat,aes(x=meanDay,y=inc/N_H/incScale,col="Reported ZIKV"),linetype="longdash")
        if(forecast){
            p <- p +
                geom_ribbon(data=predict_bounds,aes(ymin=lower,ymax=upper,x=time,fill="Model microcephaly (single wave)"),alpha=0.5) +
                geom_line(data=predict_bounds,aes(x=time,y=best,colour="Model microcephaly (single wave)"))
        }
    } else {
        p <- p +
            geom_rect(data=peakTimes,
                      aes(xmin=peakUpper,xmax=peakLower,ymin=0,ymax=Inf,group=local),
                      alpha=0.5,fill="red") +
            geom_vline(data=peakTimes,
                       aes(xintercept=peakTime,group=local),col="black",lty="dashed")
    }
    p <- p + geom_point(data=microDat, aes(x=meanDay, y = microCeph/births,col="Reported microcephaly"), size=0.6) +
        facet_wrap(~local,scales="free_y",ncol=1) +
        scale_y_continuous(limits=c(0,ylim),breaks=seq(0,ylim,by=ylim/5),expand=c(0,0),sec.axis=sec_axis(~.*(incScale),name="Reported per capita\nZIKV infection incidence (red)"))+
        theme_classic()
    if(standalone){
        p <- p + theme(axis.text.y=element_text(size=8,family="Arial"),
                       axis.title=element_text(size=8,family="Arial"),
                       strip.text=element_blank(),
                       axis.text.x=element_text(size=8, family="Arial",hjust=1,angle=45),
                       panel.grid.minor = element_blank())
        if(bot){
            p <- p + theme(
                         axis.text.x=element_blank(),
                         axis.title.x=element_blank(),
                         axis.ticks.x = element_blank())
            }


    } else {
        if(bot){
            p <- p +
                theme(axis.text.y=element_text(size=8,family="Arial"),
                      axis.title.x=element_blank(),
                      axis.title=element_text(size=10,family="Arial"),
                      strip.text=element_text(size=8,family="Arial"),
                      axis.title.y=element_blank(),
                      axis.text.x=element_text(angle=90,hjust=0.5,size=8,family="Arial"), 
                      axis.ticks.x = element_blank(),
                      panel.grid.minor = element_blank())
        } else {
            p <- p +        
                theme(axis.text.y=element_text(size=8,family="Arial"),
                      axis.title.x=element_blank(),
                      axis.title=element_text(size=10,family="Arial"),
                      strip.text=element_text(size=8,family="Arial"),
                      axis.title.y=element_blank(),
                      axis.text.x=element_blank(), 
                      axis.ticks.x = element_blank(),
                      panel.grid.minor = element_blank())
        }
    }
    x_lower <- 365
    if(!is.null(xlim)) x_lower <- xlim
    p <- p +
        scale_x_continuous(limits=c(x_lower,max(labels)),breaks=labels,labels=labels_names)+#,expand=c(0,0))+
        ylab("Reported per birth\nmicrocephaly incidence (black)") +
        xlab("") +
        geom_blank(data=data.frame(x=rep(0,5),
                                   group=c("Model microcephaly (two waves)",
                                           "Model ZIKV",
                                           "Reported ZIKV",
                                           "Model microcephaly (single wave)",
                                           "Reported microcephaly")),
                   aes(y=x,color=group,fill=group))+
        scale_colour_manual(name="",values=c("Model microcephaly (two waves)"="purple",
                                             "Model ZIKV"="green",
                                             "Reported ZIKV"="red",
                                             "Model microcephaly (single wave)"="blue",
                                             "Reported microcephaly"="black"))+
        scale_fill_manual(name="",values=c("Model microcephaly (two waves)"="purple",
                                           "Model ZIKV"="green",
                                           "Reported ZIKV"=NA,
                                           "Model microcephaly (single wave)"="blue",
                                           "Reported microcephaly"=NA)) +
        theme(legend.position=c(1,1),
              legend.justification = c(1,1),
              legend.text=element_text(size=8,family="Arial"),
              legend.background = element_blank())
    return(p)
}



