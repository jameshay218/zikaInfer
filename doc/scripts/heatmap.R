library(ggplot2)
library(reshape2)
library(reldist)
#pars[[1]] <- seq(0,1,by=1/365)
bs <- c(50,100,200)
allPlots <- NULL


for(iii in 1:length(bs)){
  print(iii)
  pars[[3]][16] <- bs[iii]
  y <- solveModel(pars)
  y <- y[y$times < 6,]
  mu_I <- pars[[3]][3] # Define distributions for head sizes
  mu_N <- pars[[3]][5]
  sd_I <- pars[[3]][4]
  sd_N <- pars[[3]][6]
  
  probMicro <- 1 # Probability that infection leads to microcephaly
  sampFreq <- 7
  
  # For each time point, generate the mixture model with weighting depending on the probability of generating 
  # a microcephaly infant vs. normal infant
  alphas <- calculate_alphas(as.matrix(unname(y[, c("If", "Sf","Ef", "Rf")])), probMicro, sampFreq)
  distributions <- matrix(nrow=length(alphas),ncol=501)
  for(i in 1:length(alphas)){
    distributions[i,] <- alphas[i]*dnorm(seq(0,50,by=0.1),mu_I,sd_I) + (1-alphas[i])*dnorm(seq(0,50,by=0.1),mu_N,sd_N)
  }
  
  # In the inference framework, we would compare this distribution over time to known head circumference distributions.
  
  
  # Calculate mean head sizes on each day
  means <- NULL
  for(i in 1:nrow(distributions)){
    x1 <- seq(0,ncol(distributions)-1,by=1)
    weights <- distributions[i,]
    #quantiles <- wtd.quantile(x1,probs=c(0.025,0.975),weight=weights)
    tmp <- NULL
    for(j in 1:ncol(distributions)){
      tmp[j] <- distributions[i,j] * j
    }
    means[i] <- sum(tmp)/sum(distributions[i,1:ncol(distributions)])/10
  }
  
  
  tmp <- melt(distributions)
  colnames(tmp) <- c("Week","Size","Proportion")
  meanDat <- data.frame(y=means,x=seq(1,max(tmp$Week),by=1))
  tmp$Size <- tmp$Size/10
                                          #labels=seq(-1,max(tmp$Size)-1,by=1),
  
  plot <- ggplot(tmp) + geom_raster(aes(x=Week,y=Size,fill=Proportion),interpolate=TRUE) +
      geom_line(data=meanDat,aes(y=y,x=x),linetype=2,colour="white",size=1)+
      scale_fill_gradientn(colours=c("darkblue","red")) +
      scale_y_continuous(expand=c(0,0),breaks=seq(0,max(tmp$Size),by=5),limits=c(25,40))+
      scale_x_continuous(expand=c(0,0))+
    ylab("")+
    xlab("")+
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
         text=element_text(size=16,colour="gray20"),
          axis.line=element_line(colour="gray20"),
          legend.direction="horizontal",
         legend.title=element_blank(),
          axis.line.x = element_line(colour = "gray20"),
          axis.line.y=element_line(colour="gray20")
      )
  if(iii == 1){
    plot <- plot + theme(legend.position="none")
    }
  if(iii == 2){
    g_legend<-function(a.gplot){
     tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)}
    
    mylegend<-g_legend(plot)
    plot <- plot + ylab("Head Circumference (cm)") + theme(legend.position="none")
  }
  if(iii == 3){
    plot <- plot + xlab("Week") + theme(legend.position="none")
  }
  allPlots[[iii]] <- plot
}

final <- grid.arrange(arrangeGrob(allPlots[[1]],allPlots[[2]],allPlots[[3]],nrow=3),mylegend,nrow=2,heights=c(10,1))
ggsave("circumferences.png",final,height=8,width=8)
