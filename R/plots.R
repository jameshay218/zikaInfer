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
