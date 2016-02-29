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
    cols <- c("times","Sm","Em","Im","Sc","Sa","Sf","Ec","Ea","Ef","Ic","Ia","If","Rc","Ra","Rf","IfA","fB")
    y <- y[,1:length(cols)]
    colnames(y) <- cols

    if(!is.null(file.name)){
        png(file.name)
    }
    par(mfrow=c(2,2))
    plot(y$Sm~y$times,col="green",ylim=c(0,N_M),type='l',main="Mosquito Dynamics",xlab="Time",ylab="Incidence")
    lines(y$Em~y$times,col="red")
    lines(y$Im~y$times,col="blue")

    plot(y$Sc~y$times,col="green",ylim=c(0,0.3*N_H),type='l',main="Children Dynamics",xlab="Time",ylab="Incidence")
    lines(y$Ec~y$times,col="red")
    lines(y$Ic~y$times,col="blue")
    lines(y$Rc~y$times,col="purple")

    plot(y$Sa~y$times,col="green",ylim=c(0,0.8*N_H),type='l',main="Adult Dynamics",xlab="Time",ylab="Incidence")
    lines(y$Ea~y$times,col="red")
    lines(y$Ia~y$times,col="blue")
    lines(y$Ra~y$times,col="purple")

    plot(y$Sf~y$times,col="green",ylim=c(0,0.004*N_H),type='l',main="First Trimester Dynamics",xlab="Time",ylab="Incidence")
    lines(y$Ef~y$times,col="red")
    lines(y$If~y$times,col="blue")
    lines(y$Rf~y$times,col="purple")
    if(!is.null(file.name)){
        dev.off()
    }
    
}

#' Head circumference heatmap
#'
#' Given a matrix or data frame of head sizes over time (rows represent sampling times), plots a heatmap showing distribution and mean head sizes over time.
#' @param
#'
