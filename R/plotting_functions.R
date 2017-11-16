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
