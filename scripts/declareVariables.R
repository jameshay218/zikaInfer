plot_chain <- function(place, burnin=25000){
  filename <- paste(place,"_chain.csv",sep="")
  greb <- fread(filename,data.table=FALSE)
  greb <- greb[burnin:nrow(greb),c("probMicro","baselineProb","r0","b")]
  plot(as.mcmc(greb))
}

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}





states <- c(
  "Acre"="acre",
  "Alagoas"="alagoas",
  "Amapá"="amapa",
  "Amazonas"="amazonas",
  "Bahia"="bahia",
  "Ceará"="ceara",
  "Distrito Federal"="distritofederal",
  "Espírito Santo"="espiritosanto",      
  "Goiás"="goias",
  "Maranhão"="maranhao",
  "Mato Grosso do Sul"="matogrossdosul",
  "Mato Grosso"="matogrosso",
  "Minas Gerais"="minasgerais",
  "Pará"="para",
  "Paraíba"="paraiba",
  "Paraná"="parana",
  "Pernambuco"="pernambuco",
  "Piauí"="piaui",
  "Rio de Janeiro"="riodejaneiro",
  "Rio Grande do Norte"="riograndedonorte",
  "Rio Grande do Sul"="riograndedosul",
  "Rondônia"="rondonia",
  "Roraima"="roraima",
  "São Paulo"="saopaulo",
  "Santa Catarina"="santacatarina",
  "Sergipe"="sergipe",
  "Tocantins"="tocantins"  
)

states_with_data <- c("pernambuco",
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

omgwow <- function(){
  pars[[1]][1] <- 100
  ys <- NULL
  index <- 1
  for(state in sort(states_with_data)){
    print(state)
    chains <- NULL
   for(i in 1:3) chains[i] <- paste(state,"_",i,"_chain.csv",sep="")
    tmpY <- get_best_trajectory(chains,30000)
    ys[[index]] <- cbind(tmpY[,"times"],rowSums(tmpY[,c("I_F","I_A","I_C")])/rowSums(tmpY[,5:ncol(tmpY)]))
    index <- index + 1
    }
  plot(ys[[1]][,2]~ys[[1]][,1],type='l',col=c25[1],xlab="Time (years)",ylab="Incidence per capita",ylim=c(0,0.15),main="Estimated time to zika resurgence")
  abline(v=ys[[1]][ys[[1]][,1] > 10,][which.max(ys[[1]][ys[[1]][,1] > 10,2]),1],col=c25[1],lty=3)
  for(i in 1:length(ys)){
    lines(ys[[i]][,2]~ys[[i]][,1],type='l',col=c25[i])
    abline(v=ys[[i]][ys[[i]][,1] > 10,][which.max(ys[[i]][ys[[i]][,1] > 10,2]),1],col=c25[i],lty=3)
  }
  legend(75,0.15,sort(states_with_data),
  lty=c(1,1),
  col=c25)
}