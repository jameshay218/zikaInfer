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

c25 <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")
overallTheme <- theme_bw() 
mycolours <-  scale_colour_manual(values=c25)
mycolours2 <- scale_fill_manual(values=c25)

mytheme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.ticks.length=unit(0.3, "lines"),
  axis.ticks.margin=unit(0.5, "lines"),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  legend.background=element_rect(fill="white", colour=NA),
  legend.key=element_rect(colour="white"),
  legend.key.size=unit(1.5, "lines"),
  legend.position="right",
  legend.text=element_text(size=rel(1.2)),
  legend.title=element_text(size=rel(1.4), face="bold", hjust=0),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.margin=unit(0, "lines"),
  plot.background=element_rect(fill="white"),
  plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"),
  plot.title=element_text(size=rel(1.8), face="bold", hjust=0.5),
  strip.background=element_rect(fill="grey90", colour="grey50"),
  strip.text.x=element_text(size=rel(0.8)),
  strip.text.y=element_text(size=rel(0.8), angle=-90) 
)


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