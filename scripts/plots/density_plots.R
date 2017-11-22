library(ggplot2)
library(gridExtra)
library(extrafont)
library(zikaProj)

## See scripts/utility/combine_all_chains.R
## This should be a reshape2::melt chain combining all of the runs.
meltedChain <- data.table::fread("~/Documents/Zika/combinedChains.csv",data.table=FALSE)

parNames <- c("maxWeek"="Peak risk week",
              "r0"="R0",
              "t0"="Epidemic seed time",
              "lower"="First risk week","upper"="Last risk week",
              "range"="Number of risk weeks",
              "tr1"="Mean 1st trimester risk","tr2"="Mean 2nd trimester risk",
              "tr3"="Mean 3rd trimester risk",
              "baselineProb"="Baseline microcephaly proportion")


####################
## Microceph parameters
####################
## Risk timings
tmpChain <- meltedChain[meltedChain$variable %in% c("maxWeek",
                                                    "lower","upper","range") & 
                          meltedChain$version == "model_1",]
tmpChain$variable <- parNames[tmpChain$variable]
tmpChain$value <- tmpChain$value/7
tmpChain$runName <- as.factor(tmpChain$runName)
p2 <- ggplot(tmpChain) + geom_violin(aes(x=runName,fill=runName,y=value),draw_quantiles=c(0.5),col="black",scale="width") +
  facet_wrap(~variable) + 
  theme_bw() +
  theme(axis.title.y=element_text(size=10,family="Arial"),
        axis.text.y=element_text(size=8,family="Arial"),
        strip.text=element_text(size=10,family="Arial"),
        axis.title.x=element_blank(),legend.position="none",
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,family="Arial",size=8)) +
  geom_hline(yintercept=c(14,28),linetype="dashed",col="red",alpha=0.5) +
  ylab("Estimated value (weeks)")


## Risk values
tmpChain1 <- meltedChain[meltedChain$variable %in% c("tr1","tr2","tr3","baselineProb") & 
                           meltedChain$version == "model_1",]
tmpChain1$variable <- parNames[tmpChain1$variable]
tmpChain1$runName <- as.factor(tmpChain1$runName)
p1 <- ggplot(tmpChain1) + geom_violin(aes(x=runName,fill=runName,y=value),draw_quantiles=c(0.5),col="black",scale="width") +
  facet_wrap(~variable,scales="free_y") + 
  theme_bw() +
  theme(axis.title.y=element_text(size=10,family="Arial"),
        axis.text.y=element_text(size=8,family="Arial"),
        strip.text=element_text(size=10,family="Arial"),
        axis.title.x=element_blank(),legend.position="none",
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1,family="Arial",size=8)) +
  ylab("Estimated value")


cairo_ps("risk_densities.eps",width=7,height=6,family="Arial")
p1
dev.off()


cairo_ps("week_densities.eps",width=7,height=6,family="Arial")
p2
dev.off()


####################
## R0 plot
####################
tmpChain <- meltedChain[meltedChain$variable =="R0" & 
                           meltedChain$version == "model_1",]

convertState <- c("bahia"="Bahia, Brazil",
                  "colombia"="Colombia",
                  "northeast"="Northeast Brazil",
                  "riograndedonorte"="Rio Grande do Norte, Brazil",
                  "pernambuco"="Pernambuco")
tmpChain$state <- convertState[tmpChain$state]
p3 <- ggplot(tmpChain) +
    geom_violin(aes(x=runName, y=value, fill=runName),scale="width",draw_quantiles=c(0.5),col="black") +
    facet_wrap(~state) +
    theme_bw() +
    geom_hline(yintercept=1,col="red",linetype="dashed")+
    scale_y_continuous(limits=c(0,5))+
    theme(axis.title.y=element_text(size=10,family="Arial"),
          axis.text.y=element_text(size=8,family="Arial"),
          strip.text=element_text(size=10,family="Arial"),
          axis.title.x=element_blank(),legend.position="none",
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0,family="Arial",size=8))

cairo_ps("r0_densities.eps",width=7,height=6,family="Arial")
p3
dev.off()
