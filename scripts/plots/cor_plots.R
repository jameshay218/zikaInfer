#############################
## Forecasting fit analysis
#############################
library(zikaInfer)
library(ggplot2)
library(extrafont)
library(cowplot)
library(hexbin)
library(RColorBrewer)

chainWD <- "~/net/home/zika/outputs/"

# Heatmaps ----------------------------------------------------------------
heatmap_theme <- theme(axis.text.x=element_text(size=8,family="Arial",color="black"),
                       plot.title=element_text(size=8,family="Arial",color="black"),
                       axis.text.y=element_text(size=8,family="Arial",color="black"),
                       axis.title.x=element_text(size=10,family="Arial",color="black"),
                       axis.title.y=element_text(size=10,family="Arial",color="black"),
                       legend.text=element_text(size=8,family="Arial",color="black"),
                       legend.title=element_text(size=8,family="Arial",color="black"),
                       strip.text=element_text(size=8,family="Arial")) + theme_bw()


# Reading in data ---------------------------------------------------------
bahia_chain <- zikaInfer::load_mcmc_chains(location = paste0(chainWD,"bahia_forecast"),
                                          asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 1,burnin = 50000)
nejm_chain <- zikaInfer::load_mcmc_chains(location = paste0(chainWD,"northeast_forecast"),
                                          asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 1,burnin = 50000)
bahia_chain$propn_increase <- bahia_chain$incPropn2/bahia_chain$incPropn
nejm_chain$propn_increase <- nejm_chain$incPropn2/nejm_chain$incPropn
                                          
bahia_dat<- make_micro_dat(bahia_chain,10000)
bahia_dat$Model <- "Bahia, Brazil (reports)"
nejm_dat<- make_micro_dat(nejm_chain,10000)
nejm_dat$Model <- "Northeast Brazil"
dat <- rbind(bahia_dat, nejm_dat)
dat$Model <- as.factor(dat$Model)

# Microceph risk curves ---------------------------------------------------
p <- ggplot(data=dat) +
  geom_vline(xintercept=c(14*7, 28*7),col="grey40",lty="dashed",size=1) +
  geom_ribbon(aes(x=weeks,ymax=upper,ymin=lower,fill=Model,group=Model),alpha=0.6) +
  geom_line(aes(x=weeks,y=means,group=Model),col="black",lty="dashed") +
  facet_wrap(~Model,ncol=1) +
  ylab("Relative probability of \nmicrocephaly given ZIKV infection")+
  xlab("Week of gestation at time of infection") +
  theme_bw() +
  theme(axis.text.x=element_text(size=8,family="Arial",color="black"),
        panel.grid.minor=element_blank(),
        plot.title=element_text(size=8,family="Arial",color="black"),
        axis.text.y=element_text(size=8,family="Arial",color="black"),
        axis.title.x=element_text(size=10,family="Arial",color="black"),
        axis.title.y=element_text(size=10,family="Arial",color="black"),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=8,family="Arial",color="black"),
        legend.justification=c(0,0)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.05)) +
  scale_x_continuous(expand=c(0,0),limits=c(0,40*7),breaks=seq(0,280,by=7*4), labels=seq(0,40,by=4))
p


# Densities
bahia_chain$version <- "Bahia, Brazil (reports)"
nejm_chain$version <- "Northeast Brazil"

b_chain <- reshape2::melt(bahia_chain, id.vars="version")
n_chain <- reshape2::melt(nejm_chain, id.vars="version")

chains <- rbind(b_chain,n_chain)
chains$version <- as.factor(chains$version)
use_vars <- c("birth_reduction","propn_increase","incPropn","incPropn2","abortion_rate")
tmp_chain <- chains[chains$variable %in% use_vars,]

fake_dat <- data.frame(version=c(rep("Bahia - reports",12),
                                 rep("Northeast Brazil NEJM", 12)),
                       variable=rep(c("propn2","propn2","birth_reduction","birth_reduction","propn_increase",
                                      "propn_increase","incPropn","incPropn","incPropn2","incPropn2",
                                      "abortion_rate","abortion_rate"),2),
                       value=rep(c(0,2,0,1,0,500,0,0.01,0,1,0,1),2)
                                 )
fake_dat <- fake_dat[fake_dat$variable %in% use_vars,]
convert_names <- c("propn2"="Microcephaly reporting rate 2016",
                   "birth_reduction"="Proportion reduction in affected births",
                   "propn_increase"="Increase in ZIKV reporting rate",
                   "incPropn"="ZIKV incidence reporting rate 2015",
                   "incPropn2"="ZIKV incidence reporting rate 2016",
                   "abortion_rate"="Proportion of microcephaly affected births aborted")
tmp_chain$variable <- convert_names[as.character(tmp_chain$variable)]
fake_dat$variable <- convert_names[as.character(fake_dat$variable)]

use_vars2 <- c("birth_reduction","propn_increase","abortion_rate","propn")
convert_names <- c("propn"="Microcephaly reporting\nrate 2015/2016",
                   "birth_reduction"="Proportion reduction\nin affected births",
                   "propn_increase"="Increase in ZIKV reporting rate",
                   "incPropn"="ZIKV incidence reporting rate 2015",
                   "incPropn2"="ZIKV incidence reporting rate 2016",
                   "abortion_rate"="Proportion of microcephaly\naffected births aborted")

tmp_chain <- b_chain[b_chain$variable %in% use_vars2,]
tmp_chain$variable <- convert_names[as.character(tmp_chain$variable)]

quantiles <- reshape2::melt(plyr::ddply(tmp_chain, "variable",function(x) quantile(x$value,c(0.025,0.5,0.975))))
colnames(quantiles) <- c("variable","quantile","value")
densities_bahia <- ggplot(tmp_chain) +
  geom_density(aes(x=value,y=..density..,fill=variable)) +
  geom_vline(data=quantiles[quantiles$quantile %in% c("2.5%","97.5%"),],aes(xintercept=value),col="black",lty="dashed")+
  geom_vline(data=quantiles[quantiles$quantile %in% c("50%"),],aes(xintercept=value),col="black")+
  facet_wrap(~variable,scale="free") +
  heatmap_theme +
  ylab("Density") +
  xlab("Value") +
  theme_bw() +
  theme(legend.position="none")
densities_bahia

tmp_chain <- n_chain[n_chain$variable %in% use_vars2,]
tmp_chain$variable <- convert_names[as.character(tmp_chain$variable)]

quantiles <- reshape2::melt(plyr::ddply(tmp_chain, "variable",function(x) quantile(x$value,c(0.025,0.5,0.975))))
colnames(quantiles) <- c("variable","quantile","value")
densities_nejm <- ggplot(tmp_chain) +
  geom_density(aes(x=value,y=..density..,fill=variable)) +
  geom_vline(data=quantiles[quantiles$quantile %in% c("2.5%","97.5%"),],aes(xintercept=value),col="black",lty="dashed")+
  geom_vline(data=quantiles[quantiles$quantile %in% c("50%"),],aes(xintercept=value),col="black")+
  facet_wrap(~variable,scale="free") +
  heatmap_theme +
  ylab("Density") +
  xlab("Value") +
  theme_bw() +
  theme(legend.position="none")
densities_nejm



# Corrplots Bahia ---------------------------------------------------------------
heatmap_theme <- heatmap_theme + theme(
  axis.text.x=element_text(size=8,family="Arial",color="black"),
  axis.text.y=element_text(size=8,family="Arial",color="black"),
  legend.text=element_text(size=6,family="Arial",color="black")
  )

p1 <- ggplot() + geom_hex(data=bahia_chain,aes(x=birth_reduction,y=propn_increase,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.005) +
  ylab("Relative increase in\n ZIKV reporting rate") +
  xlab("Proportion decrease in\nZIKV affected births") +
  labs(fill="Density") + heatmap_theme
p2 <- ggplot() + geom_hex(data=bahia_chain,aes(x=birth_reduction,y=propn,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.005) +
  ylab("Microcephaly reporting\nrate in first wave") +
  xlab("Proportion decrease in\nZIKV affected births") +
  labs(fill="Density") + heatmap_theme
p3 <- ggplot() + geom_hex(data=bahia_chain,aes(x=birth_reduction,y=abortion_rate,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.005) +
  ylab("Proportion of microcephaly\naffected births aborted") +
  xlab("Proportion decrease in\nZIKV affected births") +
  labs(fill="Density") + heatmap_theme
p4 <- ggplot() + geom_hex(data=bahia_chain,aes(x=propn_increase,y=propn,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.007) +
  ylab("Microcephaly reporting\nrate in first wave") +
  xlab("Relative increase in\n ZIKV reporting rate") +
  labs(fill="Density") + heatmap_theme
p5 <- ggplot() + geom_hex(data=bahia_chain,aes(x=propn_increase,y=abortion_rate,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.006) +
  ylab("Proportion of microcephaly\naffected births aborted") +
  xlab("Relative increase in\nZIKV reporting rate") +
  labs(fill="Density") + heatmap_theme
p6 <- ggplot() + geom_hex(data=bahia_chain,aes(x=propn,y=abortion_rate,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.006) +
  ylab("Proportion of microcephaly\naffected births aborted") +
  xlab("Microcephaly reporting\nrate in first wave") +
  labs(fill="Density") + heatmap_theme
cor_bahia <- plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)
cor_bahia

# Corrplots NEJM ----------------------------------------------------------
p1 <- ggplot() + geom_hex(data=nejm_chain,aes(x=birth_reduction,y=propn_increase,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.005) +
  ylab("Relative increase in\n ZIKV reporting rate") +
  xlab("Proportion decrease in\nZIKV affected births") +
  labs(fill="Density") + heatmap_theme
p2 <- ggplot() + geom_hex(data=nejm_chain,aes(x=birth_reduction,y=propn,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.0035) +
  ylab("Microcephaly reporting\nrate in first wave") +
  xlab("Proportion decrease in\nZIKV affected births") +
  labs(fill="Density") + heatmap_theme
p3 <- ggplot() + geom_hex(data=nejm_chain,aes(x=birth_reduction,y=abortion_rate,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.002) +
  ylab("Proportion of microcephaly\naffected births aborted") +
  xlab("Proportion decrease in\nZIKV affected births") +
  labs(fill="Density") + heatmap_theme
p4 <- ggplot() + geom_hex(data=nejm_chain,aes(x=propn_increase,y=propn,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.0075) +
  ylab("Microcephaly reporting\nrate in first wave") +
  xlab("Relative increase in\n ZIKV reporting rate") +
  labs(fill="Density") + heatmap_theme
p5 <- ggplot() + geom_hex(data=nejm_chain,aes(x=propn_increase,y=abortion_rate,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.003) +
  ylab("Proportion of microcephaly\naffected births aborted") +
  xlab("Relative increase in\nZIKV reporting rate") +
  labs(fill="Density") + heatmap_theme
p6 <- ggplot() + geom_hex(data=nejm_chain,aes(x=propn,y=abortion_rate,fill=..density..))+ 
  scale_fill_gradient2(low="#5E4FA2",mid="#FAFDB8",high="#9E0142",midpoint= 0.0025) +
  ylab("Proportion of microcephaly\naffected births aborted") +
  xlab("Microcephaly reporting\nrate in first wave") +
  labs(fill="Density") + heatmap_theme
cor_nejm <- plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)
cor_nejm
