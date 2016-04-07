library(ggplot2)
greb <- read.csv("~/Dropbox/Zika/results4/finalstats2.csv")
countries <- c("pernambuco","bahia","saopaulo","paraiba","maranhao","ceara","sergipe","riodejaneiro","piaui","riograndenorte","minasgerais","matogrosso","algoas","para","acre","espiritosanto","tocantins")
greb <- greb[greb$local %in% countries,]
greb[grep("r0",greb$X),"X"] <- "r0"
greb[grep("probMicro",greb$X),"X"] <- "probMicro"
greb[grep("baselineProb",greb$X),"X"] <- "baselineProb"
greb[grep("epiStart",greb$X),"X"] <- "epiStart"
#plot(greb[greb$X=="r0","X50."]~greb[greb$X=="probMicro","X50."],ylim=c(0,5))
dat <- greb[greb$X=="r0" | greb$X == "probMicro",c("X","X50.","local")]
colnames(dat) <- c("variable","mean","local")
dat <- data.frame("r0"=dat[dat$variable=="r0","mean"],"probMicro"=dat[dat$variable=="probMicro","mean"],"local"=dat[dat$variable=="r0","local"])
plot <- ggplot(data=dat,aes(x=probMicro,y=r0))+ geom_smooth(method='lm',formula=y~x) + geom_point(aes(colour=local)) + ggtitle("Correlation between mean R0 and mean probMicro")
print(plot)