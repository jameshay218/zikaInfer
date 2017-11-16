library(zikaProj)
library(ggplot2)

## Change this to where google trend data is saved
trends <- read.csv("~/Documents/Zika/Data/brazil/google_trends_brazil.csv")

# Trend plot --------------------------------------------------------------
#event_names <- c("PAHO Epidemiological Alert",
#                 "WHO/PAHO alert",
#                 "Brazil national public health emergency",
#                 "WHO PHEIC")
event_names <- c("A","C","B","D","E","F")
event_dates <- c(as.numeric(as.Date("2015-05-07",origin="2015-01-01")), 
                 as.numeric(as.Date("2015-12-01",origin="2015-01-01")), 
                 as.numeric(as.Date("2015-11-11",origin="2015-01-01")), 
                 as.numeric(as.Date("2016-02-01",origin="2015-01-01")),
                 as.numeric(as.Date("2016-03-02",origin="2015-01-01")),
                 as.numeric(as.Date("2016-08-17",origin="2015-01-01")))

## Shift to fit on plot
adjust <- c(25,-25,25,-25,25,25)
event_dates <- event_dates - as.numeric(as.Date("2015-01-01",origin="2015-01-01"))
events <- data.frame(Date=event_dates,names=event_names,adjust=adjust)


## Clean up trend data
trends <- reshape2::melt(trends,id.vars="Week")
trends$Week <- as.numeric(as.Date(as.character(trends$Week),origin="2015-01-01"))
trends$Week <- as.numeric(trends$Week) - as.numeric(as.Date("2015-01-01",origin="2015-01-01"))
colnames(trends) <- c("Week","Search term","Volume")
convert <- c("zika"="zika","zika.microcefalia"="zika microcefalia","microcefalia"="microcefalia")
trends$`Search term` <- convert[trends$`Search term`]

## Make axis labels
breaks <- cumsum(c(0,rep(getDaysPerMonth(breaks=2),2)))
breaks <- c(breaks,events$Date)
labels <- as.Date(breaks,origin="2015-01-01")
cols <- c(rep("black",5),rep("red",6))

## The actual plot
trend_plot <- ggplot(trends) + 
  geom_line(aes_string(x="Week",y="Volume",col="`Search term`"),size=0.6) + 
  theme_classic() + 
  ylab("Relative search volume") +
  geom_vline(data=events, aes(xintercept=Date),col="red",linetype="dashed",size=0.6) +
  scale_colour_manual(values=c("blue","green","black")) +
  geom_text(data=events,aes(x=Date,y=Inf,label=names),vjust=-0.5,size=4,angle=0,fontface="bold") +
  xlab("") +
  scale_x_continuous(limits=c(0,730),breaks=breaks,labels=labels) +
  scale_y_continuous(expand=c(0,0),limits=c(-1,100)) +
  theme(text=element_text(size=8,family = "Arial"),
        axis.text.y=element_text(size=8,family="Arial"),
        axis.text.x=element_text(colour=cols,size=8,family="Arial",angle=45,hjust=1),
        axis.title=element_text(size=10,family="Arial"),
        legend.position=c(1,1),
        legend.justification=c(1.1,1.1),
        plot.margin=unit(c(0.9,0.5,0,0.5),"lines"),
        panel.background=element_blank())
gt <- ggplot_gtable(ggplot_build(trend_plot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
trend_plot <- gt
plot(trend_plot)
