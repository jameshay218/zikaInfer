## Requires the combined and melted MCMC chains
## See scripts/utility/combine_all_chains.R

res <- data.table::fread("~/Documents/Zika/results_20112017/all_chains.csv",stringsAsFactors=FALSE,data.table=FALSE)
res <- res[res$runName %in% c("Colombia (suspected)","Northeast Brazil, Lancet","Bahia, Brazil (reports)","Rio Grande do Norte, Brazil (reports)",
                              "Colombia (confirmed)", "Pernambuco (confirmed, reports)","Rio Grande do Norte, Brazil (confirmed)","Salvador, Bahia") &
          res$version == "model_1",]

vars <- c("maxWeek","lower","upper","range","tr1","tr2","tr3","R0","t0","AR")
week_vars <- c("maxWeek","lower","upper","range")
varNames <- c("maxWeek"="Peak risk week","lower"="First risk week","upper"="Last risk week","range"="Number of weeks at risk",
              "tr1"="Mean first trimester risk","tr2"="Mean second trimester risk","tr3"="Mean third trimester risk",
              "R0"="Basic reproductive number, R0","t0"="Epidemic seed time","AR"="Attack rate")
res <- res[res$variable %in% vars,]


final <- plyr::ddply(res,c("runName","variable"),function(x){
   
  quants <- signif(quantile(x$value, c(0.025,0.975)),3)
  tmpmean <- signif(mean(x$value),3)
     if(unique(x$variable) %in% week_vars){
      quants <- signif(quants/7,3)
      tmpmean <- signif(tmpmean/7,3)
  }
  if(x$variable == "t0"){
    quants <- as.Date(quants,origin="2013-01-01")
    tmpmean <- as.Date(tmpmean, origin="2013-01-01")
  }
  
  paste0(tmpmean," (",quants[1],"-",quants[2],")")
}
)
final$variable <- varNames[final$variable]
final <- final[order(final$variable,final$runName),]
colnames(final) <- c("Analysis","Parameter","Posterior mean (95% CI)")

write.table(final,"TableS1.csv",sep=",",row.names=FALSE)
