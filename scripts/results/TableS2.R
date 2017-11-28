## Requires the combined and melted MCMC chains
## See scripts/utility/combine_all_chains_forecasting.R

res <- data.table::fread("~/Documents/Zika/results_20112017/combinedChains_forecasts.csv",stringsAsFactors=FALSE,data.table=FALSE)

vars <- c("birth_reduction","propn","abortion_rate","zikv_report_rate_change","abortions","avoided_births")
varNames <- c("birth_reduction"="Proportion of pregnant women avoiding infection","propn"="Microcephaly reporting rate 2015","abortion_rate"="Abortion rate",
              "zikv_report_rate_change"="ZIKV reporting rate change","abortions"="Number of aborted births","avoided_births"="Births avoided")
res <- res[res$variable %in% vars,]

final <- plyr::ddply(res,c("runName","variable"),function(x){
  quants <- signif(quantile(x$value, c(0.025,0.975)),3)
  tmpmean <- signif(mean(x$value),3)
  paste0(tmpmean," (",quants[1],"-",quants[2],")")
}
)
final$variable <- varNames[final$variable]
final <- final[order(final$variable,final$runName),]
colnames(final) <- c("Analysis","Parameter","Posterior mean (95% CI)")

write.table(final,"TableS2.csv",sep=",",row.names=FALSE)