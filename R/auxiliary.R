
createCounts <- function(dat){
    all <- NULL
    for(i in 1:nrow(dat)){
        tmp <- as.data.frame(table(dat[i,]))
        tmp[,2] <- tmp[,2]/sum(tmp[,2])
        tmp <- cbind(tmp, i)
        all <- rbind(all, tmp)
    }
    return(all)
}
