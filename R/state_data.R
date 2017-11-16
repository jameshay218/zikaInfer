#' Location conversion
#'
#' Converts lower-case simplified names to actual strings
#' @param name the vector or scalar of lower case location names
#' @return the same vector but with actual strings
#' @export
convert_name_to_location <- function(name){
    data(locationInfo)
    locations <- NULL
    for(i in name){
        location <- i
        if(i %in% locationInfo$rawName){
            location <- as.character(locationInfo[locationInfo$rawName == i, "fullName"])
        }
        locations <- c(locations, location)
    }
    return(locations)
}

#' Location order
#'
#' Gets the correct location order for all plots, sorted by number of microcephaly cases
#' @param birth_dat the data frame of data
#' @param named if TRUE, converts the location names to full names
#' @param normalised if TRUE, orders by per birth microcephaly rate
#' @param filePath the full file locatoin of the microcephaly data
#' @return a vector of ordered location names
#' @export
get_correct_order <- function(named=TRUE,normalised=FALSE, filePath="~/net/home/zika/data/microcephaly_data.csv"){
    birth_dat <- read.csv(filePath,stringsAsFactors=FALSE)
    if(named) birth_dat$local <- convert_name_to_location(birth_dat$local)
    if(normalised){
        order <- names(sort(sapply(unique(birth_dat$local),function(x) sum(birth_dat[birth_dat$local==x,"microCeph"]/birth_dat[birth_dat$local==x,"births"])),decreasing = TRUE))
    } else {
        order <- names(sort(sapply(unique(birth_dat$local),function(x) sum(birth_dat[birth_dat$local==x,"microCeph"])),decreasing = TRUE))
    }
    return(order)
}


#' Read all data
#'
#' Function to read in microcephaly and incidence data. Converts strings to factors, and orders them by microcephaly cases
#' @param microCephFile full filepath to microcephaly data
#' @param incDatFile full filepath ot incidence data
#' @param epidemicLocations if TRUE, gets only epidemic locations. Otherwise gives all locations
#' @param order if TRUE, orders the locations by per capita microcephaly rate
#' @param convert if TRUE, converts strings to factors
#' @return a list of microcephaly and incidence data
#' @export
read_all_data <- function(microCephFile="~/net/home/zika/data/microcephaly_data.csv",
                          incDatFile="~/net/home/zika/data/inc_data.csv",
                          epidemicLocations=FALSE, order=TRUE, convert=TRUE){
    ## Read in microcephaly data
    allDat <- read.csv(microCephFile,stringsAsFactors=FALSE)

    ## Read in incidence data
    incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)

    if(epidemicLocations){
        locations <- get_epidemic_locations()
        allDat <- allDat[allDat$local %in% locations$include,]
        incDat <- incDat[incDat$local %in% locations$include,]
    }

    ## Convert to characters and factors
    if(convert){
        allDat$local <- as.character(allDat$local)
        allDat$local <- convert_name_to_location_factor(allDat$local)
        incDat$local <- as.character(incDat$local)
        incDat$local <- convert_name_to_location_factor(incDat$local)
    }
    if(order){
        order <- get_correct_order(named=convert,normalised=TRUE)
        allDat$local <- factor(allDat$local,levels=order)
        incDat$local <- factor(incDat$local,levels=order)
    }
    
    return(list(microDat=allDat,incDat=incDat))
}


#' Convert name to factor
#'
#' Function to convert location strings to factors (with correct accents etc)
#' @param locations the vector of lower case location names to convert
#' @return the converted locations with proper naming and factored
#' @export
convert_name_to_location_factor <- function(locations){
    data(locationInfo)
    order <- locationInfo$incidenceOrder
    country_names <- locationInfo$fullName
    final <- factor(locationInfo[locationInfo$rawName %in% locations, "fullName"], levels=locationInfo[order,"fullName"])
    return(final)
}

#' Subset epidemic locations
#'
#' Gets those locations that fulfill the epidemic criteria of 2sd above mean since July for 3 months
#' @param dat the data frame of data to look through (microcephaly or incidence)
#' @param micro if TRUE, looks at microcephaly data. Incidence data otherwise
#' @param lim in days, the start of the period of consideration
#' @return a list with a vector of epidemic and non-epidemic locations
#' @export
get_epidemic_locations <- function(dat,micro=TRUE,lim=575){
    ## If looking at microcephaly data
    if(micro){
        include <- NULL
        exclude <- NULL
        ## For each state, see if 3 consecutive measurements >2sd above mean
        for(local in unique(dat$local)){
            tmpDat <- dat[dat$local==local,c("startDay","microCeph")]
            tmp <- tmpDat[tmpDat$startDay < lim,2]
            mean <- mean(tmp)
            sd <- sd(tmp)
            a <- tmpDat[tmpDat$startDay >= lim,2] >= (mean + 2*sd)
            if(any(cumul_true(a) >= 3)){
                include <- c(include, local)
            } else {
                exclude <- c(exclude, local)
            }
        }
        ## If looking at ZIKV incidence data
    } else  {
        include <- NULL
        exclude <- NULL
        for(local in unique(dat$local)){
            tmpDat <- dat[dat$local==local,c("startDay","inc")]
            tmp <- tmpDat[tmpDat$startDay < lim,2]
            mean <- mean(tmp)
            sd <- sd(tmp)
            a <- tmpDat[tmpDat$startDay >= lim,2] >= (mean + 2*sd)
            if(any(cumul_true(a) >= 3)){
                include <- c(include, local)
            } else {
                exclude <- c(exclude, local)
            }
        }
    }
    return(list("include"=include,"exclude"=exclude))
}
