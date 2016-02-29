scaletuning <- function(step, popt,pcur){
    if(pcur ==1) pcur <- 0.99
    if(pcur == 0) pcur <- 0.01
    step = (step*qnorm(popt/2))/qnorm(pcur/2)
    if(step > 1) step <- 1
    return(step)
}


proposalfunction1 <- function(param,param_table,index){
    # 4th index is upper bound, 3rd is lower
    # 1st and 2nd index used in other functions
    mn <- param_table[index,"lower_bounds"]
    mx <- param_table[index,"upper_bounds"]

    rtn <- param
    
    x <- rtn[index]
    
    x <- toUnitScale(param[index],mn,mx)

    # 5th index is step size
    stp <- param_table[index,"steps"]

    rv <- runif(1)
    rv <- (rv-0.5)*stp
    x <- x + rv

    # Bouncing boundary condition
    if (x < 0) x <- -x
    if (x > 1) x <- 2-x

    # Cyclical boundary conditions
    #if (x < 0) x <- 1 + x	
    #if (x > 1) x <- x - 1
    
    if(x < 0 | x > 1) print("Stepped outside of unit scale. Something went wrong...")

    rtn[index] <- fromUnitScale(x,mn,mx)
    rtn
}


run_metropolis_MCMC <- function(startvalue,
                                iterations=1000,
                                data,
                                ts,
                                y0s,
                                param_table,
                                popt=0.44,
                                opt_freq=50,
                                thin=1,
                                burnin=100,
                                adaptive_period=1,
                                filename,
                                save_block = 500,
                                VERBOSE=FALSE
                                ){
    TUNING_ERROR<- 0.1

    if(opt_freq ==0 && VERBOSE){ print("Not running adaptive MCMC - opt_freq set to 0")}
    else if(VERBOSE){ print("Adaptive MCMC - will adapt step size during specified burnin period")}

    # Set up quicker tables
    # Turns the data frame table into a matrix that can allow faster indexing
    param_transform_table <- as.matrix(param_table[,c("use_log","lower_bounds","upper_bounds","steps","log_proposal")])
    
                                        # Get those parameters which should be optimised
    non_fixed_params <- which(param_table$fixed==0)
    non_fixed_params_length <- length(non_fixed_params)
    all_param_length <- length(startvalue)

    # Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste(filename,"_chain.csv",sep="")
    chain_colnames <- c("sampno",param_table$names,"lnlike")
  
    # Arrays to store acceptance rates
    tempaccepted <- tempiter <- reset <- integer(all_param_length)
    reset[] <- 0

    # Create empty chain to store "save_block" iterations at a time
    empty_chain <- chain <- matrix(nrow=save_block,ncol=all_param_length+2)
    
    # Set starting value and params
    current_params <- startvalue
    
    # Create array to store values
    empty_values <- values <- sample <- numeric(save_block)

    probab <- posterior(ts, y0s, startvalue, data)

    # Set up initial csv file
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1,] <- c(1,startvalue,probab)
    colnames(tmp_table) <- chain_colnames
    
    # Write starting conditions to file
    write.table(tmp_table,file=mcmc_chain_file,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)

    no_recorded <- 1
    sampno <- 2
    
                                        # Go through chain
    for (i in 1:(iterations+adaptive_period+burnin)){
                                        # For each parameter (Gibbs)
        j <- sample(non_fixed_params,1)
     #   for(j in non_fixed_params){
                                        # Propose new parameters and calculate posterior
            #'if(i > adaptive_period + burnin) browser()
        proposal <- proposalfunction1(current_params,param_transform_table,j)
        #print(proposal)
        #proposal <- current_params
        #proposal[j] <- proposal_function(current_params[j],param_transform_table[j,"steps"])
        newprobab <- posterior(ts, y0s, proposal, data)

        #print(newprobab)
                                        # Calculate log difference in posteriors and accept/reject
        difflike <- newprobab - probab
        
        if ((!is.nan(difflike) & !is.infinite(newprobab)) & (runif(1) < exp(difflike) |  difflike > 0)){
            if(proposal[j] < param_transform_table[j,"upper_bounds"] & proposal[j] > param_transform_table[j,"lower_bounds"]){
                current_params <- proposal
            probab <- newprobab
            tempaccepted[j] <- tempaccepted[j] + 1
            }
        }
            tempiter[j] <- tempiter[j] + 1

           # If current iteration matches with recording frequency, store in the chain. If we are at the limit of the save block,
           # save this block of chain to file and reset chain
        if(sampno %% thin ==0){

            chain[no_recorded,1] <- sampno
            chain[no_recorded,2:(ncol(chain)-1)] <- current_params
                chain[no_recorded,ncol(chain)] <- probab
                no_recorded <- no_recorded + 1
                
                if(no_recorded > save_block){
                    print(i)
                    write.table(chain[1:(no_recorded-1),],file=mcmc_chain_file,col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
                    chain <- empty_chain
                    no_recorded <- 1
                }
            }
            sampno <- sampno + 1
#        }
        
                                        # Update step sizes based on acceptance rate
                                        # Note that if opt_freq is 0, then no tuning will take place
        if(opt_freq != 0 & i <= adaptive_period & i%%opt_freq== 0) {
            pcur <- tempaccepted/tempiter
            print(pcur[non_fixed_params])
            tempaccepted <- tempiter <- reset
            tmp_transform <- param_transform_table[,"steps"]
            for(x in non_fixed_params){
                if(pcur[x] < popt - (TUNING_ERROR*popt) | pcur[x] > popt + (TUNING_ERROR*popt)){
                    tmp_transform[x] <- scaletuning(tmp_transform[x],popt,pcur[x])
                    #tmp_transform[x] <- scaletuning2(tmp_transform[x],popt,pcur[x])
                }
            }
            print("Step sizes:")
            print(tmp_transform[non_fixed_params])
            param_transform_table[,"steps"] <- tmp_transform
            
        }
    }

    # If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
    # that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
    # rather than an array
    if(no_recorded > 2){
        write.table(chain[1:(no_recorded-1),],file=mcmc_chain_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
    }
    
    if(VERBOSE){
        print("Final step sizes:")
        print(param_table$step)
    }
    return(mcmc_chain_file)
}
