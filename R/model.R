#' Simulates head circumference data from an SEIR vector borne model
#'
#' Takes a list of necessary arguments to solve a system of ODEs that generates head circumferences of new born infants
#' @param allPars A list of parameter vectors. First vector should be an array of time points over which to solve the SEIR model. Second vector should be the vector of starting conditions. Third vector should be all of the parameters.
#' @return a matrix of simulated head circumferences over time.
#' @export
#' @useDynLib zikaProj
zika.sim <- function(allPars){
  ts <- allPars[[1]]
  y0s <- allPars[[2]]
  pars <- allPars[[3]]
  sampFreq <- pars[1]
  sampPropn <- pars[2]
  mu_I <- pars[3]
  sd_I <- pars[4]
  mu_N <- pars[5]
  sd_N <- pars[6]
  probMicro <- pars[7]
  
  pars <- allPars[[3]][8:length(allPars[[3]])]
  y <- ode(y0s,ts,func="derivs",parms=pars,dllname="zikaProj",initfunc="initmod",maxsteps=100000,atol=1e-10,reltol=1e-10,hmax=1e-4,nout=0)

                                         # y <- lsoda(y0s,ts,zika.ode,pars)
  y <- as.data.frame(y)
  colnames(y) <- c("times","Sm","Em","Im","Sc","Sa","Sf","Ec","Ea","Ef","Ic","Ia","If","Rc","Ra","Rf","RatePregnantI","RateInfected","RatePregnantAll","CumInc")
                  #,"leavingIf","recoverIf","allBirths","inc")
  daysPerYear <- nrow(y)/max(y$times)
  y <- y[y$times > pars[11],]

                                        #  return(y)
  tmp <- numeric(ceiling(max(y$times)*daysPerYear/sampFreq))
  tmpN <- numeric(ceiling(max(y$times)*daysPerYear/sampFreq))
  tmp[1] <- y[1,"If"]
  tmpN[1] <- sum(y[1,c("Sf","Ef","If","Rf")])
  
  i <- 1 + sampFreq
  index <- 1
  all <- NULL
  birthsPerYear <- sum(y0s[4:6])/pars[3]
  birthsPerDay <- ceiling(birthsPerYear/daysPerYear)
  
  while(i <= nrow(y)){
#      tmp[index] <- y[(i-sampFreq):i, "If"]
      is <- y[(i-sampFreq):i, "If"]
      ns <- rowSums(y[(i-sampFreq):i,c("Sf","Ef","If","Rf")])
      propns <- mean(na.omit(is/ns))
#      tmp[index] <- y$IfA[i] - y$IfA[i-sampFreq]
 #     tmpN[index] <- y$fB[i] - y$fB[i-sampFreq]
      #tmpN[index] <- sum(y[(i-sampFreq):i,c("Sf","Ef","If","Rf")])
      
#      alpha_I <- probMicro*tmp[index]/tmpN[index]
      alpha_I <- probMicro * propns
      
      N <- sampPropn*(birthsPerDay*sampFreq)
#      N <- 13500/(12*4) * sampPropn
#      print(N)
      components <- sample(1:2,c(alpha_I,1-alpha_I),size=N,replace=TRUE)
      mus <- c(mu_I,mu_N)
      sds <- c(sd_I,sd_N)
      
      distribution <- round(rnorm(n=N,mean=mus[components],sd=sds[components]),digits=1)
      all[[index]] <- distribution
      index <- index + 1
      i <- i + sampFreq
  }
  return(as.matrix(rbind.fill(lapply(all,function(x) as.data.frame(t(x))))))
}

#' R0 calculation
#'
#' Calculates the R0 of the SEIR model given a vector of parameters and human/mosquito population sizes. R0 defined as number of expected human cases given introduction of 1 infected human into a totally naive population of humans and mosquitoes.
#' @param params Vector of parameters matching those returned by \code{\link{setupListPars}}
#' @param NH integer value for human population size
#' @param NM integer value for mosquito population size
#' @return A single value for R0
#' @export
#' @seealso \code{\link{b.calc}}
r0.calc <- function(params,NH,NM){
    pars <- params[8:length(params)]
    muM <- 1/pars[1]
    sigmaM <- 1/pars[2]

    muH <- 1/pars[3]
    gammaH <- 1/pars[7]

    b <- pars[8]
    pHM <- pars[9]
    pMH <- pars[10]
    
        
    #'    first <- (b * pHM * NM)/((gammaH + muH)*NH)
    #'  second <- (sigmaM/(sigmaM+muM))*(b * pMH )/muM

    #' R0 <- first*second
    R0 <- (b^2*pHM*pMH*NM*sigmaM)/((sigmaM+muM)*muM*(gammaH+muH)*NH)
    return(R0)
}

#' Bite rate calculation
#'
#' Calculates the bite rate needed to generate a given R0 value, assuming that all other parameters are fixed.
#' @param params Vector of parameters matching those returned by \code{\link{setupListPars}}
#' @param NH integer value for human population size
#' @param NM integer value for mosquito population size
#' @param R0 desired R0 value
#' @return A single value for bite rate
#' @export
#' @seealso \code{\link{r0.calc}}
b.calc <- function(params,NH,NM,R0){    pars <- params[8:length(params)]
    muM <- 1/pars[1]
    sigmaM <- 1/pars[2]

    muH <- 1/pars[3]
    gammaH <- 1/pars[7]

    pHM <- pars[9]
    pMH <- pars[10]

    b <- sqrt(((sigmaM+muM)*muM*(gammaH+muH)*NH)/((pHM * pMH * NM * sigmaM))*R0)
    return(b)
    
}

#' Defaults ODE parameters
#'
#' Generates the default vector of parameters for the ODE model
#' @return vector of named parameter values
#' @export
#' @useDynLib zikaProj
setupParsODE <- function(){
    D_EM=10.5/360
    D_EH=4/360
    D_IH=5/360
    L_M=14/360
    L=73.6
    D_C=18
    D_F=3/12
    b=100
    p_MH=0.5
    p_HM=0.5
    t_seed <- 2
    I0 <- 10
    offset <- 0
    
    pars <- c("L_M"=L_M,"D_EM"=D_EM,"L_H"=L,"D_C"=D_C,"D_F"=D_F,"D_EH"=D_EH,"D_IH"=D_IH,"b"=b,"p_HM"=p_HM,"p_MH"=p_MH,"t_seed"=t_seed,"I0"=I0,"constSeed"=offset)
    return(pars)
}

#' Default simulation parameters
#'
#' Generates the default vector of parameters for the simulation
#' @return vector of names parameter values
#' @export
#' @seealso \code{\link{setupParsODE}}
#' @useDynLib zikaProj
setupParsLong <- function(){
    pars <- setupParsODE()
    sampFreq <- 7
    sampPropn <- 0.9
    mu_I <- 28
    sd_I <- 2
    mu_N <- 35
    sd_N <- 2
    probMicro <- 1
    pars <- c("sampFreq"=sampFreq,"sampPropn"=sampPropn,"mu_I"=mu_I,"sd_I"=sd_I,"mu_N"=mu_N,"sd_N"=sd_N,"probMicro"=probMicro,pars)
    return(pars)   
}

#' Default parameter list setup
#'
#' Sets up the default list of parameters and initial conditions to be used by \code{\link{zika.sim}}
#' @param duration Duration in years of the simulation. Defaults to 10 years
#' @return list of times, initial conditions, and ODE parameters
#' @export
#' @seealso \code{\link{setupParsLong}}
#' @useDynLib zikaProj
setupListPars <- function(duration=10){
    pars <- setupParsLong()
    D_C <- pars["D_C"]
    L <- pars["L_H"]
    N_M <- 3000000
    N_H <- 1000000
    S_M = 1*(N_M)
    E_M = 0
    I_M = 0
    S_F = 0.001*N_H
    S_C = D_C*N_H/(L+D_C)
    S_A = (1-0.001-D_C/(L+D_C))*N_H
    E_C = 0
    E_A = 0
    E_F = 0
    I_C = 0
    I_A = 0
    I_F = 0
    R_C = 0
    R_A = 0
    R_F = 0

    y0 <- c(S_M, E_M,I_M,S_C,S_A,S_F,E_C,E_A,E_F,I_C,I_A,I_F,R_C,R_A,R_F, 0, 0, 0,0)
    ts <- 0:(360*duration) / 360
    
    return(list(ts,y0,unname(pars)))
}

#' Obsolete - the R implementation of the ODE model
#'
#' Set of ODEs for the zika SEIR model to be used by deSolve
#' @param t the vector of times to be solved over
#' @param y the compartment states
#' @param pars the parameters used to solve the ODE model
#' @return set of derivatives at the given time point
#' @export
#' @useDynLib zikaProj
zika.ode <-function(t, y,pars){
  L_M <- pars[1]
  D_EM <- pars[2]
   
  L_H <- pars[3]
  D_C <- pars[4]
  D_F <- pars[5]
  D_EH <- pars[6]
  D_IH <- pars[7]
  B_H <- L_H - D_C
 
  b <- pars[8]
  P_HM <- pars[9]
  P_MH <- pars[10]
  
  t_seed <- pars[11]
  I0 <- pars[12]
  epsilon <- 0.0001
  seed <- 0
  offset <- 0
  
  if(t >= t_seed){
      offset <- pars[13]
      if(t < t_seed + epsilon){
          seed <- I0
      }
  }
  
  S_M <- y[1]
  E_M <- y[2]
  I_M <- y[3]
  S_C <- y[4]
  S_A <- y[5]
  S_F <- y[6]
  
  R_C <- y[13]
  R_A <- y[14]
  R_F <- y[15]
  
  N_H <- sum(y[4:15])
  N_M <- sum(y[1:3])
  
  lambda_M <- b*(N_M/N_H)*P_HM*(I_C + I_A + I_F + offset + seed)/N_H
  lambda_H <- b*(N_M/N_H)*P_MH*I_M/N_H
  
  dS_M <- N_M/L_M - S_M/L_M - lambda_M*S_M
  dE_M <- lambda_M*S_M - E_M/L_M - E_M/D_EM
  dI_M <- E_M/D_EM - I_M/L_M
  
  dS_C <- N_H/L_H - S_C/L_H - S_C/D_C - lambda_H*S_C
  dS_A <- S_C/D_C - S_A/L_H - S_A/B_H + S_F/D_F - lambda_H*S_A
  dS_F <- S_A/B_H - S_F/D_F - S_F/L_H - lambda_H*S_F
  
  dE_C <- lambda_H*S_C - E_C/L_H           - E_C/D_EH           - E_C/D_C
  dE_A <- lambda_H*S_A - E_A/L_H - E_A/B_H - E_A/D_EH + E_F/D_F + E_C/D_C
  dE_F <- lambda_H*S_F - E_F/L_H + E_A/B_H - E_F/D_EH - E_F/D_F
  
  dI_C <- E_C/D_EH - I_C/L_H - I_C/D_IH           - I_C/D_C
  dI_A <- E_A/D_EH - I_A/L_H - I_A/D_IH - I_A/B_H + I_C/D_C + I_F/D_F
  dI_F <- E_F/D_EH - I_F/L_H - I_F/D_IH + I_A/B_H           - I_F/D_F 
  
  dR_C <- I_C/D_IH - R_C/L_H - R_C/D_C
  dR_A <- I_A/D_IH + R_C/D_C - R_A/L_H - R_A/B_H + R_F/D_F
  dR_F <- I_F/D_IH + R_A/B_H - R_F/D_F - R_F/L_H
  
  dIf_A <- I_F/D_F
  df_B <- S_F/D_F + E_F/D_F + R_F/D_F + I_F/D_F
  dE <- E_C/D_EH + E_A/D_EH + E_F/D_EH
  
  return(list(c(dS_M,dE_M,dI_M,dS_C,dS_A,dS_F,dE_C,dE_A,dE_F,dI_C,dI_A,dI_F,dR_C,dR_A,dR_F,dIf_A,df_B,dE)))
}
