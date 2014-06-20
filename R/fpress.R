#' @title HCR - Read control file
#'
#' @description Reads a simple text file and returns a list
#' 
#' @export
#' 
#' @param file filename

hcr_read_ctr <- function (file) 
{
  x <- read.table(file)
  ctr <- list()
  for (i in 1:nrow(x)) ctr[[i]] <- x[i,]
  names(ctr) <- c("a1","a2","y1","y2","iter","f1","f2","nR","tac_y1","tac_y2",
                  "y1Bias","r_cv","r_rho","r_model","r_mean","ssb_break",
                  "a_cv","a_rho","a_error","a_bias","w_cv","w_rho","w_error",
                  "w_refB","h_alpha","h_beta","b_trigger","delay","h_number",
                  "i_number","b2","b2")
  return(ctr)
}


#' @title hcr_set_dimensions
#' 
#' @description Sets upp the dimensions and the error structure for the simulation.
#' 
#' NOTE: The current script assumes that the cv and rho in the stock are the same
#' as in the catches. 
#' @param ctr A list containing the various dimentions as well as cv's for the
#' some of the parameters.
hcr_set_dimensions <- function(ctr) {
  
  d <- list()
  
  HRATE <- ctr$HRATE
  
  nR <- ctr$nR
  
  a1 <- ctr$a1
  a2 <- ctr$a2
  n_ages <- a2 - a1 + 1
  y1 <- ctr$y1
  y2 <- ctr$y2
  n_years <- y2 - y1 + 1
  
  iter <- ctr$iter
  
  # array with age as a dimention
  x <- array(-1,dim=c(n_ages,n_years,length(HRATE),iter),
             dimnames=list(age=a1:a2,
                           year=y1:y2,
                           hrate=HRATE,
                           iter=1:iter))
  d$N <- d$tF <- d$M <- d$pM <- d$pF <- d$selF <- d$selD <- d$selB <- d$C <- d$cW  <- 
    d$sW <- d$cvcW <- d$cvsW <- d$mat  <- x
  
  # arrays without an age dimention
  x <- array(-1,dim=c(n_years,length(HRATE),iter),
             dimnames=list(year=y1:y2,
                           hrate=HRATE,
                           iter=1:iter))
  d$TAC <- d$assError <- x
  
  # array for recruit is based on year classes
  startYC <- y1-nR     # most recent year class with measurement
  nYC <- y2-startYC+1    # number of year classes that need to be
  # simulated
  x <- array(-1,dim=c(nYC,length(HRATE),iter),
             dimnames=list(yearclass=startYC:y2,
                           hrate=HRATE,
                           iter=1:iter))
  d$cvR  <- x
  
  #
  d$cvcW <- hcr_set_wgtErrors(d$cvcW,ctr)
  
  d$cvsW <- d$cvcW
  #
  d$cvR <- hcr_set_recErrors(d$cvR,ctr)
  #
  d$assError <- hcr_set_assErrors(d$assError,ctr)
  
  #X <<- d
  return(d)
}

#' @title HCR: Setup of assessment errors
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param d XXX
#' @param ctr XXX
#' 
hcr_set_assErrors <- function(d,ctr) 
  {

  n_years  <- dim(d)[1]
  n_hrates <- dim(d)[2]
  n_iters  <- dim(d)[3]

  x <- array(rnorm(n_years  * n_iters),
             dim=c(n_years ,  n_iters))

  for (i in 2:n_years) x[i,] <- ctr$a_rho * x[i-1,] + sqrt(1-ctr$a_rho^2) * x[i,]
  x <- ctr$a_cv * x
  
  # take a subset of samples, ensures that there is potential a bias
  # in the assessment year (does not matter if looking at long term
  # equilibrium). This this is not an issue, can ignore coding like this.
  
  # k <- 100:(n_years + 1000 - 100)    # ignore the 1st 100
  # firstSample <- sample(k,1)
  # lastSample  <- firstSample + n_years -1
  # x <- x[firstSample:lastSample,]
  for (h in 1:n_hrates) d[,h,] <- x
  
  # CHECK THIS:
  ## ad hoc error setup in the 1st year, fixed for iCod age range
  #d$N[(nR+1):n_ages,1,,] <- (1/ctr$Year1Bias)*N1[(nR+1):n_ages]*rep(exp(d$assError[1,,]),14)
  return(d)
}

#' @title hcr_set_wgtErrors
#' 
#' @description Sets up the weight error structure
#' 
#' @param d XXX
#' @param ctr A list with control values, here use the weight cv and weight rho.
hcr_set_wgtErrors <- function(d,ctr) 
{
  # Weight error - note for first year
  #                no cv in this implementation, need to be added
  # NOTE - same error is applied to all ages - should include an option
  #        for white noise accross ages.
  #        This is now updated in the fishvise.R script within the mac-project
  n_ages <- dim(d)[1]
  n_years  <- dim(d)[2]
  n_hrates <- dim(d)[3]
  n_iters  <- dim(d)[4]
  
  x <- array(rnorm(n_years * n_iters),
             dim=c(n_years,  n_iters))
  
  for (y in 2:n_years) x[y,] <- ctr$w_rho * x[y-1,] + sqrt(1 - ctr$w_rho^2) * x[y,]
  
  for (a in 1:n_ages) {
    for (h in 1:n_hrates) d[a,,h,] <- x * ctr$cW_cv[a]
  }
  
  return(d)
  
}


#' @title HCR: Setup of recruitment error structure
#' 
#' @description This function is modified from the original in the fishvise package.
#' In that package there was only one parameter estimate (a, b and cv). Here we have 
#' different parameters for each iteration, and even different recruitment model.
#' 
#' @export
#' 
#' @param d XXX
#' @param ctr XXX
#' 
hcr_set_recErrors <- function(d,ctr) 
{
  
  n_years  <- dim(d)[1]
  n_hrates <- dim(d)[2]
  n_iters  <- dim(d)[3]
  
  # Recruitment: cv & autocorrelation
  x <- array(rnorm(n_years * n_iters),
             dim=c(n_years , n_iters))
  
  for (y in 2:n_years) x[y,] <- ctr$r_rho * x[y-1,] + sqrt(1 - ctr$r_rho^2) * x[y,]
  
  x <- exp(x*ctr$ssb_pars$cv)
  for (h in 1:n_hrates) d[,h,] <- x
  
  return(d)
}

#' @title HCR: Reading starting years input from file
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param file XXX
#' 
hcr_read_startfile <- function(file) {
  
  d <- list()
  
  indata    <- matrix(scan(file,0,quiet=TRUE),ncol=18,byrow=TRUE)
  d$age     <- indata[,1]     # age classes
  d$N       <- indata[,2]     # population (000s)
  d$N_cv    <- indata[,3]     # population cv
  d$sW      <- indata[,4]     # spawning weight (kg)
  d$sW_cv   <- indata[,5]     # spawning weight cv
  d$cW      <- indata[,6]     # catch weight (kg)
  d$cW_cv   <- indata[,7]     # catch weight cb
  d$mat     <- indata[,8]     # maturity
  d$mat_cv  <- indata[,9]     # maturity cv
  d$selF    <- indata[,10]    # fishing mortality (selection pattern)
  d$selF_cv <- indata[,11]    # fishing mortality cv
  d$pF      <- indata[,12]    # proportion of fishing mort. bf. spawning
  d$selD    <- indata[,13]    # discard mortality
  d$selD_cv <- indata[,14]    # discard mortality cv
  d$M       <- indata[,15]    # natural mortality
  d$M_cv    <- indata[,16]    # natural mortality cv
  d$pM      <- indata[,17]    # proportio of natural mort. bf. spawning
  d$selB    <- indata[,18]    # selection pattern of the "fishable biomass"
  
  return(d)
}

#' @title HCR: Set starting condition for stock
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param dat_y1 XXX
#' @param d XXX
#' @param ctr XXX
#' 

hcr_set_starting_conditions <- function(dat_y1, d, ctr) 
  {
  
  d$N[,1,,]   <- dat_y1$N
  d$sW[,,,]   <- dat_y1$sW   # spawning weight (kg)
  d$cW[,,,]   <- dat_y1$cW   # catch weight (kg)
  d$mat[,,,]  <- dat_y1$mat  # maturity
  d$selF[,,,] <- dat_y1$selF # fishing mortality (selection pattern)
  d$pF[,,,]   <- dat_y1$pF   # proportion of fishing mort. bf. spawning
  d$selD[,,,] <- dat_y1$selD # discard mortality
  d$M[,,,]    <- dat_y1$M    # natural mortality
  d$pM[,,,]   <- dat_y1$pM   # proportio of natural mort. bf. spawning
  d$selB[,,,] <- dat_y1$selB # selection pattern of the "fishable biomass"
  
  n_ages <- length(dat_y1$age)
  n_noRec <- sum(dat_y1$N == 0)
  ## ad hoc error setup in the 1st year
  #d$N[(n_noRec+1):n_ages,1,,] <- (1 / ctr$y1Bias) * d$N[(n_noRec+1):n_ages] * rep(exp(d$assError[1,,]),n_ages-1)
    
  d$TAC[1,,]  <- ctr$tac_y1 # Already set TAC in the assessment year (year 1)
  d$TAC[2,,]  <- ctr$tac_y2 # Already set TAC in the advisory year (year 2)
  
  n_years <- dim(d$cW)[2]
  if(ctr$w_error == 1) {
    d$cW[,2:n_years,,] <- d$cW[,2:n_years,,] * exp(d$cvcW[,2:n_years,,])
    d$sW[,2:n_years,,] <- d$sW[,2:n_years,,] * exp(d$cvsW[,2:n_years,,])
  }
  
  if(ctr$w_error == 2) {
    d$cW[,2:n_years,,] <- d$cW[,2:n_years,,] * (1 + d$cvcW[,2:n_years,,])
    d$sW[,2:n_years,,] <- d$sW[,2:n_years,,] * (1 + d$cvsW[,2:n_years,,])
  }
  
  if(ctr$w_refB == 0) d$bW <- d$sW   # use stock weights to calculate ref bio
  if(ctr$w_refB == 1) d$bW <- d$cW   # use catch weights to calculate ref bio
  
  #X <<- list()
  X <<- d # Pass as a global variable - this stuff will be updated in the loop
  #return(d)
}




#' @title hcr_TAC_to_Fmult
#' 
#' @description XXX
#' 
#' @export 
#' 
#' @param y XXX
#' @param h XXX
#' 
hcr_TAC_to_Fmult <- function(y,h) {

  TAC <-  X$TAC[y,h,]
  Na  <-  X$N[,y,h,]
  Sa  <-  X$selF[,y,h,]
  Da  <-  X$selD[,y,h,]
  Ma  <-  X$M[,y,h,]
  Wa  <-  X$cW[,y,h,]
  
  epsilon <- 1e-04
  Ba <- Na * Wa

  B <- colSums(Ba)

  TAC <- ifelse(TAC > 0.9 * B, 0.9 * B, TAC)
  Fmult <- TAC/colSums(Ba * Sa * exp(-Ma)) + 0.05
  for (i in 1:5) {
    Fa <- t(Fmult * t(Sa))
    Za <- Fa + Ma + epsilon  #added on iSaithe, but why worked on iCod?
    Y1 <- colSums(Fa/Za * (1 - exp(-Za)) * Ba)
    Fa <- t((Fmult + epsilon) * t(Sa))
    Za <- Fa + Ma + epsilon #added on iSaithe, but why worked on iCod?
    Y2 <- colSums(Fa/Za * (1 - exp(-Za)) * Ba)
    dY <- (Y2 - Y1)/epsilon
    Fmult <- Fmult - (Y2 - TAC)/dY
  }
  Fmult <- ifelse(TAC == 0, 0, Fmult)
  Fmult <- ifelse(Fmult < 0, 0, Fmult)
  Fmult <- ifelse(Fmult > 1.5,1.5,Fmult)  # Quick fix
  return(Fmult)
}






#' @title Add implementation error to TAC
#' 
#' @description Dummy function
#' 
#' @export
#' 

hcr_implementation_model_1 <- function () #(TACy2, TACy1, beta1, beta2) 
  {
  #i <- TACy2 < TACy1
  #if (sum(i) > 0) 
  #  TACy2[i] <- TACy2[i] * (TACy1[i]/TACy2[i])^(rbeta(1,beta1, beta2))
  #return(TACy2)
}

#' @title Operating model
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param y XXX
#' @param h XXX
#' @param ctr ctr_rec
#' @param Fmult XXX
#' @param nR XXX

hcr_operating_model <- function(y, h, ctr, Fmult, nR=1) {
  
  n_ages  <- dim(X$N)[1]
  n_years <- dim(X$N)[2]
  
  N   <- X$N[,y ,h,]
  cW  <- X$cW[,y,h,]
  sW  <- X$sW[,y,h,]
  xN <- X$mat[,y,h,]
  M   <- X$M[,y,h,]
  pF <- X$pF[,y,h,]
  pM <- X$pM[,y,h,]
  
  tF <- t(Fmult*t(X$selF[,y,h,]))
  dF <- t(Fmult*t(X$selD[,y,h,]))
  X$tF[,y,h,] <<- tF
  
  # Conventional ssb
  ssb <- colSums(N * exp( -(pM * M + pF * (tF + dF))) * xN * sW)
  ## Mean age in the spawning stock
  # mAge <- colSums((SSBay*c(1:n_ages)))/ssb 
  ## Egg mass
  # ssb <- colSums(Ny*exp(-(My*pMy+(Fy+Dy)*pFy))*maty*sWy * (0.005*sWy))
  
  
  X$N[1,y,h,] <<- hcr_recruitment_model(ssb = ssb, reccv=X$cvR[y,h, ] ,ctr = ctr)
  
  N[1,] <- X$N[1,y,h,] # update the recruits in the current year
  # TAKE THE CATCH
  X$C[,y,h,] <<- N * tF/ (tF + dF + M + 0.00001)*(1-exp(- (tF + dF + M)))
  # NEXT YEARS STOCK

  NyEnd <- N * exp( -( tF + dF + M))
  if(y < n_years ) {
    X$N[2:n_ages, y+1, h,] <<- NyEnd[1:(n_ages-1),]
    # plus group calculation
    X$N[n_ages,y+1,h,] <<- X$N[n_ages,y+1,h,] +
      X$N[n_ages,y,h,] * exp(-X$tF[n_ages,y,h,] - X$M[n_ages,y,h,])
  }
  #return(d)
}


#' @title hcr_recruitment_model
#' 
#' @description Model to predict the recruitment.
#' 
#' NOTE: Only ricker model and vs. one or the other of segreg or bevholt
#' reccv
#' @export
#' 
#' @param ssb The true spawning stock biomass
#' @param reccv XXX
#' @param ctr The control file, containing the parameters
hcr_recruitment_model <- function(ssb,reccv,ctr) 
{
  #nsamp <- ctr$iter
  #fit <- ctr$ssb_pars
  #pR <- t(sapply(seq(nsamp), function(j) exp(match.fun(fit $ model[j]) (fit[j,], ssb)) ))
  #return(pR)
  #hcr_recruitment_model <- function(ssb,ctr) 
  #{
  
  #fit <- ctr$ssb_pars
  #rec <- ifelse(fit$model %in% "segreg",
  #              exp(log(ifelse(ssb >= fit$b,fit$a*fit$b,fit$a*ssb))) * reccv,
  #              exp(log(fit$a) + log(ssb) - fit$b * ssb) * reccv)
  #return(rec)
  ssb <- ssb * 1e6
  fit <- ctr$ssb_pars
  rec <- ifelse(fit$model %in% "segreg",
                exp(log(ifelse(ssb >= fit$b,fit$a*fit$b,fit$a*ssb))) * reccv,
                exp(log(fit$a) + log(ssb) - fit$b * ssb) * reccv)
  return(rec/1e6)
  
}



#' @title Observation error model
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param y XXX
#' @param h XXX
#' @param Fmult XXX
#' @param ctr XXX

hcr_observation_error <- function(y, h, Fmult,ctr) {
  
  hrate <- ctr$HRATE[h]
  delay <- ctr$delay
  
  assError <- X$assError[ y + delay, h,]
  
  N   <-    X$N[,y + delay, h,]
  bW  <-   X$bW[,y + delay, h,]
  sB  <- X$selB[,y + delay, h,]
  bio <- colSums(N * bW * sB)
  
  sW  <-   X$sW[,y + delay ,h,]
  xN  <-  X$mat[,y + delay, h,]
  M   <-    X$M[,y + delay, h,]
  pF  <-   X$pF[,y + delay, h,]
  pM  <-   X$pM[,y + delay, h,]
  selF <- X$selF[, y + delay, h,]
  selD <- X$selD[, y + delay, h,]
  
  totalF <- t(Fmult * t(selF))
  
  dF <- t(Fmult * t(selD))
  
  ssb <- colSums(N * exp( -( pM * M + pF * (totalF + dF))) * xN * sW)
  
  n_iters <- ctr$iter
  hrate    <- rep(hrate, n_iters)
    
  ## A. The assessment error model
  if (ctr$a_error == 1) {
    bio_hat     <- bio   * ctr$a_bias * exp(assError)
    ssb_hat     <- ssb   * ctr$a_bias * exp(assError)
    hrate_hat   <- hrate * ctr$a_bias * exp(assError)
  }
  
  if (ctr$a_error == 2) {
    bio_hat     <- bio   * ctr$a_bias * (1 + assError)
    ssb_hat     <- ssb   * ctr$a_bias * (1 + assError)
    hrate_hat   <- hrate * ctr$a_bias * (1 + assError)
  }
  
  return(list(hrate=hrate_hat,bio=bio_hat,ssb=ssb_hat))
}


#' @title Effort type management
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param hrate_hat XXX
#' @param ssb_hat XXX
#' @param Btrigger XXX
hcr_management_effort <- function(hrate_hat,ssb_hat,Btrigger) 
  {
  i <- ssb_hat < Btrigger  
  hrate_hat[i] <- hrate_hat[i] * ssb_hat[i]/Btrigger
  
  return(hrate_hat)
}

#' @title F-based Harvest Control Rule
#' 
#' @description The F-based rule is the conventional ICES decision rule. Here it
#' is implemented such that the TAC next year is calculated from the true
#' stock in numbers based on a fishing mortality that includes observation error.
#' 
#' If the Btrigger is set in the rule (Btrigger > 0) then linear reductions of 
#' fishing mortality is done relative to observed spawning stock biomass (i.e.
#' that includes observation errrors).
#' 
#' @export
#' 
#' @param y XXX
#' @param h XXX
#' @param hrate Harvest rate - with error
#' @param ssb Spawning stock biomass - with error
#' @param ctr Control file
#' @note Need to check is ssb-hat is calculated according to the correct delay
#' specification. 
#'

hcr_management_fmort <- function(y,h,hrate,ssb,ctr)
{
  selF     <- X$selF[,y + ctr$delay,h,]
  selD     <- X$selD[,y + ctr$delay,h,]
  Na       <- X$N[,y + ctr$delay,h,]
  cWa      <- X$cW[,y + ctr$delay,h,]
  selF     <- X$selF[,y + ctr$delay,h,]
  selD     <- X$selD[,y + ctr$delay,h,]
  Ma       <- X$M[,y + ctr$delay,h,] 
  
  # adjust harvest rate
  i <- ssb < ctr$b_trigger  
  hrate[i] <- hrate[i] * ssb[i]/ctr$b_trigger  
  
  
  Fa <- t(hrate * t(selF))
  Da <- t(hrate * t(selD))
  
  tac_next_year <- colSums(Na * Fa/(Fa + Da + Ma + 1e-05) * 
                             (1 - exp(-(Fa + Da + Ma))) * cWa)
  
  # 20% TAC constraint
  #tac_next_year <- ctr$h_alpha * tac_this_year + (1 - ctr$h_alpha) * tac_next_year
  
  # TAC buffer
  if(ctr$h_alpha > 0) {
    tac_this_year <- X$TAC[y,h,] # This years TAC
    i <- tac_next_year > tac_this_year * (1 + ctr$h_alpha)
    if(any(i)) tac_next_year[i] <- tac_this_year[i] * (1 + ctr$h_alpha)
    i <- tac_next_year < tac_this_year * (1 - ctr$h_alpha)
    if(any(i)) tac_next_year[i] <- tac_this_year[i] * (1 - ctr$h_alpha)
  }
  
  
  X$TAC[y+1,h,] <<- tac_next_year
  
  if(ctr$h_beta > 0) {
    # now need to add in the constraint that the F does not deviate by more than
    #  10% from the target
    # Calculate the current Fmultiplier
    xF <- hcr_TAC_to_Fmult(y+1,h)
    # Now find out which of the xF deviate by more than 0.9 and 1.1 of the hrate
    i <- xF > hrate * (1 + ctr$h_beta)
    if(any(i)) hrate[i] <- hrate[i] * (1 + ctr$h_beta)
    i <- xF < hrate * (1 - ctr$h_beta)
    if(any(i)) hrate[i] <- hrate[i] * (1 - ctr$h_beta)
    # now we have adjusted the realized harvest rate the 10% F-constraint
    # Need to update the whole calucation for the TAC
    
    # NEED TO DOUBLE CHECK NEXT STEP, THIS WAS ALREADY FIDDLED WITH ABOVE
    #  No need to worry if NO Btrigger
    #i <- ssb < ctr$b_trigger  
    #hrate[i] <- hrate[i] * ssb[i]/ctr$b_trigger  
    Fa <- t(hrate * t(selF))
    Da <- t(hrate * t(selD))
    
    tac_next_year <- colSums(Na * Fa/(Fa + Da + Ma + 1e-05) * 
                               (1 - exp(-(Fa + Da + Ma))) * cWa)
    X$TAC[y+1,h,] <<- tac_next_year
  }
}

#' @title Biomass-based Harvest Control Rule
#' 
#' @description The Biomass-based HCR is used in the case of the Icelandic cod
#' and saithe. Here it is implemented such that the TAC next year is a multiplier
#' of the target harvest rate and the observed reference biomass (i.e. that 
#' includes observation errrors).
#' 
#' If the Btrigger is set in the rule (Btrigger > 0) then linear reductions of 
#' fishing mortality is done relative to observed spawning stock biomass (i.e.
#' that includes observation errrors).
#' 
#' @export
#' 
#' @param y XXX
#' @param h XXX
#' @param bio vector Observed reference biomass.
#' @param ssb vector Observed spawning stock biomass.
#' @param ctr XXX
#' 
#' 
#' @note To do: Modify function so that buffer is not active below Btrigger and
#' also to a EU type TAC-constraint.
#' 
hcr_management_bio <- function(y,h,bio,ssb,ctr)
  {
  hrate <- rep(ctr$HRATE[h],ctr$iter)
  Btrigger <- ctr$b_trigger
  tac_this_year <- X$TAC[y,h,] # This years TAC
  
  i <- ssb < Btrigger
  hrate[i] <- hrate[i] * ssb[i]/Btrigger
  
  tac_next_year <- hrate * bio  # Next years TAC
  
  # Consider buffer
  tac_next_year <- ctr$h_alpha * tac_this_year + (1 - ctr$h_alpha) * tac_next_year
  
  X$TAC[y+1,h,] <<- tac_next_year

}

#' @title hcr_summarise_data
#' 
#' @description A wrapper function for summarising the simulation results for
#' a given target rate
#' 
#' @export
#' 
#' @param TARGET The target harvest rate to display for the dynamic results
#' @param YEARS The year range over which to summarise the equilibrium results
#' @param ctr Control file
hcr_summarise_data <- function(TARGET=0.22,YEARS=c(2037:2039,ctr=ctr)) {
  sY <- melt(colSums(X$C * X$cW))
  sS <- melt(colSums(X$N * exp(- (X$pM * X$M + X$pF * X$tF)) * X$sW * X$mat))
  sB <- melt(colSums(X$N * X$bW * X$selB))
  R <- melt(drop(X$N[1,,,]))
  Fbar <- melt(colMeans(X$tF[(ctr$f1+1):(ctr$f2+1),,,]))
  d <- data.frame(year=sY$year,iter=sY$iter,target=sY$hrate,
                  r=R$value,
                  bio=sB$value,
                  ssb=sS$value,
                  tF=Fbar$value,
                  hr=sY$value/sB$value,
                  yield=sY$value)
  i <- d$year %in% YEARS
  cn <- c("year","target","r","ssb","tF","yield")
  eq <- d[i,cn]
  eq <- melt(eq,id.vars = c("year","target"))
  eq <- ddply(eq,c("target","variable"),summarise,
              q05=quantile(value,0.05),
              q10=quantile(value,0.10),
              q25=quantile(value,0.25),
              q50=quantile(value,0.50),
              q75=quantile(value,0.75),
              q90=quantile(value,0.90),
              q95=quantile(value,0.95),
              m=mean(value))
  eq_plot <- ggplot(eq,aes(target)) + 
    geom_ribbon(aes(ymin=q05,ymax=q95),fill="red",alpha=0.2) +
    geom_ribbon(aes(ymin=q10,ymax=q90),fill="red",alpha=0.2) +
    geom_ribbon(aes(ymin=q25,ymax=q75),fill="red",alpha=0.2) +
    geom_line(aes(y=q50),col="red") +
    geom_line(aes(y=m),col="blue") +
    facet_wrap(~ variable,scales="free_y")
  i <- d$target %in% TARGET
  cn <- c("year","r","ssb","tF","yield")
  dyn <- d[i,cn]
  dyn <- melt(dyn,id.vars = c("year"))
  dyn <- ddply(dyn,c("year","variable"),summarise,
               q05=quantile(value,0.05),
               q10=quantile(value,0.10),
               q25=quantile(value,0.25),
               q50=quantile(value,0.50),
               q75=quantile(value,0.75),
               q90=quantile(value,0.90),
               q95=quantile(value,0.95),
               m=mean(value))
  dyn_plot <- ggplot(dyn,aes(year)) + 
    geom_ribbon(aes(ymin=q05,ymax=q95),fill="red",alpha=0.2) +
    geom_ribbon(aes(ymin=q10,ymax=q90),fill="red",alpha=0.2) +
    geom_ribbon(aes(ymin=q25,ymax=q75),fill="red",alpha=0.2) +
    geom_line(aes(y=q50),col="red") +
    geom_line(aes(y=m),col="blue") +
    facet_wrap(~ variable,scales="free_y")
  return(list(eq=eq,eq_plot=eq_plot,dyn=dyn,dyn_plot=dyn_plot,ctr=ctr))
}