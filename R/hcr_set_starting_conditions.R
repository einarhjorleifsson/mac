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
