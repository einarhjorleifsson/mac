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
