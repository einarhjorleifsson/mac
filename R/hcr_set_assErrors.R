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
