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
