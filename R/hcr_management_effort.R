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
