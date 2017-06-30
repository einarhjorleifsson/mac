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
