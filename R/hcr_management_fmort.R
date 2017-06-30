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

hcr_management_fmort <- function(y, h, hrate,ssb,ctr)
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
