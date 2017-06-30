#' @title B-based Harvest Control Rule
#' 
#' @description The B-based rule is the conventional iCod decision rule. Here it
#' is implemented such that the TAC next year is calculated from the observed
#' reference biomass.
#' 
#' If the Btrigger is set in the rule (Btrigger > 0) then linear reductions of 
#' fishing mortality is done relative to observed spawning stock biomass (i.e.
#' that includes observation errrors).
#' 
#' NOTE: The model currently only works on delay=0, i.e. we are basing the advise
#'       on the biomass in the assessment year.
#'       
#' NOTE: The 10% contraint on the harvest rate not implemented yet.
#' 
#' @export
#' 
#' @param y XXX
#' @param h XXX
#' @param ctr Control file
#' @note XXX
#'

hcr_management_bio_pred <- function(y, h, ctr)
{
  hrate <- rep(ctr$HRATE[h],ctr$iter)
  delay <- ctr$delay
  #delay <- 0  # NOTE - CURRENTLY HARD WIRED
  # SHORT TERM PREDICTION MODEL
  # All values in this block are actually true numbers
  # Use last three years
  # Special case: Need to take into account that the short term mean
  # can only be counted calculated in simulation year 3.
  # I.e.:
  #       when y=2, only have one historical year
  #       when y=3 have two years of history
  #       when y>3 have three years of history
  y1 <- ifelse(y < 4,1,y-3)
  
  bWa  <- X$bW[,c(y1:(y-1)), h,]
  cWa  <- X$cW[,c(y1:(y-1)), h,]
  sWa  <- X$sW[,c(y1:(y-1)), h,]
  
  # Below is currently fixed
  selF <- X$selF[,y + delay, h,]
  selB <- X$selB[,y + delay,h,] 
  selD <- X$selD[,y + delay, h,]   # Not used for mackerel
  Ma   <- X$M[,y + delay, h,]      # Actually a constant
  pF   <- X$pF[,y + delay, h,]
  pM   <- X$pM[,y + delay, h,]
  mat  <- X$mat[,y + delay, h,]
  
  # Only take mean if there are more than 1 year
  if(y > 2) {
    bWa  <- apply(bWa,c(1,3),mean)
    cWa  <- apply(cWa,c(1,3),mean)
    sWa  <- apply(sWa,c(1,3),mean)
    
    # Currently fixed in the implementation
    #selF <- apply(selF,c(1,3),mean)
    #selD <- apply(selD,c(1,3),mean)
    #Ma   <- apply(Ma,c(1,3),mean)
    #pF   <- apply(pF,c(1,3),mean)
    #pM   <- apply(pM,c(1,3),mean)
    #mat  <- apply(mat,c(1,3),mean)
  }
  
  # True stock in numbers
  Na <- X$N[,y + ctr$delay,h,] 
  Fa <- t(hrate * t(selF))
  
  # RECRUITS short term predictions
  # 1. Geometric mean assumptions applied to age 0 in the assessment year
  Na[1,] <- ctr$r_mean
  # 2. RCT3 emulation used on age 1
  #rct3_alpha <- 0.8  # the GM weight in the rct3
  #Na[2,] <- rct3_alpha * ctr$r_mean * exp(-(X$tF[1,y + delay-1,h,]+X$M[1,y + delay-1,h,])) + 
  #  (1-rct3_alpha) * X$N[1,y + ctr$delay -1,h,]
  
  # 2. Geometric mean assumptions applied to age 1 if the delay is 1, i.e.
  #    if we are using the stock in numbers in the advisory year
  #    Since pF (and pM) may vary need to use the short term prediction values
  #if(ctr$delay > 0) {    
  #  Na[2,] <- ctr$r_mean * exp(-(Fa[1,]+Ma[1,]))
  #}
  
  # Predicted biomass based only on true short term mean and recruitment assumption
  bio <- colSums(Na * bWa * selB)
  ssb <- colSums(Na * exp( -( pM * Ma + pF * Fa)) * mat * sWa)
  
  # ASSESSMENT ERROR
  assError <- X$assError[y + ctr$delay, h,]
  
  bio_hat   <- bio   * ctr$a_bias * exp(assError)
  ssb_hat   <- ssb   * ctr$a_bias * exp(assError)
  hrate <- hrate * ctr$a_bias * exp(assError)
  
  # adjust harvest rate according to the trigger
  i <- ssb_hat < ctr$b_trigger  # Note this operates on the ssb with error
  hrate[i] <- hrate[i] * ssb_hat[i]/ctr$b_trigger 
  
  #Fa_hat <- t(hrate_hat * t(selF))
  #Da_hat <- t(hrate_hat * t(selD))
  
  tac_next_year <- hrate * ssb_hat  # Next years TAC
  
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
  
  #if(ctr$h_beta > 0) {
  # now need to add in the constraint that the F does not deviate by more than
  #  10% from the target
  # Calculate the current Fmultiplier
  #  xF <- hcr_TAC_to_Fmult(y+1,h)
  # Now find out which of the xF deviate by more than 0.9 and 1.1 of the hrate
  #  i <- xF > hrate_hat * (1 + ctr$h_beta)
  #  if(any(i)) hrate_hat[i] <- hrate_hat[i] * (1 + ctr$h_beta)
  #  i <- xF < hrate_hat * (1 - ctr$h_beta)
  #  if(any(i)) hrate_hat[i] <- hrate_hat[i] * (1 - ctr$h_beta)
  # now we have adjusted the realized harvest rate the 10% F-constraint
  # Need to update the whole calucation for the TAC
  
  # NEED TO DOUBLE CHECK NEXT STEP, THIS WAS ALREADY FIDDLED WITH ABOVE
  #  No need to worry if NO Btrigger
  #i <- ssb < ctr$b_trigger  
  #hrate[i] <- hrate[i] * ssb[i]/ctr$b_trigger  
  #  Fa_hat <- t(hrate_hat * t(selF))
  #  Da_hat <- t(hrate_hat * t(selD))
  
  #  tac_next_year <- colSums(Na * Fa_hat/(Fa_hat + Da_hat + Ma + 1e-05) * 
  #                             (1 - exp(-(Fa_hat + Da_hat + Ma))) * cWa)
  #  X$TAC[y+1,h,] <<- tac_next_year
  #}
}
