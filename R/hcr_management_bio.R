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
