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
  Fmult <- ifelse(Fmult > 3,3,Fmult)  # NOTE: This was put to 1.5 in earlier
  #       versions (prior to 22.06.2014).
  # This is not equivalent to setting a 
  # max on the fishing mortality, like maxF=3
  return(Fmult)
}
