#' @title Operating model
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param y XXX
#' @param h XXX
#' @param ctr ctr_rec
#' @param Fmult XXX
#' @param nR Not used

hcr_operating_model <- function(y, h, ctr, Fmult, nR=1) {
  
  n_ages  <- dim(X$N)[1]
  n_years <- dim(X$N)[2]
  
  N   <- X$N[,y ,h,]
  cW  <- X$cW[,y,h,]
  sW  <- X$sW[,y,h,]
  xN <- X$mat[,y,h,]
  M   <- X$M[,y,h,]
  pF <- X$pF[,y,h,]
  pM <- X$pM[,y,h,]
  
  tF <- t(Fmult*t(X$selF[,y,h,]))
  dF <- t(Fmult*t(X$selD[,y,h,]))
  X$tF[,y,h,] <<- tF
  
  # Conventional ssb
  ssb <- colSums(N * exp( -(pM * M + pF * (tF + dF))) * xN * sW)
  ## Mean age in the spawning stock
  # mAge <- colSums((SSBay*c(1:n_ages)))/ssb 
  ## Egg mass
  # ssb <- colSums(Ny*exp(-(My*pMy+(Fy+Dy)*pFy))*maty*sWy * (0.005*sWy))
  
  
  X$N[1,y,h,] <<- hcr_recruitment_model(ssb = ssb, reccv=X$cvR[y,h, ] ,ctr = ctr)
  
  N[1,] <- X$N[1,y,h,] # update the recruits in the current year
  # TAKE THE CATCH
  X$C[,y,h,] <<- N * tF/ (tF + dF + M + 0.00001)*(1-exp(- (tF + dF + M)))
  # NEXT YEARS STOCK
  
  NyEnd <- N * exp( -( tF + dF + M))
  if(y < n_years ) {
    X$N[2:n_ages, y+1, h,] <<- NyEnd[1:(n_ages-1),]
    # plus group calculation
    X$N[n_ages,y+1,h,] <<- X$N[n_ages,y+1,h,] +
      X$N[n_ages,y,h,] * exp(-X$tF[n_ages,y,h,] - X$M[n_ages,y,h,])
  }
  #return(d)
}
