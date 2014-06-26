#' @title ypr
#' 
#' @description A yield per recruit model
#' 
#' @export
#' 
#' @param fmort A vector containing the fishing mortality over which to run the ypr
#' @param cW Catch weights
#' @param selC Selection pattern of the catch
#' @param ages A vector containing the age classes
#' @param plusgroup Boolean If FALSE (default) no plusgroup bookkeeping
#' @param a1 A value containing lower reference age
#' @param a2 A value containing upper reference age
#' @param dW Discard weights
#' @param ssbW Spawning stock weights
#' @param selL Selection pattern of the landings
#' @param M Natural mortality, can be a single value
#' @param mat A vector of proportion mature
#' @param pM Partial natural mortality prior to spawning
#' @param pF Partial fishing mortality prior to spawning
#' @param ssb A boolean (Not used)
#' 
ypr <- function(fmort, cW, selC, ages=3:14,plusgroup=F,a1=5,a2=10, dW, ssbW, selL, M,mat,pM=0,pF=0,ssb=F){
  
  if(missing(cW)) stop("Catch weights need to be specified")
  if(missing(selC))   stop("Selection at age for catches must be specified")
  if(missing(ages)) stop("Age vector must be specified")
  
  if(missing(dW))    dW <- rep(0,length(ages))
  if(missing(ssbW))  ssbW <- cW
  
  if(missing(mat))    mat <- rep(1,length(ages)) 
  if(length(M) == 1)  M <- rep(M,length(ages))
  if(length(pM) == 1) pM <- rep(pM,length(ages))
  if(length(pF) == 1) pF <- rep(pF,length(ages))

  if(missing(selL))   selL <- selC
  
  age <- a1:a2
  i <- !is.na(match(ages,age))
  
  selL <- selL/mean(selC[i]) 
  selC <- selC/mean(selC[i]) 
  
  res <- data.frame(fmort=fmort,yield=rep(0,length(fmort)))
  res$SSB <- res$discardN <-res$discardWt <- res$yield
  for(i in 1:nrow(res)) {
    tmp <- ypr_single_target(res$fmort[i], cW, selC,selL,M,ages=ages,plusgroup=plusgroup,mat=mat,pM=pM,pF=pF,dW=dW,ssbW=ssbW)
    res$yield[i] <- tmp$Catch
    res$SSB[i] <- tmp$SSB
    res$discardN[i] <- tmp$Discard
    res$discardWt[i] <- tmp$Discardwt
    
  }
  x <- spline(res$fmort,res$yield,n=10*nrow(res))
  maxy <- max(x$y)
  fmax <- x$x[x$y==max(x$y)]
  d <- diff(x$y)
  d <- d/d[1]
  d1 <- d[1:(length(d)-1)]
  d2 <- d[2:length(d)]
  i <- 1:length(d1)
  i <- i[d1 > 0.1 & d2 < 0.1]
  f01 <- x$x[i]
  
  x <- spline(res$fmort,res$SSB,n=10*nrow(res))
  x1  <- abs(x$y-0.35*max(res$SSB))
  fssb35 <- x$x[x1==min(x1)]
  ssb35 <- x$y[x1==min(x1)]
  refs <- data.frame(f01=f01,fmax=fmax[1],fssb35=fssb35,ssb35=ssb35,maxy=maxy)
  #referencept <- rep(0,4)
  #names(referencept) <- c("f01","fmax","maxy","ssb35")
  #referencept["f01"] <- f01
  #referencept["fmax"] <- fmax[1]
  #referencept["maxy"] <- maxy
  #referencept["ssb35"] <- ssb35
  
  #attributes(res)$refpt <- referencept
  return(list(res=res,refs=refs))
}

ypr_single_target <- function(fmort, cW, meanselcatch,meanselland,M,ages=c(2:9),plusgroup=plusgroup,mat,pM=0,pF=0,dW,ssbW)
{
  if(missing(mat)) mat <- rep(1,length(ages)) 
  if( missing(ssbW)) ssbW <- cW
  if(length(M) == 1) M <- rep(M,length(ages))
  if(length(pM) == 1) pM <- rep(pM,length(ages))
  if(length(pF) == 1) pF <- rep(pF,length(ages))
  
  Fishmort <- meanselcatch * fmort
  Fishmortlan <- meanselland * fmort
  Fishmortdisc <-  Fishmort-Fishmortlan
  Z <- Fishmort + M
  N <- 1000000.
  Catch <- SSB <-Discard <- Discardwt <-   0
  for(i in 1:length(ages)) {
    SSB <- SSB + N*exp(-(Fishmort[i]*pF[i]+M[i]*pM[i]))*ssbW[i]*mat[i]
    Catch <- Catch + Fishmortlan[i]/Z[i] * N * (1 - exp( - Z[i])) *
      cW[i]
    Discard <- Discard +  Fishmortdisc[i]/Z[i] * N * (1 - exp( - Z[i])) # Numbers
    Discardwt <- Discardwt +  Fishmortdisc[i]/Z[i] * N * (1 - exp( - Z[i]))*dW[i] # Biomass   
    N <- N * exp( - Z[i])
  }
  if(plusgroup) { # 10 umfer?ir enn ekki brottkast ? plusgruppu
    n <- length(ages)
    for(i in 1:10) {
      Catch <- Catch + Fishmortlan[n]/Z[n] * N * (1 - exp( - Z[n])) *  # ??ur var sm? villa h?r cW[i]
        cW[n]
      SSB <- SSB + N*exp(-(Fishmort[n]*pF[n]+M[n]*pM[n]))*ssbW[n]*mat[n]
      
      N <- N * exp( - Z[n])
    }
  }
  return(list(Catch=Catch/1e6,SSB=SSB/1e6,Discard=Discard/1e6,Discardwt=Discardwt/1e6))
}


