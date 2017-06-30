#' @title hcr_recruitment_model
#' 
#' @description Model to predict the recruitment.
#' 
#' NOTE: Only ricker model and vs. one or the other of segreg or bevholt
#' reccv
#' @export
#' 
#' @param ssb The true spawning stock biomass
#' @param reccv XXX
#' @param ctr The control file, containing the parameters
hcr_recruitment_model <- function(ssb,reccv,ctr) 
{
  #nsamp <- ctr$iter
  #fit <- ctr$ssb_pars
  #pR <- t(sapply(seq(nsamp), function(j) exp(match.fun(fit $ model[j]) (fit[j,], ssb)) ))
  #return(pR)
  #hcr_recruitment_model <- function(ssb,ctr) 
  #{
  
  #fit <- ctr$ssb_pars
  #rec <- ifelse(fit$model %in% "segreg",
  #              exp(log(ifelse(ssb >= fit$b,fit$a*fit$b,fit$a*ssb))) * reccv,
  #              exp(log(fit$a) + log(ssb) - fit$b * ssb) * reccv)
  #return(rec)
  ssb <- ssb * 1e6
  fit <- ctr$ssb_pars
  rec <- ifelse(fit$model %in% "segreg",
                exp(log(ifelse(ssb >= fit$b,fit$a*fit$b,fit$a*ssb))) * reccv,
                exp(log(fit$a) + log(ssb) - fit$b * ssb) * reccv)
  return(rec/1e6)
  
}
