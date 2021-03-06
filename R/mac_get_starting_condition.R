#' @title mac_get_starting_condition
#' 
#' @description Calculates mean's and cv's for weights, maturity, mortality,
#' partial F and M and selection pattern from rbya-objects.
#' 
#' @export
#' 
#' @param rbya data.frame
#' @param years_mean vector containing the years over which to take the mean
#' @param years_cv vector containing the years over which to take the cv

mac_get_starting_condition <- function(rbya,years_mean,years_cv) 
{
  cn <- c("year","age","sW","cW","mat","m","pf","pm","sel")
  
  d_mean <- rbya[rbya$year %in% years_mean,cn]
  names(d_mean) <- c("year","age","sW","cW","mat","M","pF","pM","selF")
  d_mean <- melt(d_mean,c("year","age"))
  d_mean <- ddply(d_mean,c("variable","age"),
                  summarise,
                  value=mean(value))
  
  d_cv <- rbya[rbya$year %in% years_cv,cn]
  names(d_cv) <- c("year","age","sW_cv","cW_cv","mat_cv","M_cv","pF_cv","pM_cv","selF_cv")
  d_cv <- melt(d_cv,c("year","age"))
  d_cv <- ddply(d_cv,c("variable","age"),
                summarise,
                value=sd(value)/mean(value))
  dat <- cbind(dcast(d_mean,age ~ variable,value.var="value"),dcast(d_cv,age ~ variable,value.var="value"))
  dat <- dat[,sort(names(dat))]
  i <- dat$age == 0
  dat$mat_cv[i] <- dat$sW_cv[i] <- dat$pF_cv[i] <- 0
  dat$selD <- 0
  dat$selB <- c(0,0,rep(1,11))
  dat$N <- rbya$n[rbya$year == 2013]
  return(dat)
}
