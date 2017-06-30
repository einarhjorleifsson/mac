#' @title hcr_summarise_data
#' 
#' @description A wrapper function for summarising the simulation results for
#' a given target rate
#' 
#' @export
#' 
#' @param TARGET The target harvest rate to display for the dynamic results
#' @param YEARS The year range over which to summarise the equilibrium results
#' @param ctr Control file
hcr_summarise_data <- function(TARGET=0.22,YEARS=c(2037:2039,ctr=ctr)) {
  sY <- melt(colSums(X$C * X$cW))
  sS <- melt(colSums(X$N * exp(- (X$pM * X$M + X$pF * X$tF)) * X$sW * X$mat))
  sB <- melt(colSums(X$N * X$bW * X$selB))
  R <- melt(drop(X$N[1,,,]))
  Fbar <- melt(colMeans(X$tF[(ctr$f1+1):(ctr$f2+1),,,]))
  d <- data.frame(year=sY$year,iter=sY$iter,target=sY$hrate,
                  r=R$value,
                  bio=sB$value,
                  ssb=sS$value,
                  tF=Fbar$value,
                  hr=sY$value/sB$value,
                  yield=sY$value)
  i <- d$year %in% YEARS
  cn <- c("year","target","r","ssb","tF","yield")
  eq <- d[i,cn]
  eq <- melt(eq,id.vars = c("year","target"))
  eq <- ddply(eq,c("target","variable"),summarise,
              q05=quantile(value,0.05),
              q10=quantile(value,0.10),
              q25=quantile(value,0.25),
              q50=quantile(value,0.50),
              q75=quantile(value,0.75),
              q90=quantile(value,0.90),
              q95=quantile(value,0.95),
              m=mean(value))
  eq_plot <- ggplot(eq,aes(target)) + 
    geom_ribbon(aes(ymin=q05,ymax=q95),fill="red",alpha=0.2) +
    geom_ribbon(aes(ymin=q10,ymax=q90),fill="red",alpha=0.2) +
    geom_ribbon(aes(ymin=q25,ymax=q75),fill="red",alpha=0.2) +
    geom_line(aes(y=q50),col="red") +
    geom_line(aes(y=m),col="blue") +
    facet_wrap(~ variable,scales="free_y")
  i <- d$target %in% TARGET
  cn <- c("year","r","ssb","tF","yield")
  dyn <- d[i,cn]
  dyn <- melt(dyn,id.vars = c("year"))
  dyn <- ddply(dyn,c("year","variable"),summarise,
               q05=quantile(value,0.05),
               q10=quantile(value,0.10),
               q25=quantile(value,0.25),
               q50=quantile(value,0.50),
               q75=quantile(value,0.75),
               q90=quantile(value,0.90),
               q95=quantile(value,0.95),
               m=mean(value))
  dyn_plot <- ggplot(dyn,aes(year)) + 
    geom_ribbon(aes(ymin=q05,ymax=q95),fill="red",alpha=0.2) +
    geom_ribbon(aes(ymin=q10,ymax=q90),fill="red",alpha=0.2) +
    geom_ribbon(aes(ymin=q25,ymax=q75),fill="red",alpha=0.2) +
    geom_line(aes(y=q50),col="red") +
    geom_line(aes(y=m),col="blue") +
    facet_wrap(~ variable,scales="free_y")
  return(list(eq=eq,eq_plot=eq_plot,dyn=dyn,dyn_plot=dyn_plot,ctr=ctr))
}