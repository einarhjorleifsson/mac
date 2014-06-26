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

#' @title mac_run_the_stuff
#' 
#' @description A wrapper function for running the mackerel simulation
#' 
#' @export
#' 
#' @param dat A data.frame containing the mean values
#' @param ctr The control file

mac_run_the_stuff <- function(dat,ctr) {
  monsterStructure <- hcr_set_dimensions(ctr=ctr)
  hcr_set_starting_conditions(dat=dat, d=monsterStructure, ctr=ctr)
  # Put in starting stock numbers from sam
  X$N[,1,,] <<- ctr$start_N[,1:ctr$iter]/1e6
  HRATE <- ctr$HRATE
  # the first year
  y <- 1
  for (h in 1:length(HRATE)) {
    Fmult <- hcr_TAC_to_Fmult(y, h)
    hcr_operating_model(y, h, ctr, Fmult, nR = 1)
  }
  
  # the other years 
  for (h in 1:length(HRATE)) {
    for(y in 2:(ctr$y2 - ctr$y1)) {
      # Outcommented: Implementation error
      # switch(ctr$i_number,
      #       dummy <- hcr_implementation_model_1(),
      #       stop("implementation model 2 is available for developement"),
      #       stop(" and whatever more ..."))
      
      # Get the realized Fmultiplier given the TAC decision made in the year
      # before
      Fmult <- hcr_TAC_to_Fmult(y, h)
      # Apply the Fmult to the real population numbers
      hcr_operating_model(y, h, ctr, Fmult, nR=1)
      
      # Next years TAC
      #
      # Below formerly had:
      # hat <- hcr_observation_error(y,h,Fmult,ctr)
      # The above is a bug
      hat <- hcr_observation_error(y,h,Fmult,ctr)
      # Do predictions for next year and calculate the TAC
      switch(ctr$h_number,
             stop("Not implemented yet"),
             hcr_management_fmort(y,h,hat$hrate,hat$ssb,ctr),
             hcr_management_bio(y,h,hat$bio,hat$ssb,ctr))
    }
  }
  
  # Complete the last year
  y <- ctr$y2 - ctr$y1 + 1
  for (h in 1:length(HRATE)) {
    Fmult <- hcr_TAC_to_Fmult(y, h)
    hcr_operating_model(y, h, ctr, Fmult, nR = 1)
  }
}

#' @title mac_run_the_stuff2
#' 
#' @description A wrapper function for running the mackerel simulation. Here
#' use the function \code{hcr_management_fmort_pred} function.
#' 
#' @export
#' 
#' @param dat A data.frame containing the mean values
#' @param ctr The control file

mac_run_the_stuff2 <- function(dat,ctr) {
  monsterStructure <- hcr_set_dimensions(ctr=ctr)
  hcr_set_starting_conditions(dat=dat, d=monsterStructure, ctr=ctr)
  # Put in starting stock numbers from sam
  X$N[,1,,] <<- ctr$start_N[,1:ctr$iter]/1e6
  HRATE <- ctr$HRATE
  # the first year
  y <- 1
  for (h in 1:length(HRATE)) {
    Fmult <- hcr_TAC_to_Fmult(y, h)
    hcr_operating_model(y, h, ctr, Fmult, nR = 1)
  }
  
  # the other years 
  for (h in 1:length(HRATE)) {
    for(y in 2:(ctr$y2 - ctr$y1)) {
      # Outcommented: Implementation error
      # switch(ctr$i_number,
      #       dummy <- hcr_implementation_model_1(),
      #       stop("implementation model 2 is available for developement"),
      #       stop(" and whatever more ..."))
      
      # Get the realized Fmultiplier given the TAC decision made in the year
      # before
      Fmult <- hcr_TAC_to_Fmult(y, h)
      # Apply the Fmult to the real population numbers
      hcr_operating_model(y, h, ctr, Fmult, nR=1)
      
     # Do predictions for next year and calculate the TAC
      switch(ctr$h_number,
             stop("Not implemented yet"),
             hcr_management_fmort_pred(y,h,ctr),
             hcr_management_bio_pred(y,h,ctr))
    }
  }
  
  # Complete the last year
  y <- ctr$y2 - ctr$y1 + 1
  for (h in 1:length(HRATE)) {
    Fmult <- hcr_TAC_to_Fmult(y, h)
    hcr_operating_model(y, h, ctr, Fmult, nR = 1)
  }
}

