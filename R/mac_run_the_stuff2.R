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
