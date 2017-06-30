#' @title HCR - Read control file
#'
#' @description Reads a simple text file and returns a list
#' 
#' @export
#' 
#' @param file filename

hcr_read_ctr <- function (file) 
{
  x <- read.table(file)
  ctr <- list()
  for (i in 1:nrow(x)) ctr[[i]] <- x[i,]
  names(ctr) <- c("a1","a2","y1","y2","iter","f1","f2","nR","tac_y1","tac_y2",
                  "y1Bias","r_cv","r_rho","r_model","r_mean","ssb_break",
                  "a_cv","a_rho","a_error","a_bias","w_cv","w_rho","w_error",
                  "w_refB","h_alpha","h_beta","b_trigger","delay","h_number",
                  "i_number","b2","b2")
  return(ctr)
}
