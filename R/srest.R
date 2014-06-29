#' @title Write srest.dat to disk
#' 
#' @description The function is used to write srest.dat to the disk. That file
#' is an input to the \code{srest}-model coded in ADMB by Jose Olivera. 
#' 
#' @author Einar Hjorleifsson
#' 
#' @export
#' 
#' @param name_stock Name of the stock. May later down the line be used in 
#' table and plot outputs
#' @param y1 First year of ssb and recruitment data
#' @param y2 Last year of ssb and recruitment data
#' @param aR Recruitment age
#' @param aP Plus group age
#' @param opt_sr_model Recruitment model number (1: Ricker, 2: Beverton-Holt,
#' 3: Segmented regression).
#' @param opt_pen 0=no SR constraints, 1=apply SR constrain (IS THIS ONLY FOR
#' THE SEGMENTED REGRESSION??).
#' @param r A vector containing recruitment
#' @param ssb A vector containing spawning stock biomass
#' @param year A vector containg years (just used in comments)
#' 
#' @return A file named \code{srest.dat} in the current working directory

mac_srest_cat <- function(name_stock, y1, y2, aR, aP,
                                opt_sr_model, opt_pen,
                                r, ssb, year) {
  tmpfile <- file('srest.dat',open='w')
  cat('# Header: ',name_stock,'\n',file=tmpfile,append=TRUE)
  #cat(name_stock, '# stkname: Name of the stock\n',file=tmpfile,append=TRUE)
  #cat(filename_age,'# filname: Name of the option file (2nd file\n',file=tmpfile,append=TRUE)
  cat(y1,    ' # ybeg:   First year\n',file=tmpfile,append=TRUE)
  cat(y2,    ' # yend:   Last year\n',file=tmpfile,append=TRUE)
  cat(aR,    ' # r:      Recruitment age\n',file=tmpfile,append=TRUE)
  cat(aP,    ' # A:      Plus group age\n',file=tmpfile,append=TRUE)
  cat(opt_sr_model,  ' # Ropt:   S-R function type\n',file=tmpfile,append=TRUE)
  #cat(opt_sim,' # simopt: 0=no simulation, 1=simulation                                (ifelse(nits==0,0,1)) \n',file=tmpfile,append=TRUE)
  #cat(opt_age,' # senopt: 0=error only in recr, 1=error in recr & steady-state vectors (ifelse(varybiodata,1,0))\n',file=tmpfile,append=TRUE)
  cat(opt_pen,' # penopt: 0=no SR constraints, 1=apply SR constraints\n',file=tmpfile,append=TRUE)
  cat('# r ssb\n', file=tmpfile, append=TRUE)
  cat(paste(r,ssb,'#',year), file = tmpfile,append = TRUE,sep="\n")
  close(tmpfile)
}
