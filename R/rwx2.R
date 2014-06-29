# Amended fishvise functions

#' @title Reads the stock in numbers and fishing mortality
#' 
#' @description Reads fishing mortality and stock in numbers by year and age
#' from \code{sam}.
#' 
#' @export
#' 
#' @return Returns a \code{data.frame} containing the following items:
#' 
#' @param x Object from read.fit
#' @param Scale A numerical value. If 1 (default) returns the stock in numbers
#' as in the input file.
#' 
read.rbya <- function(x,Scale=1) {
  
  minAge <- min(x$res[,3])  
  maxAge <- max(x$res[,3][x$res[,3]<98.5])  
  noN <- maxAge - minAge+1
  noFleet <- max(x$res[,2])
  
  N <- exp(x$stateEst[,1:noN])/Scale
  colnames(N) <- c(minAge:maxAge)
  rownames(N) <- x$years
  N <- melt(N,factorsAsStrings = FALSE)
  
  mort <- exp(x$stateEst[,-c(1:noN)])[,x$keys$keyLogFsta[1,]]
  colnames(mort) <- c(minAge:maxAge)
  rownames(mort) <- x$years
  mort <- melt(mort,factorsAsStrings = FALSE)
  
  res <- cbind(N,mort[,3])
  names(res) <- c("year","age","n","f")
  
  
  return(res)
}

#' @title Reads the stock summary values
#' 
#' @description Reads the stock summary numbers by year from \code{sam}.
#' 
#' 
#' @export
#'
#' @return A \code{data.frame} containing the following:
#'  
#' @param x Object from function \code{read.fit}.
#' @param range Boolean Not used
#' @param Scale A value 
read.rby <- function(x,range=FALSE,Scale=1) {  
  rby <- cbind(x$fbar,x$ssb,x$tsb)
  rby[,5:12] <- rby[,5:12]/Scale
  colnames(rby) <- paste(c(rep("f",4),rep("ssb",4),rep("tsb",4)),colnames(rby),sep="")
  rby <- data.frame(rby)
  # add recruitment
  rby$r <- x$R[,1]/Scale
  
  # add yield 
  rby$yield <- c(exp(x$logCatch[,1]),NA) / Scale
  
  rby$year <- x$years
  rownames(rby) <- NULL
  return(rby)
}

#' @title read.lowestoft
#' 
#' @description Read lowestoft files
#' 
#' @export
#' 
#' @param filename Filename
#' @param val.name Character Not used
#' @param Format Character, "List" or "Wide"
#' 
read.lowestoft <- function(filename, val.name,Format="List")
{
  y <- scan(filename, skip = 2, nlines = 1, quiet = TRUE)
  a <- scan(filename, skip = 3, nlines = 1, quiet = TRUE)
  tab <- read.delim(filename, header = FALSE, sep = "", skip = 5)
  names(tab) <- c(a[1]:a[2])
  rownames(tab) <- c(y[1]:y[2])
  if(Format == "List") return(list(y = y, a = a, tab = tab))
  if(Format == "Wide") return(tab)
  tab$year <- as.integer(rownames(tab))
  tab <- melt(tab,id.vars="year",factorsAsStrings = FALSE)
  names(tab) <- c("year","age",val.name)
  tab$age <- as.integer(as.character(tab$age))
  return(tab)
}

#' @title read.ibya
#' 
#' @description read input by year and age
#' 
#' @export
#' 
#' @param path Character containing the path to the directory
#' @param Scale A value 
#' 
read.ibya <- function(path,Scale=1) {
  oc <-  read.lowestoft(paste(path,"cn.dat",sep="/"),val.name="oC",Format = "Long")
  oc$oC <- oc$oC/Scale
  cw <-  read.lowestoft(paste(path,"cw.dat",sep="/"),val.name="cW",Format = "Long")
  sw <-  read.lowestoft(paste(path,"sw.dat",sep="/"),val.name="sW",Format = "Long")
  mat <- read.lowestoft(paste(path,"mo.dat",sep="/"),val.name="mat",Format = "Long")
  nat <- read.lowestoft(paste(path,"nm.dat",sep="/"),val.name="m",Format = "Long")
  pf <-  read.lowestoft(paste(path,"pf.dat",sep="/"),val.name="pf",Format = "Long")
  pm <-  read.lowestoft(paste(path,"pm.dat",sep="/"),val.name="pm",Format = "Long")
  
  res <- join(oc,sw,by=c("year","age"))
  res <- join(res,cw,by=c("year","age"))
  res <- join(res,mat,by=c("year","age"))
  res <- join(res,nat,by=c("year","age"))
  res <- join(res,pf,by=c("year","age"))
  res <- join(res,pm,by=c("year","age"))
  return(res)
}


