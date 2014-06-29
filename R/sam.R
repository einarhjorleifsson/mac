# roxygen2-ized sam-function

#' @title get_sam_run
#' 
#' @description Fetches the results from the run directory on stockassessment.org
#' and stored it locally.
#' 
#' @note Currently the directory structure is also copied, see datadisk/stock...
#' 
#' @export
#' 
#' @param dir Name of the directory on stockassessment.org

get_sam_run <- function(dir) {
  cmd_prefix <- "wget -r - l1 -nH --no-parent --user='guest' --password=guest https://www.stockassessment.org/datadisk/stockassessment/userdirs/user3/"
  cmd <- paste(cmd_prefix,dir,sep="")
  system(cmd)
}


#' @title read.conf
#' 
#' @description Reads sam configuration file. Used by read.fit function
#' (see below). This function originally in scr/common.R
#' 
#' @export
#' 
#' @param clonefile XXX
read.conf <- function(clonefile){
  lin <- readLines(clonefile)
  idxNam <- grep("^[[:alpha:]]",lin)
  doone <- function(i){
    idxNam <- c(idxNam,length(lin)+1)
    x <- read.table(textConnection(lin[(idxNam[i]+1):(idxNam[i+1]-1)]))
    names(x) <- NULL
    as.matrix(x)
  }
  ret<-lapply(1:length(idxNam),doone)
  names(ret)<-sub(' =','',lin[idxNam])
  ret
}


#' @title read.fit
#' 
#' @description Seems like the big mother.
#' 
#' This function originally in scr/common.R
#' 
#' @export
#' 
#' @param file Name of file
#' @param reduced Boolean 
read.fit <- function(file, reduced=FALSE){
  # Function to read a basic fit
  ret<-list()
  parfile<-as.numeric(scan(paste(file,'.par', sep=''), 
                           what='', n=16, quiet=TRUE)[c(6,11,16)])
  ret$nopar<-as.integer(parfile[1])
  ret$nlogl<-parfile[2]
  ret$maxgrad<-parfile[3]
  rep<-scan(paste(file,'.rep', sep=''), quiet=TRUE)
  ret$res<-read.table(paste(file,'.res', sep=''),header=FALSE)
  ret$stateDim<-rep[1]
  ret$years<-rep[-1]
  file<-paste(file,'.cor', sep='')
  lin<-readLines(file)
  ret$npar<-length(lin)-2
  ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2])
  sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!=''])
  ret$names<-unlist(lapply(sublin,function(x)x[2]))
  ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3])))
  ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4])))
  
  ret$cor<-matrix(NA, ret$npar, ret$npar)
  for(i in 1:ret$npar){
    ret$cor[1:i,i]<-as.numeric(unlist(lapply(sublin[i],
                                             function(x)x[5:(4+i)])))
    ret$cor[i,1:i]<-ret$cor[1:i,i]
  }
  ret$cov<-ret$cor*(ret$std%o%ret$std)
  
  mslh <- function(name){
    idx<-which(ret$names==name)
    x<-cbind(ret$est[idx], ret$std[idx], ret$est[idx]-2*ret$std[idx], 
             ret$est[idx]+2*ret$std[idx])
    colnames(x)<-c('est', 'std', 'low', 'hig')
    return(x)
  }
  
  ret$ssb<-mslh('ssb')
  ret$fbar<-mslh('fbar')
  ret$tsb<-mslh('tsb')
  ret$logssb<-mslh('logssb')
  ret$logfbar<-mslh('logfbar')
  ret$logtsb<-mslh('logtsb')
  ret$logscale<-mslh('logScale')
  ret$logFpar<-mslh('logFpar')
  ret$logCatch<-mslh('logCatch')
  x<-mslh('U')
  ret$stateEst<-matrix(x[,1],ncol=ret$stateDim, byrow=TRUE)
  ret$stateStd<-matrix(x[,2],ncol=ret$stateDim, byrow=TRUE)
  ret$stateLow<-matrix(x[,3],ncol=ret$stateDim, byrow=TRUE)
  ret$stateHig<-matrix(x[,4],ncol=ret$stateDim, byrow=TRUE)
  ret$R<-cbind(exp(ret$stateEst[,1]), NA, exp(ret$stateLow[,1]), 
               exp(ret$stateHig[,1]))
  if(reduced){
    ret <- ret[which(!names(ret)%in%c('cov','cor'))]
  }
  
  file <- sub('[[:alpha:]]+\\.cor$','confclone.log',file)
  if(file.exists(file)){
    ret$keys<-read.conf(file)
  }
  return(ret)
}



