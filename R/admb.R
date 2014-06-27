#' @title Compile  ADMB code
#' 
#' @description The \code{plotMSY} approach uses a ADModel Builder code. ...
#' 
#' @export
#' 
#' @param tpl_name Name of the admb code
#' @param windose Whether in windows (default FALSE) or not
#' @param compile Compile binaries from source. TRUE requires ADMB to be installed, FALSE only works on Windows

compile_admb2 <- function(tpl_name="srest.tpl",windose=FALSE,compile=TRUE) {
  if(!(compile | windose)) stop("Must compile in non-Windows operating systems")
  if(compile)
  {
    file.copy(paste(path.package("mac"),"/extdata/",tpl_name,sep=''),tpl_name)
    tpl <- str_replace(tpl_name,".tpl","")
    compile_admb(tpl)
    clean_admb(tpl)
  } else {
    stop(".exe file not yet available")
    #file.copy(paste(path.package("mac"),"/bin/",tpl,".exe",sep=""),paste(tpl,".exe","")
  }
}