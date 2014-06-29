#' @title Compile admb code
#' 
#' @description The function compiles the \code{srest.tpl} admb-code that is bundle in the
#' mac package in the current working directory.
#' 
#' The function is a wrapper around the \code{compile_admb} function from
#' the \code{R2admb}-package.
#' 
#' The function can be used to compile \code{.tpl} codes that come with other R-packages.
#' 
#' @export
#' 
#' @return Nothing is returned to the R-console but leaves an executable file
#' in the current working directory.
#' 
#' @param tpl_name Name of the admb code
#' @param windose Whether in windows (default FALSE) or not
#' @param compile Compile binaries from source. TRUE requires ADMB to be installed, FALSE only works on Windows

compile_admb2 <- function(tpl_name="srest.tpl",windose=FALSE,compile=TRUE) {
  #if(!(compile | windose)) stop("Must compile in non-Windows operating systems")
  #if(compile)
  #{
    file.copy(paste(path.package("mac"),"/extdata/",tpl_name,sep=''),tpl_name)
    tpl <- str_replace(tpl_name,".tpl","")
    compile_admb(tpl)
    clean_admb(tpl)
  #} else {
  #  stop(".exe file not yet available")
    #file.copy(paste(path.package("mac"),"/bin/",tpl,".exe",sep=""),paste(tpl,".exe","")
  #}
}