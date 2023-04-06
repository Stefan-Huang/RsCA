#' Title
#'
#' @param scRNA_EXP
#'
#' @return
#' @export
#'
#' @examples
Get.data.tpm<-function(scRNA_EXP){
  library(dplyr)
  library(data.table)
  a<-as.matrix(scRNA_EXP)
  b<-row.names.data.frame(scRNA_EXP)
  c<-colnames(a)
  r<-list(a,a,b,c)
  names(r)<-c("tpm","cd","genes","cells")
  return(r)
}
