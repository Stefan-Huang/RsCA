#' Title
#'
#' @param scRNA_EXP
#' @param length
#'
#' @return
#' @export
#'
#' @examples
Get.tpm<-function(scRNA_EXP,length){
  library(dplyr)
  library(data.table)
  a<-as.matrix(scRNA_EXP)
  scRNA_EXP$Gene<-row.names.data.frame(scRNA_EXP)
  merge<-merge(scRNA_EXP,length,by="Gene")
  merge <- na.omit(merge)
  rownames(merge)<-merge[,1]
  b<-merge[,1]
  merge<-merge[,-1]
  n<-ncol(merge)
  c<-colnames(merge)
  #TPM
  kb <- merge$Length / 1000
  countdata <- merge[,1:n]
  rpk <- countdata / kb
  tpm <- t(t(rpk)/colSums(rpk) * 1000000)
  tpm<-tpm[,-n]
  return(tpm)
}
