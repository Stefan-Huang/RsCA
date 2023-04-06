#' Title
#'
#' @param r
#' @param gene
#'
#' @return
#' @export
#'
#' @examples
Get.mean1.scRNA <- function(r,gene){
  library(dplyr)
  library(data.table)
  X<-10*((2^r$tpm)-1)
  r$genes.dist<-log2(rowMeans(X,na.rm = T)+1)
  c<-as.data.frame(r$genes.dist)
  c[gene,]->a
  return(a)
}
