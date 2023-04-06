#' Title
#'
#' @param r
#' @param gene.sign
#' @param num.rounds
#' @param mat.return.flag
#'
#' @return
#' @export
#'
#' @examples
Get.Mean2.scRNA <- function(r,gene.sign,num.rounds = 1000,mat.return.flag = T){
  library(dplyr)
  library(data.table)
  r$genes.mean<-rowMeans(r$tpm)
  r$zscores<-sweep(r$tpm,1,r$genes.mean,FUN = '-')
  if(any(r$tpm<0)){
    print("Using counts to bin genes!!!")
    r$genes.dist<-rowMeans(r$cd>0)
  }else{
    X<-10*((2^r$tpm)-1)
    r$genes.dist<-log2(rowMeans(X,na.rm = T)+1)
  }
  r$genes.dist.q<-discretize(r$genes.dist,n.cat = 50)
  r$sign.scores<-matrix(data = 0,nrow = ncol(r$tpm),ncol = length(gene.sign))
  colnames(r$sign.scores)<-gene.sign
  r$sign.scores.raw<-r$sign.scores
  c.sign<-is.element(r$genes,gene.sign[[1]])
  rand.scores<-Get.random1.score(r,r$genes.dist.q,c.sign,num.rounds = num.rounds)
  d<-as.data.frame(r$zscores)
  raw.scores<-d[gene.sign,]
  final.scores<-raw.scores-rand.scores
  X<-10*((2^final.scores)-1)
  a<-log2(rowMeans(X,na.rm = T)+1)
  return(a)
}
