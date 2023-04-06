#' Title
#'
#' @param r
#' @param genes.dist.q
#' @param c.sign
#' @param num.rounds
#' @param full.flag
#'
#' @return
#' @export
#'
#' @examples
Get.random.score <- function(r,genes.dist.q,c.sign,num.rounds = 1000,full.flag = F){
  library(dplyr)
  library(data.table)
  sign.q<-as.matrix(table(genes.dist.q[c.sign]))
  q<-rownames(sign.q)
  idx.all<-c()
  B<-matrix(data = F,nrow = length(genes.dist.q),ncol = num.rounds)
  Q<-matrix(data = 0,nrow = length(genes.dist.q),ncol = num.rounds)
  for (i in 1:nrow(sign.q)){
    num.genes<-sign.q[i]
    if(num.genes>0){
      idx<-which(is.element(genes.dist.q,q[i]))
      for (j in 1:num.rounds){
        idxj<-sample(idx,num.genes)
        B[idxj,j]<-T
        Q[i,j]<-sum(B[idxj,j])
      }
    }
  }
  rand.scores<-apply(B,2,function(x) colMeans(r$zscores[x,]))
  if(full.flag){return(rand.scores)}
  rand.scores<-rowMeans(rand.scores)
  return(rand.scores)
}
