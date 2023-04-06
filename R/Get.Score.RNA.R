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
Get.Score.RNA <- function(r,gene.sign,num.rounds = 1000,mat.return.flag = T){
  library(dplyr)
  library(data.table)
  flag1=F
  r$genes.mean<-rowMeans(r$tpm)
  r$zscores<-sweep(r$tpm,1,r$genes.mean,FUN = '-')
  if(any(r$tpm<0)){
    print("Using counts to bin genes!!!")
    r$genes.dist<-rowMeans(r$cd>0)
  }else{
    r$genes.dist<-rowMeans(r$tpm)
  }
  r$genes.dist.q<-discretize(r$genes.dist,n.cat = 50)
  r$sign.scores<-matrix(data = 0,nrow = ncol(r$tpm),ncol = length(gene.sign))
  sign.names<-names(gene.sign)
  colnames(r$sign.scores)<-sign.names
  sig.names<-names(gene.sign)
  r$sign.scores.raw<-r$sign.scores
  for (i in 1:length(gene.sign)){
    b.sign<-is.element(r$genes,gene.sign[[i]])
    if(!flag1){
      flag1=T
      c.sign <- b.sign
    }
    else{
      c.sign[which(b.sign==T)] <- T
    }
  }
  for (i in 1:length(gene.sign)){
    if (sum(c.sign) > 1){
      rand.scores<-Get.random.score(r,r$genes.dist.q,c.sign,num.rounds = num.rounds)
      raw.scores<-colMeans(r$zscores[c.sign,])
      final.scores<-raw.scores-rand.scores
      r$sign.scores[,i]<-final.scores
      r$sign.scores.raw[,i]<-raw.scores
    }
  }
  if(mat.return.flag){
    sign.scores<-r$sign.scores
    sign.scores<-cbind(r$cells,rowMeans(sign.scores))
    return(sign.scores)
  }else{
    return(r)
  }
}
