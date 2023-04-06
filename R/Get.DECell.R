#' Title
#'
#' @param Sign_scores
#' @param i
#' @param j
#'
#' @return
#' @export
#'
#' @examples
Get.DECell<-function(Sign_scores,i,j){
  library(dplyr)
  library(data.table)
  library(Seurat)
  library(SeuratObject)
  library(limma)
  as.data.frame(Sign_scores)->Sign_scores
  as.numeric(Sign_scores$V2)->Sign_scores$V2
  Sorted <- Sign_scores[order(Sign_scores$V2),]
  n <- nrow(Sorted)
  ni <- round(0.01*i * nrow(Sorted))
  nj <- round(0.01*j * nrow(Sorted))
  df_i <- Sorted[1:ni, 1]
  df_j <- Sorted[nj:n, 1]
  exp<-matrix(r$tpm,nrow=nrow(r$tpm),ncol=ncol(r$tpm),dimnames =list(r$genes,r$cells))
  exp<-cbind(exp[,df_j],exp[,df_i])
  group<-c(rep('a',times=(n-nj+1)),rep('b',times=ni))
  seuratobject <- CreateSeuratObject(counts = exp,project = "Leukocyte_sample", min.cells = 3, min.features = 200)
  Idents(seuratobject)<-group
  seuratobject[["RNA"]]@data=as(as.matrix(log(exp + 1)), "dgCMatrix")
  results<- FindMarkers(seuratobject, ident.1="a", ident.2="b", slot="data", logfc.threshold=0.25, min.pct=0.1, test.use = "wilcox")
  result<-data.frame(row.names(results),results$avg_log2FC)
  row.names(result)<-row.names(results)
  return(result)
}
