#' Title
#'
#' @param DeCell_result
#' @param RNA_EXP
#'
#' @return
#' @export
#'
#' @examples
Get.ptScore.RNA<-function(DeCell_result,RNA_EXP){
  library(utils)
  library(estimate)
  library(dplyr)
  write.table(RNA_EXP, "RNA_EXP.txt",sep = "\t",
              row.names = T,col.names = NA,quote = F)
  filterCommonGenes(input.f ="RNA_EXP.txt",   #输入文件名
                    output.f = "RNA_EXP.gct",   #输出文件名
                    id = "GeneSymbol")   #行名为gene symbol
  estimateScore("RNA_EXP.gct",
                "Estimate_score.txt",
                platform="affymetrix")   #默认平台
  est <- read.table("Estimate_score.txt",
                    sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  est <- est[,-1]   #移除第一列
  colnames(est) <- est[1,]   #设置列名
  est <- as.data.frame(t(est[-1,]))
  rownames(est) <- colnames(RNA_EXP)
  sub <- subset(est,TumorPurity>=0.60)
  common<-row.names(sub)
  RNA_EXP<-as.matrix(RNA_EXP)
  RNA_EXP<-RNA_EXP[,common]
  n<-ncol(RNA_EXP)
  common<-intersect(row.names(DeCell_result),row.names(RNA_EXP))
  Ref<-DeCell_result[common,]
  Test<-RNA_EXP[common,]
  Sig_Score<-data_frame(1)
  for(i in 1:n){
    Sig_Score[i,1]<-cor(Test[,i],Ref[,2],method="pearson")
  }
  Sig_Score<-data.frame(Sig_Score)
  row.names(Sig_Score)<-colnames(RNA_EXP)
  return(Sig_Score)
}
