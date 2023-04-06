#' Title
#'
#' @param DeCell_result
#' @param RNA_EXP
#'
#' @return
#' @export
#'
#' @examples
Get.pScore.RNA<-function(DeCell_result,RNA_EXP){
  library(dplyr)
  RNA_EXP<-as.matrix(RNA_EXP)
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
