#' Title
#'
#' @param Sign_scores
#'
#' @return
#' @export
#'
#' @examples
Get.prop<-function(Sign_scores){
  library(dplyr)
  library(data.table)
  Sign_scores<-as.data.frame(Sign_scores)
  length(Sign_scores$V2[Sign_scores$V2>0])/length(Sign_scores$V2)
}
