#' Title
#'
#' @param v
#' @param n.cat
#'
#' @return
#' @export
#'
#' @examples
discretize<-function(v,n.cat){
  library(dplyr)
  library(data.table)
  q1<-quantile(v,seq(from = (1/n.cat),to = 1,by = (1/n.cat)))
  u<-matrix(nrow = length(v))
  for(i in 2:n.cat){
    u[(v>=q1[i-1])&(v<q1[i])]<-i
  }
  return(u)
}
