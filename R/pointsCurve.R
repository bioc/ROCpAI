#' @title Points of the ROC curve
#' @description It calculates the coordinates (fpr, sen) of the ROC curve.
#' This function sorts the scores of a model test and generates
#' the points which will be used to plot its the ROC curve
#' @param x It is the vector of the status (gold standar)
#' @param y It is the vector with the values of a predictor
#' variable or clasificator
#' @return return a matrix with the points of 1-specificity and
#' sensibility that will be used to generate a ROC curve
#' @export pointsCurve
#' @examples
#'library(fission)
#'data("fission")
#'strain <- fission@colData@listData$strain
#'pointsCurve<- pointsCurve(strain, t(assay(fission))[,"SPNCRNA.1080"])

##1st col<-x<-gold standard, 2nd col<-y<-sample
pointsCurve<-function(x, y){
  stopifnot(is.numeric(y)||is.integer(y))
  xsample <- NULL; ysample <- NULL
  xsample <- cbind(x[which(is.na(x)==FALSE & is.na(y)==FALSE)])
  ysample <- y[which(is.na(x)==FALSE & is.na(y)==FALSE)]
  points<-NULL; pre.point<-NULL; fpr.point<-NULL; sen.point<-NULL; xy<-NULL
  points<-sort(ysample)
  points<-append(points[-length(points)]+diff(points)/2, min(points)-1, 0)
  points<-append(points, max(ysample)+1, length(points))
  for (point in points) {
    pre.point<-(ysample>point)*1
    fpr.point[which(points == point)]<-sum((pre.point == 1)*(xsample == 1))/sum(xsample == 1)
    sen.point[which(points == point)]<-sum((pre.point == 1)*(xsample == 2))/sum(xsample == 2)

  }
  if (is.unsorted(sen.point)) {
    sen.point<-rev(sen.point)
    fpr.point<-rev(fpr.point)}
  xy<-cbind(fpr.point, sen.point)
  return(xy)
}



