#' @title points.curve function to calculate the coordinates (fpr,sen)
#' @description Function that calculates a serie of coordinates that will be used to do the curve ROC.
#' This function sorts a vector resulted of a model test and generates several points that will be used to divide the ROC curve and analyse it
#' @param xsample sample etiquetes
#' @param ysample vector with the values of sensibility of a model
#' @return return a matrix with a serie of coordenates that will be used to generate a ROC curve
#' @export points.curve
#' @examples
#' install.packages("golubEsets")
#' data(Golub_Merge)
#' dataset <-data.frame(as.factor(dataset[,1])-1,Golub_Merge@assayData$exprs[,1])
#' datasample <- sample(1:length(dataset[,1]),length(dataset:3)
#' pointsCurve<- points_curve(dataset[,1], dataset[,2])

########## PRELIMINARIES
### Define data, variables, initial values, functions...
### Limites (lw; up) de fpr de la regi?n parcial de inter?s
lw.extreme<-NULL; up.extreme<-NULL; n.extreme<-NULL; i.extreme<-NULL
datafiles<-NULL; i.datafiles<-NULL; n.datafiles<-NULL
datasets<-NULL; dataset<-NULL
##########
##1st col<-xsample<-gold standard, 2nd col<-ysample<-sample
points_curve<-function(xsample, ysample){
  points<-NULL; prepoint<-NULL; fpr.point<-NULL; sen.point<-NULL; xy<-NULL
  points<-sort(ysample)
  points<-append(points[-length(points)]+diff(points)/2, min(points)-1, 0)
  points<-append(points, max(ysample)+1, length(points))
  for (point in points) {
    pre.point<-(ysample>point)*1
    fpr.point[which(points == point)]<-sum((pre.point == 1)*(xsample == 0))/sum(xsample == 0)
    sen.point[which(points == point)]<-sum((pre.point == 1)*(xsample == 1))/sum(xsample == 1)
  }
  if (is.unsorted(sen.point)) {
    sen.point<-rev(sen.point)
    fpr.point<-rev(fpr.point)}
  xy<-cbind(fpr.point, sen.point)
  return(xy)
}
