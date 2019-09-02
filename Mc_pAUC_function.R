#' @title Clasification of area under ROC curve following McClish method
#' @description Calculate the area under the ROC curve following McClish methodologic from a dataset and a sample from that dataset.
#' @param bsdataset Dataframe of the complete information of the samples
#' @param bsdatasample vector that contains the samples that are going to be analysed
#' @return
#' @export Mc.pAUC
#' @examples
#' install.packages("golubEsets")
#' data(Golub_Merge)

#' dataset <- data.frame(as.factor(dataset[,1])-1,Golub_Merge@assayData$exprs[,1])
#' datasample <- sample(1:length(dataset[,1]),length(dataset:3)
#' Mc_pAUC <- Mc_pAUC_function(dataset, datasample)



   Mc.pAUC <- function(bsdataset, bsdatasample,  valor.inferior=NULL, valor.superior=NULL) {
    ## Variables and initial values for each sample<-ROC curve
    dataset.boot<-bsdataset[bsdatasample, ]
    setsample<-dataset.boot[,1]
    datasample<-dataset.boot[, 2]

    fpr.roc<-points_curve(setsample, datasample)[, 1] #fpr es especificidad VP/VP+FN
    sen.roc<-points_curve(setsample, datasample)[, 2] #sen es sensibilidad VN/VN+FR
    ### ORDER fpr & sen from 0 to 1
    if (is.unsorted(sen.roc)) {
      sen.roc<-rev(sen.roc)
      fpr.roc<-rev(fpr.roc)}
    ## Variables and initial values for the partial area of each ROC curve
    fpr.proc<-NULL; sen.proc<-NULL; proper.proc<-NULL; plr.proc<-NULL; plr.proc.bounded<-NULL
    i.low<-NULL; i.up<-NULL; j.low<-NULL; j.up<-NULL; limite.superior <- NULL; limite.inferior <- NULL
    Ap.roc<-NULL; Mc.pAUC.min.roc<-NULL; TpAUC.roc<-NULL; TpAUC.min.roc<-NULL; TpAUC.max.roc<-NULL
    ### PARTIAL ROC curve (fpr.proc; sen.proc) on [lower.fp <= e <= upper.fp]

    if (!is.null(valor.superior)){limite.superior <- valor.superior}else{limite.superior <- fpr.roc[length(fpr.roc)]}
    if (!is.null(valor.inferior)){limite.inferior <- valor.inferior}else{limite.inferior <- fpr.roc[1]}

    i.low<-min(which(fpr.roc >= limite.inferior))  ## me dice que no existe lower.fp, puesto que es el fp más bajo lo cambio
    j.low<-max(i.low-1, 1)
    i.up<-max(which(fpr.roc <= limite.superior)) ##aqui iría upper.fp, pero es su lugar pongo fpr.roc[length(fpr.roc)]
    j.up<-min(1+i.up, length(fpr.roc))
    fpr.proc<-fpr.roc[i.low:i.up]
    sen.proc<-sen.roc[i.low:i.up]
    if (fpr.roc[i.low] > limite.inferior) {
      fpr.proc<-append(fpr.proc, limite.inferior, 0)
      sen.proc<-append(sen.proc, sen.roc[j.low]+(sen.roc[i.low]-sen.roc[j.low])*(fpr.proc[1]-fpr.roc[j.low])/(fpr.roc[i.low]-fpr.roc[j.low]), 0)}
    if (fpr.roc[i.up] < limite.superior) {
      fpr.proc<-append(fpr.proc, limite.superior, length(fpr.proc))
      sen.proc<-append(sen.proc, sen.roc[j.up]-(sen.roc[j.up]-sen.roc[i.up])*(fpr.roc[j.up]-fpr.proc[length(fpr.proc)])/(fpr.roc[j.up]-fpr.roc[i.up]), length(sen.proc))}
    #fpr.proc es la especificidad
    #sen.proc es la sensibilidad
    Ap.roc<-sum(diff(fpr.proc)*apply(cbind(sen.proc[-1], sen.proc[-length(sen.proc)]), 1, mean))

  ### TYPE OF PORTION ROC curve on [lower.fp <= e <= upper.fp]
  if (all(sen.proc >= fpr.proc)) {proper.proc<-TRUE} else {proper.proc<-FALSE}
  #warning("Improper ROC curve on partial region: min.pauc<-(e2^2-e1^2)/2>pAUC")
  ### BOUNDS FOR PARTIAL AREA INDEXES: MIN & MAX boundaries for Mc.pAUC index
  Mc.pAUC.max.roc<-sum(diff(fpr.proc))
  if (Mc.pAUC.max.roc < Ap.roc) {warning("Impossible: Nonincreasing partial ROC curve")}
  Mc.pAUC.min.roc<-sum(diff(fpr.proc^2))/2
  if (proper.proc) {
    if (Mc.pAUC.max.roc == Ap.roc) {Mc.pAUC.roc<-1} else {
      if (Mc.pAUC.max.roc == Mc.pAUC.min.roc) {Mc.pAUC.roc<-1
      warning("Mc.pAUC index: Constant partial ROC curve")} else {
        Mc.pAUC.roc<-(1+((Ap.roc-Mc.pAUC.min.roc)/(Mc.pAUC.max.roc-Mc.pAUC.min.roc)))/2}
    }} else {Mc.pAUC.roc<-NA}
  #warning("Improper partial ROC curve: McClish's pAUC index is not well defined")
  objet <- list(Area.Parcial.Estandarizada =Mc.pAUC.roc ,Area.Parcial.Sin.Estandarizar=Ap.roc, Sensitivities =sen.proc, Specifities.Complementary = fpr.proc,Mc.pAUC.min.roc=Mc.pAUC.min.roc,Mc.pAUC.max.roc=Mc.pAUC.max.roc)
  return(objet)

}
