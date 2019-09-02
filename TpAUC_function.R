#' @title Clasification of partial ROC curve
#' @description Calculates the partial ROC curve and analyses it so it can be classified
#' @param bsdataset Dataframe of the complete information of the samples
#' @param bsdatasample vector that contains the samples that are going to be analysed
#' @return
#' @export TpAUC.roc
#' @examples
#' install.packages("golubEsets")
#' data(Golub_Merge)
#' dataset<- data.frame(as.factor(dataset[,1])-1,Golub_Merge@assayData$exprs[,1])
#' datasample <- sample(1:length(dataset[,1]),length(dataset:3)
#' TpAUC <- TpAUC.function(dataset, datasample)





TpAUC <- function(bsdataset, bsdatasample,  valor.inferior=NULL, valor.superior=NULL) {
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
  if (all(sen.proc>=fpr.proc)) {proper.proc<-TRUE} else {proper.proc<-FALSE}
  #warning("Improper ROC curve on partial region: min.pauc<-(e2^2-e1^2)/2>pAUC")
  plr.proc<-(sen.proc-sen.proc[1])/(fpr.proc-fpr.proc[1])
  plr.proc<-plr.proc[is.finite(plr.proc)]
  if (all(plr.proc>=plr.proc[length(plr.proc)])) {plr.proc.bounded<-TRUE} else {plr.proc.bounded<-FALSE}
  #warning("PORTION ROC curve over partial plr")
  ### BOUNDS FOR PARTIAL AREA INDEXES: MIN & MAX boundaries for TpAUC index
  Mc.pAUC.min.roc<-sum(diff(fpr.proc^2))/2
  TpAUC.max.roc<-sum(diff(fpr.proc))*max(sen.proc)

  if (TpAUC.max.roc < Ap.roc) {warning("Impossible: Nonincreasing partial ROC curve")
    TpAUC.max.roc<-sum(diff(fpr.proc))}
  if (plr.proc.bounded) {TpAUC.min.roc<-sum(diff(fpr.proc))*mean(c(min(sen.proc), max(sen.proc)))
  } else {if (proper.proc) {
    TpAUC.min.roc<-max(sum(diff(fpr.proc))*min(sen.proc), Mc.pAUC.min.roc)} else {
      TpAUC.min.roc<-sum(diff(fpr.proc))*min(sen.proc)}}

  if (TpAUC.min.roc > Ap.roc) {warning("Impossible: Improper & nonincreasing partial ROC curve")}
  if (TpAUC.max.roc == Ap.roc) {TpAUC.roc<-1} else {
    if (TpAUC.max.roc == TpAUC.min.roc) {#max=min => max=Ap => this case is impossible
      TpAUC.roc<-1
      warning("TpAUC index: Constant partial ROC curve")} else {
        TpAUC.roc<-(1+((Ap.roc-TpAUC.min.roc)/(TpAUC.max.roc-TpAUC.min.roc)))/2}} #esto calcula la diferencia entre la curva sin pasar y cuando cruza
  #correccion de McClish para aucROC impropias (curtan la daigonal.)


  object <- list(Area.Parcial.Estandarizada =TpAUC.roc ,Area.Parcial.Sin.Estandarizar=Ap.roc, Sensitivities =sen.proc, Specifities.Complementary = fpr.proc,TpAUC.max.roc=TpAUC.max.roc,TpAUC.min.roc=TpAUC.min.roc)
  return(object)
}
