
#to clasificate the portion of the ROC curve
#portion ROC over a fixed interval
portion_ROC <- function(up.limit, low.limit, fpr.roc,sen.roc){

  i.low<-min(which(fpr.roc >= low.limit))
  j.low<-max(i.low-1, 1)
  i.up<-max(which(fpr.roc <= up.limit))
  j.up<-min(1+i.up, length(fpr.roc))
  fpr.proc<-fpr.roc[i.low:i.up]
  sen.proc<-sen.roc[i.low:i.up]
  if (fpr.roc[i.low] > low.limit) {
    fpr.proc<-append(fpr.proc, low.limit, 0)
    sen.proc<-append(sen.proc, sen.roc[j.low]+(sen.roc[i.low]-sen.roc[j.low])*(fpr.proc[1]-fpr.roc[j.low])/(fpr.roc[i.low]-fpr.roc[j.low]), 0)}
  if (fpr.roc[i.up] < up.limit) {
    fpr.proc<-append(fpr.proc, up.limit, length(fpr.proc))
    sen.proc<-append(sen.proc, sen.roc[j.up]-(sen.roc[j.up]-sen.roc[i.up])*(fpr.roc[j.up]-fpr.proc[length(fpr.proc)])/(fpr.roc[j.up]-fpr.roc[i.up]), length(sen.proc))}

  return(cbind(fpr.proc, sen.proc))
}

classification_Tp <- function(fpr.proc,sen.proc){

  if (all(sen.proc>=fpr.proc)) {
    proper.proc<-TRUE} else {proper.proc<-FALSE}
  plr.proc<-(sen.proc-sen.proc[1])/(fpr.proc-fpr.proc[1])
  plr.proc<-plr.proc[is.finite(plr.proc)]
  if (all(plr.proc>=plr.proc[length(plr.proc)])) {
    plr.proc.bounded<-TRUE} else {plr.proc.bounded<-FALSE}

  classification <- c(plr.proc.bounded,proper.proc)

  return(classification)
}

#tighter partial area under a portion
TpA <- function(fpr.proc, sen.proc){
  pA.roc <- pA(fpr.proc, sen.proc)
  type_roc <- classification_Tp(fpr.proc, sen.proc)

  min.pAUC<-sum(diff(fpr.proc^2))/2
  max.pAUC<-sum(diff(fpr.proc))
  TpAUC.max.roc<-max.pAUC*max(sen.proc)

  if (min(sen.proc) == max(sen.proc) ) {
    TpAUC.min.roc =0} else {
      if (type_roc[1]) {
        TpAUC.min.roc<-sum(diff(fpr.proc))*mean(c(min(sen.proc), max(sen.proc)))
        } else {

        if (type_roc[2]) {
            TpAUC.min.roc<-max(max.pAUC*min(sen.proc), min.pAUC)
          } else {TpAUC.min.roc<-max.pAUC*min(sen.proc)}
        }
      }
  if (min(fpr.proc) == max(fpr.proc) ) { TpAUC.max.roc = 1}

  if (max(sen.proc)!=0 )
  {TpA.roc<-(1+((pA.roc-TpAUC.min.roc)/(TpAUC.max.roc-TpAUC.min.roc)))/2} else {TpA.roc=0}

  return(TpA.roc)

}
#Calculate the partial area under a portion
pA <- function(fpr.proc,sen.proc){
  aux <- sum(diff(fpr.proc)*
              apply(cbind(sen.proc[-1],
                          sen.proc[-length(sen.proc)]), 1, mean))
  return(aux)
}

#McClish partial area under a portion

MCpA <- function(sen.proc,fpr.proc){
  pA.roc <- pA(fpr.proc, sen.proc)
  type_roc <- classification_Tp(fpr.proc, sen.proc)

  max.pAUC <- sum(diff(fpr.proc))
  min.pAUC <- sum(diff(fpr.proc^2))/2

  if (type_roc[2]) {
    if (min(fpr.proc) == max(fpr.proc) ) { max.pAUC = 1}

            MCpA.roc<-(1+((pA.roc-min.pAUC)/(max.pAUC-min.pAUC)))/2
    } else {
      MCpA.roc<-NA}

  #warning("Improper partial ROC curve: McClish's pAUC index is not well defined")
  return(MCpA.roc)
}

fbootT<- function(dataset,bssample, low.limit, up.limit){
  SpAUC <- NULL
  bsdata <- dataset[bssample,]
  for (i in 2:dim(bsdata)[2]) {
    bsdata_temp <- cbind(bsdata[,1],bsdata[,i])
    sen.roc<- points_curve(bsdata_temp[,1],bsdata_temp[,2])[,2]
    fpr.roc<- points_curve(bsdata_temp[,1],bsdata_temp[,2])[,1]

    fpr.proc <- portion_ROC(up.limit, low.limit, fpr.roc,sen.roc)[,1]

    sen.proc <- portion_ROC(up.limit, low.limit, fpr.roc,sen.roc)[,2]
    SpAUC[i-1] <- TpA(fpr.proc,sen.proc)


  }
  return(SpAUC)
}

fbootM<- function(dataset,bssample, low.limit, up.limit){
  SpAUC <- NULL
  bsdata <- dataset[bssample,]
  for (i in 2:dim(bsdata)[2]) {
    bsdata_temp <- cbind(bsdata[,1],bsdata[,i])
    sen.roc<- points_curve(bsdata_temp[,1],bsdata_temp[,2])[,2]
    fpr.roc<- points_curve(bsdata_temp[,1],bsdata_temp[,2])[,1]

    fpr.proc <- portion_ROC(up.limit, low.limit, fpr.roc,sen.roc)[,1]

    sen.proc <- portion_ROC(up.limit, low.limit, fpr.roc,sen.roc)[,2]
    SpAUC[i-1] <- MCpA(sen.proc,fpr.proc)


  }
  return(SpAUC)
}

createSE <- function(object, names){
  names(object) <- names
  names <- names

  data.matrix <- as.matrix(object)
  #cambiar
  se=SummarizedExperiment::SummarizedExperiment(assays=data.matrix,
                          colData<-data.frame(metrics = (names)))
  names(se@assays@data@listData[[1]]) <- names
  return(se)
}


