#' @title Tigher partial area under the ROC curve
#' @description It standarizes the partial area under the ROC curve by the tigher index
#' @param dataset Dataframe of the complete information of the samples
#' @param low.value inferior limit
#' @param up.value inferior limit
#' @param plot ROC plot
#' @param low.value lower false positive rate value that the function will use to calculate the pAUC
#' @param up.value upper false positive rate value that the function will use to calculate the pAUC
#' @param selection  vector that will only be used if the parameter "dataset" is a RangedSummarizedExperiment object.
#' This parameter is used to select the variables that will be analysed
#' @return RangedSummarizedExperiment object with the pAUC and the TpAUC scores,and the TPR and FPR values for each ROC curve generated
#' @export TpAUC.roc
#' @examples
#'library(fission)
#'data("fission")
#'resultsT <- TpAUC(fission, low.value = 0, up.value = 0.25, plot = TRUE, selection = c("SPNCRNA.1080","SPAC186.08c"))




TpAUC <- function(dataset,  low.value = NULL, up.value = NULL, plot = FALSE, selection = NULL ) {
  St_pAUC <-NULL; pAUC <- NULL; sensitivity <- NULL; FPR <- NULL;
  fpr.proc<-NULL; sen.proc<-NULL;  up.limit <- NULL; low.limit <- NULL
  Ap.roc<-NULL;   object <- NULL;  par <- NULL; legend <- NULL; abline <- NULL;
    ## Variables and initial values for each sample<-ROC curve

  if (class(dataset)=="RangedSummarizedExperiment") {
    strain <- dataset@colData@listData$strain
    dataset <- as.data.frame(SummarizedExperiment::assay(dataset))
    dataset <- scale(t(as.matrix(dataset[selection,])), center=TRUE, scale = TRUE)
    name.variable <- colnames(dataset)
    dataset <- as.data.frame(cbind(strain,dataset))
  }  else {  dataset <- as.data.frame(dataset)
  name.variable <- colnames(dataset)
  }
  dimension <- dim(dataset)

  if(dimension[2]<2) {stop("database has to have at least 2 colums")}


  for (i in 2:dimension[2]) {
    dataset_temporal <- cbind(dataset[,1],dataset[i])
    sen.roc<- points_curve(dataset_temporal[,1],dataset_temporal[,2])[,2]
    fpr.roc<- points_curve(dataset_temporal[,1],dataset_temporal[,2])[,1]

    if (!is.null(up.value)){up.limit <- up.value}else{up.limit <- 1}
    if (!is.null(low.value)){low.limit <- low.value}else{low.limit <- 0}
    fpr.proc <- portion_ROC(up.limit, low.limit, fpr.roc,sen.roc)[,1]
    sen.proc <- portion_ROC(up.limit, low.limit, fpr.roc,sen.roc)[,2]

    St_pAUC[[i-1]] <- TpA(fpr.proc,sen.proc)
    pAUC[[i-1]] <- as.vector(pA(fpr.proc,sen.proc))
    sensitivity[[i-1]] <-as.vector(sen.proc)
    FPR[[i-1]] <- as.vector(fpr.proc)
    if (isTRUE(plot)) {plot(sen.roc~fpr.roc, type="l", col=i, ylab="TPR", xlab="FPR")
      legend(x= "bottomright",legend = name.variable[1:i-1], fill = 2:i, cex = 0.8)
      abline(a=c(0,1),lwd=1, col="grey")
      abline(v=low.value,col="black")
      abline(v=up.value,col="black")
      par(new=TRUE)

      }
  }
  object <- list(St_pAUC, pAUC, sensitivity, FPR)
  names <- c("St_pAUC","pAUC", "Sensitivity", "FPR")
  se <- createSE(object, names)
  return(se)

}

