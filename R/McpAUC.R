#' @title Clasification of area under ROC curve following McClish method
#' @description Calculate the area under the ROC curve following McClish methodologic from a dataset and a sample from that dataset.
#' @param dataset Dataframe of the complete information of the samples
#' @param low.value lower false positive rate value that the function will use to calculate the pAUC
#' @param up.value upper false positive rate value that the function will use to calculate the pAUC
#' @param selection  vector that will only be used if the parameter "dataset" is a RangedSummarizedExperiment object.
#' This parameter is used to select the variables that will be analysed
#' @param variable in case that dataset is a  SummarizedExperiment, indicate the Gold Standard
#' @param plot ROC plot
#' @return RangedSummarizedExperiment object with the pAUC and the McpAUC scores,and the TPR and FPR values for each ROC curve generated
#' @export Mc.pAUC
#' @examples
#'library(fission)
#'data("fission")
#'resultsMC <- McpAUC(fission, low.value = 0, up.value = 0.25, plot = TRUE, selection = c("SPNCRNA.1080","SPAC186.08c"), varidable="strain")




McpAUC <- function(dataset,  low.value=NULL, up.value=NULL, plot=FALSE, selection = NULL, variable=NULL) {
  St_pAUC <-NULL; pAUC <- NULL; sensitivity <- NULL; FPR <- NULL;
  fpr.proc<-NULL; sen.proc<-NULL;  up.limit <- NULL; low.limit <- NULL;
  Ap.roc<-NULL;   object <- NULL;par(new=FALSE);  par <- NULL; legend <- NULL; abline <- NULL;

  par(new=FALSE)
  name.variable <- colnames(dataset)
  if (class(dataset)=="RangedSummarizedExperiment") {
    strain <- dataset@colData@listData
    strain <- strain[variable][[1]]
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

    ## Variables and initial values for the partial area of each ROC curve
    ### PARTIAL ROC curve (fpr.proc; sen.proc) on [lower.fp <= e <= upper.fp]

    if (!is.null(up.value)){up.limit <- up.value}else{up.limit <- 1}
    if (!is.null(low.value)){low.limit <- low.value}else{low.limit <- 0}
    fpr.proc <- portion_ROC(up.limit, low.limit, fpr.roc,sen.roc)[,1]
    sen.proc <- portion_ROC(up.limit, low.limit, fpr.roc,sen.roc)[,2]

         ### TYPE OF PORTION ROC curve on [lower.fp <= e <= upper.fp]


         #warning("Improper ROC curve on partial region: min.pauc<-(e2^2-e1^2)/2>pAUC")
         #warning("PORTION ROC curve over partial plr")
         ### BOUNDS FOR PARTIAL AREA INDEXES: MIN & MAX boundaries for TpAUC index
    St_pAUC[[i-1]] <- MCpA(sen.proc,fpr.proc)
    pAUC[[i-1]] <- as.vector(pA(fpr.proc,sen.proc))
    sensitivity[[i-1]] <-as.vector(sen.proc)
    FPR[[i-1]] <- as.vector(fpr.proc)
    if (isTRUE(plot)) {plot(sen.roc~fpr.roc, type="l", col=i, ylab="TPR", xlab="FPR")
      legend(x= "bottomright",legend = name.variable[2:i], fill = 2:i, cex = 0.8)
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

