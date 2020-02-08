#' @title TpAUCboot
#' @description Calculates the confidence interval using a boot analysis
#' @param dataset dataframe or RangedSummarizedExperiment objetc
#' @param r number of iterations.
#' @param type.interval String that represent the type of intervals required.
#' The value should be any subset of the values c("norm","basic", "stud", "perc", "bca")
#' or simply "all" which will compute all five types of intervals.
#' @param seed Seed
#' @param low.value lower false positive rate value that the function will use to calculate the pAUC
#' @param up.value upper false positive rate value that the function will use to calculate the pAUC
#' @param selection  vector that will only be used if the parameter "dataset" is a RangedSummarizedExperiment object.
#' This parameter is used to select the variables that will be analysed
#' @param level confidence level
#' @param variable in case that dataset is a  SummarizedExperiment, indicate the Gold Standard
#' @return SummarizedExperiment object with the Tp_AUC, the standard desviation, and the lower and upper limits of the confidence interval
#' @export TpAUCboot
#' @import boot
#' @examples
#'library(fission)
#'data("fission")
#'resultstboot<- TpAUCboot(fission,low.value = 0, up.value = 0.25, seed = 1234, selection = c("SPNCRNA.1080","SPAC186.08c"), varidable="strain")


TpAUCboot <- function(dataset,  low.value = NULL, up.value = NULL,
                      r=50, seed=NULL, level = 0.95, type.interval="perc", selection = NULL, variable=NULL) {

  ci_TpAUC <- NULL; CpA=NULL;ci_MCpAUC <- NULL; sd <- NULL; par <- NULL; legend <- NULL; abline <- NULL;

  if (class(dataset)=="RangedSummarizedExperiment") {
    strain <- dataset@colData@listData
    strain <- strain[variable][[1]]
    dataset <- as.data.frame(SummarizedExperiment::assay(dataset))
    dataset <- scale(t(as.matrix(dataset[selection,])), center=TRUE, scale = TRUE)
    name.variable <- colnames(dataset)
    dataset <- as.data.frame(cbind(strain,dataset))
  }  else {  dataset <- as.data.frame(dataset)  }

  if(dim(dataset)[2]<2) {stop("database has to have at least 2 colums")}
  if (!is.null(up.value)){
    up.limit <- up.value
  }else{
    up.limit <- 1
  }
  if (!is.null(low.value)){
    low.limit <- low.value
  }else{
    low.limit <- 0
  }
  if(is.null(seed)){
    old.seed <- .Random.seed
    on.exit({.Random.seed <- old.seed})
  }

  result_boot <- boot::boot(dataset, statistic = fbootT, R=r, low.limit = low.limit, up.limit =up.limit)

  intervalo_confianza <- matrix(0,nrow=dim(result_boot$t)[2],ncol=4)
  for (i in 1:dim(result_boot$t)[2]) {
    ci_TpAUC <- boot::boot.ci(result_boot, type <- type.interval, conf = level, index = i )
    p_max=length(ci_TpAUC[[4]])
    intervalo_confianza[i,]=c(result_boot$t0[i],ci_TpAUC[[4]][c(p_max-1,p_max)],sd(result_boot$t[,i]))
  }
  colnames(intervalo_confianza)=c("Tp_AUC","lwr","upr","sd")
  names <- c("Tp_AUC","lwr","upr","sd")
  obj <- list(intervalo_confianza[,1],intervalo_confianza[,2],intervalo_confianza[,3],intervalo_confianza[,4])
  x <- createSE(obj, names)

  return(x)

}
