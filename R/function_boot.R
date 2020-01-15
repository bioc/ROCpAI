#' @title funcion.boot
#' @description Calculates the confidence interval using a boot analysis
#' @param bsdatabase Dataframe of the complete information of the samples
#' @param r number of iterations.
#' @param funcion function
#' @param type.interval String that represent the type of intervals required. The value should be any subset of the values c("norm","basic", "stud", "perc", "bca") or simply "all" which will compute all five types of intervals.
#' @param seed Seed
#' @return
#' @export funcion.boot
#' @import boot
#' @examples
#'
fboot<- function(dataset,bssample, low.limit, up.limit){
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
TpAUCboot <- function(dataset,  low.value = NULL, up.value = NULL,
                      r=50, seed=NULL, level = 0.95, type.interval="perc") {

  ci_TpAUC <- NULL; CpA=NULL;

  CpA <- TpA

  if (class(dataset)=="RangedSummarizedExperiment") {
    dataset <- as.data.frame(SummarizedExperiment::assay(dataset))
  } else {
    dataset <- as.data.frame(dataset)}
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
  }else{
    set.seed(as.integer(seed))
  }

  result_boot <- boot::boot(dataset, statistic = fboot, R=r, low.limit = low.limit, up.limit =up.limit)

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


MCpAUCboot <- function(dataset,  low.value = NULL, up.value = NULL,
                      r=50, seed=NULL, level = 0.95, type.interval="perc") {

  ci_TpAUC <- NULL; CpA=NULL;

  CpA <- MCpA

  if (class(dataset)=="RangedSummarizedExperiment") {
    dataset <- as.data.frame(SummarizedExperiment::assay(dataset))
  } else {
    dataset <- as.data.frame(dataset)}
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
  }else{
    set.seed(as.integer(seed))
  }

  result_boot <- boot::boot(dataset, statistic = fboot, R=r, low.limit = low.limit, up.limit =up.limit)

  intervalo_confianza <- matrix(0,nrow=dim(result_boot$t)[2],ncol=4)
  for (i in 1:dim(result_boot$t)[2]) {
    ci_TpAUC <- boot::boot.ci(result_boot, type <- type.interval, conf = level, index = i )
    p_max=length(ci_TpAUC[[4]])
    intervalo_confianza[i,]=c(result_boot$t0[i],ci_TpAUC[[4]][c(p_max-1,p_max)],sd(result_boot$t[,i]))
  }
  colnames(intervalo_confianza)=c("MCp_AUC","lwr","upr","sd")
  names <- c("MCp_AUC","lwr","upr","sd")
  obj <- list(intervalo_confianza[,1],intervalo_confianza[,2],intervalo_confianza[,3],intervalo_confianza[,4])
  x <- createSE(obj, names)

  return(x)

}
