#' @title mcpAUCboot
#' @description Calculates the confidence interval using a boot analysis
#' @param dataset dataframe or RangedSummarizedExperiment objetc
#' @param r number of iterations.
#' @param type.interval String that represent the type of intervals required.
#' The value should be any subset of the values
#' c("norm","basic", "stud", "perc", "bca")
#' or simply "all" which will compute all five types of intervals.
#' @param low.value lower false positive rate value that the function
#'  will use to calculate the pAUC
#' @param up.value upper false positive rate value that the function
#'  will use to calculate the pAUC
#' @param selection  vector that will only be used if the parameter
#' "dataset" is a RangedSummarizedExperiment object.
#' This parameter is used to select the variables that will be analysed
#' @param variable in case that dataset is a  SummarizedExperiment,
#' indicate the Gold Standard
#' @param level confidence level
#' @return SummarizedExperiment object with the mcpAUC, the standard
#'  desviation, and the lower and upper limits of the confidence interval.
#' @export mcpAUCboot
#' @import boot
#' @examples
#'library(fission)
#'data("fission")
#'resultsMCboot <- mcpAUCboot(fission,low.value = 0, up.value = 0.25,
#' selection = c("SPNCRNA.1080","SPAC186.08c"), variable="strain")


mcpAUCboot <- function(dataset,  low.value = NULL, up.value = NULL,
                       r=50, level = 0.95, type.interval="perc", selection = NULL, variable=NULL) {

  ci_mcpAUC <- NULL; sd <- NULL; par <- NULL; legend <- NULL; abline <- NULL;

  stopifnot(is.data.frame(dataset) || is(dataset, "SummarizedExperiment"),
            is.numeric(low.value), low.value>=0 && low.value <=1,
            is.numeric(up.value), up.value>=0 && up.value <=1,
            is.numeric(level), level<=1 && level>=0,
            is.numeric(r), type.interval == "norm" || type.interval == "perc" || type.interval == "basic" ||
              type.interval == "stud" || type.interval == "bca")
  if (is(dataset, "SummarizedExperiment")) {
    stopifnot(is.character(selection), is.character(variable))
    strain <- dataset@colData@listData
    strain <- strain[variable][[1]]
    dataset <- as.data.frame(SummarizedExperiment::assay(dataset))
    dataset <- scale(t(as.matrix(dataset[selection,])), center=TRUE, scale = TRUE)
    name.variable <- colnames(dataset)
    dataset <- as.data.frame(cbind(strain,dataset))
  } else {  dataset <- as.data.frame(dataset)  }

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

  value <- mcpAUC(dataset, low.value=low.value, up.value = up.value)
  McpNA <-assay(value)$St_pAUC

  result_boot <- boot::boot(dataset, statistic = fbootM, R=r, low.limit = low.limit, up.limit =up.limit)
  intervalo_confianza <- matrix(0,nrow=dim(result_boot$t)[2],ncol=4)
  for (i in seq_along(dim(result_boot$t)[2])) {
    if (is.na(McpNA[i])){intervalo_confianza[i,] =rep(NA,4)} else {
      ci_mcpAUC <- boot::boot.ci(result_boot, type <- type.interval, conf = level, index = i )
      p_max=length(ci_mcpAUC[[4]])
      intervalo_confianza[i,]=c(result_boot$t0[i],ci_mcpAUC[[4]][c(p_max-1,p_max)],sd(result_boot$t[,i]))
    }
  }
  colnames(intervalo_confianza)=c("MCp_AUC","lwr","upr","sd")
  names <- c("MCp_AUC","lwr","upr","sd")
  obj <- list(intervalo_confianza[,1],intervalo_confianza[,2],intervalo_confianza[,3],intervalo_confianza[,4])
  x <- createSE(obj, names)

  return(x)

}

