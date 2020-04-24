# !/usr/bin/env Rscript
# 04/23/2020
#
#' @title CUE
#' 
#' @description This function CUE.impute is designed to perform CpG impUtation Ensemble (CUE), which leverages multiple classical statistical and modern machine learning methods, to impute from the Illumina HumanMethylation450 (HM450) BeadChip to the Illumina HumanMethylationEPIC (HM850) BeadChip.
#' It takes as input the beta methylation matrix for the probse of the Illumina HumanMethylation450 (HM450) BeadChip.
#' It outputs the methylation matrix of for imputed probes of the Illumina HumanMethylationEPIC (HM850) BeadChip.
#' @usage CUE.impute(HM450)
#' @param HM450 is matrix of XXXK HM450 probes.
#'
#' @return SMNNcorrect returns the HM850 imputed probes
#' @author Gang Li <franklee@live.unc.edu>, Laura Raffield, Yun Li <yunli@med.unc.edu>
#' @examples 
#' # Load the example data data_CUE
#' data("data_CUE")
#' 
#' # Perform imputation using CUE.impute
#' HM850.impute <- CUE.impute(HM450)
#' @import 
#' @export
CUE.impute <- function(HM450){

  return (HM850.impute)
}
