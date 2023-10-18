#' @title Perform the GEE-based power calculation for stepped wedge trials
#' and multilevel cluster randomized trials. Wrapper
#'
#' @description
#'

CalcGEEPower <- function(design_mat_list = NULL,
                         param_vec = NULL,
                         n_clust_per_seq = NULL,
                         clust_trt_effect = NULL,
                         ind_trt_effect = NULL,
                         interaction_effect = NULL,
                         period_effect_vec = NULL,
                         incr_effect_scalar = NULL,
                         cor_mat_list = NULL,
                         incidence_mat_list = NULL,
dist = "normal",
                         link = "identity",
dispersion,
                         alpha = .05,
test_df){

}

#' @title Internal Power Calculation Function
#'
#' @description
#' Perform the GEE-based power calculation for stepped wedge trials and
#' multilevel cluster randomized trials. This function is the internal function
#' for performing the power calculation. It is not exported.
#'
#' @param design_mat_list A list of design matrices, one for each cluster.
#' @param incidence_mat_list A list of incidence matrices, one for each cluster.
#' @param cor_mat_list A list of correlation matrices, one for each cluster.
#' @param trt_param_col_vec A vector of treatment parameters as
#' a matrix with one column.
#' @param dispersion The dispersion parameter for the GEE model.
#' @param link The link function for the GEE model.
#' @param dist The distribution for the GEE model.
#' @param alpha The desired type I error rate.
#' @param test_df The degrees of freedom for the test statistic.
#'
#' @return A list with the following components:
#' \item{power}{The power of the test.}
#' \item{sample_size}{The sample size required to achieve the desired power.}
#' \item{test_stat}{The test statistic.}
#' \item{p_value}{The p-value for the test statistic.}
#' \item{test_df}{The degrees of freedom for the test statistic.}
#'
#'
.CalcPower <- function(design_mat_list,
                       incidence_mat_list,
                       cor_mat_list,
                       trt_param_col_vec,
                       dispersion,
                       link,
                       dist,
                       alpha,
                       test_df){


}

.CreateEEMeanComponents <- function(){

}

.CreatEEVarComponents <- function (){

}