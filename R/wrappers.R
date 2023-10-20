##############################################################################################
## GEE Power Calculation Wrapper
##
##############################################################################################


#' @title Perform the GEE-based power calculation for stepped wedge trials
#' and multilevel cluster randomized trials. Wrapper
#'
#' @description
#'

CalcSWTPower <- function(design_mat_list = NULL,
                         dgm_param_vec = NULL,
                         n_clust_per_seq = NULL,
                         n_ind_per_clust = NULL,
                         clust_trt_effect = NULL,
                         ind_trt_effect = NULL,
                         interaction_effect = NULL,
                         period_effect_vec = NULL,
                         incr_effect_scalar = NULL,
                         var_fun,
                         cor_mat_list = NULL,
                         incidence_mat_list = NULL,
                         dist = "normal",
                         link = "identity",
                         dispersion_scalar_scalar = NULL,
                         alpha = .05,
                         test_df,
                         contrast_mat = NULL) {

  if(is.null(design_mat_list)) design_mat_list <- CreateClusterCompleteDesignMatrixList(n_study_time_periods,
                                                                                       n_clust_per_sequence,
                                                                                       n_ind_per_clust,
                                                                                       time_effects_type,
                                                                                       trt_effects_type,
                                                                                       linear_trt_scale_factor = 1,
                                                                                       incidence_matrix_list = NULL)

    if(is.null(cor_mat_list)) cor_mat_list <- CreateClusterCompleteCorMatList(n_study_time_periods,
                                                                              n_clust_per_sequence,
                                                                              n_ind_per_clust,
                                                                              within_period_corr,
                                                                              between_period_corr,
                                                                              within_subject_correlation)

  swt_power <- .CalcPower(design_mat_list,
             working_cor_mat_list,
             incidence_mat_list = NULL,
             trt_param_col_vec,
             dispersion_scalar,
             var_fun = var_fun,
             link,
             response_type,
             alpha,
             test_df,
             contrast_mat,
             null_val_vec,
             power_only_flag = FALSE)
  stop("This function is not yet implemented")
}
