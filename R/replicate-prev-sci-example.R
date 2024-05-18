#' @title Replicate Example Published in Prevention Science
#'
#' @description Replicates the example published in Prevention Science.
#'
#' @export
ReplicatePrevSciPaperExample <- function(){

  # Design matrices
  n_study_periods <- 6
  n_clust_trt_seqs <- 5
  n_clust_per_seq <- 10
  crossover_time_period_vec <- .CreateCrossoverTimePeriodVec(n_study_periods = n_study_periods,
                                                             n_clust_trt_seqs = n_clust_trt_seqs)
  n_ind_per_clust <- 9
  time_model_type <- "linear"
  trt_model_type <- "linear"
  linear_trt_scale_factor <- 1/3
  mli_study_flag <- TRUE

  # Working correlation matrices
  sample_type <- "cohort"
  within_period_cor <-  .05
  between_period_cor <-  .025
  within_subject_cor <-  .4
  cor_decay_rate <- NULL

  # Data Generating Model
  time_intercept_param <- 2
  time_trend_param <-  -.05
  period_effect_param_vec <-  NULL
  clust_trt_param <- .15
  ind_trt_param <-  .1
  interaction_param <-  -.05

  # Power Calculation
  dispersion_scalar <-  1
  family <- gaussian()
  alpha <- .05

  ##

  design_mat_list <- CreateClusterCompleteDesignMatrixList(n_study_periods = n_study_periods,
                                                           n_clust_trt_seqs = n_clust_trt_seqs,
                                                           n_clust_per_seq = n_clust_per_seq,
                                                           crossover_time_period_vec = crossover_time_period_vec,
                                                           n_ind_per_clust = n_ind_per_clust,
                                                           time_model_type = time_model_type,
                                                           trt_model_type = trt_model_type,
                                                           linear_trt_scale_factor = linear_trt_scale_factor,
                                                           mli_study_flag = mli_study_flag)

  n_mean_model_params <- ncol(design_mat_list[[1]])

  test_df <- .CalculateDefaultTestDf(n_clusters = n_clust_trt_seqs*n_clust_per_seq,
                                     n_mean_model_params = n_mean_model_params)
  contrast_mat <- matrix(c(0, 0, 1, 1, 1), nrow = 1, ncol = 5)
  null_val_vec <- .CreateDefaultNullValVec(mli_study_flag = mli_study_flag)


  working_cor_mat_list <- CreateCompleteWorkingCorMatList(
    n_study_periods = n_study_periods,
    n_clusters = n_clust_trt_seqs*n_clust_per_seq,
    n_ind_per_clust = n_ind_per_clust,
    sample_type = "cohort",
    within_period_cor = within_period_cor,
    between_period_cor = between_period_cor,
    within_subject_cor = within_subject_cor,
    cor_decay_rate = cor_decay_rate
  )

  set.seed(42)

  dgm_param_vec <- SpecifyDataGeneratingModel(n_study_periods = n_study_periods,
                                              mli_study_flag = mli_study_flag,
                                              time_model_type = time_model_type,
                                              time_intercept_param = time_intercept_param,
                                              time_trend_param = time_trend_param,
                                              period_effect_param_vec = period_effect_param_vec,
                                              trt_model_type = trt_model_type,
                                              clust_trt_param = clust_trt_param,
                                              ind_trt_param = ind_trt_param,
                                              interaction_param = interaction_param)

  dgm_param_col_vec <- matrix(dgm_param_vec, ncol = 1)

  power_list <- .CalcPower(design_mat_list = design_mat_list,
                           working_cor_mat_list = working_cor_mat_list,
                           dgm_param_col_vec = dgm_param_col_vec,
                           dispersion_scalar = dispersion_scalar,
                           family = family,
                           alpha = alpha,
                           test_df = 1,
                           contrast_mat = contrast_mat,
                           null_val_vec = null_val_vec,
                           power_only_flag = FALSE
  )

  return(power_list)
}


ReplicatePrevSciPaperExample_CalcSampleSize <- function(n_clust_per_seq=10, n_ind_per_clust=9){
  
  # Design matrices
  n_study_periods <- 6
  n_clust_trt_seqs <- 5
  n_clust_per_seq <- n_clust_per_seq
  crossover_time_period_vec <- .CreateCrossoverTimePeriodVec(n_study_periods = n_study_periods,
                                                             n_clust_trt_seqs = n_clust_trt_seqs)
  n_ind_per_clust <- n_ind_per_clust
  time_model_type <- "linear"
  trt_model_type <- "linear"
  linear_trt_scale_factor <- 1/3
  mli_study_flag <- TRUE
  
  # Working correlation matrices
  sample_type <- "cohort"
  within_period_cor <-  .05
  between_period_cor <-  .025
  within_subject_cor <-  .4
  cor_decay_rate <- NULL
  
  # Data Generating Model
  time_intercept_param <- 2
  time_trend_param <-  -.05
  period_effect_param_vec <-  NULL
  clust_trt_param <- .15
  ind_trt_param <-  .1
  interaction_param <-  -.05
  
  # Power Calculation
  dispersion_scalar <-  1
  family <- gaussian()
  alpha <- .05
  
  ##
  
  design_mat_list <- CreateClusterCompleteDesignMatrixList(n_study_periods = n_study_periods,
                                                           n_clust_trt_seqs = n_clust_trt_seqs,
                                                           n_clust_per_seq = n_clust_per_seq,
                                                           crossover_time_period_vec = crossover_time_period_vec,
                                                           n_ind_per_clust = n_ind_per_clust,
                                                           time_model_type = time_model_type,
                                                           trt_model_type = trt_model_type,
                                                           linear_trt_scale_factor = linear_trt_scale_factor,
                                                           mli_study_flag = mli_study_flag)
  
  n_mean_model_params <- ncol(design_mat_list[[1]])
  
  test_df <- .CalculateDefaultTestDf(n_clusters = n_clust_trt_seqs*n_clust_per_seq,
                                     n_mean_model_params = n_mean_model_params)
  contrast_mat <- matrix(c(0, 0, 1, 1, 1), nrow = 1, ncol = 5)
  null_val_vec <- .CreateDefaultNullValVec(mli_study_flag = mli_study_flag)
  
  
  working_cor_mat_list <- CreateCompleteWorkingCorMatList(
    n_study_periods = n_study_periods,
    n_clusters = n_clust_trt_seqs*n_clust_per_seq,
    n_ind_per_clust = n_ind_per_clust,
    sample_type = "cohort",
    within_period_cor = within_period_cor,
    between_period_cor = between_period_cor,
    within_subject_cor = within_subject_cor,
    cor_decay_rate = cor_decay_rate
  )
  
  set.seed(42)
  
  dgm_param_vec <- SpecifyDataGeneratingModel(n_study_periods = n_study_periods,
                                              mli_study_flag = mli_study_flag,
                                              time_model_type = time_model_type,
                                              time_intercept_param = time_intercept_param,
                                              time_trend_param = time_trend_param,
                                              period_effect_param_vec = period_effect_param_vec,
                                              trt_model_type = trt_model_type,
                                              clust_trt_param = clust_trt_param,
                                              ind_trt_param = ind_trt_param,
                                              interaction_param = interaction_param)
  
  dgm_param_col_vec <- matrix(dgm_param_vec, ncol = 1)
  
  power_list <- .CalcPower(design_mat_list = design_mat_list,
                           working_cor_mat_list = working_cor_mat_list,
                           dgm_param_col_vec = dgm_param_col_vec,
                           dispersion_scalar = dispersion_scalar,
                           family = family,
                           alpha = alpha,
                           test_df = 1,
                           contrast_mat = contrast_mat,
                           null_val_vec = null_val_vec,
                           power_only_flag = FALSE
  )
  
  return(power_list)
}