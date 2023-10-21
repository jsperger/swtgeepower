#' @title Create study design matrix for a complete study where every cluster
#' observed in every period.
#'
#' @description
#' r lifecycle::badge("experimental")`
#'
#' Simple function that assumes that each cluster is observed in each time
#' period, and that every sequence has at least one period on control before
#' crossing over to treatment.
#'
#' @param n_study_periods The total number of time periods in the study.
#' @param n_obs_per_sequence
#' @param time_effects_type
#' @param trt_effects_type
#' @param linear_trt_scale_factor
#'
#' @return A matrix
#'
#' @export

CreateClusterCompleteDesignMatrixList <- function(
  n_study_periods,
  n_clust_trt_seqs,
  n_clust_per_seq,
  crossover_time_period,
  n_ind_per_clust,
  time_model_type,
  treatment_model_type,
  linear_trt_scale_factor = 1,
  mli_study_flag = FALSE){
  
  design_mat_list <- vector("list", length = n_clust_trt_seqs*n_clust_per_seq
  
  for(i in 1:n_clust_trt_seqs){
    seq_lower_clust_index <- (i - 1)*n_clust_per_seq + 1
    seq_upper_clust_index <- i*n_clust_per_seq
    design_mat_list[[seq_lower_clust_index:seq_upper_clust_index]] <- CreateClusterCompleteDesignMatrix(n_study_periods = n_study_periods,
                                                              mli_study_flag = mli_study_flag,
                                                              crossover_time_period = crossover_time_period[i],
                                                              n_ind_per_clust = n_ind_per_clust,
                                                              time_model_type = time_model_type,
                                                              treatment_model_type = treatment_model_type,
                                                              linear_trt_scale_factor = linear_trt_scale_factor)
  }
  
  return(design_mat_list)

}

.CreateClusterCompleteDesignMatrix <- function(n_study_periods,
                                               mli_study_flag,
                                               crossover_time_period,
                                        n_ind_per_clust,
                                        time_effects_type,
                                        trt_effects_type,
                                        linear_trt_scale_factor = 1) {

  dim_time_params <- switch(time_effects_type,
                            linear = 2,
                            categorical = n_study_periods,
                            quadratic = stop("Not implemented yet [will be 3]"),
  )

  dim_trt_params <- switch(trt_effects_type,
                           average = 1,
                           linear = 1,
                           quadratic = stop("Not implemented yet [will be 2]")
  )

  if(mli_study_flag == TRUE) dim_trt_params <- 3*dim_trt_params

  n_params <- dim_time_params + dim_trt_params

  time_period_vec <- 1:n_study_periods

  trt_time_vec <- c(rep(0, times = crossover_time_period - 1),
                    seq(from = 1, to = n_study_periods - crossover_time_period + 1, by = 1))

  trt_time_mat <- matrix(1, nrow = dim_trt_params, ncol = 1) %x% matrix(trt_time_vec, ncol = 1)

  # Matrix Dimension will depend on the type of time effects
  transformed_time_mat <- switch(time_effects_type,
                                 linear = matrix(time_period_vec, ncol = 1),
                                 quadratic = matrix(c(time_period_vec, time_period_vec^2), ncol = 2),
                                 period = diag(n_study_periods)
  )

  # Matrix Dimension will depend on the type of study design
  # Note for quadratic that x^2 for a matrix in R does element-wise squaring
  transformed_trt_mat <- switch(trt_effects_type,
                                average = pmin(trt_time_mat, 1),
                                linear = linear_trt_scale_factor * trt_time_mat,
                                quadratic = cbind(trt_time_mat, trt_time_mat^2)
  )

  # Create the design matrix

  single_ind_design_mat <- cbind(rep(1, times = n_study_periods),
                         transformed_time_mat,
                         transformed_trt_mat)

  cluster_design_mat <- matrix(1, nrow = n_ind_per_clust, ncol = 1) %x%
    single_ind_design_mat

  stopifnot(all.equal(dim(cluster_design_mat), c(n_ind_per_clust, n_params + 1)))

  colnames(cluster_design_mat) <- c(.NameTimeEffects(time_effects_type), .NameTrtEffects(mli_study_flag))

  return(cluster_design_mat)
}