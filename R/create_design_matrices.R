#' @title Create Cluster Complete Design Matrix List
#'
#' @description
#' Creates a list of design matrices for multiple clusters over a given number of study periods.
#' These matrices account for the study design, treatment sequences, time effects, and other relevant factors.
#'
#' @param n_study_periods Integer. The number of study periods.
#' @param n_clust_trt_seqs Integer. The number of cluster treatment sequences.
#' @param n_clust_per_seq Integer. The number of clusters per sequence.
#' @param crossover_time_period_vec Numeric vector. The time periods for crossover.
#' @param n_ind_per_clust Integer. The number of individuals per cluster.
#' @param time_model_type Character. The type of time model ("linear" or "categorical").
#' @param trt_model_type Character. The type of treatment model ("average" or "linear").
#' @param linear_trt_scale_factor Numeric. The scaling factor for linear treatment effects. Default is 1.
#' @param mli_study_flag Logical. Flag indicating if the study is a multi-level intervention study. Default is FALSE.
#'
#' @return A list of design matrices.
CreateClusterCompleteDesignMatrixList <- function(
  n_study_periods,
  n_clust_trt_seqs,
  n_clust_per_seq,
  crossover_time_period_vec,
  n_ind_per_clust,
  time_model_type,
  trt_model_type,
  linear_trt_scale_factor = 1,
  mli_study_flag = FALSE){

  .CheckClusterCompleteDesignMatrixInputs(
    n_study_periods = n_study_periods,
    n_clust_trt_seqs = n_clust_trt_seqs,
    n_clust_per_seq = n_clust_per_seq,
    crossover_time_period_vec = crossover_time_period_vec,
    n_ind_per_clust = n_ind_per_clust,
    time_model_type = time_model_type,
    trt_model_type = trt_model_type,
    linear_trt_scale_factor = linear_trt_scale_factor,
    mli_study_flag = mli_study_flag)
  
  design_mat_list <- vector("list", length = n_clust_trt_seqs*n_clust_per_seq)
  
  for (i in 1:n_clust_trt_seqs) {
    seq_lower_clust_index <- (i - 1)*n_clust_per_seq + 1
    seq_upper_clust_index <- i*n_clust_per_seq
    cur_seq_indices <- seq(seq_lower_clust_index, seq_upper_clust_index)

    sequence_design_mat <- .CreateClusterCompleteDesignMatrix(
      n_study_periods = n_study_periods,
      mli_study_flag = mli_study_flag,
      crossover_time_period = crossover_time_period_vec[i],
      n_ind_per_clust = n_ind_per_clust,
      time_model_type = time_model_type,
      trt_model_type = trt_model_type,
      linear_trt_scale_factor = linear_trt_scale_factor)

    design_mat_list[c(cur_seq_indices)] <- list(sequence_design_mat)
  }

  .CheckClusterCompleteDesignMatrixOutput(design_mat_list)

  return(design_mat_list)

}

#' @title Create Cluster Complete Design Matrix
#'
#' @description
#' Creates a design matrix for a single cluster over a given number of study periods.
#' This matrix accounts for the study design, treatment sequences, time effects, and other relevant factors.
#'
#' @param n_study_periods Integer. The number of study periods.
#' @param mli_study_flag Logical. Flag indicating if the study is a multi-level intervention study.
#' @param crossover_time_period Integer. The time period for crossover.
#' @param n_ind_per_clust Integer. The number of individuals per cluster.
#' @param time_model_type Character. The type of time model ("linear" or "categorical").
#' @param trt_model_type Character. The type of treatment model ("average" or "linear").
#' @param linear_trt_scale_factor Numeric. The scaling factor for linear treatment effects. Default is 1.
#'
#' @return A matrix that serves as the design matrix for a single cluster.
#'
#' @keywords internal
.CreateClusterCompleteDesignMatrix <- function(n_study_periods,
                                               mli_study_flag,
                                               crossover_time_period,
                                        n_ind_per_clust,
                                        time_model_type,
                                        trt_model_type,
                                        linear_trt_scale_factor = 1) {

  checkmate::expect_int(crossover_time_period)

  # Note: this includes the intercept
  dim_time_params <- .DefaultDimTimeParams(time_model_type = time_model_type,
                                           n_study_periods = n_study_periods)

  # Assumes the intercept is included as part of the time parameters
  dim_trt_params <- .DefaultDimTrtParams(mli_study_flag = mli_study_flag,
                                         trt_model_type = trt_model_type)

  # Define the number of observations and parameters (design matrix row and column dimensions)
  n_observations <- n_study_periods*n_ind_per_clust
  n_params <- dim_time_params + dim_trt_params

  # Create the time and treatment matrices
  time_period_vec <- 1:n_study_periods

  trt_time_vec <- c(rep(0, times = (crossover_time_period - 1)),
                    seq(from = 1, to = n_study_periods - crossover_time_period + 1, by = 1))

  trt_time_mat <- matrix(1, nrow = dim_trt_params, ncol = 1) %x% matrix(trt_time_vec, ncol = 1)

  # Matrix Dimension will depend on the type of time effects
  # There will be an intercept so the first period is not included in categorical
  transformed_time_mat <- switch(time_model_type,
                                 linear = matrix(time_period_vec, ncol = 1),
                                 quadratic = matrix(c(time_period_vec, time_period_vec^2), ncol = 2),
                                 categorical = diag(n_study_periods)[,-1],
                                 stop("time_model_type must be one of 'linear', 'quadratic', or 'categorical'")
  )

  # Matrix Dimension will depend on the type of study design
  # Note for quadratic that x^2 for a matrix in R does element-wise squaring
  transformed_trt_mat <- switch(trt_model_type,
                                average = pmin(trt_time_mat, 1),
                                linear = linear_trt_scale_factor * trt_time_mat,
                                quadratic = cbind(trt_time_mat, trt_time_mat^2),
                                stop("trt_model_type must be one of 'average', 'linear', or 'quadratic'")
  )

  # Create the design matrix

  single_ind_design_mat <- cbind(rep(1, times = n_study_periods),
                         transformed_time_mat,
                         transformed_trt_mat)

  cluster_design_mat <- matrix(1, nrow = n_ind_per_clust, ncol = 1) %x%
    single_ind_design_mat

  expect_equal(nrow(cluster_design_mat), n_observations)
  expect_equal(ncol(cluster_design_mat), n_params)

  colnames(cluster_design_mat) <- c(.NameTimeEffects(time_model_type, n_study_periods), .NameTrtEffects(mli_study_flag))

  return(cluster_design_mat)
}

###############################################################################
## Design Matrix Input Checks
##
###############################################################################
#' @title Check Inputs for CreateClusterCompleteDesignMatrixList
#'
#' @description
#' This function checks the validity of the inputs for the `CreateClusterCompleteDesignMatrixList` function.
#' It uses the `checkmate` package to perform these checks.
#'
#' @details
#' For more information on the parameters, see \code{\link{CreateClusterCompleteDesignMatrixList}}.
#'
#' @return NULL. The function performs checks and will return an error if any check fails.
#' @keywords internal,check
.CheckClusterCompleteDesignMatrixInputs <- function(
  n_study_periods,
  n_clust_trt_seqs,
  n_clust_per_seq,
  crossover_time_period_vec,
  n_ind_per_clust,
  time_model_type,
  trt_model_type,
  linear_trt_scale_factor = 1,
  mli_study_flag = FALSE) {

  checkmate::expect_int(n_study_periods)
  checkmate::expect_int(n_clust_trt_seqs)
  checkmate::expect_int(n_clust_per_seq)
  checkmate::expect_int(n_ind_per_clust)

  checkmate::expect_numeric(crossover_time_period_vec,
                            any.missing = FALSE, null.ok = FALSE)

  checkmate::expect_choice(time_model_type, choices = c("categorical", "linear"))
  checkmate::expect_choice(trt_model_type, choices = c("average", "linear"))

  checkmate::expect_numeric(linear_trt_scale_factor)

  checkmate::expect_flag(mli_study_flag)

  return(NULL)
}

#' @title Check Output for CreateClusterCompleteDesignMatrixList
#'
#' @description
#' This function checks the validity of the output for the `CreateClusterCompleteDesignMatrixList` function.
#' Specifically, it checks the number of columns in each matrix in the list to ensure they are consistent.
#'
#' @details
#' For more information on the expected output structure, see \code{\link{CreateClusterCompleteDesignMatrixList}}.
#'
#' @return NULL. The function performs checks and will return an error if any check fails.
#' @keywords internal,check
.CheckClusterCompleteDesignMatrixOutput <- function(design_mat_list){
  columns_per_matrix <- sapply(design_mat_list, ncol, simplify = TRUE, USE.NAMES = FALSE)

  expect_equal(length(unique(columns_per_matrix)), 1)

  return(NULL)
}