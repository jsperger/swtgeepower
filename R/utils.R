#' Description: Utility functions for the project

###############################################################################
# Utility Functions - Defaults
#
###############################################################################

# Default Names
###############################################################################
#' @title Names for the Time Effects Parameters
#'
#' @description Gets its own function because it's used in the model parameter
#' and design matrix functions.
.NameTimeEffects <- function(time_model_type, n_study_periods){
  time_effect_names <- switch(time_model_type,
                              linear = c("Intercept", "Period"),
                              categorical = c("Intercept", paste0("Period", 2:n_study_periods))
  )

  return(time_effect_names)
}

.NameTrtEffects <- function(mli_study_flag){
  trt_effect_names <- ifelse(mli_study_flag == TRUE,
                             yes = c("TrtClust", "TrtInd", "TrtClustInd"),
                             no = "TrtClust")

  return(trt_effect_names)
}

# Default Parameter Dimensions
###############################################################################

.DefaultDimTimeParams <- function(time_model_type, n_study_periods){
  # Note: this includes the intercept

  dim_time_params <- switch(time_model_type,
                            linear = 2,
                            categorical = n_study_periods,
                            quadratic = stop("Not implemented yet [will be 3]"),
  )

  return(dim_time_params)
}

.DefaultDimTrtParams <- function(mli_study_flag, trt_model_type){
  # Assumes the intercept is included as part of the time parameters
  dim_trt_params <- switch(trt_model_type,
                           average = 1,
                           linear = 1,
                           quadratic = stop("Not implemented yet [will be 2]")
  )

  if(mli_study_flag == TRUE) dim_trt_params <- 3*dim_trt_params

  return(dim_trt_params)
}

#' @title Create an integer vector whose entries are the unique calendar time
#' periods at which the intervention is introduced to a cluster(s).
#'
#' @description
#' Tries to guess the calendar time periods at which the intervention is
#' introduced based on the number of study periods and the number of randomization
#' sequences.
#'
#' There are three options currently implemented:
#' - If the number of study periods is equal to the number of randomization sequences,
#' then every time period is in the vector of crossover times.
#' - If the number of study periods is one more than the number of randomization sequences,
#' then every time period except the first is in the vector of crossover times.
#' - If the number of study periods is greater than the number of randomization sequences
#' by more than one, then the crossover times are equally spaced starting in period two.
#' Throws an error if the number of study periods minus one is not evenly divisible by
#' the number of randomization sequences.
#'
#' @param n_study_periods The number of study periods
#' @param n_clust_trt_seqs The number of randomization sequences
.CreateCrossoverTimePeriodVec <- function(n_study_periods,
                                          n_clust_trt_seqs){
  # There is one sequence that only receives the intervention and never the control condition
  if(n_study_periods == n_clust_trt_seqs) return(1:n_study_periods)

  # Assumes that the first period is control only
  if(n_study_periods == n_clust_trt_seqs + 1) return(2:n_study_periods)

  # Assumes that the first period is always control
  equally_spaced_flag <- ((n_study_periods - 1) %% n_clust_trt_seqs) == 0

  if(equally_spaced_flag == TRUE){
    crossover_time_period_vec <- seq(from = 2,
                                     to = n_study_periods,
                                     by = (n_study_periods - 1)/n_clust_trt_seqs)

    return(crossover_time_period_vec)
  }

  stop("n_clust_trt_seqs must be equal to one of \n n_study_periods, \n n_clust_trt_seqs + 1, \n or a factor of n_study_periods - 1")

  return(NULL)
}


.CreateDefaultNullValVec <- function(mli_study_flag){
  contrast_dim <- ifelse(mli_study_flag == TRUE, 3, 1)

  return(rep(0, times = contrast_dim))
}

#' @title Create a default contrast matrix
#'
.CreateDefaultContrastMat <- function(n_time_model_params,
                                      n_trt_model_params) {
  contrast_row_dim <- n_trt_model_params

  contrast_mat <- cbind( matrix(0, nrow = contrast_row_dim, ncol = n_time_model_params),
                         diag(1, nrow = n_trt_model_params))

  return(contrast_mat)
}

.CalculateDefaultTestDf <- function(n_clusters,
                                 n_mean_model_params){
  return(n_clusters - n_mean_model_params)
}




###############################################################################
# Utility Functions - Common Checks
#
###############################################################################

.CheckFrechetBound <- function(){
  stop("Not implemented yet")
}