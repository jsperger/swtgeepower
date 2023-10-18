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
#' @param n_study_time_periods The total number of time periods in the study.
#' @param n_obs_per_sequence
#' @param time_effects_type
#' @param trt_effects_type
#' @param linear_trt_scale_factor
#'
#' @return A matrix
#'
#' @export

CreateCompleteStudyDesignMatrix <- function(n_study_time_periods,
                                     n_obs_per_sequence,
                                     time_effects_type,
                                     trt_effects_type,
                                     linear_trt_scale_factor = 1,
                                     incidence_matrix_list = NULL){
  n_sequences <- n_study_time_periods - 1



}