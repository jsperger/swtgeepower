#' swtgeepower
#'
#' Power and sample size calculations for stepped wedge cluster randomized trials and multilevel cluster randomized trials
#'
#' @docType package
#' @name swtgeepower


##############################################################################################
## Wrapper to Specify All Incidence Matrices
##
##############################################################################################

#' @title Create Incidence Matrices that define the missingness/observation structure of the data.
#'
#' @description
#' r lifecycle::badge("experimental")`
#'
#' @return list of incidence matrices
CreateIncidenceMatricesList <- function (n_study_periods,
                                         n_obs_periods){
  n_unique_sequences <- n_study_periods - n_obs_periods + 1

  incidence_mats_list <- vector(mode = "list", length = n_unique_sequences)

  for(i in 1:n_unique_sequences){
    incidence_mats_list[[i]] <- CreateSequenceIncidenceMat(start_period = i,
                                                 end_period = i + n_obs_periods - 1,
                                                 n_study_periods = n_study_periods)
  }

  return(incidence_mats_list)
}

##############################################################################################
## Specify individual Incidence Matrices
##
##############################################################################################

#' @title Generate the incidence matrix for a specific sequence (applies to all clusters in that sequence)
#'
#' @description
#' r lifecycle::badge("experimental")`
#' The incidence matrix defines the time periods in which each cluster is observed.
#' Each row corresponds to a time period the cluster is observed, and each column corresponds to a time period in the study.
#'
#'
#' @param start_period The first time period in which the cluster is observed
#' @param end_period The last time period in which the cluster is observed
#' @param n_study_periods The total number of time periods in the study.
#' @param unobserved_periods An optional vector of time periods in which the cluster is not observed. Only necessary if the cluster is not observed in time periods between the start and end period.
#' Not required if unobserved periods are at the beginning or end of the study.
#'
#' @return A matrix with one row for each time period the cluster is observed and .
#'
#' @examples
#' CreateIncidenceMat(start_period = 1,
#'                   end_period = 3,
#'                  n_study_periods = 5)
#'
#' CreateIncidenceMat(start_period = 3,
#'                  end_period = 7,
#'                 n_study_periods = 10,
#'                unobserved_periods = c(4, 5))
#'
#' @export
CreateSequenceIncidenceMat <- function(start_period,
                               end_period,
                               n_study_periods,
                               unobserved_periods = NULL){
  stopifnot(length(start_period) == 1)
  stopifnot(length(end_period) == 1)

  start_period <- as.integer(start_period)
  end_period <- as.integer(end_period)

  if(!is.null(unobserved_periods)) unobserved_periods <- as.integer(unobserved_periods)

  # Check that the start and end periods are valid
  if (start_period > end_period) {
    stop("ERROR: start_period must be less than or equal to end_period")
  }
  if (start_period < 1 || start_period > n_study_periods) {
    stop("ERROR: start_period must be between 1 and n_study_periods")
  }
  if (end_period < 1 || end_period > n_study_periods) {
    stop("ERROR: end_period must be between 1 and n_study_periods")
  }

  # Check that the unobserved periods are valid
  if (!is.null(unobserved_periods)) {
    if (any(unobserved_periods < 1) || any(unobserved_periods > n_study_periods)) {
      stop("ERROR: unobserved_periods must be between 1 and n_study_periods")
    }
    if (any(unobserved_periods < start_period) || any(unobserved_periods > end_period)) {
      warning("Warning: unobserved_periods are unnecessary except for periods unobserved between start_period and end_period. Verify that the input periods are correct.")
    }
  }

  # Create the incidence matrix
  n_obs_periods <- 1 + (end_period - start_period) - length(unobserved_periods)

  observed_periods <- setdiff(start_period:end_period, unobserved_periods) # If NULL setdiff returns the first argument

  incidence_mat <- matrix(0, nrow = n_obs_periods, ncol = n_study_periods)

  # Fill in the observed periods with 1s
  # Takes unobserved periods into account

  incidence_mat[1:n_obs_periods, observed_periods] <- diag(1, n_obs_periods)

  return(incidence_mat)
}
