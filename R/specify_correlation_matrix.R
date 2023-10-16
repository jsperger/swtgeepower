#' swtgeepower
#'
#' Power and sample size calculations for stepped wedge cluster randomized trials and multilevel cluster randomized trials
#'
#' @docType package
#' @name swtgeepower

#TODO: explain why there isn't an option for exchangeable correlation structure
#TODO: Add a note about the correlation when there is an implementation period.
#TODO: Add a way to accomodate an implementation period with decay models so that the decay is calculated using calendar time period equivalents

##############################################################################################
## Correlation Matrix Wrapper Function
##
##############################################################################################

#' Create a Correlation Matrix Specification
#'
#' This function generates a specification for a correlation matrix based on
#' a variety of correlation structures, such as exchangeable, exponential decay,
#' and block exchangeable. The function is intended to be used as a utility
#' function to specify the correlation structure for GEE-based power calculations
#' or sample size calculations in stepped wedge trials. The type of correlation structure that is created depends
#' on which arguments are provided.
#'
#' @param n_individuals Integer. The number of individuals within each cluster.
#' @param n_time_periods Integer. The number of time periods for the study.
#' @param design_type String. Specifies whether the design is cross-sectional or cohort. Valid values are "cross" and "cohort".
#'
#' @return An (n_individuals*n_time_periods)-dimensional matrix specifying the correlation structure between observations
#' in the same cluster over time.
#'
#' @examples
#' # For an exchangeable structure
#' corr_matrix_spec(10, 5, "exchangeable", 0.5)
#'
#' # For an exponential decay structure
#' corr_matrix_spec(10, 5, "exponential_decay", 0.5, time_decay_rate = 0.8)
#'
#' @export
#'
SpecifyCorrelationMatrix <- function(
  n_subj_per_period,              # Number of individuals in a cluster
  n_obs_periods,              # Number of observed time periods
  design_type,                # Type of design
  within_period_cor, # Correlation between observations on different subjects in the same cluster in the same time period
  cor_decay_rate = NULL,      # Correlation decary rate over time cor_decay_rate^|time difference|
  between_period_cor = NULL,# Correlation between observations on different subjects in the same cluster in different time periods
  within_subject_cor = NULL # Correlation between observations on the same subject in different time periods
  ) {

  design_type_sanitized <- tolower(design_type) %>% stringr::str_trim(.)

  # Check inputs
  .CheckCorMatInputs(n_obs_periods = n_obs_periods,
                     n_subj_per_period = n_subj_per_period,
                     design_type_sanitized = design_type_sanitized,
                     within_period_cor = within_period_cor,
                     between_period_cor = between_period_cor,
                     within_subject_cor = within_subject_cor,
                     cor_decay_rate = cor_decay_rate)
  
  # Not logically exhaustive, but all other cases are errors that should be caught by input checks
  cor_structure <- dplyr::case_when(
    design_type_sanitized == "cross" & !is.null(cor_decay_rate) ~ "exponential_decay",
    design_type_sanitized == "cross" & !is.null(between_period_cor) & is.null(within_subject_cor) ~ "nested_exchangeable",
    design_type_sanitized == "cohort" & !is.null(cor_decay_rate) ~ "proportional_decay",
    design_type_sanitized == "cohort" & !any(is.null(between_period_cor), is.null(within_subject_cor)) ~ "block_exchangeable",
    TRUE ~ NA_character_
  )

  correlation_matrix <- switch(cor_structure,
                               exponential_decay = .CreateExponentialDecayCorMat(n_obs_periods = n_obs_periods,
                                                                                  n_subj_per_period = n_subj_per_period,
                                                                                  within_period_cor = within_period_cor,
                                                                                  cor_decay_rate = cor_decay_rate),
                               nested_exchangeable = .CreateNestedExchangeableCorMat(n_obs_periods = n_obs_periods,
                                                                                   n_subj_per_period = n_subj_per_period,
                                                                                   within_period_cor = within_period_cor,
                                                                                   between_period_cor = between_period_cor),
                               proportional_decay = .CreateProportionalDecayCorMat(n_obs_periods = n_obs_periods,
                                                                                   n_subj_per_period = n_subj_per_period,
                                                                                   within_period_cor = within_period_cor,
                                                                                   cor_decay_rate = cor_decay_rate),
                               block_exchangeable = .CreateBlockExchangeableCorMat(n_obs_periods = n_obs_periods,
                                                                                   n_subj_per_period = n_subj_per_period,
                                                                                   within_period_cor = within_period_cor,
                                                                                   between_period_cor = between_period_cor,
                                                                                   within_subject_cor = within_subject_cor))

  return(correlation_matrix)
}


##############################################################################################
## Correlation Matrix Creation Functions
##
##############################################################################################

## Cross-sectional Design
###############################################################################################
#' @title Create Nested Exchangeable Correlation Matrix
#'
#' @description
#' Generates a nested-exchangeable correlation matrix based on the given parameters. This function assumes that input validation
#' has been performed by its parent function, 'SpecifyCorrelationMatrix'.
#'
#' @param n_obs_periods Integer. The number of observation periods.
#' @param n_subj_per_period Integer. The number of subjects per observation period.
#' @param within_period_cor Numeric. The correlation coefficient between observations from different subjects
#' in the same cluster within the same time period.
#' @param between_period_cor Numeric. The correlation coefficient between observations from different subjects
#' in the same cluster but in different time periods.
#'
#' @return
#' A nested-exchangeable correlation matrix with dimensions [n_obs_periods * n_subj_per_period, n_obs_periods * n_subj_per_period].
#'
#' @examples
#' \dontrun{
#' .CreateNestedExchangeableCorMat(n_obs_periods = 2, n_subj_per_period = 3, within_period_cor = 0.1, between_period_cor = 0.2)
#' }
#'
.CreateNestedExchangeableCorMat <- function(n_obs_periods,
                                           n_subj_per_period,
                                           within_period_cor,
                                           between_period_cor) {
#TODO: Maybe just replace this with a call to .CreateBlockExchangeableCorMat with within_subject_cor = between_period_cor

  # Create sub-components of the correlation matrix
  diagonal_mat <- diag(1 - within_period_cor,
                       nrow = n_obs_periods * n_subj_per_period)

  within_period_component <- (within_period_cor - between_period_cor) *
    diag(1, n_obs_periods)  %x% matrix(data = 1, nrow = n_subj_per_period, ncol = n_subj_per_period)

  between_period_component <- between_period_cor * matrix(data = 1,
                                                          nrow = n_obs_periods * n_subj_per_period,
                                                          ncol = n_obs_periods * n_subj_per_period)

  # Create nested exchangeable correlation matrix
  nested_exchangeable_cor_mat <-   diagonal_mat + within_period_component + between_period_component

  return(nested_exchangeable_cor_mat)

}

#' @title Create Exponential Decay Correlation Matrix
#'
#' @description This function constructs a correlation matrix based on an exponential decay model. The function
#' assumes that inputs have been validated in the parent function 'SpecifyCorrelationMatrix'.
#'
#' @param n_obs_periods Integer. The number of observation periods.
#' @param n_subj_per_period Integer. The number of subjects per observation period.
#' @param within_period_cor Numeric. Baseline correlation between observations on different subjects within the same time period.
#' @param cor_decay_rate Numeric. Decay rate applied to the correlation between observations across different time periods.
#' The correlation between observations in different periods is given by within_period_cor * cor_decay_rate^|period difference|.
#'
#' @return A correlation matrix following the exponential decay model with row- (and column-) dimension n_obs_periods*n_subj_per_period.
#'
#' @examples
#' .CreateExponentialDecayCorMat(n_obs_periods = 3,
#'                               n_subj_per_period = 2,
#'                               within_period_cor = 0.2,
#'                               cor_decay_rate = 0.9)

.CreateExponentialDecayCorMat <- function (n_obs_periods,
                              n_subj_per_period,
                              within_period_cor,
                              cor_decay_rate){

  cor_dim <- n_obs_periods * n_subj_per_period

  cor_mat <- matrix(0, nrow = cor_dim,
                     ncol = cor_dim)

  for(i in 1:cor_dim){
    for(j in i:cor_dim){
      if(i == j){
        cor_mat[i,j] <- 1
      } else {
        i_period <- ceiling(i / n_subj_per_period)
        j_period <- ceiling(j / n_subj_per_period)
        period_dif <- abs(i_period - j_period)

        cor_mat[i,j] <- within_period_cor * cor_decay_rate^period_dif
        cor_mat[j,i] <- cor_mat[i,j]
      }
    }
  }

  return(cor_mat)
}

## Cohort Designs
##############################################################################################
#' @title Create Block Exchangeable Correlation Matrix
#'
#' @description This function creates a block-exchangeable correlation matrix using the provided parameters.
#' Assumes inputs were checked in the parent function 'SpecifyCorrelationMatrix'.
#'
#' @param n_obs_periods Number of observation periods
#' @param n_subj_per_period Number of subjects per observation period
#' @param within_period_cor Correlation between observations on different subjects
#' in the same cluster in the same time period
#' @param between_period_cor Correlation between observations on different subjects
#' in the same cluster in different time periods
#' @param within_subject_cor Correlation between observations on the same subject
#' in different time periods
#'
#' @return A block-exchangeable correlation matrix with row- (and column-) dimension n_obs_periods*n_subj_per_period
#'
#' @examples
#' .CreateBlockExchangeableCorMat(n_obs_periods = 2,
#'                               n_subj_per_period = 3,
#'                               within_period_cor = 0.1,
#'                               between_period_cor = 0.2,
#'                               within_subject_cor = 0.3)
.CreateBlockExchangeableCorMat <- function(n_obs_periods,
                                          n_subj_per_period,
                                          within_period_cor,
                                          between_period_cor,
                                          within_subject_cor){

  # Create sub-components of the correlation matrix
  diagonal_mat <- diag(1 - within_period_cor + between_period_cor - within_subject_cor,
                       nrow = n_obs_periods * n_subj_per_period)

  within_subj_component <- (within_subject_cor - between_period_cor) *
    matrix(data = 1, nrow = n_obs_periods, ncol = n_obs_periods) %x%
      diag(1, nrow = n_subj_per_period)

  within_period_component <- (within_period_cor - between_period_cor) *
    diag(1, n_obs_periods)  %x% matrix(data = 1, nrow = n_subj_per_period, ncol = n_subj_per_period)

  between_period_component <- between_period_cor * matrix(data = 1,
                                                          nrow = n_obs_periods * n_subj_per_period,
                                                          ncol = n_obs_periods * n_subj_per_period)

  # Create block exchangeable correlation matrix
  block_exchangeable_cor_mat <-   diagonal_mat + within_period_component + within_subj_component + between_period_component

  return(block_exchangeable_cor_mat)

}
#' @title Create Proportional Decay Correlation Matrix
#'
#' @description This function constructs a correlation matrix based on a proportional decay model. The function
#' assumes that inputs have been validated in the parent function 'SpecifyCorrelationMatrix'.
#'
#' @param n_obs_periods Integer. The number of observation periods.
#' @param n_subj_per_period Integer. The number of subjects per observation period.
#' @param within_period_cor Numeric. Baseline correlation between observations on different subjects within the same time period.
#' @param cor_decay_rate Numeric. Decay rate applied to the correlation between observations across different time periods.
#' The correlation between observations on different subjects in different periods is given by within_period_cor * cor_decay_rate^|period difference|.
#' The correlation between observations on the same subject in different periods is given by cor_decay_rate^|period difference|.
#'
#' @return A correlation matrix following the proportional decay model with row- (and column-) dimension n_obs_periods*n_subj_per_period.
#'
#' @examples
#' .CreateProportionalDecayCorMat(n_obs_periods = 3,
#'                                n_subj_per_period = 2,
#'                                within_period_cor = 0.2,
#'                                cor_decay_rate = 0.9)
.CreateProportionalDecayCorMat <- function (n_obs_periods,
                                           n_subj_per_period,
                                           within_period_cor,
                                           cor_decay_rate){
  cor_dim <- n_obs_periods * n_subj_per_period

  cor_mat <- matrix(0, nrow = cor_dim,
                    ncol = cor_dim)

  for(i in 1:cor_dim){
    for(j in i:cor_dim){
      if(i == j){
        cor_mat[i,j] <- 1
      } else {
        i_period <- ceiling(i / n_subj_per_period)
        j_period <- ceiling(j / n_subj_per_period)
        period_dif <- abs(i_period - j_period)

        same_individual_indicator <- (i %% n_subj_per_period == j %% n_subj_per_period)

        cor_mat[i,j] <- ifelse(same_individual_indicator == TRUE,
                                 cor_decay_rate^period_dif,
                                 within_period_cor * cor_decay_rate^period_dif)

        cor_mat[j,i] <- cor_mat[i,j]
      }
    }
  }

  return(cor_mat)
}
##############################################################################################
## Helper Functions
##
##############################################################################################

#' @title Check inputs for correlation matrix functions
#'
#' @description Internal helper function to check the inputs for functions used to create correlation matrices.
#' Called from the wrapper function \code{SpecifyCorrelationMatrix} rather than individual correlation-specific helper functions.
#'
#' @param n_obs_periods Number of observation periods (positive integer).
#' @param n_subj_per_period Number of subjects per observation period (positive integer).
#' @param design_type_sanitized A string specifying the type of study design, cohort or cross-sectional (character).
#' @param within_period_cor Correlation between observations on different subjects in the same cluster in the same time period. Required for all correlation structures (numeric between 0 and 1 exclusive).
#' @param between_period_cor Correlation between observations on different subjects in the same cluster in different time periods. Required for "block_exchangeable" and "nested_exchangeable" structures (numeric between 0 and 1 exclusive).
#' @param within_subject_cor Correlation between observations on the same subject in different time periods. Required for "block_exchangeable" structure (numeric between 0 and 1 exclusive).
#' @param cor_decay_rate The rate of correlation decay. Required for proportional and exponential decay structures; its exact meaning depends on the correlation structure (numeric between 0 and 1 exclusive).
#'
#' @return NULL if the inputs are correct; otherwise, it raises an error with a descriptive message.
#'
#' @examples
#'
#' # Example usage:
#'

.CheckCorMatInputs <- function(n_obs_periods,
                               n_subj_per_period,
                               design_type_sanitized,
                               within_period_cor = NULL,
                               between_period_cor = NULL,
                               within_subject_cor = NULL,
                               cor_decay_rate = NULL) {

  .CheckCorMatInputTypes(n_obs_periods = n_obs_periods,
                         n_subj_per_period = n_subj_per_period,
                         design_type = design_type_sanitized,
                         within_period_cor = within_period_cor,
                         between_period_cor = between_period_cor,
                         within_subject_cor = within_subject_cor,
                         cor_decay_rate = cor_decay_rate)

  .CheckCorMatInputsNumericVals(n_obs_periods = n_obs_periods,
                                n_subj_per_period = n_subj_per_period,
                                within_period_cor = within_period_cor,
                                between_period_cor = between_period_cor,
                                within_subject_cor = within_subject_cor,
                                cor_decay_rate = cor_decay_rate)

  .CheckCorMatInputsLogic(design_type_sanitized = design_type_sanitized,
                          between_period_cor = between_period_cor,
                          cor_decay_rate = cor_decay_rate)

  return(NULL)
}

#' @title Check logic of correlation matrix inputs
#'
#' @description Internal helper function to check the logical consistency of correlation matrix inputs.
#'
#' @param design_type_sanitized A string specifying the type of study design, cohort or cross-sectional (character).
#' @param between_period_cor Correlation between observations on different subjects in the same cluster in different time periods. Required for "block_exchangeable" and "nested_exchangeable" structures (numeric between 0 and 1 exclusive).
#' @param cor_decay_rate The rate of correlation decay. Required for proportional and exponential decay structures; its exact meaning depends on the correlation structure (numeric between 0 and 1 exclusive).
#'
#' @return NULL if the inputs are logically consistent; otherwise, it raises an error with a descriptive message.
#'
#' @examples
#'
#' # Example usage:
#' .CheckCorMatInputsLogic(design_type_sanitized = "cohort", between_period_cor = 0.3, cor_decay_rate = 0.8)
#'

.CheckCorMatInputsLogic <- function(design_type_sanitized,
                                    between_period_cor = NULL,
                                    cor_decay_rate = NULL){

  if((design_type_sanitized %in% c("cross", "cohort")) == FALSE){
    stop("ERROR: design_type must be either 'cross' or 'cohort'")
  }

  # Check that the inputs only contain parameters for one correlation structure

  if(any(is.null(c(between_period_cor, within_subject_cor)) == FALSE) & !is.null(cor_decay_rate)) {
    stop("ERROR: Parameters for both exchangeable and decay correlation structures were provided. Please provide parameters for only one correlation structure.")
  }
  
  if(all(is.null(c(between_period_cor, within_subject_cor, cor_decay_rate)))) {
    stop("ERROR: No parameters for any correlation structure were provided. Please provide parameters for one correlation structure.")
  }

  # Check that the necessary parameters for each correlation structure type are provided
  if(is.null(cor_decay_rate)){
    if (design_type_sanitized == "cross" && any(is.null(c(within_period_cor, between_period_cor)))){
      stop("ERROR: For the nested exchangeable structure, you must provide within_period_cor and between_period_cor")
    }

    if (design_type_sanitized == "cross" && !is.null(within_subject_cor)){
      stop("ERROR: within_subject_cor must not be speficied in cross-sectional designs")
    }
    
    if (design_type_sanitized == "cohort" && any(is.null(c(within_period_cor, between_period_cor, within_subject_cor)))){
      stop("ERROR: For the block exchangeable structure, you must provide within_period_cor, between_period_cor, and within_subject_cor")
    }
  }
  
  return(NULL)
}

.CheckCorMatInputsNumericVals <- function(n_obs_periods,
                                          n_subj_per_period,
                                          within_period_cor = NULL,
                                          between_period_cor = NULL,
                                          within_subject_cor = NULL,
                                          cor_decay_rate = NULL){

  # Check that correlations are between 0 and 1
  if (!is.null(within_period_cor) && (within_period_cor <= 0 || within_period_cor >= 1)){
    stop("ERROR: within_period_cor must be between 0 and 1 exclusive")
  }
  if (!is.null(between_period_cor) && (between_period_cor <= 0 || between_period_cor >= 1)){
    stop("ERROR: between_period_cor must be between 0 and 1 exclusive")
  }
  if (!is.null(within_subject_cor) && (within_subject_cor <= 0 || within_subject_cor >= 1)){
    stop("ERROR: within_subject_cor must be between 0 and 1 exclusive")
  }
  if (!is.null(cor_decay_rate) && (cor_decay_rate <= 0 || cor_decay_rate >= 1)){
    stop("ERROR: cor_decay_rate must be between 0 and 1 exclusive")
  }

  # Check that n_obs_periods and n_subj_per_period are positive integers
  if (n_obs_periods < 1 || !all.equal(n_obs_periods %% 1, 0)){
    stop("ERROR: n_obs_periods must be a positive integer")
  }
  if (n_subj_per_period < 1 || !all.equal(n_subj_per_period %% 1, 0)){
    stop("ERROR: n_subj_per_period must be a positive integer")
  }

  return(NULL)
}

.CheckCorMatInputTypes <- function(n_obs_periods,
                                  n_subj_per_period,
                                  design_type,
                                  within_period_cor,
                                  between_period_cor = NULL,
                                  within_subject_cor = NULL,
                                  cor_decay_rate = NULL){

  if(is.character(design_type) == FALSE){
    stop("ERROR: design_type must be a character string")
  }

  # Check that the required inputs are numeric
  if(all(is.numeric(c(n_obs_periods, n_subj_per_period, within_period_cor))) == FALSE){
  stop("ERROR: n_obs_periods, n_subj_per_period, and within_period_cor must be numeric")
  }

  # Check that the potentially null correlations are numeric
  if (!is.null(between_period_cor) && !is.numeric(between_period_cor)){
  stop("ERROR: between_period_cor must be numeric")
  }
  if (!is.null(within_subject_cor) && !is.numeric(within_subject_cor)){
  stop("ERROR: within_subject_cor must be numeric")
  }
  if (!is.null(cor_decay_rate) && !is.numeric(cor_decay_rate)){
  stop("ERROR: cor_decay_rate must be numeric")
  }

  return(NULL)
}