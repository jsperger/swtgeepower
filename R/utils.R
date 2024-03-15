#' Description: Utility functions for the project

# Family wrangling
###############################################################################
#' @title Input processing for family objects
#' @description
#' Ensures that the input `family` is a family object (list). If it is not,
#' `.WrangleFamily` will attempt to convert it to a family object.
#'
#' @param family A family object (list), string, or function.
#' @return A family object (list)
.WrangleFamily <- function (family){

  # String input
  if(is.character(family) == TRUE){
    if(family %in% c("quasi", "quasibinomial", "quasipoisson")){
      stop("quasi families are not supported. If you are only looking for overdispersion,
       the package supports a separate overdispersion parameter")
    }
    checkmate::expect_choice(family, c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian"))

    family <- stats::family(family)
  }
  # Function input
  if(is.function(family) == TRUE){
        family <- family()
    }

  stopifnot(is(family) == "family")
  # General error checking
  return(family)
}
###############################################################################
# Utility Functions
#
###############################################################################

# Default Names
###############################################################################
#' @title Defaults - Default Names for the Time Effects Parameters
#'
#' @description
#' Obtains the names for time effect parameters based on the specified time model type
#' and the number of study periods. This function is crucial for naming conventions
#' in model parameters and the construction of the design matrix.
#'
#' @param time_model_type A character string specifying the type of time model.
#' Can be "linear" or "categorical".
#' @param n_study_periods An integer indicating the number of study periods.
#'
#' @return
#' A character vector containing the names of the time effects parameters,
#' including an "Intercept" and depending on the model type, "Period" for linear
#' models or "Period2", "Period3", ..., for categorical models.
#'
#' @examples
#' .NameTimeEffects("linear", 4)
#' .NameTimeEffects("categorical", 3)
.NameTimeEffects <- function(time_model_type, n_study_periods){
  time_effect_names <- switch(time_model_type,
                              linear = c("Intercept", "Period"),
                              categorical = c("Intercept", paste0("Period", 2:n_study_periods))
  )

  return(time_effect_names)
}

#' @title Defaults - Default Names for the Treatment Effects Parameters
#'
#' @description
#' Determines the names for treatment effect parameters based on whether the study
#' is a multi-level intervention study. This function aids in maintaining consistent
#' naming conventions across treatment parameters.
#'
#' @param mli_study_flag A logical flag indicating if the study is a multi-level
#' intervention study (`TRUE`) or not (`FALSE`).
#'
#' @return
#' A character vector containing the names of the treatment effects parameters.
#' In the case of a multi-level intervention study, this includes "TrtClust",
#' "TrtInd", and "TrtClustInd". Otherwise, only "TrtClust" is returned.
#'
#' @examples
#' .NameTrtEffects(TRUE)
#' .NameTrtEffects(FALSE)
.NameTrtEffects <- function(mli_study_flag){
  if (mli_study_flag == TRUE) return(c("TrtClust", "TrtInd", "TrtClustInd"))

  return("TrtClust")
}

# Default Parameter Dimensions
###############################################################################

#' @title Defaults - Default Dimension for Time Parameters
#'
#' @description
#' Calculates the dimension (number of parameters) for time effects based on the
#' model type and the number of study periods. This function is essential for
#' dimensionality management in models incorporating time effects.
#'
#' @param time_model_type A character string indicating the type of time model
#' being used. Supported values are "linear" and "categorical".
#' @param n_study_periods An integer denoting the number of study periods.
#'
#' @return
#' An integer representing the dimension of the time parameters. This includes
#' both the intercept and the effects (linear or categorical) associated with
#' the time model.
#'
#' @examples
#' .DefaultDimTimeParams("linear", 4)
#' .DefaultDimTimeParams("categorical", 3)
.DefaultDimTimeParams <- function(time_model_type, n_study_periods){
  # Note: this includes the intercept

  dim_time_params <- switch(time_model_type,
                            linear = 2,
                            categorical = n_study_periods,
                            quadratic = stop("Not implemented yet [will be 3]"),
  )

  return(as.integer(dim_time_params))
}


#' @title Defaults - Default Dimension for Treatment Parameters
#'
#' @description
#' Computes the dimension (number of parameters) for treatment effects, considering
#' the study's design (single-level or multi-level intervention) and the treatment
#' model type. It is vital for accurately sizing treatment effect vectors in models.
#'
#' @param mli_study_flag A logical flag indicating if the study involves a
#' multi-level intervention (`TRUE`) or not (`FALSE`).
#' @param trt_model_type A character string specifying the type of treatment model.
#' Supported values are "average" and "linear".
#' @return
#' An integer representing the dimension of the treatment parameters. If the study
#' is a multi-level intervention, this value is adjusted accordingly.
#'
#' @examples
#' .DefaultDimTrtParams(FALSE, "average")
#' .DefaultDimTrtParams(TRUE, "linear")
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

#' @title Defaults - Create Vector of Default Null Hypothesis Values
#'
#' @description
#' Generates a zero vector to use as the default values for hypothesis testing.
#' The dimension of the zero vector depends on whether the study is a
#' multi-level intervention study.
#'
#' @param mli_study_flag A logical flag indicating if the study is a multi-level
#' intervention study (`TRUE`) or not (`FALSE`).
#'
#' @return
#' A column matrix (vector) of zeros with the number of rows determined by the
#' study's intervention level: 5 rows if the study is a multi-level intervention,
#' and 3 rows otherwise.
#'
#' @examples
#' .CreateDefaultNullValVec(TRUE)
#' .CreateDefaultNullValVec(FALSE)
.CreateDefaultNullValVec <- function(mli_study_flag){
  if(mli_study_flag == TRUE) contrast_dim <- 5
  if(mli_study_flag == FALSE) contrast_dim <- 3

  return(matrix(rep(0, times = contrast_dim), ncol = 1))
}


#' @title Defaults - Create a Default Contrast Matrix
#'
#' @description
#' Constructs a default contrast matrix for testing treatment effect parameters.
#' The default contrast matrix is a block diagonal matrix with a block of zeroes
#' for the time-related parameters, and an identity matrix for the block
#' corresponding to the treatment effect parameters.
#'
#' @param n_time_model_params An integer indicating the number of parameters in
#' the time model.
#' @param n_trt_model_params An integer indicating the number of parameters in
#' the treatment model.
#'
#' @return
#' A matrix where the first set of columns (corresponding to the time model
#' parameters) are filled with zeros, and the subsequent columns form an identity
#' matrix based on the number of treatment model parameters.
#' @examples
#' .CreateDefaultContrastMat(2, 3)
#' .CreateDefaultContrastMat(3, 1)
.CreateDefaultContrastMat <- function(n_time_model_params,
                                      n_trt_model_params) {
        contrast_row_dim <- n_trt_model_params

        contrast_mat <- cbind(
                matrix(0, nrow = contrast_row_dim, ncol = n_time_model_params),
                diag(1, nrow = n_trt_model_params)
        )

        return(contrast_mat)
}

# Crossover Time Helpers
###############################################################################

#' @title Utility - Create Vector of Crossover Time Guesses
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

#' @title Utility - Calculate Test Degrees of Freedom
#'
#' @description
#' Calculates the default degrees of freedom for statistical tests of the
#' treatment effects based on the number of clusters and the number of parameters in the mean model.
#'
#' @param n_clusters An integer representing the total number of clusters in the
#' study.
#' @param n_mean_model_params An integer representing the number of parameters
#' in the mean model.
#'
#' @return
#' An integer representing the degrees of freedom for the test, calculated as
#' the number of clusters minus the number of mean model parameters.
#'
#' @examples
#' .CalculateDefaultTestDf(10, 5)
#' .CalculateDefaultTestDf(20, 3)
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
