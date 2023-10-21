
##############################################################################################
## Data Generating Model Specification Function
##
##############################################################################################

#' @title Specify Data Generating Model
#'
#' @description
#' Wrapper function for specifying the data generating model for a given study.
#' This function allows for the specification of both calendar time and treatment effects.
#'
#' @param n_study_periods Integer. The number of study periods. Must be at least 2.
#' @param mli_study_flag Logical. Flag to indicate if the study is investigating a multi-level intervention.
#' @param time_model_type Character. The type of time effect model to be used. Valid values are "linear" and "categorical".
#' @param time_intercept_param Numeric. The intercept parameter for the time model. Required if `time_model_type` is "linear".
#' @param time_trend_param Numeric. The trend parameter for the time model. Required if `time_model_type` is "linear".
#' @param period_effect_param_vec Numeric vector. The period effect parameters for the time model. Required if `time_model_type` is "categorical".
#' @param trt_model_type Character. The type of treatment effect model to be used. Valid values are "linear" and "average".
#' @param clust_trt_param Numeric. The treatment effect parameter at the cluster level.
#' @param ind_trt_param Numeric. The treatment effect parameter at the individual level. Only required if `mli_study_flag` is TRUE.
#' @param interaction_param Numeric. The interaction parameter for the treatment effect. Only required if `mli_study_flag` is TRUE.
#'
#' @return A named numeric vector of model parameters.
#'
#'
#' @examples
#' \dontrun{
#' SpecifyDataGeneratingModel(n_study_periods = 5,
#'                            mli_study_flag = FALSE,
#'                            time_model_type = "linear",
#'                            time_intercept_param = 0.5,
#'                            time_trend_param = 0.1,
#'                            trt_model_type = "linear",
#'                            clust_trt_param = 0.3)
#' }
#'
#' @export
SpecifyDataGeneratingModel <- function (n_study_periods = NULL,
                                        mli_study_flag = FALSE,
                                        time_model_type = "linear",
                                        time_intercept_param = NULL,
                                        time_trend_param = NULL,
                                        period_effect_param_vec = NULL,
                                        trt_model_type = "linear",
                                        clust_trt_param,
                                        ind_trt_param = NULL,
                                        interaction_param = NULL){
  # Perform input checks
  .CheckDGMInputs(n_study_periods = n_study_periods,
                  mli_study_flag = mli_study_flag,
                  time_model_type = time_model_type,
                  time_intercept_param = time_intercept_param,
                  time_trend_param = time_trend_param,
                  period_effect_param_vec = period_effect_param_vec,
                  trt_model_type = trt_model_type,
                  clust_trt_param = clust_trt_param,
                  ind_trt_param = ind_trt_param,
                  interaction_param = interaction_param)

  # Specify model parameters in the order that corresponds to the design matrix
  named_time_params <- switch(time_model_type,
                           "linear" = .SpecifyLinearTimeEffect(time_intercept_param,
                                                               time_trend_param),
                           "categorical" = .SpecifyCategoricalTimeEffect(n_study_periods,
                                                                         period_effect_param_vec))

  named_trt_params <- switch(trt_model_type,
                             "linear" = .SpecifyAvgOrLinearTrtEffect(mli_study_flag,
                                                                    clust_trt_param,
                                                                    ind_trt_param,
                                                                    interaction_param),
                             "average" = .SpecifyAvgOrLinearTrtEffect(mli_study_flag,
                                                                     clust_trt_param,
                                                                     ind_trt_param,
                                                                     interaction_param))
  # Provide informative names for the parameters
  names(named_time_params) <- .NameTimeEffects(time_model_type, n_study_periods = n_study_periods)
  names(named_trt_params) <- .NameTrtEffects(mli_study_flag = mli_study_flag)

  dgm_param_vec <- c(named_time_params, named_trt_params, use.names = TRUE)

  .CheckDGMOutput(dgm_param_vec)

  return(dgm_param_vec)
}

##############################################################################################
## Individual Effect-type Specification Functions
##
##############################################################################################
#' @title Specify Linear Time Effect Parameters
#'
#' @description
#' Internal function to specify the parameters for a linear time effect in the data generating model.
#'
#' @param time_intercept_param Numeric. The intercept parameter for the time model.
#' @param time_trend_param Numeric. The trend parameter for the time model.
#'
#' @return A named numeric vector with the time effect parameters.
#'
#' @keywords internal
.SpecifyLinearTimeEffect <- function(time_intercept_param,
                                     time_trend_param){
  checkmate::expect_number(time_intercept_param, null.ok = FALSE)
  checkmate::expect_number(time_trend_param, null.ok = FALSE)

  time_effect_param_vec <- c(time_intercept_param, time_trend_param)

  return(time_effect_param_vec)
}

#' @title Specify Categorical Time Effect Parameters
#'
#' @description
#' Internal function to specify the parameters for a categorical time effect in the data generating model.
#'
#' @param n_study_periods Numeric. The number of study periods, must be at least 2.
#' @param period_effect_param_vec Numeric vector. The period effect parameters for the time model.
#'
#' @return A named numeric vector with the time effect parameters.
#'
#' @keywords internal
.SpecifyCategoricalTimeEffect <- function(n_study_periods, period_effect_param_vec){
  checkmate::expect_number(n_study_periods, lower = 2, null.ok = FALSE)
  checkmate::expect_numeric(period_effect_param_vec, any.missing = FALSE, null.ok = FALSE)
  testthat::expect_length(period_effect_param_vec, n = n_study_periods)

  return(period_effect_param_vec)
}

#' @title Specify Average or Linear Treatment Effect Parameters
#'
#' @description
#' Internal function to specify the parameters for average or linear treatment effects in the data generating model.
#'
#' @param mli_study_flag Logical. Flag to indicate if the study is investigating a multi-level intervention.
#' @param clust_trt_param Numeric. The treatment effect parameter at the cluster level.
#' @param ind_trt_param Numeric. The treatment effect parameter at the individual level. Optional if `mli_study_flag` is FALSE.
#' @param interaction_param Numeric. The interaction parameter for the treatment effect. Optional if `mli_study_flag` is FALSE.
#'
#' @return A named numeric vector with the treatment effect parameters.
#'
#' @keywords internal
.SpecifyAvgOrLinearTrtEffect <- function(mli_study_flag,
                                   clust_trt_param,
                                   ind_trt_param,
                                   interaction_param){
  checkmate::expect_flag(mli_study_flag, null.ok = FALSE)
  checkmate::expect_number(clust_trt_param, null.ok = FALSE)
  checkmate::expect_number(ind_trt_param, null.ok = (mli_study_flag == FALSE))
  checkmate::expect_number(interaction_param, null.ok = (mli_study_flag == FALSE))

  trt_param_vec <- c(clust_trt_param, ind_trt_param, interaction_param)
  trt_param_vec <- trt_param_vec[is.null(trt_param_vec) == FALSE]

  return(trt_param_vec)
}




##############################################################################################
## Input Checks
##
##############################################################################################
  #' @title Check Inputs for Data Generating Model Specification
  #'
  #' @description
  #' Internal function to validate the inputs for the `SpecifyDataGeneratingModel` function.
  #' For details on the parameters, see the documentation for `SpecifyDataGeneratingModel`.
  #'
  #' @return NULL. Function returns NULL if all checks pass.
  #'
  #' @seealso \code{\link{SpecifyDataGeneratingModel}}
  #' @keywords internal
.CheckDGMInputs <- function(n_study_periods,
                            mli_study_flag,
                            time_model_type,
                            time_intercept_param,
                            time_trend_param,
                            period_effect_param_vec,
                            trt_model_type,
                            clust_trt_param,
                            ind_trt_param,
                            interaction_param) {
  # Should be at least two study periods
  checkmate::expect_number(n_study_periods, lower = 2)

  # Is the study investigating a multi-level intervention?
  checkmate::expect_flag(mli_study_flag)

  # Current time effect models are 'linear' and 'categorical'
  checkmate::expect_choice(time_model_type, choices = c("linear", "categorical"))


  checkmate::expect_number(time_intercept_param, null.ok = (time_model_type != "linear"))
  checkmate::expect_number(time_trend_param, null.ok = (time_model_type != "linear"))

  # Numeric vector
  checkmate::expect_numeric(period_effect_param_vec, any.missing = FALSE,
                            null.ok = (time_model_type != "categorical"))

  # Should be either 'linear' or 'average'
  checkmate::expect_choice(trt_model_type, choices = c("linear", "average"))

  # Numeric
  checkmate::expect_number(clust_trt_param)

  # Numeric, allowed to be NULL if mli_study_flag is FALSE
  checkmate::expect_number(ind_trt_param, null.ok = (mli_study_flag == FALSE))

  # Numeric, allowed to be NULL if mli_study_flag is FALSE
  checkmate::expect_number(interaction_param, null.ok = (mli_study_flag == FALSE))

  return(NULL)
}


#' @title Check Inputs for Data Generating Model Specification
#'
#' @description
#' Internal function to validate the output for the `SpecifyDataGeneratingModel` function.
#' Checks that the output is a named numeric vector.
#'
#' @return NULL. Function returns NULL if all checks pass.
#'
#' @seealso \code{\link{SpecifyDataGeneratingModel}}
#' @keywords internal
.CheckDGMOutput <- function(dgm_param_vec){
  checkmate::expect_numeric(dgm_param_vec, any.missing = FALSE)
  checkmate::expect_names(names(dgm_param_vec))
  return(NULL)
}