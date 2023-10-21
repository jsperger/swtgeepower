
##############################################################################################
## Data Generating Model Specification Function
##
##############################################################################################

#' @title Specify Data Generating Model
#'
#' @description
#' Wrapper function for specifying the data generating model including calendar time and treatment effects.
#' @return a numeric vector of model parameters
SpecifyDataGeneratingModel <- function (n_study_periods = NULL,
                                        mli_study_flag = FALSE,
                                        time_model_type = "linear",
                                        time_intercept_param = NULL,
                                        time_trend_param = NULL,
                                        period_effect_param_vec = NULL,
                                        treatment_model_type = "linear",
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
                  treatment_model_type = treatment_model_type,
                  clust_trt_param = clust_trt_param,
                  ind_trt_param = ind_trt_param,
                  interaction_param = interaction_param)

  # Specify model parameters in the order that corresponds to the design matrix
  named_time_params <- switch(time_model_type,
                           "linear" = .SpecifyLinearTimeEffect(time_intercept_param,
                                                               time_trend_param),
                           "categorical" = .SpecifyCategoricalTimeEffect(n_study_periods,
                                                                         period_effect_param_vec))

  named_trt_params <- switch(treatment_model_type,
                             "linear" = .SpecifyAvgOrLinearTrtEffect(mli_study_flag,
                                                                    clust_trt_param,
                                                                    ind_trt_param,
                                                                    interaction_param),
                             "average" = .SpecifyAvgOrLinearTrtEffect(mli_study_flag,
                                                                     clust_trt_param,
                                                                     ind_trt_param,
                                                                     interaction_param))
  # Provide informative names for the parameters
  names(named_time_params) <- .NameTimeEffects(time_model_type)
  names(named_trt_params) <- .NameTrtEffects(mli_study_flag = mli_study_flag)

  dgm_param_vec <- c(named_time_params, named_trt_params, use.names = TRUE)

  .CheckDGMOutput(dgm_param_vec)

  return(dgm_param_vec)
}

##############################################################################################
## Individual Effect-type Specification Functions
##
##############################################################################################

.SpecifyLinearTimeEffect <- function(time_intercept_param,
                                     time_trend_param){
  checkmate::expect_number(time_intercept_param, null.ok = FALSE)
  checkmate::expect_number(time_trend_param, null.ok = FALSE)

  time_effect_param_vec <- c(time_intercept_param, time_trend_param)

  return(time_effect_param_vec)
}

.SpecifyCategoricalTimeEffect <- function(n_study_periods, period_effect_param_vec){
  checkmate::expect_number(n_study_periods, lower = 2, null.ok = FALSE)
  checkmate::expect_numeric(period_effect_param_vec, any.missing = FALSE, null.ok = FALSE)
  testthat::expect_length(period_effect_param_vec, n = n_study_periods)

  return(period_effect_param_vec)
}


.SpecifyAvgOrLinearTrtEffect <- function(mli_study_flag,
                                   clust_trt_param,
                                   ind_trt_param,
                                   interaction_param){
  checkmate::expect_flag(mli_study_flag, null.ok = FALSE)
  checkmate::expect_number(clust_trt_param, null.ok = FALSE)
  checkmate::expect_number(ind_trt_param, null.ok = (mli_study_flag == FALSE))
  checkmate::expect_number(interaction_param, null.ok = (mli_study_flag == FALSE))

  trt_param_vec <- c(clust_trt_param, ind_trt_param, interaction_param)
  trt_param_vec <- trt_param_vec[is.null(trt_param_vec) == FALSE)]

  return(trt_param_vec)
  }
}



##############################################################################################
## Input Checks
##
##############################################################################################

.CheckDGMInputs <- function(n_study_periods,
                            mli_study_flag,
                            time_model_type,
                            time_intercept_param,
                            time_trend_param,
                            period_effect_param_vec,
                            treatment_model_type,
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
  checkmate::expect_choice(treatment_model_type, choices = c("linear", "average"))

  # Numeric
  checkmate::expect_number(clust_trt_param)

  # Numeric, allowed to be NULL if mli_study_flag is FALSE
  checkmate::expect_number(ind_trt_param, null.ok = (mli_study_flag == FALSE))

  # Numeric, allowed to be NULL if mli_study_flag is FALSE
  checkmate::expect_number(interaction_param, null.ok = (mli_study_flag == FALSE))

  return(NULL)
}

#' Check output of SpecifyDataGeneratingModel
.CheckDGMOutput <- function(dgm_param_vec){
  checkmate::expect_numeric(dgm_param_vec, any.missing = FALSE)
  checkmate::expect_names(names(dgm_param_vec))
  return(NULL)
}