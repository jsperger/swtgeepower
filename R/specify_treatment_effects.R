
##############################################################################################
## Data Generating Model Specification Function
##
##############################################################################################

#' @title Specify Data Generating Model
#'
#' @description
#' Wrapper function for specifying the data generating model including calendar time and treatment effects.
#' @return a numeric vector of model parameters
SpecifyDataGeneratingModel <- function (time_model_type,
                                        time_model_params,
                                        treatment_model_type,
                                        treatment_model_params){
  stopifnot(time_model_type %in% c("linear", "per-period"))
  stopifnot(treatment_model_type %in% c("linear", "average"))

  param_vec <- c(time_model_params, treatment_model_params)
}

##############################################################################################
## Individual Effect-type Specification Functions
##
##############################################################################################

SpecifyCalendarTimeEffects <- function(time_effects_model = "linear",
                                       time_effect_parameters = NULL){

}

SpecifyTreatmentEffects <- function(treatment_effects_model = "linear",
                                    treatment_effect_parameters = NULL){

}



##############################################################################################
## Input Checks
##
##############################################################################################