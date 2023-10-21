#' Description: Utility functions for the project

#' @title Names for the Time Effects Parameters
#'
#' @description Gets its own function because it's used in the model parameter
#' and design matrix functions.
.NameTimeEffects <- function(time_model_type){
  time_effect_names <- switch(time_model_type,
                              linear = c("Intercept", "LinearTimeTrend"),
                              categorical = c("Intercept", paste0("Period", 2:n_study_time_periods))
  )

  return(time_effect_names)
}

.NameTrtEffects <- function(mli_study_flag){
  trt_effect_names <- ifelse(mli_study_flag == TRUE,
                             yes = c("ClustTrtEffect", "IndTrtEffect", "InteractionEffect"),
                             no = "ClustTrtEffect")

  return(trt_effect_names)
}

}

.CheckFrechetBound <- function(){
  stop("Not implemented yet")
}