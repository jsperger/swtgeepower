
if(FALSE){
#' @title Create a `working_cor` (Working Correlation Matrix) object
#' @noRd
#' @description Under development
#'
#' @keywords internal,s3,class
  NewWorkingCorrelationMatrix <- function(){
    new_instance <- structure(list(), class = "working_correlation_matrix")
  }

  #' @title Create a `swtdesign` object
  #' @noRd
  #' @description Under development
  #'
  #' @keywords internal,s3,class
  NewSwtdesign <- function(n_time_periods,
                           n_obs_periods_per_clust,
                           n_clust,
                           n_seq,
                           n_subj_per_clust,
                           n_enroll_periods,
                           n_rand_periods,
                           n_obs_periods,
                           n_per_clust_period){
    new_instance <- structure(list(), class = "swtdesign")
  }

}
