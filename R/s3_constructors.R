
if(FALSE){
#' @title Create a `working_cor` (Working Correlation Matrix) object
#' @noRd
#' @description Under development
#'
#' @keywords internal,s3,class
  NewWorkingCorrelationMatrix <- function(correlation_matrix, n_subj_per_period, n_obs_periods, sample_type,
                            within_period_cor, between_period_cor = NULL, within_subject_cor = NULL,
                            cor_decay_rate = NULL, cor_structure_type) {
    if (!is.matrix(correlation_matrix)) {
      stop("correlation_matrix must be a matrix.")
    }

    new_cor_mat <- structure(
      list(
        correlation_matrix = correlation_matrix,
        n_subj_per_period = n_subj_per_period,
        n_obs_periods = n_obs_periods,
        sample_type = sample_type,
        within_period_cor = within_period_cor,
        between_period_cor = between_period_cor,
        within_subject_cor = within_subject_cor,
        cor_decay_rate = cor_decay_rate,
        cor_structure_type = cor_structure_type,
        is_valid = is_valid
      ),
      class = "working_cor_mat"
    )

    return(new_cor_mat)
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
