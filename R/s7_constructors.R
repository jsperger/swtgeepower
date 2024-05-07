#### Constructors for S7 Classes

                                        # Base Classes

#' Study Design Class
#'
#' Only incorporates the elements of the study. Does not include the data generating model.
#' @noRd
study_design <- new_class("study_design",
                          properties = list(
                              clust_sample_size = class_numeric,
                              ind_per_clust = class_numeric,
                              n_time_periods = class_numeric,
                              sample_type = class_character,
                              ind_sample_size = new_property(
                                getter = function(self) {
                                  n_ind <- ifelse(length(self@ind_per_clust) == 1,
                                          self@clust_sample_size * self@ind_per_clust,
                                          sum(self@ind_per_clust))
                                return(n_ind)
                                }
                              ),
                              multilevel_study = class_logical
                          ))





#' @title Correlation Structure Class
#' @noRd
#'
#' @keywords internal, s7, class
cor_struct <- new_class("correlation_structure",
                     properties = list(
                         sample_type = class_character,
                         structure_type = class_character,
                         within_period_cor = class_numeric
))

data_gen_model <- new_class("data_generating_model", properties = list(
  time_model_type = class_character,
  trt_model_type = class_character,
  time_model_params = class_numeric,
  trt_model_params = class_numeric,
  family = class_character,
  dispersion_scalar = class_numeric
))

power_summary <- new_class("power_summary", properties = list(
  study_design = class_list,
  data_gen_model = class_list,
  cor_struct = class_list,
  power = class_numeric,
  alpha = class_numeric,
  test_df = class_numeric,
  contrast_mat = class_matrix,
  null_val_vec = class_numeric
))

# Child Classes

#TODO: Consider whether a treatment assignment matrix should be included in the study design class
# or if the child classes should have their simplified versions of the treatment assignment matrix.
# E.g. a parallel group design only needs to know the number of treatment groups and the number of
# subjects in each group. A crossover design needs to know the number of treatment groups and the
# number of subjects in each group, but also the crossover time period.

parallel_group_design <- new_class("parallel_group_design",
                                   parent = "study_design",
                                  properties = list(
                                      baseline_period = class_logical,
                                      final_all_trt_period = class_logical
                                  ))

# This assumes only two treatment groups and a single crossover time period
crossover_design <- new_class("crossover_design",
                                    parent = "study_design",
                              properties = list(
                                  crossover_time_period = class_numeric
                              ))


stepped_wedge_design <- new_class("stepped_wedge_design", parent = "study_design",
                                  properties = list(
                                      n_seq = class_numeric,
                                      open_cohort = class_logical,
                                      n_enroll_periods = class_numeric,
                                      n_rand_periods = class_numeric,
                                      n_obs_periods = class_numeric,
                                      n_per_clust_period = class_numeric
                                  ))

autoregressive_cor_struct <- new_class("autoregressive_correlation_structure",
                                       parent = "correlation_structure",
                                       properties = list(
                                           cor_decay_rate = class_numeric
                                       ))

fixed_cor_struct <- new_class("fixed_correlation_structure",
                              parent = "correlation_structure",
                              properties = list(
                                  between_period_cor = class_numeric,
                                  within_subject_cor = class_numeric
                              ))