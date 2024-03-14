library(devtools)

load_all()


## Specify All the relevant parameters

# Design matrices
n_study_periods <- 4
n_clust_trt_seqs <- 3
n_clust_per_seq <- 10
crossover_time_period_vec <- .CreateCrossoverTimePeriodVec(n_study_periods = n_study_periods,
                                                           n_clust_trt_seqs = n_clust_trt_seqs)
n_ind_per_clust <- 1
time_model_type <- "linear"
trt_model_type <- "average"
linear_trt_scale_factor <- 1
mli_study_flag <- FALSE

# Working correlation matrices
sample_type <- "cohort"
within_period_cor <-  0
between_period_cor <-  0
within_subject_cor <-  0
cor_decay_rate <- NULL

# Data Generating Model
time_intercept_param <- 0
time_trend_param <-  0
period_effect_param_vec <-  NULL
clust_trt_param <- .25
ind_trt_param <-  NULL
interaction_param <-  NULL

# Power Calculation
dispersion_scalar <-  1
family <- gaussian()
alpha <- .05

##

design_mat_list <- CreateClusterCompleteDesignMatrixList(n_study_periods = n_study_periods,
                                                         n_clust_trt_seqs = n_clust_trt_seqs,
                                                         n_clust_per_seq = n_clust_per_seq,
                                                         crossover_time_period_vec = crossover_time_period_vec,
                                                         n_ind_per_clust = n_ind_per_clust,
                                                         time_model_type = time_model_type,
                                                         trt_model_type = trt_model_type,
                                                         linear_trt_scale_factor = linear_trt_scale_factor,
                                                         mli_study_flag = mli_study_flag)

n_mean_model_params <- ncol(design_mat_list[[1]])

test_df <- .CalculateDefaultTestDf(n_clusters = n_clust_trt_seqs*n_clust_per_seq,
                                   n_mean_model_params = n_mean_model_params)
contrast_mat <- .CreateDefaultContrastMat(n_time_model_params = 2,
                                          n_trt_model_params = 1)
null_val_vec <- .CreateDefaultNullValVec(mli_study_flag = mli_study_flag)


working_cor_mat_list <- CreateCompleteWorkingCorMatList(
  n_study_periods = n_study_periods,
  n_clusters = n_clust_trt_seqs*n_clust_per_seq,
  n_ind_per_clust = n_ind_per_clust,
  sample_type = "cohort",
  within_period_cor = within_period_cor,
  between_period_cor = between_period_cor,
  within_subject_cor = within_subject_cor,
  cor_decay_rate = cor_decay_rate
)


dgm_param_vec <- SpecifyDataGeneratingModel(n_study_periods = n_study_periods,
                                            mli_study_flag = mli_study_flag,
                                            time_model_type = time_model_type,
                                            time_intercept_param = time_intercept_param,
                                            time_trend_param = time_trend_param,
                                            period_effect_param_vec = period_effect_param_vec,
                                            trt_model_type = trt_model_type,
                                            clust_trt_param = clust_trt_param,
                                            ind_trt_param = ind_trt_param,
                                            interaction_param = interaction_param)

dgm_param_col_vec <- matrix(dgm_param_vec, ncol = 1)

power_list <- .CalcPower(design_mat_list = design_mat_list,
                         working_cor_mat_list = working_cor_mat_list,
                         dgm_param_col_vec = dgm_param_col_vec,
                         dispersion_scalar = dispersion_scalar,
                         family = family,
                         alpha = alpha,
                         test_df = test_df,
                         contrast_mat = contrast_mat,
                         null_val_vec = null_val_vec,
                         power_only_flag = FALSE
)

ex_power <- power_list$power
noncentrality_param <- power_list$noncentrality_param
model_based_var_mat <- round(power_list$model_based_var_mat, 5)
estimated_param_vec <- round(power_list$estimated_param_vec, 5)
bi_list <- power_list$bi_var_matrix_list
vi_list <- power_list$vi_var_matrix_list
di_list <- power_list$di_mean_jacobian_mat_list

library(tidyverse)

n_reps <- 10000
sim_lms <- vector(mode = "list", length = n_reps)
sim_t_tests <- vector(mode = "list", length = n_reps)
sim_p <- vector(mode = "numeric", length = n_reps)
x <- bind_rows(lapply(design_mat_list, as_tibble))
xmat <- as.matrix(x)
y <- xmat %*% dgm_param_col_vec
n_obs <- nrow(xmat)

trt_group <- x$TrtClust
set.seed(42)
for(i in 1:length(sim_res)){
  obs <- y + rnorm(n = n_obs, mean = 0, sd = 1)
  rep_lm <- lm(obs ~ -1 + xmat, x = FALSE, y = FALSE)
  sim_res[[i]] <- rep_lm
  sim_p[[i]] <- summary(rep_lm)$coefficients[3,4]
  sim_t_tests[[i]] <- t.test(obs ~ trt_group)
}
