##############################################################################################
## GEE Power Calculation Wrapper
##
##############################################################################################


#' @title Perform the GEE-based power calculation for stepped wedge trials
#' and multilevel cluster randomized trials. Wrapper
#'
#' @description
#'

CalcGEEPower <- function(design_mat_list = NULL,
                         param_vec = NULL,
                         n_clust_per_seq = NULL,
                         n_ind_per_clust = NULL,
                         clust_trt_effect = NULL,
                         ind_trt_effect = NULL,
                         interaction_effect = NULL,
                         period_effect_vec = NULL,
                         incr_effect_scalar = NULL,
                         var_fun,
                         cor_mat_list = NULL,
                         incidence_mat_list = NULL,
dist = "normal",
                         link = "identity",
dispersion_scalar_scalar = NULL,
                         alpha = .05,
test_df,
                         contrast_mat = NULL){

}


##############################################################################################
## GEE Power Calculation Internal Functions
##
##############################################################################################

#' @title Internal Power Calculation Function
#'
#' @description
#' Perform the GEE-based power calculation for stepped wedge trials and
#' multilevel cluster randomized trials. This function is the internal function
#' for performing the power calculation. It is not exported.
#'
#' @param design_mat_list A list of design matrices, one for each cluster.
#' @param incidence_mat_list A list of incidence matrices, one for each cluster.
#' @param cor_mat_list A list of correlation matrices, one for each cluster.
#' @param trt_param_col_vec A vector of treatment parameters as
#' a matrix with one column.
#' @param dispersion_scalar The dispersion parameter for the GEE model.
#' @param link The link function for the GEE model.
#' @param dist The distribution for the GEE model.
#' @param alpha The desired type I error rate.
#' @param test_df The degrees of freedom for the test statistic.
#' @param power_only_flag A logical flag indicating whether to return only the
#' calculated power. Used internally for sample size calculations
#'
#' @return A list with the following components:
#' \item{power}{The power of the test.}
#' \item{sample_size}{The sample size required to achieve the desired power.}
#' \item{test_stat}{The test statistic.}
#' \item{p_value}{The p-value for the test statistic.}
#' \item{test_df}{The degrees of freedom for the test statistic.}
#'
#'
.CalcPower <- function(design_mat_list,
                       working_cor_mat_list,
                       incidence_mat_list,
                       cor_mat_list,
                       trt_param_col_vec,
                       dispersion_scalar,
                       var_fun,
                       link,
                       response_type,
                       alpha,
                       test_df,
                       contrast_mat,
                       null_val_vec,
                       power_only_flag = FALSE){

  n_clust <- length(design_mat_list)

  # Create placeholder lists
  mui_col_vec_list <- vector(mode = "list", length = n_clust)
  di_mean_jacobian_mat_list <- vector(mode = "list", length = n_clust)
  bi_var_matrix_list <- vector(mode = "list", length = n_clust)
  vi_var_matrix_list <- vector(mode = "list", length = n_clust)

  # Calculate estimating equation components for each cluster
  for(i in 1:n_clust){
    clust_mu_vec <- .CreateMuiVector(design_mat = design_mat_list[[i]],
                                              trt_param_col_vec = trt_param_col_vec,
                                              response_type = response_type,
                                              link = link)

    clust_mean_jacobian_mat_list <- .CreateMeanJacobianMatDi(mu_col_vec = clust_mu_vec,
                                                               link = link)

    clust_bi_var_matrix <- .CreateVarParamMatrixBi(mu_col_vec = clust_mu_vec,
                                                   dispersion_scalar = dispersion_scalar,
                                                   var_fun = var_fun)

    clust_vi_var_matrix <- .CalcVarMatrixVi(var_mat_bi = clust_bi_var_matrix,
                                            working_cor_mat_ri = working_cor_mat_list[[i]])

    # Store the results in the placeholder lists
    mui_col_vec_list[[i]] <- clust_mu_vec
    di_mean_jacobian_mat_list[[i]] <- clust_mean_jacobian_mat_list
    bi_var_matrix_list[[i]] <- clust_bi_var_matrix
    vi_var_matrix_list[[i]] <- clust_vi_var_matrix
  }

# Calculate the model-based variance matrix

  model_based_var_mat <- .CalcModelBasedVariance(
    mean_jacobian_mat_list = di_mean_jacobian_mat_list,
    var_mat_vi_list = vi_var_matrix_list)

  # Calculate the estimated parameter vector
  estimated_param_vec <- .CalcEstimatedParamVec(
    design_mat_list = design_mat_list,
    mean_jacobian_mat_list = di_mean_jacobian_mat_list,
    var_mat_vi_list = vi_var_matrix_list,
    mui_vec_list = mui_col_vec_list,
    invlink_fun = invlink_fun)

  # Calculate the noncentrality parameter for the Wald test
  noncentrality_param <- .CalcWaldNcp(n_clust = n_clust,
                                      contrast_mat = contrast_mat,
                                      null_val_vec = null_val_vec,
                                      estimated_param_vec = estimated_param_vec,
                                      model_based_var_mat = model_based_var_mat)

    # Calculate the critical value for the Wald test
    null_crit_val <- .CalcWaldCritVal(alpha = alpha,
                                 test_df = test_df)

    # Calculate the power of the test
    power <- 1 - pchisq(q = crit_val,
                        df = test_df,
                        ncp = noncentrality_param)

  if(power_only_flag == TRUE) return(power)

  power_calc_results <- list(power = power,
                             cluster_sample_size = n_clust,
                             test_stat = noncentrality_param,
                             alpha = alpha,
                             test_df = test_df)

  return(power_calc_results)
}

.CreateMuiVector <- function (design_mat,
                             trt_param_vec,
                             response_type,
                             link){

  # Expected response on the scale of the link function
  linear_pred_vec <- design_mat %*% trt_param_vec

  mu_vec <- switch(link,
                   identity = linear_pred_vec,
                   logit = exp(linear_pred_vec)/(1 + exp(linear_pred_vec)),
                   log = exp(linear_pred_vec))

  if(response_type == "count" & min(mu_vec) <= 0 ){
    stop("The mean of the response variable must be positive")
  }

  mu_col_vec <- matrix(data = mu_vec, ncol = 1)

  return(mu_col_vec)
}

.CreateMeanJacobianMatDi <- function(mu_col_vec, link){

  mean_jacobian_mat_di <- switch(link,
                                 identity = diag(1, nrow = mu_col_vec),
                                 logit = diag(mu_col_vec) %*% diag(1 - mu_col_vec),
                                 log = diag(mu_col_vec))

  return(mean_jacobian_mat_di)
}

#' @title Create the variance matrix for the GEE model
#'
#' @description
#' Create the variance matrix for the GEE model. This function is not exported.
#' $$B_i = \phi \diag(v(\mu_{i11}), \ldots, v(\mu_{IJK}) )$$
#' where $\phi$ is the dispersion parameter and $v(\mu_{ijk})$ is the variance
.CreateVarParamMatrixBi <- function(mu_col_vec,
                                    dispersion_scalar,
                                    var_fun){

  #concatenate matrix to vector
  var_vec <- var_fun(c(mu_col_vec))

  var_mat <- diag(var_vec)

  var_mat_bi <- dispersion_scalar * var_mat

  return(var_mat_bi)
}

.CalcVarMatrixVi <- function(var_mat_bi,
                             working_cor_mat_ri){

  root_bi <- sqrt(var_mat_bi)

  var_mat_vi <- root_bi %*% working_cor_mat_ri %*% root_bi

  return(var_mat_vi)
}

.CalcModelBasedVariance <- function(mean_jacobian_mat_list,
                                    var_mat_vi_list){

  param_dim <- ncol(mean_jacobian_mat_list[[1]])

  vinv_list <- lapply(var_mat_vi_list, solve)

  model_based_var_mat <- matrix(0, nrow = param_dim, ncol = param_dim)

  for(i in 1:length(mean_jacobian_mat_list)){
    model_based_var_mat <- model_based_var_mat +
      t(mean_jacobian_mat_list[[i]]) %*%
      vinv_list[[i]] %*%
      mean_jacobian_mat_list[[i]]
  }

  return(model_based_var_mat)
}

.CalcEstimatedParamVec <- function(design_mat_list,
                                   mean_jacobian_mat_list,
                                   var_mat_vi_list,
                                   mui_vec_list,
                                   invlink_fun){
  n_clust <- length(design_mat_list)

  shared_term_list <- lapply(1:n_clust,
                             function(i){
                               t(design_mat_list[[i]]) %*%
                                 mean_jacobian_mat_list[[i]] %*%
                                 var_mat_vi_list[[i]] %*%
                                 t(mean_jacobian_mat_list[[i]])
                             })

  #TODO: Come up with a more descriptive name for these variables
left_part_list <- vector(mode = "list", length = n_clust)
  right_part_list <- vector(mode = "list", length = n_clust)
  for(i in 1:n_clust){
    left_part_list[[i]] <- shared_term_list[[i]] %*% design_mat_list[[i]]
    right_part_list[[i]] <- shared_term_list[[i]] %*% invlink_fun(mui_vec_list[[i]])
  }

  left_part_term <- solve(Reduce("+", left_part_list))

  right_part_term <- Reduce("+", right_part_list)

  estimated_param_vec <- left_part_term %*% right_part_term

  return(estimated_param_vec)
}

#' @title Calculate the noncentrality parameter for the Wald test
.CalcWaldNcp <- function(n_clust,
                         contrast_mat,
                         null_val_vec,
                         estimated_param_vec,
                         model_based_var_mat){
 model_based_var_inv <- solve(model_based_var_mat)

  noncentrality_param <- n_clust *
    t(contrast_mat %*% estimated_param_vec - null_val_vec) %*%
    model_based_var_inv %*%
    (contrast_mat %*% estimated_param_vec - null_val_vec)

    return(noncentrality_param)

}

.CalcWaldNullCritVal <- function(alpha, test_df){
  null_crit_val <- qchisq(q = alpha,
                     df =  test_df)

  return(null_crit_val)
}

##############################################################################################
## GEE Power Validity Checks
##
##############################################################################################

#' @title Check the validity of the inputs to the GEE power calculation