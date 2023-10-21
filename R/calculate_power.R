# Internal functions for calculating power. These functions are not intended to be called directly by the user.
# The user-facing functions are in the file wrappers.R

##############################################################################################
## GEE Power Calculation Internal Functions
##
##############################################################################################

#' @title Internal Power Calculation Function
#'
#' Internal Power Calculation Function: .CalcPower
#'
#' @description
#' The function \code{.CalcPower} performs the core power calculation for stepped wedge trials and
#' multilevel cluster randomized trials using Generalized Estimating Equations (GEE).
#' This function is not intended to be directly called by the user, but rather serves as
#' an internal function for the main power calculation function.
#'
#' @param design_mat_list A list of design matrices, one for each cluster.
#' @param working_cor_mat_list A list of working correlation matrices, one for each cluster.
#' @param incidence_mat_list A list of incidence matrices, one for each cluster.
#' @param dgm_param_col_vec A column vector of all data generating model parameters.
#' @param dispersion_scalar The dispersion parameter for the GEE model.
#' @param var_fun A function defining the variance as a function of the mean.
#' @param link The link function for the GEE model.
#' @param response_type Type of the response variable (e.g., continuous, count, binary).
#' @param alpha The desired Type I error rate.
#' @param test_df The degrees of freedom for the test statistic.
#' @param contrast_mat A matrix specifying the contrasts for hypothesis testing.
#' @param null_val_vec A vector specifying the null hypothesis values for the contrasts.
#' @param power_only_flag A logical flag indicating whether to return only the calculated power.
#' Used internally for sample size calculations.
#'
#' @details
#' The function performs the following calculations:
#' \enumerate{
#'   \item Calculates the components for the estimating equation for each cluster: \( \boldsymbol{\mu}_i \),
#'   \( \textbf{D}_{i} \), \( \textbf{B}_{i} \), \( \textbf{R}_i(\boldsymbol{\alpha}) \), and \( \textbf{V}_{i} \).
#'   \item Solves for \( \boldsymbol{\Sigma}_1^{-1} \) and \( \widehat{\boldsymbol{\theta}} \).
#'   \item Uses the Wald test to calculate the power based on the noncentrality parameter \( \lambda \),
#'   the critical value, and the degrees of freedom \( q \).
#' }
#'
#' @return
#' The function returns a list with the following components:
#' \itemize{
#'   \item \code{power}: The power of the test.
#'   \item \code{cluster_sample_size}: The sample size for each cluster.
#'   \item \code{test_stat}: The test statistic.
#'   \item \code{alpha}: The Type I error rate.
#'   \item \code{test_df}: The degrees of freedom for the test statistic.
#' }
#'
#' @keywords internal

.CalcPower <- function(design_mat_list,
                       working_cor_mat_list,
                       incidence_mat_list = NULL,
                       dgm_param_col_vec,
                       dispersion_scalar,
                       var_fun,
                       link,
                       response_type,
                       alpha,
                       test_df,
                       contrast_mat,
                       null_val_vec,
                       power_only_flag = FALSE) {

  .CheckCalcPowerInputs(design_mat_list = design_mat_list,
                        working_cor_mat_list = working_cor_mat_list,
                        incidence_mat_list = incidence_mat_list,
                        dgm_param_col_vec = dgm_param_col_vec,
                        dispersion_scalar = dispersion_scalar,
                        var_fun = var_fun,
                        link = link,
                        response_type = response_type,
                        alpha = alpha,
                        test_df = test_df,
                        contrast_mat = contrast_mat,
                        null_val_vec = null_val_vec,
                        power_only_flag = power_only_flag)

  n_clust <- length(design_mat_list)

  # Create placeholder lists
  mui_col_vec_list <- vector(mode = "list", length = n_clust)
  di_mean_jacobian_mat_list <- vector(mode = "list", length = n_clust)
  bi_var_matrix_list <- vector(mode = "list", length = n_clust)
  vi_var_matrix_list <- vector(mode = "list", length = n_clust)

  # Calculate estimating equation components for each cluster
  for (i in 1:n_clust) {
    clust_mu_vec <- .CreateMuiVector(
      design_mat = design_mat_list[[i]],
      dgm_param_col_vec = dgm_param_col_vec,
      response_type = response_type,
      link = link
    )

    clust_mean_jacobian_mat_list <- .CreateMeanJacobianMatDi(
      mu_col_vec = clust_mu_vec,
      link = link
    )

    clust_bi_var_matrix <- .CreateVarParamMatrixBi(
      mu_col_vec = clust_mu_vec,
      dispersion_scalar = dispersion_scalar,
      var_fun = var_fun
    )

    clust_vi_var_matrix <- .CalcVarMatrixVi(
      var_mat_bi = clust_bi_var_matrix,
      working_cor_mat_ri = working_cor_mat_list[[i]]
    )

    # Store the results in the placeholder lists
    mui_col_vec_list[[i]] <- clust_mu_vec
    di_mean_jacobian_mat_list[[i]] <- clust_mean_jacobian_mat_list
    bi_var_matrix_list[[i]] <- clust_bi_var_matrix
    vi_var_matrix_list[[i]] <- clust_vi_var_matrix
  }

  # Calculate the model-based variance matrix

  model_based_var_mat <- .CalcModelBasedVariance(
    mean_jacobian_mat_list = di_mean_jacobian_mat_list,
    var_mat_vi_list = vi_var_matrix_list
  )

  # Calculate the estimated parameter vector
  estimated_param_vec <- .CalcEstimatedParamVec(
    design_mat_list = design_mat_list,
    mean_jacobian_mat_list = di_mean_jacobian_mat_list,
    var_mat_vi_list = vi_var_matrix_list,
    mui_vec_list = mui_col_vec_list,
    invlink_fun = invlink_fun
  )

  # Calculate the noncentrality parameter for the Wald test
  noncentrality_param <- .CalcWaldNcp(
    n_clust = n_clust,
    contrast_mat = contrast_mat,
    null_val_vec = null_val_vec,
    estimated_param_vec = estimated_param_vec,
    model_based_var_mat = model_based_var_mat
  )

  # Calculate the critical value for the Wald test
  null_crit_val <- .CalcWaldNullCritVal(
    alpha = alpha,
    test_df = test_df
  )

  # Calculate the power of the test
  power <- 1 - pchisq(
    q = null_crit_val,
    df = test_df,
    ncp = noncentrality_param
  )

  if (power_only_flag == TRUE) return(power)

  power_calc_results <- list(
    power = power,
    cluster_sample_size = n_clust,
    test_stat = noncentrality_param,
    alpha = alpha,
    test_df = test_df
  )

  return(power_calc_results)
}

#' Create the Expected Response Vector (\( \mu \))
#'
#' @description
#' This internal function calculates the expected response vector \( \mu \) based on the design matrix,
#' treatment parameters, response type, and link function.
#'
#' @param design_mat A numeric matrix representing the design of the experiment for a particular cluster.
#' @param dgm_param_col_vec A numeric vector containing the treatment parameters.
#' @param response_type A string indicating the type of the response variable (e.g., 'count', 'binary').
#' @param link A string specifying the link function to use (e.g., 'identity', 'logit', 'log').
#'
#' @details
#' The expected response \( \mu \) is calculated as:
#' \[
#' \mu = g^{-1}(X \beta)
#' \]
#' where \( g^{-1} \) is the inverse of the link function, \( X \) is the design matrix, and \( \beta \) is the treatment parameters.
#'
#' @return
#' Returns a numeric column vector \( \mu \) containing the expected responses.
#'
#' @keywords internal
.CreateMuiVector <- function(design_mat,
                             dgm_param_col_vec,
                             response_type,
                             link) {
  # Expected response on the scale of the link function
  linear_pred_vec <- design_mat %*% dgm_param_col_vec

  mu_vec <- switch(link,
                   identity = linear_pred_vec,
                   logit = exp(linear_pred_vec) / (1 + exp(linear_pred_vec)),
                   log = exp(linear_pred_vec)
  )

  if (response_type == "count" & min(mu_vec) <= 0) {
    stop("The mean of the response variable must be positive")
  }

  mu_col_vec <- matrix(data = mu_vec, ncol = 1)

  return(mu_col_vec)
}

#' Create the Mean Jacobian Matrix (\( D_i \))
#'
#' @description
#' This internal function calculates the mean Jacobian matrix \( D_i \), which is the derivative of the expected response \( \mu \)
#' with respect to the treatment parameters \( \beta \), based on the link function.
#'
#' @param mu_col_vec A numeric column vector containing the expected responses \( \mu \).
#' @param link A string specifying the link function to use (e.g., 'identity', 'logit', 'log').
#'
#' @details
#' The mean Jacobian matrix \( D_i \) is calculated based on the link function as follows:
#' \[
#' D_i = \frac{\partial \mu}{\partial \beta^T}
#' \]
#'
#' @return
#' Returns a numeric matrix, the mean Jacobian matrix \( D_i \).
#'
#' @keywords internal
.CreateMeanJacobianMatDi <- function(mu_col_vec, link) {
  mean_jacobian_mat_di <- switch(link,
                                 identity = diag(1, nrow = mu_col_vec),
                                 logit = diag(mu_col_vec) %*% diag(1 - mu_col_vec),
                                 log = diag(mu_col_vec)
  )

  return(mean_jacobian_mat_di)
}

#' Create the Variance Parameter Matrix (\( B_i \))
#'
#' @description
#' This internal function calculates the variance parameter matrix \( B_i \) based on the expected responses \( \mu \),
#' the dispersion parameter \( \phi \), and a variance function \( v \).
#'
#' @param mu_col_vec A numeric column vector containing the expected responses \( \mu \).
#' @param dispersion_scalar A numeric scalar representing the dispersion parameter \( \phi \).
#' @param var_fun A function defining the variance as a function of the mean.
#'
#' @details
#' The variance parameter matrix \( B_i \) is calculated as:
#' \[
#' B_i = \phi \, \text{diag}(v(\mu_{i11}), \ldots, v(\mu_{IJK}))
#' \]
#' where \( \phi \) is the dispersion parameter and \( v(\mu_{ijk}) \) is the variance.
#'
#' @return
#' Returns a numeric matrix \( B_i \), the variance parameter matrix.
#'
#' @keywords internal
.CreateVarParamMatrixBi <- function(mu_col_vec,
                                    dispersion_scalar,
                                    var_fun) {
  # concatenate matrix to vector
  var_vec <- var_fun(c(mu_col_vec))

  var_mat <- diag(var_vec)

  var_mat_bi <- dispersion_scalar * var_mat

  return(var_mat_bi)
}

#' Calculate the Covariance Matrix (\( V_i \))
#'
#' @description
#' This internal function calculates the covariance matrix \( V_i \) based on the variance parameter matrix \( B_i \)
#' and the working correlation matrix \( R_i \).
#'
#' @param var_mat_bi A numeric matrix representing the variance parameter matrix \( B_i \).
#' @param working_cor_mat_ri A numeric matrix representing the working correlation matrix \( R_i \).
#'
#' @details
#' The covariance matrix \( V_i \) is calculated as:
#' \[
#' V_i = B_i^{1/2} R_i B_i^{1/2}
#' \]
#'
#' @return
#' Returns a numeric matrix \( V_i \), the covariance matrix.
#'
#' @keywords internal
.CalcVarMatrixVi <- function(var_mat_bi,
                             working_cor_mat_ri) {
  root_bi <- sqrt(var_mat_bi)

  var_mat_vi <- root_bi %*% working_cor_mat_ri %*% root_bi

  return(var_mat_vi)
}

#' Calculate the Model-Based Variance Matrix
#'
#' @description
#' This internal function calculates the model-based variance matrix based on the list of mean Jacobian matrices \( D_i \)
#' and the list of covariance matrices \( V_i \).
#'
#' @param mean_jacobian_mat_list A list of numeric matrices representing the mean Jacobian matrices \( D_i \) for each cluster.
#' @param var_mat_vi_list A list of numeric matrices representing the covariance matrices \( V_i \) for each cluster.
#'
#' @details
#' The model-based variance matrix is calculated as:
#' \[
#' \Sigma_1^{-1} = \left[ \sum_{i = 1}^I D_{i}^T V_{i}^{-1}D_{i} \right]^{-1}
#' \]
#'
#' @return
#' Returns a numeric matrix representing the model-based variance matrix.
#'
#' @keywords internal

.CalcModelBasedVariance <- function(mean_jacobian_mat_list,
                                    var_mat_vi_list) {

  param_dim <- ncol(mean_jacobian_mat_list[[1]])

  vinv_list <- lapply(var_mat_vi_list, solve)

  model_based_var_mat <- matrix(0, nrow = param_dim, ncol = param_dim)

  for (i in 1:length(mean_jacobian_mat_list)) {
    model_based_var_mat <- model_based_var_mat +
      t(mean_jacobian_mat_list[[i]]) %*%
        vinv_list[[i]] %*%
        mean_jacobian_mat_list[[i]]
  }

  return(model_based_var_mat)
}


#' Calculate the Estimated Parameter Vector (\( \hat{\theta} \))
#'
#' @description
#' This internal function calculates the estimated parameter vector \( \hat{\theta} \)
#' based on the design matrices, mean Jacobian matrices, covariance matrices, expected responses,
#' and the inverse link function.
#'
#' @param design_mat_list A list of numeric matrices representing the design matrices \( X_i \) for each cluster.
#' @param mean_jacobian_mat_list A list of numeric matrices representing the mean Jacobian matrices \( D_i \) for each cluster.
#' @param var_mat_vi_list A list of numeric matrices representing the covariance matrices \( V_i \) for each cluster.
#' @param mui_vec_list A list of numeric vectors representing the expected responses \( \mu_i \) for each cluster.
#' @param invlink_fun A function representing the inverse link function.
#'
#' @details
#' The estimated parameter vector \( \hat{\theta} \) is calculated as:
#' \[
#' \hat{\theta} = \left( \sum_{i = 1}^I X_i^T D_i V_i D_i^T X_i \right)^{-1} \sum_{i = 1}^I X_i^T D_i V_i D_i^T g(\mu_i)
#' \]
#' where \( X_i \) are the design matrices, \( D_i \) are the mean Jacobian matrices, \( V_i \) are the covariance matrices,
#' and \( g(\mu_i) \) are the expected responses transformed by the inverse link function.
#'
#' @return
#' Returns a numeric vector \( \hat{\theta} \), the estimated parameter vector.
#'
#' @keywords internal
.CalcEstimatedParamVec <- function(design_mat_list,
                                   mean_jacobian_mat_list,
                                   var_mat_vi_list,
                                   mui_vec_list,
                                   invlink_fun) {
  num_clusters <- length(design_mat_list)

  # Calculate the shared term for each cluster
  shared_term_for_clusters <- lapply(
    1:num_clusters,
    function(i) {
      t(design_mat_list[[i]]) %*%
        mean_jacobian_mat_list[[i]] %*%
        var_mat_vi_list[[i]] %*%
        t(mean_jacobian_mat_list[[i]])
    }
  )

  # Initialize lists to hold the left and right parts of the equation
  left_equation_terms <- vector(mode = "list", length = num_clusters)
  right_equation_terms <- vector(mode = "list", length = num_clusters)

  # Calculate the left and right parts for each cluster
  for (i in 1:num_clusters) {
    left_equation_terms[[i]] <- shared_term_for_clusters[[i]] %*% design_mat_list[[i]]
    right_equation_terms[[i]] <- shared_term_for_clusters[[i]] %*% invlink_fun(mui_vec_list[[i]])
  }

  # Solve for the final left and right terms
  solved_left_term <- solve(Reduce("+", left_equation_terms))
  summed_right_term <- Reduce("+", right_equation_terms)

  # Calculate the estimated parameter vector by multiplying the left and right terms
  estimated_param_vector <- solved_left_term %*% summed_right_term

  return(estimated_param_vector)
}


#' @title Calculate the noncentrality parameter for the Wald test
#'
#' @description
#' This internal function calculates the noncentrality parameter \( \lambda \) for the Wald test.
#'
#' @param n_clust An integer representing the number of clusters.
#' @param contrast_mat A numeric matrix representing the contrast matrix \( L \).
#' @param null_val_vec A numeric vector representing the null value vector \( \ell \).
#' @param estimated_param_vec A numeric vector representing the estimated parameter vector \( \hat{\theta} \).
#' @param model_based_var_mat A numeric matrix representing the model-based variance matrix \( \Sigma_1 \).
#'
#' @details
#' The noncentrality parameter \( \lambda \) is calculated as:
#' \[
#' \lambda = I \left( L \hat{\theta} - \ell \right)^T \Sigma_1^{-1} \left( L \hat{\theta} - \ell \right)
#' \]
#'
#' @return
#' Returns a numeric value representing the noncentrality parameter \( \lambda \).
#'
#' @keywords internal
.CalcWaldNcp <- function(n_clust,
                         contrast_mat,
                         null_val_vec,
                         estimated_param_vec,
                         model_based_var_mat) {
  model_based_var_inv <- solve(model_based_var_mat)

  noncentrality_param <- n_clust *
    t(contrast_mat %*% estimated_param_vec - null_val_vec) %*%
      model_based_var_inv %*%
      (contrast_mat %*% estimated_param_vec - null_val_vec)

  return(noncentrality_param)
}


#' Calculate the Critical Value under the Null Hypothesis for the Wald Test
#'
#' @description
#' This internal function calculates the critical value under the null hypothesis for the Wald test.
#'
#' @param alpha A numeric value representing the type I error rate \( \alpha \).
#' @param test_df An integer representing the degrees of freedom.
#'
#' @details
#' The critical value under the null hypothesis is calculated using the chi-squared distribution quantile function:
#' \[
#' \chi^2_{1-\alpha, \text{df}}
#' \]
#'
#' @return
#' Returns a numeric value representing the critical value under the null hypothesis.
#'
#' @keywords internal
.CalcWaldNullCritVal <- function(alpha, test_df) {
  null_crit_val <- qchisq(
    q = 1 - alpha,
    df = test_df
  )

  return(null_crit_val)
}

###############################################################################
## GEE Power Validity Checks
##
###############################################################################
#' @keywords internal,check
.CheckCalcPowerInputs <- function(
  design_mat_list, working_cor_mat_list, incidence_mat_list, dgm_param_col_vec,
  dispersion_scalar, var_fun, link, response_type, alpha, test_df, contrast_mat,
  null_val_vec, power_only_flag) {
  checkmate::expect_list(design_mat_list,
                         any.missing = FALSE,
                         types = "numeric",
                         min.len = 1,
                         null.ok = FALSE
  )

  checkmate::expect_list(working_cor_mat_list,
                         any.missing = FALSE,
                         types = "numeric",
                         min.len = 1,
                         null.ok = FALSE
  )

  checkmate::expect_list(incidence_mat_list,
                         any.missing = FALSE,
                         types = "numeric",
                         min.len = 1,
                         null.ok = TRUE
  )

  # Check if lists are of the same length
  checkmate::expect_equal(length(design_mat_list), length(working_cor_mat_list))

  if (!is.null(incidence_mat_list)) {
    checkmate::expect_equal(length(design_mat_list), length(incidence_mat_list))
  }


  checkmate::expect_numeric(dgm_param_col_vec, any.missing = FALSE)
  checkmate::expect_number(dispersion_scalar, lower = 0)
  checkmate::expect_function(var_fun)
  checkmate::expect_function(link)
  checkmate::expect_choice(response_type, choices = c("continuous", "count", "binary")) # Replace with actual choices
  checkmate::expect_number(alpha, lower = 0, upper = 1)
  checkmate::expect_number(test_df, lower = 1)
  checkmate::expect_matrix(contrast_mat, any.missing = FALSE)
  checkmate::expect_numeric(null_val_vec, any.missing = FALSE)
  checkmate::expect_flag(power_only_flag)

  return(NULL)
}