% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_power.R
\name{.CalcModelBasedVariance}
\alias{.CalcModelBasedVariance}
\title{Calculate the Model-Based Variance Matrix}
\usage{
.CalcModelBasedVariance(n_clust, mean_jacobian_mat_list, var_mat_vi_list)
}
\arguments{
\item{mean_jacobian_mat_list}{A list of numeric matrices representing the mean Jacobian matrices \( D_i \) for each cluster.}

\item{var_mat_vi_list}{A list of numeric matrices representing the covariance matrices \( V_i \) for each cluster.}
}
\value{
Returns a numeric matrix representing the model-based variance matrix.
}
\description{
This internal function calculates the model-based variance matrix based on the list of mean Jacobian matrices \( D_i \)
and the list of covariance matrices \( V_i \).
}
\details{
The model-based variance matrix is calculated as:
[
\Sigma_1^{-1} = \left\link[= \\sum_{i = 1}^I D_{i}^T V_{i}^{-1}D_{i} \\right]{ \sum_{i = 1}^I D_{i}^T V_{i}^{-1}D_{i} \right}^{-1}
]
}
\keyword{internal}
