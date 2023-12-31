% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_power.R
\name{.CreateMuiVector}
\alias{.CreateMuiVector}
\title{Create the Expected Response Vector (\( \mu \))}
\usage{
.CreateMuiVector(linear_pred_vec, family)
}
\arguments{
\item{linear_pred_vec}{}

\item{family}{A family object (from base R \code{stats} class)}

\item{dgm_param_col_vec}{A numeric vector containing the treatment parameters.}
}
\value{
Returns a numeric column vector \( \mu \) containing the expected responses.
}
\description{
This internal function calculates the expected response vector \( \mu \) based on the design matrix,
treatment parameters, response type, and (inverse) link function.
}
\details{
The expected response \( \mu \) is calculated as:
[
\mu = g^{-1}(X \beta)
]
where \( g^{-1} \) is the inverse of the link function, \( X \) is the design matrix, and \( \beta \) is the treatment parameters.
}
\keyword{internal}
