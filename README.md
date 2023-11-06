
<!-- README.md is generated from README.Rmd. Please edit that file -->

# swtgeepower

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/jsperger/swtgeepower/branch/main/graph/badge.svg)](https://app.codecov.io/gh/jsperger/swtgeepower?branch=main)
<!-- badges: end -->

The goal of swtgeepower is to …

## Design philosophy

All internal computing functions are written using clusters as the basic
unit of computation. The exported user-facing functions are for defining
a study do not yet support the full range of functionality this package
is capable of. The user-facing functions are designed to be simple and
intuitive, and to support the most common use cases.

The internal functions are designed to be flexible and scalable, and
power users who can create their own cluster-level design matrices and
working correlation matrices can directly call the `foo` function and
provide .

- Supports cross-sectional, closed cohort, and open cohort sampling
  designs.
- Supports incomplete designs where clusters are not observed in every
  study period.
- Supports using different working correlation structures for clusters.

# Using `swtgeepower`

## Installation

You can install the development version of swtgeepower from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jsperger/swtgeepower")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(swtgeepower)
## basic example code
```

# Development Roadmap

- shiny interface
- 
