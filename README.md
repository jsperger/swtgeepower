---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# swtgeepower

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/jsperger/swtgeepower/branch/main/graph/badge.svg)](https://app.codecov.io/gh/jsperger/swtgeepower?branch=main)
<!-- badges: end -->

The goal of swtgeepower is to provide a unified package for GEE-based power calculations
for stepped wedge cluster randomized trials (SWTs) and multilevel intervention stepped wedge trials
(MLI-SWTs). The package supports closed cohort, open cohort, and cross-sectional sampling designs, and
it supports incomplete designs where clusters are not observed in every study period.

Warning: this package is in alpha, and breaking changes may potentially occur.

## Design philosophy
All internal computing functions are written using clusters as the basic unit of computation.
The exported user-facing functions
are for defining a study do not yet support the full range of functionality this package is capable of.
The user-facing functions are designed to be simple and intuitive, and to support the most common use cases.


The internal functions are designed to be flexible and scalable, and power users
who can create their own cluster-level design matrices and working correlation
matrices can directly call the `foo` function and provide .

 - Supports cross-sectional, closed cohort, and open cohort sampling designs.
 - Supports incomplete designs where clusters are not observed in every study period.
 - Supports using different working correlation structures for clusters.



# Using `swtgeepower`
## Installation

You can install the development version of swtgeepower from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jsperger/swtgeepower")
```

## Example

This is a basic example which shows you how to solve a common problem:


```r
library(swtgeepower)

## Example Usage

## Prevention Science Paper Replication

prev_sci_ex <- swtgeepower::ReplicatePrevSciPaperExample()
paste0("Example power:",round(100*print(prev_sci_ex$power), 1),"%")
#> [1] 0.8092754
#> [1] "Example power:80.9%"
```
```


# Development Roadmap
 - Calculation speed improvements
 - Create classes for study design objects and working correlation matrices
 - Shiny interface
