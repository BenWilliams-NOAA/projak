
<!-- README.md is generated from README.Rmd. Please edit that file -->

# projak

<!-- badges: start -->

<!-- badges: end -->

The goal of projak is to make an R package that can generate the
projection module in `RTMB`

## Installation

You can install the development version of projak from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("BenWilliams-NOAA/projak")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(projak)
# pull in results from a RTMB assessment
rpt = projak::m22.1 # 2024 goa northern rockfish assessment results
result = run_projections(rpt, future_catch = c("2024" = 1170.036, "2025" = 1685, "2026" = 1565))
format_output(result, var = "catch")
format_output(result, var = "f")
format_output(result, var = "ssb")
```
