
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RGremlinsConjoint

<!-- badges: start -->

[![R-CMD-check](https://github.com/statuser/RGremlinsConjoint/workflows/R-CMD-check/badge.svg)](https://github.com/statuser/RGremlinsConjoint/actions)
[![R-CMD-check](https://github.com/statuser/RGremlinsConjoint/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/statuser/RGremlinsConjoint/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The gremlins package provides the tools and utilities to estimate a the
model described in “Gremlins in the Data: Identifying the Information
Content of Research Subjects”
(\[<https://doi.org/10.1177/0022243720965930>\]) using conjoint analysis
data such as that collected in Sawtooth Software’s Lighthouse or
Discover Products. The packages also contains utility functions for
formatting the input data and extracting the relevant results.

## Installation

<!-- You can install the released version of RGremlinsConjoint from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("RGremlinsConjoint") -->
<!-- ``` -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("statuser/RGremlinsConjoint")
```

## Example

The package exposes basically one function You can use it like:

``` r
library(RGremlinsConjoint)

# Read in the data
truck_design_file <- system.file("extdata", "simTruckDesign.csv", package = "RGremlinsConjoint")
truck_data_file <- system.file("extdata", "simTruckData.csv", package = "RGremlinsConjoint")
truckDesign <- read.csv(truck_design_file)
truckData <- read.csv(truck_data_file)

# Covert the design file to be dummy coded is necessary
# The simulated data is already coded
# codedTruck <- code_sawtooth_design(truckDesign, c(4:9), include_none_option=TRUE)

outputSimData_burn <- estimateGremlinsModel(truckData,
                                            truckDesign,
                                            R = 10,
                                            keepEvery = 1,
                                            num_lambda_segments = 2)
#> Finding Starting Values
#> Beginning MCMC Routine
#> Completing iteration :  1 
#> Accept rate slopes:  0 
#> Accept rate lambda:  0 
#> Mu_adapt lambda:     50 
#> Gamma_adapt lambda:  10 
#> metstd lambda:       10 
#> current lambda:      50
```
