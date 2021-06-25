# INTRODUCTION
# The files in this directory constitute the R package behind the 'Gremlins in the Data'
# JMR paper. The scripts use MCMC to estimate parameters of the model proposed in the paper.
# Authors:
#   John R. Howell (jrhowell@byu.edu)
#   Peter Ebbes (ebbes@hec.fr)
#   John C. Liechty (jcl12@psu.edu)


#' 'RGremlinsConjoint': A package for estimating the "Gremlins in the Data" model
#'
#' The tools and utilities to estimate the model described in "Gremlin's in the
#' Data: Identifying the Information Content of Research Subjects" (Howell et al.
#' (2021) <doi:10.1177/0022243720965930>)  using conjoint analysis data such as
#' that collected in Sawtooth Software's 'Lighthouse' or 'Discover' products.
#' Additional utilities are included for formatting the input data.
#'
#' @docType package
#' @name gremlins
NULL


#' Set global options for the gremlins models.  These options are not expected to be modified by the user
#' but are extracted from the functions to simplify the coding.
gremlinsEnv <- new.env()

# TODO create functions to set these values
gremlinsEnv$jumpSizes <- c(0.05, 0.1, 0.5)
gremlinsEnv$jumpSizeProbs <- c(0.65, 0.25, 0.10)
gremlinsEnv$totalConstraintTries <- 100
gremlinsEnv$num_lambda_segments <- 2


#' Simulated data for the "Gremlins in the Data Model"
#'
#' A dataset containing simulated choices from a CBC study where some of the
#' respondents are information poor or 'Gremlins'.  The data is simulated data
#' and does not reflect actual preferences.
#'
#' @format A data frame with 32000 rows and 10 variables:
#' \describe{
#'   \item{resp.id}{A respondent identifier}
#'   \item{ques}{The question or task number}
#'   \item{alt}{The choice alternative}
#'   \item{choice}{An indicator that takes on a value of 1 if the alternative was chosen.  (Default is 0.)}
#'   \item{brandFord}{A dummy coded variable indicating the brand is Ford}
#'   \item{brandGM}{A dummy coded variable indicating the brand is GM}
#'   \item{brandDodge}{A dummy coded variable indicating the brand is Dodge}
#'   \item{enghyb}{A dummy coded variable indicating the engine is a hybrid}
#'   \item{engelec}{A dummy coded variable indicating the engine is electric}
#'   \item{price}{A continuous variable for the relative price of the individual offerings.}
#' }
#'
"cbc.df"
