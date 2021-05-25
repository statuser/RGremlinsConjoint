# INTRODUCTION
# The files in this directory constitute the r package behind the 'Gremlins in the Data'
# JMR paper. The scripts use MCMC to estimate parameters of the model proposed in the paper.
# Authors:
#   John R. Howell (jrhowell@byu.edu)
#   Peter Ebbes (ebbes@hec.fr)
#   John C. Liechty (jcl12@psu.edu)


#' gremlins: A package for estimating the "Gremlins in the Data" model
#'
#' The gremlins package provides the tools and utilities to estimate a the model described in
#' "Gremlin's in the Data: Identifying the Information Content of Research Subjects" using
#' conjoint analysis data such as that collected in Sawtooth Software's Lighthouse or Discover
#' Products.  The packages also contains utility functions for formatting the input data and
#' extracting the relevant results.
#'
#'
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
#'   \item{price}{A continuopus variable for the relative price of the individual offerings.}
#' }
#'
"cbc.df"

#' Estimate Gremlin's Model - Hierarchical MNL
#'
#' The function estimates the model described in "Gremlin's in the Data: Identifying the
#' Information Content of Research Subjects" using a hierarchical multinomial logit model
#'
#'
#' @param data A matrix containing the raw data.
#'   The first column the respondent identifier, followed by the design number, the remaining
#'   columns indicate the choices for the tasks that coincide to the design file.
#' @param design A matrix representing the coded (dummy of effects) design file.  The design file should be
#'   formatted as a matrix with number of versions X number of tasks X number of alternatives
#'   rows and number of parameters + 3 columns.  The first column contains the version number,
#'   the second columns contains the task number, the third column contains the alternative,
#'   and the remaining columns contain the coded design.  A generic Sawtooth Software design
#'   file can be converted to this format using the \code{\link{code_sawtooth_design}} function.
#' @param Priors A data structure that contains the priors for to the model.  Can be null indicating
#'   the use of default priors or must contain a full prior specification.
#' @param R The number of repetitions in the chain
#' @param keepEvery saves every keepEvery-th draw for output
#' #' Estimate Gremlin's Model - Hierarchical MNL
#'
#' The function estimates the model described in "Gremlin's in the Data: Identifying the
#' Information Content of Research Subjects" using a hierarchical multinomial logit model
#'
#'
#' @param data A matrix containing the raw data.
#'   The first column the respondent identifier, followed by the design number, the remaining
#'   columns indicate the choices for the tasks that coincide to the design file.
#' @param design A matrix representing the coded (dummy of effects) design file.  The design file should be
#'   formatted as a matrix with number of versions X number of tasks X number of alternatives
#'   rows and number of parameters + 3 columns.  The first column contains the version number,
#'   the second columns contains the task number, the third column contains the alternative,
#'   and the remaining columns contain the coded design.  A generic Sawtooth Software design
#'   file can be converted to this format using the \code{\link{code_sawtooth_design}} function.
#' @param Priors A data structure that contains the priors for to the model.  Can be null indicating
#'   the use of default priors or must contain a full prior specification.
#' @param R The number of repetitions in the chain
#' @param keepEvery saves every keepEvery-th draw for output
#' @param nSegments (Default = 2) The number of segments for the scale factor
#' @param covariates (Optional) A matrix of covariates for the model.  One row per respondent with the
#'   respondent identified in the first column.
#' @param constraints (Optional) a vector of length n-param specifying the constraints
#'    to impose on the parameters or NULL.  a 1 indicates the parameter is constrained to be positive
#'    a -1 constrains to be  negative, and a 0 indicates no constraint.
#' @param segmentCovariates (Optional) a matrix of covariates that influence the probability that each
#'    individual belongs to the gremlins or reference group.  One row per respondent with the respondent
#'    identified in the first column.
#' @param startingValues (Optional) starting values to use for the MCMC algorithm.  This is a list of
#'    containing: slope = a nRespondent by nParamter matrix of slopes for the respondent
#'                slopeBar = a nParameter vector of the slopeBar parameter
#'                slopeCov = a nParameter by nParameter matrix containing the variance covariance matrix
#'                           for the slopeBar parameter
#'                lambda = a nSegment vector containing the starting values for the lambda parameters.  The first
#'                         element in the vector should be 1.
#'                segMembership = a nRespondent vector containing the segment membership for each respondent.
#'                phi_lambda = a nParameter vector containing the base probabilities that an individual belongs
#'                             to each segment.  Should sum to 1.
#' @param Atch_mcmc_cnt_in Parameters for the Atchade Algorithm
#' @param Atch_starting_values_slopes_in Parameters for the Atchade Algorithm
#' @param Atch_starting_values_lambda_in Parameters for the Atchade Algorithm
#' @return A data structure containing the draws from the complete MCMC chain
#'
#' @export
#' @seealso \code{\link{code_sawtooth_design}}

RGremlinsConjoint <- function() {

}
