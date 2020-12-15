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
#'   file can be converted to this format using the \code{codeSawtoothDesign} function.
#'  @param Priors A data structure that contains the priors for to the model.  Can be null indicating
#'   the use of default priors or must contain a full prior specification.
#'  @param numSegments (Default = 2) The number of segments for the scale factor
#'  @param R The number of repitions in the chain
#'  @param covariates (Optional) A matrix of covariates for the model.  One row per respondent with the
#'   respondent identified in the first column.
#'  @param constraints (Optional) a vector of length n-param sepecifying the constraints
#'    to impose on the parameters or NULL.  a 1 indicates the parameter is constrained to be positive
#'    a -1 constrains to be  negative, and a 0 indicates no constraint.
#'  @return A data structure containing the draws from the complete MCMC chain
#'
#'
estimateGremlinsModel <- function(data, design, Priors = NULL, R = NULL, nSegments = 2, covariates = NULL, constraints = NULL, segmentCovariates = NULL, startingValues = NULL) {

  # Validate and Format Data
  # Set Parameters and Environment variable
  # Estimate bayesm model to find starting values
  # Estimate Atchade code for full model Peter uses two steps to reduce storage 1. Burnin 2. Saved Results  There is probablya better way to find these.
  # Outstanding Questions:
  # - How do you tune the Atchade steps.  This does not seem to be automatic
  # -



}
