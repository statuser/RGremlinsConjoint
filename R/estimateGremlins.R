#' Estimate Gremlin's Model - Hierarchical MNL
#'
#' The function estimates the model described in "Gremlin's in the Data: Identifying the
#' Information Content of Research Subjects" (Howell et al.
#' (2021) <doi:10.1177/0022243720965930>) using a hierarchical multinomial logit model
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
#' @param verbose Print intermediate results to the screen (default = TRUE)
#' @param num_lambda_segments (Default = 2) The number of segments for the scale factor
#' @param constraints (Optional) a vector of length n-param specifying the constraints
#'    to impose on the parameters or NULL.  a 1 indicates the parameter is constrained to be positive
#'    a -1 constrains to be  negative, and a 0 indicates no constraint.
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
#' @param previous_iterations The number of previous iterations run.  This parameter is used for the Atchade
#'      adaptive MCMC step size algorithm.  This is used since the Atchade update does not happen for less
#'      than 1000 iterations.  (Default = 0)
#' @param Atchade_slope_tuning The tuning factor to use for Atchade step for the slopes parameter.  Larger values
#'      decrease the acceptance rate.  (Default = 0.01)
#' @param Atchade_lambda_tuning The tuning factor to use for the Atchade step for the lambda parameter.  Larger
#'      values decrease the acceptance rate. (Default = 10)
#' @return A data structure containing the draws from the complete MCMC chain
#'
#' @examples
#'
#' truck_design_file <- system.file("extdata", "simTruckDesign.csv", package = "RGremlinsConjoint")
#' truck_data_file <- system.file("extdata", "simTruckData.csv", package = "RGremlinsConjoint")
#'
#' truckDesign <- read.csv(truck_design_file)
#' truckData <- read.csv(truck_data_file)
#' outputSimData_burn <- estimateGremlinsModel(truckData,
#'                                             truckDesign,
#'                                             R = 10,
#'                                             keepEvery = 1,
#'                                             num_lambda_segments = 2)
#'
#' @export
#' @seealso \code{\link{code_sawtooth_design}}
estimateGremlinsModel <- function(data,
                                  design,
                                  Priors = NULL,
                                  R = NULL,
                                  keepEvery = 1,
                                  verbose = TRUE,
                                  num_lambda_segments = 2,
                                  constraints = NULL,
                                  startingValues = NULL,
                                  previous_iterations = 0,
                                  Atchade_slope_tuning = 0.1,
                                  Atchade_lambda_tuning = 10) {

  #Set global parameters used by many of the functions
  gremlinsEnv$num_parameters <- ncol(design) - 3
  gremlinsEnv$num_respondents <- nrow(data)
  gremlinsEnv$num_tasks <- max(design[,2])
  gremlinsEnv$num_concepts <- max(design[,3])
  gremlinsEnv$num_lambda_segments <- num_lambda_segments

  # Used in Atchade steps
  gremlinsEnv$eps1 <- 10^(-7)
  gremlinsEnv$eps2 <- 10^(-6)
  gremlinsEnv$A1 <- 10^7
  gremlinsEnv$tau_bar <- 0.28		    # Targeted acceptance rate MH for beta_i's
  gremlinsEnv$tau_bar_lambda <- 0.28  # Targeted acceptance rate MH for lambda's
  gremlinsEnv$acceptanceRate_lambda <- rep(0,times=num_lambda_segments)

  if (is.null(constraints)){
    gremlinsEnv$acceptanceRate_slopes <- matrix(0, nrow=nrow(data), ncol=1)       # Monitor accept rates for each type of constraint
  }else{
    gremlinsEnv$acceptanceRate_slopes <- matrix(0, nrow=nrow(data), ncol=length(unique(constraints)))       # Monitor accept rates for each type of constraint 0, -1, 1
  }



  # Ensure design and data are actually matrices instead of data frames
  design <- as.matrix(design)
  data <- as.matrix(data)

  #Set up covariates to just an intercept if not supplied
  covariates = as.matrix(rep(1, gremlinsEnv$num_respondents))
  gremlinsEnv$nCovariates = 1

  if(is.null(constraints)) {
    constraints = double(gremlinsEnv$num_parameters)
  } else {
    if(length(constraints) != gremlinsEnv$num_parameters || !all(constraints %in% c(-1, 0, 1))) {
      stop("constraints must contain a value fo reach parameter and
           can only be -1, 0, or 1 for negative, no constraint, or positive")
    }
  }

  #validate and Set priors
  if(is.null(Priors)) {
    Priors = list(
      mu_not = matrix(0, ncol = gremlinsEnv$num_parameters, nrow = gremlinsEnv$nCovariates),
      Ainv = solve(1000*diag(ncol(covariates))),

      nu_not = gremlinsEnv$num_parameters + 2,
      V_not = 1*diag(gremlinsEnv$num_parameters),

      lambdaScale = c(1, rep(5, gremlinsEnv$num_lambda_segments - 1)),
      lambdaShape = c(1, rep(5, gremlinsEnv$num_lambda_segments - 1)),

      psi_k = rep(1, gremlinsEnv$num_lambda_segments)

    )
  } else {
    validatePriors(Priors)
  }



  respDesign <- list()
  respData <- list()

  for(ind in 1:gremlinsEnv$num_respondents) {
    respData[[ind]] <- as.matrix(data[ind, -c(1:2)]) # Drop the RespID and Design Version
    respDesign[[ind]] <- as.matrix(design[design[,1] == data[ind, 2], -c(1:3)]) # Drop the Version, Task, and Concept Columns
  }

  if(is.null(startingValues)) {
    if(verbose) {
      cat("Finding Starting Values\n")
    }
    startingValues <- find_starting_values(respData, respDesign)
  } else {
    # TODO: Validate starting values
  }

  # Initialization Atchade parameters; for each respondent a list
  Atch_MCMC_list_slopes <- generate_starting_atchade_slopes(startingValues$slopeBar,0.1)


  # Initialization Atchade parameters for lambda
  Atch_MCMC_list_lambda <- generate_starting_atchade_lambdas(10, startingValues$lambda)

  # Initialize Storage and MCMC
  slope <- array(0, dim=c(R%/%keepEvery, gremlinsEnv$num_respondents, gremlinsEnv$num_parameters))
  slope_MCMC <- startingValues$slope

  slopeBar <- array(0, dim=c(R%/%keepEvery, gremlinsEnv$num_parameters))
  slopeBar_MCMC <- matrix(startingValues$slopeBar, nrow = 1)

  slopeCov <- array(0, dim=c(R%/%keepEvery, gremlinsEnv$num_parameters, gremlinsEnv$num_parameters))
  slopeCov_MCMC <- startingValues$slopeCov

  lambda <- matrix(0, nrow = R%/%keepEvery, ncol = gremlinsEnv$num_lambda_segments)
  lambda_MCMC <- startingValues$lambda

  K <- matrix(0, nrow = R%/%keepEvery, ncol = gremlinsEnv$num_respondents)
  K_MCMC <- startingValues$segMembership

  phi_lambda <- matrix(0, nrow = R%/%keepEvery, ncol = gremlinsEnv$num_lambda_segments)
  phi_lambda_MCMC <- startingValues$phi_lambda



  if(verbose) {
    cat("Beginning MCMC Routine\n")
  }

  # Run MCMC
  for(rep in 1:R) {

      #; Provide status update on code
      if(rep %% 1000 == 1 && verbose) {
        cat(paste("Completing iteration : ", rep, "\n"))
        cat(paste("Accept rate slopes: ", colMeans(gremlinsEnv$acceptanceRate_slopes/rep), "\n"))
        cat(paste("Accept rate lambda: ", gremlinsEnv$acceptanceRate_lambda[2:num_lambda_segments]/rep, "\n"))
        cat(paste("Mu_adapt lambda:    ", Atch_MCMC_list_lambda$mu_adapt[2:num_lambda_segments], "\n"))
        cat(paste("Gamma_adapt lambda: ", Atch_MCMC_list_lambda$Gamma_adapt[2:num_lambda_segments], "\n"))
        cat(paste("metstd lambda:      ", Atch_MCMC_list_lambda$metstd[2:num_lambda_segments], "\n"))
        cat(paste("current lambda:     ", lambda_MCMC[2:num_lambda_segments], "\n"))

      }

      for(ind in 1:gremlinsEnv$num_respondents) {

        slope_atch_list_rep_ind <- generateIndividualSlope_ATCH(respData[[ind]],
                                                                respDesign[[ind]],
                                                                slope_MCMC[ind, ],
                                                                lambda_MCMC[K_MCMC[ind]],
                                                                slopeBar_MCMC,
                                                                slopeCov_MCMC, constraints,
                                                                ind_in = ind,
                                                                cur_mcmc_iter = (previous_iterations+rep),
                                                                metstd_in = Atch_MCMC_list_slopes[[ind]]$metstd,
                                                                Vmet_in = Atch_MCMC_list_slopes[[ind]]$Vmet,
                                                                Gamma_adapt_in = Atch_MCMC_list_slopes[[ind]]$Gamma_adapt,
                                                                mu_adapt_in = Atch_MCMC_list_slopes[[ind]]$mu_adapt
        )

        slope_MCMC[ind, ] <- slope_atch_list_rep_ind$currentSlope_out
        Atch_MCMC_list_slopes[[ind]] <- slope_atch_list_rep_ind$Atch_out_respi_list


        K_MCMC[ind] <- generateSegmentMembership(respData[[ind]], respDesign[[ind]], slope_MCMC[ind, ], lambda_MCMC, phi_lambda_MCMC, Priors)
      }

      multiregOut <- multivariateRegression(slope_MCMC, as.matrix(covariates), Priors)
      slopeBar_MCMC <- multiregOut$Beta
      slopeCov_MCMC <- multiregOut$Sigma

      lambda_atch_list_rep <- generateLambdaMix_ATCH(lambda_MCMC, respData, respDesign, slope_MCMC, K_MCMC, Priors,
                                                     cur_mcmc_iter = (previous_iterations+rep),
                                                     metstd_in = Atch_MCMC_list_lambda$metstd,
                                                     Vmet_in = Atch_MCMC_list_lambda$Vmet,
                                                     Gamma_adapt_in = Atch_MCMC_list_lambda$Gamma_adapt,
                                                     mu_adapt_in = Atch_MCMC_list_lambda$mu_adapt)


      lambda_MCMC <- lambda_atch_list_rep$lambda_out
      Atch_MCMC_list_lambda <- lambda_atch_list_rep$Atch_out_lambda_list


      phi_lambda_MCMC <- drawPhi(K_MCMC, Priors$psi_k)





    # Store only those draws that we intend to keep
    if(rep %% keepEvery == 0) {
      savedIterNum = rep %/% keepEvery

      slope[savedIterNum,,] <- slope_MCMC
      slopeBar[savedIterNum,] <- slopeBar_MCMC
      slopeCov[savedIterNum,,] <- slopeCov_MCMC
      lambda[savedIterNum,] <- lambda_MCMC
      K[savedIterNum,] <- K_MCMC
      phi_lambda[savedIterNum, ] <- phi_lambda_MCMC
    }


  }


  #; Save the output

  return(list(slope = slope,
              slopeBar = slopeBar,
              slopeCov = slopeCov,
              lambda = lambda,
              segmentMembership = K,
              lambdaSegProb = phi_lambda,
              Atch_MCMC_list_slopes = Atch_MCMC_list_slopes,
              Atch_MCMC_list_lambda = Atch_MCMC_list_lambda,
              nmcmc_tot = R))

}


validatePriors <- function(Priors) {

  if(is.null(Priors$mu_not) ||
     ncol(Priors$mu_not) != gremlinsEnv$num_parameters ||
     nrow(Priors$mu_not) != gremlinsEnv$nCovariates ||
     typeof(Priors$mu_not) != "double") {
    stop("mu_not should be a matrix of numerics values
         with columns equal to the number of parameters
         and rows equal to the number of covariates.")
  }

  if(is.null(Priors$Ainv) ||
     !is.matrix(Priors$Ainv) ||
     ncol(Priors$Ainv) != nrow(Priors$Ainv) ||
     ncol(Priors$Ainv) != gremlinsEnv$nCovariates ||
     typeof(Priors$Ainv) != "double") {
    stop("Ainv should be a square matrix of numeric values
         with dimensions equal to the number of covariates")
  }



  #TODO Finish validating the priors

}

drawPhi <- function(K, priors) {
  N_k <- double(length(priors))
  for(i in 1:length(priors)) {
    N_k[i] <- sum(K == i)
  }
  return(rdirichlet(N_k + priors))
}

#' @importFrom stats rgamma
rdirichlet <- function(alpha) {
  y <- rgamma(length(alpha), alpha)
  return(y/sum(y))
}
