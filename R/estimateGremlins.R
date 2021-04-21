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
#'  @param R The number of repetitions in the chain
#'  @param keepEvery saves every keepEvery-th draw for output
#'  @param numSegments (Default = 2) The number of segments for the scale factor
#'  @param covariates (Optional) A matrix of covariates for the model.  One row per respondent with the
#'   respondent identified in the first column.
#'  @param constraints (Optional) a vector of length n-param specifying the constraints
#'    to impose on the parameters or NULL.  a 1 indicates the parameter is constrained to be positive
#'    a -1 constrains to be  negative, and a 0 indicates no constraint.
#'  @param segmentCovariates (Optional) a matrix of covariates that influence the probability that each
#'    individual belongs to the gremlins or reference group.  One row per respondent with the respondent
#'    identified in the first column.
#'  @param startingValues (Optional) starting values to use for the MCMC algorithm.  This is a list of
#'    containing: slope = a nRespondent by nParamter matrix of slopes for the respondent
#'                slopeBar = a nParameter vector of the slopeBar parameter
#'                slopeCov = a nParameter by nParameter matrix containing the variance covariance matrix
#'                           for the slopeBar parameter
#'                lambda = a nSegment vector containing the starting values for the lambda parameters.  The first
#'                         element in the vector should be 1.
#'                segMembership = a nRespondent vector containing the segment membership for each respondent.
#'                phi_lambda = a nParameter vector containing the base probabilities that an individual belongs
#'                             to each segment.  Should sum to 1.
#'  @param Atch_mcmc_cnt_in ...
#'  @param Atch_starting_values_slopes_in ...
#'  @param Atch_starting_values_lambda_in
#'  @return A data structure containing the draws from the complete MCMC chain
#'
#'
estimateGremlinsModel <- function(data,
                                  design,
                                  Priors = NULL,
                                  R = NULL,
                                  keepEvery = 1,
                                  nSegments = 2,
                                  covariates = NULL,
                                  constraints = NULL,
                                  segmentCovariates = NULL,
                                  startingValues = NULL,
                                  Atch_mcmc_cnt_in,
                                  Atch_starting_values_slopes_in,
                                  Atch_starting_values_lambda_in) {

  #Set global parameters used by many of the functions
  gremlinsEnv$nParams <- ncol(design) - 3
  gremlinsEnv$nUnits <- nrow(data)
  gremlinsEnv$nTasks <- max(design[,2])
  gremlinsEnv$nConcepts <- max(design[,3])
  gremlinsEnv$nSegments <- nSegments

  # Used in Atchade steps
  gremlinsEnv$eps1 <- 10^(-7)
  gremlinsEnv$eps2 <- 10^(-6)
  gremlinsEnv$A1 <- 10^7
  gremlinsEnv$tau_bar <- 0.28		    # Targeted acceptance rate MH for beta_i's
  gremlinsEnv$tau_bar_lambda <- 0.28  # Targeted acceptance rate MH for lambda's
  gremlinsEnv$acceptanceRate_lambda <- rep(0,times=nSegments)

  if (is.null(constraints)){
    gremlinsEnv$acceptanceRate_slopes <- matrix(0, nrow=nrow(data), ncol=1)       # Monitor accept rates for each type of constraint
  }else{
    gremlinsEnv$acceptanceRate_slopes <- matrix(0, nrow=nrow(data), ncol=length(unique(constraints)))       # Monitor accept rates for each type of constraint 0, -1, 1
  }



  # Ensure design and data are actually matrices instead of data frames
  design <- as.matrix(design)
  data <- as.matrix(data)

  #Check for covariates on segment Variables
  if(!is.null(segmentCovariates)) {
    gremlinsEnv$useDelta = TRUE
    gremlinsEnv$nDeltaParams = ncol(segmentCovariates)
  } else {
    gremlinsEnv$useDelta = FALSE
  }

  #Set up covariates to just an intercept if not supplied
  if(is.null(covariates)) {
    covariates = as.matrix(rep(1, gremlinsEnv$nUnits))
    gremlinsEnv$nCovariates = 1
  } else {
    if(!is.matrix(covariates) || nrow(covariates) != gremlinsEnv$nUnits) {
      stop("covariates must be a matrix with
           nrows equal to the number of units")
    }
    gremlinsEnv$nCovariates = ncol(covariates)
  }

  if(is.null(constraints)) {
    constraints = double(gremlinsEnv$nParams)
  } else {
    if(length(constraints) != gremlinsEnv$nParams || !all(constraints %in% c(-1, 0, 1))) {
      stop("constraints must contain a value fo reach parameter and
           can only be -1, 0, or 1 for negative, no constraint, or positive")
    }
  }

  #validate and Set priors
  if(is.null(Priors)) {
    Priors = list(
      mu_not = matrix(0, ncol = gremlinsEnv$nParams, nrow = gremlinsEnv$nCovariates),
      Ainv = solve(1000*diag(ncol(covariates))),

      nu_not = gremlinsEnv$nParams + 2,
      V_not = 1*diag(gremlinsEnv$nParams),

      lambdaScale = c(NA, 5),
      lambdaShape = c(NA, 5),

      psi_k = rep(1, gremlinsEnv$nSegments)

    )
  } else {
    validatePriors(Priors, useDelta = useDelta)
  }

  if(is.null(startingValues)) {
    # TODO: We have a method using bayesm.  We need to convert this to the function and
    # print appropriate timing and progress
    # startingValues = findStartingValues()

    startingValues <- list(slope = matrix(0, nrow = gremlinsEnv$nUnits, ncol = gremlinsEnv$nParams),
                           slopeBar = double(gremlinsEnv$nParams),
                           slopeCov = diag(gremlinsEnv$nParams),
                           lambda = seq(1, by=25, length.out = gremlinsEnv$nSegments),
                           segMembership = sample(1:gremlinsEnv$nSegments, size = gremlinsEnv$nUnits, replace = TRUE),
                           phi_lambda = rep(1/gremlinsEnv$nSegments, gremlinsEnv$nSegments))
  } else {
    # TODO: Validate starting values
  }

  # Initialization Atchade parameters; for each respondent a list
  Atch_MCMC_list_slopes <- Atch_starting_values_slopes_in


  # Initialization Atchade parameters for lambda
  Atch_MCMC_list_lambda <- Atch_starting_values_lambda_in

  # Initialize Storage and MCMC
  slope <- array(0, dim=c(R%/%keepEvery, gremlinsEnv$nUnits, gremlinsEnv$nParams))
  slope_MCMC <- startingValues$slope

  slopeBar <- array(0, dim=c(R%/%keepEvery, gremlinsEnv$nParams))
  slopeBar_MCMC <- startingValues$slopeBar

  slopeCov <- array(0, dim=c(R%/%keepEvery, gremlinsEnv$nParams, gremlinsEnv$nParams))
  slopeCov_MCMC <- startingValues$slopeCov

  lambda <- matrix(0, nrow = R%/%keepEvery, ncol = gremlinsEnv$nSegments)
  lambda_MCMC <- startingValues$lambda

  K <- matrix(0, nrow = R%/%keepEvery, ncol = gremlinsEnv$nUnits)
  K_MCMC <- startingValues$segMembership

  phi_lambda <- matrix(0, nrow = R%/%keepEvery, ncol = gremlinsEnv$nSegments)
  phi_lambda_MCMC <- startingValues$phi_lambda

  respDesign <- list()
  respData <- list()

  for(ind in 1:gremlinsEnv$nUnits) {
    respData[[ind]] <- data[ind, -c(1:2)] # Drop the RespID and Design Version
    respDesign[[ind]] <- design[design[,1] == data[ind, 2], -c(1:3)] # Drop the Version, Task, and Concept Columns
  }

  cat("Beginning MCMC Routine\n")

  # Run MCMC
  for(rep in 1:R) {

      #; Provide status update on code
      if(rep %% 1000 == 0) {
        cat(paste("Completing iteration : ", R, "\n"))
        cat(paste("Accept rate slopes: ", colMeans(gremlinsEnv$acceptanceRate_slopes/rep), "\n"))
        cat(paste("Accept rate lambda: ", gremlinsEnv$acceptanceRate_lambda[2:nSegments]/rep, "\n"))
        cat(paste("Mu_adapt lambda:    ", Atch_MCMC_list_lambda$mu_adapt[2:nSegments], "\n"))
        cat(paste("Gamma_adapt lambda: ", Atch_MCMC_list_lambda$Gamma_adapt[2:nSegments], "\n"))
        cat(paste("metstd lambda:      ", Atch_MCMC_list_lambda$metstd[2:nSegments], "\n"))
        cat(paste("current lambda:     ", lambda_MCMC[2:nSegments], "\n"))

      }

      for(ind in 1:gremlinsEnv$nUnits) {

        slope_atch_list_rep_ind <- generateIndividualSlope_ATCH(respData[[ind]],
                                                                respDesign[[ind]],
                                                                slope_MCMC[ind, ],
                                                                lambda_MCMC[K_MCMC[ind]],
                                                                slopeBar_MCMC[L_MCMC[ind],],
                                                                slopeCov_MCMC[L_MCMC[ind],,], constraints,
                                                                ind_in = ind,
                                                                cur_mcmc_iter = (Atch_mcmc_cnt_in+rep),
                                                                metstd_in = Atch_MCMC_list_slopes[[ind]]$metstd,
                                                                Vmet_in = Atch_MCMC_list_slopes[[ind]]$Vmet,
                                                                Gamma_adapt_in = Atch_MCMC_list_slopes[[ind]]$Gamma_adapt,
                                                                mu_adapt_in = Atch_MCMC_list_slopes[[ind]]$mu_adapt
        )

        slope_MCMC[ind, ] <- slope_atch_list_rep_ind$currentSlope_out
        Atch_MCMC_list_slopes[[ind]] <- slope_atch_list_rep_ind$Atch_out_respi_list


        K_MCMC[ind] <- generateSegmentMembership(respData[[ind]], respDesign[[ind]], slope_MCMC[ind, ], lambda_MCMC, phi_lambda_MCMC, Priors)
        L_MCMC[ind] <- generateBetaSegmentMembership(slope_MCMC[ind, ], slopeBar_MCMC, slopeCov_MCMC, phi_beta_MCMC)
      }

      for(l in 1:gremlinsEnv$nBetaSegments) {
        multiregOut <- multivariateRegression(slope_MCMC[L_MCMC == l,], as.matrix(covariates[L_MCMC == l]), Priors)
        slopeBar_MCMC[l,] <- multiregOut$Beta
        slopeCov_MCMC[l,,] <- multiregOut$Sigma
      }

      lambda_atch_list_rep <- generateLambdaMix_ATCH(lambda_MCMC, respData, respDesign, slope_MCMC, K_MCMC, Priors,
                                                     cur_mcmc_iter = (Atch_mcmc_cnt_in+rep),
                                                     metstd_in = Atch_MCMC_list_lambda$metstd,
                                                     Vmet_in = Atch_MCMC_list_lambda$Vmet,
                                                     Gamma_adapt_in = Atch_MCMC_list_lambda$Gamma_adapt,
                                                     mu_adapt_in = Atch_MCMC_list_lambda$mu_adapt)


      lambda_MCMC <- lambda_atch_list_rep$lambda_out
      Atch_MCMC_list_lambda <- lambda_atch_list_rep$Atch_out_lambda_list


      phi_lambda_MCMC <- drawPhi(K_MCMC, Priors$psi_k)
      phi_beta_MCMC <- drawPhi(L_MCMC, Priors$psi_l)





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
              betaSegmentMembership = L,
              lambdaSegProb = phi_lambda,
              betaSegProb = phi_beta,
              Atch_MCMC_list_slopes = Atch_MCMC_list_slopes,
              Atch_MCMC_list_lambda = Atch_MCMC_list_lambda,
              nmcmc_tot = nreps_total))

}
