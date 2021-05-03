library(bayesm)

#  Load the example Data File
data(cbc.df)

# convert to bayesm format
cbc_bayesm.df <- cbc.df


# Read in the Sawtooth Formatted data
cameraDesign <- read.csv("data-raw/CameraDesign.csv")
cameraData <- read.csv("data-raw/CameraFullData.csv")

## Covert the design file to be dummy coded
price_list = c(0.79, 1.29, 1.79, 2.29, 2.79)

cameraDesign$price_lin <- price_list[cameraDesign$price_lin]
codedCamera <- code_sawtooth_design(cameraDesign, c(4:9), include_none_option=TRUE)

#  Run bayesm with default values to set starting values.  This does not need to
#  be that accurate since we are just using this to start the chain at a
#  reasonable spot and speed up convergence

gremlinsEnv$num_lambda_segments <- 2

# Format the datafile for bayesm
bayesm_data <- convert_to_bayesm(cameraData, codedCamera)

bayesm_prior <- list(ncomp=1)
starting_value_iters <- 2000
starting_value_thin <- 10
starting_value_burnin <- 1750
bayesm_mcmc <- list(R=starting_value_iters*starting_value_thin, keep=starting_value_thin, nprint = 0)

invisible(capture.output(bayesm_results <- rhierMnlRwMixture(Data=bayesm_data, Prior=bayesm_prior, Mcmc=bayesm_mcmc)))

# Set Starting Values
slope_start <- apply(bayesm_results$betadraw[,,starting_value_burnin:starting_value_iters], c(1,2), mean)
slopeBar_start <- colMeans(slope_start)
slopeCov_start <- diag(length(slopeBar_start))
lambda_start <- seq(1, 50, length.out = num_lambda_segments)
segment_membership_start <- calculate_segment_memberships(bayesm_results$betadraw[,,starting_value_burnin:starting_value_iters], bayesm_data)
phi_lambda <- as.numeric(table(segment_membership_start)/length(segment_membership_start))

startingValues <- list(slope = slope_start,
                       slopeBar = slopeBar_start,
                       slopeCov = slopeCov_start,
                       lambda = lambda_start,
                       segMembership = segment_membership_start,
                       phi_lambda = phi_lambda)

# Set Priors as needed
nParams <- ncol(codedCamera) - 3
nUnits <- nrow(cameraData)
nCovariates <- 1
Set_Priors = list(
  mu_not = matrix(0, ncol=nParams, nrow=nCovariates),
  Ainv = solve(1000*diag(nCovariates)),
  nu_not = nParams + 2,
  V_not = 1*diag(nParams),
  lambdaScale = c(NA, 5),
  lambdaShape = c(NA, 5),
  psi_k = rep(1, 2)
)

Atch_starting_values_slopes<-generate_starting_atchade_slopes(startingValues$slopeBar,0.1,length(startingValues$slopeBar), nrow(startingValues$slope))
Atch_starting_values_lambda<-generate_starting_atchade_lambdas(15, 2, c(NA,20), c(NA,15))
Atch_cnt <- 0
niter_atchade_past <- 0


# Tune Atchade algorithm; monitor accept rates; doesn't update if total MCMC < 1000
# Adjust Atch_tau_tune_slopes and Atch_tau_tune_lambda until accept rates are about 20-80%
# Put those numbers in the function input of the burn in

outputSimData_burn <- estimateGremlinsModel(cameraData,
                                            codedCamera,
                                            startingValues = startingValues,
                                            Priors = Set_Priors,
                                            R = 4000,
                                            keepEvery = 1,
                                            nSegments = 2,
                                            Atch_mcmc_cnt_in = Atch_cnt,
                                            Atch_starting_values_slopes_in = Atch_starting_values_slopes,
                                            Atch_starting_values_lambda_in = Atch_starting_values_lambda)

