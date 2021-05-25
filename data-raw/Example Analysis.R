# To install the package you will need to run the following command
# devtools::install_github("statuser/RGremlinsConjoint")
library(RGremlinsConjoint)

camera_design_file <- system.file("extdata", "CameraDesign.csv", package = "RGremlinsConjoint")
camera_data_file <- system.file("extdata", "CameraFullData.csv", package = "RGremlinsConjoint")

# Read in the Sawtooth Formatted data
cameraDesign <- read.csv(camera_design_file)
cameraData <- read.csv(camera_data_file)

## Covert the design file to be dummy coded
price_list = c(0.79, 1.29, 1.79, 2.29, 2.79)

cameraDesign$price_lin <- price_list[cameraDesign$price_lin]
codedCamera <- code_sawtooth_design(cameraDesign, c(4:9), include_none_option=TRUE)


num_parameters <- ncol(codedCamera) - 3
num_covariates <- 1
# Set Priors as needed
Set_Priors = list(
  mu_not = matrix(0, ncol=num_parameters, nrow=num_covariates),
  Ainv = solve(1000*diag(num_covariates)),
  nu_not = num_parameters + 2,
  V_not = 1*diag(num_parameters),
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

# TODO - Check and update the setting or priors into the package code
# Figure out how to pass the appropriate values to the Atchade tuning steps and package that code


outputSimData_burn <- estimateGremlinsModel(cameraData,
                                            codedCamera,
                                            Priors = Set_Priors,
                                            R = 4000,
                                            keepEvery = 1,
                                            num_lambda_segments = 2,
                                            Atch_mcmc_cnt_in = Atch_cnt,
                                            Atch_starting_values_slopes_in = Atch_starting_values_slopes,
                                            Atch_starting_values_lambda_in = Atch_starting_values_lambda)

