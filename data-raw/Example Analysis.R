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

num_lambda_segments <- 2

# Format the datafile for bayesm
bayesm_data <- convert_to_bayesm(cameraData, codedCamera)

bayesm_prior <- list(ncomp=1)
starting_value_iters <- 2000
starting_value_thin <- 10
starting_values_burnin <- 1750
bayesm_mcmc <- list(R=starting_value_iters*starting_value_thin, keep=starting_value_thin, nprint = 0)

invisible(capture.output(bayesm_results <- rhierMnlRwMixture(Data=bayesm_data, Prior=bayesm_prior, Mcmc=bayesm_mcmc)))


slope_start <- apply(bayesm_results$betadraw[,,starting_values_burnin:starting_value_iters], c(1,2), mean)
slopeBar_start <- colMeans(betai_start)
slopeCov_start <- diag(length(slopeBar_start))
lambda_start <- seq(1, 50, length.out = num_lambda_segments)
segment_membership_start <- calculate_segment_memberships(bayesm_results$betadraw[,,starting_values_burnin:starting_value_iters], bayesm_data)
phi_lambda <- as.numeric(table(segment_membership_start)/length(segment_membership_start))



