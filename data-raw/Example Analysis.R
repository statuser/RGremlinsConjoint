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


# Tune Atchade algorithm; monitor accept rates; doesn't update if total MCMC < 1000
# Adjust Atch_tau_tune_slopes and Atch_tau_tune_lambda until accept rates are about 20-80%
# Put those numbers in the function input of the burn in


outputSimData_burn <- estimateGremlinsModel(cameraData,
                                            codedCamera,
                                            R = 4000,
                                            keepEvery = 1,
                                            num_lambda_segments = 2,
                                            Atchade_lambda_tuning = 1)

