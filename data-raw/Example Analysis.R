# To install the package you will need to run the following command
# devtools::install_github("statuser/RGremlinsConjoint")
library(RGremlinsConjoint)

truck_design_file <- system.file("extdata", "simTruckDesign.csv", package = "RGremlinsConjoint")
truck_data_file <- system.file("extdata", "simTruckData.csv", package = "RGremlinsConjoint")

# Read in the Sawtooth Formatted data
truckDesign <- read.csv(truck_design_file)
truckData <- read.csv(truck_data_file)

## Covert the design file to be dummy coded
# the truck data is already coded
## codedTruck <- code_sawtooth_design(truckDesign, c(4:9), include_none_option=TRUE)


# Tune Atchade algorithm; monitor accept rates; doesn't update if total MCMC < 1000
# Adjust Atch_tau_tune_slopes and Atch_tau_tune_lambda until accept rates are about 20-80%
# Put those numbers in the function input of the burn in


outputSimData_burn <- estimateGremlinsModel(truckData,
                                            truckDesign,
                                            R = 1010,
                                            keepEvery = 1,
                                            num_lambda_segments = 2,
                                            Atchade_lambda_tuning = 1)

