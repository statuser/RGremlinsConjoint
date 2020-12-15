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
codedCamera <- code_sawtooth_design(cameraDesign, c(4:9))

#  Run bayesm with default values to set starting values.  This does not need to
#  be that accurate since we are just using this to start the chain at a
#  reasonable spot and speed up convergence

# Format the datafile for bayesm
bayesm_lgtData <- convert_to_bayesm(cameraData, codedCamera)
