## code to prepare `DATASET` dataset goes here

loadstatemcmc_str <- "./data-raw/cbc_sim_data_logit.rds"
cbc.df <- readRDS(loadstatemcmc_str)


# Convert data to sawtooth format

design <- cbc.df[,-4]
data <- cbc.df[,c(1,2,3,4)]
data <- data[data$choice != 0, ]
data.temp <- matrix(data$alt, nrow=length(unique(data$resp.id)), byrow=TRUE)
data <- data.frame(cbind(unique(data$resp.id), unique(data$resp.id), data.temp))
names(data) <- c("RespID", "Version", paste0("Task", 1:ncol(data.temp)))

write.csv(data, "./inst/extdata/simTruckData.csv", row.names = FALSE)
write.csv(design, "./inst/extdata/simTruckDesign.csv", row.names = FALSE)
