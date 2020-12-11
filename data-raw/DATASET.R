## code to prepare `DATASET` dataset goes here

loadstatemcmc_str <- "./data-raw/cbc_sim_data_logit.rds"
cbc.df <- readRDS(loadstatemcmc_str)

usethis::use_data(cbc.df, overwrite = TRUE)
