# Internal data

load("inst/scripts/standardization_sd.rda")
load("inst/scripts/standardization_mean.rda")
load("inst/scripts/standardization_mean_all.rda")

usethis::use_data(SD, MEAN, MEAN_all, internal = TRUE, overwrite = TRUE)
