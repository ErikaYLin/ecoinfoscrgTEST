# For planar data
sf::sf_use_s2(FALSE)

# Load packages
library(sf)
library(dplyr)

# Requires having run `wranglingCleaning_kriging.R` then `standardize_kriging.R`
## Can be sourced below
source("inst/scripts/wranglingCleaning_kriging.R")
source("inst/scripts/standardize_kriging.R")

############################
# SITEWIDE INTERPOLATION
############################

# `automap::autoKrige()` used for automated interpolation at specific sites

# Prepare spatial data
crs <- sp::CRS("EPSG:6346")  # retrieve coordinate reference system

## NUMERICAL VARS ----

# Function for interpolating numerical and binary vars
site_num <- function(site, var, formula, duplicates = TRUE) {

  # Select variable to krige for
  num.var <- as.data.frame(dplyr::select(site, var))  # automatically selects geometry
  newdat <- num.var[is.na(num.var[var]),]  # store rows where data is missing
  num.var <- num.var[!is.na(num.var[var]),]  # remove rows missing data

  # Convert to spatial object
  num.var <- sp::SpatialPointsDataFrame(coords = sf::st_coordinates(sf::st_as_sf(num.var)),
                                        data = num.var,
                                        proj4string = crs)

  # # Empirical variogram
  # variogram <- automap::autofitVariogram(var ~ 1, bin.var)
  # variogram

  # Krige for numerical/binary variable
  krige <- automap::autoKrige(as.formula(formula), num.var,
                              new_data = sp::SpatialPointsDataFrame(sf::st_coordinates(sf::st_as_sf(newdat)),
                                                                    data = as.data.frame(newdat),
                                                                    proj4string = crs),
                              remove_duplicates = duplicates, verbose = TRUE)
  return(krige)  # kriged output
}


# IRL Hab_W4
IRL_habw4_krige <- site_num(IRL, "Hab_W4", formula = "Hab_W4 ~ 1")  # kriging for unknowns
# saveRDS(IRL_habw4_krige, file = "../Data/Kriging_outputs/IRL_Hab_W4_krige.rds")
IRL_habw4_krige <- readRDS(file = "../Data/Kriging_outputs/IRL_Hab_W4_krige.rds")
IRL$Hab_W4[is.na(IRL$Hab_W4)] <- IRL_habw4_krige$krige_output$var1.pred

# Add interpolated data to statewide predictors
pred$Hab_W4[pred$study == "IRL"] <- IRL$Hab_W4


## BINARY VARS ----

# IRL Rest_Opp
IRL_restopp_krige <- site_num(IRL, var = "Rest_Opp", formula = "Rest_Opp ~ 1")
# saveRDS(IRL_restopp_krige, file = "../Data/Kriging_outputs/IRL_Rest_Opp_krige.rds")
IRL_restopp_krige <- readRDS(file = "../Data/Kriging_outputs/IRL_Rest_Opp_krige.rds")
IRL$Rest_Opp[is.na(IRL$Rest_Opp)] <- as.factor(ifelse(IRL_restopp_krige$krige_output$var1.pred >= 1.5, 1, 0))

# Add interpolated data to statewide predictors
pred$Rest_Opp[pred$study == "IRL"] <- IRL$Rest_Opp

# IRL X0yster_Pr
IRL_oysterpr_krige <- site_num(IRL, var = "X0yster_Pr", formula = "X0yster_Pr ~ 1")
# saveRDS(IRL_oysterpr_krige, file = "../Data/Kriging_outputs/IRL_X0yster_Pr_krige.rds")
IRL_oysterpr_krige <- readRDS(file = "../Data/Kriging_outputs/IRL_X0yster_Pr_krige.rds")
IRL$X0yster_Pr[is.na(IRL$X0yster_Pr)] <- as.factor(ifelse(IRL_oysterpr_krige$krige_output$var1.pred >= 1.5, 1, 0))

# Add interpolated data to statewide predictors
pred$X0yster_Pr[pred$study == "IRL"] <- IRL$X0yster_Pr

# IRL Seagrass_P
IRL_seagrass_krige <- site_num(IRL, var = "Seagrass_P", formula = "Seagrass_P ~ 1")
# saveRDS(IRL_seagrass_krige, file = "../Data/Kriging_outputs/IRL_Seagrass_P_krige.rds")
IRL_seagrass_krige <- readRDS(file = "../Data/Kriging_outputs/IRL_Seagrass_P_krige.rds")
IRL$Seagrass_P[is.na(IRL$Seagrass_P)] <- as.factor(ifelse(IRL_seagrass_krige$krige_output$var1.pred >= 1.5, 1, 0))

# Add interpolated data to statewide predictors
pred$Seagrass_P[pred$study == "IRL"] <- IRL$Seagrass_P


# tampa WideBeach
# tampa_wbeach_krige <- site_num(tampa, var = "WideBeach", formula = "WideBeach ~ 1")
## CANNOT INTERPOLATE: All observations identical and equal to 1
tampa$WideBeach[is.na(tampa$WideBeach)] <- 1

# Add interpolated data to statewide predictors
pred$WideBeach[pred$study == "tampa"] <- tampa$WideBeach

# tampa PublicRamp
# tampa_pramp_krige <- site_num(tampa, var = "PublicRamp", formula = "PublicRamp ~ 1")
## CANNOT INTERPOLATE: All observations identical and equal to 1
tampa$PublicRamp[is.na(tampa$PublicRamp)] <- 1

# Add interpolated data to statewide predictors
pred$PublicRamp[pred$study == "tampa"] <- tampa$PublicRamp

# tampa canal
# tampa_canal_krige <- site_num(tampa, var = "canal", formula = "canal ~ 1")
## CANNOT INTERPOLATE: All observations identical and equal to 1
tampa$canal[is.na(tampa$canal)] <- 1

# Add interpolated data to statewide predictors
pred$canal[pred$study == "tampa"] <- tampa$canal

# tampa SandSpit
# tampa_sandspit_krige <- site_num(tampa, var = "SandSpit", formula = "SandSpit ~ 1")
## CANNOT INTERPOLATE: All observations identical and equal to 1
tampa$SandSpit[is.na(tampa$SandSpit)] <- 1

# Add interpolated data to statewide predictors
pred$SandSpit[pred$study == "tampa"] <- tampa$SandSpit

# tampa SAV
tampa_sav_krige <- site_num(tampa, var = "SAV", formula = "SAV ~ 1")
## CANNOT INTERPOLATE: All observations identical and equal to 1
tampa$SAV[is.na(tampa$SAV)] <- 1

# Add interpolated data to statewide predictors
pred$SAV[pred$study == "tampa"] <- tampa$SAV

# tampa defended
tampa_def_krige <- site_num(tampa, var = "defended", formula = "defended ~ 1")
## CANNOT INTERPOLATE: All observations identical and equal to 1
tampa$defended[is.na(tampa$defended)] <- 1

# Add interpolated data to statewide predictors
pred$defended[pred$study == "tampa"] <- tampa$defended


# # choc canal
# choc_canal_krige <- site_num(choc, var = "canal", formula = "canal ~ 1")

# # choc SandSpit
# choc_sandspit_krige <- site_num(choc, var = "SandSpit", formula = "SandSpit ~ 1")

# choc SAV
choc_sav_krige <- site_num(choc, var = "SAV", formula = "SAV ~ 1")
# saveRDS(choc_sav_krige, file = "../Data/Kriging_outputs/choc_SAV_krige.rds")
choc$SAV[is.na(choc$SAV)] <- as.factor(ifelse(choc_sav_krige$krige_output$var1.pred >= 0.5, 1, 0))
### GENERATES ONLY NAs

# choc selectThis
choc_select_krige <- site_num(choc, var = "selectThis", formula = "selectThis ~ 1")
# saveRDS(choc_select_krige, file = "../Data/Kriging_outputs/choc_selectThis_krige.rds")
choc$selectThis[is.na(choc$selectThis)] <- as.factor(ifelse(choc_select_krige$krige_output$var1.pred >= 1.5, 1, 0))

# Add interpolated data to statewide predictors
pred$selectThis[pred$study == "choc"] <- choc$selectThis

# Save predictors
# saveRDS(pred, file = "../Data/Kriging_outputs/predictors_kriged.rds")


# Statewide angle transformation
# Transform angles to 0,90
ang <- pred %>%
  dplyr::select(angle) %>%
  as.data.frame()  # convert from sf to dataframe

# Extract and transform angles >90
ang[which(ang$angle > 90), "angle"] <- 90-abs((ang[which(ang$angle > 90), "angle"]%%180)-90)
pred$angle <- ang$angle
any(pred$angle > 90)  # check transformation

# Save predictors
# saveRDS(pred, file = "data/predictors_kriged_fix.rds")
pred <- readRDS(file = "data/predictors_kriged_fix.rds")

# Update site-specific predictors
choc <- pred[pred$study == "choc",]
IRL <- pred[pred$study == "IRL",]
pens <- pred[pred$study == "pens",]
tampa <- pred[pred$study == "tampa",]


## CATEGORICAL VARS

# # Function for interpolating categorical variables
# site_cat <- function(site, var, duplicates.rm = TRUE) {
#
#   # Select variable
#   cat.var <- as.data.frame(dplyr::select(site, starts_with(var)))  # automatically selects geometry
#
#   # Convert to spatial object
#   cat.var <- sp::SpatialPointsDataFrame(coords = sf::st_coordinates(sf::st_as_sf(cat.var)),
#                                         data = cat.var,
#                                         proj4string = crs)
#
#   # Locations to krige for
#   name <- deparse(substitute(site))
#   newdat <- miss_data %>%
#     filter(study == name) %>%
#     dplyr::select(starts_with(var))
#   newdat <- newdat[newdat[[1]] == 1,]
#
#   # Empty list
#   krige <- list()
#
#   # Krige for categorical variable (if possible)
#   for (i in 1:(ncol(cat.var)-1)) {
#     print(paste(i, "out of", ncol(cat.var)-1))  # to track progress
#     if (all(cat.var[[i]] == 1)) {
#       krige[i] <- list(paste("All identical and equal to 1"),
#                        var1.pred = 1)
#     } else if (all(cat.var[[i]] == 0)) {
#       krige[i] <- list(paste("All identical and equal to 0"),
#                        var1.pred = 0)
#     } else {
#       # Kriging formula
#       formula <- as.formula(paste(names(cat.var)[i], " ~ 1"))
#       krige[i] <- automap::autoKrige(formula = formula, cat.var,
#                                      new_data = sp::SpatialPointsDataFrame(sf::st_coordinates(sf::st_as_sf(newdat)),
#                                                                            data = as.data.frame(newdat),
#                                                                            proj4string = crs),
#                                      remove_duplicates = duplicates.rm, verbose = TRUE)
#     }
#   }
#
#   # Rename list elements to match categories
#   names(krige) <- names(cat.var)[-ncol(cat.var)]
#
#   return(krige)  # kriged output
# }

# # Example for choc PermStruc variable
# choc_Perm_krige <- site_cat(choc, "PermStruc", duplicates.rm = FALSE)
#
# # Combine all kriging outputs
# new_perm_choc <- new_perm_choc %>%
#   mutate(PermStruc_1 = choc_Perm_krige[[1]]$var1.pred,
#          PermStruc_2 = choc_Perm_krige[[2]]$var1.pred,
#          PermStruc_3 = choc_Perm_krige[[3]]$var1.pred)
#
# for (i in 1:nrow(new_perm_choc)) {
#   # Select categories with highest probabilities
#   max_cat <- which(max(c(new_perm_choc$PermStruc_1,
#                          new_perm_choc$PermStruc_2,
#                          new_perm_choc$PermStruc_3)))
#
#   if (all(new_perm_choc[i,] < 0.5)) {
#     choc$PermStruc_1[choc$geometry == new_perm_choc$geometry[i]] <- 0
#   } else if (max_cat == 1) {
#     choc$PermStruc_1[choc$geometry == new_perm_choc$geometry[i]] <- 1
#   } else if (max_cat == 2) {
#     choc$PermStruc_2[choc$geometry == new_perm_choc$geometry[i]] <- 1
#   } else if (max_cat == 3) {
#     choc$PermStruc_3[choc$geometry == new_perm_choc$geometry[i]] <- 1
#   }
# }
#
# # Add to statewide predictors
# pred$PermStruc_1[pred$study == "choc"] <- choc$PermStruc_1
# pred$PermStruc_2[pred$study == "choc"] <- choc$PermStruc_2
# pred$PermStruc_3[pred$study == "choc"] <- choc$PermStruc_3

## Kriging returned all NAs for indicator variables (Unable to interpolate)


## Categorical variables took extensive time to Krige and locations were insufficiently autocorrelated.
## As such, interpolation results have not been incorporated into the data input for model selection


## Updated predictors with spatial interpolation then fed through `pruning_kriging.R` for model selection


