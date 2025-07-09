#############################################
# GET RESPONSE AND STANDARDIZED PREDICTORS
#############################################
resp <- as.data.frame(cbind(state$Response, state$study))
colnames(resp) <- c("Response", "study")
# NEW!
# resp <- data.frame(Response = state$Response)
resp$Response <- as.factor(resp$Response)
resp$Response <- ordered(resp$Response)


pred <- as.data.frame(pred)

# Replace NAs with means for numerical variables
pred <- pred %>%
  mutate(across(all_of(numerical_vars), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# Standardize numeric vars  ## MUST USE THIS EXACT MEAN AND SD EVERYTIME A NEW DATASET IS FED THROUGH
pred <- pred %>%
  mutate(across(all_of(numerical_vars), ~ (.-mean(., na.rm = TRUE))/sd(., na.rm = TRUE)))

# List all predictors (AKA column names of known predictors)
numeric_pred <- pred %>%
  select_if(is.numeric)
factor_pred <- pred %>%
  select_if(is.factor)
# factor_pred <- factor_pred %>%
  # mutate(across(all_of(colnames(factor_pred)), as.character)) %>%
  # mutate(across(all_of(colnames(factor_pred)), as.numeric))
predictors <- colnames(cbind(factor_pred, numeric_pred))
predictors <- predictors[-5]  # remove SandSpit
# predictors <- predictors[-c(247,248)]  # try removing WTLD_VEG...
## Inclusion of SandSpit causes Error in str2lang(x) : <text>:2:0:
##                              unexpected end of input 1: Response ~ ^

# pred <- pred %>%
#   mutate(across(all_of(colnames(factor_pred)), as.character)) %>%
#   mutate(across(all_of(colnames(factor_pred)), as.numeric))  # convert factors to numeric

# # Find binary predictors still coded as numeric
# factor_fix <- c()
# for (i in colnames(numeric_pred)) {
#   if (length(unique(numeric_pred[[i]])) <= 2) {
#     factor_fix[i] <- i
#   }
#   factor_fix <- na.omit(factor_fix)
# }
#
# # Convert to binary vars to factor if not already
# for (var in factor_fix) {
#   if (1 %in% numeric_pred[,var] && 0 %in% numeric_pred[,var]) {  # convert if binary
#     pred[,var] <- as.factor(pred[,var])
#   }
# }

# Save standardized predictors
# saveRDS(pred, file = "data/predictors_kriged_standardized.RDS")
# pred <- readRDS("data/predictors_kriged_standardized.RDS")

# Correct for any colinear predictors

# Filter by study
resp_choc <- resp %>% filter(study == "choc")
resp_pens <- resp %>% filter(study == "pens")
resp_tampa <- resp %>% filter(study == "tampa")
resp_IRL <- resp %>% filter(study == "IRL")

pred_choc <- pred %>% filter(study == "choc")
pred_pens <- pred %>% filter(study == "pens")
pred_tampa <- pred %>% filter(study == "tampa")
pred_IRL <- pred %>% filter(study == "IRL")

##### CHOSE STUDY HERE #####
# combine response and pred
data <- cbind(resp_choc, pred_choc) # choc example

# Specify a short name of the model
name <- "choc_ordinal_bounded"
############################

# Define the response variable
response_var <- "Response"

study <- data.frame(study = state$study)
input <- cbind(resp_choc, pred_choc)
input$SMMv5Def <- NULL
input$study <- as.factor(input$study)

# Run build-up/pair-down R scripts
start_time <- Sys.time()
# source("inst/scripts/BUPD_nonparallel.R")
source("inst/scripts/BUPD.R")
end_time <- Sys.time()


## ERROR for tampa/IRL ##
# Warning in ctmm:::cov.loglike(DERIV$hess, DERIV$grad) :
#   MLE is near a boundary or optimizer failed.
# No further improvement, final model selected.
# Warning in ctmm:::cov.loglike(DERIV$hess, DERIV$grad) :
#   MLE is near a boundary or optimizer failed.

# Warning in Q[Y] - LE :
#   longer object length is not a multiple of shorter object length


#######################
# RETRIEVE BETAS
#######################

# Beta estimates were retrieved from each model
vcov.ordinal <- final_model$COV  # extract variance-covariance matrix
# Build matrix to summarize model coefficients
coefs <- matrix(NA, length(final_model$est), 2,
                dimnames = list(names(final_model$est),
                                c("Estimate", "Std. Error")))
coefs[,1] <- final_model$est  # beta estimates
coefs[,2] <- sqrt(diag(vcov.ordinal))  # standard error

# # Recreate coefficients table in model summary
# vcov.clm <- final_model$vcov  # extract VAR-COV matrix
# # Build matrix to summarize model coefficients
# coefs <- matrix(NA, length(final_model$coefficients), 4,
#                 dimnames = list(names(final_model$coefficients),
#                                 c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
# coefs[, 1] <- final_model$coefficients
# alias.clm <- unlist(final_model$aliased)
# coefs[!alias.clm, 2] <- sd <- sqrt(diag(vcov.clm))  # standard error
# coefs[!alias.clm, 3] <- coefs[!alias.clm, 1]/coefs[!alias.clm, 2]  # Z-value
# coefs[!alias.clm, 4] <- 2 * pnorm(abs(coefs[!alias.clm, 3]),
#                               lower.tail=FALSE)

# # Betas were retrieved from the model summary for each site
# average_betas <- as.data.frame(summary(final_model)[1])
# colnames(average_betas) <- c("Estimate", "Std. Error", "t value")
# average_betas <- as.matrix(average_betas)

# Save output
output_directory <- "data"
saveRDS(coefs, file = file.path(output_directory, paste0(name, "_betas.rds")))


#######################
# COMPARE MODELS
#######################

# old = prior to Kriging
# fix = after Kriging and angle transformation (missing 7 predictors)
# new = final set of predictors

# choc (~91 mins)
choc_ordinal <- readRDS("data/choc_ordinal_bounded_final_model.rds")
# Response ~ selectThis + rd_pstruc_1 + Structure_7 + Fetch__2 +
#            StrucList_2 + bathymetry_1 + offshorest_1 + StrucList_4 +
#            marsh_all_1 + RiparianLU_3 + marsh_all_7 + RiparianLU_17 +
#            RiparianLU_18 + marsh_all_5 + roads_3 + tribs_2 + RiparianLU_8 +
#            StrucList_8 + roads_1 + Structure_2 + marsh_all_6
## AIC = 7452.343

choc_mod_old <- readRDS("data/old/choc_parallel_test_final_model.rds")
choc_mod_fix <- readRDS("data/choc_non_parallel_fix_final_model.rds")
# Response ~ marsh_all_4 + rd_pstruc_1 + RiparianLU_3 + Structure_1 + RiparianLU_18 +
#           RiparianLU_7 + StrucList_4 + Exposure_2 + tribs_2 + Structure_3 + bathymetry_3 +
#           RiparianLU_8 + marsh_all_1 + Exposure_3 + marsh_all_6 + RiparianLU_13 +
#           Structure_2 + roads_3 + RiparianLU_16 + offshorest_1 + roads_1 + RiparianLU_1 +
#           RiparianLU_15 + RiparianLU_11 + PermStruc_3 + StrucList_8
## AIC = 7949.623
choc_mod_new <- readRDS("data/choc_non_parallel_new_final_model.rds")
# Response ~ selectThis + rd_pstruc_1 + Structure_7 + Exposure_2 + RiparianLU_18 +
#           bathymetry_1 + StrucList_4 + marsh_all_4 + tribs_2 + marsh_all_1 + RiparianLU_17 +
#           offshorest_1 + roads_3 + RiparianLU_3 + roads_1 + StrucList_8 + StrucList_2 +
#           RiparianLU_7 + marsh_all_6 + RiparianLU_8 + Structure_2 + RiparianLU_13 +
#           PermStruc_3
## AIC = 5705.783

# pens (~166 mins)
pens_ordinal <- readRDS("data/pens_ordinal_bounded_final_model.rds")
# Response ~ RiparianLU_7 + Beach + SAV + canal + bathymetry_1 +
#            RiparianLU_3 + SandSpit + RiparianLU_17 + marsh_all_4 + Structure_1 +
#            Exposure_1 + PermStruc_3 + marsh_all_5 + RiparianLU_15 +
#            RiparianLU_9 + offshorest_2 + tribs_1 + PublicRamp + RiparianLU_11 +
#            Structure_7 + angle + roads_1 + WideBeach + bnk_height_1
## AIC = 19212.63

pens_mod_old <- readRDS("data/old/pensTest_final_model.rds")
pens_mod_fix <- readRDS("data/pens_non_parallel_fix_final_model.rds")
# Response ~ RiparianLU_7 + bathymetry_1 + Structure_1 + Exposure_1 + PermStruc_3 +
#           RiparianLU_17 + marsh_all_4 + RiparianLU_3 + offshorest_2 + marsh_all_6 +
#           roads_3 + RiparianLU_9 + Exposure_2 + bnk_height_2 + angle + tribs_2 +
#           RiparianLU_11 + bathymetry_2 + Structure_9
## AIC = 7544.411
pens_mod_new <- readRDS("data/pens_non_parallel_new_final_model.rds")
# Response ~ RiparianLU_7 + bathymetry_1 + Beach + SAV + canal + PermStruc_2 +
#           Exposure_1 + marsh_all_5 + RiparianLU_17 + Structure_7 + RiparianLU_3 +
#           offshorest_2 + marsh_all_4 + roads_1 + defended + tribs_2 + Exposure_2 +
#           RiparianLU_9 + PublicRamp + angle
## AIC = 6510.424

# tampa (~449 mins)
tampa_ordinal <- readRDS("data/tampa_ordinal_bounded_final_model.rds")
# Response ~ tribs_2 + marsh_all_6 + forestshl_1 + RiparianLU_6 +
#            marsh_all_2 + Exposure_1 + bathymetry_1 + RiparianLU_14 +
#            RiparianLU_7 + RiparianLU_18 + bnk_height_1 + RiparianLU_3 +
#            tribs_3 + roads_3 + marsh_all_5 + RiparianLU_4 + RiparianLU_12 +
#            marsh_all_3 + RiparianLU_8 + Structure_1 + Beach + PermStruc_3 +
#            RiparianLU_5 + Exposure_3 + Structure_7 + Structure_4 + RiparianLU_2 +
#            angle + offshorest_2 + bathymetry_2
## AIC = 46389.31

tampa_mod_old <- readRDS("data/old/tampaTest_final_model.rds")
tampa_mod_fix <- readRDS("data/tampa_non_parallel_fix_final_model.rds")
# Response ~ Structure_1 + RiparianLU_8 + Structure_7 + Exposure_3 + Structure_4 +
#           Structure_5 + RiparianLU_19 + marsh_all_5 + roads_3 + Structure_10 +
#           RiparianLU_7 + marsh_all_2 + forestshl_1 + bnk_height_2 + tribs_2 +
#           RiparianLU_6 + RiparianLU_14 + RiparianLU_4 + RiparianLU_5 + bathymetry_1 +
#           tribs_3 + RiparianLU_3 + RiparianLU_10 + PermStruc_3 + RiparianLU_15 +
#           offshorest_3 + RiparianLU_11 + RiparianLU_2 + RiparianLU_17 + RiparianLU_1 +
#           marsh_all_3
## AIC = 14984.59
tampa_mod_new <- readRDS("data/tampa_non_parallel_fix_final_model.rds")
# Response ~ Structure_1 + RiparianLU_8 + Structure_7 + Exposure_3 + Structure_4 +
#           Structure_5 + RiparianLU_19 + marsh_all_5 + roads_3 + Structure_10 +
#           RiparianLU_7 + marsh_all_2 + forestshl_1 + bnk_height_2 + tribs_2 + Beach +
#           RiparianLU_6 + RiparianLU_14 + RiparianLU_4 + RiparianLU_5 + bathymetry_1 +
#           RiparianLU_18 + RiparianLU_9 + tribs_3 + PermStruc_3 + offshorest_3 +
#           RiparianLU_10 + RiparianLU_1 + Exposure_2 + RiparianLU_15 + offshorest_1 +
#           offshorest_2
## AIC = 14943.37

# IRL (~11 mins)
IRL_ordinal <- readRDS("data/IRL_ordinal_bounded_final_model.rds")
# Response ~ 1
## AIC = 2

IRL_mod_old <- readRDS("data/old/IRLTestNonParallel_final_model.rds")
IRL_mod_fix <- readRDS("data/IRL_non_parallel_fix_final_model.rds")
# Response ~ Hardened_1 + WTLD_VEG_3_2 + Slope_4 + City_5 + Erosion_1_2 + Adj_LU_7 +
#           Rest_Opp + Adj_H1_6 + City_6
## AIC = 54.28279
IRL_mod_new <- readRDS("data/IRL_non_parallel_fix_final_model.rds")
# Response ~ Hardened_1 + WTLD_VEG_3_2 + Slope_4 + City_5 + Erosion_1_2 + Adj_LU_7 +
#           Rest_Opp + Adj_H1_6 + City_6
## AIC = 54.28279




# # NEW average betas
# choc_betas <- as.data.frame(summary(choc_mod_fix)[1])
# colnames(choc_betas) <- c("Estimate", "Std. Error", "t value")
# choc_betas <- as.matrix(choc_betas)
# saveRDS(choc_betas, file = "output/choc_non_parallel_fix_average_betas.rds")
#
# pens_betas <- as.data.frame(summary(pens_mod_fix)[1])
# colnames(pens_betas) <- c("Estimate", "Std. Error", "t value")
# pens_betas <- as.matrix(pens_betas)
# saveRDS(pens_betas, file = "output/pens_non_parallel_fix_average_betas.rds")
#
# tampa_betas <- as.data.frame(summary(tampa_mod_fix)[1])
# colnames(tampa_betas) <- c("Estimate", "Std. Error", "t value")
# tampa_betas <- as.matrix(tampa_betas)
# saveRDS(tampa_betas, file = "output/tampa_non_parallel_fix_average_betas.rds")
#
# IRL_betas <- as.data.frame(summary(IRL_mod_fix)[1])
# colnames(IRL_betas) <- c("Estimate", "Std. Error", "t value")
# IRL_betas <- as.matrix(IRL_betas)
# saveRDS(IRL_betas, file = "output/IRL_non_parallel_fix_average_betas.rds")

# # OLD average betas
# chocBetas <- readRDS("output/chocContinuous_average_betas.rds")
# pensBetas <- readRDS("output/pensContinuous_average_betas.rds")
# IRLBetas <- readRDS("output/IRLContinuous_average_betas.rds")
# tampaBetas <- readRDS("output/tampaContinuous_average_betas.rds")



### ERRORS ###

# In update.uC(rho) : iteration limit reached when updating the random effects at iteration 1280


