library(tidyr)

#########################
# NUMERICAL VARS
#########################
# List numerical vars
numerical_vars <- c("angle", "IT_Width", "Hab_W1",
                    "Hab_W2", "Hab_W3", "Hab_W4", "Slope", "X3_m_depth", "X5_m_depth", "Slope_4",
                     "X10th", "X20th", "X30th", "X40th", "X50th", "X60th", "X70th", "X80th",
                     "X90th", "X99th", "MANGROVE")
pred <- pred %>%
    mutate(across(all_of(numerical_vars), as.numeric)) # convert them to numeric if not already

# Replace NAs with means
pred <- pred %>%
  mutate(across(all_of(numerical_vars), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# Standardize numeric vars
pred <- pred %>%
  mutate(across(all_of(numerical_vars), ~ (.-mean(., na.rm = TRUE))/sd(., na.rm = TRUE)))  ## MUST USE THIS EXACT MEAN AND SD EVERY TIME FOR NEW DATA USING THIS MODE

# Checks
## make sure sd is close to 1 and mean is close to 0
summary(pred[numerical_vars]) # summary stats
sapply(pred[numerical_vars], sd, na.rm = TRUE) # stdev

#########################
# CATEGORICAL VARS
#########################
# List cat vars
categorical_vars <- c("bnk_height", "Beach", "WideBeach", "Exposure", "bathymetry",
                      "roads", "PermStruc", "PublicRamp", "RiparianLU", "canal",
                      "SandSpit", "Structure", "offshorest", "SAV", "marsh_all",
                      "tribs", "defended", "rd_pstruc", "lowBnkStrc", "ShlType",
                      "Fetch_", "selectThis", "StrucList", "forestshl",
                      "City", "Point_Type", "Edge_Type", "Hard_Mater", "Adj_LU",
                      "Erosion", "Erosion_2", "Owner", "Adj_H1", "Adj_H2", "Adj_H3", "Adj_H4", "V_Type1",
                      "V_Type2", "V_Type3", "V_Type4", "Rest_Opp", "X0yster_Pre", "Seagrass_P",
                      "Hardened_1", "WTLD_VEG_3")
                      # note that study column is excluded here for easier processing later
pred <- pred %>%
    mutate(across(all_of(categorical_vars), as.factor))

##### DUMMY VARIABLES #####
# Replace NA with "Missing" for dummy vars
pred <- pred %>%
  mutate(across(all_of(categorical_vars), ~ factor(ifelse(is.na(.), "Missing", .), levels = unique(c(.,"Missing")))))

# Make dummy vars
for (var in categorical_vars) {
  dummies <- model.matrix(~ . - 1, data = pred[var])  # suggested to avoid itercept
  colnames(dummies) <- paste(var, levels(pred[[var]]), sep = "_")
  pred <- cbind(pred, as.data.frame(dummies))
}

# Remove OG categorical columns
pred <- pred %>%
  select(-all_of(setdiff(categorical_vars, "study")))

# Remove "Length" columns
pred <- pred %>%
  select(-contains(c("Length", "Lgth")))

# Check data
# str(pred)
