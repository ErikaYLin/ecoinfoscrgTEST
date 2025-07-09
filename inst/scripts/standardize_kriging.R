# For planar data
sf::sf_use_s2(FALSE)

# Load packages
library(dplyr)
library(sf)

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

# # Replace NAs with means
# pred <- pred %>%
#   mutate(across(all_of(numerical_vars), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
#
# # Standardize numeric vars
# pred <- pred %>%
#   mutate(across(all_of(numerical_vars), ~ (.-mean(., na.rm = TRUE))/sd(., na.rm = TRUE)))

# Checks
## make sure sd is close to 1 and mean is close to 0
summary(pred[numerical_vars]) # summary stats
# sapply(pred[numerical_vars], sd, na.rm = TRUE) # stdev

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
                      "Erosion_1", "Erosion_2", "Owner", "Adj_H1", "Adj_H2", "Adj_H3", "Adj_H4", "V_Type1",
                      "V_Type2", "V_Type3", "V_Type4", "Rest_Opp", "X0yster_Pr", "Seagrass_P",
                      "Hardened_1", "WTLD_VEG_3")
# note that study column is excluded here for easier processing later

# List binary variables
binary_vars <- c("Beach", "WideBeach", "PublicRamp", "canal", "SandSpit", "SAV",
                 "defended", "selectThis", "Erosion_2", "Rest_Opp", "X0yster_Pr", "Seagrass_P", "Hardened_1")

# Update list of cat vars (exclude binary)
categorical_vars2 <- setdiff(categorical_vars, binary_vars)

# Fix lists of categories to separate columns
## Identify categorical variables with string of categories
for (var in categorical_vars2) {
  print(unique(pred[[var]]))
}
categorical_vars2  # Edge_Type, Hard_Mater, Adj_LU, Adj_H1, Adj_H2, Adj_H3, Adj_H4, V_Type4

# Indicator variables to convert to dummy vars
dummy_vars <- c("Edge_Type", "Hard_Mater", "Adj_LU", "Adj_H1",
                "Adj_H2", "Adj_H3", "Adj_H4", "V_Type4")

# JUST IN CASE
pred2 <- pred
dum <- pred[, dummy_vars]


# Edge_Type
edgetype <- matrix(data = 0, nrow = nrow(pred),  # matrix of unique categories
                   ncol = length(unique(unlist(strsplit(as.character(pred$Edge_Type), ",")))),
                   dimnames = list(NULL, c("Edge_Type_Missing", unique(unlist(strsplit(as.character(pred$Edge_Type), ",")))[-1])))
edgetype <- as.data.frame(edgetype)  # convert matrix to data frame
# Convert categories to binaries
for (i in 2:ncol(edgetype)) {
  # Identify categories present at each observation
  edgetype[which(grepl(colnames(edgetype)[i], pred$Edge_Type)), i] <- 1
  colnames(edgetype)[i] <- paste0("Edge_Type_", i-1)
}
edgetype[which(is.na(pred$Edge_Type)), 1] <- 1  # find NAs
edgetype <- edgetype %>%  # factorize all columns
  mutate(across(all_of(colnames(edgetype)), as.factor))
pred <- cbind(pred, edgetype)  # add to predictors

# Hard_Mater
for (i in 1:nrow(pred)) {
  if (!is.na(pred$Hard_Mater[i])) {
    pred$Hard_Mater[i] <- paste0(pred$Hard_Mater[i], ".")  # add "." to each string
  }
}
pred$Hard_Mater <- gsub(",", ".,", pred$Hard_Mater)  # add "." to each category in string

hardmater <- matrix(data = 0, nrow = nrow(pred),  # matrix of unique categories
                    ncol = length(unique(unlist(strsplit(as.character(pred$Hard_Mater), ", ")))),
                    dimnames = list(NULL, c("Hard_Mater_Missing", unique(unlist(strsplit(as.character(pred$Hard_Mater), ", ")))[-1])))
hardmater <- as.data.frame(hardmater)  # convert matrix to data frame
# Convert categories to binaries
for (i in 2:ncol(hardmater)) {
  # Identify categories present at each observation
  hardmater[which(grepl(colnames(hardmater)[i], pred$Hard_Mater, fixed = TRUE)), i] <- 1
  colnames(hardmater)[i] <- paste0("Hard_Mater_", i-1)
}
hardmater[which(is.na(pred$Hard_Mater)), 1] <- 1  # find NAs
hardmater <- hardmater %>%  # factorize all columns
  mutate(across(all_of(colnames(hardmater)), as.factor))
pred <- cbind(pred, hardmater)  # add to predictors

# Adj_LU
for (i in 1:nrow(pred)) {
  if (!is.na(pred$Adj_LU[i])) {
    pred$Adj_LU[i] <- paste0(pred$Adj_LU[i], ".")  # add "." to each string
  }
}
pred$Adj_LU <- gsub(",", ".,", pred$Adj_LU)  # add "." to each category in string

adj_lu <- matrix(data = 0, nrow = nrow(pred),  # matrix of unique categories
                 ncol = length(unique(unlist(strsplit(as.character(pred$Adj_LU), ", ")))),
                 dimnames = list(NULL, c("Adj_LU_Missing", unique(unlist(strsplit(as.character(pred$Adj_LU), ", ")))[-1])))
adj_lu <- as.data.frame(adj_lu)  # convert matrix to data frame
# Convert categories to binaries
for (i in 2:ncol(adj_lu)) {
  # Identify categories present at each observation
  adj_lu[which(grepl(colnames(adj_lu)[i], pred$Adj_LU, fixed = TRUE)), i] <- 1
  colnames(adj_lu)[i] <- paste0("Adj_LU_", i-1)
}
adj_lu[which(is.na(pred$Adj_LU)), 1] <- 1  # find NAs
adj_lu <- adj_lu %>%  # factorize all columns
  mutate(across(all_of(colnames(adj_lu)), as.factor))
pred <- cbind(pred, adj_lu)  # add to predictors

# Adj_H1
for (i in 1:nrow(pred)) {
  if (!is.na(pred$Adj_H1[i])) {
    pred$Adj_H1[i] <- paste0(pred$Adj_H1[i], ".")  # add "." to each string
  }
}
pred$Adj_H1 <- gsub(",", ".,", pred$Adj_H1)  # add "." to each category in string

adj_h1 <- matrix(data = 0, nrow = nrow(pred),  # matrix of unique categories
                 ncol = length(unique(unlist(strsplit(as.character(pred$Adj_H1), ", ")))),
                 dimnames = list(NULL, c("Adj_H1_Missing", unique(unlist(strsplit(as.character(pred$Adj_H1), ", ")))[-1])))
adj_h1 <- as.data.frame(adj_h1)  # convert matrix to data frame
# Convert categories to binaries
for (i in 2:ncol(adj_h1)) {
  # Identify categories present at each observation
  adj_h1[which(grepl(colnames(adj_h1)[i], pred$Adj_H1, fixed = TRUE)), i] <- 1
  colnames(adj_h1)[i] <- paste0("Adj_H1_", i-1)
}
adj_h1[which(is.na(pred$Adj_H1)), 1] <- 1  # find NAs
adj_h1 <- adj_h1 %>%  # factorize all columns
  mutate(across(all_of(colnames(adj_h1)), as.factor))
pred <- cbind(pred, adj_h1)  # add to predictors

# Adj_H2
for (i in 1:nrow(pred)) {
  if (!is.na(pred$Adj_H2[i])) {
    pred$Adj_H2[i] <- paste0(pred$Adj_H2[i], ".")  # add "." to each string
  }
}
pred$Adj_H2 <- gsub(",", ".,", pred$Adj_H2)  # add "." to each category in string

adj_h2 <- matrix(data = 0, nrow = nrow(pred),  # matrix of unique categories
                 ncol = length(unique(unlist(strsplit(as.character(pred$Adj_H2), ", ")))),
                 dimnames = list(NULL, c("Adj_H2_Missing", unique(unlist(strsplit(as.character(pred$Adj_H2), ", ")))[-1])))
adj_h2 <- as.data.frame(adj_h2)  # convert matrix to data frame
# Convert categories to binaries
for (i in 2:ncol(adj_h2)) {
  # Identify categories present at each observation
  adj_h2[which(grepl(colnames(adj_h2)[i], pred$Adj_H2, fixed = TRUE)), i] <- 1
  colnames(adj_h2)[i] <- paste0("Adj_H2_", i-1)
}
adj_h2[which(is.na(pred$Adj_H2)), 1] <- 1  # find NAs
adj_h2 <- adj_h2 %>%  # factorize all columns
  mutate(across(all_of(colnames(adj_h2)), as.factor))
pred <- cbind(pred, adj_h2)  # add to predictors

# Adj_H3
for (i in 1:nrow(pred)) {
  if (!is.na(pred$Adj_H3[i])) {
    pred$Adj_H3[i] <- paste0(pred$Adj_H3[i], ".")  # add "." to each string
  }
}
pred$Adj_H3 <- gsub(",", ".,", pred$Adj_H3)  # add "." to each category in string

adj_h3 <- matrix(data = 0, nrow = nrow(pred),  # matrix of unique categories
                 ncol = length(unique(unlist(strsplit(as.character(pred$Adj_H3), ", ")))),
                 dimnames = list(NULL, c("Adj_H3_Missing", unique(unlist(strsplit(as.character(pred$Adj_H3), ", ")))[-1])))
adj_h3 <- as.data.frame(adj_h3)  # convert matrix to data frame
# Convert categories to binaries
for (i in 2:ncol(adj_h3)) {
  # Identify categories present at each observation
  adj_h3[which(grepl(colnames(adj_h3)[i], pred$Adj_H3, fixed = TRUE)), i] <- 1
  colnames(adj_h3)[i] <- paste0("Adj_H3_", i-1)
}
adj_h3[which(is.na(pred$Adj_H3)), 1] <- 1  # find NAs
adj_h3 <- adj_h3 %>%  # factorize all columns
  mutate(across(all_of(colnames(adj_h3)), as.factor))
pred <- cbind(pred, adj_h3)  # add to predictors

# Adj_H4
for (i in 1:nrow(pred)) {
  if (!is.na(pred$Adj_H4[i])) {
    pred$Adj_H4[i] <- paste0(pred$Adj_H4[i], ".")  # add "." to each string
  }
}
pred$Adj_H4 <- gsub(",", ".,", pred$Adj_H4)  # add "." to each category in string

adj_h4 <- matrix(data = 0, nrow = nrow(pred),  # matrix of unique categories
                 ncol = length(unique(unlist(strsplit(as.character(pred$Adj_H4), ", ")))),
                 dimnames = list(NULL, c("Adj_H4_Missing", unique(unlist(strsplit(as.character(pred$Adj_H4), ", ")))[-1])))
adj_h4 <- as.data.frame(adj_h4)  # convert matrix to data frame
# Convert categories to binaries
for (i in 2:ncol(adj_h4)) {
  # Identify categories present at each observation
  adj_h4[which(grepl(colnames(adj_h4)[i], pred$Adj_H4, fixed = TRUE)), i] <- 1
  colnames(adj_h4)[i] <- paste0("Adj_H4_", i-1)
}
adj_h4[which(is.na(pred$Adj_H4)), 1] <- 1  # find NAs
adj_h4 <- adj_h4 %>%  # factorize all columns
  mutate(across(all_of(colnames(adj_h4)), as.factor))
pred <- cbind(pred, adj_h4)  # add to predictors

# V_Type4
for (i in 1:nrow(pred)) {
  if (!is.na(pred$V_Type4[i])) {
    pred$V_Type4[i] <- paste0(pred$V_Type4[i], ".")  # add "." to each string
  }
}
pred$V_Type4 <- gsub(",", ".,", pred$V_Type4)  # add "." to each category in string

v_type4 <- matrix(data = 0, nrow = nrow(pred),  # matrix of unique categories
                  ncol = length(unique(unlist(strsplit(as.character(pred$V_Type4), ", ")))),
                  dimnames = list(NULL, c("V_Type4_Missing", unique(unlist(strsplit(as.character(pred$V_Type4), ", ")))[-1])))
v_type4 <- as.data.frame(v_type4)  # convert matrix to data frame
# Convert categories to binaries
for (i in 2:ncol(v_type4)) {
  # Identify categories present at each observation
  v_type4[which(grepl(colnames(v_type4)[i], pred$V_Type4, fixed = TRUE)), i] <- 1
  colnames(v_type4)[i] <- paste0("V_Type4_", i-1)
}
v_type4[which(is.na(pred$V_Type4)), 1] <- 1  # find NAs
v_type4 <- v_type4 %>%  # factorize all columns
  mutate(across(all_of(colnames(v_type4)), as.factor))
pred <- cbind(pred, v_type4)  # add to predictors


pred <- pred %>%
  mutate(across(all_of(categorical_vars), as.factor)) # convert them to factor if not already


# VARIABLES WITH MISSING INFORMATION ----
## Auto-kriging will be used to interpolate for variables where possible

# Site names
site <- c("choc", "IRL", "pens", "tampa")
pred2 <- as.data.frame(pred)

# List to store numerical variables with missing data
miss_num <- list()

# For each numerical variable, check for NAs
for (i in site) {
  for (var in numerical_vars) {
    # No info for entire site
    if (all(is.na(pred2[pred2$study == i, var]))) {
      miss_num[[i]]$statewide[[var]] <- var
      # Some info missing at site
    } else if (any(is.na(pred2[pred2$study == i, var]))) {
      miss_num[[i]]$sitewide[[var]] <- var
    }
  }
}

# List to store binary variables with missing data
miss_bin <- list()

# For each categorical variable, check for NAs
for (var in binary_vars) {
  for (i in site) {
    # No info for entire site
    if (all(is.na(pred2[pred2$study == i, var]))) {
      miss_bin[[i]]$statewide[[var]] <- var
      # Some info missing at site
    } else if (any(is.na(pred2[pred2$study == i, var]))) {
      miss_bin[[i]]$sitewide[[var]] <- var
    }
  }
}

# List to store categorical variables with missing data
miss_cat <- list()

# For each categorical variable, check for NAs
for (var in categorical_vars2) {
  for (i in site) {
    # No info for entire site
    if (all(is.na(pred2[pred2$study == i, var]))) {
      miss_cat[[i]]$statewide[[var]] <- var
      # Some info missing at site
    } else if (any(is.na(pred2[pred2$study == i, var]))) {
      miss_cat[[i]]$sitewide[[var]] <- var
    }
  }
}


#########################
# DUMMY VARS
#########################

# Variables to convert to dummy vars
categorical_vars3 <- setdiff(categorical_vars2, dummy_vars)

# Replace NA with "Missing" for dummy vars
pred <- pred %>%
  mutate(across(all_of(categorical_vars3), ~ factor(ifelse(is.na(.), "Missing", .),
                                                    levels = unique(c(.,"Missing")))))

# Drop geometry for dummy encoding
geom <- pred$geometry  # save geometry separately
pred <- st_drop_geometry(pred)

# Make dummy vars for all
for (var in categorical_vars3) {
  dummies <- model.matrix(~ . - 1, data = pred[var])  # suggested to avoid intercept
  colnames(dummies) <- paste(var, levels(pred[[var]]), sep = "_")
  pred <- cbind(pred, as.data.frame(dummies))
}

# Save columns with "Missing"
miss_data <- pred %>%
  dplyr::select(c(contains("Missing"), study)) %>%
  mutate(across(all_of(contains("Missing")), as.factor)) %>%
  mutate(geometry = geom)
miss_data <- st_as_sf(miss_data)  # convert to spatial object

# Separate sites for missing data
miss_choc <- miss_data[miss_data$study == "choc",]
miss_IRL <- miss_data[miss_data$study == "IRL",]
miss_pens <- miss_data[miss_data$study == "pens",]
miss_tampa <- miss_data[miss_data$study == "tampa",]

# Remove OG categorical columns
pred <- pred %>%
  dplyr::select(-all_of(categorical_vars2)) %>%
  dplyr::select(-contains(c("Length", "Lgth"))) %>%  # remove "Length" columns
  dplyr::select(-contains("Missing"))  # remove columns with "Missing"

# Remove spatial data artifacts
pred <- pred %>%
  dplyr::select(-c("feature_x", "feature_y", "nearest_x", "nearest_y", "Shape__Len",
                   "field_84", "distance_2", "field_83"))

# Reset geometry for predictors
st_geometry(pred) <- geom
pred <- st_as_sf(pred)

# Check data
str(pred)

# Separate data by study site
choc <- pred[pred$study == "choc",]
IRL <- pred[pred$study == "IRL",]
pens <- pred[pred$study == "pens",]
tampa <- pred[pred$study == "tampa",]

# Data prepared for `kriging_predictors.R`
