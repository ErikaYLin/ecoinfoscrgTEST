# Cleaning data from predictors for Kriging

# For planar data
sf::sf_use_s2(FALSE)

# Load packages
library(dplyr)
# library(sp)
library(sf)
library(stringr)

# Import statewide predictors for each site
## Requires DBF files to load properly (too large for GitHub)
choc <- st_transform(st_read("../output/Final_Shapefile_all_data/Choctawatchee Bay/choc_predicted.shp"))
IRL <- st_transform(st_read("../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted.shp"))
pens <- st_transform(st_read("../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted.shp"))
tampa <- st_transform(st_read("../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted.shp"))

# Rename Erosion column due to duplicate naming
colnames(choc)[45] <- "Erosion_1"
colnames(IRL)[45] <- "Erosion_1"
colnames(pens)[45] <- "Erosion_1"
colnames(tampa)[44] <- "Erosion_1"

# Keep full data set for convenient access
choc_full <- choc
IRL_full <- IRL
pens_full <- pens
tampa_full <- tampa


####################
# REFORMAT DATA
####################

# List numerical vars
numerical_vars <- c("angle", "IT_Width", "Hab_W1",
                    "Hab_W2", "Hab_W3", "Hab_W4", "Slope", "X3_m_depth", "X5_m_depth", "Slope_4",
                    "X10th", "X20th", "X30th", "X40th", "X50th", "X60th", "X70th", "X80th",
                    "X90th", "X99th", "MANGROVE")

# List categorical vars
categorical_vars <- c("bnk_height", "Beach", "WideBeach", "Exposure", "bathymetry",
                      "roads", "PermStruc", "PublicRamp", "RiparianLU", "canal",
                      "SandSpit", "Structure", "offshorest", "SAV", "marsh_all",
                      "tribs", "defended", "rd_pstruc", "lowBnkStrc", "ShlType",
                      "Fetch_", "selectThis", "StrucList", "forestshl",
                      "City", "Point_Type", "Edge_Type", "Hard_Mater", "Adj_LU",
                      "Erosion_1", "Erosion_2", "Owner", "Adj_H1", "Adj_H2", "Adj_H3",
                      "Adj_H4", "V_Type1", "V_Type2", "V_Type3", "V_Type4",
                      "Rest_Opp", "X0yster_Pr", "Seagrass_P", "Hardened_1", "WTLD_VEG_3")
# note that study column is excluded here for easier processing later

# List binary variables
binary_vars <- c("Beach", "WideBeach", "PublicRamp", "canal", "SandSpit", "SAV",
                 "defended", "selectThis", "Erosion_2", "Rest_Opp", "X0yster_Pr", "Seagrass_P", "Hardened_1")

# Update list of categorical variables (no binary)
categorical_vars2 <- setdiff(categorical_vars, binary_vars)

# Combine data into list
pred <- list(choc = choc, IRL = IRL, pens = pens, tampa = tampa)

# Factorize categorical variables
for (i in 1:length(pred)) {
  for (j in 1:ncol(pred[[i]])) {
    pred[[i]][[j]] <- gsub("NA", NA, pred[[i]][[j]])
  }
}

# Reconvert to site-specific sf objects
choc <- pred$choc %>%
  st_set_geometry(choc_full$geometry) %>%
  st_as_sf()
IRL <- pred$IRL %>%
  st_set_geometry(IRL_full$geometry) %>%
  st_as_sf()
pens <- pred$pens %>%
  st_set_geometry(pens_full$geometry) %>%
  st_as_sf()
tampa <- pred$tampa %>%
  st_set_geometry(tampa_full$geometry) %>%
  st_as_sf()

# Combine Data
state <- dplyr::bind_rows(choc, pens, tampa, IRL)
pred <- state  %>%
  dplyr::select(-"Response") # Remove response variables

pred <- pred %>%
  mutate(across(all_of(numerical_vars), as.numeric)) # convert them to numeric if not already
# pred <- pred %>%
#   mutate(across(all_of(categorical_vars), as.factor)) # convert them to factor if not already


##################
# SPELL CHECK
##################

# Spelling and capitalization corrections (words needs to be chosen manually)
corrections <- data.frame(
  incorrect = c("YEs", "no", "NO", "No'", "RIprap", "riprap", "Permament",
                "Permanenet", "Permanenent", "Bulkead", "Bulkhea", "BUlkhead",
                "Bulkkhead", "moderate", "Scrub-shurb", "toe",
                "S, Shell, v", "R, s, Shell", "high", "low"),
  correct = c("Yes", "No", "No", "No", "Riprap", "Riprap", "Permanent",
              "Permanent", "Permanent", "Bulkhead", "Bulkhead", "Bulkhead",
              "Bulkhead", "Moderate", "Scrub-shrub", "Toe",
              "S, Shell, V", "R, S, Shell", "High", "Low")
)

# Fix misspelled words
pred <- pred %>%
  mutate(across(where(is.character), ~{
    column <- .
    for (i in 1:nrow(corrections)) {
      # for (i in seq_along(corrections$incorrect)) {
      # Use regex to match case-insensitively
      pattern <- str_c("(?i)\\b", corrections$incorrect[i], "\\b")
      column <- str_replace_all(column, regex(pattern), corrections$correct[i])
    }
    column
  }))

# Check for additional inconsistencies in data entry
for (var in categorical_vars) {
  print(unique(pred[[var]]))
}

# Fix unique cases of misspellings
pred$Adj_H2[44059] <- "V"
pred$WideBeach <- gsub("Yes", "1", pred$WideBeach)
pred$PublicRamp <- gsub("Yes", "1", pred$PublicRamp)
pred$canal <- gsub("Canal", "1", pred$canal)
pred$canal <- gsub("No", "0", pred$canal)
pred$SandSpit <- gsub("Yes", "1", pred$SandSpit)
pred$SAV <- gsub("Yes", "1", pred$SAV)
pred$SAV <- gsub("No", "0", pred$SAV)
pred$marsh_all <- gsub("No'", "No", pred$marsh_all)
pred$defended <- gsub("Yes", "1", pred$defended)
pred$selectThis <- gsub("Yes", "1", pred$selectThis)
pred$selectThis <- gsub("No", "0", pred$selectThis)
pred$Rest_Opp <- gsub("Y", "1", pred$Rest_Opp)
pred$Rest_Opp <- gsub("N", "0", pred$Rest_Opp)
pred$X0yster_Pr <- gsub("Y", "1", pred$X0yster_Pr)
pred$X0yster_Pr <- gsub("N", "0", pred$X0yster_Pr)
pred$Seagrass_P <- gsub("Y", "1", pred$Seagrass_P)
pred$Seagrass_P <- gsub("N", "0", pred$Seagrass_P)

# Data cleaned for `standardize_kriging.R`
