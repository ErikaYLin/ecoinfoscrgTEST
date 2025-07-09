### Data wrangling and cleaning functions

## Cleaning and pruning data ----

# Initial data cleaning and combining
#' @export
wranglingCleaning <- function(data, response) {

  # `data` should be a named list of shape (shp) files, 1 file per study
  # `response` is a vector of shoreline suitability response variable names, 1 per study

  # For planar data
  sf::sf_use_s2(FALSE)
  # library(dplyr, verbose = FALSE)

  studies <- list()  # empty list to store shp files
  for (i in 1:length(data)) {

    NAME <- names(data)[i]  # keep name of each study
    studies[[NAME]] <- sf::st_transform(data[[i]], crs = 6346)  # match crs
  }

  # Keep original data for easy access
  studies_og <- studies

  # Convert to dataframe to avoid selecting "geometry"
  studies <- lapply(studies, as.data.frame)

  # Convert response variable to numeric levels 1-3
  for (i in 1:length(studies)) {

    if (!is.null(response)) {

      # Rename columns that should have consistent names
      ## Exposure
      if ("MxQExpCode" %in% colnames(studies[[i]])) {
        studies[[i]] <- dplyr::rename(studies[[i]], Exposure = "MxQExpCode")
      } else if ("exposure" %in% colnames(studies[[i]])) {
        studies[[i]] <- dplyr::rename(studies[[i]], Exposure = "exposure")
      }
      # WideBeach
      if ("widebeach" %in% colnames(studies[[i]])) {
        studies[[i]] <- dplyr::rename(studies[[i]], WideBeach = "widebeach")
      }
      # Beach
      if ("beach" %in% names(studies[[i]])) {
        studies[[i]] <- dplyr::rename(studies[[i]], Beach = "beach")
      }

      # Rename response column
      studies[[i]] <- dplyr::rename(studies[[i]], Response = paste0(response[i]))

      # Convert LSSM evaluation to levels 1-3
      studies[[i]] <- studies[[i]] %>%
        dplyr::mutate(Response = as.numeric(dplyr::case_when(
          Response %in% c("Maintain Beach or Offshore Breakwater with Beach Nourishment",
                          "Non-Structural Living Shoreline",
                          "Plant Marsh with Sill", "Existing Marsh Sill", "Existing Breakwater",
                          "Maintain/Enhance/Restore Riparian Buffer<br>Maintain Beach OR Offshore Breakwaters with Beach Nourishment",
                          "Maintain/Enhance/Restore Riparian Buffer<br>Maintain/Enhance/Create Marsh",
                          "Land Use Management<br>Maintain Beach OR Offshore Breakwaters with Beach Nourishment",
                          "Land Use Management<br>Maintain/Enhance/Create Marsh",
                          "Maintain/Enhance/Restore Riparian Buffer<br>Plant Marsh with Sill",
                          "Land Use Management<br>Plant Marsh with Sill",
                          "Option 2 or 5", "Option 1", "Maintain/Enhance/Restore Riparian Buffer<br>", "3") ~ "3",
          Response %in% c("Ecological Conflicts. Seek regulatory advice.",
                          "Highly Modified Area. Seek expert advice.",
                          "Land Use Management", "Land Use Management<br>",
                          "No Action Needed",
                          "Special Geomorphic Feature. Seek expert advice.",
                          "Land Use Management<br>Ecological Conflicts. Seek regulatory advice.",
                          "Maintain/Enhance/Restore Riparian Buffer<br>Ecological Conflicts. Seek regulatory advice.", "2") ~ "2",
          Response %in% c("Groin Field with Beach Nourishment",
                          "Maintain/Enhance/Restore Riparian Buffer<br>Groin Field with Beach Nourishment",
                          "Revetment", "Maintain/Enhance/Restore Riparian Buffer<br>Revetment",
                          "Revetment/Bulkhead Toe Revetment",
                          "Revetment/Bulkhead Toe Revetment Replacement",
                          "Option B3 or B4", "Option B8 or B9", "Option R7 or R8",
                          "Option B7", "Option R3 or R4", "Option 6", "1") ~ "1",
          TRUE ~ "1"  # NAs to 1
        )))
      # NOTE: other forms of living shoreline suitability evaluations will need additional steps
      ## See `Meta-Analysis_Model_Full_Workflow.Rmd` for examples of pre-processing/modifications
    }

    # Add "study" column
    studies[[i]]$study <- names(studies)[i]
  }

  # Combine Data
  state <- dplyr::bind_rows(studies)
  if (!is.null(response)) {
    pred <- state %>%
      dplyr::select(-"Response") # Remove response variables
  } else { pred <- state }

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
                        "Rest_Opp", "X0yster_Pre", "Seagrass_P", "Hardened_1", "WTLD_VEG_3")
  # NOTE: study column is excluded here for easier processing later

  # List binary variables
  binary_vars <- c("Beach", "WideBeach", "PublicRamp", "canal", "SandSpit", "SAV",
                   "defended", "selectThis", "Erosion_2", "Rest_Opp", "X0yster_Pre",
                   "Seagrass_P", "Hardened_1")

  # Update list of categorical variables (no binary)
  categorical_vars2 <- setdiff(categorical_vars, binary_vars)

  # Convert numeric variables to numeric if not already
  pred <- pred %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(numerical_vars), as.numeric))

  # Spelling and capitalization corrections (words must be chosen manually)
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
    dplyr::mutate(dplyr::across(dplyr::where(is.character), ~{
      column <- .
      for (i in 1:nrow(corrections)) {
        # Use regex to match case-insensitively
        pattern <- stringr::str_c("(?i)\\b", corrections$incorrect[i], "\\b")
        column <- stringr::str_replace_all(column, stringr::regex(pattern), corrections$correct[i])
      }
      column
    }))

  # Convert binary variables to 0/1
  for (var in binary_vars) {
    if ("Yes" %in% pred[[var]]) {
      pred[[var]] <- gsub("Yes", 1, pred[[var]])
    } else if ("Y" %in% pred[[var]]) {
      pred[[var]] <- gsub("Y", 1, pred[[var]])
    }
    if ("No" %in% pred[[var]]) {
      pred[[var]] <- gsub("No", 0, pred[[var]])
    } else if ("N" %in% pred[[var]]) {
      pred[[var]] <- gsub("N", 0, pred[[var]])
    }
  }

  # Return list of outputs
  if (!is.null(response)) {
    return(list(state = as.data.frame(state[,c("Response", "study")]),
                predictors = pred,
                numerical_vars = numerical_vars,
                categorical_vars = categorical_vars2,
                binary_vars = binary_vars))
  } else {
    return(list(state = as.data.frame(state$study),
                predictors = pred,
                numerical_vars = numerical_vars,
                categorical_vars = categorical_vars2,
                binary_vars = binary_vars))
  }
}


# Standardizing and pre-processing variables
#' @export
standardize <- function(data,
                        site.method = c("krige", "medianImpute", "meanImpute"),
                        state.method = c("meanImpute", "medianImpute"),  # kriging only available within study sites
                        duplicates = TRUE) {  # `duplicates` argument only used for kriging

  # Load dplyr quietly
  # library(dplyr, verbose = FALSE)
  sf::sf_use_s2(FALSE)  # for planar data

  # if (class(data) != "list") {
  if (!is.list(data)) {
    stop("'data' argument is not formatted properly as a list.
    Please use the same format as the output from 'wranglingCleaning()':
    list with 'state', 'predictors', 'numerical_vars', 'categorical_vars', 'binary_vars'.")
  }

  pred <- data$predictors  # extract predictors
  numerical_vars <- data$numerical_vars  # extract variables
  categorical_vars <- data$categorical_vars
  binary_vars <- data$binary_vars

  # Convert variables to proper format
  pred <- pred %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(numerical_vars), as.numeric)) %>%  # numeric if not already
    dplyr::mutate(dplyr::across(dplyr::all_of(binary_vars), as.numeric)) %>%  # numeric if not already
    dplyr::mutate(dplyr::across(dplyr::all_of(categorical_vars), as.character))  # character if not already


  ## Dummy encoding

  message("Converting categorical data to binary dummy variables...")

  # Check for categorical variables that were recorded as character strings
  ## We want these indicators to be binary dummy variables (0/1)
  dummy_df <- data.frame(row.names = seq(1:nrow(pred)))
  for (var in categorical_vars) {
    if (any(grepl(",", unique(pred[[var]])))) {
      dummy_df[[var]] <- pred[[var]]
    }
  }
  dummy_vars <- colnames(dummy_df)  # save variable names

  for (var in dummy_vars) {
    for (i in 1:nrow(pred)) {
      if (!is.na(pred[[var]][i])) {
        pred[[var]][i] <- paste0(pred[[var]][i], ".")  # add "." to each string
      }
    }
    pred[[var]] <- gsub(",", ".,", pred[[var]])  # add "." to each category in string

    dummy <- matrix(data = 0, nrow = nrow(pred),  # matrix of unique categories
                    ncol = length(unique(unlist(strsplit(as.character(pred[[var]]), ", ")))),
                    dimnames = list(NULL, c(paste0(var, "_Missing"),
                                            unique(unlist(strsplit(as.character(pred[[var]]), ", ")))[-1])))
    dummy <- as.data.frame(dummy)  # convert matrix to data frame

    # Convert categories to binaries
    for (i in 2:ncol(dummy)) {
      # Identify categories present at each observation
      dummy[which(grepl(colnames(dummy)[i], pred[[var]], fixed = TRUE)), i] <- 1
      colnames(dummy)[i] <- paste0(var, "_", i-1)
    }

    dummy[which(is.na(pred[[var]])), 1] <- 1  # find NAs
    dummy <- dummy %>%  # factor all columns
      dplyr::mutate(dplyr::across(dplyr::all_of(colnames(dummy)), as.factor))
    pred <- cbind(pred, dummy)  # add to predictors
  }

  pred <- pred %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(categorical_vars), as.factor)) # convert them to factor if not already

  # Remaining variables to convert to dummy vars
  categorical_vars2 <- setdiff(categorical_vars, dummy_vars)

  # Replace NA with "Missing" for dummy vars
  pred <- pred %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(categorical_vars2), ~ factor(ifelse(is.na(.), "Missing", .),
                                                      levels = unique(c(.,"Missing")))))

  # # Drop geometry for dummy encoding
  geom <- pred$geometry  # save geometry separately
  pred <- sf::st_drop_geometry(pred)

  # Make dummy variables for all
  for (var in categorical_vars2) {
    dummies <- stats::model.matrix(~ . - 1, data = pred[var])  # suggested to avoid intercept
    colnames(dummies) <- paste(var, levels(pred[[var]]), sep = "_")
    pred <- cbind(pred, as.data.frame(dummies))  # merge with original predictors
  }

  # Save columns with "Missing"
  miss_data <- pred %>%
    dplyr::select(c(dplyr::contains("Missing"), "study")) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(dplyr::contains("Missing")), as.factor)) #%>%
    # mutate(geometry = geom)
  # miss_data <- st_as_sf(miss_data)  # convert to spatial object

  # # Separate sites for missing data
  # miss_choc <- miss_data[miss_data$study == "choc",]
  # miss_IRL <- miss_data[miss_data$study == "IRL",]
  # miss_pens <- miss_data[miss_data$study == "pens",]
  # miss_tampa <- miss_data[miss_data$study == "tampa",]

  # Remove OG categorical columns
  pred <- pred %>%
    dplyr::select(-dplyr::all_of(categorical_vars)) %>%
    dplyr::select(-dplyr::contains(c("Length", "Lgth"))) %>%  # remove "Length" columns
    dplyr::select(-dplyr::contains("Missing"))  # remove columns with "Missing"

  # Remove spatial data artifacts
  pred <- pred %>%
    dplyr::select(-c("feature_x", "feature_y", "nearest_x", "nearest_y", "Shape__Len", "distance_2"))
                     # "field_84", "distance_2", "field_83"))

  # # Check data
  # str(pred)

  # # Separate data by study site
  # choc <- pred[pred$study == "choc",]
  # IRL <- pred[pred$study == "IRL",]
  # pens <- pred[pred$study == "pens",]
  # tampa <- pred[pred$study == "tampa",]

  # Vector of study names
  sites <- unique(pred$study)

  ## Check for missing numerical variables

  # List to store numerical variables with missing data
  miss_num <- list()

  # For each numeric/binary variable, check for NAs
  for (i in unique(pred$study)) {
    for (var in c(numerical_vars, binary_vars)) {
      # No info for entire site
      if (all(is.na(pred[pred$study == i, var]))) {
        miss_num[[i]]$statewide[[var]] <- var
        # Some info missing at site
      } else if (any(is.na(pred[pred$study == i, var]))) {
        miss_num[[i]]$sitewide[[var]] <- var
      }
    }
  }

  # Reset geometry for predictors
  sf::st_geometry(pred) <- geom
  pred <- sf::st_as_sf(pred)


  ### Site-wide Imputation

  message("Site-wide Imputation...")

  ## Auto-Kriging

  if (site.method == "krige" | is.null(site.method)) {  # krige by default

    message("Kriging selected for site-wide imputation")

    crs <- sp::CRS("EPSG:6346")  # retrieve coordinate reference system

    # Check for missing within-site (site-wide) data
    site.krige <- list()
    for (i in sites) {

      message(paste0("Sites ", which(sites == i), "/", length(sites)))

      # List missing numeric variables to krige for (by study)
      if ("sitewide" %in% names(miss_num[[i]])) {
        site.krige[[i]] <- miss_num[[i]]$sitewide
      }

      # Krige each predictor at each site
      for (var in site.krige[[i]]) {

        message(paste0("Predictor ", which(site.krige[[i]] == var), "/", length(site.krige[[i]])))

        # Check if binary (would need additional post-processing)
        if (length(unique(stats::na.omit(site.krige[var]))) <= 2 && (0 %in% site.krige[var] | 1 %in% site.krige[var])) {
          BINARY <- "YES"
        } else { BINARY <- "NO" }

        site.full <- pred[pred$study == i,]
        krige.var <- try(krigePredictors(site = pred[pred$study == i,], var = var,
                                         formula = stats::as.formula(paste0(var, " ~ 1")),
                                         duplicates = duplicates))

        if(!inherits(krige.var,'try-error')) {

          # Store kriging results
          if (BINARY == "YES") {

            # Store binary results based on a threshold
            site.full[[var]][is.na(site.full[[var]])] <-
              as.factor(ifelse(krige.var$krige_output$var1.pred >= 1.5, 1, 0))
          } else {

            # Store numeric results as is
            site.full[[var]][is.na(site.full[[var]])] <- krige.var$krige_output$var1.pred
          }

          # Add interpolated data to statewide predictors
          pred[[var]][pred$study == i] <- site.full[[var]]
        }
      }
    }
    pred <- as.data.frame(pred)  # convert back to data frame

    ## Mean Imputation

  } else if (site.method == "meanImpute") {

    message("Mean site-wide imputation selected")

    pred <- as.data.frame(pred)

    # Check for missing within-site (site-wide) data
    site.impute <- list()
    for (i in sites) {

      message(paste0("Sites ", which(sites == i), "/", length(sites)))

      # List missing numeric variables to krige for (by study)
      if ("sitewide" %in% names(miss_num[[i]])) {
        site.impute[[i]] <- miss_num[[i]]$sitewide
      }

      # Impute each predictor at each site using sitewide mean
      for (var in site.impute[[i]]) {
        site.full <- pred[pred$study == i,]
        impute.var <- meanImpute(site = pred[pred$study == i,], var = var)

        # Store imputation results
        site.full[[var]] <- impute.var

        # Add interpolated data to statewide predictors
        pred[[var]][pred$study == i] <- site.full[[var]]
      }
    }

    ## Median Imputation

  } else if (site.method == "medianImpute") {

    message("Median site-wide imputation selected")

    pred <- as.data.frame(pred)

    # Check for missing within-site (site-wide) data
    site.impute <- list()
    for (i in sites) {

      message(paste0("Sites ", which(sites == i), "/", length(sites)))

      # List missing numeric variables to krige for (by study)
      if ("sitewide" %in% names(miss_num[[i]])) {
        site.impute[[i]] <- miss_num[[i]]$sitewide
      }

      # Impute each predictor at each site using sitewide median
      for (var in site.impute[[i]]) {

        # message(paste0("Predictor ", i, "/", length(site.impute[[i]])))

        site.full <- pred[pred$study == i,]
        impute.var <- medianImpute(site = pred[pred$study == i,], var = var)

        # Store imputation results
        site.full[[var]] <- impute.var

        # Add interpolated data to statewide predictors
        pred[[var]][pred$study == i] <- site.full[[var]]
      }
    }
  }


  ### State-wide Imputation
  ## Mean Imputation

  if (state.method == "meanImpute" | is.null(state.method)) {  # mean by default

    message("Mean state-wide imputation selected")

    # Replace remaining NAs in numeric variables with state-wide mean
    pred <- pred %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(numerical_vars), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

    ## Median Imputation

  } else if (state.method == "medianImpute") {  # mean by default

    message("Median state-wide imputation selected")

    # Replace remaining NAs in numeric variables with state-wide median
    pred <- pred %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(numerical_vars), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
  }


  # # Statewide angle transformation
  # # Transform angles to 0,90
  # ang <- pred %>%
  #   dplyr::select(angle) %>%
  #   as.data.frame()  # convert from sf to dataframe
  #
  # # Extract and transform angles >90
  # ang[which(ang$angle > 90), "angle"] <- 90-abs((ang[which(ang$angle > 90), "angle"]%%180)-90)
  # pred$angle <- ang$angle


  ## Standardization

  message("Standardizing numeric predictors")

  # Save standardization MEAN and SD
  pred_num <- pred %>%
    dplyr::select(dplyr::all_of(numerical_vars))  # store numeric variables
  # MEAN <<- colMeans(pred_num, na.rm = TRUE)
  #
  # SD <- apply(pred_num, 2, sd, na.rm = TRUE)
  # names(SD) <- names(MEAN)
  # SD <<- SD

  # Load mean and sd to use for standardization
  # load("R/standardization_mean.rda")
  # load("R/standardization_sd.rda")

  # Standardize numeric predictors
  # pred <- pred %>%
  #   mutate(across(all_of(numerical_vars), ~ (.-MEAN)/SD))

  for(i in 1:ncol(pred_num)) {
    pred_num[,i] <- (pred_num[,i] - MEAN[i])/SD[i]
  }
  pred <- pred %>%
    dplyr::select(-dplyr::all_of(numerical_vars))  # remove OG numeric vars
  pred <- cbind(pred, pred_num)  # add standardized numeric vars
  ## Must use same mean and sd to standardize (transform) any new data

  # Update list of items to return
  data$predictors <- pred

  return(data)
}


## Handling missing values ----

# `automap::autoKrige()` used for automated spatial interpolation at specific sites
#' @noRd
#' @keywords internal
krigePredictors <- function(site, var, formula, duplicates = TRUE) {

  # Prepare spatial data
  crs <- sp::CRS("EPSG:6346")  # retrieve coordinate reference system

  # Select variable to krige for
  num.var <- as.data.frame(dplyr::select(site, c(var, "geometry")))  # automatically selects geometry
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
  krige <- automap::autoKrige(stats::as.formula(formula), num.var,
                              new_data = sp::SpatialPointsDataFrame(sf::st_coordinates(sf::st_as_sf(newdat)),
                                                                    data = as.data.frame(newdat),
                                                                    proj4string = crs),
                              remove_duplicates = duplicates, verbose = TRUE)
  return(krige)  # kriged output
}


# `caret::preProcess()` for median, K-nearest neighbors, or bagged imputation
#' @noRd
#' @keywords internal
medianImpute <- function(site, var, method = "medianImpute", k = 5) {

  # Select variables to impute for
  num.var <- as.data.frame(dplyr::select(site, var))  # automatically selects geometry

  # Impute
  # if (method != "bagImpute") {
  pre.params <- caret::preProcess(num.var, method = "medianImpute", verbose = T)
  # } else {
  #   pre.params <- caret::preProcess(num.var, method = method, k = k, verbose = T)
  # }

  # Create pre-processed data
  var.impute <- predict(pre.params, num.var)

  return(var.impute)  # imputed output
}


# Mean imputation with `dplyr`
#' @noRd
#' @keywords internal
meanImpute <- function(site, var) {

  # library(dplyr, verbose = FALSE)

  # Select variables to impute for
  num.var <- as.data.frame(dplyr::select(site, var))  # automatically selects geometry

  var.impute <- num.var %>%
    dplyr::mutate(dplyr::across(var, ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

  return(var.impute)
}


