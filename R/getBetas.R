# Function for organizing beta estimates
#' @export
getBetas <- function(data, predictors, Betas) {

  # library(dplyr, quietly = TRUE)

  ############################
  # GRAB MODEL OUTPUT
  ############################
  # Original
  # source("inst/scripts/wranglingCleaning.R")
  # source("inst/scripts/standardize.R")
  # Updated with Kriging
  # source("inst/scripts/wranglingCleaning_kriging.R")
  # source("inst/scripts/standardize_kriging.R")

  # Retrieve cleaned predictor data and study names
  pred <- data$predictors
  studies <- unique(data$state$study)

  # Retrieve beta estimates
  # Check if `Betas` exists
  if(!exists("Betas")) {
    stop("No betas found. Please assign beta coefficient estimates from selected models to a list named `Betas`.")
  } else { message(paste("Identified beta coefficient estimates from", length(Betas), "studies.")) }

  # Standardize betas (Z-score)
  for (i in 1:length(Betas)) {
    Betas[[i]][,1] <- (Betas[[i]][,1] + mean(Betas[[i]][,1]))/stats::sd(Betas[[i]][,1])
  }

  # Get predictors (excluding study and definitions)
  numeric_pred_cols <- predictors

  # pred <- pred %>%
  #   dplyr::mutate(dplyr::across(c("OBJECTID", "ID", "bmpCountv5", "n", "distance", "X", "Y"), as.character))
  # ## keep non predictor data as characters to avoid being selected as predictors
  # numeric_pred <- pred %>%
  #   dplyr::select_if(is.numeric)
  # factor_pred <- pred %>%
  #   dplyr::select_if(is.factor)
  # numeric_pred_cols <- colnames(cbind(factor_pred, numeric_pred))

  # Function to convert model output to dataframe (only keep "Estimate")
  prepare_df <- function(matrix, source, stat = c("beta", "se")) {
    df <- as.data.frame(t(matrix))
    if(stat == "beta") {
      df <- df[1, , drop = FALSE] # Keep "Estimate"
    } else if (stat == "se") {
      df <- df[2, , drop = FALSE] # Keep "Estimate"
    }
    colnames(df) <- rownames(matrix)
    missing_cols <- setdiff(numeric_pred_cols, rownames(matrix)) # missing predictors as NA
    df[missing_cols] <- NA
    df <- df[, numeric_pred_cols] # order columns
    df$study <- source
    return(df)
  }

  DF <- list()  # list to store dataframes of betas from each study
  for (i in 1:length(Betas)) {
    DF[[i]] <- prepare_df(Betas[[i]], studies[i], stat = "beta")
  }

  # Combine all dataframes
  combined_betas <- do.call(rbind, DF)
  combined_betas[] <- lapply(combined_betas, function(x) as.numeric(as.character(x))) # will probably get NAs

  # Get column names from pred, set any other missing columns as NA
  missing_cols <- setdiff(numeric_pred_cols, colnames(combined_betas))
  combined_betas[missing_cols] <- NA


  ############################
  # GENERATE AVERAGES AND REPLACE NAs WITH THEM
  ############################
  # Get average for each row ignore NA
  # combined_betas[] <- lapply(combined_betas, function(x) as.numeric(as.character(x)))
  row_averages <- apply(combined_betas, 1, function(row) mean(row, na.rm = TRUE))

  # Replace NAs
  for (i in 1:nrow(combined_betas)) {
    combined_betas[i, ][is.na(combined_betas[i, ])] <- row_averages[i]
  }

  # Add study column (optional)
  combined_betas$study <- studies
  # Remove study column
  combined_betas_only <- combined_betas[, !colnames(combined_betas) %in% "study"]

  message("Betas estimates combined. See objects `combined_betas` and `combined_betas_only`.")


  ################################
  # EXTRA: SE (for effect size)
  ################################

  # Standardize SE (Z-score)
  for (i in 1:length(Betas)) {
    Betas[[i]][,2] <- (Betas[[i]][,2] + mean(Betas[[i]][,2]))/stats::sd(Betas[[i]][,2])
  }

  DF_se <- list()
  for (i in 1:length(Betas)) {
    DF_se[[i]] <- prepare_df(Betas[[i]], studies[i], stat = "se")
  }

  combined_se <- do.call(rbind, DF_se)
  combined_se[] <- lapply(combined_se, function(x) as.numeric(as.character(x)))

  row_averages_se <- apply(combined_se, 1, function(row) mean(row, na.rm = TRUE))

  # choc_avg_se <- row_averages_se[1]
  # pens_avg_se <- row_averages_se[2]
  # tampa_avg_se <- row_averages_se[3]
  # IRL_avg_se <- row_averages_se[4]

  # Replace NAs
  for (i in 1:nrow(combined_se)) {
    combined_se[i, ][is.na(combined_se[i, ])] <- row_averages_se[i]
  }
  # combined_se[1, ][is.na(combined_se[1, ])] <- row_averages_se[1]
  # combined_se[2, ][is.na(combined_se[2, ])] <- row_averages_se[2]
  # combined_se[3, ][is.na(combined_se[3, ])] <- row_averages_se[3]
  # combined_se[4, ][is.na(combined_se[4, ])] <- row_averages_se[4]

  combined_betas$study <- studies
  # View(combined_se)

  # remove study column
  combined_se_only <- combined_se[, !colnames(combined_se) %in% "study"]

  message("Standard error estimates combined. See objects `combined_se` and `combined_se_only`.")

  return(list(combined_betas = combined_betas,
              combined_betas_only = combined_betas_only,
              combined_se = combined_se,
              combined_se_only = combined_se_only))

}





################################
# SCALE BETAS TO REFERENCE STUDY
################################
# # Tampa as reference
# reference_study <- combined_betas[3, ]
# # str(reference_study, list.len=ncol(reference_study))
#
# for (i in 2:(ncol(combined_betas))) { # the last column is the study
#   reference_value <- as.numeric(reference_study[i])
#
#   # Scale columns by reference study's beta value
#   combined_betas[, i] <- combined_betas[, i] / reference_value
# }
#
# # add study column (optional)
# combined_betas$study <- c("choc", "pens", "tampa", "IRL")
# # Remove study column
# combined_betas_only <- combined_betas[, !colnames(combined_betas) %in% "study"]
# # View(combined_betas_only)


# # Function to convert model out put to dataframe (only keep "Estimate")
# prepare_se_df <- function(matrix, source) {
#   df <- as.data.frame(t(matrix))
#   df <- df[2, , drop = FALSE] # Keep "Estimate"
#   colnames(df) <- rownames(matrix)
#   missing_cols <- setdiff(numeric_pred_cols, rownames(matrix)) # missing predictors as NA
#   df[missing_cols] <- NA
#   df <- df[, numeric_pred_cols] # order columns
#   df$study <- source
#   return(df)
# }

