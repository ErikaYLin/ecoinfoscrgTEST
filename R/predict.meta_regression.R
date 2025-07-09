# Predicting from the meta-analysis model
#' @export
predict.meta_regression <- function(data, meta_regression) {

  # Retrieve predictor data and variable names
  pred <- data$predictors
  predictors <- unique(meta_regression$data$predictor)  # retrieve predictors

  # Remove non-predictor columns
  pred <- pred %>%
    dplyr::select(dplyr::all_of(predictors))
  predictor_data <- pred[, match(predictors, colnames(pred))]  # reorder columns

  # Reformat predictor columns
  predictor_data <- predictor_data %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(predictors), as.character)) %>% # convert predictors to character first (needed for factorized columns)
    dplyr::mutate(dplyr::across(dplyr::all_of(predictors), as.numeric)) # convert all predictors to numeric
    # dplyr::mutate(dplyr::across(dplyr::all_of(predictors), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))  # sub NAs for statewide means
    # dplyr::mutate(dplyr::across(dplyr::all_of(predictors), ~ ifelse(is.na(.), MEAN_all, .)))

  for (i in 1:length(predictors)) {
    if(anyNA(predictor_data[,i])){
      predictor_data[which(is.na(predictor_data[,i])),i] <- MEAN_all[i]
    }
  }
  ## NAs in the data can result in predicted values to be NA

  # Extract coefficients and intercept from meta-analysis model
  betas <- meta_regression$beta

  # Remove "predictor" from predictor names
  rownames(betas) <- gsub("predictor", "", rownames(betas))  # match predictor names
  setdiff(colnames(predictor_data), rownames(betas))  # which predictor used as the intercept?
  rownames(betas)[1] <- setdiff(colnames(predictor_data), rownames(betas))  # replace with predictor name
  betas <- as.vector(betas)

  # Multiply predictors by the meta-analysis model betas
  expected_outcome <- as.matrix(predictor_data) %*% betas

  # Inverse probit for predicted values
  predicted <- pnorm(expected_outcome)

  return(predicted)

}



