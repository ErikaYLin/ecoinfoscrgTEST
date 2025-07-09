# Function for model selection
#' @export
BUPD <- function(data, predictors, parallel = FALSE) {

  # Retrieve response data
  resp <- data$state
  resp$Response <- as.integer(resp$Response)  # must be an integer sequence from 1:K

  # Define the response variable
  response_var <- "Response"

  # Retrieve predictor data
  pred <- data$predictors
  pred <- pred %>%
    dplyr::select(dplyr::all_of(predictors))
  pred$study <- resp$study

  # Identify the studies
  studies <- unique(resp$study)  # vector of study names

  # Filter by study
  resp_all <- list()
  pred_all <- list()
  data_all <- list()
  input <- list()
  for (s in 1:length(studies)){
    resp_all[[s]] <- resp %>% dplyr::filter(study == studies[[s]])
    pred_all[[s]] <- pred %>% dplyr::filter(study == studies[[s]])

    # Merge response and predictors
    data_all[[s]] <- cbind(resp_all[[s]], pred_all[[s]])
    input[[s]] <- cbind(resp_all[[s]], pred_all[[s]])
    input[[s]]$study <- as.factor(input[[s]]$study)
  }

  predictors_all <- predictors  # save predictor names

  RESULTS <- list()

  if (parallel) {

    # Use total cores - 1 to parallelize
    numCores <- parallel::detectCores() - 1

    # Repeat for every study
    for (s in 1:length(studies)) {

      # start_time <- Sys.time()  # to track run time per study
      tictoc::tic()  # to track run time per study

      # Initial formula (intercept only)
      initial_formula <- stats::as.formula(paste(response_var, "~ 1"))
      best_formula <- initial_formula

      # Initial model
      best_model <- ordinal(best_formula, data = input[[s]])  # ordinal probit regression
      best_aic <- best_model$AIC

      ############################################
      # BUILD-UP PHASE
      ############################################
      # Function to fit and evaluate models
      fit_model <- function(formula, data) {
        tryCatch({
          model <- ordinal(formula, data = data)
          aic <- model$AIC
          return(list(formula_str = paste(deparse(formula, width.cutoff = 500), collapse = ""),
                      model = model, aic = aic, error = NULL))
        }, error = function(e) {
          return(list(formula_str = NULL, model = NULL, aic = Inf, error = e$message))
        })
      }

      # Remove predictors with zero variance
      ZV <- c()
      for (var in predictors_all) {
        if (length(unique(input[[s]][[var]])) == 1 && is.na(unique(input[[s]][[var]]))) {
          ZV[var] <- var
        }
      }
      predictors_all <- setdiff(predictors_all, ZV)

      # Build-up phase
      for (i in 1:length(predictors_all)) {
        current_predictors <- all.vars(best_formula)[-1]
        remaining_predictors <- setdiff(predictors_all, current_predictors)

        # Use mclapply to parallelize the model fitting with forked processes
        results <- parallel::mclapply(remaining_predictors, function(predictor) {
          new_formula <- stats::update(best_formula, paste(". ~ . +", predictor))
          fit <- fit_model(new_formula, input[[s]])

          if (!is.null(fit$error)) {
            cat("Error for predictor", predictor, ":", fit$error, "\n")
            return(NULL)
          }

          formula_str <- paste(deparse(new_formula, width.cutoff = 500), collapse = "")
          cat("Tested formula:", formula_str, "AIC:", fit$aic, "\n")
          list(formula_str = formula_str, aic = fit$aic)
        }, mc.cores = numCores)

        # Filter out results with errors
        results <- Filter(function(x) !is.null(x$formula_str), results)

        if (length(results) > 0) {
          best_candidate <- which.min(sapply(results, function(x) x$aic))
          best_candidate_formula <- results[[best_candidate]]$formula_str
          best_candidate_aic <- results[[best_candidate]]$aic

          if (best_candidate_aic < best_aic) {
            best_formula <- stats::as.formula(best_candidate_formula)
            best_aic <- best_candidate_aic
            cat("New best model:", best_candidate_formula, "AIC:", best_aic, "\n")
          } else {
            cat("No further improvement, stopping build-up.\n")
            break
          }
        } else {
          cat("No more predictors to test, stopping build-up.\n")
          break
        }
      }

      # Final model
      best_model <- ordinal(best_formula, data = input[[s]])
      best_aic <- best_model$AIC

      ############################################
      # PAIR-DOWN PHASE
      ############################################
      current_formula <- best_formula  # start with best FORMULA from build-up phase
      current_aic <- best_model$AIC

      # >>>>> WORKS TO THIS POINT (~66 preds, ~ 12 hours) <<<<<
      # 08/09/2024: Seems to be working for pairdown; outputting formulas as it tests. Waiting for it to finish. If need to, can save outut from buildup for later testing.

      repeat {
        predictors_in_model <- all.vars(current_formula)[-1]

        # Use mclapply to parallelize the pair-down model fitting with forked processes
        results <- parallel::mclapply(predictors_in_model, function(predictor) {
          if(length(predictors_in_model <= 1)) {
            pairdown_formula <- stats::as.formula(paste(response_var, "~", 1))
          } else {
            pairdown_formula <- stats::as.formula(paste(response_var, "~", paste(setdiff(predictors_in_model, predictor), collapse = "+")))
          }
          if (length(all.vars(pairdown_formula)[-1]) == 0) {
            return(NULL)
          }

          fit <- fit_model(pairdown_formula, input[[s]])
          formula_str <- paste(deparse(pairdown_formula, width.cutoff = 500), collapse = "")
          list(formula_str = formula_str, aic = fit$aic)
        }, mc.cores = numCores)

        # Filter out results with errors
        results <- Filter(function(x) !is.null(x$formula_str), results)

        if (length(results) > 0) {
          best_pairdown_aic <- min(sapply(results, function(x) x$aic))
          if (best_pairdown_aic < current_aic) {
            best_pairdown_formula <- results[[which.min(sapply(results, function(x) x$aic))]]$formula_str
            current_formula <- stats::as.formula(best_pairdown_formula)
            current_aic <- best_pairdown_aic
            cat("New best pairdown model:", best_pairdown_formula, "AIC:", best_pairdown_aic, "\n")
          } else {
            cat("No further improvement, final model selected.\n")
            break
          }
        } else {
          cat("No improvement, stopping pair-down phase.\n")
          break
        }
      }

      # end_time <- Sys.time()  # to track run time per study
      run_time <- tictoc::toc()  # to track run time per study

      # Final pairdown model
      final_model <- ordinal(current_formula, data = input[[s]])
      RESULTS[[s]]$final_model <- final_model
      RESULTS[[s]]$final_form <- current_formula

      # Useful info for meta-analysis
      RESULTS[[s]]$odds_ratios <- exp(final_model$est)

      # Retrieve beta estimates
      vcov.ordinal <- final_model$COV  # extract variance-covariance matrix
      # Build matrix to summarize model coefficients
      RESULTS[[s]]$coefs <- matrix(NA, length(final_model$est), 2,
                      dimnames = list(names(final_model$est),
                                      c("Estimate", "Std. Error")))
      RESULTS[[s]]$coefs[,1] <- final_model$est  # beta estimates
      RESULTS[[s]]$coefs[,2] <- sqrt(diag(vcov.ordinal))  # standard error
      # results[[s]]$run_time <- end_time-start_time  # total run time for study
      RESULTS[[s]]$run_time <- run_time$callback_msg  # total run time for study

    }  # end for loop of studies

  } else {

    # Non-parallel format for non-Unix systems
    for (s in 1:length(studies)) {

      # start_time <- Sys.time()  # for run time of each study
      tictoc::tic()  # to track run time per study

      # Initial formula (intercept only)
      initial_formula <- stats::as.formula(paste(response_var, "~ 1"))
      best_formula <- initial_formula

      # Initial model
      best_model <- ordinal(best_formula, data = input[[s]])  # ordinal probit regression
      best_aic <- best_model$AIC

      ############################################
      # BUILD-UP PHASE
      ############################################
      fit_model <- function(formula, data) {
        tryCatch({
          model <- ordinal(formula, data = data)
          aic <- model$AIC
          return(list(model = model, aic = aic, error = NULL))
        }, error = function(e) {
          return(list(model = NULL, aic = Inf, error = e$message))
        })
      }

      # Remove predictors with zero variance
      ZV <- c()
      for (var in predictors) {
        if (length(unique(input[[s]][[var]])) == 1 && is.na(unique(input[[s]][[var]]))) {
          ZV[var] <- var
        }
      }
      predictors <- setdiff(predictors, ZV)

      for (i in 1:length(predictors)) {

        current_predictors <- all.vars(best_formula)[-1]
        remaining_predictors <- setdiff(predictors, current_predictors)

        candidate_models <- list()  # for storing model and their AIC

        # Build models iteratively based on current best model
        for (predictor in remaining_predictors) {
          new_formula <- stats::update(best_formula, paste0("~ . +", predictor))
          # Note that the width.cutoff and collapse are EXTREMELY important or deparse will split your string into two lines by default
          formula_str <- paste0(deparse(new_formula, width.cutoff = 500), collapse = "")
          fit <- fit_model(new_formula, input[[s]])

          if (!is.null(fit$error)) {
            # Log error if predictor is bad
            cat("Error for predictor", predictor, ":", fit$error, "\n")
            next  # Skip to next iteration if error
          }

          candidate_models[[formula_str]] <- fit$aic # store the new model and AIC in cadidate model list
          print(paste("Testing build-up formula:", formula_str, "with AIC:", fit$aic)) # helpful output
        }

        # ID the best model
        if (length(candidate_models) > 0) {
          best_candidate <- which.min(unlist(candidate_models)) # ID model with lowest AIC
          best_candidate_formula <- names(candidate_models)[best_candidate] # Get formula of best candidate model (string)
          best_candidate_aic <- unlist(candidate_models[best_candidate_formula])
          # Update to the new best model if it improves AIC
          if (best_candidate_aic < best_aic) { # If the AIC of the best candidate model is lower than the current best model's AIC, it replaces the model with the new model
            best_formula <- stats::as.formula(best_candidate_formula) # Converts the model formula string back into a regular formula object
            best_aic <- best_candidate_aic # Updates best_aic to the AIC of the new best model
            print(paste("New best model:", best_candidate_formula, "AIC:", best_aic)) # helpful output
          } else {
            print("No further improvement, stopping build-up.")
            break
          }
        } else {
          print("No more predictors to test, stopping build-up.")
          break
        }
      }

      # Final model
      best_model <- ordinal(best_formula, data = input[[s]])

      ############################################
      # PAIR-DOWN PHASE
      ############################################
      current_formula <- best_formula  # start with best FORMULA from build-up phase
      # Iteratively remove predictors
      repeat {
        predictors_in_model <- all.vars(current_formula)[-1]  # get all predictors currently in the best model
        candidate_models <- list() # list to store models and AIC
        current_model <- ordinal(current_formula, data = input[[s]])
        current_AIC <- current_model$AIC

        for (predictor in predictors_in_model) { # Loop through each predictors to test their removal
          if(length(predictors_in_model <= 1)) {
            pairdown_formula <- stats::as.formula(paste(response_var, "~", 1))
          } else {
            pairdown_formula <- stats::as.formula(paste(response_var, "~", paste(setdiff(predictors_in_model, predictor), collapse = "+")))
          }

          if (length(all.vars(pairdown_formula)[-1]) == 0) { # check if model is empty (no predictors)
            next  # skip iteration if no predictors are left
          }

          fit <- fit_model(pairdown_formula, input[[s]])
          formula_str <- deparse(pairdown_formula, width.cutoff = 500)
          candidate_models[formula_str] <- fit$aic # Store AIC and formula of pairdown model
          print(paste("Testing pairdown formula:", deparse(pairdown_formula), "with AIC:", fit$aic)) # Helpful output
        }

        # See if any pairdown model is better than current best model
        if (length(candidate_models) > 0) {
          best_pairdown_aic <- min(sapply(candidate_models, identity)) # find smallest AIC among pairdown models

          if (best_pairdown_aic < current_aic) { # If a pairdown model has lower AIC, update current best model
            best_pairdown_formula <- names(candidate_models)[which.min(sapply(candidate_models, identity))]
            current_formula <- stats::as.formula(best_pairdown_formula)
            current_aic <- best_pairdown_aic
            print(paste("New best pairdown model:", best_pairdown_formula, "AIC:", best_pairdown_aic)) # Output new model's formula and AIC
          } else {
            print("No further improvement, final model selected.")
            break  # Exit loop if no improvement
          }
        } else {
          print("No improvement, stopping pair-down phase.") # If no pairdown models found with lower AIC, then stops the process
          break
        }
      }

      # end_time <- Sys.time()  # for run time of each study
      run_time <- tictoc::toc()  # to track run time per study

      # Save selected model
      final_model <- ordinal(current_formula, data = input[[s]])
      RESULTS[[s]]$final_model <- final_model
      RESULTS[[s]]$final_form <- current_formula  # selected model formula
      RESULTS[[s]]$odds_ratios <- exp(final_model$est) # Calculate OR

      # Retrieve beta estimates
      vcov.ordinal <- final_model$COV  # extract variance-covariance matrix
      # Build matrix to summarize model coefficients
      RESULTS[[s]]$coefs <- matrix(NA, length(final_model$est), 2,
                                   dimnames = list(names(final_model$est),
                                                   c("Estimate", "Std. Error")))

      RESULTS[[s]]$coefs[,1] <- final_model$est  # beta estimates
      RESULTS[[s]]$coefs[,2] <- sqrt(diag(vcov.ordinal))  # standard error
      # results[[s]]$run_time <- end_time-start_time  # total run time for study
      RESULTS[[s]]$run_time <- run_time$callback_msg  # total run time for study

    }  # end for loop of studies
  }

  # Rename results list according to study names
  names(RESULTS) <- studies

  return(RESULTS)

}

