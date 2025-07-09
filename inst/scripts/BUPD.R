############################################
# INITIALIZE
############################################
library(parallel)
# library(MASS)
source("R/ordinal.R")
numCores <- detectCores() - 1

predictors <- predictors_full  # save predictors

# Initial empty model
initial_formula <- as.formula(paste(response_var, "~ 1"))
best_formula <- initial_formula

best_model <- ordinal(best_formula, data = input)
best_aic <- best_model$AIC
# best_model <- polr(best_formula, data = input, Hess = TRUE, method = "probit")
# best_aic <- AIC(best_model)
name_prefix <- gsub(" ", "_", name) # add underscores

############################################
# BUILD-UP PHASE
############################################
# Function to fit and evaluate models
fit_model <- function(formula, data) {
  tryCatch({
    # model <- polr(formula, data = data, Hess = TRUE, method = "probit")
    # aic <- AIC(model)
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
for (var in predictors) {
  if (length(unique(input[[var]])) == 1 && is.na(unique(input[[var]]))) {
    ZV[var] <- var
  }
}
predictors <- setdiff(predictors, ZV)

# Build-up phase
for (i in 1:length(predictors)) {
    current_predictors <- all.vars(best_formula)[-1]
    remaining_predictors <- setdiff(predictors, current_predictors)

    # Use mclapply to parallelize the model fitting with forked processes
    results <- mclapply(remaining_predictors, function(predictor) {
        new_formula <- update(best_formula, paste(". ~ . +", predictor))
        fit <- fit_model(new_formula, input)

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
            best_formula <- as.formula(best_candidate_formula)
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
# best_model <- polr(best_formula, data = input, Hess = TRUE, method = "probit")
best_model <- ordinal(best_formula, data = input)
best_aic <- best_model$AIC

############################################
# PAIR-DOWN PHASE
############################################
current_formula <- best_formula  # start with best FORMULA from build-up phase
# current_aic <- AIC(best_model)
current_aic <- best_model$AIC

# >>>>> WORKS TO THIS POINT (~66 preds, ~ 12 hours) <<<<<
# 08/09/2024: Seems to be working for pairdown; outputting formulas as it tests. Waiting for it to finish. If need to, can save outut from buildup for later testing.

repeat {
    predictors_in_model <- all.vars(current_formula)[-1]

    # Use mclapply to parallelize the pair-down model fitting with forked processes
    results <- mclapply(predictors_in_model, function(predictor) {
      if(length(predictors_in_model <= 1)) {
        pairdown_formula <- as.formula(paste(response_var, "~", 1))
      } else {
        pairdown_formula <- as.formula(paste(response_var, "~", paste(setdiff(predictors_in_model, predictor), collapse = "+")))
      }
      if (length(all.vars(pairdown_formula)[-1]) == 0) {
        return(NULL)
      }

      fit <- fit_model(pairdown_formula, input)
      formula_str <- paste(deparse(pairdown_formula, width.cutoff = 500), collapse = "")
      list(formula_str = formula_str, aic = fit$aic)
    }, mc.cores = numCores)

    # Filter out results with errors
    results <- Filter(function(x) !is.null(x$formula_str), results)

    if (length(results) > 0) {
      best_pairdown_aic <- min(sapply(results, function(x) x$aic))
      if (best_pairdown_aic < current_aic) {
          best_pairdown_formula <- results[[which.min(sapply(results, function(x) x$aic))]]$formula_str
          current_formula <- as.formula(best_pairdown_formula)
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

# Final pairdown model
# final_model <- polr(current_formula, data = input, Hess = TRUE, method = "probit")
final_model <- ordinal(current_formula, data = input)
# final_form <- formula(final_model)
final_form <- current_formula

# Useful info for meta-analysis
# coeff <- coef(final_model)
coeff <- final_model$est
odds_ratios <- exp(coeff)

# `assign` to new variables based on name prefix chosen
assign(paste0(name_prefix, "_final_model"), final_model)
assign(paste0(name_prefix, "_final_form"), final_form)
assign(paste0(name_prefix, "_odds_ratios"), odds_ratios)

# Save for later
output_directory <- "data"
saveRDS(get(paste0(name_prefix, "_final_model")), file = file.path(output_directory, paste0(name_prefix, "_final_model.rds")))
saveRDS(get(paste0(name_prefix, "_final_form")), file = file.path(output_directory, paste0(name_prefix, "_final_form.rds")))
saveRDS(get(paste0(name_prefix, "_odds_ratios")), file = file.path(output_directory, paste0(name_prefix, "_odds_ratios.rds")))
