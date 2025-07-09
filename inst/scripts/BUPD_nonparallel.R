############################################
# INITIALIZE
############################################
# library(nnet)
# library(lme4)
# library(MASS)
# library(ordinal)
# library(glmmTMB)
# library(ordbetareg)
source("R/ordinal.R")

predictors <- predictors_full  # save predictors

# Initial empty model
# Note that `best_formula` starts with `initial_formula` as baseline
# initial_formula <- as.formula(paste(response_var, "~ 1 + (1|study)"))
# initial_formula <- as.formula(paste(response_var, "~ (1|study)"))
initial_formula <- as.formula(paste(response_var, "~ 1"))
best_formula <- initial_formula

# Save number of predictors
# k <- length(strsplit(as.character(best_formula[3]), split = " + ", fixed = TRUE)) - 1
# k <- length(all.vars(best_formula)) - 1  # this one works

# best_model <- multinom(best_formula, data = data, MaxNWts = 5000) #nnet
# best_model <- lme4::lmer(best_formula, data = input) # lme4
# best_model <- polr(best_formula, data = input, method = "probit") # MASS
# best_model <- polr(best_formula, data = input, Hess=TRUE, method="probit")
# best_model <- ordinal::clmm(best_formula, data = input, Hess = TRUE,
#                            link = "probit", threshold = "equidistant")
# best_model <- glmmTMB::glmmTMB(best_formula, data = input,
#                                family = glmmTMB::ordbeta(link = "probit"))
# best_model <- ordbetareg::ordbetareg(best_formula, data = input, true_bounds = c(1,3))
best_model <- ordinal(best_formula, data = input)

# best_mod_dmatrix <- ordinal::clm(formula = best_formula, data = input,
#                                  link = "probit", start = c(1/3, 2/3, rep(0, k)),
#                                  control = list(method = "design",  # try adjusting tolerances to get convergence
#                                                 lower = c(1/3 - .Machine$double.eps, 2/3 - .Machine$double.eps, rep(-Inf, k)),
#                                                 upper = c(1/3 + .Machine$double.eps, 2/3 - .Machine$double.eps, rep(+Inf, k))))
# best_mod_dmatrix$control$method <- "optim"
# best_model <- clm.fit(best_mod_dmatrix)

# best_model <- ordinal::clm(best_formula, data = input, Hess = TRUE,
#                            method = "optim", link = "probit", threshold = "symmetric")
# control_settings <- clmm.control(maxIter = 2000)  # clmm stuff
# best_model <- clmm(best_formula, data = input, control = control_settings, link = "probit") # ordinal
# best_model <- glmer(best_formula, data = input, family=binomial(link="probit"))
# best_model <- glm(best_formula, data = input, family=binomial(link="probit"))
# glmer(Response ~ (1 | study) + Beach_1, data = input, family=binomial(link="probit"))

best_aic <- best_model$AIC
# best_aic <- AIC(best_model)
# best_aic <- (2*length(best_model$coefficients)) + (2*best_model$logLik)
name_prefix <- gsub(" ", "_", name) # add underscores
############################################
# BUILD-UP PHASE
############################################
# Function to fit and evaluate models (similar to Chris' only with a multinomial logistic regression for this specific dataset)
fit_model <- function(formula, data) {
  # model <- lme4::lmer(formula, data = data)
  # model <- clmm(formula, data = data, link = "probit")
  # model <- polr(formula, data = data, method = "probit")
  tryCatch({ # NEW LINE
    # model <- lme4::lmer(formula, data = data)
    # model <- polr(formula, data = data, Hess=TRUE, method="probit")
    # model <- ordinal::clm(formula, data = data, Hess = TRUE, method = "optim", link = "probit")
    # k <- length(strsplit(as.character(formula[3]), split = " + ", fixed = TRUE))
    # model <- ordbetareg::ordbetareg(formula, data = data, true_bounds = c(1,3))
    model <- ordinal(formula, data = data)

    # model <- glmmTMB::glmmTMB(formula, data = data, family = glmmTMB::ordbeta(link = "probit"))

    # k <- length(all.vars(formula)) - 1
    # mod_dmatrix <- ordinal::clm(formula, data = data,
    #                             link = "probit", start = c(1/3, 2/3, rep(0, k)),
    #                             control = list(method = "design",  # try adjusting tolerances to get convergence
    #                                            lower = c(1/3 - .Machine$double.eps, 2/3 - .Machine$double.eps, rep(-Inf, k)),
    #                                            upper = c(1/3 + .Machine$double.eps, 2/3 - .Machine$double.eps, rep(+Inf, k))))
    # mod_dmatrix$control$method <- "optim"
    # model <- clm.fit(mod_dmatrix)

    # model <- ordinal::clmm(formula, data = input, Hess = TRUE,
    #                        link = "probit", threshold = "equidistant")
    aic <- model$AIC
    # aic <- AIC(model) # NEW LINE
    # aic <- (2*length(model$coefficients)) + (2*model$logLik)
    return(list(model = model, aic = aic, error = NULL))
  }, error = function(e) { # NEW LINE
    return(list(model = NULL, aic = Inf, error = e$message))
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

# Similar to Chris' only with a multinomial logistic regression (which can be switched out)
for (i in 1:length(predictors)) {
    current_predictors <- all.vars(best_formula)[-1]
    remaining_predictors <- setdiff(predictors, current_predictors)
    # k <- length(current_predictors)  # number of parameters for each model to test

    candidate_models <- list()  # for storing model and their AIC
    # Build models iteratively based on current best model
    for (predictor in remaining_predictors) {
        # new_formula <- update(best_formula, paste(". ~ . +", predictor))
        new_formula <- update(best_formula, paste0("~ . +", predictor))
        formula_str <- paste0(deparse(new_formula, width.cutoff = 500), collapse = "") # converts model formula to string. Note that the width.cutoff and collapse are EXTREMELY important or deparse will split your string into two lines by default
        fit <- fit_model(new_formula, input)
        if (!is.null(fit$error)) { # NEW LINE
        # Log error if predictor is bad # NEW LINE
          cat("Error for predictor", predictor, ":", fit$error, "\n") # NEW LINE
          next  # Skip to next iteration if error # NEW LINE
        } # NEW LINE
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
            best_formula <- as.formula(best_candidate_formula) # Converts the model formula string back into a regular formula object
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

#####
# Final model
# best_model <- multinom(best_formula, data = data, MaxNWts = 5000, trace = FALSE)
# best_model <- lme4::lmer(best_formula, data = input)
# best_model <- clmm(best_formula, data = input, link = "probit")
# best_model <- polr(best_formula, data = input, Hess=TRUE, method="probit")
# best_model <- ordinal::clm(best_formula, data = input, Hess = TRUE,
#                            method = "optim", link = "probit", threshold = "symmetric")
# best_model <- ordinal::clmm(best_formula, data = input, Hess = TRUE,
#                             link = "probit", threshold = "equidistant")
# best_model <- ordbetareg::ordbetareg(best_formula, data = input, true_bounds = c(1,3))
best_model <- ordinal(best_formula, data = input)

# best_model <- glmmTMB::glmmTMB(best_formula, data = input,
#                                family = glmmTMB::ordbeta(link = "probit"))

# k <- length(all.vars(best_formula)) - 1
# best_mod_dmatrix <- ordinal::clm(formula = best_formula, data = input,
#                                  link = "probit", start = c(1/3, 2/3, rep(0, k)),
#                                  control = list(method = "design",  # try adjusting tolerances to get convergence
#                                                 lower = c(1/3 - .Machine$double.eps, 2/3 - .Machine$double.eps, rep(-Inf, k)),
#                                                 upper = c(1/3 + .Machine$double.eps, 2/3 - .Machine$double.eps, rep(+Inf, k))))
# best_mod_dmatrix$control$method <- "optim"
# best_model <- clm.fit(best_mod_dmatrix)

# summary(best_model)
############################################
# PAIR-DOWN PHASE
############################################
current_formula <- best_formula  # start with best FORMULA from build-up phase
# Iteratively remove predictors
repeat {
  predictors_in_model <- all.vars(current_formula)[-1]  # get all predictors currently in the best model
  candidate_models <- list() # list to store models and AIC
  # current_model <- lme4::lmer(current_formula, data = input)
  # current_model <- clmm(current_formula, data = input, link = "probit")
  # LM_PD <- lm(current_formula, data = input)
  # current_model <- polr(current_formula, data = input, start = c(LM_PD$coefficients[-1], qnorm(c(1/3,2/3))),
  #                       Hess=TRUE, method="probit")
  # current_model <- polr(current_formula, data = input, Hess=TRUE, method="probit")
  # current_model <- ordinal::clm(current_formula, data = input, Hess = TRUE,
  #                               method = "optim", link = "probit", threshold = "symmetric")
  # current_model <- ordinal::clmm(current_formula, data = input, Hess = TRUE,
  #                               link = "probit", threshold = "equidistant")
  # current_model <- ordbetareg::ordbetareg(current_formula, data = input, true_bounds = c(1,3))
  current_model <- ordinal(current_formula, data = input)

  # current_model <- glmmTMB::glmmTMB(current_formula, data = input,
  #                                   family = glmmTMB::ordbeta(link = "probit"))

  # k <- length(all.vars(current_formula)) - 1
  # current_mod_dmatrix <- ordinal::clm(current_formula, data = input,
  #                             link = "probit", start = c(1/3, 2/3, rep(0, k)),
  #                             control = list(method = "design",  # try adjusting tolerances to get convergence
  #                                            lower = c(1/3 - .Machine$double.eps, 2/3 - .Machine$double.eps, rep(-Inf, k)),
  #                                            upper = c(1/3 + .Machine$double.eps, 2/3 - .Machine$double.eps, rep(+Inf, k))))
  # current_mod_dmatrix$control$method <- "optim"
  # current_model <- clm.fit(current_mod_dmatrix)

  # current_aic <- AIC(multinom(current_formula, data = data, MaxNWts = 5000, trace = FALSE)) # calculate AIC of current best model
  current_AIC <- current_model$AIC
  # current_aic <- AIC(current_model)
  # current_aic <- (2*length(current_model$coefficients)) + (2*current_model$logLik)

  for (predictor in predictors_in_model) { # Loop through each predictors to test their removal
    if(length(predictors_in_model <= 1)) {
      pairdown_formula <- as.formula(paste(response_var, "~", 1))
    } else {
      pairdown_formula <- as.formula(paste(response_var, "~", paste(setdiff(predictors_in_model, predictor), collapse = "+")))
    }
    # pairdown_formula <- as.formula(paste(response_var, "~", paste(setdiff(predictors_in_model, predictor), collapse = "+"))) # New formula without current predictor
    if (length(all.vars(pairdown_formula)[-1]) == 0) { # check if model is empty (no predictors)
      next  # skip iteration if no predictors are left
    }
    # pairdown_aic <- AIC(multinom(pairdown_formula, data = data, MaxNWts = 5000, trace = FALSE)) # fit pairdown model, get AIC
    fit <- fit_model(pairdown_formula, input)
    # formula_str <- paste(deparse(pairdown_formula, width.cutoff = 500), collapse = "") # store the name/string of the model properly
    formula_str <- deparse(pairdown_formula, width.cutoff = 500)
    candidate_models[formula_str] <- fit$aic # Store AIC and formula of pairdown model
    print(paste("Testing pairdown formula:", deparse(pairdown_formula), "with AIC:", fit$aic)) # Helpful output
  }

  # See if any pairdown model is better than current best model
  if (length(candidate_models) > 0) {
    best_pairdown_aic <- min(sapply(candidate_models, identity)) # find smallest AIC among pairdown models
    if (best_pairdown_aic < current_aic) { # If a pairdown model has lower AIC, update current best model
      best_pairdown_formula <- names(candidate_models)[which.min(sapply(candidate_models, identity))]
      current_formula <- as.formula(best_pairdown_formula)
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
# Final pairdown model
# final_model <- multinom(current_formula, data = data, MaxNWts = 5000, trace = FALSE)
# final_model <- lme4::lmer(current_formula, data = input)
# final_model <- clmm(current_formula, data = input, link = "probit")
# final_model <- polr(current_formula, data = input, Hess=TRUE, method="probit")
# final_model <- ordinal::clm(current_formula, data = input, Hess = TRUE,
#                             method = "optim", link = "probit", threshold = "symmetric")
# final_model <- ordinal::clmm(current_formula, data = input, Hess = TRUE,
#                            link = "probit", threshold = "equidistant")
# final_model <- ordbetareg::ordbetareg(current_formula, data = input, true_bounds = c(1,3))
final_model <- ordinal(current_formula, data = input)

# final_model <- glmmTMB::glmmTMB(current_formula, data = input,
#                                 family = glmmTMB::ordbeta(link = "probit"))

# k <- length(all.vars(current_formula)) - 1
# final_mod_dmatrix <- ordinal::clm(current_formula, data = input,
#                                  link = "probit", start = c(1/3, 2/3, rep(0, k)),
#                                  control = list(method = "design",  # try adjusting tolerances to get convergence
#                                                 lower = c(1/3 - .Machine$double.eps, 2/3 - .Machine$double.eps, rep(-Inf, k)),
#                                                 upper = c(1/3 + .Machine$double.eps, 2/3 - .Machine$double.eps, rep(+Inf, k))))
# final_mod_dmatrix$control$method <- "optim"
# final_model <- clm.fit(final_mod_dmatrix)

# summary(final_model)
# final formula
# final_form <- formula(final_model)  # final_model
final_form <- current_formula
# print(final_form)
# # Useful info for meta-analysis
# coeff <- coef(final_model) # grab coefficients
coeff <- final_model$est
# standard_err <- sqrt(diag(vcov(final_model))) # Calculate SE (method OK?)
# confidence_intervals <- confint(final_model, level = 0.95) # Calculate CI
odds_ratios <- exp(coeff) # Calculate OR
# `assign` to new variables based on name prefix chosen
assign(paste0(name_prefix, "_final_model"), final_model)
assign(paste0(name_prefix, "_final_form"), final_form)
assign(paste0(name_prefix, "_odds_ratios"), odds_ratios)
# save for later
output_directory <- "data" # output directory
saveRDS(get(paste0(name_prefix, "_final_model")), file = file.path(output_directory, paste0(name_prefix, "_final_model.rds")))
saveRDS(get(paste0(name_prefix, "_final_form")), file = file.path(output_directory, paste0(name_prefix, "_final_form.rds")))
saveRDS(get(paste0(name_prefix, "_odds_ratios")), file = file.path(output_directory, paste0(name_prefix, "_odds_ratios.rds")))
