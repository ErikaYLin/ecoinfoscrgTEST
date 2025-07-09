# QUANTILE REGRESSION

load("data/data2_standardize_impute.rda")

# Scale response to be on [0,1] scale for glmmTMB
resp <- data2$state
resp$Response <- as.integer(resp$Response)
# resp$Response <- (resp$Response-1)/(3-1)  # scale response for ordbeta family
# resp$Response <- as.factor(resp$Response)
# resp$Response <- ordered(resp$Response)

# Define the response variable
response_var <- "Response"

pred <- data2$predictors  # retrieve predictor data
pred <- pred %>%
  mutate(across(c(OBJECTID, ID, bmpCountv5, n, distance, X, Y), as.character))

# Filter by study
resp_choc <- resp %>% filter(study == "choc")
resp_pens <- resp %>% filter(study == "pens")
resp_tampa <- resp %>% filter(study == "tampa")
resp_IRL <- resp %>% filter(study == "IRL")

pred_choc <- pred %>% filter(study == "choc")
pred_pens <- pred %>% filter(study == "pens")
pred_tampa <- pred %>% filter(study == "tampa")
pred_IRL <- pred %>% filter(study == "IRL")

# List all predictors (AKA column names of known predictors)
numeric_pred <- pred %>%
  select_if(is.numeric)
factor_pred <- pred %>%
  select_if(is.factor)
# factor_pred <- factor_pred %>%
# mutate(across(all_of(colnames(factor_pred)), as.character)) %>%
# mutate(across(all_of(colnames(factor_pred)), as.numeric))
predictors <- colnames(cbind(factor_pred, numeric_pred))
# predictors <- predictors[-5]  # remove SandSpit

# data3 is pre-standardized data

# Dataset
full_data <- cbind(resp, pred)
choc_data <- cbind(resp_choc, pred_choc)
pens_data <- cbind(resp_pens, pred_pens)
tampa_data <- cbind(resp_tampa, pred_tampa)
IRL_data <- cbind(resp_IRL, pred_IRL)

# Test formula (5 categorical vars, 5 numeric vars)
FORMULA <- as.formula(paste0(paste(response_var, "~"), paste(predictors[c(1:5,85:90)], collapse = "+")))
FORMULA_full <- as.formula(paste("Response ~", paste(predictors, collapse = "+")))

## Ordinal probit regression via custom function

source("R/ordinal.R")

mod <- ordinal(formula = FORMULA_full, data = full_data)


## Ordered beta regression with `brms` wrapper (`ordbetareg`)

library(ordbetareg)

mod <- ordbetareg::ordbetareg(formula = Response ~ 1, data = choc_data, true_bounds = c(1,3))


# # Test quantile regression on full dataset
# mod <- quantreg::rq(FORMULA, data = full_data, tau = 1/3, method = "br")
## Error in rq.fit.br(x, y, tau = tau, ...) : Singular design matrix
## likely because of factored response, can't include a probit link

# for linear regression, truncate response to 1-3 (0-100%)

## Ordered beta regression with glmmTMB

library(glmmTMB)

# Unlike the implementation in the `ordbeta` package (wrapper for `brms`), this family will not
# automatically scale the data. If your response variable is defined on the closed interval [a,b],
# transform it to [0,1] via y_scaled <- (y-a)/(b-a)

# Test formula (5 categorical vars, 5 numeric vars)
FORMULA2 <- as.formula(paste("Response ~", paste(predictors[c(1:5,85:90)], collapse = "+"), paste0(" + (1|study)")))

mod <- glmmTMB::glmmTMB(formula = Response ~ angle + (1|study), data = full_data,
                        family = glmmTMB::ordbeta(link = "probit"),
                        start = rep(0.5, 10))



## Ordinal probit regression with CLM

library(ordinal)

# Ordinal cumulative link model w/ probit link
mod <- ordinal::clm(formula = FORMULA, data = cbind(resp,pred),
                    link = "probit",
                    start = c(1/3, 2/3, rep(2, 5)),
                    control = list(method = "design",  # try adjusting tolerances to get convergence
                                   lower = c((1/3 - .Machine$double.eps), (2/3 - .Machine$double.eps), rep(-Inf, 5)),
                                   upper = c((1/3 + .Machine$double.eps), (2/3 - .Machine$double.eps), rep(+Inf, 5))))
mod$control$method <- "optim"
mod2 <- ordinal::clm.fit(mod)
## Can't specify threshold or pass arguments to `optim`


# CALCULATE LOG-LIKELIHOOD
## get model matrix multiplied by betas to get linear effect
## calculate the pdf under the link function and sum to get log-likelihoods

# Prob density of normal distribution
NLL_norm <- function(p) {

  # Get model matrix (design matrix) for linear effect
  mod.mat <- mod$X[,-1] %*% mod2$beta
  mu <- pnorm(mod.mat)

  -sum(dnorm(x, mean = mu, sd = ..., log = TRUE)) # negative log-likelihood
}


# Try to optimize and use fixed parameters (TRY AFTER FIXING LOGLIK ISSUE)
fn <- function(PRED, verbose=FALSE)
{
  PRED <- c(mod$alpha[1], mod$alpha[2], mod$beta[1], mod$beta[2], mod$beta[3], mod$beta[4],
            mod$beta[5], mod$beta[6])
#   FORMULA <- as.formula(paste(response_var, "~",
#                               paste(PRED[3], PRED[4]*predictors[1], PRED[5]*predictors[2],
#                                     PRED[6]*predictors[3], PRED[7]*predictors[4], PRED[8]*predictors[5],
#                                     pred[9]*predictors[6], PRED[10]*predictors[7],
#                                     PRED[11]*predictors[8], PRED[11]*predictors[9], PRED[12]*predictors[10], collapse = "+")))

  FORMULA <- as.formula(paste(response_var, "~ Beach + WideBeach + PublicRamp + canal + SAV"))
  mod <- ordinal::clm(formula = FORMULA, data = cbind(resp,pred),
                      link = "probit",
                      control = list(method = "optim"))
                      # start = c(1/3, 2/3, rep(0, 5)),
                      # control = list(method = "design",  # try adjusting tolerances to get convergence
                      #                lower = c((1/3 - .Machine$double.eps), (2/3 - .Machine$double.eps), rep(-Inf, 5)),
                      #                upper = c((1/3 + .Machine$double.eps), (2/3 - .Machine$double.eps), rep(+Inf, 5))))
  # mod$control$method <- "optim"
  # mod2 <- ordinal::clm.fit(mod)
  NLL <- -mod$logLik
  if(verbose)
  { return(mod) }  # can give option to return model or neg loglik
  else
  { return(NLL) }
}

start <- c(1/3, 2/3, rep(0, 5))
OPT <- optim(par = start, fn, method = "L-BFGS-B",
             lower = c((1/3 - .Machine$double.eps), (2/3 - .Machine$double.eps), rep(-Inf, 5)),
             upper = c((1/3 + .Machine$double.eps), (2/3 + .Machine$double.eps), rep(+Inf, 5)))
OPT  # most likely parameters

FIT <- fn(OPT$minimum,verbose=TRUE)
FIT  # model with best period

## lower and upper bounds as arrays in order of parameters, betas would be (-inf,inf), thetas would be (1/3 - epsilon) and (2/3 - epsilon), starting values as 1/3, 2/3

# # Ordinal probit regression (CAN pass arguments to `optim`)
# mod <- MASS::polr(FORMULA, data = full_data, Hess = TRUE, method = "probit",
#                   lower = (1/3), upper = (2/3))
## In optim(s0, fmin, gmin, method = "BFGS", ...): bounds can only be used with method L-BFGS-B (or Brent)

# Try feeding rq results from optim into polr model
start <- c(1/3, 2/3)  # ???
thresholds <- optim(start, )



mod <- MASS::polr(FORMULA, data = full_data, Hess = TRUE, method = "probit")






