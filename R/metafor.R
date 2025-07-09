# library(metafor)
# library(Matrix)

# source("ecoinfoscrg/R/getBetas.R")
# source("ecoinfoscrg/R/varCov.R")
# source("R/getBetas.R")
# source("R/varCov.R")

##############################
# Meta-Analytic Regression
##############################
# Note: Testing with betas but need effect sizes (e.g. odd ratios) from BUPD.R model output
#' @export
meta_regression <- function(data, varCov) {

  # Retrieve necessary components
  combined_betas_only <- varCov$combined_betas_only
  VAR <- varCov$VAR
  studies <- data$state$study

  # Flatten beta coefficients
  betas_vector <- as.vector(t(combined_betas_only))

  # Get study vector corresponding to each beta (required for study random effect)
  # study_labels <- c("choc", "pens", "tampa", "IRL") # will need to change labels as more data is available
  study_vector <- rep(studies, each = ncol(combined_betas_only))

  # specify input data
  beta = betas_vector
  study = factor(study_vector)
  predictor = rep(colnames(combined_betas_only), times = length(studies))
  # variance <- as.vector(t(adjusted_VAR))
  variance <- as.vector(t(VAR))
  # variance <- adjusted_cov_matrix
  # variance <- diag(adjusted_cov_matrix)

  # FUNCTION FOR CALLING rma.mv (# `meta_regression` will run with default options )
  meta_model <- function(
    beta = NULL,
    variance = NULL,
    predictor = NULL,
    study = NULL,
    method = "REML",
    struct = "UN",
    verbose = FALSE
  ) {
    if (is.null(beta)) beta <- get("beta", envir = .GlobalEnv) # default option
    if (is.null(variance)) variance <- get("variance", envir = .GlobalEnv) # default option
    if (is.null(predictor)) predictor <- get("predictor", envir = .GlobalEnv) # default option
    if (is.null(study)) study <- get("study", envir = .GlobalEnv) # default option

    result <- metafor::rma.mv(
      yi = beta,          # Vector of all beta coefficients
      V = variance,       # Vector of sampling variances or variance-covariance matrix
      method = method,    # Method for the meta-analysis
      # random = ~ 1 | study,  # Random effects for studies (contact authors to understand if they had random effects or separate regressions)
      # struct = struct,    # Structure of the variance-covariance matrix
      mods = ~ predictor, # Including predictors as fixed effects
      verbose = verbose   # Verbosity of the output
    )

    return(result)
  }

  unscaled.meta <- meta_regression()

  return(unscaled.meta)

}





#########################
# FOR FURTHER TESTING
#########################

# saveRDS(unscaled.meta, file = "../output/unscaled_meta.rds")  # file.path("output", "unscaled_meta.rds"))
# saveRDS(unscaled.meta, file = "../output/unscaled_meta3.rds")  # file.path("output", "unscaled_meta.rds"))
# unscaled.meta <- readRDS("../output/unscaled_meta.rds")


# output <- meta_regression()
# print(output)

# png("forest_plot.png")
# forest(output)  # effect sizes and CIs
# dev.off()

# png("funnel_plot.png")
# funnel(output) # pub bias
# dev.off()

################################
# SCALE BETAS TO REFERENCE STUDY
################################

# # Tampa as reference
# scale_reference <- 3
# # Could set it as index with matching study name for including as arguments in functions
# # scale_reference <- which(combined_betas$study=="tampa")

# # Function to scale betas
# # scale_betas function scales each beta coef by theta/GM(theta)
# # covs are scaled by theta%*%t(theta)/(GM(theta)^2), but due to r code structure, we use thetas%*%thetas/(GM(theta)^2)
# # thetas is a vector of scaling parameters
# # keep is a vector of column indices to hold constant (not scale), if any
# scale.betas <- function (thetas, adjust = TRUE) {
#   results <- list()
#   scaled_betas <- combined_betas_only
#   scaled_cov <- cov_matrix
#   scaled_var <- VAR
#   gm_thetas <- prod(thetas)^(1/length(thetas))
#   scalars.normalized <- thetas/gm_thetas
#   for (i in 1:length(thetas)){
#     scaled_betas[i, ] <- scaled_betas[i, ]*scalars.normalized[i]
#   }
#   # normalize it and then square
#   # need to fix in order to accommodate for per-site cov matrices
#   # this can be easily implemented by taking each normalized scalar (theta/gm_thetas)
#   # and then squaring them to obtain the cov matrix scalar for each site
#   cov.scalar <- scalars.normalized^2
#   # cov.scalar <- prod(scalars.normalized^2)
#
#   # if (adjust) {
#   #   scaled_cov_adjusted <- list()
#   #   scaled_cov_new <- list()
#   #   for (i in 1:length(scaled_cov)) {
#   #     cov.names <- dimnames(scaled_cov[[i]])  # save dimnames
#   #
#   #     # Scale covariance matrices
#   #     scaled_cov_new[[i]] <- cov.scalar[i]*scaled_cov[[i]]
#   #
#   #     ####### cov matrix adjustments
#   #     eigen_decomp <- eigen(scaled_cov_new[[i]])
#   #     eigenvalues <- eigen_decomp$values
#   #     eigenvectors <- eigen_decomp$vectors
#   #
#   #     # Find the smallest positive eigenvalue
#   #     smallest_eigenvalue <- min(eigenvalues[eigenvalues > 0])
#   #
#   #     # # Define the maximum allowed variance
#   #     # max_allowed_variance <- sqrt(1 / .Machine$double.eps) * smallest_eigenvalue
#   #
#   #     # Adjust eigenvalues
#   #     adjusted_eigenvalues <- pmax(eigenvalues, (.Machine$double.eps)^(1/3))
#   #     adjusted_eigenvalues <- pmin(adjusted_eigenvalues, (.Machine$double.eps)^(-1/3))
#   #     # adjusted_eigenvalues <- pmin(eigenvalues, max_allowed_variance)
#   #
#   #     # adjusted_eigenvalues <- pmax(adjusted_eigenvalues, (.Machine$double.eps)^(1/3))
#   #     # adjusted_eigenvalues <- pmin(adjusted_eigenvalues, (.Machine$double.eps)^(-1/3))
#   #     scaled_cov_adjusted[[i]] <- eigenvectors %*% diag(adjusted_eigenvalues) %*% t(eigenvectors)
#   #
#   #     # another check
#   #     smallest_eigenvalue <- min(adjusted_eigenvalues[adjusted_eigenvalues > 0])
#   #     # smallest_eigenvalue
#   #
#   #     # Set very large variances to very large value
#   #     large_value <- (.Machine$double.eps)^(-1/3)
#   #     diag(scaled_cov_adjusted[[i]])[which(diag(scaled_cov_adjusted[[i]])>large_value)] <- large_value
#   #
#   #     #####
#   #     # Recommendation: Add a small jitter to the diagonal to ensure positive definiteness
#   #     epsilon <- 1e-6
#   #     scaled_cov_adjusted[[i]] <- scaled_cov_adjusted[[i]] + diag(epsilon, nrow(scaled_cov_adjusted[[i]]))
#   #
#   #     dimnames(scaled_cov_adjusted[[i]]) <- cov.names
#   #   }
#   #   #####
#   #   #######
#   #
#   #   # VARIANCES
#   #   scaled_var_adjusted <- scaled_var
#   #   scaled_var_new <- scaled_var
#   #   for (i in 1:nrow(scaled_var)) {
#   #
#   #     # Scale variances
#   #     scaled_var_new[i,] <- cov.scalar[i]*scaled_var[i,]
#   #
#   #     # Adjustments to variance
#   #     scaled_var_adjusted[i,] <- pmax(scaled_var_new[i,], (.Machine$double.eps)^(1/3))
#   #     scaled_var_adjusted[i,] <- pmin(scaled_var_adjusted[i,], (.Machine$double.eps)^(-1/3))
#   #
#   #   }
#   #
#   #   # DEBUG <<- list(thetas = thetas, GM = gm_thetas, scalars = scalars.normalized,
#   #   #                beta = scaled_betas, variance = scaled_var_adjusted)
#   #
#   #   #####
#   #   #######
#   #   results[[1]] <- scaled_betas
#   #   results[[2]] <- scaled_cov_adjusted
#   #   results[[3]] <- scaled_var_adjusted
#   #
#   # } else {
#
#   if (adjust) {
#
#     scaled_cov_new <- list()
#     scaled_cov_adjusted <- list()
#     for (i in 1:length(scaled_cov)) {
#       cov.names <- dimnames(scaled_cov[[i]])  # save dimnames
#
#       # Scale covariance matrices
#       scaled_cov_new[[i]] <- cov.scalar[i]*scaled_cov[[i]]
#       # Clean each covariance matrix
#       scaled_cov_adjusted[[i]] <- clean_matrix(scaled_cov_new[[i]])
#     }
#
#     # Grab variances from each matrix
#     scaled_var_adjusted <- scaled_var  # simulate structure of combined betas dataframe
#     for (i in 1:nrow(scaled_var_adjusted)) {
#       for (j in 1:ncol(scaled_var_adjusted)) {
#
#         # Add variances to corresponding predictor and study
#         # VAR[i,j] <- diag(cov_matrix[[i]])[which(names(diag(cov_matrix[[i]])) == colnames(VAR)[j])]
#         scaled_var_adjusted[i,j] <- diag(scaled_cov_adjusted[[i]])[which(names(diag(scaled_cov_adjusted[[i]])) == colnames(scaled_var_adjusted)[j])]
#       }
#     }
#
#     results[[1]] <- scaled_betas
#     results[[2]] <- scaled_cov_adjusted
#     results[[3]] <- scaled_var_adjusted
#
#   } else {
#
#     # No adjustments
#     scaled_cov_new <- list()
#     for (i in 1:length(scaled_cov)) {
#       cov.names <- dimnames(scaled_cov[[i]])  # save dimnames
#
#       # Scale covariance matrices
#       scaled_cov_new[[i]] <- cov.scalar[i]*scaled_cov[[i]]
#     }
#
#     # VARIANCES
#     scaled_var_new <- scaled_var
#     for (i in 1:nrow(scaled_var)) {
#
#       # Scale variances
#       scaled_var_new[i,] <- cov.scalar[i]*scaled_var[i,]
#     }
#
#     results[[1]] <- scaled_betas
#     results[[2]] <- scaled_cov_new
#     results[[3]] <- scaled_var_new
#   }
#   return (results)
# }

### IGNORE

# ## VARIANCES ONLY ##
# ## scale.betas for variances
# scale.betas <- function (thetas) {
#   results <- list()
#   scaled_betas <- combined_betas_only
#   scaled_var <- VAR
#   scaled_cov <- cov_matrix
#   gm_thetas <- prod(thetas)^(1/length(thetas))
#   scalars.normalized <- thetas/gm_thetas
#   for (i in 1:length(thetas)){
#     scaled_betas[i, ] <- scaled_betas[i, ]*scalars.normalized[i]
#   }
#   # normalize it and then square
#   # need to fix in order to accommodate for per-site cov matrices
#   # this can be easily implemented by taking each normalized scalar (theta/gm_thetas)
#   # and then squaring them to obtain the cov matrix scalar for each site
#   cov.scalar <- prod(scalars.normalized^2)
#
#   scaled_var_adjusted <- scaled_var
#   scaled_var_new <- scaled_var
#   for (i in 1:nrow(scaled_var)) {
#
#     # Scale variances
#     scaled_var_new[i,] <- cov.scalar[i]*scaled_var[i,]
#
#     # Adjustments to variance
#     scaled_var_adjusted[i,] <- pmax(scaled_var_new[i,], (.Machine$double.eps)^(1/3))
#     scaled_var_adjusted[i,] <- pmin(scaled_var_new[i,], (.Machine$double.eps)^(-1/3))
#
#   }
#   #####
#   #######
#   results[[1]] <- scaled_betas
#   results[[2]] <- scaled_cov_adjusted
#   results[[3]] <- scaled_var_adjusted
#   return (results)
# }

### END IGNORE


# # Function to calculate log-likelihood of meta-analytic regression with given scaling parameters for beta
# # Takes in nr_thetas, meaning non-reference thetas: a vector of all thetas in order of study site, excluding the reference site
# meta.ll <- function (log.thetas, adjust = FALSE) {
#   nr_thetas <- exp(log.thetas)  # exponentiate after optimizing log.thetas
#   thetas <- c(nr_thetas[0:(scale_reference-1)], 1, nr_thetas[scale_reference:length(nr_thetas)])
#   scaled_betas <- scale.betas(thetas, adjust = adjust)
#   #-----------------
#   # Formatting betas to a format the meta_regression can take
#   # Not sure if there is a more efficient way
#   betas_vec <- as.vector(t(scaled_betas[[1]]))
#   #-----------------
#   # Retrieve variances
#   covs <- scaled_betas[[2]]  # extract scaled covariance matrices
#   scaled_VAR <- scaled_betas[[3]]  # simulate structure of combined betas dataframe
#   # for (i in 1:length(thetas)) {
#   #   for (j in 1:ncol(scaled_VAR)) {
#   #     # Add variances to corresponding predictor and study
#   #     scaled_VAR[i,j] <- diag(covs[[i]])[which(names(diag(covs[[i]])) == colnames(scaled_VAR)[j])]
#   #   }
#   # }
#
#   model <- meta_regression(beta = betas_vec, variance = as.vector(t(scaled_VAR)))
#   ll <- model[["fit.stats"]]["ll", "ML"]
#   return(ll)
# }
#
#
# # Estimate scaling parameters theta by maximizing meta analytic regression likelihood (log likelihood is maximized here)
#
# # try L-BFGS-B
#
# # optim.results <- optim(c(1,1,1), fn=meta.ll, method="L-BFGS-B", control = list(fnscale = -1))
# optim.results <- optim(c(0,0,0), fn=meta.ll, control=list(fnscale=-1))
#
#
# # saveRDS(optim.results, file = "../output/scaling_thetas2.rds")
# # saveRDS(optim.results, file=file.path("output", "scaling_thetas.rds"))
# # saveRDS(optim.results, file=file.path("output", "scaling_thetas_LBFGSB.rds"))
# #
# # thetas.maxll <- readRDS("output/scaling_thetas.rds")
# # thetas.maxll <- readRDS("output/scaling_thetas_LBFGSB.rds")
# thetas.maxll <- readRDS("../output/scaling_thetas2.rds")
#
# #-----------------
# # Fitting model with new thetas
# exp.thetas <- exp(thetas.maxll[["par"]])
# thetas <- c(exp.thetas[0:(scale_reference-1)], 1, exp.thetas[scale_reference:length(exp.thetas)])
# # thetas.normalized <- thetas/(prod(thetas)^(1/length(thetas)))  # for checking thetas
# scaled_betas <- scale.betas(thetas)
# betas_vec <- as.vector(t(scaled_betas[[1]]))
# covs <- scaled_betas[[2]]
# scaled_VAR <- scaled_betas[[3]]
# # Scaled meta-analytic regression model
# scaled.meta <- meta_regression(beta = betas_vec, variance = as.vector(t(scaled_VAR)))
#
# # saveRDS(scaled.meta, file = "../output/scaled_meta2.rds")
# # saveRDS(scaled.meta, file=file.path("output", "scaled_meta.rds"))
# # scaled.meta <- readRDS("output/scaled_meta.rds")
# scaled.meta <- readRDS("../output/scaled_meta.rds")
#
# #-----------------
# unscaled.meta[["fit.stats"]]["ll", "ML"]
# scaled.meta[["fit.stats"]]["ll", "ML"]
#
# unscaled.meta[["fit.stats"]]["AICc", "ML"]
# scaled.meta[["fit.stats"]]["AICc", "ML"]
#
# unscaled.meta[["fit.stats"]]
# scaled.meta[["fit.stats"]]
#
# # Check scaling parameters
# thetas  # 1.547785e+02 3.734221e-06 1.000000e+00 2.502109e+02
# thetas.normalized  # 2.509900e+02 6.055441e-06 1.621607e+00 4.057439e+02




#############################
# OLD COV ADJUSTMENT CODE
#############################
# Model output for each study can be found here:
# chocBetas <- readRDS("data/choc_non_parallel_new_average_betas.rds")
# pensBetas <- readRDS("data/pens_non_parallel_new_average_betas.rds")
# IRLBetas <- readRDS("data/IRL_non_parallel_new_average_betas.rds")
# tampaBetas <- readRDS("data/tampa_non_parallel_new_average_betas.rds")

###############
# Define predictor columns
# predictor_columns <- colnames(combined_betas_only)

# Create the formula string dynamically
# formula_string <- paste("~", paste(predictor_columns, collapse = " + "))
# formula <- as.formula(formula_string)

# Prepare the beta estimates
# betas <- combined_betas_only


###########
#
# eigen_decomp <- eigen(cov_matrix)
# eigenvalues <- eigen_decomp$values
# eigenvectors <- eigen_decomp$vectors
#
#
# # Find the smallest positive eigenvalue
# smallest_eigenvalue <- min(eigenvalues[eigenvalues > 0])
#
# # # Define the maximum allowed variance
# # max_allowed_variance <- sqrt(1 / .Machine$double.eps) * smallest_eigenvalue
#
# # Adjust eigenvalues
# adjusted_eigenvalues <- pmax(eigenvalues, (.Machine$double.eps)^(1/3))
# adjusted_eigenvalues <- pmin(eigenvalues, (.Machine$double.eps)^(-1/3))
# # adjusted_eigenvalues <- pmin(eigenvalues, max_allowed_variance)
#
# # adjusted_eigenvalues <- pmax(adjusted_eigenvalues, (.Machine$double.eps)^(1/3))
# # adjusted_eigenvalues <- pmin(adjusted_eigenvalues, (.Machine$double.eps)^(-1/3))
# adjusted_cov_matrix <- eigenvectors %*% diag(adjusted_eigenvalues) %*% t(eigenvectors)
#
# # another check
# smallest_eigenvalue <- min(adjusted_eigenvalues[adjusted_eigenvalues > 0])
# # smallest_eigenvalue
#
# #####
# # Recommendation: Add a small jitter to the diagonal to ensure positive definiteness
# epsilon <- 1e-6
# adjusted_cov_matrix <- adjusted_cov_matrix + diag(epsilon, nrow(adjusted_cov_matrix))


## FULL MATRIX ADJUSTMENT ##

# adjusted_cov_matrix <- list()
# # Adjustments to covariance matrices
# for (i in 1:length(cov_matrix)) {
#
#   cov_names <- dimnames(cov_matrix[[i]])  # save dimnames
#
#   eigen_decomp <- eigen(cov_matrix[[i]])
#   eigenvalues <- eigen_decomp$values
#   eigenvectors <- eigen_decomp$vectors
#
#
#   # Find the smallest positive eigenvalue
#   smallest_eigenvalue <- min(eigenvalues[eigenvalues > 0])
#
#   # # Define the maximum allowed variance
#   # max_allowed_variance <- sqrt(1 / .Machine$double.eps) * smallest_eigenvalue
#
#   # Adjust eigenvalues
#   adjusted_eigenvalues <- pmax(eigenvalues, (.Machine$double.eps)^(1/3))
#   adjusted_eigenvalues <- pmin(adjusted_eigenvalues, (.Machine$double.eps)^(-1/3))
#   # adjusted_eigenvalues <- pmin(eigenvalues, max_allowed_variance)
#
#   # adjusted_eigenvalues <- pmax(adjusted_eigenvalues, (.Machine$double.eps)^(1/3))
#   # adjusted_eigenvalues <- pmin(adjusted_eigenvalues, (.Machine$double.eps)^(-1/3))
#   adjusted_cov_matrix[[i]] <- eigenvectors %*% diag(adjusted_eigenvalues) %*% t(eigenvectors)
#
#   # another check
#   smallest_eigenvalue <- min(adjusted_eigenvalues[adjusted_eigenvalues > 0])
#   # smallest_eigenvalue
#
#   # Set very large variances to very large value
#   large_value <- (.Machine$double.eps)^(-1/3)
#   diag(adjusted_cov_matrix[[i]])[which(diag(adjusted_cov_matrix[[i]])>large_value)] <- large_value
#
#   #####
#   # Recommendation: Add a small jitter to the diagonal to ensure positive definiteness
#   epsilon <- 1e-6
#   adjusted_cov_matrix[[i]] <- adjusted_cov_matrix[[i]] + diag(epsilon, nrow(adjusted_cov_matrix[[i]]))
# ####
#
#   dimnames(adjusted_cov_matrix[[i]]) <- cov_names
# }
#
# # Grab variances from each matrix
# VAR <- combined_betas_only  # simulate structure of combined betas dataframe
# for (i in 1:nrow(VAR)) {
#   for (j in 1:ncol(VAR)) {
#
#     # Add variances to corresponding predictor and study
#     VAR[i,j] <- diag(adjusted_cov_matrix[[i]])[which(names(diag(adjusted_cov_matrix[[i]])) == colnames(VAR)[j])]
#   }
# }


## VARIANCES ONLY ##

# # Grab variances from each unadjusted matrix
# VAR <- combined_betas_only  # simulate structure of combined betas dataframe
# adjusted_VAR <- VAR
# for (i in 1:nrow(VAR)) {
#   for (j in 1:ncol(VAR)) {
#
#     # Add variances to corresponding predictor and study
#     VAR[i,j] <- diag(cov_matrix[[i]])[which(names(diag(cov_matrix[[i]])) == colnames(VAR)[j])]
#   }
#
#   adjusted_VAR[i,] <- pmax(VAR[i,], (.Machine$double.eps)^(1/3))
#   adjusted_VAR[i,] <- pmin(adjusted_VAR[i,], (.Machine$double.eps)^(-1/3))
# }

# adjusted_VAR <- pmax(VAR, (.Machine$double.eps)^(1/3))
# adjusted_VAR <- pmin(VAR, (.Machine$double.eps)^(-1/3))




# If keeping some columns constant when scaling betas
#
# scale_betas <- function (
#     thetas,
#     keep = NULL
# ) {
#   results <- list()
#   scaled_betas <- combined_betas_only
#   scaled_cov <- cov_matrix
#   gm_thetas <- prod(thetas)^(1/length(thetas))
#   if (keep == NULL) {
#     for (i in 1:length(thetas)) {
#       scaled_betas[i, ] <- scaled_betas[i, ]*thetas[i]/gm_thetas
#     }
#   }
#   else {
#     for (i in 1:length(thetas)) {
#       scaled_betas[i, -keep] <- scaled_betas[i, ]*thetas[i]/gm(thetas)
#     }
#   }
#   scaled_cov <- (thetas%*%thetas)/(gm_thetas^2)*scaled_cov
#   results[[1]] <- scaled_betas
#   results[[2]] <- scaled_cov
#   return (results)
# }



#########################################
# SANDBOX
#########################################
# run metafor meta-analytic regression
# Y_i ~ mu + RE_i + E_i
# WORKING, INTERPRETABLE
# result <- rma.mv(
#   yi = beta, # Vector of all beta coefficients (may need effect size e.g. odds ratios from BUPD.R output instead of scaled betas here)
#   V = variance, # CHECK AND INPROVE THIS!!! vector of length k with the corresponding sampling variances or a k x k variance-covariance matrix of the sampling errors.
#   method = "REML", # default
#   random = ~ 1 | study, # Random effects for studies
#   struct =  "UN", # unstructured varcov matrix
#   mods = ~ predictor, # Including predictors as fixed effects without an intercept (AKA is predictors in formula)
#   verbose = TRUE
# )
# summary(result)

#
# scale.betas <- function (thetas) {
#   results <- list()
#   scaled_betas <- combined_betas_only
#   gm_thetas <- prod(thetas)^(1/length(thetas))
#   for (i in 1:length(thetas)){
#     scaled_betas[i, ] <- scaled_betas[i, ]*thetas[i]/gm_thetas
#   }
#   # Generate covariance matrix
#   scaled_cov <- cov(scaled_betas, use = "pairwise.complete.obs") # calculates the correlation between each pair of variables using all complete pairs of observations for those variables
#
#   # Find where there's missing values
#   missing_values <- is.na(scaled_cov)
#
#   # Set missing off-diagonals to zero
#   scaled_cov[missing_values & !row(scaled_cov) == col(scaled_cov)] <- 0
#
#   # Set missing variances to very large value
#   large_value <- 10000
#   diag(scaled_cov)[missing_values[diag(TRUE, nrow(scaled_cov))]] <- large_value
#   # View(cov_matrix)
#   ####### cov matrix adjustments
#   eigen_decomp <- eigen(scaled_cov)
#   eigenvalues <- eigen_decomp$values
#   eigenvectors <- eigen_decomp$vectors
#
#   # Find the smallest positive eigenvalue
#   smallest_eigenvalue <- min(eigenvalues[eigenvalues > 0])
#
#   # # Define the maximum allowed variance
#   # max_allowed_variance <- sqrt(1 / .Machine$double.eps) * smallest_eigenvalue
#
#   # Adjust eigenvalues
#   adjusted_eigenvalues <- pmax(eigenvalues, (.Machine$double.eps)^(1/3))
#   adjusted_eigenvalues <- pmin(eigenvalues, (.Machine$double.eps)^(-1/3))
#   # adjusted_eigenvalues <- pmin(eigenvalues, max_allowed_variance)
#
#   # adjusted_eigenvalues <- pmax(adjusted_eigenvalues, (.Machine$double.eps)^(1/3))
#   # adjusted_eigenvalues <- pmin(adjusted_eigenvalues, (.Machine$double.eps)^(-1/3))
#   scaled_cov_adjusted <- eigenvectors %*% diag(adjusted_eigenvalues) %*% t(eigenvectors)
#
#   # another check
#   smallest_eigenvalue <- min(adjusted_eigenvalues[adjusted_eigenvalues > 0])
#   # smallest_eigenvalue
#
#   #####
#   # Recommendation: Add a small jitter to the diagonal to ensure positive definiteness
#   epsilon <- 1e-6
#   scaled_cov_adjusted <- scaled_cov_adjusted + diag(epsilon, nrow(scaled_cov_adjusted))
#   #####
#   #######
#   results[[1]] <- scaled_betas
#   results[[2]] <- scaled_cov_adjusted
#   return (results)
# }
