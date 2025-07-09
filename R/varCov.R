# Function for retrieving/cleaning variance-covariance matrices from model selection
#' @export
varCov <- function(data, mods, Betas) {

  # Check if `mods` exists
  if(!exists("mods")) {
    stop("No models found. Please assign the selected model from each study to a list named `mods`.")
  } else { message(paste("Identified selected models from", length(mods), "studies.")) }

  # Check if `Betas` exists
  if(!exists("Betas")) {
    stop("No betas found. Please assign beta coefficient estimates from selected models to a list named `Betas`.")
  } else { message(paste("Identified beta coefficient estimates from", length(Betas), "studies.")) }

  # Retrieve cleaned predictor data and study names
  pred <- data$predictors
  studies <- unique(data$state$study)

  # Get predictors (excluding study and definitions)
  pred <- pred %>%
    dplyr::mutate(dplyr::across(c("OBJECTID", "ID", "bmpCountv5", "n", "distance", "X", "Y"), as.character))
  ## keep non predictor data as characters to avoid being selected as predictors
  numeric_pred <- pred %>%
    dplyr::select_if(is.numeric)
  factor_pred <- pred %>%
    dplyr::select_if(is.factor)
  numeric_pred_cols <- colnames(cbind(factor_pred, numeric_pred))

  # Retrieve covariance matrices from each model
  cov_matrix <- list(mods[[1]]$COV,
                     mods[[2]]$COV,
                     mods[[3]]$COV,
                     mods[[4]]$COV)


  # Source getBetas.R script for predictors and to combine betas
  # source("R/getBetas.R")

  # Retrieve and combine betas (using getBetas.R)
  GB <- getBetas(data, Betas)

  # Re-assign outputs for easy access
  combined_betas <- GB$combined_betas
  combined_betas_only <- GB$combined_betas_only
  combined_se <- GB$combined_se
  combined_se_only <- GB$combined_se_only

  for (i in 1:length(cov_matrix)) {

    # # Symmetrization of matrices
    # cov_matrix[[i]] <- (cov_matrix[[i]] + t(cov_matrix[[i]]))/2

    # Remove intercepts
    cov_matrix[[i]] <- cov_matrix[[i]][2:(dim(cov_matrix[[i]])[[1]]),
                                       2:(dim(cov_matrix[[i]])[[2]])]
  }

  # Find where there are missing values
  missing_values <- list()
  for (i in 1:length(cov_matrix)) {

    # Find non-selected predictors
    selected <- dimnames(cov_matrix[[i]])
    missing_preds <- setdiff(numeric_pred_cols, selected[[1]])
    missing_col <- matrix(NA, length(selected[[1]]), length(missing_preds))
    missing_row <- matrix(NA, length(missing_preds), length(c(selected[[1]], missing_preds)))
    cov_matrix[[i]] <- cbind(cov_matrix[[i]], missing_col)
    cov_matrix[[i]] <- rbind(cov_matrix[[i]], missing_row)
    dimnames(cov_matrix[[i]]) <- list(c(selected[[1]], missing_preds),
                                      c(selected[[1]], missing_preds))

    # Predictors with missing values
    missing_values[[i]] <- is.na(cov_matrix[[i]])

    # Set missing off-diagonals to zero
    cov_matrix[[i]][missing_values[[i]] & !row(cov_matrix[[i]]) == col(cov_matrix[[i]])] <- 0

    # Set missing variances to very large value
    # large_value <- 10000
    large_value <- (.Machine$double.eps)^(-1/3)
    diag(cov_matrix[[i]])[missing_values[[i]][diag(TRUE, nrow(cov_matrix[[i]]))]] <- large_value
    # # Set infinite variances to very large value
    # diag(cov_matrix[[i]])[which(is.infinite(diag(cov_matrix[[i]])))] <- large_value
    # # Set very large variances to very large value
    # diag(cov_matrix[[i]])[which(diag(cov_matrix[[i]])>large_value)] <- large_value

  }

  # Source matrix cleaning script
  # source("inst/scripts/clean_matrix.R")

  # Clean each covariance matrix
  adjusted_cov_matrix <- list()
  for (i in 1:length(cov_matrix)) {
    adjusted_cov_matrix[[i]] <- clean_matrix(cov_matrix[[i]])
  }


  # Grab variances from each matrix
  VAR <- combined_betas_only  # simulate structure of combined betas dataframe
  for (i in 1:nrow(VAR)) {
    for (j in 1:ncol(VAR)) {

      # Add variances to corresponding predictor and study
      # VAR[i,j] <- diag(cov_matrix[[i]])[which(names(diag(cov_matrix[[i]])) == colnames(VAR)[j])]
      VAR[i,j] <- diag(adjusted_cov_matrix[[i]])[which(names(diag(adjusted_cov_matrix[[i]])) == colnames(VAR)[j])]
    }
  }

  return(list(VAR = VAR,
              adjusted_cov_matrix = adjusted_cov_matrix,
              combined_betas = combined_betas,
              combined_betas_only = combined_betas_only,
              combined_se = combined_se,
              combined_se_only = combined_se_only))

}







## OLD CODE ##

# # Only using variances for now
# cov_matrix <- list(choc = ctmm::pd.solve(choc_mod_new$Hessian),
#                    pens = ctmm::pd.solve(pens_mod_new$Hessian),
#                    tampa = ctmm::pd.solve(tampa_mod_new$Hessian),
#                    IRL = ctmm::pd.solve(IRL_mod_new$Hessian))

# Function to update thresholds from probit regression (EQUIDISTANT)
# zetas <- function(model, betas) {
#
#   ## CHECK
#   avg_zeta <- mean(model$zeta)  # average zetas
#   d_zeta <- diff(model$zeta)  # model$zeta[2]-model$zeta[1]  # difference in zetas
#   int <- -avg_zeta  # new intercept
#   # adjust <- (2*qnorm(0.975))/d_zeta  # threshold adjustment
#   adjust <- (d_zeta/diff(pnorm(model$zeta))) * ((0.975-0.025)/(2*qnorm(0.975)))
#   betas_new <- betas
#   betas_new[,1] <- adjust*betas[,1]  # new betas
#
#   # Adjust covariance matrix
#   matrix_initial <- ctmm::pd.solve(model$Hessian)  # retrieve original COV
#   matrix_new <- (adjust^2)*matrix_initial
#
#   return(list(zetas = c(avg_zeta, d_zeta),
#               intercept = int,
#               betas_new = betas_new,
#               cov_matrix = matrix_new))
#
# }
#
# # Threshold adjustments
# updated_mods <- list()
# cov_matrix <- list()  # store all COVs together
# for (i in 1:length(studies)) {
#   updated_mods[[i]] <- zetas(studies[[i]], studies_betas[[i]])
#   cov_matrix[[i]] <- updated_mods[[i]]$cov_matrix
# }

# Reassign betas for getBetas.R
# chocBetas <- updated_mods[[1]]$betas_new
# pensBetas <- updated_mods[[2]]$betas_new
# tampaBetas <- updated_mods[[3]]$betas_new
# IRLBetas <- updated_mods[[4]]$betas_new



# # Generate covariance matrix
# # cov_matrix <- cov(combined_betas_only, use = "pairwise.complete.obs") # calculates the correlation between each pair of variables using all complete pairs of observations for those variables
#
# # Predictors with missing values
# missing_values <- is.na(cov_matrix)
#
# # Set missing off-diagonals to zero
# cov_matrix[missing_values & !row(cov_matrix) == col(cov_matrix)] <- 0
#
# # Set missing variances to very large value
# # large_value <- 10000
# large_value <- (.Machine$double.eps)^(-1/3)
# diag(cov_matrix)[missing_values[diag(TRUE, nrow(cov_matrix))]] <- large_value
# # Set infinite variances to very large value
# # diag(cov_matrix)[which(is.infinite(diag(cov_matrix[[3]])))] <- large_value

# # Retrieve covariance matrices from model outputs
# cov_matrix <- list(choc = vcov(choc_mod_new), pens = vcov(pens_mod_new),
#                    tampa = vcov(tampa_mod_new), IRL = vcov(IRL_mod_new))



