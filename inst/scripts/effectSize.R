# # Number of columns
# num_columns <- ncol(combined_betas_only)

# # Initialize vectors to store overall effect sizes and their standard errors
# overall_effect <- numeric(num_columns)
# overall_standard_error <- numeric(num_columns)

# # Loop over each column
# for (i in 1:num_columns) {
#   estimates <- combined_betas_only[, i]
#   standard_errors <- combined_se_only[, i]

#   # Calculate variances
#   variances <- standard_errors^2

#   # Calculate weights
#   weights <- 1 / variances

#   # Calculate the weighted mean effect size
#   weighted_mean <- sum(weights * estimates) / sum(weights)

#   # Calculate the variance of the combined effect size
#   combined_variance <- 1 / sum(weights)

#   # Calculate the standard error of the combined effect size
#   combined_standard_error <- sqrt(combined_variance)

#   # Store the results
#   overall_effect[i] <- weighted_mean # effect of all four studies averaged for each predictor.
#   overall_standard_error[i] <- combined_standard_error
# }




# #####
# col_names <- colnames(combined_betas_only)
# overall_effect_df <- as.data.frame(t(overall_effect))  
# colnames(overall_effect_df) <- col_names 
# View(overall_effect_df)

# ##### average betas #####
# average_betas <- data.frame(t(colMeans(combined_betas_only, na.rm = TRUE)))
# view(average_betas)

