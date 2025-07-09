# PCA for predictors to test for collinearity

library(FactoMineR)
library(dplyr)

# Retrieve predictors
pred <- readRDS("data/predictors_kriged_standardized.rds")
pred <- as.data.frame(pred)  # convert to dataframe so "geometry" is not selected
pred_test <- pred %>%
  select(all_of(predictors)) %>%
  mutate(across(all_of(predictors), as.character)) %>%
  mutate(across(all_of(predictors), as.numeric)) %>%
  mutate(study = pred$study)
# Separate study sites
choc <- pred_test[pred_test$study == "choc", predictors]
pens <- pred_test[pred_test$study == "pens", predictors]
tampa <- pred_test[pred_test$study == "tampa", predictors]
IRL <- pred_test[pred_test$study == "IRL", predictors]
pred_test <- pred_test[, predictors]  # remove study column

sites <- list(choc, pens, tampa, IRL)  # combine studies

# Remove columns with no information
for (m in 1:length(sites)) {
  missed <- c()
  for (i in 1:length(predictors)) {
    if (all(is.na(sites[[m]][,predictors[i]]))) {
      missed[i] <- predictors[i]  # track predictors with no info (all NA)
      missed <- na.omit(missed)
    }
  }
  sites[[m]] <- sites[[m]] %>%
    select(-all_of(missed))
}

# Combined PCA
pca.pred <- FactoMineR::PCA(pred_test)

# Study-specific PCA
pca.choc <- FactoMineR::PCA(sites[[1]])
pca.pens <- FactoMineR::PCA(sites[[2]])
pca.tampa <- FactoMineR::PCA(sites[[3]])
pca.IRL <- FactoMineR::PCA(sites[[4]])

# Check correlation between predictors for each site
cor.choc <- na.omit(pca.choc$var$cor)
cor.pens <- na.omit(pca.pens$var$cor)
cor.tampa <- na.omit(pca.tampa$var$cor)
cor.IRL <- na.omit(pca.IRL$var$cor)

