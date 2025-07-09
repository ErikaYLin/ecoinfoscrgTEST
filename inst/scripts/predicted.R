## NEW MODEL ## ----

# For planar data
sf::sf_use_s2(FALSE)

# Load requisite packages
library(dplyr)
library(sf)
library(ggplot2)
basemaps::set_defaults(map_service = "esri", map_type = "world_imagery")

# Import meta-analytic regression model
meta_analysis <- readRDS("../output/scaled_meta2.rds")  # scaled model
meta_analysis <- readRDS("../output/unscaled_meta3.rds")  # unscaled model

# Retrieve studies for predictions
pred <- readRDS("data/predictors_kriged_standardized.rds")
colnames(pred)[41] <- "predicted_OLD"  # rename to compare to new predictions
predictor_data <- as.data.frame(pred)
predictors <- meta_analysis$data$predictor[1:248]  # retrieve predictors

# Remove non-predictor columns
predictor_data <- predictor_data %>%
  select(all_of(predictors)) %>%
  select(order(colnames(.)))
# Reformat predictor columns
predictor_data <- predictor_data %>%
  mutate(across(all_of(predictors), as.character)) %>% # convert predictors to character first (needed for factorized columns)
  mutate(across(all_of(predictors), as.numeric)) %>% # convert all predictors to numeric
  mutate(across(all_of(predictors), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))  # sub NAs for statewide means
## NAs in the data can result in predicted values to be NA

## CALCULATE AND STORE STATEWIDE MEANS FOR PREDICTORS AS FILL-INS FOR EMPTY/MISSING PREDICTORS

# Extract coefficients and intercept from final meta-analysis model
betas <- meta_analysis$beta  # IS THIS WHERE THE NEW INTERCEPT GOES???
# Remove "predictor" from predictor names
rownames(betas) <- gsub("predictor", "", rownames(betas))  # match predictor names
setdiff(colnames(predictor_data), rownames(betas))  # which predictor used as the intercept?
rownames(betas)[1] <- setdiff(colnames(predictor_data), rownames(betas))  # replace with predictor name
betas <- as.vector(betas)

# # Create the new model matrix and remove the intercept
# predgrid <- model.matrix(~predictors, data = predictor_data)[,-1]  ## TEST ##
# colnames(predgrid) <- gsub("predictors", "predictor", colnames(predgrid))
#
# test <- metafor::predict.rma(meta_analysis, newmods = predictor_data)

# Multiply predictors by the meta-analysis model betas
expected_outcome <- as.matrix(predictor_data) %*% betas

# Inverse probit for predicted values
predicted <- pnorm(expected_outcome)

# Add predictions to predictors data
pred$predicted <- predicted



# Load in original predictions and data for each site
choc <- st_transform(st_read("../Data/choctawatchee_bay/choctawatchee_bay_lssm_POINTS_0.001deg.shp"), crs = 6346) # Response:
pens <- st_transform(st_read("../Data/santa_rosa_bay/Santa_Rosa_Bay_Living_Shoreline_POINTS_0.001deg.shp"), crs = 6346) # Response:
tampa <- st_transform(st_read("../Data/tampa_bay/Tampa_Bay_Living_Shoreline_Suitability_Model_Results_POINTS_0.001deg.shp"), crs = 6346) # Response:
IRL <- st_transform(st_read("../Data/indian_river_lagoon/UCF_livingshorelinemodels_MosquitoNorthIRL_111m.shp"), crs = 6346) # Response:

# Full data to fix transformed geometry
choc_full <- st_transform(st_read("data/choctawatchee_bay/choctawatchee_bay_lssm_POINTS_0.001deg.shp"))
IRL_full <- st_transform(st_read("../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted.shp"))
pens_full <- st_transform(st_read("../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted.shp"))
tampa_full <- st_transform(st_read("../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted.shp"))

coo <- c("OBJECTID", "ID", "geometry", "feature_x", "feature_y", "nearest_x", "nearest_y", "shape__len", "Shape__Len", "distance", "distance_2", "n", "x", "y", "X", "Y")
pred_dat <- as.data.frame(pred)

# Add predictions to each site
## Choctawatchee
choc_available_columns <- coo[coo %in% colnames(choc)]
choc_coo <- choc %>% select(all_of(choc_available_columns))
choc_pred <- pred_dat[pred_dat$study == "choc",]
choc_pred$predicted <- pred_dat[pred_dat$study == "choc", "predicted"]
choc_cols <- intersect(colnames(choc_pred), colnames(choc_coo))
choc_pred <- choc_pred %>%
  select(!all_of(choc_cols))

choc_predicted <- cbind(choc_coo, choc_pred)
# st_write(choc_predicted, "../output/Final_Shapefile_all_data/Choctawatchee Bay/choc_predicted_new.shp", append = FALSE)
choc_test <- choc_predicted
st_geometry(choc_test) <- choc_full$geometry
choc_test <- choc_test %>%
  select(all_of(c("feature_x", "feature_y", "predicted")))  # automatically selects geometry
colnames(choc_test)[c(1,2)] <- c("X", "Y")
st_write(choc_test, "../output/Final_Shapefile_all_data/Choctwatchee Bay/choc_predicted_shape.shp")

plot_choc <- st_transform(choc_test, crs = 3857)   # st_transform(st_read("../output/Final_Shapefile_all_data/Choctawatchee Bay/choc_predicted_shape.shp"), crs = 3857)
choc_bb <- st_bbox(c(xmin = -9644651-15000,
                     ymin = 3551143-15000,
                     xmax = -9584357+15000,
                     ymax = 3570713+15000), crs = 3857)
choc_coords <- as.data.frame(st_coordinates(plot_choc))
plot_choc2 <- cbind(plot_choc[,-c(1,2)], choc_coords)

choc_map <- basemaps::basemap_ggplot(ext = choc_bb) +  # ESRI World Imagery basemap
  geom_sf(data = plot_choc2, aes(x = X, y = Y, color = predicted), size = 1) +
  labs(x = "Latitude", y = "Longitude", title = "Choctawatchee Bay") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.08,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  # scale_color_viridis_c(name = "Predicted", option = "C")
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = choc_map, file = "../QGIS/choc_predicted_map.png", width = 6, height = 6)


## Pensacola Bay
pens_available_columns <- coo[coo %in% colnames(pens)]
pens_coo <- pens %>% select(all_of(pens_available_columns))
pens_pred <- pred_dat[pred_dat$study == "pens",]
pens_pred$predicted <- pred_dat[pred_dat$study == "pens", "predicted"]
pens_cols <- intersect(colnames(pens_pred), colnames(pens_coo))
pens_pred <- pens_pred %>%
  select(!all_of(pens_cols))

pens_predicted <- cbind(pens_coo, pens_pred)
write.csv(pens_predicted, "../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted_new.csv", row.names = FALSE)
pens_test <- pens_predicted
st_geometry(pens_test) <- pens_full$geometry
pens_test <- pens_test %>%
  select(all_of(c("feature_x", "feature_y", "predicted")))  # automatically selects geometry
colnames(pens_test)[c(1,2)] <- c("X", "Y")
write.csv(pens_test, "../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted_test.csv", row.names = FALSE)
st_write(pens_test, "../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted_shape.shp")


# plot_pens <- st_read("../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted_shape.shp")
plot_pens <- st_transform(st_read("../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted_shape.shp"), crs = 3857)
# plot_pens <- st_crop(plot_pens, c(xmin = -9751391, ymin = 3536051, xmax = -9640820, ymax = 3594841))
pens_bb <- st_bbox(c(xmin = -9732366-13000,
                     ymin = 3542349-13000,
                     xmax = -9659068+13000,
                     ymax = 3589286+13000), crs = 3857)
# plot_pens2 <-st_crop(plot_pens, pens_bb)
# plot_pens2 <- st_intersection(plot_pens, pens_bb)
pens_coords <- as.data.frame(st_coordinates(plot_pens))
plot_pens2 <- cbind(plot_pens[,-c(1,2)], pens_coords)

# ggplot(plot_pens2, aes(x = X, y = Y)) +
#   geom_sf(aes(color = predicted)) +
#   basemaps::basemap_gglayer(ext = pens_bb)
# ggplot(plot_pens2, aes(x = X, y = Y)) +
#   geom_sf(aes(color = predicted))

pens_map <- basemaps::basemap_ggplot(ext = pens_bb) +
  geom_sf(data = plot_pens2, aes(x = X, y = Y, color = predicted), size = 1) +
  labs(x = "Latitude", y = "Longitude", title = "Pensacola Bay") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.5,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  # scale_color_viridis_c(name = "Predicted", option = "C")
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = pens_map, file = "../QGIS/pens_predicted_map.png", width = 6, height = 6)


# Tampa Bay
tampa_available_columns <- coo[coo %in% colnames(tampa)]
tampa_coo <- tampa %>% select(all_of(tampa_available_columns))
tampa_pred <- pred_dat[pred_dat$study == "tampa",]
tampa_pred$predicted <- pred_dat[pred_dat$study == "tampa", "predicted"]
tampa_cols <- intersect(colnames(tampa_pred), colnames(tampa_coo))
tampa_pred <- tampa_pred %>%
  select(!all_of(tampa_cols))

tampa_predicted <- cbind(tampa_coo, tampa_pred)
# write.csv(tampa_predicted, "../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted_new.csv", row.names = FALSE)
tampa_test <- tampa_predicted
st_geometry(tampa_test) <- tampa_full$geometry
# write.csv(tampa_test, "../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted_test.csv", row.names = FALSE)
tampa_test <- tampa_test %>%
  select(all_of(c("feature_x", "feature_y", "predicted")))  # automatically selects geometry
colnames(tampa_test)[c(1,2)] <- c("X", "Y")
st_write(tampa_test, "../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted_shape.shp")

plot_tampa <- st_transform(tampa_test, crs = 3857)   # st_transform(st_read("../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted_shape.shp"), crs = 3857)
tampa_bb <- st_bbox(c(xmin = -9222837-15000,
                      ymin = 3178783-15000,
                      xmax = -9161061+15000,
                      ymax = 3253985+15000), crs = 3857)
tampa_coords <- as.data.frame(st_coordinates(plot_tampa))
plot_tampa2 <- cbind(plot_tampa[,-c(1,2)], tampa_coords)

tampa_map <- basemaps::basemap_ggplot(ext = tampa_bb) +
  geom_sf(data = plot_tampa2, aes(x = X, y = Y, color = predicted), size = 1) +
  labs(x = "Latitude", y = "Longitude", title = "Tampa Bay") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.00013,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  # scale_color_viridis_c(name = "Predicted", option = "C")
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = tampa_map, file = "../QGIS/tampa_predicted_map.png", width = 6, height = 6)


# IRL
IRL_available_columns <- coo[coo %in% colnames(IRL)]
IRL_coo <- IRL %>% select(all_of(IRL_available_columns))
IRL_pred <- pred_dat[pred_dat$study == "IRL",]
IRL_pred$predicted <- pred_dat[pred_dat$study == "IRL", "predicted"]
IRL_cols <- intersect(colnames(IRL_pred), colnames(IRL_coo))
IRL_pred <- IRL_pred %>%
  select(!all_of(IRL_cols))

IRL_predicted <- cbind(IRL_coo, IRL_pred)
write.csv(IRL_predicted, "../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted_new.csv", row.names = FALSE)
IRL_test <- IRL_predicted
st_geometry(IRL_test) <- IRL_full$geometry
# write.csv(IRL_test, "../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted_test.csv", row.names = FALSE)
IRL_test <- IRL_test %>%
  select(all_of(c("X", "Y", "predicted")))  # automatically selects geometry
st_write(IRL_test, "../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted_shape.shp")

plot_IRL <- st_transform(IRL_test, crs = 3857)   # st_transform(st_read("../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted_shape.shp"), crs = 3857)
IRL_bb <- st_bbox(c(xmin = -9008088-15000,
                    ymin = 3315585-15000,
                    xmax = -8976838+15000,
                    ymax = 3385010+15000), crs = 3857)
IRL_coords <- as.data.frame(st_coordinates(plot_IRL))
plot_IRL2 <- cbind(plot_IRL[,-c(1,2)], IRL_coords)

IRL_map <- basemaps::basemap_ggplot(ext = IRL_bb) +
  geom_sf(data = plot_IRL2, aes(x = X, y = Y, color = predicted), size = 1) +
  labs(x = "Latitude", y = "Longitude", title = "Indian River Lagoon") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.08,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  # scale_color_viridis_c(name = "Predicted", option = "C")
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = IRL_map, file = "../QGIS/IRL_predicted_map.png", width = 6, height = 6)





## LOCAL BETAS ## ----

# Full data to fix corrupted geometry
choc_full <- st_transform(st_read("data/choctawatchee_bay/choctawatchee_bay_lssm_POINTS_0.001deg.shp"))
IRL_full <- st_transform(st_read("../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted.shp"))
pens_full <- st_transform(st_read("../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted.shp"))
tampa_full <- st_transform(st_read("../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted.shp"))

# Import meta-analytic regression model
# meta_analysis <- readRDS("../output/scaled_meta.rds")  # scaled model
meta_analysis <- readRDS("../output/unscaled_meta.rds")  # unscaled model

# Retrieve studies for predictions
pred <- readRDS("data/predictors_kriged_standardized.rds")
colnames(pred)[41] <- "predicted_OLD"  # rename to compare to new predictions
predictor_data <- as.data.frame(pred)
predictors <- meta_analysis$data$predictor[1:248]  # retrieve predictors

# Remove non-predictor columns
predictor_data <- predictor_data %>%
  select(all_of(predictors)) %>%
  select(order(colnames(.)))
# Reformat predictor columns
predictor_data <- predictor_data %>%
  mutate(across(all_of(predictors), as.character)) %>% # convert predictors to character first (needed for factorized columns)
  mutate(across(all_of(predictors), as.numeric)) %>% # convert all predictors to numeric
  mutate(across(all_of(predictors), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))  # sub NAs for statewide means
## NAs in the data can result in predicted values to be NA

# retrieve vector of studies
predictor_data$study <- pred$study

# coo <- c("OBJECTID", "ID", "geometry", "feature_x", "feature_y", "nearest_x", "nearest_y", "shape__len", "Shape__Len", "distance", "distance_2", "n", "x", "y", "X", "Y")


# Choctawatchee
choc_dat <- predictor_data[predictor_data$study == "choc",]  # subset to choc
choc_dat <- choc_dat[,-ncol(choc_dat)]

betas_choc <- t(combined_betas_only[1,])  # retrieve local betas

# Multiply predictors by the local betas
expected_choc <- as.matrix(choc_dat) %*% as.vector(betas_choc)

# Inverse probit for predicted values
predicted_choc <- pnorm(expected_choc)

# Add predictions to predictors data
choc_dat$predicted <- predicted_choc

# Save predictions as spatial data
choc_test <- st_as_sf(data.frame(X = choc_full$feature_x,
                                 Y = choc_full$feature_y,
                                 predicted = predicted_choc,
                                 geometry = choc_full$geometry))
# Plot predictions
plot_choc <- st_transform(choc_test, crs = 3857)
choc_bb <- st_bbox(c(xmin = -9644651-15000,
                     ymin = 3551143-15000,
                     xmax = -9584357+15000,
                     ymax = 3570713+15000), crs = 3857)
choc_coords <- as.data.frame(st_coordinates(plot_choc))
plot_choc2 <- cbind(plot_choc[,-c(1,2)], choc_coords)

choc_map <- basemaps::basemap_ggplot(ext = choc_bb) +
  geom_sf(data = plot_choc2, aes(color = predicted), size = 0.5) +
  labs(x = "Latitude", y = "Longitude", title = "Choctawatchee Bay") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.5,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = choc_map, file = "../QGIS/choc_predicted_map_localB.png", width = 6, height = 6)


# Pensacola
pens_dat <- predictor_data[predictor_data$study == "pens",]  # subset to pens
pens_dat <- pens_dat[,-ncol(pens_dat)]

betas_pens <- t(combined_betas_only[2,])  # retrieve local betas

# Multiply predictors by the local betas
expected_pens <- as.matrix(pens_dat) %*% as.vector(betas_pens)

# Inverse probit for predicted values
predicted_pens <- pnorm(expected_pens)

# Add predictions to predictors data
pens_dat$predicted <- predicted_pens

# Save predictions as spatial data
pens_test <- st_as_sf(data.frame(X = pens_full$feature_x,
                                 Y = pens_full$feature_y,
                                 predicted = predicted_pens,
                                 geometry = pens_full$geometry))
# Plot predictions
plot_pens <- st_transform(pens_test, crs = 3857)
pens_bb <- st_bbox(c(xmin = -9732366-13000,
                     ymin = 3542349-13000,
                     xmax = -9659068+13000,
                     ymax = 3589286+13000), crs = 3857)
pens_coords <- as.data.frame(st_coordinates(plot_pens))
plot_pens2 <- cbind(plot_pens[,-c(1,2)], pens_coords)

pens_map <- basemaps::basemap_ggplot(ext = pens_bb) +
  geom_sf(data = plot_pens2, aes(color = predicted), size = 0.5) +
  labs(x = "Latitude", y = "Longitude", title = "Pensacola Bay") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.5,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = pens_map, file = "../QGIS/pens_predicted_map_localB.png", width = 6, height = 6)


# Tampa
tampa_dat <- predictor_data[predictor_data$study == "tampa",]  # subset to tampa
tampa_dat <- tampa_dat[,-ncol(tampa_dat)]

betas_tampa <- t(combined_betas_only[3,])  # retrieve local betas

# Multiply predictors by the local betas
expected_tampa <- as.matrix(tampa_dat) %*% as.vector(betas_tampa)

# Inverse probit for predicted values
predicted_tampa <- pnorm(expected_tampa)

# Add predictions to predictors data
tampa_dat$predicted <- predicted_tampa

# Save predictions as spatial data
tampa_test <- st_as_sf(data.frame(X = tampa_full$feature_x,
                                 Y = tampa_full$feature_y,
                                 predicted = predicted_tampa,
                                 geometry = tampa_full$geometry))
# Plot predictions
plot_tampa <- st_transform(tampa_test, crs = 3857)
tampa_bb <- st_bbox(c(xmin = -9222837-15000,
                      ymin = 3178783-15000,
                      xmax = -9161061+15000,
                      ymax = 3253985+15000), crs = 3857)
tampa_coords <- as.data.frame(st_coordinates(plot_tampa))
plot_tampa2 <- cbind(plot_tampa[,-c(1,2)], tampa_coords)

tampa_map <- basemaps::basemap_ggplot(ext = tampa_bb) +
  geom_sf(data = plot_tampa2, aes(color = predicted), size = 0.5) +
  labs(x = "Latitude", y = "Longitude", title = "Tampa Bay") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.5,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = tampa_map, file = "../QGIS/tampa_predicted_map_localB.png", width = 6, height = 6)


# Indian River Lagoon
IRL_dat <- predictor_data[predictor_data$study == "IRL",]  # subset to IRL
IRL_dat <- IRL_dat[,-ncol(IRL_dat)]

betas_IRL <- t(combined_betas_only[4,])  # retrieve local betas

# Multiply predictors by the local betas
expected_IRL <- as.matrix(IRL_dat) %*% as.vector(betas_IRL)

# Inverse probit for predicted values
predicted_IRL <- pnorm(expected_IRL)

# Add predictions to predictors data
IRL_dat$predicted <- predicted_IRL

# Save predictions as spatial data
IRL_test <- st_as_sf(data.frame(X = IRL_full$feature_x,
                                Y = IRL_full$feature_y,
                                predicted = predicted_IRL,
                                geometry = IRL_full$geometry))
# Plot predictions
plot_IRL <- st_transform(IRL_test, crs = 3857)
IRL_bb <- st_bbox(c(xmin = -9008088-15000,
                    ymin = 3315585-15000,
                    xmax = -8976838+15000,
                    ymax = 3385010+15000), crs = 3857)
IRL_coords <- as.data.frame(st_coordinates(plot_IRL))
plot_IRL2 <- cbind(plot_IRL[,-c(1,2)], IRL_coords)

IRL_map <- basemaps::basemap_ggplot(ext = IRL_bb) +
  geom_sf(data = plot_IRL2, aes(color = predicted), size = 0.5) +
  labs(x = "Latitude", y = "Longitude", title = "Indian River Lagoon") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.5,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = IRL_map, file = "../QGIS/IRL_predicted_map_localB.png", width = 6, height = 6)



# Add predictions to each site
# # Choctawatchee
# choc_available_columns <- coo[coo %in% colnames(choc)]
# choc_coo <- choc %>% select(all_of(choc_available_columns))
# choc_pred <- pred_dat[pred_dat$study == "choc",]
# choc_pred$predicted <- pred_dat[pred_dat$study == "choc", "predicted"]
# choc_cols <- intersect(colnames(choc_pred), colnames(choc_coo))
# choc_pred <- choc_pred %>%
#   select(!all_of(choc_cols))
#
# choc_predicted <- cbind(choc_coo, choc_pred)








## AVERAGE BETAS ## ----

average_betas <- data.frame(t(colMeans(combined_betas_only, na.rm = TRUE)))
# View(average_betas)
# str(average_betas)
# max_index <- which.max(average_betas)
# max_predictor <- predictor_columns[max_index]
# max_predictor

# Retrieve studies for predictions
pred <- readRDS("data/predictors_kriged_standardized.rds")
colnames(pred)[41] <- "predicted_OLD"  # rename to compare to new predictions
predictor_data <- as.data.frame(pred)
predictors <- meta_analysis$data$predictor[1:248]  # retrieve predictors

# Remove non-predictor columns
predictor_data <- predictor_data %>%
  select(all_of(predictors)) %>%
  select(order(colnames(.)))
# Reformat predictor columns
predictor_data <- predictor_data %>%
  mutate(across(all_of(predictors), as.character)) %>% # convert predictors to character first (needed for factorized columns)
  mutate(across(all_of(predictors), as.numeric)) %>% # convert all predictors to numeric
  mutate(across(all_of(predictors), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))  # sub NAs for statewide means
## NAs in the data can result in predicted values to be NA

# Calculate expected outcome
# expected_outcome <- predictor_matrix %*% average_betas_vector
expected_outcome <- as.matrix(predictor_data) %*% average_betas_vector

# Inverse probit for predicted values
predicted <- pnorm(expected_outcome)

# Add predictions to predictors data
pred$predicted <- predicted



# Load in original predictions and data for each site
choc <- st_transform(st_read("../Data/choctawatchee_bay/choctawatchee_bay_lssm_POINTS_0.001deg.shp"), crs = 6346) # Response:
pens <- st_transform(st_read("../Data/santa_rosa_bay/Santa_Rosa_Bay_Living_Shoreline_POINTS_0.001deg.shp"), crs = 6346) # Response:
tampa <- st_transform(st_read("../Data/tampa_bay/Tampa_Bay_Living_Shoreline_Suitability_Model_Results_POINTS_0.001deg.shp"), crs = 6346) # Response:
IRL <- st_transform(st_read("../Data/indian_river_lagoon/UCF_livingshorelinemodels_MosquitoNorthIRL_111m.shp"), crs = 6346) # Response:

# Full data to fix transformed geometry
choc_full <- st_transform(st_read("data/choctawatchee_bay/choctawatchee_bay_lssm_POINTS_0.001deg.shp"))
IRL_full <- st_transform(st_read("../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted.shp"))
pens_full <- st_transform(st_read("../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted.shp"))
tampa_full <- st_transform(st_read("../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted.shp"))

coo <- c("OBJECTID", "ID", "geometry", "feature_x", "feature_y", "nearest_x", "nearest_y", "shape__len", "Shape__Len", "distance", "distance_2", "n", "x", "y", "X", "Y")
pred_dat <- as.data.frame(pred)

# Add predictions to each site
## Choctawatchee
choc_available_columns <- coo[coo %in% colnames(choc)]
choc_coo <- choc %>% select(all_of(choc_available_columns))
choc_pred <- pred_dat[pred_dat$study == "choc",]
choc_pred$predicted <- pred_dat[pred_dat$study == "choc", "predicted"]
choc_cols <- intersect(colnames(choc_pred), colnames(choc_coo))
choc_pred <- choc_pred %>%
  select(!all_of(choc_cols))

choc_predicted <- cbind(choc_coo, choc_pred)
# st_write(choc_predicted, "../output/Final_Shapefile_all_data/Choctawatchee Bay/choc_predicted_new.shp", append = FALSE)
choc_test <- choc_predicted
st_geometry(choc_test) <- choc_full$geometry
choc_test <- choc_test %>%
  select(all_of(c("feature_x", "feature_y", "predicted")))  # automatically selects geometry
colnames(choc_test)[c(1,2)] <- c("X", "Y")
# st_write(choc_test, "../output/Final_Shapefile_all_data/Choctwatchee Bay/choc_predicted_shape.shp")

plot_choc <- st_transform(choc_test, crs = 3857)   # st_transform(st_read("../output/Final_Shapefile_all_data/Choctawatchee Bay/choc_predicted_shape.shp"), crs = 3857)
choc_bb <- st_bbox(c(xmin = -9644651-15000,
                     ymin = 3551143-15000,
                     xmax = -9584357+15000,
                     ymax = 3570713+15000), crs = 3857)
choc_coords <- as.data.frame(st_coordinates(plot_choc))
plot_choc2 <- cbind(plot_choc[,-c(1,2)], choc_coords)

choc_map <- basemaps::basemap_ggplot(ext = choc_bb) +  # ESRI World Imagery basemap
  geom_sf(data = plot_choc2, aes(color = predicted), size = 0.5) +
  labs(x = "Latitude", y = "Longitude", title = "Choctawatchee Bay") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.5,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = choc_map, file = "../QGIS/choc_predicted_map_avgB.png", width = 6, height = 6)


## Pensacola Bay
pens_available_columns <- coo[coo %in% colnames(pens)]
pens_coo <- pens %>% select(all_of(pens_available_columns))
pens_pred <- pred_dat[pred_dat$study == "pens",]
pens_pred$predicted <- pred_dat[pred_dat$study == "pens", "predicted"]
pens_cols <- intersect(colnames(pens_pred), colnames(pens_coo))
pens_pred <- pens_pred %>%
  select(!all_of(pens_cols))

pens_predicted <- cbind(pens_coo, pens_pred)
# write.csv(pens_predicted, "../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted_new.csv", row.names = FALSE)
pens_test <- pens_predicted
st_geometry(pens_test) <- pens_full$geometry
pens_test <- pens_test %>%
  select(all_of(c("feature_x", "feature_y", "predicted")))  # automatically selects geometry
colnames(pens_test)[c(1,2)] <- c("X", "Y")
# write.csv(pens_test, "../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted_test.csv", row.names = FALSE)
# st_write(pens_test, "../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted_shape.shp")


# plot_pens <- st_read("../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted_shape.shp")
# plot_pens <- st_transform(st_read("../output/Final_Shapefile_all_data/Pensacola Bay/pens_predicted_shape.shp"), crs = 3857)
# plot_pens <- st_crop(plot_pens, c(xmin = -9751391, ymin = 3536051, xmax = -9640820, ymax = 3594841))
pens_bb <- st_bbox(c(xmin = -9732366-13000,
                     ymin = 3542349-13000,
                     xmax = -9659068+13000,
                     ymax = 3589286+13000), crs = 3857)
# plot_pens2 <-st_crop(plot_pens, pens_bb)
# plot_pens2 <- st_intersection(plot_pens, pens_bb)
pens_coords <- as.data.frame(st_coordinates(plot_pens))
plot_pens2 <- cbind(plot_pens[,-c(1,2)], pens_coords)

# ggplot(plot_pens2, aes(x = X, y = Y)) +
#   geom_sf(aes(color = predicted)) +
#   basemaps::basemap_gglayer(ext = pens_bb)
# ggplot(plot_pens2, aes(x = X, y = Y)) +
#   geom_sf(aes(color = predicted))

pens_map <- basemaps::basemap_ggplot(ext = pens_bb) +
  geom_sf(data = plot_pens2, aes(color = predicted), size = 0.5) +
  labs(x = "Latitude", y = "Longitude", title = "Pensacola Bay") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.5,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = pens_map, file = "../QGIS/pens_predicted_map_avgB.png", width = 6, height = 6)


# Tampa Bay
tampa_available_columns <- coo[coo %in% colnames(tampa)]
tampa_coo <- tampa %>% select(all_of(tampa_available_columns))
tampa_pred <- pred_dat[pred_dat$study == "tampa",]
tampa_pred$predicted <- pred_dat[pred_dat$study == "tampa", "predicted"]
tampa_cols <- intersect(colnames(tampa_pred), colnames(tampa_coo))
tampa_pred <- tampa_pred %>%
  select(!all_of(tampa_cols))

tampa_predicted <- cbind(tampa_coo, tampa_pred)
# write.csv(tampa_predicted, "../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted_new.csv", row.names = FALSE)
tampa_test <- tampa_predicted
st_geometry(tampa_test) <- tampa_full$geometry
# write.csv(tampa_test, "../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted_test.csv", row.names = FALSE)
tampa_test <- tampa_test %>%
  select(all_of(c("feature_x", "feature_y", "predicted")))  # automatically selects geometry
colnames(tampa_test)[c(1,2)] <- c("X", "Y")
# st_write(tampa_test, "../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted_shape.shp")

plot_tampa <- st_transform(tampa_test, crs = 3857)   # st_transform(st_read("../output/Final_Shapefile_all_data/Tampa Bay/tampa_predicted_shape.shp"), crs = 3857)
tampa_bb <- st_bbox(c(xmin = -9222837-15000,
                      ymin = 3178783-15000,
                      xmax = -9161061+15000,
                      ymax = 3253985+15000), crs = 3857)
tampa_coords <- as.data.frame(st_coordinates(plot_tampa))
plot_tampa2 <- cbind(plot_tampa[,-c(1,2)], tampa_coords)

tampa_map <- basemaps::basemap_ggplot(ext = tampa_bb) +
  geom_sf(data = plot_tampa2, aes(color = predicted), size = 0.5) +
  labs(x = "Latitude", y = "Longitude", title = "Tampa Bay") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.5,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = tampa_map, file = "../QGIS/tampa_predicted_map_avgB.png", width = 6, height = 6)


# IRL
IRL_available_columns <- coo[coo %in% colnames(IRL)]
IRL_coo <- IRL %>% select(all_of(IRL_available_columns))
IRL_pred <- pred_dat[pred_dat$study == "IRL",]
IRL_pred$predicted <- pred_dat[pred_dat$study == "IRL", "predicted"]
IRL_cols <- intersect(colnames(IRL_pred), colnames(IRL_coo))
IRL_pred <- IRL_pred %>%
  select(!all_of(IRL_cols))

IRL_predicted <- cbind(IRL_coo, IRL_pred)
# write.csv(IRL_predicted, "../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted_new.csv", row.names = FALSE)
IRL_test <- IRL_predicted
st_geometry(IRL_test) <- IRL_full$geometry
# write.csv(IRL_test, "../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted_test.csv", row.names = FALSE)
IRL_test <- IRL_test %>%
  select(all_of(c("X", "Y", "predicted")))  # automatically selects geometry
# st_write(IRL_test, "../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted_shape.shp")

plot_IRL <- st_transform(IRL_test, crs = 3857)   # st_transform(st_read("../output/Final_Shapefile_all_data/Indian River Lagoon/IRL_predicted_shape.shp"), crs = 3857)
IRL_bb <- st_bbox(c(xmin = -9008088-15000,
                    ymin = 3315585-15000,
                    xmax = -8976838+15000,
                    ymax = 3385010+15000), crs = 3857)
IRL_coords <- as.data.frame(st_coordinates(plot_IRL))
plot_IRL2 <- cbind(plot_IRL[,-c(1,2)], IRL_coords)

IRL_map <- basemaps::basemap_ggplot(ext = IRL_bb) +
  geom_sf(data = plot_IRL2, aes(color = predicted), size = 0.5) +
  labs(x = "Latitude", y = "Longitude", title = "Indian River Lagoon") +
  scale_color_gradient2(name = "Predicted", midpoint = 0.5,
                        low = "#B31B1B", mid = "#FCF75E", high = "#0BDA51") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank())
ggsave(plot = IRL_map, file = "../QGIS/IRL_predicted_map_avgB.png", width = 6, height = 6)






############### TESTING #############
# coo <- c("OBJECTID", "ID", "geometry", "feature_x", "feature_y", "nearest_x", "nearest_y", "shape__len", "Shape__Len", "distance", "distance_2", "n", "x", "y", "X", "Y")
# state$predicted <- expected_outcome
#
# ##### TAMPA #####
# tampa_test <- subset(state, study == "tampa")
# tampa <- st_transform(st_read("../Data/tampa_bay/Tampa_Bay_Living_Shoreline_Suitability_Model_Results_POINTS_0.001deg.shp"), crs = 6346) # Response:
# tampa_available_columns <- coo[coo %in% colnames(tampa)]
# tampa_coo <- tampa %>% select(all_of(tampa_available_columns))
# tampa_predicted <- cbind(tampa_coo, tampa_test)
# write.csv(tampa_predicted, "../output/tampa_predicted.csv", row.names = FALSE)
#
# ##### CHOC #####
# choc_test <- subset(state, study == "choc")
# choc <- st_transform(st_read("../Data/choctawatchee_bay/choctawatchee_bay_lssm_POINTS_0.001deg.shp"), crs = 6346) # Response:
# choc_available_columns <- coo[coo %in% colnames(choc)]
# choc_coo <- choc %>% select(all_of(choc_available_columns))
# choc_predicted <- cbind(choc_coo, choc_test)
# write.csv(choc_predicted, "../output/choc_predicted.csv", row.names = FALSE)
#
# ##### PENS #####
# pens_test <- subset(state, study == "pens")
# pens <- st_transform(st_read("../Data/santa_rosa_bay/Santa_Rosa_Bay_Living_Shoreline_POINTS_0.001deg.shp"), crs = 6346) # Response:
# pens_available_columns <- coo[coo %in% colnames(pens)]
# pens_coo <- pens %>% select(all_of(pens_available_columns))
# pens_predicted <- cbind(pens_coo, pens_test)
# write.csv(pens_predicted, "../output/pens_predicted.csv", row.names = FALSE)
#
# ##### IRL #####
# IRL_test <- subset(state, study == "IRL")
# IRL <- st_transform(st_read("../Data/indian_river_lagoon/UCF_livingshorelinemodels_MosquitoNorthIRL_111m.shp"), crs = 6346) # Response:
# IRL_available_columns <- coo[coo %in% colnames(IRL)]
# IRL_coo <- IRL %>% select(all_of(IRL_available_columns))
# IRL_predicted <- cbind(IRL_coo, IRL_test)
# IRL_predicted$feature_x <- IRL$feature_x
# IRL_predicted$feature_y <- IRL$feature_y
# # note that coordinates need to be adjusted at some point before mapping (either in R or GIS)
# ## note that exports at Florida 17N
# write.csv(IRL_predicted, "../output/Final_Shapefile_all_data/IRL_predicted.csv", row.names = FALSE)
#
# # ##### LSSM4 #####
# # # All 4 LSSMs together
# # # LSSM4_predicted <- rbind(choc_predicted, pens_predicted, tampa_predicted, IRL_predicted)
# # LSSM4_predicted <- bind_rows(choc_predicted, pens_predicted, tampa_predicted, IRL_predicted)
# # write.csv(LSSM4_predicted, "output/LSSM4_predicted.csv", row.names = FALSE)
#
