pacman::p_load(tidyverse, tidymodels, readr, skimr, knitr,
               dplyr, janitor, carplot, vip, forcats, recipes, 
               parsnips, glmnet, tinytex, ggplot2, ggplot,
               lme4, lmerTest, glmmTMB, nlme, betareg, haven, 
               randomForest, caret, DHARMa, RColorBrewer,MuMIn)

#Loading the dataset
project_data <- read.delim("MODS_Rohrlach_EAGER.tsv")

#Processing Categorical Predictor 'sample'

# Define unique strings
sample_unique <- unique(project_data$sample)

# Convert strings to numeric using match
project_data$sample_new <- match(project_data$sample, sample_unique)

#normalizing some numeric variables
project_data$nr_of_input_reads <- scale(project_data$nr_of_input_reads)
project_data$covered_sn_ps_on_1240k <- scale(project_data$covered_sn_ps_on_1240k)
project_data$nr_of_mapped_reads_over_30bp <- scale(project_data$nr_of_mapped_reads_over_30bp)
project_data$mt_to_nuclear_read_ratio <- scale(project_data$mt_to_nuclear_read_ratio)
project_data$nr_of_unique_mapped_reads <- scale(project_data$nr_of_unique_mapped_reads)
project_data$mean_fold_coverage <- scale(project_data$mean_fold_coverage)
project_data$nr_sn_ps_used_in_contamination_estimation <- scale(project_data$nr_sn_ps_used_in_contamination_estimation)

## VARIABLE NAME: covered_sn_ps_on_1240k ##

#now after doing this, numeric = sample
ggplot(data  = project_data,
       aes(x = nr_of_input_reads,
           y = covered_sn_ps_on_1240k)) +
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter") + # to add some random noise for plotting purposes
  theme_bw() +
  labs(title = "Covered SNPs on 1240k vs. No. of Input Reads",
       x = "No. of Input Reads",
       y = "Covered SNPs on 1240k") +
  theme(plot.title = element_text(size = 18), # Adjust the size as needed
        plot.subtitle = element_text(size = 14)) # Adjust the size as needed

ggplot(data = project_data, 
     aes(x = nr_of_input_reads,
         y = covered_sn_ps_on_1240k, 
         col = as.factor(source))) +
  geom_point(size = 1, 
           alpha = .7, 
           position = "jitter") +
  geom_smooth(method = lm,
            se = TRUE, 
            size = 1.5, 
            linetype = 1, 
            alpha = .7) +
  theme_bw() +
  labs(title = "'Covered SNPs on 1240k' and 'No. of Input Reads'", 
     subtitle = "Linear Relationship",
     x = "No. of Input Reads",
     y = "Covered SNPs on 1240k") +
  scale_color_manual(name = "Sources",
                   labels = c("Twist", "1240K"),
                   values = c("blue", "orange")) +
  theme(plot.title = element_text(size = 18, face = "bold"), 
        plot.subtitle = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),  
        axis.title.y = element_text(size = 14, face = "bold"))

ggplot(data    = project_data,
       aes(x   = nr_of_input_reads,
           y   = covered_sn_ps_on_1240k,
           col = sample_new))+ #to add the colours for different samples
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(25))+
  labs(title    = "'Covered SNPs on 1240k' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different samples",
       x = "No. of Input Reads",
       y = "Covered SNPs on 1240k")

ggplot(data      = project_data,
       aes(x     = nr_of_input_reads,
           y     = covered_sn_ps_on_1240k,
           col   = sample,
           group = sample))+ #to add the colours for different classes
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "right")+
  #scale_color_gradientn(colours = rainbow(25))+
  geom_smooth(method = lm,
              se     = FALSE,
              size   = .5, 
              alpha  = .8)+ # to add regression line
  labs(title    = "'Covered SNPs on 1240k' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different sample and regression lines",
       x = "No. of Input Reads",
       y = "Covered SNPs on 1240k")+
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) + 
  geom_smooth(col='black',group=1,method='lm')

#### ME Model fitting ----
lm_interactive <- lm(formula = covered_sn_ps_on_1240k ~ nr_of_input_reads * source, 
                     data    = project_data)
summary(lm_interactive)

lm_additive <- lm(formula = covered_sn_ps_on_1240k ~ source + nr_of_input_reads, 
                  data    = project_data)
summary(lm_additive)

anova(lm_interactive,
      lm_additive)

#### Is sample important?
lm_interactive_sample <- lm(formula = covered_sn_ps_on_1240k ~ source + nr_of_input_reads + factor(sample_new), 
                           data    = project_data)
summary(lm_interactive_sample)

#### Mix Effects Model fitting: Random intercept----
lmer_interactive <- lmer(formula = covered_sn_ps_on_1240k ~ source * nr_of_input_reads + (1|sample_new), 
                         data    = project_data)
summary(lmer_interactive)

lmer_additive <- lmer(formula = covered_sn_ps_on_1240k ~ source + nr_of_input_reads + (1|sample_new), 
                      data    = project_data)
summary(lmer_additive)

anova(lmer_interactive,
      lmer_additive)

AIC(lmer_interactive,
    lmer_additive)

#### But is the random effect required? ----
ranova(lmer_additive)

# Calculate the R-squared values
r_squared_values.covsnp.add <- r.squaredGLMM(lmer_additive)
r_squared_values.covsnp.inter <- r.squaredGLMM(lmer_interactive)

# Print the R-squared values
print(r_squared_values.covsnp.add)
print(r_squared_values.covsnp.inter)

#### Still need to check diagnostics! ----
par(mfrow=c(2,2))
plot(fitted(lmer_additive), resid(lmer_additive, type = "pearson"))# this will create the plot
abline(0,0, col="red")
qqnorm(resid(lmer_additive),main='Fixed Effects residuals') 
qqline(resid(lmer_additive), col = "red")
qqnorm(ranef(lmer_additive)$sample_new[,1],main='Random Effect Intercept residuals')
qqline(ranef(lmer_additive)$sample_new[,1], col = "red")
par(mfrow=c(1,1))

# ANALYSIS USING RANDOM FOREST

# Convert sample and source to factors
project_data$sample <- as.factor(project_data$sample)
project_data$source <- as.factor(project_data$source)

# Set seed for reproducibility
set.seed(1871448)

# Create a random train-test split (80% train, 20% test)
train_indices <- sample(seq_len(nrow(project_data)), size = 0.8 * nrow(project_data))
train_data <- project_data[train_indices, ]
test_data <- project_data[-train_indices, ]

# Fit the Random Forest model with the encoded data
rf_model_covsnp <- randomForest(covered_sn_ps_on_1240k ~ nr_of_input_reads + source + sample, 
                                data = train_data, 
                                ntree = 500, 
                                importance = TRUE)

# Print the model summary
print(rf_model_covsnp)

# Predict on the test set
predictions_covsnp <- predict(rf_model_covsnp, newdata = test_data)

# Calculate Mean Squared Error (MSE)
mse_covsnp <- mean((test_data$covered_sn_ps_on_1240k - predictions_covsnp)^2)
print(paste("Mean Squared Error:", mse_covsnp))

# Calculate R-squared
r_squared_covsnp <- 1 - (sum((test_data$covered_sn_ps_on_1240k - predictions_covsnp)^2) / 
                           sum((test_data$covered_sn_ps_on_1240k - mean(test_data$covered_sn_ps_on_1240k))^2))
print(paste("R-squared:", r_squared_covsnp))

# Plot variable importance
varImpPlot(rf_model_covsnp)

#Hyperparameter Tuning

# Define the grid of hyperparameters
tune_grid <- expand.grid(
  mtry = c(1, 2, 3)  # Number of variables randomly sampled as candidates at each split
)

# Set up cross-validation
train_control <- trainControl(
  method = "cv",     # Cross-validation
  number = 5,        # Number of folds
  search = "grid"    # Perform a grid search
)

# Train the model using cross-validation
rf_covsnp_tuned <- train(
  covered_sn_ps_on_1240k ~ nr_of_input_reads + source + sample,
  data = train_data,
  method = "rf",                 # Random Forest
  trControl = train_control,     # Cross-validation settings
  tuneGrid = tune_grid,          # Hyperparameter grid
  ntree = 500,
  importance = TRUE
)

# Print the best model
print(rf_covsnp_tuned)

# Best hyperparameters
best_params_covsnps <- rf_covsnp_tuned$bestTune
print(best_params_covsnps)

# Predict on the test set using the best model
predictions_tuned_covsnp <- predict(rf_covsnp_tuned, newdata = test_data)

# Calculate Mean Squared Error (MSE) on test set
mse_test_covsnp <- mean((test_data$covered_sn_ps_on_1240k - predictions_tuned_covsnp)^2)
print(paste("Test Set Mean Squared Error:", mse_test_covsnp))

# Calculate R-squared on test set
r_squared_covsnps_test <- 1 - (sum((test_data$covered_sn_ps_on_1240k - predictions_tuned_covsnp)^2) / sum((test_data$covered_sn_ps_on_1240k - mean(test_data$covered_sn_ps_on_1240k))^2))
print(paste("Test Set R-squared:", r_squared_covsnps_test))
# ------------------------------------------------------------------

## VARIABLE NAME: nr_of_mapped_reads_over_30bp ##

#Having a first glance
ggplot(data  = project_data,
       aes(x = nr_of_input_reads,
           y = nr_of_mapped_reads_over_30bp)) +
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter") + # to add some random noise for plotting purposes
  theme_bw() +
  labs(title = "Mapped Reads over 30bp vs. No. of Input Reads",
       x = "No. of Input Reads",
       y = "No. of Mapped Reads over 30bp") +
  theme(plot.title = element_text(size = 18), # Adjust the size as needed
        plot.subtitle = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size =15)) # Adjust the size as needed

ggplot(data = project_data, 
       aes(x = nr_of_input_reads,
           y = nr_of_mapped_reads_over_30bp, 
           col = as.factor(source))) +
  geom_point(size = 1, 
             alpha = .7, 
             position = "jitter") +
  geom_smooth(method = lm,
              se = TRUE, 
              size = 1.5, 
              linetype = 1, 
              alpha = .7) +
  theme_bw() +
  labs(title = "'Mapped Reads over 30bp' and 'No. of Input Reads'", 
       subtitle = "Linear relationship",
       x = "No. of Input Reads",
       y = "No. of Mapped Reads over 30bp") +
  scale_color_manual(name = "Sources",
                     labels = c("Twist", "1240K"),
                     values = c("blue", "orange")) +
  theme(plot.title = element_text(size = 18, face = "bold"), # Adjust the size as needed
        plot.subtitle = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")) # Adjust the size as needed

ggplot(data    = project_data,
       aes(x   = nr_of_input_reads,
           y   = nr_of_mapped_reads_over_30bp,
           col = sample_new))+ #to add the colours for different samples
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(25))+
  labs(title    = "'Mapped Reads over 30bp' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different samples",
       x = "No. of Input Reads",
       y = "No. of Mapped Reads over 30bp")

ggplot(data      = project_data,
       aes(x     = nr_of_input_reads,
           y     = nr_of_mapped_reads_over_30bp,
           col   = sample,
           group = sample))+ #to add the colours for different classes
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "right")+
  #scale_color_gradientn(colours = rainbow(25))+
  geom_smooth(method = lm,
              se     = FALSE,
              size   = .5, 
              alpha  = .8)+ # to add regression line
  labs(title    = "'Mapped Reads over 30bp' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different sample and regression lines",
       x = "No. of Input Reads",
       y = "No.of Mapped Reads over 30bp")+
  theme(plot.title = element_text(size = 18, face="bold"), # Adjust the size as needed
        plot.subtitle = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")) + 
  geom_smooth(col='black',group=1,method='lm')

#### ME Model fitting ----
lm_mapreadsov30bp_interactive <- lm(formula = nr_of_mapped_reads_over_30bp ~ nr_of_input_reads * source, 
                     data    = project_data)
summary(lm_mapreadsov30bp_interactive)

lm_mapreadsov30bp_additive <- lm(formula = nr_of_mapped_reads_over_30bp ~ source + nr_of_input_reads, 
                  data    = project_data)
summary(lm_mapreadsov30bp_additive)

anova(lm_mapreadsov30bp_interactive,
      lm_mapreadsov30bp_additive)

AIC(lm_mapreadsov30bp_interactive,
    lm_mapreadsov30bp_additive)

#### Is sample important?
lm_mapreadsov30bp_interactive_sample <- lm(formula = nr_of_mapped_reads_over_30bp ~ source + nr_of_input_reads + factor(sample_new), 
                            data    = project_data)
summary(lm_mapreadsov30bp_interactive_sample)

#### Mix Effects Model fitting: Random intercept----
lmer_mapreadsov30bp_interactive <- lmer(formula = nr_of_mapped_reads_over_30bp ~ source * nr_of_input_reads + (1|sample_new), 
                         data    = project_data)
summary(lmer_mapreadsov30bp_interactive)

lmer_mapreadsov30bp_additive <- lmer(formula = nr_of_mapped_reads_over_30bp ~ source + nr_of_input_reads + (1|sample_new), 
                      data    = project_data)
summary(lmer_mapreadsov30bp_additive)

anova(lmer_mapreadsov30bp_interactive,
      lmer_mapreadsov30bp_additive)

AIC(lmer_mapreadsov30bp_interactive,
    lmer_mapreadsov30bp_additive)

#### But is the random effect required? ----
ranova(lmer_mapreadsov30bp_additive)

# Calculate the R-squared values
r_squared_values.add <- r.squaredGLMM(lmer_mapreadsov30bp_additive)
r_squared_values.inter <- r.squaredGLMM(lmer_mapreadsov30bp_interactive)

# Print the R-squared values
print(r_squared_values.add)
print(r_squared_values.inter)

#### Still need to check diagnostics! ----
par(mfrow=c(2,2))
plot(fitted(lmer_mapreadsov30bp_additive), resid(lmer_mapreadsov30bp_additive, type = "pearson"))# this will create the plot
abline(0,0, col="red")
qqnorm(resid(lmer_mapreadsov30bp_additive),main='Fixed Effects residuals') 
qqline(resid(lmer_mapreadsov30bp_additive), col = "red")
qqnorm(ranef(lmer_mapreadsov30bp_additive)$sample_new[,1],main='Random Effect Intercept residuals')
qqline(ranef(lmer_mapreadsov30bp_additive)$sample_new[,1], col = "red")
par(mfrow=c(1,1))

# ANALYSIS USING RANDOM FOREST

# Convert sample and source to factors
project_data$sample <- as.factor(project_data$sample)
project_data$source <- as.factor(project_data$source)

# Set seed for reproducibility
set.seed(1871448)

# Create a random train-test split (80% train, 20% test)
train_indices <- sample(seq_len(nrow(project_data)), size = 0.8 * nrow(project_data))
train_data <- project_data[train_indices, ]
test_data <- project_data[-train_indices, ]

# Fit the Random Forest model with the encoded data
rf_model_mapreadsov30bp <- randomForest(nr_of_mapped_reads_over_30bp ~ nr_of_input_reads + source + sample, 
                                data = train_data, 
                                ntree = 500, 
                                importance = TRUE)

# Print the model summary
print(rf_model_mapreadsov30bp)

# Predict on the test set
predictions_mapreadsov30bp <- predict(rf_model_mapreadsov30bp, newdata = test_data)

# Calculate Mean Squared Error (MSE)
mse_mapreadsov30bp <- mean((test_data$nr_of_mapped_reads_over_30bp - predictions_mapreadsov30bp)^2)
print(paste("Mean Squared Error:", mse_mapreadsov30bp))

# Calculate R-squared
r_squared_mapreadsov30bp <- 1 - (sum((test_data$nr_of_mapped_reads_over_30bp - predictions_mapreadsov30bp)^2) / 
                           sum((test_data$nr_of_mapped_reads_over_30bp - mean(test_data$nr_of_mapped_reads_over_30bp))^2))
print(paste("R-squared:", r_squared_mapreadsov30bp))

# Plot variable importance
varImpPlot(rf_model_mapreadsov30bp)

#Hyperparameter Tuning

# Define the grid of hyperparameters
tune_grid <- expand.grid(
  mtry = c(1, 2, 3)  # Number of variables randomly sampled as candidates at each split
)

# Set up cross-validation
train_control <- trainControl(
  method = "cv",     # Cross-validation
  number = 5,        # Number of folds
  search = "grid"    # Perform a grid search
)

# Train the model using cross-validation
rf_mapreadsov30bp_tuned <- train(
  nr_of_mapped_reads_over_30bp ~ nr_of_input_reads + source + sample,
  data = train_data,
  method = "rf",                 # Random Forest
  trControl = train_control,     # Cross-validation settings
  tuneGrid = tune_grid,          # Hyperparameter grid
  ntree = 500,
  importance = TRUE
)

# Print the best model
print(rf_mapreadsov30bp_tuned)

# Best hyperparameters
best_params_mapreadsov30bp <- rf_mapreadsov30bp_tuned$bestTune
print(best_params_mapreadsov30bp)

# Predict on the test set using the best model
predictions_tuned_mapreadsov30bp <- predict(rf_mapreadsov30bp_tuned, newdata = test_data)

# Calculate Mean Squared Error (MSE) on test set
mse_test_mapreadsov30bp <- mean((test_data$nr_of_mapped_reads_over_30bp - predictions_tuned_mapreadsov30bp)^2)
print(paste("Test Set Mean Squared Error:", mse_test_mapreadsov30bp))

# Calculate R-squared on test set 
r_squared_mapreadsov30bp_test <- 1 - (sum((test_data$nr_of_mapped_reads_over_30bp - predictions_tuned_mapreadsov30bp)^2) / sum((test_data$nr_of_mapped_reads_over_30bp - mean(test_data$nr_of_mapped_reads_over_30bp))^2))
print(paste("Test Set R-squared:", r_squared_mapreadsov30bp_test))
# ----------------------------------------------------------------------------------

## VARIABLE NAME: nr_of_unique_mapped_reads ##

#Having a first glance
ggplot(data  = project_data,
       aes(x = nr_of_input_reads,
           y = nr_of_unique_mapped_reads)) +
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter") + # to add some random noise for plotting purposes
  theme_bw() +
  labs(title = "Unique Mapped Reads vs. No. of Input Reads",
       x = "No. of Input Reads",
       y = "No. of Unique Mapped Reads") +
  theme(plot.title = element_text(size = 18), # Adjust the size as needed
        plot.subtitle = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) # Adjust the size as needed

ggplot(data = project_data, 
       aes(x = nr_of_input_reads,
           y = nr_of_unique_mapped_reads, 
           col = as.factor(source))) +
  geom_point(size = 1, 
             alpha = .7, 
             position = "jitter") +
  geom_smooth(method = lm,
              se = TRUE, 
              size = 1.5, 
              linetype = 1, 
              alpha = .7) +
  theme_bw() +
  labs(title = "'Unique Mapped Reads' and 'No. of Input Reads'", 
       subtitle = "Linear Relationship",
       x = "No. of Input Reads",
       y = "No. of Unique Mapped Reads") +
  scale_color_manual(name = "Sources",
                     labels = c("Twist", "1240K"),
                     values = c("blue", "orange")) +
  theme(plot.title = element_text(size = 18, face="bold"), # Adjust the size as needed
        plot.subtitle = element_text(size = 16, face="bold"),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold")) # Adjust the size as needed

ggplot(data    = project_data,
       aes(x   = nr_of_input_reads,
           y   = nr_of_unique_mapped_reads,
           col = sample_new))+ #to add the colours for different samples
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(25))+
  labs(title    = "'Unique Mapped Reads' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different samples",
       x = "No. of Input Reads",
       y = "Unique Mapped Reads")

ggplot(data      = project_data,
       aes(x     = nr_of_input_reads,
           y     = nr_of_unique_mapped_reads,
           col   = sample,
           group = sample))+ #to add the colours for different classes
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "right")+
  #scale_color_gradientn(colours = rainbow(25))+
  geom_smooth(method = lm,
              se     = FALSE,
              size   = .5, 
              alpha  = .8)+ # to add regression line
  labs(title    = "'Unique Mapped Reads' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different sample and regression lines",
       x = "No. of Input Reads",
       y = "Unique Mapped Reads")+
  theme(plot.title = element_text(size = 18, face = "bold"), # Adjust the size as needed
        plot.subtitle = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")) + 
  geom_smooth(col='black',group=1,method='lm')

#### ME Model fitting ----
lm_unqmapreads_interactive <- lm(formula = nr_of_unique_mapped_reads ~ nr_of_input_reads * source, 
                                    data    = project_data)
summary(lm_unqmapreads_interactive)

lm_unqmapreads_additive <- lm(formula = nr_of_unique_mapped_reads ~ source + nr_of_input_reads, 
                                 data    = project_data)
summary(lm_unqmapreads_additive)

anova(lm_unqmapreads_interactive,
      lm_unqmapreads_additive)

#### Is sample important?
lm_unqmapreads_interactive_sample <- lm(formula = nr_of_unique_mapped_reads ~ source + nr_of_input_reads + factor(sample_new), 
                                           data    = project_data)
summary(lm_unqmapreads_interactive_sample)

#### Mix Effects Model fitting: Random intercept----
lmer_unqmapreads_interactive <- lmer(formula = nr_of_unique_mapped_reads ~ source * nr_of_input_reads + (1|sample_new), 
                                        data    = project_data)
summary(lmer_unqmapreads_interactive)

lmer_unqmapreads_additive <- lmer(formula = nr_of_unique_mapped_reads ~ source + nr_of_input_reads + (1|sample_new), 
                                     data    = project_data)
summary(lmer_unqmapreads_additive)

anova(lmer_unqmapreads_interactive,
      lmer_unqmapreads_additive)
#Here result says p-value is not significant, so both models explain the relationship equally well. So, we go for the simpler model which is the additive model.

AIC(lmer_unqmapreads_additive,
    lmer_unqmapreads_interactive)


#### But is the random effect required? ----
ranova(lmer_unqmapreads_additive)

# Calculate the R-squared values
r_squared_values.unqmapreads.add <- r.squaredGLMM(lmer_unqmapreads_interactive)
r_squared_values.unqmapreads.inter <- r.squaredGLMM(lmer_unqmapreads_interactive)

# Print the R-squared values
print(r_squared_values.unqmapreads.add)
print(r_squared_values.unqmapreads.inter)

# ________________________________________________________________________________

#### Still need to check diagnostics! ----
par(mfrow=c(2,2))
plot(fitted(lmer_unqmapreads_additive), resid(lmer_unqmapreads_additive, type = "pearson"))# this will create the plot
abline(0,0, col="red")
qqnorm(resid(lmer_unqmapreads_additive),main='Fixed Effects residuals') 
qqline(resid(lmer_unqmapreads_additive), col = "red")
qqnorm(ranef(lmer_unqmapreads_additive)$sample_new[,1],main='Random Effect Intercept residuals')
qqline(ranef(lmer_unqmapreads_additive)$sample_new[,1], col = "red")
par(mfrow=c(1,1))

# ANALYSIS USING RANDOM FOREST

# Fit the Random Forest model with the encoded data
rf_model_unqmapreads <- randomForest(nr_of_unique_mapped_reads ~ nr_of_input_reads + source + sample, 
                                        data = train_data, 
                                        ntree = 500, 
                                        importance = TRUE)

# Print the model summary
print(rf_model_unqmapreads)

# Predict on the test set
predictions_unqmapreads <- predict(rf_model_unqmapreads, newdata = test_data)

# Calculate Mean Squared Error (MSE)
mse_unqmapreads <- mean((test_data$nr_of_unique_mapped_reads - predictions_unqmapreads)^2)
print(paste("Mean Squared Error:", mse_unqmapreads))

# Calculate R-squared
r_squared_unqmapreads <- 1 - (sum((test_data$nr_of_unique_mapped_reads - predictions_unqmapreads)^2) / 
                                   sum((test_data$nr_of_unique_mapped_reads - mean(test_data$nr_of_unique_mapped_reads))^2))
print(paste("R-squared:", r_squared_unqmapreads))

# Plot variable importance
varImpPlot(rf_model_unqmapreads)

#Hyperparameter Tuning

# Define the grid of hyperparameters
tune_grid <- expand.grid(
  mtry = c(1, 2, 3)  # Number of variables randomly sampled as candidates at each split
)

# Set up cross-validation
train_control <- trainControl(
  method = "cv",     # Cross-validation
  number = 5,        # Number of folds
  search = "grid"    # Perform a grid search
)

# Train the model using cross-validation
rf_unqmapreads_tuned <- train(
  nr_of_unique_mapped_reads ~ nr_of_input_reads + source + sample,
  data = train_data,
  method = "rf",                 # Random Forest
  trControl = train_control,     # Cross-validation settings
  tuneGrid = tune_grid,          # Hyperparameter grid
  ntree = 500,
  importance = TRUE
)

# Print the best model
print(rf_unqmapreads_tuned)

# Best hyperparameters
best_params_unqmapreads <- rf_unqmapreads_tuned$bestTune
print(best_params_unqmapreads)

# Predict on the test set using the best model
predictions_tuned_unqmapreads <- predict(rf_unqmapreads_tuned, newdata = test_data)

# Calculate Mean Squared Error (MSE) on test set
mse_test_unqmapreads <- mean((test_data$nr_of_unique_mapped_reads - predictions_tuned_unqmapreads)^2)
print(paste("Test Set Mean Squared Error:", mse_test_unqmapreads))

# Calculate R-squared on test set
r_squared_unqmapreads_test <- 1 - (sum((test_data$nr_of_unique_mapped_reads - predictions_tuned_unqmapreads)^2) / sum((test_data$nr_of_unique_mapped_reads - mean(test_data$nr_of_unique_mapped_reads))^2))
print(paste("Test Set R-squared:", r_squared_unqmapreads_test))

#-----------------------------------------------------------------------------------

## VARIABLE NAME: mean_fold_coverage ##

#Having a first glance
ggplot(data  = project_data,
       aes(x = nr_of_input_reads,
           y = mean_fold_coverage)) +
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter") + # to add some random noise for plotting purposes
  theme_bw() +
  labs(title = "Mean Fold Coverage vs. No. of Input Reads",
       x = "No. of Input Reads",
       y = "Mean Fold Coverage") +
  theme(plot.title = element_text(size = 18), # Adjust the size as needed
        plot.subtitle = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) # Adjust the size as needed

ggplot(data = project_data, 
       aes(x = nr_of_input_reads,
           y = mean_fold_coverage, 
           col = as.factor(source))) +
  geom_point(size = 1, 
             alpha = .7, 
             position = "jitter") +
  geom_smooth(method = lm,
              se = TRUE, 
              size = 1.5, 
              linetype = 1, 
              alpha = .7) +
  theme_bw() +
  labs(title = "'Mean Fold Coverage' and 'No. of Input Reads'", 
       subtitle = "Linear relationship",
       x = "No. of Input Reads",
       y = "Mean Fold Coverage") +
  scale_color_manual(name = "Sources",
                     labels = c("Twist", "1240K"),
                     values = c("blue", "orange")) +
  theme(plot.title = element_text(size = 18, face = "bold"), # Adjust the size as needed
        plot.subtitle = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold")) # Adjust the size as needed

ggplot(data    = project_data,
       aes(x   = nr_of_input_reads,
           y   = mean_fold_coverage,
           col = sample))+ #to add the colours for different samples
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(25))+
  labs(title    = "'Mean Fold Coverage' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different samples",
       x = "No. of Input Reads",
       y = "Mean Fold Coverage")

#project_data$sample_new <- as.character(project_data$sample_new)
#project_data$source <- as.character(project_data$source)

ggplot(data      = project_data,
       aes(x     = nr_of_input_reads,
           y     = mean_fold_coverage,
           col   = sample,
           group = sample))+ #to add the colours for different classes
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "right")+
  #scale_color_gradientn(colours = rainbow(25))+
  geom_smooth(method = lm,
              se     = FALSE,
              size   = .5, 
              alpha  = .8)+ # to add regression line
  labs(title    = "'Mean Fold Coverage' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different sample and regression lines",
       x = "No. of Input Reads",
       y = "Mean Fold Coverage")+
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold")) +
geom_smooth(col='black',group=1,method='lm')


#### ME Model fitting ----
lm_foldcov_interactive <- lm(formula = mean_fold_coverage ~ nr_of_input_reads * source, 
                                 data    = project_data)
summary(lm_foldcov_interactive)

lm_foldcov_additive <- lm(formula = mean_fold_coverage ~ source + nr_of_input_reads, 
                              data    = project_data)
summary(lm_foldcov_additive)

anova(lm_foldcov_interactive,
      lm_foldcov_additive)

#### Is sample important?
lm_foldcov_interactive_sample <- lm(formula = mean_fold_coverage ~ source + nr_of_input_reads + factor(sample_new), 
                                        data    = project_data)
summary(lm_foldcov_interactive_sample)

#### Mix Effects Model fitting: Random intercept----
lmer_foldcov_interactive <- lmer(formula = mean_fold_coverage ~ source * nr_of_input_reads + (1|sample_new), 
                                     data    = project_data)
summary(lmer_foldcov_interactive)

lmer_foldcov_additive <- lmer(formula = mean_fold_coverage ~ source + nr_of_input_reads + (1|sample_new), 
                                  data    = project_data)
summary(lmer_foldcov_additive)

anova(lmer_foldcov_interactive,
      lmer_foldcov_additive)

AIC(lmer_foldcov_interactive,
      lmer_foldcov_additive)

#### But is the random effect required? ----
ranova(lmer_foldcov_additive)

# Calculate the R-squared values
r_squared_values.foldcov.add <- r.squaredGLMM(lmer_foldcov_additive)
r_squared_values.foldcov.inter <- r.squaredGLMM(lmer_foldcov_interactive)

# Print the R-squared values
print(r_squared_values.foldcov.add)
print(r_squared_values.foldcov.inter)
# ________________________________________________________________________________

#### Still need to check diagnostics! ----
par(mfrow=c(2,2))
plot(fitted(lmer_foldcov_additive), resid(lmer_foldcov_additive, type = "pearson"))# this will create the plot
abline(0,0, col="red")
qqnorm(resid(lmer_foldcov_additive),main='Fixed Effects residuals') 
qqline(resid(lmer_foldcov_additive), col = "red")
qqnorm(ranef(lmer_foldcov_additive)$sample_new[,1],main='Random Effect Intercept residuals')
qqline(ranef(lmer_foldcov_additive)$sample_new[,1], col = "red")
par(mfrow=c(1,1))

# ANALYSIS USING RANDOM FOREST

# Fit the Random Forest model with the encoded data
rf_model_foldcov <- randomForest(mean_fold_coverage ~ nr_of_input_reads + source + sample, 
                                     data = train_data, 
                                     ntree = 500, 
                                     importance = TRUE)

# Print the model summary
print(rf_model_foldcov)

# Predict on the test set
predictions_foldcov <- predict(rf_model_foldcov, newdata = test_data)

# Calculate Mean Squared Error (MSE)
mse_foldcov <- mean((test_data$mean_fold_coverage - predictions_foldcov)^2)
print(paste("Mean Squared Error:", mse_foldcov))

# Calculate R-squared
r_squared_foldcov <- 1 - (sum((test_data$mean_fold_coverage - predictions_foldcov)^2) / 
                                sum((test_data$mean_fold_coverage - mean(test_data$mean_fold_coverage))^2))
print(paste("R-squared:", r_squared_foldcov))

# Plot variable importance
varImpPlot(rf_model_foldcov)

#Hyperparameter Tuning

# Define the grid of hyperparameters
tune_grid <- expand.grid(
  mtry = c(1, 2, 3)  # Number of variables randomly sampled as candidates at each split
)

# Set up cross-validation
train_control <- trainControl(
  method = "cv",     # Cross-validation
  number = 5,        # Number of folds
  search = "grid"    # Perform a grid search
)

# Train the model using cross-validation
rf_foldcov_tuned <- train(
  mean_fold_coverage ~ nr_of_input_reads + source + sample,
  data = train_data,
  method = "rf",                 # Random Forest
  trControl = train_control,     # Cross-validation settings
  tuneGrid = tune_grid,          # Hyperparameter grid
  ntree = 500,
  importance = TRUE
)

# Print the best model
print(rf_foldcov_tuned)

# Best hyperparameters
best_params_foldcov <- rf_foldcov_tuned$bestTune
print(best_params_foldcov)

# Predict on the test set using the best model
predictions_tuned_foldcov <- predict(rf_foldcov_tuned, newdata = test_data)

# Calculate Mean Squared Error (MSE) on test set
mse_test_foldcov <- mean((test_data$mean_fold_coverage - predictions_tuned_foldcov)^2)
print(paste("Test Set Mean Squared Error:", mse_test_foldcov))

# Calculate R-squared on test set
r_squared_foldcov_test <- 1 - (sum((test_data$mean_fold_coverage - predictions_tuned_foldcov)^2) / sum((test_data$mean_fold_coverage - mean(test_data$mean_fold_coverage))^2))
print(paste("Test Set R-squared:", r_squared_foldcov_test))
#-----------------------------------------------------------------------------------

## VARIABLE NAME: mt_to_nuclear_read_ratio ##

#Having a first glance
ggplot(data  = project_data,
       aes(x = nr_of_input_reads,
           y = mt_to_nuclear_read_ratio)) +
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter") + # to add some random noise for plotting purposes
  theme_bw() +
  labs(title = "'Mitochondrial-to-Nuclear Read Ratio' vs. 'No. of Input Reads'",
       x = "No. of Input Reads",
       y = "Mitochondrial-to-Nuclear Read Ratio") +
  theme(plot.title = element_text(size = 18), # Adjust the size as needed
        plot.subtitle = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) # Adjust the size as needed

ggplot(data = project_data, 
       aes(x = nr_of_input_reads,
           y = mt_to_nuclear_read_ratio, 
           col = as.factor(source))) +
  geom_point(size = 1, 
             alpha = .7, 
             position = "jitter") +
  geom_smooth(method = lm,
              se = TRUE, 
              size = 1.5, 
              linetype = 1, 
              alpha = .7) +
  theme_bw() +
  labs(title = "'Mitochondria-to-Nuclear Read Ratio' vs 'No. of Input Reads'", 
       subtitle = "Linear Relationship",
       x = "No. of Input Reads",
       y = "Mitochondrial-to-Nuclear Read Ratio") +
  scale_color_manual(name = "Sources",
                     labels = c("Twist", "1240K"),
                     values = c("blue", "orange")) +
  theme(plot.title = element_text(size = 18, face="bold"), # Adjust the size as needed
        plot.subtitle = element_text(size = 16, face="bold"),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold")) # Adjust the size as needed

ggplot(data    = project_data,
       aes(x   = nr_of_input_reads,
           y   = mt_to_nuclear_read_ratio,
           col = sample_new))+ #to add the colours for different samples
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(25))+
  labs(title    = "'Mitochondria-to-Nuclear Read Ratio' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different samples",
       x = "No. of Input Reads",
       y = "Mitochondrial-to-Nuclear Read Ratio")

ggplot(data      = project_data,
       aes(x     = nr_of_input_reads,
           y     = mt_to_nuclear_read_ratio,
           col   = sample,
           group = sample))+ #to add the colours for different classes
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "right")+
  #scale_color_gradientn(colours = rainbow(25))+
  geom_smooth(method = lm,
              se     = FALSE,
              size   = .5, 
              alpha  = .8)+ # to add regression line
  labs(title    = "'Mitochondrial-to-Nuclear Read Ratio' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different sample and regression lines",
       x = "No. of Input Reads",
       y = "Mitochondrial-to-Nuclear Read Ratio")+
  theme(plot.title = element_text(size = 18, face="bold"),
        plot.subtitle = element_text(size = 16, face="bold"),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold")) + 
  geom_smooth(col='black',group=1,method='lm')

#### ME Model fitting ----
lm_mtnuclread_interactive <- lm(formula = mt_to_nuclear_read_ratio ~ nr_of_input_reads * source, 
                             data    = project_data)
summary(lm_mtnuclread_interactive)

lm_mtnuclread_additive <- lm(formula = mt_to_nuclear_read_ratio ~ source + nr_of_input_reads, 
                          data    = project_data)
summary(lm_mtnuclread_additive)

lm_mtnuclread_additive_simpler <- lm(formula = mt_to_nuclear_read_ratio ~ source, 
                             data    = project_data)
summary(lm_mtnuclread_additive_simpler)


anova(lm_mtnuclread_interactive,
      lm_mtnuclread_additive,
      lm_mtnuclread_additive_simpler)

#### Is sample important?
lm_mtnuclread_interactive_sample <- lm(formula = mt_to_nuclear_read_ratio ~ source + factor(sample_new), 
                                    data    = project_data)
summary(lm_mtnuclread_interactive_sample)

#### Mix Effects Model fitting: Random intercept----
lmer_mtnuclread_interactive <- lmer(formula = mt_to_nuclear_read_ratio ~ source * nr_of_input_reads + (1|sample_new), 
                                 data    = project_data)
summary(lmer_mtnuclread_interactive)

lmer_mtnuclread_additive <- lmer(formula = mt_to_nuclear_read_ratio ~ source + nr_of_input_reads + (1|sample_new), 
                              data    = project_data)
summary(lmer_mtnuclread_additive)

lmer_mtnuclread_additive_simpler <- lmer(formula = mt_to_nuclear_read_ratio ~ source + (1|sample_new), 
                              data    = project_data)
summary(lmer_mtnuclread_additive_simpler)

anova(lmer_mtnuclread_interactive,
      lmer_mtnuclread_additive,
      lmer_mtnuclread_additive_simpler)

AIC(lmer_mtnuclread_interactive,
    lmer_mtnuclread_additive_simpler)

#### But is the random effect required? ----
ranova(lmer_mtnuclread_additive_simpler)

#sample had no impact and it makes sense because. mitchondrial DNA carries no information that could be significant in the contex to DNA sequencing.
#so, lm results needed to be added sometime

# Calculate the R-squared values
r_squared_values.mtnuclread.add.simpler <- r.squaredGLMM(lmer_mtnuclread_additive_simpler)
r_squared_values.mtnuclread.add <- r.squaredGLMM(lmer_mtnuclread_additive)
r_squared_values.mtnuclread.inter <- r.squaredGLMM(lmer_mtnuclread_interactive)

# Print the R-squared values
r_squared_values.mtnuclread.add.simpler
print(r_squared_values.mtnuclread.add)
print(r_squared_values.mtnuclread.inter)
# ________________________________________________________________________________
#### Still need to check diagnostics! ----
par(mfrow=c(2,2))
plot(fitted(lmer_mtnuclread_additive_simpler), resid(lmer_mtnuclread_additive_simpler, type = "pearson"))# this will create the plot
abline(0,0, col="red")
qqnorm(resid(lmer_mtnuclread_additive_simpler),main='Fixed Effects residuals') 
qqline(resid(lmer_mtnuclread_additive_simpler), col = "red")
qqnorm(ranef(lmer_mtnuclread_additive_simpler)$sample_new[,1],main='Random Effect Intercept residuals')
qqline(ranef(lmer_mtnuclread_additive_simpler)$sample_new[,1], col = "red")
par(mfrow=c(1,1))

# ANALYSIS USING RANDOM FOREST

# Fit the Random Forest model with the encoded data
rf_model_mtnuclread <- randomForest(mt_to_nuclear_read_ratio ~ nr_of_input_reads + source + sample, 
                                 data = train_data, 
                                 ntree = 500, 
                                 importance = TRUE)

# Print the model summary
print(rf_model_mtnuclread)

# Predict on the test set
predictions_mtnuclread <- predict(rf_model_mtnuclread, newdata = test_data)

# Calculate Mean Squared Error (MSE)
mse_mtnuclread <- mean((test_data$mt_to_nuclear_read_ratio - predictions_mtnuclread)^2)
print(paste("Mean Squared Error:", mse_mtnuclread))

# Calculate R-squared
r_squared_mtnuclread <- 1 - (sum((test_data$mt_to_nuclear_read_ratio - predictions_mtnuclread)^2) / 
                            sum((test_data$mt_to_nuclear_read_ratio - mean(test_data$mt_to_nuclear_read_ratio))^2))
print(paste("R-squared:", r_squared_mtnuclread))

# Plot variable importance
varImpPlot(rf_model_mtnuclread)

#Hyperparameter Tuning

# Define the grid of hyperparameters
tune_grid <- expand.grid(
  mtry = c(1, 2, 3)  # Number of variables randomly sampled as candidates at each split
)

# Set up cross-validation
train_control <- trainControl(
  method = "cv",     # Cross-validation
  number = 5,        # Number of folds
  search = "grid"    # Perform a grid search
)

# Train the model using cross-validation
rf_mtnuclread_tuned <- train(
  mt_to_nuclear_read_ratio ~ nr_of_input_reads + source + sample,
  data = train_data,
  method = "rf",                 # Random Forest
  trControl = train_control,     # Cross-validation settings
  tuneGrid = tune_grid,          # Hyperparameter grid
  ntree = 500,
  importance = TRUE
)

# Print the best model
print(rf_mtnuclread_tuned)

# Best hyperparameters
best_params_mtnuclread <- rf_mtnuclread_tuned$bestTune
print(best_params_mtnuclread)

# Predict on the test set using the best model
predictions_tuned_mtnuclread <- predict(rf_mtnuclread_tuned, newdata = test_data)

# Calculate Mean Squared Error (MSE) on test set
mse_test_mtnuclread <- mean((test_data$mt_to_nuclear_read_ratio - predictions_tuned_mtnuclread)^2)
print(paste("Test Set Mean Squared Error:", mse_test_mtnuclread))

# Calculate R-squared on test set
r_squared_mtnuclread_test <- 1 - (sum((test_data$mt_to_nuclear_read_ratio - predictions_tuned_mtnuclread)^2) / sum((test_data$mt_to_nuclear_read_ratio - mean(test_data$mt_to_nuclear_read_ratio))^2))
print(paste("Test Set R-squared:", r_squared_mtnuclread_test))
#-----------------------------------------------------------------------------------

## VARIABLE NAME: nr_sn_ps_used_in_contamination_estimation ##

#Having a first glance
ggplot(data  = project_data,
       aes(x = nr_of_input_reads,
           y = nr_sn_ps_used_in_contamination_estimation)) +
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter") + # to add some random noise for plotting purposes
  theme_bw() +
  labs(title = "No. of SNPs Used for Contamination Estimation vs. No. of Input Reads",
       x = "No. of Input Reads",
       y = "No. of SNPs Used for Contamination Estimation") +
  theme(plot.title = element_text(size = 18, face = "bold"), # Adjust the size as needed
        plot.subtitle = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) # Adjust the size as needed

ggplot(data = project_data, 
       aes(x = nr_of_input_reads,
           y = nr_sn_ps_used_in_contamination_estimation, 
           col = as.factor(source))) +
  geom_point(size = 1, 
             alpha = .7, 
             position = "jitter") +
  geom_smooth(method = lm,
              se = TRUE, 
              size = 1.5, 
              linetype = 1, 
              alpha = .7) +
  theme_bw() +
  labs(title = "'SNPs used in Contaminaiton Estimation' and 'No. of Input Reads'", 
       subtitle = "Linear relationship",
       x = "No. of Input Reads",
       y = "No. of SNPs Used in Contamination Estimation") +
  scale_color_manual(name = "Sources",
                     labels = c("Twist", "1240K"),
                     values = c("blue", "orange")) +
  theme(plot.title = element_text(size = 18, face = "bold"), # Adjust the size as needed
        plot.subtitle = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) # Adjust the size as needed

ggplot(data    = project_data,
       aes(x   = nr_of_input_reads,
           y   = nr_sn_ps_used_in_contamination_estimation,
           col = sample_new))+ #to add the colours for different samples
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(25))+
  labs(title    = "'No. of SNPs Used in Contamination Estimation' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different samples",
       x = "No. of Input Reads",
       y = "No. of SNPs Used in Contamination Estimation")

ggplot(data      = project_data,
       aes(x     = nr_of_input_reads,
           y     = nr_sn_ps_used_in_contamination_estimation,
           col   = sample,
           group = sample))+ #to add the colours for different classes
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_bw()+
  theme(legend.position = "right")+
  #scale_color_gradientn(colours = rainbow(25))+
  geom_smooth(method = lm,
              se     = FALSE,
              size   = .5, 
              alpha  = .8)+ # to add regression line
  labs(title    = "'No. of SNPs Used in Contamination Estimation' vs. 'No. of Input Reads'",
       subtitle = "Added colours for different sample and regression lines",
       x = "No. of Input Reads",
       y = "No. of SNPs Used in Contamination Estimation")+
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")) +
  geom_smooth(col='black',group=1,method='lm')




#### ME Model fitting ----
lm_snpsincontamestim_interactive <- lm(formula = nr_sn_ps_used_in_contamination_estimation ~ nr_of_input_reads * source, 
                                data    = project_data)
summary(lm_snpsincontamestim_interactive)

lm_snpsincontamestim_additive <- lm(formula = nr_sn_ps_used_in_contamination_estimation ~ source + nr_of_input_reads, 
                             data    = project_data)
summary(lm_snpsincontamestim_additive)

anova(lm_snpsincontamestim_interactive,
      lm_snpsincontamestim_additive)

#### Is sample important?
lm_snpsincontamestim_interactive_sample <- lm(formula = nr_sn_ps_used_in_contamination_estimation ~ source + factor(sample_new), 
                                       data    = project_data)
summary(lm_snpsincontamestim_interactive_sample)

#### Mix Effects Model fitting: Random intercept----
lmer_snpsincontamestim_interactive <- lmer(formula = nr_sn_ps_used_in_contamination_estimation ~ source * nr_of_input_reads + (1|sample_new), 
                                    data    = project_data)
summary(lmer_snpsincontamestim_interactive)

#we cannot remove input reads and keep interaction term, so we cannot keep the interaction

lmer_snpsincontamestim_additive <- lmer(formula = nr_sn_ps_used_in_contamination_estimation ~ source + nr_of_input_reads + (1|sample_new), 
                                 data    = project_data)
summary(lmer_snpsincontamestim_additive)

anova(lmer_snpsincontamestim_interactive,
      lmer_snpsincontamestim_additive
      )

AIC(lmer_snpsincontamestim_additive,
    lmer_snpsincontamestim_interactive)

#### But is the random effect required? ----
ranova(lmer_snpsincontamestim_additive)

# Calculate the R-squared values
r_squared_values.snpsincontamestim.add <- r.squaredGLMM(lmer_snpsincontamestim_additive)
r_squared_values.snpsincontamestim.inter <- r.squaredGLMM(lmer_snpsincontamestim_interactive)

# Print the R-squared values
print(r_squared_values.snpsincontamestim.add)
print(r_squared_values.snpsincontamestim.inter)

# ________________________________________________________________________________

#### Still need to check diagnostics! ----
par(mfrow=c(2,2))
plot(fitted(lmer_snpsincontamestim_additive), resid(lmer_snpsincontamestim_additive, type = "pearson"))# this will create the plot
abline(0,0, col="red")
qqnorm(resid(lmer_snpsincontamestim_additive),main='Fixed Effects residuals') 
qqline(resid(lmer_snpsincontamestim_additive), col = "red")
qqnorm(ranef(lmer_snpsincontamestim_additive)$sample_new[,1],main='Random Effect Intercept residuals')
qqline(ranef(lmer_snpsincontamestim_additive)$sample_new[,1], col = "red")
par(mfrow=c(1,1))

# ANALYSIS USING RANDOM FOREST

# Fit the Random Forest model with the encoded data
rf_model_snpsincontamestim <- randomForest(nr_sn_ps_used_in_contamination_estimation ~ nr_of_input_reads + source + sample, 
                                    data = train_data, 
                                    ntree = 500, 
                                    importance = TRUE)

# Print the model summary
print(rf_model_snpsincontamestim)

# Predict on the test set
predictions_snpsincontamestim <- predict(rf_model_snpsincontamestim, newdata = test_data)

# Calculate Mean Squared Error (MSE)
mse_snpsincontamestim <- mean((test_data$nr_sn_ps_used_in_contamination_estimation - predictions_snpsincontamestim)^2)
print(paste("Mean Squared Error:", mse_snpsincontamestim))

# Calculate R-squared
r_squared_snpsincontamestim <- 1 - (sum((test_data$nr_sn_ps_used_in_contamination_estimation - predictions_snpsincontamestim)^2) / 
                               sum((test_data$nr_sn_ps_used_in_contamination_estimation - mean(test_data$nr_sn_ps_used_in_contamination_estimation))^2))
print(paste("R-squared:", r_squared_snpsincontamestim))

# Plot variable importance
varImpPlot(rf_model_snpsincontamestim)

#Hyperparameter Tuning

# Define the grid of hyperparameters
tune_grid <- expand.grid(
  mtry = c(1, 2, 3)  # Number of variables randomly sampled as candidates at each split
)

# Set up cross-validation
train_control <- trainControl(
  method = "cv",     # Cross-validation
  number = 5,        # Number of folds
  search = "grid"    # Perform a grid search
)

# Train the model using cross-validation
rf_snpsincontamestim_tuned <- train(
  nr_sn_ps_used_in_contamination_estimation ~ nr_of_input_reads + source + sample,
  data = train_data,
  method = "rf",                 # Random Forest
  trControl = train_control,     # Cross-validation settings
  tuneGrid = tune_grid,          # Hyperparameter grid
  ntree = 500,
  importance = TRUE
)

# Print the best model
print(rf_snpsincontamestim_tuned)

# Best hyperparameters
best_params_snpsincontamestim <- rf_snpsincontamestim_tuned$bestTune
print(best_params_snpsincontamestim)

# Predict on the test set using the best model
predictions_tuned_snpsincontamestim <- predict(rf_snpsincontamestim_tuned, newdata = test_data)

# Calculate Mean Squared Error (MSE) on test set
mse_test_snpsincontamestim <- mean((test_data$nr_sn_ps_used_in_contamination_estimation - predictions_tuned_snpsincontamestim)^2)
print(paste("Test Set Mean Squared Error:", mse_test_snpsincontamestim))

# Calculate R-squared on test set
r_squared_snpsincontamestim_test <- 1 - (sum((test_data$nr_sn_ps_used_in_contamination_estimation - predictions_tuned_snpsincontamestim)^2) / sum((test_data$nr_sn_ps_used_in_contamination_estimation - mean(test_data$nr_sn_ps_used_in_contamination_estimation))^2))
print(paste("Test Set R-squared:", r_squared_snpsincontamestim_test))

#----------------------------------------------------------------------------------- 
## VARIABLE NAME: percent_endogenous_dna_over_30bp ##
 
# Beta Regression
summary(project_data$percent_endogenous_dna_over_30bp)

# Transform the response variable to be within (0, 1)
#project_data <- project_data %>%
#  mutate(percent_endogenous_dna_over_30bp_scaled = (percent_endogenous_dna_over_30bp - min(percent_endogenous_dna_over_30bp) + 0.00001) /
#           (max(percent_endogenous_dna_over_30bp) - min(percent_endogenous_dna_over_30bp) + 0.00002))

# Transform the response variable to be within (0, 1)
project_data <- project_data %>%
  mutate(percent_endogenous_dna_over_30bp_new = (percent_endogenous_dna_over_30bp/100))

# Applying the same transformation to 'train_data' and 'test_data' too
train_data <- train_data %>%
  mutate(percent_endogenous_dna_over_30bp_new = percent_endogenous_dna_over_30bp / 100)

test_data <- test_data %>%
  mutate(percent_endogenous_dna_over_30bp_new = percent_endogenous_dna_over_30bp / 100)


# Check the summary of the transformed variable
summary(project_data$percent_endogenous_dna_over_30bp_new)

project_data %>%
  ggplot(aes(x=nr_of_input_reads,
             y=percent_endogenous_dna_over_30bp_new,
             shape=source))+
  theme_bw()+
  xlab("No. of Input Reads") +
  ylab("Endogenous DNA over 300bp(%)") +
  geom_point(size=3)+
  scale_shape_manual(values=c(21,22))


project_data %>%
  ggplot(aes(x=nr_of_input_reads,
             y=percent_endogenous_dna_over_30bp_new,
             fill=sample,
             shape=source)) +
  theme_bw() +
  xlab("No. of Input Reads") +
  ylab("Endogenous DNA over 30bp(%)") +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,22)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              aes(col=source, group=source),
              linetype= 'solid',
              se = TRUE) +
  labs(
    title = "'Endogenous DNA over 30bp(%)' vs 'No. of Input Reads'",
    subtitle = "Scatterplot: Added colors for different samples"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"), # Title font size
    plot.subtitle = element_text(size = 14),             # Subtitle font size
    axis.title.x = element_text(size = 12),              # X-axis label font size
    axis.title.y = element_text(size = 12)               # Y-axis label font size
  ) 


# Complex Model with Interaction Term
model_percentendoover30bp.beta.int <- betareg(percent_endogenous_dna_over_30bp_new~nr_of_input_reads*source, 
                          data=project_data)
summary(model_percentendoover30bp.beta.int)

#only source is significant

# Additive Model without Interaction Term
model_percentendoover30bp.beta.add <- betareg(percent_endogenous_dna_over_30bp_new~nr_of_input_reads+source, 
                          data=project_data)
summary(model_percentendoover30bp.beta.add)

AIC(model_percentendoover30bp.beta.int,
    model_percentendoover30bp.beta.add) %>% 
  dplyr::arrange(AIC)

par(mfrow = c(3, 2))
plot(model_percentendoover30bp.beta.add, which = 1:6)
par(mfrow = c(1, 1))

# Mixed effects?
#model_percentendoover30bp.beta.mef.int <- glmmTMB(percent_endogenous_dna_over_30bp_new~nr_of_input_reads*source+(1|sample),
#                                              data=project_data, 
#                                              family=beta_family(link="logit"))
#summary(model_percentendoover30bp.beta.mef.int)


model_percentendoover30bp.beta.mef.add <- glmmTMB(percent_endogenous_dna_over_30bp_new~nr_of_input_reads+source+(1|sample),
                          data=project_data, 
                          family=beta_family(link="logit"))
summary(model_percentendoover30bp.beta.mef.add)

model_percentendoover30bp.beta.fef <- glmmTMB(percent_endogenous_dna_over_30bp_new~nr_of_input_reads+source,
                          data=project_data, 
                          family=beta_family(link="logit"))
summary(model_percentendoover30bp.beta.fef)

anova(model_percentendoover30bp.beta.mef.add,
      model_percentendoover30bp.beta.fef)
AIC(model_percentendoover30bp.beta.mef.add,
    model_percentendoover30bp.beta.fef) %>% 
  dplyr::arrange(AIC)

# Calculate R-squared values for the mixed effect models
# r_squared_values_mefint_percentendoover30bp <- r.squaredGLMM(model_percentendoover30bp.beta.mef.int)
r_squared_values_mefadd_percentendoover30bp <- r.squaredGLMM(model_percentendoover30bp.beta.mef.add)
r_squared_values_fef_percentendoover30bp <- r.squaredGLMM(model_percentendoover30bp.beta.fef)

# Print the R-squared values
# print(r_squared_values_mefint_percentendoover30bp)
print(r_squared_values_mefadd_percentendoover30bp)
print(r_squared_values_fef_percentendoover30bp)

# Simulate residuals
sim_res <- simulateResiduals(fittedModel = model_percentendoover30bp.beta.mef.add)

# Plot QQ plot of the residuals
plotQQunif(sim_res)


## RANDOM FOREST ##

# Train the Random Forest model
rf_model_percentendoover30bp <- randomForest(percent_endogenous_dna_over_30bp_new ~ nr_of_input_reads + source + sample, 
                         data = train_data, 
                         ntree = 500, 
                         importance = TRUE)

# Print the model summary
print(rf_model_percentendoover30bp)

# Make predictions on the test set
predictions_percentendoover30bp <- predict(rf_model_percentendoover30bp, newdata = test_data)

# Calculate performance metrics
mse_rf_percentendoover30bp <- mean((test_data$percent_endogenous_dna_over_30bp_new - predictions_percentendoover30bp)^2)
r_squared_rf_percentendoover30bp <- 1 - (sum((test_data$percent_endogenous_dna_over_30bp_new - predictions_percentendoover30bp)^2) / 
                    sum((test_data$percent_endogenous_dna_over_30bp_new - mean(test_data$percent_endogenous_dna_over_30bp_new))^2))

print(paste("Test Set Mean Squared Error:", mse_rf_percentendoover30bp))
print(paste("Test Set R-squared:", r_squared_rf_percentendoover30bp))

# Plot variable importance
varImpPlot(rf_model_percentendoover30bp)

## HYPERPARAMETER TUNING ##

# Define the control method for cross-validation
train_control <- trainControl(method = "cv", number = 5)

# Define the hyperparameter grid
tune_grid <- expand.grid(mtry = c(1, 2, 3, 4))

# Train the model using cross-validation
rf_tuned_percentendoover30bp <- train(
  percent_endogenous_dna_over_30bp_new ~ nr_of_input_reads + source + sample,
  data = train_data,
  method = "rf",
  trControl = train_control,
  tuneGrid = tune_grid,
  ntree = 500,
  importance = TRUE
)

# Print the best model
print(rf_tuned_percentendoover30bp$bestTune)

# Make predictions on the test set using the tuned model
predictions_tuned_percentendoover30bp <- predict(rf_tuned_percentendoover30bp, newdata = test_data)

# Calculate performance metrics for the tuned model
mse_tuned_rf_percentendoover30bp <- mean((test_data$percent_endogenous_dna_over_30bp_new - predictions_tuned_percentendoover30bp)^2)
r_squared_tuned_rf_percentendoover30bp <- 1 - (sum((test_data$percent_endogenous_dna_over_30bp_new - predictions_tuned_percentendoover30bp)^2) /
                          sum((test_data$percent_endogenous_dna_over_30bp_new - mean(test_data$percent_endogenous_dna_over_30bp_new))^2))

print(paste("Tuned Test Set Mean Squared Error:", mse_tuned_rf_percentendoover30bp))
print(paste("Tuned Test Set R-squared:", r_squared_tuned_rf_percentendoover30bp))

#-------------------------------------------------------------------------------------

## VARIABLE NAME: proportion_of_duplicate_reads ##

# Beta Regression

# Transform the response variable to be within (0, 1)
#project_data <- project_data %>%
#  mutate(proportion_of_duplicate_reads_scaled = (proportion_of_duplicate_reads - min(proportion_of_duplicate_reads) + 0.00001) /
#           (max(proportion_of_duplicate_reads) - min(proportion_of_duplicate_reads) + 0.00002))

# Check the summary of the transformed variable
summary(project_data$proportion_of_duplicate_reads)

project_data %>%
  ggplot(aes(x=nr_of_input_reads,
             y=proportion_of_duplicate_reads,
             shape=source))+
  theme_bw()+
  xlab("No. of Input Reads") +
  ylab("Proportion of Duplicate Reads") +
  geom_point(size=3)+
  scale_shape_manual(values=c(21,22))


project_data %>%
  ggplot(aes(x=nr_of_input_reads,
             y=proportion_of_duplicate_reads,
             fill=sample,
             shape=source)) +
  theme_bw() +
  xlab("No. of Input Reads") +
  ylab("Proportion of Duplicate Reads") +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,22)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              aes(col=source, group=source),
              linetype= 'solid',
              se = TRUE) +
  labs(
    title = "'Proportion of Duplicate Reads' vs 'No. of Input Reads'",
    subtitle = "Scatterplot: Added colors for different samples"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"), # Title font size
    plot.subtitle = element_text(size = 14, face= "bold"),             # Subtitle font size
    axis.title.x = element_text(size = 12, face = "bold"),              # X-axis label font size
    axis.title.y = element_text(size = 12, face = "bold")               # Y-axis label font size
  ) 


# Complex Model with Interaction Term
propdupread.beta.int <- betareg(proportion_of_duplicate_reads~nr_of_input_reads*source, 
                          data=project_data)
summary(propdupread.beta.int)

# Additive Model without Interaction Term
propdupread.beta.add <- betareg(proportion_of_duplicate_reads~nr_of_input_reads+source, 
                          data=project_data)
summary(propdupread.beta.add)

propdupread.beta.add.simpler <- betareg(proportion_of_duplicate_reads~source, 
                                data=project_data)
summary(propdupread.beta.add.simpler)


AIC(propdupread.beta.int,
    propdupread.beta.add.simpler) %>% 
  dplyr::arrange(AIC)

par(mfrow = c(3, 2))
plot(propdupread.beta.add.simpler, which = 1:6)
par(mfrow = c(1, 1))

# Mixed effects?

propdupread.beta.mef.add.simpler <- glmmTMB(proportion_of_duplicate_reads~source+(1|sample),
                          data=project_data, 
                          family=beta_family(link="logit"))
summary(propdupread.beta.mef.add.simpler)

propdupread.beta.fef <- glmmTMB(proportion_of_duplicate_reads~source,
                          data=project_data, 
                          family=beta_family(link="logit"))
summary(propdupread.beta.fef)

anova(propdupread.beta.mef.add.simpler,
      propdupread.beta.fef)
AIC(propdupread.beta.mef.add.simpler,
    propdupread.beta.fef)

# Calculate R-squared values for the mixed effect models
r_squared_values_mefadd_propdupread <- r.squaredGLMM(propdupread.beta.mef.add.simpler)
r_squared_values_fef_propdupread <- r.squaredGLMM(propdupread.beta.fef)

# Print the R-squared values
print(r_squared_values_mefadd_propdupread)
print(r_squared_values_fef_propdupread)

# Simulate residuals
sim_res_propdupread <- simulateResiduals(fittedModel = propdupread.beta.mef.add.simpler)

# Plot QQ plot of the residuals
plotQQunif(sim_res_propdupread)

## RANDOM FOREST ##

# Train the Random Forest model
rf_model_propdupread <- randomForest(proportion_of_duplicate_reads ~ nr_of_input_reads + source + sample, 
                         data = train_data, 
                         ntree = 500, 
                         importance = TRUE)

# Print the model summary
print(rf_model_propdupread)

# Make predictions on the test set
predictions_propdupread <- predict(rf_model_propdupread, newdata = test_data)

# Calculate performance metrics
mse_propdupread <- mean((test_data$proportion_of_duplicate_reads - predictions_propdupread)^2)
r_squared_propdupread <- 1 - (sum((test_data$proportion_of_duplicate_reads - predictions_propdupread)^2) / 
                    sum((test_data$proportion_of_duplicate_reads - mean(test_data$proportion_of_duplicate_reads))^2))

print(paste("Test Set Mean Squared Error:", mse_propdupread))
print(paste("Test Set R-squared:", r_squared_propdupread))

# Plot variable importance
varImpPlot(rf_model_propdupread)

## HYPERPARAMETER TUNING ##

# Define the control method for cross-validation
train_control <- trainControl(method = "cv", number = 5)

# Define the hyperparameter grid
tune_grid <- expand.grid(mtry = c(1, 2, 3, 4))

# Train the model using cross-validation
rf_tuned_propdupread <- train(
  proportion_of_duplicate_reads ~ nr_of_input_reads + source + sample,
  data = train_data,
  method = "rf",
  trControl = train_control,
  tuneGrid = tune_grid,
  ntree = 500,
  importance = TRUE
)

# Print the best model
print(rf_tuned_propdupread$bestTune)

# Make predictions on the test set using the tuned model
predictions_tuned_propdupread <- predict(rf_tuned_propdupread, newdata = test_data)

# Calculate performance metrics for the tuned model
mse_tuned_propdupread <- mean((test_data$proportion_of_duplicate_reads - predictions_tuned_propdupread)^2)
r_squared_tuned_propdupread <- 1 - (sum((test_data$proportion_of_duplicate_reads - predictions_tuned_propdupread)^2) /
                          sum((test_data$proportion_of_duplicate_reads - mean(test_data$proportion_of_duplicate_reads))^2))

print(paste("Tuned Test Set Mean Squared Error:", mse_tuned_propdupread))
print(paste("Tuned Test Set R-squared:", r_squared_tuned_propdupread))

#-----------------------------------------------------------------------------------------

## VARIABLE NAME: percent_gc_of_unique_reads ##

# Beta Regression
summary(project_data$percent_gc_of_unique_reads)

# Transform the response variable to be within (0, 1)
project_data <- project_data %>%
  mutate(percent_gc_of_unique_reads_new = (percent_gc_of_unique_reads/100))

# Applying the same transformation to 'train_data' and 'test_data' too
train_data <- train_data %>%
  mutate(percent_gc_of_unique_reads_new = percent_gc_of_unique_reads / 100)

test_data <- test_data %>%
  mutate(percent_gc_of_unique_reads_new = percent_gc_of_unique_reads / 100)


# Check the summary of the transformed variable
summary(project_data$percent_gc_of_unique_reads_new)

project_data %>%
  ggplot(aes(x=nr_of_input_reads,
             y=percent_gc_of_unique_reads_new,
             shape=source))+
  theme_bw()+
  xlab("No. of Input Reads") +
  ylab("Covered sites % which are G or C in Reference") +
  geom_point(size=3)+
  scale_shape_manual(values=c(21,22))

project_data %>%
  ggplot(aes(x=nr_of_input_reads,
             y=percent_gc_of_unique_reads_new,
             fill=sample,
             shape=source)) +
  theme_bw() +
  xlab("No. of Input Reads") +
  ylab("Covered site % which are G or C in Reference") +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,22)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              aes(col=source, group=source),
              linetype= 'solid',
              se = TRUE) +
  labs(
    title = "'% GC of Unique Reads' vs 'No. of Input Reads'",
    subtitle = "Scatterplot: Added colors for different samples"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"), # Title font size
    plot.subtitle = element_text(size = 14, face= "bold"),             # Subtitle font size
    axis.title.x = element_text(size = 12, face = "bold"),              # X-axis label font size
    axis.title.y = element_text(size = 12, face = "bold")               # Y-axis label font size
  ) 



# Complex Model with Interaction Term
coveredgc.beta.int <- betareg(percent_gc_of_unique_reads_new~nr_of_input_reads*source, 
                                data=project_data)
summary(coveredgc.beta.int)

# Additive Model without Interaction Term
coveredgc.beta.add <- betareg(percent_gc_of_unique_reads_new~nr_of_input_reads+source, 
                                data=project_data)
summary(coveredgc.beta.add)

AIC(coveredgc.beta.int,
    coveredgc.beta.add) %>% 
  dplyr::arrange(AIC)

par(mfrow = c(3, 2))
plot(coveredgc.beta.add, which = 1:6)
par(mfrow = c(1, 1))

# Mixed effects?

coveredgc.beta.mef.add <- glmmTMB(percent_gc_of_unique_reads_new~nr_of_input_reads+source+(1|sample),
                                data=project_data, 
                                family=beta_family(link="logit"))
summary(coveredgc.beta.mef.add)


coveredgc.beta.fef <- glmmTMB(percent_gc_of_unique_reads_new~nr_of_input_reads+source,
                                data=project_data, 
                                family=beta_family(link="logit"))
summary(coveredgc.beta.fef)

anova(coveredgc.beta.mef.add,
      coveredgc.beta.fef)
AIC(coveredgc.beta.mef.add,
    coveredgc.beta.fef)

# Calculate R-squared values for the mixed effect models
r_squared_values_mefadd_coveredgc <- r.squaredGLMM(coveredgc.beta.mef.add)
r_squared_values_fef_coveredgc <- r.squaredGLMM(coveredgc.beta.fef)

# Print the R-squared values
print(r_squared_values_mefadd_coveredgc)
print(r_squared_values_fef_coveredgc)

# Simulate residuals
sim_res_coveredgc <- simulateResiduals(fittedModel = coveredgc.beta.mef.add)

# Plot QQ plot of the residuals
plotQQunif(sim_res_coveredgc)

## RANDOM FOREST ##

# Train the Random Forest model
rf_model_coveredgc <- randomForest(percent_gc_of_unique_reads_new ~ nr_of_input_reads + source + sample, 
                                     data = train_data, 
                                     ntree = 500, 
                                     importance = TRUE)

# Print the model summary
print(rf_model_coveredgc)

# Make predictions on the test set
predictions_coveredgc <- predict(rf_model_coveredgc, newdata = test_data)

# Calculate performance metrics
mse_coveredgc <- mean((test_data$percent_gc_of_unique_reads_new - predictions_coveredgc)^2)
r_squared_coveredgc <- 1 - (sum((test_data$percent_gc_of_unique_reads_new - predictions_coveredgc)^2) / 
                                sum((test_data$percent_gc_of_unique_reads_new - mean(test_data$percent_gc_of_unique_reads_new))^2))

print(paste("Test Set Mean Squared Error:", mse_coveredgc))
print(paste("Test Set R-squared:", r_squared_coveredgc))

# Plot variable importance
varImpPlot(rf_model_coveredgc)

# Min-Max scaling function
min_max_scaling <- function(x, new_min, new_max) {
  (x - min(x)) / (max(x) - min(x)) * (new_max - new_min) + new_min
}

# Normalize the predictions to the range of the actual values
predictions_scaled <- min_max_scaling(predictions_coveredgc, min(test_data$percent_gc_of_unique_reads_new), max(test_data$percent_gc_of_unique_reads_new))

# Calculate performance metrics with scaled predictions
mse_coveredgc <- mean((test_data$percent_gc_of_unique_reads_new - predictions_scaled)^2)
r_squared_coveredgc <- 1 - (sum((test_data$percent_gc_of_unique_reads_new - predictions_scaled)^2) / 
                              sum((test_data$percent_gc_of_unique_reads_new - mean(test_data$percent_gc_of_unique_reads_new))^2))

print(paste("Test Set Mean Squared Error:", mse_coveredgc))
print(paste("Test Set R-squared:", r_squared_coveredgc))


## HYPERPARAMETER TUNING ##

# Define the control method for cross-validation
train_control <- trainControl(method = "cv", number = 5)

# Define the hyperparameter grid
tune_grid <- expand.grid(mtry = c(1, 2, 3, 4))

# Train the model using cross-validation
#set.seed(1871448)
rf_tuned_coveredgc <- train(
  percent_gc_of_unique_reads_new ~ nr_of_input_reads + source + sample,
  data = train_data,
  method = "rf",
  trControl = train_control,
  tuneGrid = tune_grid,
  ntree = 500,
  importance = TRUE
)

# Print the best model
print(rf_tuned_coveredgc$bestTune)

# Make predictions on the test set using the tuned model
predictions_tuned_coveredgc <- predict(rf_tuned_coveredgc, newdata = test_data)

# Calculate performance metrics for the tuned model
mse_tuned_coveredgc <- mean((test_data$percent_gc_of_unique_reads - predictions_tuned_coveredgc)^2)
r_squared_tuned_coveredgc <- 1 - (sum((test_data$percent_gc_of_unique_reads - predictions_tuned_coveredgc)^2) /
                                      sum((test_data$percent_gc_of_unique_reads - mean(test_data$percent_gc_of_unique_reads))^2))

print(paste("Tuned Test Set Mean Squared Error:", mse_tuned_coveredgc))
print(paste("Tuned Test Set R-squared:", r_squared_tuned_coveredgc))

#Min-Max scaling function
min_max_scaling_tuned <- function(x, new_min, new_max) {
 (x - min(x)) / (max(x) - min(x)) * (new_max - new_min) + new_min
}

# Normalize the predictions to the range of the actual values
predictions_tuned_coveredgc_scaled <- min_max_scaling_tuned(predictions_tuned_coveredgc, min(test_data$percent_gc_of_unique_reads_new), max(test_data$percent_gc_of_unique_reads_new))

# Calculate performance metrics with scaled predictions
mse_coveredgc_tuned <- mean((test_data$percent_gc_of_unique_reads_new - predictions_tuned_coveredgc_scaled)^2)
r_squared_coveredgc_tuned <- 1 - (sum((test_data$percent_gc_of_unique_reads_new - predictions_tuned_coveredgc_scaled)^2) / 
                              sum((test_data$percent_gc_of_unique_reads_new - mean(test_data$percent_gc_of_unique_reads_new))^2))

print(paste("Test Set Mean Squared Error:", mse_coveredgc_tuned))
print(paste("Test Set R-squared:", r_squared_coveredgc_tuned))

#---------------------------------------------------------------------------------------

## VARIABLE NAME: nuclear_contamination_m1_ml ##

# Beta Regression
summary(project_data$nuclear_contamination_m1_ml)

# Transform the response variable to be within (0, 1)
project_data <- project_data %>%
  mutate(nuclear_contamination_m1_ml_scaled = (nuclear_contamination_m1_ml - min(nuclear_contamination_m1_ml) + 0.00001) /
           (max(nuclear_contamination_m1_ml) - min(nuclear_contamination_m1_ml) + 0.00002))

# Check the summary of the transformed variable
summary(project_data$nuclear_contamination_m1_ml_scaled)

project_data %>%
  ggplot(aes(x=nr_of_input_reads,
             y=nuclear_contamination_m1_ml_scaled,
             shape=source))+
  theme_bw()+
  xlab("No. of Input Reads") +
  ylab("Nuclear Contamination M1") +
  geom_point(size=3)+
  scale_shape_manual(values=c(21,22))


project_data %>%
  ggplot(aes(x=nr_of_input_reads,
             y=nuclear_contamination_m1_ml_scaled,
             fill=sample,
             shape=source)) +
  theme_bw() +
  xlab("No. of Input Reads") +
  ylab("Nuclear Contamination M1") +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,22)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              aes(col=source, group=source),
              linetype= 'solid',
              se = TRUE) +
  labs(
    title = "'Nuclear Contamination M1' vs 'No. of Input Reads'",
    subtitle = "Scatterplot: Added colors for different samples"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"), # Title font size
    plot.subtitle = element_text(size = 14, face= "bold"),             # Subtitle font size
    axis.title.x = element_text(size = 12, face = "bold"),              # X-axis label font size
    axis.title.y = element_text(size = 12, face = "bold")               # Y-axis label font size
  ) 


# Complex Model with Interaction Term
nuclcont.beta.int <- betareg(nuclear_contamination_m1_ml_scaled~nr_of_input_reads*source, 
                              data=project_data)
summary(nuclcont.beta.int)

# Additive Model without Interaction Term
nuclcont.beta.add <- betareg(nuclear_contamination_m1_ml_scaled~nr_of_input_reads+source, 
                              data=project_data)
summary(nuclcont.beta.add)

nuclcont.beta.add.simpler1 <- betareg(nuclear_contamination_m1_ml_scaled~nr_of_input_reads, 
                             data=project_data)
summary(nuclcont.beta.add.simpler1)

nuclcont.beta.add.simpler2 <- betareg(nuclear_contamination_m1_ml_scaled~source, 
                                      data=project_data)
summary(nuclcont.beta.add.simpler2)

nuclcont.beta.add.simplest <- betareg(nuclear_contamination_m1_ml_scaled~ 1-1, 
                                      data=project_data)
summary(nuclcont.beta.add.simplest)


AIC(nuclcont.beta.int,
    nuclcont.beta.add.simplest) %>% 
  dplyr::arrange(AIC)

par(mfrow = c(3, 2))
plot(nuclcont.beta.add.simplest, which = 1:6)
par(mfrow = c(1, 1))

# Mixed effects?
nuclcont.beta.mef.add.simplest <- glmmTMB(nuclear_contamination_m1_ml_scaled~1-1+(1|sample),
                              data=project_data, 
                              family=beta_family(link="logit"))
summary(nuclcont.beta.mef.add.simplest)


nuclcont.beta.fef <- glmmTMB(nuclear_contamination_m1_ml_scaled~1-1,
                              data=project_data, 
                              family=beta_family(link="logit"))
summary(nuclcont.beta.fef)

anova(nuclcont.beta.mef.add.simplest,
      nuclcont.beta.fef)

# Calculate R-squared values for the mixed effect models
#r_squared_values_mefint_nuclcont <- r.squaredGLMM(nuclcont.beta.mef.int)
r_squared_values_mefadd_nuclcont <- r.squaredGLMM(nuclcont.beta.mef.add.simplest)
r_squared_values_fef_nuclcont <- r.squaredGLMM(nuclcont.beta.fef)

# Print the R-squared values
#print(r_squared_values_mefint_nuclcont)
print(r_squared_values_mefadd_nuclcont)
print(r_squared_values_fef_nuclcont)

# Simulate residuals
sim_res_nuclcont <- simulateResiduals(fittedModel = nuclcont.beta.mef.add.simplest)

# Plot QQ plot of the residuals
plotQQunif(sim_res_nuclcont)

## RANDOM FOREST ##

# Train the Random Forest model
rf_model_nuclcont <- randomForest(nuclear_contamination_m1_ml ~ nr_of_input_reads + source + sample, 
                                   data = train_data, 
                                   ntree = 500, 
                                   importance = TRUE)

# Print the model summary
print(rf_model_nuclcont)

# Make predictions on the test set
predictions_nuclcont <- predict(rf_model_nuclcont, newdata = test_data)

# Calculate performance metrics
mse_nuclcont <- mean((test_data$nuclear_contamination_m1_ml - predictions_nuclcont)^2)
r_squared_nuclcont <- 1 - (sum((test_data$nuclear_contamination_m1_ml - predictions_nuclcont)^2) / 
                              sum((test_data$nuclear_contamination_m1_ml - mean(test_data$nuclear_contamination_m1_ml))^2))

print(paste("Test Set Mean Squared Error:", mse_nuclcont))
print(paste("Test Set R-squared:", r_squared_nuclcont))

# Plot variable importance
varImpPlot(rf_model_nuclcont)

## HYPERPARAMETER TUNING ##

# Define the control method for cross-validation
train_control <- trainControl(method = "cv", number = 5)

# Define the hyperparameter grid
tune_grid <- expand.grid(mtry = c(1, 2, 3, 4))

# Train the model using cross-validation
#set.seed(1871448)
rf_tuned_nuclcont <- train(
  nuclear_contamination_m1_ml ~ nr_of_input_reads + source + sample,
  data = train_data,
  method = "rf",
  trControl = train_control,
  tuneGrid = tune_grid,
  ntree = 500,
  importance = TRUE
)

# Print the best model
print(rf_tuned_nuclcont$bestTune)

# Make predictions on the test set using the tuned model
predictions_tuned_nuclcont <- predict(rf_tuned_nuclcont, newdata = test_data)

# Calculate performance metrics for the tuned model
mse_tuned_nuclcont <- mean((test_data$nuclear_contamination_m1_ml - predictions_tuned_nuclcont)^2)
r_squared_tuned_nuclcont <- 1 - (sum((test_data$nuclear_contamination_m1_ml - predictions_tuned_nuclcont)^2) /
                                    sum((test_data$nuclear_contamination_m1_ml - mean(test_data$nuclear_contamination_m1_ml))^2))

print(paste("Tuned Test Set Mean Squared Error:", mse_tuned_nuclcont))
print(paste("Tuned Test Set R-squared:", r_squared_tuned_nuclcont))
#-----------------------------------------------------------------------------------------
