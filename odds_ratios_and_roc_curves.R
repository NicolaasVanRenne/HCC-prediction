###########################################################
#
# Calculate odds ratios and ROC curves
# For the training and validation dataset
# For the prognosis based on the gene signature
#
##########################################################

### Packages ###

library(dplyr)
library(gtsummary)
library(pROC)
library(gt)
library(rms)
library(flextable)
library(officer)


### Settings ###

# To save the plots
savePlots <- F


### Load data ###

# Load the data with the prediction results
data_pred_train <- read.delim("output/training_NTP_prediction_result.txt")
data_pred_val <- read.delim("output/validation_NTP_prediction_result.txt")
# Load the clinical metadata
meta_train <- read.delim("data_files/meta_data/training_full_meta.txt")
meta_val <- read.delim("data_files/meta_data/validation_full_meta.txt")


### Data preparation ###

# Adjust the prediction labels to good - intermediate - poor
get_prognosis <- function (pred.summary) {
  prognosis <- rep(0, nrow(pred.summary)) 
  for(i in seq_along(prognosis)){
    if(pred.summary[i,2] == 1 & pred.summary[i,6] < 0.05) { 
      prognosis[i] <- "poor"
    }
    if(pred.summary[i,2] == 2 & pred.summary[i,6] < 0.05) { 
      prognosis[i] <- "good"
    }
    if(pred.summary[i,6] >= 0.05){ 
      prognosis[i] <- "intermediate"
    }
  }
  return(prognosis)
}
data_pred_train$prognosis <- get_prognosis(data_pred_train)
data_pred_val$prognosis <- get_prognosis(data_pred_val)
data_pred_train$prognosis <- as.factor(data_pred_train$prognosis)
data_pred_val$prognosis <- as.factor(data_pred_val$prognosis)

# Combine the dataset of the prediction results with the metadata
data_train <- full_join(meta_train, data_pred_train, by = c("SubjectName" = "sample.names"))
data_val <- full_join(meta_val, data_pred_val, by = c("SubjectName" = "sample.names"))


### Calculate odds ratios ###

# Report odds ratios for the HCC occurrence of the good, intermediate and poor prognosis groups

## Training data

model_train <- glm(HCC_event ~ prognosis, data = data_train, family = binomial)
summary(model_train)
exp(coef(model_train))  # odds ratios

# Outcome table
results_tbl <- model_train %>%
  tbl_regression(
    exponentiate = TRUE,
  ) %>%
  bold_labels() %>%
  modify_caption("**Training model: regression results**")

results_tbl

if (savePlots == T) {
  results_tbl_docx <- results_tbl %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "model_train_oddsratio.docx")
}

## Validation data

model_val <- glm(HCC_event ~ prognosis, data = data_val, family = binomial)
summary(model_val)
exp(coef(model_val))  # odds ratios

# Outcome table
results_tbl <- model_val %>%
  tbl_regression(
    exponentiate = TRUE,
  ) %>%
  bold_labels() %>%
  modify_caption("**Validation model: regression results**")

results_tbl

if (savePlots == T) {
  results_tbl_docx <- results_tbl %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "model_val_oddsratio.docx")
}


### Create ROC curves ###

# Create ROC curves for the predicted prognosis of the gene signature compared with the observed outcome
# And calculate the corresponding AUCs

## Training

# Predicted probabilities
pred_train <- predict(model_train, type = "response")

# Compute ROC curve
roc_train <- roc(response = data_train$HCC_event, predictor = pred_train)

# AUC
auc_train <- auc(roc_train)
print(auc_train)

# Convert ROC object to data frame
roc_df_train <- data.frame(
  fpr = 1 - roc_train$specificities,
  tpr = roc_train$sensitivities
)

# Plot ROC curve
p <- ggplot(roc_df_train, aes(x = fpr, y = tpr)) +
  geom_line(color = "#1B9E77", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("ROC Curve — AUC = ", round(auc_train, 3)),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p

if (savePlots == T) {
  ggsave("roc_train.png", p, width = 8, height = 6, dpi = 600, bg = "white")
}

## Validation

# Predicted probabilities
pred_val <- predict(model_val, type = "response")

# Compute ROC curve
roc_val <- roc(response = data_val$HCC_event, predictor = pred_val)

# AUC
auc_val <- auc(roc_val)
print(auc_val)

# Convert ROC object to data frame
roc_df_val <- data.frame(
  fpr = 1 - roc_val$specificities,
  tpr = roc_val$sensitivities
)

# Plot ROC curve
p <- ggplot(roc_df_val, aes(x = fpr, y = tpr)) +
  geom_line(color = "#1B9E77", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("ROC Curve — AUC = ", round(auc_val, 3)),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p

if (savePlots == T) {
  ggsave("roc_val.png", p, width = 8, height = 6, dpi = 600, bg = "white")
}

