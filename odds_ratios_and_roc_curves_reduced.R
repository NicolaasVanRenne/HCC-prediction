
###########################################################
#
# Calculate odds ratios and ROC curves
# For the training and validation dataset
# For the prognosis based on the reduced gene signature
#
##########################################################

### Packages ###

library(dplyr)
library(gtsummary)
library(survival)
library(survminer)
library(timeROC)
library(gt)
library(rms)
library(flextable)
library(officer)
library(ggplot2)
library(pagedown)


### Settings ###

savePlots <- T


### Load data ###

# Load the data with the prediction results
data_pred_train <- read.delim("output/3_NTP/training_NTP_prediction_result.txt")
data_pred_val <- read.delim("output/3_NTP/validation_NTP_prediction_result.txt")
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

# Combine the dataset with the prediction results with the metadata
data_train <- full_join(meta_train, data_pred_train, by = c("SubjectName" = "sample.names"))
data_val <- full_join(meta_val, data_pred_val, by = c("SubjectName" = "sample.names"))


### Calculate odds ratios ###

# Report odds ratios for the HCC occurrence of the good, intermediate and poor prognosis groups

## Training data

# Fit model
model_train <- glm(HCC_event ~ prognosis, data = data_train, family = binomial)

# Results
summary(model_train)

# Results based on Wald test
exp(coef(model_train))  # odds ratios
exp(confint.default(model_train)) # confidence intervals
coef(summary(model_train))[,4] # p-values

# Results based on profile likelihood
exp(coef(model_train))  # odds ratios
exp(confint(model_train)) # confidence intervals
anova(model_train, test = "LRT") # p-value for prognosis overall
# p-values per prognosis category based on likelihoods
# - for poor vs good
data_sub <- subset(data_train, prognosis %in% c("good", "poor"))
data_sub$prognosis <- droplevels(data_sub$prognosis)
model_train_sub <- glm(HCC_event ~ prognosis, family = binomial, data = data_sub)
anova(model_train_sub, test = "LRT")
p_train_poor <- anova(model_train_sub, test = "LRT")["prognosis","Pr(>Chi)"]
exp(confint(model_train_sub))
# - for intermediate vs good
data_sub <- subset(data_train, prognosis %in% c("good", "intermediate"))
data_sub$prognosis <- droplevels(data_sub$prognosis)
model_train_sub <- glm(HCC_event ~ prognosis, family = binomial, data = data_sub)
anova(model_train_sub, test = "LRT")
exp(confint(model_train_sub))
p_train_inter <- anova(model_train_sub, test = "LRT")["prognosis","Pr(>Chi)"]

# Outcome table
results_tbl <- model_train %>%
  tbl_regression(
    exponentiate = TRUE,
    pvalue_fun = function(x) {
      vapply(x, function(p) {
        if (!is.na(p) & p < 0.001) {
          exp <- floor(log10(abs(p)))
          base <- round(p / 10^exp, 3)
          paste0(base, "*10^", exp)
        } else {
          as.character(round(p, 3))
        }
      }, character(1))
    }
  )
results_tbl[["table_body"]][["p.value"]] <- c(NA, NA, p_train_inter, p_train_poor)
results_tbl <- results_tbl %>%
  bold_labels() %>%
  modify_caption("Training model (reduced signature): regression results")
results_tbl <- results_tbl %>%
  modify_abbreviation("")
results_tbl

if (savePlots == T) {
  results_tbl_docx <- results_tbl %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "check/model_train_reduc_oddsratio.docx")
  
  ft <- results_tbl %>%
    as_flex_table()
  ft <- flextable::delete_part(ft, part = "footer")
  ft <- flextable::set_caption(
    ft,
    caption = flextable::as_paragraph(
      flextable::as_b("Training model (reduced signature): regression results")
    )
  )
  ft %>%
    flextable::save_as_html(path = "check/model_train_reduc_oddsratio.html")
  pagedown::chrome_print(
    "check/model_train_reduc_oddsratio.html",
    output = "check/model_train_reduc_oddsratio.pdf"
  )
}


## Validation data

# Fit model
model_val <- glm(HCC_event ~ prognosis, data = data_val, family = binomial)

# Results
summary(model_val)

# Results based on Wald test
exp(coef(model_val))  # odds ratios
exp(confint.default(model_val)) # confidence intervals
coef(summary(model_val))[,4] # p-values

# Results based on profile likelihood
exp(coef(model_val))  # odds ratios
exp(confint(model_val)) # confidence intervals
anova(model_val, test = "LRT") # p-value for prognosis overall
# p-values per prognosis category based on likelihoods
# - for poor vs good
data_sub <- subset(data_val, prognosis %in% c("good", "poor"))
data_sub$prognosis <- droplevels(data_sub$prognosis)
model_val_sub <- glm(HCC_event ~ prognosis, family = binomial, data = data_sub)
anova(model_val_sub, test = "LRT")
exp(confint(model_val_sub))
p_val_poor <- anova(model_val_sub, test = "LRT")["prognosis","Pr(>Chi)"]
# - for intermediate vs good
data_sub <- subset(data_val, prognosis %in% c("good", "intermediate"))
data_sub$prognosis <- droplevels(data_sub$prognosis)
model_val_sub <- glm(HCC_event ~ prognosis, family = binomial, data = data_sub)
anova(model_val_sub, test = "LRT")
exp(confint(model_val_sub))
p_val_inter <- anova(model_val_sub, test = "LRT")["prognosis","Pr(>Chi)"]

# Outcome table
results_tbl <- model_val %>%
  tbl_regression(
    exponentiate = TRUE,
    pvalue_fun = function(x) {
      vapply(x, function(p) {
        if (!is.na(p) & p < 0.001) {
          exp <- floor(log10(abs(p)))
          base <- round(p / 10^exp, 3)
          paste0(base, "*10^", exp)
        } else {
          as.character(round(p, 3))
        }
      }, character(1))
    }
  )
results_tbl[["table_body"]][["p.value"]] <- c(NA, NA, p_val_inter, p_val_poor)
results_tbl %>%
  bold_labels() %>%
  modify_caption("Validation model (reduced signature): regression results")

if (savePlots == T) {
  results_tbl_docx <- results_tbl %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = "check/model_val_reduc_oddsratio.docx")
  
  ft <- results_tbl %>%
    as_flex_table()
  ft <- flextable::delete_part(ft, part = "footer")
  ft <- flextable::set_caption(
    ft,
    caption = flextable::as_paragraph(
      flextable::as_b("Validation model (reduced signature): regression results")
    )
  )
  ft %>%
    flextable::save_as_html(path = "check/model_val_reduc_oddsratio.html")
  pagedown::chrome_print(
    "check/model_val_reduc_oddsratio.html",
    output = "check/model_val_reduc_oddsratio.pdf"
  )
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
    title = paste0("Training: ROC Curve, \n using reduced Van Renne signature \n — AUC = ", round(auc_train, 3)),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p

if (savePlots == T) {
  ggsave("check/roc_train_reduc.png", p, width = 8, height = 6, dpi = 600, bg = "white")
  ggsave("check/roc_train_reduc.pdf", p, width = 8, height = 6, dpi = 600, bg = "white")
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
    title = paste0("Validation: ROC Curve, \n using reduced Van Renne signature \n — AUC = ", round(auc_val, 3)),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p

if (savePlots == T) {
  ggsave("check/roc_val_reduc.png", p, width = 8, height = 6, dpi = 600, bg = "white")
  ggsave("check/roc_val_reduc.pdf", p, width = 8, height = 6, dpi = 600, bg = "white")
}