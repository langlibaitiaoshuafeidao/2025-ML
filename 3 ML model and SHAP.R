# 1. Set working directory
setwd("C:/Users/86139/Desktop/机器学习/20250915 6模型初筛选择3models 加KEGG DIABLO建模/HAPC VS NHAPC")

# 2. Load required packages
library(pbapply)
library(rlang)
library(tidyverse)
library(reshape2)
library(openxlsx)
library(DALEX)
library(readr)
library(dplyr)
library(caret)
library(ggplot2)
library(pROC)
library(rms)
library(rmda)
library(dcurves)
library(Hmisc)
library(ResourceSelection)
library(survey)
library(foreign)
library(plotROC)
library(shapper)
library(iml)
library(e1071)
library(ROCR)
library(corrplot)
library(lattice)
library(DMwR)
library(kernelshap)
library(shapviz)
library(glmnet)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)

# 3. Data loading and preprocessing
data <- read.csv("DIABLOtop6.csv", header = TRUE, stringsAsFactors = FALSE)
data$Result <- factor(data$Result, levels = c(0, 1), labels = c('No', 'Yes'))

preprocessParams <- preProcess(data[, -which(colnames(data)=="Result")], method = c("center", "scale"))
data_scaled_features <- predict(preprocessParams, data[, -which(colnames(data)=="Result")])
data_scaled <- cbind(Result = data$Result, data_scaled_features)

# 4. Data splitting and class balancing (SMOTE)
set.seed(52)
inTrain <- createDataPartition(y = data_scaled$Result, p = 0.7, list = FALSE)
traindata <- data_scaled[inTrain, ]
testdata <- data_scaled[-inTrain, ]

set.seed(52)
traindata_smote <- SMOTE(Result ~ ., data = traindata, k = 5, perc.over = 100, perc.under = 150)

write.csv(traindata_smote, "dev_smote.csv", row.names = FALSE)
write.csv(testdata, "vad.csv", row.names = FALSE)

# 5. Model training
models_config <- list(
  LogisticRegression = list(
    method = "glmnet",
    params = list(
      family = "binomial",
      tuneGrid = expand.grid(alpha = 0, lambda = exp(seq(-0.5, 5, length = 40)))
    )
  ),
  Lasso = list(
    method = "glmnet",
    params = list(tuneGrid = expand.grid(alpha = 1, lambda = exp(seq(-3, 3, length = 50))))
  ),
  ElasticNet = list(
    method = "glmnet",
    params = list(
      tuneGrid = expand.grid(alpha = seq(0.05,0.5, by = 0.05), 
                             lambda = exp(seq(-1, 3, length = 50))),
      maxit = 30000
    )
  )
)

set.seed(520)
train_control <- trainControl(
  method = 'repeatedcv',
  number = 5,
  repeats = 6,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  sampling = "up"
)

ML_models <- list()
train_probe <- data.frame(Result = traindata_smote$Result)
test_probe <- data.frame(Result = testdata$Result)
importance <- list()

for (model_name in names(models_config)) {
  config <- models_config[[model_name]]
  
  set.seed(52)
  fit <- do.call(train, args = c(
    list(
      x = traindata_smote[, -which(colnames(traindata_smote)=="Result")],
      y = traindata_smote$Result,
      method = config$method,
      metric = "ROC",
      trControl = train_control
    ),
    config$params
  ))
  
  ML_models[[model_name]] <- fit
  train_probe[[model_name]] <- predict(fit, newdata = traindata_smote, type = 'prob')$Yes
  test_probe[[model_name]] <- predict(fit, newdata = testdata, type = 'prob')$Yes
  importance[[model_name]] <- varImp(fit, scale = TRUE)
}

# 6. Model evaluation
evaluate_model <- function(data_type, pred_data) {
  ROC_list <- list()
  ROC_label <- list()
  Evaluation_metrics <- data.frame(
    Model = NA, Threshold = NA, Accuracy = NA, 
    Sensitivity = NA, Specificity = NA, Precision = NA, F1 = NA, AUC = NA
  )
  
  for (model_name in names(ML_models)) {
    ROC <- roc(response = pred_data$Result, predictor = pred_data[[model_name]])
    AUC <- round(auc(ROC), 3)
    
    if (AUC == 1) {
      bestp <- 0.5
    } else {
      best_coords <- coords(ROC, "best", ret = "threshold", transpose = TRUE)
      bestp <- best_coords[1]
    }
    
    if (AUC < 1) {
      CI <- ci.auc(ROC)
      ROC_label[[model_name]] <- paste0(model_name, " (AUC=", sprintf("%.3f", AUC), 
                                        ", 95%CI:", sprintf("%.3f", CI[1]), "-", sprintf("%.3f", CI[3]), ")")
    } else {
      ROC_label[[model_name]] <- paste0(model_name, " (AUC=1.000 - Potential overfitting)")
    }
    
    predlab <- factor(ifelse(pred_data[[model_name]] > bestp, "Yes", "No"), levels = c("No", "Yes"))
    cm <- confusionMatrix(data = predlab, reference = pred_data$Result, positive = "Yes", mode = "everything")
    
    conf_matrix <- table(Actual = pred_data$Result, Predicted = predlab)
    conf_df <- as.data.frame(conf_matrix) %>%
      rename(Actual = Actual, Predicted = Predicted, Count = Freq) %>%
      group_by(Actual) %>%
      mutate(Percentage = Count / sum(Count) * 100, Label = paste0(Count, "\n(", round(Percentage, 1), "%)")) %>%
      ungroup()
    
    confusion_plot <- ggplot(conf_df, aes(x = Predicted, y = Actual, fill = Count)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = Label), color = "black", size = 5) +
      scale_fill_gradient(low = "#D6EAF8", high = "#2E86C1") +
      labs(x = "Predicted", y = "Actual", title = paste0(model_name, " Confusion Matrix - ", data_type)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none"
      )
    
    ggsave(paste0(data_type, "_", model_name, "_cm_plot.pdf"), confusion_plot, width = 6, height = 5)
    
    Evaluation_metrics <- rbind(Evaluation_metrics,
                                data.frame(
                                  Model = model_name,
                                  Threshold = round(bestp, 3),
                                  Accuracy = round(cm$overall["Accuracy"], 3),
                                  Sensitivity = round(cm$byClass["Sensitivity"], 3),
                                  Specificity = round(cm$byClass["Specificity"], 3),
                                  Precision = round(cm$byClass["Precision"], 3),
                                  F1 = round(cm$byClass["F1"], 3),
                                  AUC = AUC
                                ))
    
    ROC_list[[model_name]] <- ROC
  }
  
  Evaluation_metrics <- Evaluation_metrics[-1, ]
  write.csv(Evaluation_metrics, paste0(data_type, "_Evaluation_metrics.csv"), row.names = FALSE)
  
  ROC_plot <- pROC::ggroc(ROC_list, size = 1.5, legacy.axes = TRUE) +
    theme_bw() +
    labs(title = paste0(data_type, " Set ROC Curves"), x = "1-Specificity", y = "Sensitivity") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 12, face = "bold"),
      legend.position = c(0.7, 0.25),
      legend.background = element_blank(),
      axis.title = element_text(size = 12, face = "bold"),
      panel.border = element_rect(color = "black", size = 1),
      panel.background = element_blank()
    ) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour = 'grey', linetype = 'dotdash') +
    scale_colour_discrete(breaks = names(ML_models), labels = ROC_label)
  
  ggsave(paste0(data_type, "_ROC.pdf"), ROC_plot, width = 7, height = 7)
  
  return(Evaluation_metrics)
}

train_metrics <- evaluate_model("Train", train_probe)
test_metrics <- evaluate_model("Test", test_probe)

# 7. Feature importance visualization
for (model_name in names(ML_models)) {
  imp <- importance[[model_name]]
  imp_table <- as.data.frame(imp$importance) %>%
    mutate(Features = rownames(.), .before = 1)
  
  score_col <- if ("Yes" %in% colnames(imp_table)) "Yes" else "Overall"
  if (!score_col %in% colnames(imp_table)) {
    stop(paste0(model_name, " model has no valid importance column"))
  }
  
  imp_plot <- ggplot(imp_table, aes(x = .data[[score_col]], y = reorder(Features, .data[[score_col]]))) +
    geom_bar(aes(fill = .data[[score_col]]), stat = "identity", width = 0.6) +
    scale_fill_gradient(low = "#E8F4FD", high = "#2E86AB") +
    theme_classic() +
    labs(x = "Importance Score", y = "Feature", title = paste0(model_name, " Feature Importance")) +
    theme(
      plot.title = element_text(hjust= 0.5, size = 16, face = "bold"),
      legend.position = "none",
      axis.text = element_text(size = 10, face = "bold", color = "black"),
      axis.title.x = element_text(size = 12, face = "bold", color = "black"),
      axis.title.y = element_text(size = 12, face = "bold", color = "black")
    )
  
  ggsave(paste0(model_name, "_important.pdf"), imp_plot, width = 7, height = 5, family = "serif")
}

# 8. Performance comparison
performance_comparison <- data.frame(
  Model = rep(c("LogisticRegression", "Lasso", "ElasticNet"), each = 6),
  Metric = rep(c("Accuracy", "AUC", "F1", "Precision", "Sensitivity", "Specificity"), 3),
  Value = NA
)

test_metrics <- read.csv("Test_Evaluation_metrics.csv")

for (i in 1:nrow(performance_comparison)) {
  model <- performance_comparison$Model[i]
  metric <- performance_comparison$Metric[i]
  
  if (metric == "AUC") {
    roc_obj <- roc(response = test_probe$Result, predictor = test_probe[[model]])
    performance_comparison$Value[i] <- auc(roc_obj)
  } else {
    row_idx <- which(test_metrics$Model == model)
    if (metric == "Accuracy") {
      performance_comparison$Value[i] <- as.numeric(test_metrics$Accuracy[row_idx])
    } else if (metric == "F1") {
      performance_comparison$Value[i] <- as.numeric(test_metrics$F1[row_idx])
    } else if (metric == "Precision") {
      performance_comparison$Value[i] <- as.numeric(test_metrics$Precision[row_idx])
    } else if (metric == "Sensitivity") {
      performance_comparison$Value[i] <- as.numeric(test_metrics$Sensitivity[row_idx])
    } else if (metric == "Specificity") {
      performance_comparison$Value[i] <- as.numeric(test_metrics$Specificity[row_idx])
    }
  }
}

model_labels <- c(
  "LogisticRegression" = "Logistic",
  "Lasso" = "Lasso", 
  "ElasticNet" = "ElasticNet"
)

performance_comparison$Model_Display <- factor(performance_comparison$Model, 
                                               levels = names(model_labels),
                                               labels = model_labels)

metric_levels <- c("Accuracy", "AUC", "F1", "Precision", "Sensitivity", "Specificity")
performance_comparison$Metric <- factor(performance_comparison$Metric, levels = metric_levels)

p_perf <- ggplot(performance_comparison, aes(x = Metric, y = Value, fill = Model_Display)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = c("Logistic" = "#1f77b4", 
                               "Lasso" = "#ff7f0e", 
                               "ElasticNet" = "#2ca02c")) +
  labs(title = "Model Performance Comparison (Test Set)",
       x = "", 
       y = "Score",
       fill = "Model") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "top") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  geom_hline(yintercept = seq(0, 1, 0.25), linetype = "dotted", alpha = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Value)), 
            position = position_dodge(0.8), 
            vjust = -0.5, 
            size = 4)

pdf("Model_Performance_Comparison.pdf", 12, 8)
print(p_perf)
dev.off()

# 9. Overfitting assessment
overfitting_metrics <- data.frame(
  Model = names(ML_models),
  Train_AUC = NA,
  Test_AUC = NA,
  Overfitting = NA
)

for (model_name in names(ML_models)) {
  roc_train <- roc(response = train_probe$Result, predictor = train_probe[[model_name]])
  roc_test <- roc(response = test_probe$Result, predictor = test_probe[[model_name]])
  
  overfitting_metrics$Train_AUC[overfitting_metrics$Model == model_name] <- round(auc(roc_train), 3)
  overfitting_metrics$Test_AUC[overfitting_metrics$Model == model_name] <- round(auc(roc_test), 3)
  overfitting_metrics$Overfitting[overfitting_metrics$Model == model_name] <- round(auc(roc_train) - auc(roc_test), 8)
}

write.csv(overfitting_metrics, "overfitting_metrics.csv", row.names = FALSE)

# 10. Permutation test
set.seed(123)
n_perm <- 100
models_names <- names(ML_models)
perm_results <- matrix(NA, nrow = n_perm, ncol = length(models_names))
colnames(perm_results) <- models_names

best_lasso_params <- ML_models[["Lasso"]]$bestTune
best_elastic_params <- ML_models[["ElasticNet"]]$bestTune

for (i in 1:n_perm) {
  permuted_data <- traindata_smote
  permuted_data$Result <- sample(permuted_data$Result)
  
  for (model_name in models_names) {
    if (model_name == "LogisticRegression") {
      fit <- train(Result ~ ., data = permuted_data, method = "glm", 
                   family = "binomial", trControl = trainControl(method = "none"))
    } else if (model_name == "Lasso") {
      fit <- train(Result ~ ., data = permuted_data, method = "glmnet", 
                   tuneGrid = best_lasso_params,
                   trControl = trainControl(method = "none"))
    } else {
      fit <- train(Result ~ ., data = permuted_data, method = "glmnet", 
                   tuneGrid = best_elastic_params,
                   trControl = trainControl(method = "none"))
    }
    
    pred <- predict(fit, newdata = testdata, type = "prob")$Yes
    perm_roc <- roc(response = testdata$Result, predictor = pred)
    perm_results[i, model_name] <- auc(perm_roc)
  }
}

actual_aucs <- sapply(models_names, function(m) {
  auc(roc(response = test_probe$Result, predictor = test_probe[[m]]))
})

perm_pvalues <- sapply(models_names, function(m) {
  mean(perm_results[, m] >= actual_aucs[m])
})

perm_metrics <- data.frame(Model = models_names, 
                           Actual_AUC = actual_aucs, 
                           Permutation_pvalue = perm_pvalues)
write.csv(perm_metrics, "permutation_test_results.csv", row.names = FALSE)

# 11. Final report
final_report <- data.frame(
  Model = models_names,
  Train_AUC = overfitting_metrics$Train_AUC,
  Test_AUC = overfitting_metrics$Test_AUC,
  Overfitting = overfitting_metrics$Overfitting,
  Accuracy = as.numeric(test_metrics$Accuracy),
  Sensitivity = as.numeric(test_metrics$Sensitivity),
  Specificity = as.numeric(test_metrics$Specificity),
  Precision = as.numeric(test_metrics$Precision),
  F1 = as.numeric(test_metrics$F1),
  Permutation_AUC = perm_metrics$Actual_AUC,
  Permutation_pvalue = perm_metrics$Permutation_pvalue
)

write.csv(final_report, "Final_Model_Performance_Report.csv", row.names = FALSE)

# 12. Manual best model selection
cat("Trained models list:\n")
for (model_name in names(ML_models)) {
  if (is.null(ML_models[[model_name]])) {
    cat("- ", model_name, " [Warning: Model training failed]\n")
  } else {
    cat("- ", model_name, " [Status OK]\n")
  }
}

best_model_name <- readline(prompt = "\nEnter the name of the best model (must match list above): ")

while (TRUE) {
  if (!(best_model_name %in% names(ML_models))) {
    cat("Invalid input! Model name not in available list.\n")
  } else if (is.null(ML_models[[best_model_name]])) {
    cat(paste0("Invalid input! Model", best_model_name, "training failed.\n"))
  } else {
    break
  }
  
  cat("Please select from the following models:\n")
  for (model_name in names(ML_models)) {
    if (!is.null(ML_models[[model_name]])) {
      cat("- ", model_name, "\n")
    }
  }
  best_model_name <- readline(prompt = "Enter model name: ")
}

best_model <- ML_models[[best_model_name]]
if (is.null(best_model)) {
  stop(paste0("Critical error: Model", best_model_name, "exists but is NULL. Please rerun model training."))
}

saveRDS(best_model, file = "best_model.rds")

model_selection_reason <- data.frame(
  Selection_Criteria = "Manual_Selection",
  Selected_Model = best_model_name,
  Model_Class = paste(class(best_model), collapse = ","),
  Selection_Date = Sys.Date(),
  Reason = paste0("User manually selected ", best_model_name, " as best model")
)
write.csv(model_selection_reason, "Model_Selection_Reason.csv", row.names = FALSE)

# 13. Threshold optimization
if (best_model_name %in% c("LogisticRegression", "Lasso", "ElasticNet")) {
  test_probs <- test_probe[[best_model_name]]
  
  thresholds <- seq(0.3, 0.7, by = 0.05)
  performance_at_thresholds <- data.frame()
  
  for (thresh in thresholds) {
    pred_labels <- ifelse(test_probs > thresh, "Yes", "No")
    pred_labels <- factor(pred_labels, levels = c("No", "Yes"))
    
    cm <- confusionMatrix(pred_labels, testdata$Result, positive = "Yes")
    
    performance_at_thresholds <- rbind(performance_at_thresholds, 
                                       data.frame(Threshold = thresh,
                                                  Sensitivity = cm$byClass["Sensitivity"],
                                                  Specificity = cm$byClass["Specificity"],
                                                  Accuracy = cm$overall["Accuracy"],
                                                  F1 = cm$byClass["F1"]))
  }
  
  optimal_thresh <- performance_at_thresholds$Threshold[
    which.max(performance_at_thresholds$Sensitivity + performance_at_thresholds$Specificity)
  ]
  
  write.csv(performance_at_thresholds, paste0(best_model_name, "_threshold_analysis.csv"), row.names = FALSE)
  
  final_pred_labels <- ifelse(test_probs > optimal_thresh, "Yes", "No")
  final_pred_labels <- factor(final_pred_labels, levels = c("No", "Yes"))
  final_cm <- confusionMatrix(final_pred_labels, testdata$Result, positive = "Yes")
  
  p_threshold <- ggplot(performance_at_thresholds, aes(x = Threshold)) +
    geom_line(aes(y = Sensitivity, color = "Sensitivity"), size = 1.2) +
    geom_line(aes(y = Specificity, color = "Specificity"), size = 1.2) +
    geom_line(aes(y = Accuracy, color = "Accuracy"), size = 1.2) +
    geom_vline(xintercept = optimal_thresh, linetype = "dashed", color = "red") +
    annotate("text", x = optimal_thresh, y = 0.9, 
             label = paste("Optimal Threshold =", round(optimal_thresh, 3)),
             hjust = -0.1, color = "red") +
    labs(title = paste("Threshold Analysis -", best_model_name),
         x = "Threshold", y = "Performance") +
    scale_color_manual(values = c("Sensitivity" = "blue", 
                                  "Specificity" = "green", 
                                  "Accuracy" = "purple")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  pdf(paste0(best_model_name, "_threshold_analysis.pdf"), 8, 6)
  print(p_threshold)
  dev.off()
}

# 14. SHAP explanation
if (is.null(best_model)) {
  stop("Cannot perform SHAP analysis, best model is NULL")
}

explainer_data <- traindata_smote[, -1]
background_data <- testdata[, -1]

if (!is.matrix(explainer_data)) {
  explainer_data <- as.matrix(explainer_data)
}
if (!is.matrix(background_data)) {
  background_data <- as.matrix(background_data)
}

set.seed(123)
predict_fun <- function(model, newdata) {
  predict(model, newdata = as.data.frame(newdata), type = "prob")$Yes
}

explain_kernel <- kernelshap(best_model, X = explainer_data, bg_X = background_data, 
                             predict_fun = predict_fun)
shap_value <- shapviz(explain_kernel, X_pred = as.data.frame(explainer_data))

pdf(paste0("SHAP_", best_model_name, "_force_plot.pdf"), 8, 5)
sv_force(shap_value$Yes, row_id = 1, size = 10) +
  ggtitle(label = paste0(best_model_name, " - SHAP Force Plot")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

pdf(paste0("SHAP_", best_model_name, "_beeswarm.pdf"), 8, 6)
sv_importance(shap_value$Yes, kind = "beeswarm",
              viridis_args = list(begin = 0.25, end = 0.85, option = "B"),
              show_numbers = FALSE) +
  ggtitle(label = paste0(best_model_name, " - SHAP Feature Importance")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

pdf(paste0("SHAP_", best_model_name, "_barplot.pdf"), 8, 6)
sv_importance(shap_value$Yes, kind = "bar", show_numbers = FALSE,
              fill = "#fca50a") +
  theme_bw() +
  ggtitle(label = paste0(best_model_name, " - SHAP Feature Importance")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

feature_names <- colnames(explainer_data)

for (feature in feature_names) {
  pdf(paste0("SHAP_", best_model_name, "_dependence_", feature, ".pdf"), 7, 6)
  print(
    sv_dependence(shap_value$Yes, v = feature,
                  color = "#3b528b") +
      theme_bw() +
      ggtitle(label = paste0(best_model_name, " - SHAP Dependence: ", feature)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )
  dev.off()
}

pdf(paste0("SHAP_", best_model_name, "_waterfall.pdf"), 7, 6)
sv_waterfall(shap_value$Yes, row_id = 1,
             fill_colors = c("#f7d13d", "#a52c60")) +
  theme_bw() +
  ggtitle(label = paste0(best_model_name, " - SHAP Waterfall Plot")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()