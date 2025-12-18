setwd("C:/Users/86139/Desktop/机器学习/20250913 7模型改正错误/TA R VS NR/LIP")

library(caret)
library(glmnet)
library(pls)
library(xgboost)
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(kernlab)
library(gbm)
library(naivebayes)
library(rpart)
library(party)

full_data <- read.csv('TA R vs NR LIP.csv', header = TRUE)

colnames(full_data)[1] <- "Group"
full_data$Group <- ifelse(full_data$Group == 0, "Control", "Case")
full_data$Group <- as.factor(full_data$Group)

set.seed(12345)
train_idx <- createDataPartition(full_data$Group, p = 0.7, list = FALSE)
train_data <- full_data[train_idx, ]
test_data <- full_data[-train_idx, ]

cv_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  allowParallel = FALSE,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter = TRUE
)

model_settings <- data.frame(
  AlgorithmName = c("Lasso", "ElasticNet", "LogisticRegression", 
                    "RandomForest", "XGBoost", "DecisionTree"),
  Implementation = c("glmnet", "glmnet", "glm", 
                     "rf", "xgbTree", "rpart")
)

training_times <- numeric()
feature_importance_df <- data.frame()
gene_lists <- list()
top_features_per_model <- 20

compute_feature_importance <- function(model, model_name, feature_names) {
  tryCatch({
    if (model_name %in% c("RandomForest")) {
      imp <- varImp(model)$importance
      if (!is.null(imp)) {
        imp_df <- data.frame(Feature = rownames(imp), Overall = imp$Overall, 
                             check.names = FALSE, stringsAsFactors = FALSE)
        return(imp_df)
      }
    } 
    else if (model_name %in% c("XGBoost")) {
      xgb_imp <- xgb.importance(model = model$finalModel)
      if (!is.null(xgb_imp)) {
        imp_df <- data.frame(Feature = xgb_imp$Feature, Overall = xgb_imp$Gain,
                             check.names = FALSE, stringsAsFactors = FALSE)
        return(imp_df)
      }
    }
    else if (model_name %in% c("Lasso", "ElasticNet")) {
      coefficients <- as.matrix(coef(model$finalModel, s = model$bestTune$lambda))
      coefficients <- coefficients[-1, , drop = FALSE]
      imp_df <- data.frame(Feature = rownames(coefficients), 
                           Overall = abs(coefficients[, 1]),
                           check.names = FALSE, stringsAsFactors = FALSE)
      return(imp_df)
    }
    else if (model_name %in% c("LogisticRegression")) {
      imp <- varImp(model)$importance
      if (!is.null(imp)) {
        imp_df <- data.frame(Feature = rownames(imp), Overall = imp$Overall,
                             check.names = FALSE, stringsAsFactors = FALSE)
        return(imp_df)
      }
    }
    else if (model_name %in% c("DecisionTree")) {
      imp <- varImp(model)$importance
      if (!is.null(imp)) {
        imp_df <- data.frame(Feature = rownames(imp), Overall = imp$Overall,
                             check.names = FALSE, stringsAsFactors = FALSE)
        return(imp_df)
      }
    }
    
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

total_start <- Sys.time()

feature_names <- colnames(train_data)[-1]

for (idx in 1:nrow(model_settings)) {
  model_name <- as.character(model_settings$AlgorithmName[idx])
  method_name <- as.character(model_settings$Implementation[idx])
  
  start_time <- Sys.time()
  
  model <- tryCatch({
    if (method_name == "glmnet") {
      if (model_name == "Lasso") {
        tuneGrid <- expand.grid(alpha = 1, lambda = seq(0.001, 0.1, length = 10))
      } else if (model_name == "ElasticNet") {
        tuneGrid <- expand.grid(alpha = 0.5, lambda = seq(0.001, 0.1, length = 10))
      }
      train(Group ~ ., data = train_data, method = method_name,
            trControl = cv_control, metric = "ROC", tuneGrid = tuneGrid,
            preProcess = c("center", "scale"))
    }
    else if (method_name == "glm") {
      train(Group ~ ., data = train_data, method = method_name,
            trControl = cv_control, metric = "ROC",
            preProcess = c("center", "scale"), 
            family = binomial(link = "logit"))
    }
    else if (method_name == "rf") {
      train(Group ~ ., data = train_data, method = method_name,
            trControl = cv_control, metric = "ROC",
            preProcess = c("center", "scale"), tuneLength = 3)
    }
    else if (method_name == "xgbTree") {
      tuneGrid <- expand.grid(
        nrounds = 100, max_depth = 6, eta = 0.3,
        gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8
      )
      train(Group ~ ., data = train_data, method = method_name,
            trControl = cv_control, metric = "ROC", tuneGrid = tuneGrid,
            preProcess = c("center", "scale"))
    }
    else if (method_name == "rpart") {
      train(Group ~ ., data = train_data, method = method_name,
            trControl = cv_control, metric = "ROC",
            preProcess = c("center", "scale"), tuneLength = 5)
    }
    else {
      return(NULL)
    }
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(model)) {
    end_time <- Sys.time()
    training_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    training_times[model_name] <- training_time
    
    importance_scores <- compute_feature_importance(model, model_name, feature_names)
    if (!is.null(importance_scores) && nrow(importance_scores) > 0) {
      importance_scores$Model <- model_name
      feature_importance_df <- rbind(
        feature_importance_df, 
        importance_scores[, c("Feature", "Overall", "Model")]
      )
      
      top_features <- importance_scores %>%
        arrange(desc(Overall)) %>%
        head(top_features_per_model) %>%
        pull(Feature)
      gene_lists[[model_name]] <- top_features
      
      write.table(top_features, 
                  file = sprintf("top20_metabolites_%s.txt", model_name),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
}

if (length(training_times) > 0) {
  write.csv(data.frame(
    Model = names(training_times),
    Seconds = unname(training_times)),
    "1_Model_training_times.csv", row.names = FALSE)
}

if (nrow(feature_importance_df) > 0) {
  write.csv(feature_importance_df, "2_Feature_importance_RAW.csv", row.names = FALSE)
}

if (length(gene_lists) > 0) {
  feature_counts <- table(unlist(gene_lists))
  feature_counts_df <- data.frame(
    Feature = names(feature_counts),
    Count = as.numeric(feature_counts)
  ) %>% arrange(desc(Count))
  
  write.csv(feature_counts_df, "3_Feature_selection_counts.csv", row.names = FALSE)
}

if (length(gene_lists) > 0) {
  feature_counts <- table(unlist(gene_lists))
  feature_counts_df <- data.frame(
    Feature = names(feature_counts),
    Count = as.numeric(feature_counts)
  ) %>% arrange(desc(Count))
  
  feature_counts_filtered <- feature_counts_df %>% filter(Count >= 2)
  
  if (nrow(feature_counts_filtered) > 0) {
    ggplot(feature_counts_filtered, aes(x = reorder(Feature, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "orange") +
      coord_flip() +
      labs(title = "Metabolite Selection Frequency Across Models",
           x = "Metabolites", y = "Number of Models Selecting the Metabolite") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
    ggsave("4_metabolite_selection_by_counts.pdf", width = 10, height = 8)
  }
}

if (length(gene_lists) > 0) {
  valid_models <- names(gene_lists)[sapply(gene_lists, length) > 0]
  if (length(valid_models) > 1) {
    all_features <- unique(unlist(gene_lists[valid_models]))
    binary_matrix <- matrix(0, nrow = length(all_features), ncol = length(valid_models))
    rownames(binary_matrix) <- all_features
    colnames(binary_matrix) <- valid_models
    
    for (i in 1:length(valid_models)) {
      binary_matrix[all_features %in% gene_lists[[valid_models[i]]], i] <- 1
    }
    
    binary_df <- as.data.frame(binary_matrix)
    pdf("5_metabolite_selection_by_upset_plot.pdf", width = 12, height = 8)
    print(upset(binary_df, 
                nsets = length(valid_models),
                nintersects = 20,
                sets = valid_models,
                mb.ratio = c(0.6, 0.4),
                order.by = "freq",
                decreasing = TRUE,
                mainbar.y.label = "Number of Intersected Metabolites",
                sets.x.label = "Metabolites in Each Model",
                text.scale = c(1.3, 1.3, 1, 1, 2, 0.75)))
    dev.off()
  }
}

if (nrow(feature_importance_df) > 0) {
  top_features <- feature_importance_df %>%
    group_by(Feature) %>%
    summarise(mean_imp = mean(Overall, na.rm = TRUE)) %>%
    arrange(desc(mean_imp)) %>%
    head(30) %>%
    pull(Feature)
  
  if (length(top_features) > 0) {
    imp_subset <- feature_importance_df %>% filter(Feature %in% top_features)
    
    ggplot(imp_subset, aes(x = reorder(Feature, Overall, mean), y = Overall, color = Model)) +
      geom_point(position = position_jitter(width = 0.2), alpha = 0.7) +
      coord_flip() +
      labs(title = "Feature Importance by Model (Top 30 Features)",
           x = "Metabolites", y = "Importance Score") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8),
            legend.position = "bottom")
    ggsave("6_feature_importance_by_model.pdf", width = 12, height = 10)
  }
}