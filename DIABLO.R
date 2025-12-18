library(mixOmics)
library(tidyverse)
library(doParallel)
library(foreach)
library(igraph)
library(pheatmap)

set.seed(2025)

# ---------- 1. Load data ----------
input_dir <- "C:/Users/86139/Desktop/机器学习/R DIABLO(1)/R DIABLO 20250915 加KEGG DEM/HAPC VS NHAPC"
output_dir <- file.path(input_dir, "output/DIABLO  20251014")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metab <- read.csv(file.path(input_dir, "met_data v2.csv"), 
                  check.names = FALSE, 
                  row.names = 1)
lipid <- read.csv(file.path(input_dir, "lip_data v2.csv"), 
                  check.names = FALSE, 
                  row.names = 1)
meta <- read.csv(file.path(input_dir, "sample_info.csv"), 
                 check.names = FALSE, 
                 row.names = 1)

# ---------- 2. Data preprocessing ----------
common_samps <- Reduce(intersect, list(rownames(meta), rownames(metab), rownames(lipid)))
metab  <- metab[common_samps, , drop = FALSE]
lipid <- lipid[common_samps, , drop = FALSE]
meta  <- meta[common_samps, , drop = FALSE]

Y <- factor(meta$group)

metab_scaled <- scale(metab) %>% as.data.frame()
lipid_scaled <- scale(lipid) %>% as.data.frame()

# ---------- 3. Feature pre-screening ----------
calc_cv <- function(data) {
  apply(data, 2, function(x) {
    if(sd(x) == 0) return(0)
    sd(x) / mean(abs(x), na.rm = TRUE)
  })
}

metab_cv <- calc_cv(metab_scaled)
lipid_cv <- calc_cv(lipid_scaled)

metab_keep <- names(sort(metab_cv, decreasing = TRUE)[1:floor(0.6 * length(metab_cv))])
lipid_keep <- names(sort(lipid_cv, decreasing = TRUE)[1:floor(0.6 * length(lipid_cv))])

metab_filtered <- metab_scaled[, metab_keep, drop = FALSE]
lipid_filtered <- lipid_scaled[, lipid_keep, drop = FALSE]

blocks <- list(METAB = metab_filtered, LIPID = lipid_filtered)

# ---------- 4. Model tuning ----------
design <- matrix(0.8, ncol = length(blocks), nrow = length(blocks))
diag(design) <- 0
colnames(design) <- rownames(design) <- names(blocks)

max_metab_features <- min(8, ncol(metab_filtered))
max_lipid_features <- min(8, ncol(lipid_filtered))

list.keepX <- list(
  METAB = sort(unique(c(6, 7, 8, max_metab_features))),
  LIPID = sort(unique(c(6, 7, 8, max_lipid_features)))
)

ncomp_try <- min(2, max_metab_features, max_lipid_features)

use_parallel <- TRUE
if (use_parallel) {
  cores <- min(parallel::detectCores() - 1, 4)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}

set.seed(2025)
tune.res <- tune.block.splsda(
  X = blocks,
  Y = Y,
  ncomp = ncomp_try,
  test.keepX = list.keepX,
  design = design,
  validation = 'Mfold',
  folds = 5,
  nrepeat = 5,
  dist = "centroids.dist",
  progressBar = TRUE
)

if (use_parallel) { 
  stopCluster(cl)
  registerDoSEQ() 
}

optimal.keepX <- tune.res$choice.keepX
optimal.ncomp <- length(optimal.keepX$METAB)

# ---------- 5. Train final model ----------
adjust_keepX <- function(keepX, X) {
  lapply(names(keepX), function(block_name) {
    pmin(keepX[[block_name]], ncol(X[[block_name]]))
  }) %>% setNames(names(keepX))
}

adjusted.keepX <- adjust_keepX(optimal.keepX, blocks)
adjusted.ncomp <- min(optimal.ncomp, 
                      min(sapply(blocks, ncol)),
                      min(sapply(adjusted.keepX, length)))

if (adjusted.ncomp < optimal.ncomp) {
  adjusted.keepX <- lapply(adjusted.keepX, function(x) x[1:adjusted.ncomp])
}

set.seed(2025)
diablo.fit <- block.splsda(
  X = blocks,
  Y = Y,
  ncomp = adjusted.ncomp,
  keepX = adjusted.keepX,
  design = design
)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
saveRDS(diablo.fit, file.path(output_dir, "diablo_final_model.rds"))

# ---------- 6. Feature extraction ----------
get_selected_features <- function(fit, block, comp) {
  var_info <- tryCatch({
    selectVar(fit, block = block, comp = comp)
  }, error = function(e) {
    NULL
  })
  
  build_df <- function(loading_vec, comp, block) {
    if (is.null(loading_vec) || length(loading_vec) == 0) return(NULL)
    if (is.null(names(loading_vec))) names(loading_vec) <- paste0("V", seq_along(loading_vec))
    data.frame(
      Feature = names(loading_vec),
      Loading = as.numeric(loading_vec),
      Component = paste0("comp", comp),
      Block = block,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  
  if (!is.null(var_info)) {
    if (!is.null(var_info$name)) {
      loading_vec <- NULL
      if (!is.null(var_info$value)) {
        if (is.list(var_info$value) && !is.null(var_info$value$value.var)) {
          loading_vec <- var_info$value$value.var
        } else if (is.data.frame(var_info$value) || is.matrix(var_info$value)) {
          loading_vec <- as.numeric(var_info$value[, 1])
        } else {
          loading_vec <- tryCatch(as.numeric(unlist(var_info$value)), error = function(e) NULL)
        }
      }
      if (!is.null(loading_vec) && length(loading_vec) == length(var_info$name)) {
        names(loading_vec) <- var_info$name
        return(build_df(loading_vec, comp, block))
      }
      if (!is.null(loading_vec) && length(loading_vec) != length(var_info$name)) {
        names(loading_vec) <- var_info$name[seq_len(length(loading_vec))]
        return(build_df(loading_vec, comp, block))
      }
    }
    
    if (is.matrix(var_info) || is.data.frame(var_info)) {
      if (ncol(var_info) >= 1) {
        loading_vec <- as.numeric(var_info[, 1])
        names(loading_vec) <- rownames(var_info)
        return(build_df(loading_vec, comp, block))
      }
    }
  }
  
  load_mat <- tryCatch(fit$loadings[[block]], error = function(e) NULL)
  if (!is.null(load_mat) && is.matrix(load_mat) && ncol(load_mat) >= comp) {
    loading_vec <- load_mat[, comp]
    names(loading_vec) <- rownames(load_mat)
    return(build_df(loading_vec, comp, block))
  }
  
  return(NULL)
}

sel.metab.comp1 <- get_selected_features(diablo.fit, "METAB", 1)
sel.metab.comp2 <- get_selected_features(diablo.fit, "METAB", 2)
sel.lipid.comp1 <- get_selected_features(diablo.fit, "LIPID", 1)
sel.lipid.comp2 <- get_selected_features(diablo.fit, "LIPID", 2)

all_sel_list <- list(
  sel.metab.comp1 = sel.metab.comp1,
  sel.metab.comp2 = sel.metab.comp2,
  sel.lipid.comp1 = sel.lipid.comp1,
  sel.lipid.comp2 = sel.lipid.comp2
)

sel_all_df <- do.call(rbind, lapply(all_sel_list, function(x) { if(is.null(x)) NULL else x }))

if (!is.null(sel_all_df) && nrow(sel_all_df) > 0) {
  write.csv(sel_all_df, file.path(output_dir, "selected_features_comp1_comp2.csv"), row.names = FALSE)
}

if (!is.null(sel.metab.comp1)) write.csv(sel.metab.comp1, file.path(output_dir, "selected_METAB_comp1.csv"), row.names = FALSE)
if (!is.null(sel.metab.comp2)) write.csv(sel.metab.comp2, file.path(output_dir, "selected_METAB_comp2.csv"), row.names = FALSE)
if (!is.null(sel.lipid.comp1)) write.csv(sel.lipid.comp1, file.path(output_dir, "selected_LIPID_comp1.csv"), row.names = FALSE)
if (!is.null(sel.lipid.comp2)) write.csv(sel.lipid.comp2, file.path(output_dir, "selected_LIPID_comp2.csv"), row.names = FALSE)

# ---------- 7. Performance evaluation ----------
if (use_parallel) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}

set.seed(2025)
perf.res <- perf(
  diablo.fit,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  dist = "centroids.dist",
  progressBar = TRUE
)

if (use_parallel) { 
  stopCluster(cl)
  registerDoSEQ() 
}

n_metab_features <- ncol(blocks$METAB)
n_lipid_features <- ncol(blocks$LIPID)
total_features <- n_metab_features + n_lipid_features

weight_metab <- n_metab_features / total_features
weight_lipid <- n_lipid_features / total_features

weighted_ber <- matrix(0, nrow = adjusted.ncomp, ncol = 1)
rownames(weighted_ber) <- paste0("comp", 1:adjusted.ncomp)
colnames(weighted_ber) <- "Weighted_BER"

for (comp in 1:adjusted.ncomp) {
  metab_error <- perf.res$error.rate$METAB[comp, 1]
  lipid_error <- perf.res$error.rate$LIPID[comp, 1]
  weighted_ber[comp, 1] <- metab_error * weight_metab + lipid_error * weight_lipid
}

performance_details <- data.frame(
  Component = rep(1:adjusted.ncomp, 2),
  Block = rep(c("METAB", "LIPID"), each = adjusted.ncomp),
  ErrorRate = c(perf.res$error.rate$METAB[, 1], perf.res$error.rate$LIPID[, 1]),
  ErrorRate_SD = c(perf.res$error.rate.sd$METAB[, 1], perf.res$error.rate.sd$LIPID[, 1]),
  Weighted_BER = rep(weighted_ber[, 1], 2),
  stringsAsFactors = FALSE
)

performance_details$Accuracy <- 1 - performance_details$ErrorRate
performance_details$Weighted_Accuracy <- 1 - performance_details$Weighted_BER

write.csv(performance_details, file.path(output_dir, "detailed_performance_metrics.csv"), row.names = FALSE)

weighted_ber_df <- data.frame(
  Component = 1:adjusted.ncomp,
  Weighted_BER = weighted_ber[, 1],
  Weighted_Accuracy = 1 - weighted_ber[, 1],
  stringsAsFactors = FALSE
)
write.csv(weighted_ber_df, file.path(output_dir, "weighted_performance.csv"), row.names = FALSE)

# ---------- 8. Visualization ----------
pdf(file.path(output_dir, "block_performance_comparison.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))

matplot(1:adjusted.ncomp, cbind(perf.res$error.rate$METAB[, 1], 
                                perf.res$error.rate$LIPID[, 1],
                                weighted_ber[, 1]), 
        type = "b", pch = 1:3, col = 1:3, lwd = 2,
        xlab = "Component", ylab = "Error Rate", 
        main = "Error Rate by Block and Component")
legend("topright", legend = c("METAB", "LIPID", "Weighted Average"), 
       pch = 1:3, col = 1:3, lwd = 2)
grid()

matplot(1:adjusted.ncomp, cbind(1 - perf.res$error.rate$METAB[, 1], 
                                1 - perf.res$error.rate$LIPID[, 1],
                                1 - weighted_ber[, 1]), 
        type = "b", pch = 1:3, col = 1:3, lwd = 2,
        xlab = "Component", ylab = "Accuracy", 
        main = "Accuracy by Block and Component")
legend("bottomright", legend = c("METAB", "LIPID", "Weighted Average"), 
       pch = 1:3, col = 1:3, lwd = 2)
grid()
dev.off()

pdf(file.path(output_dir, "performance_with_ci.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

metab_lower <- pmax(0, perf.res$error.rate$METAB[, 1] - perf.res$error.rate.sd$METAB[, 1])
metab_upper <- pmin(1, perf.res$error.rate$METAB[, 1] + perf.res$error.rate.sd$METAB[, 1])

plot(1:adjusted.ncomp, perf.res$error.rate$METAB[, 1], type = "b", 
     ylim = c(0, max(metab_upper) * 1.1), pch = 19, col = "red", lwd = 2,
     xlab = "Component", ylab = "Error Rate", 
     main = "METAB Error Rate with Confidence Intervals")
arrows(1:adjusted.ncomp, metab_lower, 1:adjusted.ncomp, metab_upper, 
       angle = 90, code = 3, length = 0.1, col = "red")
grid()

lipid_lower <- pmax(0, perf.res$error.rate$LIPID[, 1] - perf.res$error.rate.sd$LIPID[, 1])
lipid_upper <- pmin(1, perf.res$error.rate$LIPID[, 1] + perf.res$error.rate.sd$LIPID[, 1])

plot(1:adjusted.ncomp, perf.res$error.rate$LIPID[, 1], type = "b", 
     ylim = c(0, max(lipid_upper) * 1.1), pch = 19, col = "blue", lwd = 2,
     xlab = "Component", ylab = "Error Rate", 
     main = "LIPID Error Rate with Confidence Intervals")
arrows(1:adjusted.ncomp, lipid_lower, 1:adjusted.ncomp, lipid_upper, 
       angle = 90, code = 3, length = 0.1, col = "blue")
grid()
dev.off()

pdf(file.path(output_dir, "diablo_plotIndiv.pdf"), width = 7, height = 6)
plotIndiv(diablo.fit, comp = c(1,2), group = Y, legend = TRUE, ind.names = FALSE, 
          title = "DIABLO: sample plot (comp1 vs comp2)")
dev.off()

pdf(file.path(output_dir, "circos_plot.pdf"), width = 9, height = 9)
circosPlot(
  diablo.fit,
  cutoff = 0.3,
  legend = TRUE,
  line = TRUE,
  size.variables = 0.5,
  color.blocks = c(METAB = "#E41A1C", LIPID = "#377EB8"),
  track.height = 0.1,
  show.block.names = TRUE
)
dev.off()

plot_loadings <- function(fit, block, comp, filename) {
  pdf(file.path(output_dir, filename), width = 7, height = 5)
  plotLoadings(
    fit,
    comp = comp,
    block = block,
    method = "mean",
    contrib = "max"
  )
  dev.off()
}

plot_loadings(diablo.fit, "METAB", 1, "loadings_METAB_comp1.pdf")
plot_loadings(diablo.fit, "METAB", 2, "loadings_METAB_comp2.pdf")
plot_loadings(diablo.fit, "LIPID", 1, "loadings_LIPID_comp1.pdf")
plot_loadings(diablo.fit, "LIPID", 2, "loadings_LIPID_comp2.pdf")

# ---------- 9. Heatmap function ----------
make_heatmap_from_selected <- function(metab_scaled, lipid_scaled, sel_metab_df, sel_lipid_df, outname, Y) {
  if (is.null(sel_metab_df) || nrow(sel_metab_df) == 0 || is.null(sel_lipid_df) || nrow(sel_lipid_df) == 0) {
    return(NULL)
  }
  sel_metab_feats <- sel_metab_df$Feature
  sel_lipid_feats <- sel_lipid_df$Feature
  heat.metab <- metab_scaled[, intersect(colnames(metab_scaled), sel_metab_feats), drop = FALSE]
  heat.lipid <- lipid_scaled[, intersect(colnames(lipid_scaled), sel_lipid_feats), drop = FALSE]
  heat.data <- cbind(heat.metab, heat.lipid)
  if (ncol(heat.data) == 0) {
    return(NULL)
  }
  
  annotation_col <- data.frame(Group = Y)
  rownames(annotation_col) <- rownames(heat.data)
  annotation_row <- data.frame(Block = c(rep("Metabolite", ncol(heat.metab)), rep("Lipid", ncol(heat.lipid))))
  rownames(annotation_row) <- colnames(heat.data)
  
  pdf(file.path(output_dir, outname), width = 10, height = 8)
  pheatmap(t(heat.data),
           annotation_col = annotation_col,
           annotation_row = annotation_row,
           scale = "row",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           show_colnames = FALSE,
           fontsize_row = 8)
  dev.off()
  return(TRUE)
}

if (exists("metab_scaled") && exists("lipid_scaled")) {
  make_heatmap_from_selected(metab_scaled, lipid_scaled, sel.metab.comp1, sel.lipid.comp1, "heatmap_selected_comp1.pdf", Y)
  make_heatmap_from_selected(metab_scaled, lipid_scaled, sel.metab.comp2, sel.lipid.comp2, "heatmap_selected_comp2.pdf", Y)
}

save_loadings_table <- function(fit, block, comp, outfile) {
  load_mat <- tryCatch(fit$loadings[[block]], error = function(e) NULL)
  if (is.null(load_mat) || ncol(load_mat) < comp) return(NULL)
  df <- data.frame(Feature = rownames(load_mat), Loading = load_mat[, comp], stringsAsFactors = FALSE)
  df <- df[order(abs(df$Loading), decreasing = TRUE), ]
  write.csv(df, outfile, row.names = FALSE)
  return(df)
}

save_loadings_table(diablo.fit, "METAB", 1, file.path(output_dir, "loadings_table_METAB_comp1.csv"))
save_loadings_table(diablo.fit, "METAB", 2, file.path(output_dir, "loadings_table_METAB_comp2.csv"))
save_loadings_table(diablo.fit, "LIPID", 1, file.path(output_dir, "loadings_table_LIPID_comp1.csv"))
save_loadings_table(diablo.fit, "LIPID", 2, file.path(output_dir, "loadings_table_LIPID_comp2.csv"))

# ---------- 10. LOOCV validation ----------
loocv_diablo <- function(X, Y, design, ncomp, keepX) {
  n_samples <- nrow(X[[1]])
  predictions <- matrix(NA, nrow = n_samples, ncol = ncomp)
  rownames(predictions) <- rownames(X[[1]])
  
  for (i in 1:n_samples) {
    test_index <- i
    train_index <- setdiff(1:n_samples, test_index)
    
    X_train <- lapply(X, function(x) x[train_index, , drop = FALSE])
    X_test <- lapply(X, function(x) x[test_index, , drop = FALSE])
    Y_train <- Y[train_index]
    
    model <- block.splsda(
      X = X_train,
      Y = Y_train,
      ncomp = ncomp,
      keepX = keepX,
      design = design
    )
    
    pred <- predict(model, newdata = X_test, dist = "max.dist", ncomp = ncomp)
    
    if (!all(c("class", "max.dist") %in% c(names(pred), names(pred$class)))) {
      stop(paste("Sample", i, "prediction missing required elements"))
    }
    
    pred_results <- pred$class$max.dist
    
    if (is.list(pred_results)) {
      for (comp in 1:ncomp) {
        if (comp <= length(pred_results) && length(pred_results[[comp]]) > 0) {
          predictions[i, comp] <- pred_results[[comp]][1]
        } else {
          predictions[i, comp] <- NA
        }
      }
    } else if (is.matrix(pred_results)) {
      for (comp in 1:ncomp) {
        if (comp <= ncol(pred_results)) {
          predictions[i, comp] <- pred_results[1, comp]
        } else {
          predictions[i, comp] <- NA
        }
      }
    } else if (is.vector(pred_results)) {
      for (comp in 1:ncomp) {
        if (comp <= length(pred_results)) {
          predictions[i, comp] <- pred_results[comp]
        } else {
          predictions[i, comp] <- NA
        }
      }
    } else {
      stop(paste("Sample", i, "prediction format unknown:", class(pred_results)))
    }
  }
  
  accuracy <- sapply(1:ncomp, function(comp) {
    mean(predictions[, comp] == as.character(Y), na.rm = TRUE)
  })
  
  confusion_matrices <- lapply(1:ncomp, function(comp) {
    table(Actual = Y, Predicted = predictions[, comp])
  })
  
  return(list(predictions = predictions, 
              accuracy = accuracy, 
              confusion_matrices = confusion_matrices))
}

loocv_results <- loocv_diablo(
  X = blocks,
  Y = Y,
  design = design,
  ncomp = optimal.ncomp,
  keepX = optimal.keepX
)

write.csv(
  data.frame(
    Component = 1:optimal.ncomp,
    Accuracy = loocv_results$accuracy
  ),
  file.path(output_dir, "loocv_accuracy.csv"),
  row.names = FALSE
)

# ---------- 11. Save final results ----------
final_results <- list(
  diablo_model = diablo.fit,
  tuning = tune.res,
  performance = perf.res,
  parameters = list(
    optimal_ncomp = optimal.ncomp,
    optimal_keepX = optimal.keepX,
    adjusted_ncomp = adjusted.ncomp,
    adjusted_keepX = adjusted.keepX
  )
)

saveRDS(final_results, file.path(output_dir, "diablo_complete_results.rds"))