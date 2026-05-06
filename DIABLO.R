library(mixOmics)
library(tidyverse)
library(doParallel)
library(foreach)
library(igraph)
library(pheatmap)

set.seed(123)
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

common_samps <- Reduce(intersect, list(rownames(meta), rownames(metab), rownames(lipid)))
metab <- metab[common_samps, , drop = FALSE]
lipid <- lipid[common_samps, , drop = FALSE]
meta <- meta[common_samps, , drop = FALSE]

Y <- factor(meta$group)
cat("Group information:\n")
print(table(Y))

metab_scaled <- scale(metab) %>% as.data.frame()
lipid_scaled <- scale(lipid) %>% as.data.frame()

calc_cv <- function(data) {
  apply(data, 2, function(x) {
    if(sd(x) == 0) return(0)
    sd(x) / mean(abs(x), na.rm = TRUE)
  })
}

metab_cv <- calc_cv(metab_scaled)
lipid_cv <- calc_cv(lipid_scaled)

metab_keep <- names(sort(metab_cv, decreasing = TRUE)[1:floor(0.9 * length(metab_cv))])
lipid_keep <- names(sort(lipid_cv, decreasing = TRUE)[1:floor(0.9 * length(lipid_cv))])

metab_filtered <- metab_scaled[, metab_keep, drop = FALSE]
lipid_filtered <- lipid_scaled[, lipid_keep, drop = FALSE]

cat("Filtered features - Metabolites:", ncol(metab_filtered), "Lipids:", ncol(lipid_filtered), "\n")

blocks <- list(METAB = metab_filtered, LIPID = lipid_filtered)

design <- matrix(0.1, ncol = length(blocks), nrow = length(blocks))
diag(design) <- 0
colnames(design) <- rownames(design) <- names(blocks)

cat("\nDesign matrix:\n")
print(design)

cat("\n========== Step 1: Tune number of components ==========\n")

max_ncomp <- min(4, ncol(blocks$METAB), ncol(blocks$LIPID))
basic.model <- block.splsda(
  X = blocks,
  Y = Y,
  ncomp = max_ncomp,
  design = design
)

set.seed(123)
perf.diablo <- perf(
  basic.model,
  validation = 'Mfold',
  folds = 10,
  nrepeat = 10,
  dist = c('centroids.dist', 'mahalanobis.dist', 'max.dist'),
  progressBar = TRUE
)

pdf(file.path(output_dir, "ncomp_selection.pdf"), width = 10, height = 8)
plot(perf.diablo)
dev.off()

cat("\n=== Weighted vote recommended components ===\n")
print(perf.diablo$choice.ncomp$WeightedVote)

ncomp_optimal_centroids <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
cat("\nRecommended components based on centroids.dist BER:", ncomp_optimal_centroids, "\n")

ncomp_values <- c(
  perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"],
  perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "mahalanobis.dist"],
  perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "max.dist"]
)
ncomp_optimal <- as.numeric(names(sort(table(ncomp_values), decreasing = TRUE)[1]))
cat("Recommended components based on majority vote:", ncomp_optimal, "\n")

cat("\n=== Actual error rates ===\n")

if (!is.null(perf.diablo$error.rate) && length(perf.diablo$error.rate) > 0) {
  if (!is.null(perf.diablo$error.rate$centroids.dist)) {
    cat("\ncentroids.dist error rate:\n")
    print(perf.diablo$error.rate$centroids.dist)
  }
  
  if (!is.null(perf.diablo$error.rate$mahalanobis.dist)) {
    cat("\nmahalanobis.dist error rate:\n")
    print(perf.diablo$error.rate$mahalanobis.dist)
  }
  
  if (!is.null(perf.diablo$error.rate$max.dist)) {
    cat("\nmax.dist error rate:\n")
    print(perf.diablo$error.rate$max.dist)
  }
} else {
  cat("Note: error.rate is empty, attempting alternative extraction\n")
  
  if (!is.null(perf.diablo$Overall)) {
    cat("\nOverall performance:\n")
    print(perf.diablo$Overall)
  }
  
  if (!is.null(perf.diablo$BER)) {
    cat("\nBalanced error rate:\n")
    print(perf.diablo$BER)
  }
}

cat("\n========== Step 2: Set keepX search range ==========\n")

max_metab_keep <- min(10, ncol(blocks$METAB))
max_lipid_keep <- min(10, ncol(blocks$LIPID))

test.keepX <- list(
  METAB = sort(unique(c(2:5, seq(6, max_metab_keep, 2), max_metab_keep))),
  LIPID = sort(unique(c(2:5, seq(6, max_lipid_keep, 2), max_lipid_keep)))
)

test.keepX$METAB <- test.keepX$METAB[test.keepX$METAB <= ncol(blocks$METAB)]
test.keepX$LIPID <- test.keepX$LIPID[test.keepX$LIPID <= ncol(blocks$LIPID)]

cat("\nkeepX search range:\n")
cat("METAB:", paste(test.keepX$METAB, collapse = ", "), "\n")
cat("LIPID:", paste(test.keepX$LIPID, collapse = ", "), "\n")

cat("\n========== Step 3: Tune keepX ==========\n")

cores <- min(parallel::detectCores() - 1, 4)
cl <- makeCluster(cores)
registerDoParallel(cl)

set.seed(123)
tune.res <- tune.block.splsda(
  X = blocks,
  Y = Y,
  ncomp = ncomp_optimal,
  test.keepX = test.keepX,
  design = design,
  validation = 'Mfold',
  folds = 5,
  nrepeat = 5,
  dist = "centroids.dist",
  progressBar = TRUE
)

stopCluster(cl)
registerDoSEQ()

optimal.keepX <- tune.res$choice.keepX
cat("\nOptimal keepX:\n")
print(optimal.keepX)

optimal.keepX <- lapply(names(optimal.keepX), function(b) {
  pmin(optimal.keepX[[b]], ncol(blocks[[b]]))
}) %>% setNames(names(optimal.keepX))

cat("\n========== Step 4: Train temporary model for performance evaluation ==========\n")

adjusted.ncomp <- min(ncomp_optimal, 
                      min(sapply(blocks, ncol)),
                      min(sapply(optimal.keepX, length)))

set.seed(123)
temp.model <- block.splsda(
  X = blocks,
  Y = Y,
  ncomp = adjusted.ncomp,
  keepX = optimal.keepX,
  design = design
)

cat("\n========== Step 5: Evaluate optimal parameters using perf() ==========\n")

use_parallel <- TRUE
if (use_parallel) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}

set.seed(123)
perf.res <- perf(
  temp.model,
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

cat("\n========== Step 6: Train final model ==========\n")

set.seed(123)
diablo.fit <- block.splsda(
  X = blocks,
  Y = Y,
  ncomp = adjusted.ncomp,
  keepX = optimal.keepX,
  design = design
)

cat("\nFinal model summary:\n")
print(diablo.fit)

saveRDS(diablo.fit, file.path(output_dir, "diablo_final_model.rds"))

get_selected_features <- function(fit, block, comp) {
  selected <- selectVar(fit, block = block, comp = comp)
  
  if (!is.null(selected$value)) {
    loading_vec <- selected$value$value.var
    data.frame(
      Feature = names(loading_vec),
      Loading = as.numeric(loading_vec),
      Component = paste0("comp", comp),
      Block = block,
      stringsAsFactors = FALSE
    )
  }
}

sel.metab.comp1 <- get_selected_features(diablo.fit, "METAB", 1)
sel.metab.comp2 <- get_selected_features(diablo.fit, "METAB", 2)
sel.lipid.comp1 <- get_selected_features(diablo.fit, "LIPID", 1)
sel.lipid.comp2 <- get_selected_features(diablo.fit, "LIPID", 2)

sel_all_df <- do.call(rbind, list(
  sel.metab.comp1, sel.metab.comp2, sel.lipid.comp1, sel.lipid.comp2
))

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

final_metab_acc <- 1 - perf.res$error.rate$METAB[adjusted.ncomp, 1]
final_lipid_acc <- 1 - perf.res$error.rate$LIPID[adjusted.ncomp, 1]
final_weighted_acc <- 1 - weighted_ber[adjusted.ncomp, 1]


if (final_weighted_acc >= 0.85) {
  performance_grade <- "Excellent"
} else if (final_weighted_acc >= 0.75) {
  performance_grade <- "Good"
} else if (final_weighted_acc >= 0.65) {
  performance_grade <- "Fair"
} else {
  performance_grade <- "Needs improvement"
}

for (strategy in names(perf.res$choice.ncomp)) {
  cat(sprintf("  %s: %d components\n", strategy, perf.res$choice.ncomp[[strategy]][1, 1]))
}

if (adjusted.ncomp >= 2) {
  pdf(file.path(output_dir, "diablo_plotIndiv.pdf"), width = 7, height = 6)
  plotIndiv(diablo.fit, comp = c(1,2), group = Y, legend = TRUE, ind.names = FALSE, 
            title = "DIABLO: sample plot (comp1 vs comp2)")
  dev.off()
} else {
  pdf(file.path(output_dir, "diablo_plotIndiv_comp1.pdf"), width = 7, height = 6)
  plotIndiv(diablo.fit, comp = 1, group = Y, legend = TRUE, ind.names = FALSE,
            title = "DIABLO: sample plot (Component 1)")
  dev.off()
  cat("Note: Model has only 1 component, 1D plot created\n")
}

if (!is.null(diablo.fit)) {
  pdf(file.path(output_dir, "circos_plot.pdf"), width = 9, height = 9)
  tryCatch({
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
  }, error = function(e) {
    cat("Circos plot failed:", e$message, "\n")
  })
  dev.off()
  cat("Circos plot saved: circos_plot.pdf\n")
} else {
  cat("Model does not exist, cannot create Circos plot\n")
}

plot_loadings <- function(fit, block, comp, filename) {
  model_ncomp <- if (is.null(fit$ncomp)) {
    if (!is.null(fit$loadings[[block]])) ncol(fit$loadings[[block]]) else 0
  } else if (length(fit$ncomp) > 1) {
    if (block %in% names(fit$ncomp)) fit$ncomp[block] else fit$ncomp[1]
  } else {
    fit$ncomp
  }
  
  if (comp > model_ncomp) {
    cat("Skipping", filename, "- component", comp, "does not exist\n")
    return(NULL)
  }
  
  if (!is.null(fit$keepX[[block]])) {
    if (length(fit$keepX[[block]]) < comp || fit$keepX[[block]][comp] == 0) {
      cat("Skipping", filename, "- no features selected\n")
      return(NULL)
    }
  }
  
  pdf(file.path(output_dir, filename), width = 8, height = 6)
  tryCatch({
    plotLoadings(fit, comp = comp, block = block, method = "mean", contrib = "max")
  }, error = function(e) {
    cat("Error plotting loadings:", filename, "-", e$message, "\n")
  })
  dev.off()
  cat("Loadings plot saved:", filename, "\n")
}

plot_loadings(diablo.fit, "METAB", 1, "loadings_METAB_comp1.pdf")
plot_loadings(diablo.fit, "METAB", 2, "loadings_METAB_comp2.pdf")
plot_loadings(diablo.fit, "LIPID", 1, "loadings_LIPID_comp1.pdf")
plot_loadings(diablo.fit, "LIPID", 2, "loadings_LIPID_comp2.pdf")

pdf(file.path(output_dir, "plotDiablo.pdf"), width = 8, height = 6)
plotDiablo(diablo.fit, ncomp = 1)
dev.off()

pdf(file.path(output_dir, "plotArrow.pdf"), width = 8, height = 6)
plotArrow(diablo.fit, ind.names = FALSE, legend = TRUE)
dev.off()

pdf(file.path(output_dir, "heatmap_cimDiablo_comp1.pdf"), width = 12, height = 10)
cimDiablo(diablo.fit, comp = 1, 
          title = "DIABLO Heatmap - Component 1",
          legend.position = "right")
dev.off()

if (adjusted.ncomp >= 2) {
  pdf(file.path(output_dir, "heatmap_cimDiablo_comp2.pdf"), width = 12, height = 10)
  cimDiablo(diablo.fit, comp = 2, 
            title = "DIABLO Heatmap - Component 2",
            legend.position = "right")
  dev.off()
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

cat("Feature extraction and loadings table saved. Completed.\n")

selected_features_for_save <- list(
  METAB_comp1 = sel.metab.comp1,
  METAB_comp2 = sel.metab.comp2,
  LIPID_comp1 = sel.lipid.comp1,
  LIPID_comp2 = sel.lipid.comp2
)

final_results <- list(
  diablo_model = diablo.fit,
  tuning = tune.res,
  performance = perf.res,
  selected_features = selected_features_for_save,
  parameters = list(
    optimal_ncomp = ncomp_optimal,
    optimal_keepX = optimal.keepX,
    adjusted_ncomp = adjusted.ncomp,
    adjusted_keepX = optimal.keepX,
    design_matrix = design
  )
)

saveRDS(final_results, file.path(output_dir, "diablo_complete_results.rds"))