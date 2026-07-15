library(mixOmics)
library(tidyverse)
library(doParallel)
library(foreach)
library(igraph)
library(pheatmap)
library(doRNG)
set.seed(123)
summarize_perf_stability <- function(perf_obj) {
  stable <- perf_obj$features$stable
  if (is.null(stable)) {
    return(NULL)
  }
  unlisted_data <- unlist(stable, use.names = TRUE)
  if (length(unlisted_data) == 0) {
    return(NULL)
  }
  split_names <- strsplit(names(unlisted_data), "\\.", perl = TRUE)
  df <- data.frame(
    nrepeat = vapply(split_names, function(x) x[1], character(1)),
    block = vapply(split_names, function(x) x[2], character(1)),
    component = vapply(split_names, function(x) x[3], character(1)),
    feature = vapply(
      split_names,
      function(x) paste(x[4:length(x)], collapse = "."),
      character(1)
    ),
    stability = as.numeric(unlisted_data),
    stringsAsFactors = FALSE
  )
  agg <- aggregate(
    stability ~ block + component + feature,
    data = df,
    FUN = mean
  )
  agg <- agg[order(agg$block, agg$component, -agg$stability), ]
  rownames(agg) <- NULL
  list(raw = df, summary = agg)
}

save_feature_stability <- function(perf_obj, out_dir) {
  stab <- summarize_perf_stability(perf_obj)
  if (is.null(stab)) {
    warning("perf$features$stable is empty; stability CSV not written.")
    return(invisible(NULL))
  }
  write.csv(
    stab$raw,
    file.path(out_dir, "feature_stability_all_repeats.csv"),
    row.names = FALSE
  )
  write.csv(
    stab$summary,
    file.path(out_dir, "feature_stability_summary.csv"),
    row.names = FALSE
  )
  for (blk in unique(stab$summary$block)) {
    sub <- stab$summary[stab$summary$block == blk, , drop = FALSE]
    fname <- sprintf("feature_stability_%s.csv", blk)
    write.csv(sub, file.path(out_dir, fname), row.names = FALSE)
    cat(sprintf("\n%s stability (top 10):\n", blk))
    print(head(sub, 10))
  }
  invisible(stab)
}


input_dir <- "input_dir"
output_dir <- "output_dir"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metab <- read.csv(file.path(input_dir, "met_data_FDR.csv"),
                  check.names = FALSE,
                  row.names = 1)
lipid <- read.csv(file.path(input_dir, "lip_data_FDR.csv"),
                  check.names = FALSE,
                  row.names = 1)
meta <- read.csv(file.path(input_dir, "sample_info.csv"), 
                 check.names = FALSE, 
                 row.names = 1)

common_samps <- Reduce(intersect, list(rownames(meta), rownames(metab), rownames(lipid)))
metab  <- metab[common_samps, , drop = FALSE]
lipid <- lipid[common_samps, , drop = FALSE]
meta  <- meta[common_samps, , drop = FALSE]

Y <- factor(meta$group)
print(table(Y))

metab_scaled <- scale(metab) %>% as.data.frame()
lipid_scaled <- scale(lipid) %>% as.data.frame()
blocks <- list(METAB = metab_scaled, LIPID = lipid_scaled)
design <- matrix(0.1, ncol = length(blocks), nrow = length(blocks))
diag(design) <- 0
colnames(design) <- rownames(design) <- names(blocks)
print(design)

max_ncomp <- min(4, ncol(blocks$METAB), ncol(blocks$LIPID))
basic.model <- block.splsda(
  X = blocks,
  Y = Y,
  ncomp = max_ncomp,
  design = design
)
ncomp_optimal <- 2
set.seed(123)
perf.diablo <- perf(
  basic.model,
  validation = 'Mfold',
  folds = 5,
  nrepeat = 10,
  dist = c('centroids.dist', 'mahalanobis.dist', 'max.dist'),
  progressBar = TRUE
)

pdf(file.path(output_dir, "ncomp_selection.pdf"), width = 10, height = 8)
plot(perf.diablo)
dev.off()

print(perf.diablo$choice.ncomp$WeightedVote)

test.keepX <- list(
  METAB = c(5, 8, 10, 12),
  LIPID = c(4, 5, 6, 8, 10)
)
test.keepX$METAB <- test.keepX$METAB[test.keepX$METAB <= ncol(blocks$METAB)]
test.keepX$LIPID <- test.keepX$LIPID[test.keepX$LIPID <= ncol(blocks$LIPID)]

cat("METAB:", paste(test.keepX$METAB, collapse = ", "), "\n")
cat("LIPID:", paste(test.keepX$LIPID, collapse = ", "), "\n")

cores <- min(parallel::detectCores() - 1, 4)
cl <- makeCluster(cores)
registerDoParallel(cl)
# ==========  ==========
registerDoRNG(123)   
# =========================================
set.seed(123)
tune.res <- tune.block.splsda(
  X = blocks,
  Y = Y,
  ncomp = ncomp_optimal, 
  test.keepX = test.keepX,
  design = design,
  validation = 'Mfold',
  folds = 5,
  nrepeat = 50,
  dist = "centroids.dist",
  progressBar = TRUE
)

stopCluster(cl)
registerDoSEQ()

optimal.keepX <- tune.res$choice.keepX
print(optimal.keepX)

for (b in names(optimal.keepX)) {
  optimal.keepX[[b]] <- pmin(optimal.keepX[[b]], ncol(blocks[[b]]))
}

if (length(optimal.keepX$METAB) >= 1) {
  optimal.keepX$METAB[1] <- max(optimal.keepX$METAB[1], 5)
}
if (length(optimal.keepX$METAB) >= 2) {
  optimal.keepX$METAB[2] <- max(optimal.keepX$METAB[2], 8)
  optimal.keepX$METAB[2] <- min(optimal.keepX$METAB[2], 12)
}
if (length(optimal.keepX$LIPID) >= 2) {
  optimal.keepX$LIPID[2] <- max(optimal.keepX$LIPID[2], 5)
}
print(optimal.keepX)

keepX_df <- rbind(
  data.frame(
    block = "METAB",
    component = paste0("comp", seq_along(optimal.keepX$METAB)),
    keepX = optimal.keepX$METAB,
    stringsAsFactors = FALSE
  ),
  data.frame(
    block = "LIPID",
    component = paste0("comp", seq_along(optimal.keepX$LIPID)),
    keepX = optimal.keepX$LIPID,
    stringsAsFactors = FALSE
  )
)
write.csv(keepX_df, file.path(output_dir, "diablo_keepX_parameters.csv"), row.names = FALSE)

adjusted.ncomp <- min(ncomp_optimal, 
                      length(optimal.keepX$METAB),
                      length(optimal.keepX$LIPID))


keepX_adjusted <- list(
  METAB = optimal.keepX$METAB[1:adjusted.ncomp],
  LIPID = optimal.keepX$LIPID[1:adjusted.ncomp]
)

print(keepX_adjusted)

set.seed(123)
temp.model <- block.splsda(
  X = blocks,
  Y = Y,
  ncomp = adjusted.ncomp,
  keepX = keepX_adjusted, 
  design = design
)


use_parallel <- TRUE
if (use_parallel) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}

registerDoRNG(123)
set.seed(123)
perf.res <- perf(
  temp.model,
  validation = "Mfold",
  folds = 5,
  nrepeat = 50,
  dist = "centroids.dist",
  progressBar = TRUE
)

save_feature_stability(perf.res, output_dir)


set.seed(123)
diablo.fit <- block.splsda(
  X = blocks,
  Y = Y,
  ncomp = adjusted.ncomp,
  keepX = keepX_adjusted, 
  design = design
)

print(diablo.fit)

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
sel.metab.comp2 <- if(adjusted.ncomp >= 2) get_selected_features(diablo.fit, "METAB", 2)
sel.lipid.comp1 <- get_selected_features(diablo.fit, "LIPID", 1)
sel.lipid.comp2 <- if(adjusted.ncomp >= 2) get_selected_features(diablo.fit, "LIPID", 2)

sel_all_df <- do.call(rbind, Filter(Negate(is.null), list(
  sel.metab.comp1, sel.metab.comp2, sel.lipid.comp1, sel.lipid.comp2
)))

n_metab_features <- ncol(blocks$METAB)
n_lipid_features <- ncol(blocks$LIPID)
total_features <- n_metab_features + n_lipid_features

weight_metab <- n_metab_features / total_features
weight_lipid <- n_lipid_features / total_features

actual_ncomp <- nrow(perf.res$error.rate$METAB)

ncomp_to_use <- min(adjusted.ncomp, actual_ncomp)

weighted_ber <- matrix(0, nrow = ncomp_to_use, ncol = 1)
rownames(weighted_ber) <- paste0("comp", 1:ncomp_to_use)
colnames(weighted_ber) <- "Weighted_BER"

for (comp in 1:ncomp_to_use) {
  metab_error <- perf.res$error.rate$METAB[comp, 1]
  lipid_error <- perf.res$error.rate$LIPID[comp, 1]
  weighted_ber[comp, 1] <- metab_error * weight_metab + lipid_error * weight_lipid
}

print(weighted_ber)
print(perf.res$AveragedPredict.error.rate)
print(perf.res$WeightedVote.error.rate$centroids.dist)
print(perf.res$choice.ncomp)

pdf(file.path(output_dir, "block_performance_comparison.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))
matplot(1:ncomp_to_use, cbind(perf.res$error.rate$METAB[1:ncomp_to_use, 1], 
                              perf.res$error.rate$LIPID[1:ncomp_to_use, 1],
                              weighted_ber[1:ncomp_to_use, 1]), 
        type = "b", pch = 1:3, col = 1:3, lwd = 2,
        xlab = "Component", ylab = "Error Rate", 
        main = "Error Rate by Block and Component")
legend("topright", legend = c("METAB", "LIPID", "Weighted Average"), 
       pch = 1:3, col = 1:3, lwd = 2)
grid()

matplot(1:ncomp_to_use, cbind(1 - perf.res$error.rate$METAB[1:ncomp_to_use, 1], 
                              1 - perf.res$error.rate$LIPID[1:ncomp_to_use, 1],
                              1 - weighted_ber[1:ncomp_to_use, 1]), 
        type = "b", pch = 1:3, col = 1:3, lwd = 2,
        xlab = "Component", ylab = "Accuracy", 
        main = "Accuracy by Block and Component")
legend("bottomright", legend = c("METAB", "LIPID", "Weighted Average"), 
       pch = 1:3, col = 1:3, lwd = 2)
grid()

dev.off()

open_perf_pdf <- function(path, width = 12, height = 6) {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::CairoPDF(path, width = width, height = height, family = "Arial")
  } else {
    pdf(path, width = width, height = height, useDingbats = FALSE, family = "sans")
  }
}
open_perf_pdf(file.path(output_dir, "performance_with_ci.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2), family = "sans")

metab_lower <- pmax(0, perf.res$error.rate$METAB[1:ncomp_to_use, 1] - perf.res$error.rate.sd$METAB[1:ncomp_to_use, 1])
metab_upper <- pmin(1, perf.res$error.rate$METAB[1:ncomp_to_use, 1] + perf.res$error.rate.sd$METAB[1:ncomp_to_use, 1])

plot(1:ncomp_to_use, perf.res$error.rate$METAB[1:ncomp_to_use, 1], type = "b", 
     ylim = c(0, max(metab_upper) * 1.1), pch = 19, col = "red", lwd = 2,
     xlab = "Component", ylab = "Error Rate", 
     main = "METAB Error Rate with Confidence Intervals")
arrows(1:ncomp_to_use, metab_lower, 1:ncomp_to_use, metab_upper, 
       angle = 90, code = 3, length = 0.1, col = "red")
grid()

lipid_lower <- pmax(0, perf.res$error.rate$LIPID[1:ncomp_to_use, 1] - perf.res$error.rate.sd$LIPID[1:ncomp_to_use, 1])
lipid_upper <- pmin(1, perf.res$error.rate$LIPID[1:ncomp_to_use, 1] + perf.res$error.rate.sd$LIPID[1:ncomp_to_use, 1])

plot(1:ncomp_to_use, perf.res$error.rate$LIPID[1:ncomp_to_use, 1], type = "b", 
     ylim = c(0, max(lipid_upper) * 1.1), pch = 19, col = "blue", lwd = 2,
     xlab = "Component", ylab = "Error Rate", 
     main = "LIPID Error Rate with Confidence Intervals")
arrows(1:ncomp_to_use, lipid_lower, 1:ncomp_to_use, lipid_upper, 
       angle = 90, code = 3, length = 0.1, col = "blue")
grid()

dev.off()


final_metab_acc <- 1 - perf.res$error.rate$METAB[ncomp_to_use, 1]
final_lipid_acc <- 1 - perf.res$error.rate$LIPID[ncomp_to_use, 1]
final_weighted_acc <- 1 - weighted_ber[ncomp_to_use, 1]



for (strategy in names(perf.res$choice.ncomp)) {
  cat(sprintf("  %s: 组件 %d\n", strategy, perf.res$choice.ncomp[[strategy]][1, 1]))
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
  })
  dev.off()
} else {
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
    return(NULL)
  }
  
  if (!is.null(fit$keepX[[block]])) {
    if (length(fit$keepX[[block]]) < comp || fit$keepX[[block]][comp] == 0) {
      return(NULL)
    }
  }
  
  pdf(file.path(output_dir, filename), width = 8, height = 6)
  tryCatch({
    plotLoadings(fit, comp = comp, block = block, method = "mean", contrib = "max")
  }, error = function(e) {
  })
  dev.off()
}

plot_loadings(diablo.fit, "METAB", 1, "loadings_METAB_comp1.pdf")
if(adjusted.ncomp >= 2) plot_loadings(diablo.fit, "METAB", 2, "loadings_METAB_comp2.pdf")
plot_loadings(diablo.fit, "LIPID", 1, "loadings_LIPID_comp1.pdf")
if(adjusted.ncomp >= 2) plot_loadings(diablo.fit, "LIPID", 2, "loadings_LIPID_comp2.pdf")

pdf(file.path(output_dir, "plotDiablo.pdf"), width = 8, height = 6)
plotDiablo(diablo.fit, ncomp = 1)
dev.off()

pdf(file.path(output_dir, "plotArrow.pdf"), width = 8, height = 6)
plotArrow(diablo.fit, ind.names = FALSE, legend = TRUE)
dev.off()


pdf(file.path(output_dir, "heatmap_cimDiablo_com1.pdf"), width = 12, height = 10)
cimDiablo(diablo.fit, comp = 1, 
          title = "DIABLO Heatmap - Component 1",
          legend.position = "right")
dev.off()

if(adjusted.ncomp >= 2) {
  pdf(file.path(output_dir, "heatmap_cimDiablo_com2.pdf"), width = 12, height = 10)
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
if(adjusted.ncomp >= 2) save_loadings_table(diablo.fit, "METAB", 2, file.path(output_dir, "loadings_table_METAB_comp2.csv"))
save_loadings_table(diablo.fit, "LIPID", 1, file.path(output_dir, "loadings_table_LIPID_comp1.csv"))
if(adjusted.ncomp >= 2) save_loadings_table(diablo.fit, "LIPID", 2, file.path(output_dir, "loadings_table_LIPID_comp2.csv"))


# =======================================
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
    adjusted_keepX = keepX_adjusted,
    design_matrix = design
  )
)

saveRDS(final_results, file.path(output_dir, "diablo_complete_results.rds"))

cat("\n===================\n")

