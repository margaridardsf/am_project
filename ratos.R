options(stringsAsFactors = FALSE, scipen = 999)
if (!interactive()) grDevices::pdf(NULL)

if (!requireNamespace("readxl", quietly = TRUE)) {
  user_library <- file.path(
    Sys.getenv("LOCALAPPDATA", unset = path.expand("~")),
    "R", "win-library",
    paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
  )
  if (!dir.exists(user_library)) dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_library, .libPaths()))
  install.packages("readxl", lib = user_library, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. Install it with install.packages('readxl').")
}

data_path <- "Data_Cortex_Nuclear.xls"
if (!file.exists(data_path)) data_path <- file.choose()
read_path <- file.path(tempdir(), "Data_Cortex_Nuclear.xls")
if (!file.copy(data_path, read_path, overwrite = TRUE)) {
  stop("Could not copy the source workbook to a temporary path for reading.")
}
dados <- as.data.frame(readxl::read_excel(read_path, na = c("NA", "")))
names(dados) <- trimws(names(dados))

required_metadata <- c("MouseID", "Genotype", "Treatment", "Behavior", "class")
missing_metadata <- setdiff(required_metadata, names(dados))
if (length(missing_metadata)) {
  stop("Expected metadata columns not found: ", paste(missing_metadata, collapse = ", "))
}

dados$MeasurementID <- as.character(dados$MouseID)
dados$MouseID <- sub("_[0-9]+$", "", dados$MeasurementID)
for (v in c("Genotype", "Treatment", "Behavior", "class")) dados[[v]] <- factor(dados[[v]])
protein_names <- setdiff(names(dados), c(required_metadata, "MeasurementID"))
if (!length(protein_names)) stop("No protein columns were detected.")
non_numeric <- protein_names[!vapply(dados[protein_names], is.numeric, logical(1))]
if (length(non_numeric)) {
  stop("Protein columns must be numeric: ", paste(non_numeric, collapse = ", "))
}
protein_pairs <- combn(protein_names, 2, simplify = FALSE)
exact_duplicate_list <- Filter(function(pair) identical(dados[[pair[1]]], dados[[pair[2]]]), protein_pairs)
exact_duplicate_features <- if (length(exact_duplicate_list)) {
  as.data.frame(do.call(rbind, exact_duplicate_list), stringsAsFactors = FALSE)
} else {
  data.frame(protein_1 = character(0), protein_2 = character(0))
}
names(exact_duplicate_features) <- c("protein_1", "protein_2")

measurements_per_mouse <- table(dados$MouseID)
if (any(measurements_per_mouse != 15)) {
  warning("Not every mouse has exactly 15 measurements. Check the MouseID suffix parsing.")
}

mouse_metadata <- dados[!duplicated(dados$MouseID), required_metadata]
for (v in required_metadata[-1]) {
  n_unique <- tapply(as.character(dados[[v]]), dados$MouseID,
                     function(z) length(unique(na.omit(z))))
  if (any(n_unique > 1)) stop("Metadata variable changes within MouseID: ", v)
}
row.names(mouse_metadata) <- mouse_metadata$MouseID

mouse_means <- matrix(NA_real_, nrow = nrow(mouse_metadata), ncol = length(protein_names),
                      dimnames = list(mouse_metadata$MouseID, protein_names))
missing_counts <- matrix(0L, nrow = nrow(mouse_metadata), ncol = length(protein_names),
                         dimnames = list(mouse_metadata$MouseID, protein_names))
for (i in seq_len(nrow(mouse_metadata))) {
  idx <- dados$MouseID == mouse_metadata$MouseID[i]
  mouse_means[i, ] <- vapply(dados[idx, protein_names, drop = FALSE], function(z) {
    if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
  }, numeric(1))
  missing_counts[i, ] <- colSums(is.na(dados[idx, protein_names, drop = FALSE]))
}

class_levels <- levels(droplevels(mouse_metadata$class))
mice_by_class <- table(factor(mouse_metadata$class, levels = class_levels))
measurements_by_class <- table(factor(dados$class, levels = class_levels))
class_counts <- rbind(Mice = as.integer(mice_by_class),
                      Measurement_rows = as.integer(measurements_by_class))
colnames(class_counts) <- class_levels

factor_balance <- do.call(rbind, lapply(c("Genotype", "Treatment", "Behavior"), function(v) {
  counts <- table(mouse_metadata[[v]])
  data.frame(factor = v, level = names(counts), mice = as.integer(counts), row.names = NULL)
}))

missing_by_protein <- data.frame(
  protein = protein_names,
  missing_measurements = vapply(dados[protein_names], function(z) sum(is.na(z)), integer(1)),
  total_measurements = nrow(dados),
  missing_percent = 100 * vapply(dados[protein_names], function(z) mean(is.na(z)), numeric(1)),
  mice_missing_all_15 = vapply(protein_names, function(p) {
    sum(vapply(split(dados[[p]], dados$MouseID), function(z) all(is.na(z)), logical(1)))
  }, integer(1)),
  stringsAsFactors = FALSE
)
missing_by_protein <- missing_by_protein[order(-missing_by_protein$missing_measurements,
                                               missing_by_protein$protein), ]

missing_by_mouse <- data.frame(
  MouseID = mouse_metadata$MouseID,
  class = mouse_metadata$class,
  missing_measurements = rowSums(missing_counts),
  missing_protein_averages = rowSums(missing_counts == 15),
  partly_missing_proteins = rowSums(missing_counts > 0 & missing_counts < 15),
  proteins_missing_after_averaging = rowSums(is.na(mouse_means)),
  stringsAsFactors = FALSE
)

missingness_by_class <- do.call(rbind, lapply(class_levels, function(cl) {
  in_class <- mouse_metadata$class == cl
  missing_averages <- missing_by_mouse$proteins_missing_after_averaging[in_class]
  data.frame(
    class = cl,
    mice = sum(in_class),
    mice_with_missing_averages = sum(missing_averages > 0),
    percent_mice_with_missing_averages = 100 * mean(missing_averages > 0),
    median_missing_protein_averages = stats::median(missing_averages),
    max_missing_protein_averages = max(missing_averages),
    stringsAsFactors = FALSE
  )
}))

total_protein_values <- nrow(dados) * length(protein_names)
total_missing <- sum(is.na(as.matrix(dados[protein_names])))
proteins_with_missing <- sum(colSums(missing_counts) > 0)
measurement_rows_with_missing <- sum(rowSums(is.na(as.matrix(dados[protein_names]))) > 0)
fully_missing_mouse_protein_pairs <- sum(missing_counts == 15)
partly_missing_mouse_protein_pairs <- sum(missing_counts > 0 & missing_counts < 15)
within_sd_by_protein <- vapply(protein_names, function(p) {
  mouse_sds <- vapply(split(dados[[p]], dados$MouseID), function(z) {
    if (sum(!is.na(z)) < 2) NA_real_ else sd(z, na.rm = TRUE)
  }, numeric(1))
  if (all(is.na(mouse_sds))) NA_real_ else mean(mouse_sds, na.rm = TRUE)
}, numeric(1))
between_sd_by_protein <- apply(mouse_means, 2, sd, na.rm = TRUE)
within_between_ratio <- within_sd_by_protein / between_sd_by_protein
variability_by_protein <- data.frame(
  protein = protein_names,
  mean_within_mouse_sd = within_sd_by_protein,
  between_mouse_sd_of_averages = between_sd_by_protein,
  within_to_between_sd_ratio = within_between_ratio,
  stringsAsFactors = FALSE
)
variability_by_protein <- variability_by_protein[order(-variability_by_protein$within_to_between_sd_ratio), ]
eta_squared <- function(x, group) {
  keep <- !is.na(x) & !is.na(group)
  x <- x[keep]
  group <- droplevels(factor(group[keep]))
  if (length(x) < 2 || nlevels(group) < 2) return(NA_real_)
  total_ss <- sum((x - mean(x))^2)
  if (total_ss == 0) return(NA_real_)
  group_ss <- sum(vapply(split(x, group), function(z) {
    length(z) * (mean(z) - mean(x))^2
  }, numeric(1)))
  group_ss / total_ss
}
factor_effects <- data.frame(
  protein = protein_names,
  genotype_eta2 = vapply(protein_names, function(p) eta_squared(mouse_means[, p], mouse_metadata$Genotype), numeric(1)),
  treatment_eta2 = vapply(protein_names, function(p) eta_squared(mouse_means[, p], mouse_metadata$Treatment), numeric(1)),
  behavior_eta2 = vapply(protein_names, function(p) eta_squared(mouse_means[, p], mouse_metadata$Behavior), numeric(1)),
  eight_class_eta2 = vapply(protein_names, function(p) eta_squared(mouse_means[, p], mouse_metadata$class), numeric(1)),
  stringsAsFactors = FALSE
)
class_colors <- setNames(grDevices::hcl.colors(nlevels(droplevels(mouse_metadata$class)), "Dark 3"),
                         levels(droplevels(mouse_metadata$class)))

cat("\nMICE PROTEIN EXPRESSION - PRELIMINARY ANALYSIS\n")
cat("Source file: Data_Cortex_Nuclear.xls\n")
cat("Measurement rows:", nrow(dados), "| Mice:", nrow(mouse_metadata),
    "| Proteins:", length(protein_names), "| Classes:", nlevels(droplevels(dados$class)), "\n")
cat("Measurements per mouse (min-max):", min(measurements_per_mouse), "-",
    max(measurements_per_mouse), "\n\n")
cat("Class counts (two rows: mice and measurement rows):\n")
print(class_counts)
cat("\nMouse counts by experimental factor:\n")
print(factor_balance, row.names = FALSE)
cat("\nMissing protein readings:", total_missing, "of", total_protein_values,
    sprintf("(%.2f%%)\n", 100 * total_missing / total_protein_values))
cat("Proteins with missing values:", proteins_with_missing, "of", length(protein_names), "\n")
cat("Measurement rows with at least one missing protein:", measurement_rows_with_missing,
    "of", nrow(dados), "| complete rows:", nrow(dados) - measurement_rows_with_missing, "\n")
cat("Mouse-protein pairs fully missing (15/15):", fully_missing_mouse_protein_pairs,
    "| partly missing (1-14/15):", partly_missing_mouse_protein_pairs, "\n")
cat("Mouse-level protein averages still missing after averaging:", sum(is.na(mouse_means)), "\n")
cat("Most incomplete proteins (top six):\n")
top_missing_console <- utils::head(missing_by_protein[, c("protein", "missing_measurements",
                                                           "missing_percent", "mice_missing_all_15")], 6)
top_missing_console$missing_percent <- round(top_missing_console$missing_percent, 1)
print(top_missing_console, row.names = FALSE)
cat("Exact duplicate protein columns:", nrow(exact_duplicate_features), "\n")
if (nrow(exact_duplicate_features)) print(exact_duplicate_features, row.names = FALSE)
cat("Mice with at least one protein average missing:", sum(rowSums(is.na(mouse_means)) > 0),
    "of", nrow(mouse_means), "| complete mouse profiles:",
    sum(rowSums(is.na(mouse_means)) == 0), "\n")
cat("Mice with the most missing readings (raw 15-measurement data):\n")
print(utils::head(missing_by_mouse[order(-missing_by_mouse$missing_measurements),
                                   c("MouseID", "class", "missing_measurements",
                                     "missing_protein_averages", "partly_missing_proteins")], 5),
      row.names = FALSE)
cat("\nMissing mouse-level protein averages by class (descriptive; not a test of missingness mechanism):\n")
missingness_by_class_print <- data.frame(
  Class = missingness_by_class$class,
  Mice = missingness_by_class$mice,
  Mice_with_NA = missingness_by_class$mice_with_missing_averages,
  Pct_mice_NA = round(missingness_by_class$percent_mice_with_missing_averages, 1)
)
print(missingness_by_class_print, row.names = FALSE)
cat("Median within-mouse SD / between-mouse SD ratio:",
    sprintf("%.2f", median(within_between_ratio, na.rm = TRUE)),
    "(larger values mean repeated readings vary more within mice)\n")
cat("Proteins with highest within/between variability ratio (full table in Viewer):\n")
print(utils::head(variability_by_protein[, c("protein", "within_to_between_sd_ratio")], 5),
      row.names = FALSE)
cat("\nMedian marginal eta-squared across proteins by factor (descriptive, not causal):\n")
eta_summary <- data.frame(
  factor = c("Genotype", "Treatment", "Behavior", "Eight-class label"),
  median_eta2 = c(median(factor_effects$genotype_eta2, na.rm = TRUE),
                  median(factor_effects$treatment_eta2, na.rm = TRUE),
                  median(factor_effects$behavior_eta2, na.rm = TRUE),
                  median(factor_effects$eight_class_eta2, na.rm = TRUE))
)
eta_summary$median_eta2 <- round(eta_summary$median_eta2, 3)
print(eta_summary, row.names = FALSE)
cat("These marginal effect sizes describe association with each factor; they are not classifier accuracy or adjusted causal effects.\n")

cat("\nSlide-ready summary: 72 mice, 77 proteins, 8 classes, 15 measurements per mouse; 1.68% of readings are missing.\n")

class_factor <- droplevels(mouse_metadata$class)
observed_medians <- apply(mouse_means, 2, median, na.rm = TRUE)
usable_proteins <- names(observed_medians)[is.finite(observed_medians)]
X <- mouse_means[, usable_proteins, drop = FALSE]
for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- observed_medians[colnames(X)[j]]
variable_proteins <- apply(X, 2, stats::sd) > 0
X <- X[, variable_proteins, drop = FALSE]

protein_variances <- apply(mouse_means, 2, stats::var, na.rm = TRUE)
protein_variances <- protein_variances[is.finite(protein_variances)]
selected_proteins <- names(sort(protein_variances, decreasing = TRUE))[seq_len(min(4, length(protein_variances)))]
cat("\nHighest-variance proteins shown in the boxplots:\n")
print(selected_proteins)

# Figure 1. Proteins with the greatest amount of missing measurement data.
n_show <- min(15, nrow(missing_by_protein))
miss_plot <- missing_by_protein[seq_len(n_show), ]
par(mar = c(5, 10, 4, 2))
bp <- barplot(rev(miss_plot$missing_percent), names.arg = rev(miss_plot$protein), horiz = TRUE,
              las = 1, col = "#4477AA", border = NA, xlab = "Missing measurements (%)",
              xlim = c(0, max(miss_plot$missing_percent) * 1.18),
              main = "Missingness is concentrated in a small set of proteins")
text(rev(miss_plot$missing_percent), bp,
     labels = sprintf("%.1f%%", rev(miss_plot$missing_percent)), pos = 4, cex = 0.75)

# Figure 2. Missing reading counts for the most incomplete proteins; each cell label is out of 15.
top_missing <- head(missing_by_protein$protein[missing_by_protein$missing_measurements > 0], 10)
if (length(top_missing)) {
  missing_plot <- missing_counts[, top_missing, drop = FALSE]
  mice_to_show <- rowSums(missing_plot) > 0
  mice_selected <- mouse_metadata[mice_to_show, , drop = FALSE]
  display_order <- order(mice_selected$class, mice_selected$MouseID)
  mice_ordered <- mice_selected[display_order, , drop = FALSE]
  missing_plot <- missing_plot[mice_to_show, , drop = FALSE][display_order, , drop = FALSE]
  missing_plot_reversed <- missing_plot[nrow(missing_plot):1, , drop = FALSE]
  missing_colors <- grDevices::colorRampPalette(c("#FFF7F7", "#E58B9A", "#8F1230"))(16)
  par(mar = c(9, 8, 4, 2))
  image(x = seq_len(ncol(missing_plot_reversed)), y = seq_len(nrow(missing_plot_reversed)),
        z = t(missing_plot_reversed), col = missing_colors,
        breaks = seq(-0.5, 15.5, by = 1), axes = FALSE,
        xlab = "Protein", ylab = "Mouse (grouped by class)",
        main = "Missing readings per mouse and protein (each label is out of 15)")
  mtext("Darker red indicates more missing readings", side = 3, line = 0.3, cex = 0.75)
  axis(1, at = seq_len(ncol(missing_plot_reversed)), labels = colnames(missing_plot_reversed),
       las = 2, cex.axis = 0.75)
  axis(2, at = seq_len(nrow(missing_plot_reversed)), labels = rev(rownames(missing_plot_reversed)),
       las = 2, cex.axis = 0.6)
  missing_cells <- which(missing_plot_reversed > 0, arr.ind = TRUE)
  if (nrow(missing_cells)) {
    cell_values <- missing_plot_reversed[missing_cells]
    text(x = missing_cells[, 2], y = missing_cells[, 1], labels = cell_values,
         cex = 0.58, col = ifelse(cell_values >= 8, "white", "#4A0913"))
  }
  class_breaks <- which(as.character(mice_ordered$class[-1]) !=
                          as.character(mice_ordered$class[-nrow(mice_ordered)]))
  if (length(class_breaks)) {
    abline(h = nrow(missing_plot_reversed) - class_breaks + 0.5, col = "#555555", lwd = 0.8)
  }
}

# Figure 3. Mouse-level distributions of four high-variance proteins across the eight classes.
if (length(selected_proteins)) {
  oldpar <- par(mfrow = c(2, 2), mar = c(8, 4, 4, 1))
  for (p in selected_proteins) {
    boxplot(mouse_means[, p] ~ class_factor, col = class_colors[levels(class_factor)],
            las = 2, xlab = "Class", ylab = "Mouse-level mean protein expression",
            main = p, outline = TRUE)
  }
  par(oldpar)
}

# Figure 4. Compare the spread of univariate effect sizes across experimental factors.
par(mar = c(6, 5, 4, 2))
boxplot(list(Genotype = factor_effects$genotype_eta2,
             Treatment = factor_effects$treatment_eta2,
             Behavior = factor_effects$behavior_eta2,
             Eight_classes = factor_effects$eight_class_eta2),
        col = c("#4477AA", "#66CCEE", "#228833", "#CC6677"),
        ylab = "Marginal eta-squared across proteins", xlab = "Factor",
        main = "Exploratory per-protein association with each factor",
        las = 2, outline = TRUE)

if (nrow(X) >= 3 && ncol(X) >= 2) {
  pca_fit <- stats::prcomp(X, center = TRUE, scale. = TRUE)
  pca_variance <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  pca_cumulative <- cumsum(pca_variance)
  n_pc_80 <- which(pca_cumulative >= 0.80)[1]
  n_pc_90 <- which(pca_cumulative >= 0.90)[1]
  pca_loadings <- data.frame(protein = rownames(pca_fit$rotation), pca_fit$rotation,
                             check.names = FALSE, row.names = NULL)
  pca_scores <- data.frame(MouseID = mouse_metadata$MouseID, class = class_factor,
                           Genotype = mouse_metadata$Genotype,
                           Treatment = mouse_metadata$Treatment,
                           Behavior = mouse_metadata$Behavior,
                           PC1 = pca_fit$x[, 1], PC2 = pca_fit$x[, 2])
  cat(sprintf("\nPCA after median imputation and standardization: PC1 %.1f%%; PC2 %.1f%%; combined %.1f%%.\n",
              100 * pca_variance[1], 100 * pca_variance[2], 100 * sum(pca_variance[1:2])))
  cat("Components needed for 80% / 90% cumulative variance:", n_pc_80, "/", n_pc_90, "\n")
  cat("Largest absolute loadings for PC1 and PC2 (for component interpretation):\n")
  for (component in c("PC1", "PC2")) {
    loading_order <- order(abs(pca_loadings[[component]]), decreasing = TRUE)[1:5]
    cat(component, "\n")
    print(data.frame(protein = pca_loadings$protein[loading_order],
                     loading = round(pca_loadings[[component]][loading_order], 3)), row.names = FALSE)
  }
  cat("PCA imputation and scaling use all mice for description only; refit them within training folds for prediction.\n")

  # Figure 4. Mouse-level PCA scores show class overlap in the first two components.
  par(mar = c(5, 5, 4, 10), xpd = TRUE)
  plot(pca_scores$PC1, pca_scores$PC2, pch = 19,
       col = class_colors[as.character(pca_scores$class)],
       xlab = sprintf("PC1 (%.1f%%)", 100 * pca_variance[1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * pca_variance[2]),
       main = "PCA of mouse-level protein averages")
  legend("topright", inset = c(-0.32, 0), legend = names(class_colors),
         col = class_colors, pch = 19, bty = "n", cex = 0.85)

  # Figure 5. The same PCA scores colored by genotype compare the binary and eight-class questions.
  genotype_levels <- levels(droplevels(mouse_metadata$Genotype))
  genotype_colors <- setNames(grDevices::hcl.colors(length(genotype_levels), "Set 2"), genotype_levels)
  par(mar = c(5, 5, 4, 8), xpd = TRUE)
  plot(pca_scores$PC1, pca_scores$PC2, pch = 19,
       col = genotype_colors[as.character(pca_scores$Genotype)],
       xlab = sprintf("PC1 (%.1f%%)", 100 * pca_variance[1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * pca_variance[2]),
       main = "PCA scores by genotype")
  legend("topright", inset = c(-0.25, 0), legend = names(genotype_colors),
         col = genotype_colors, pch = 19, bty = "n", cex = 0.9)

  # Figure 6. Proportion of total variance explained by successive principal components.
  par(mar = c(5, 5, 4, 2), xpd = FALSE)
  n_scree <- min(20, length(pca_variance))
  plot(seq_len(n_scree), pca_variance[seq_len(n_scree)], type = "b", pch = 19,
       xlab = "Principal component", ylab = "Proportion of variance explained",
       main = "PCA scree plot")

  # Figure 7. Correlations among the most variable proteins at mouse level.
  corr_proteins <- names(sort(apply(X, 2, stats::var), decreasing = TRUE))[seq_len(min(20, ncol(X)))]
  corr <- stats::cor(X[, corr_proteins, drop = FALSE])
  all_correlations <- stats::cor(X)
  correlation_indices <- which(upper.tri(all_correlations), arr.ind = TRUE)
  correlation_pairs <- data.frame(
    protein_1 = rownames(all_correlations)[correlation_indices[, 1]],
    protein_2 = colnames(all_correlations)[correlation_indices[, 2]],
    correlation = all_correlations[correlation_indices],
    stringsAsFactors = FALSE
  )
  correlation_pairs <- correlation_pairs[order(-abs(correlation_pairs$correlation)), ]
  cat("\nStrongest absolute protein correlations (median-imputed mouse averages; exploratory):\n")
  correlation_pairs$correlation <- round(correlation_pairs$correlation, 3)
  print(utils::head(correlation_pairs, 6), row.names = FALSE)
  cat("\nProteins included in the correlation heatmap:\n")
  print(corr_proteins)
  par(mar = c(10, 10, 4, 2))
  pal <- grDevices::colorRampPalette(c("#4477AA", "white", "#CC6677"))(101)
  image(seq_len(ncol(corr)), seq_len(nrow(corr)), t(corr[nrow(corr):1, , drop = FALSE]),
        col = pal, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "",
        main = "Correlation heatmap: most variable proteins")
  axis(1, at = seq_len(ncol(corr)), labels = colnames(corr), las = 2, cex.axis = 0.7)
  axis(2, at = seq_len(nrow(corr)), labels = rev(rownames(corr)), las = 2, cex.axis = 0.7)

  # Figure 8. Compare average within-mouse spread with between-mouse spread for each protein.
  valid_variability <- is.finite(variability_by_protein$mean_within_mouse_sd) &
    is.finite(variability_by_protein$between_mouse_sd_of_averages) &
    variability_by_protein$mean_within_mouse_sd > 0 &
    variability_by_protein$between_mouse_sd_of_averages > 0
  variability_plot <- variability_by_protein[valid_variability, ]
  par(mar = c(5, 5, 4, 2), xpd = FALSE)
  plot(variability_plot$between_mouse_sd_of_averages, variability_plot$mean_within_mouse_sd,
       log = "xy", pch = 19, col = "#4477AA",
       xlab = "Between-mouse SD of mouse averages",
       ylab = "Mean within-mouse SD across 15 readings",
       main = "Within- versus between-mouse variability")
  abline(0, 1, lty = 2, col = "#CC6677")
  label_rows <- head(order(variability_plot$within_to_between_sd_ratio, decreasing = TRUE), 5)
  text(variability_plot$between_mouse_sd_of_averages[label_rows],
       variability_plot$mean_within_mouse_sd[label_rows],
       labels = variability_plot$protein[label_rows], pos = 4, cex = 0.65)
}

cat("\nQuestions for the group:\n")
cat("1. Keep the eight genotype-treatment-behavior classes as the primary target, with genotype-only as a secondary comparison?\n")
cat("2. Use mouse-level means as the primary analysis unit and keep each mouse within a single validation fold?\n")
cat("3. Compare removing the six most incomplete proteins against imputation fitted within each training fold?\n")
cat("\nNext: compare supervised classifiers on original and PCA-reduced data using mouse-grouped validation.\n")
cat("With 77 proteins and only 72 mice, standard LDA/QDA covariance estimates may be singular; reduce dimension or use regularization.\n")
cat("\nFive-slide priority: dataset and mouse counts; research questions; missing-count map plus PCA; validation/missingness/dimension challenges; workplan and discussion questions.\n")
cat("For the five-minute talk, prioritize the missing-count map and eight-class PCA; use the boxplots or factor-effect plot only if time permits.\n")

protein_diagnostics <- merge(missing_by_protein, variability_by_protein, by = "protein", all = TRUE, sort = FALSE)
protein_diagnostics <- merge(protein_diagnostics, factor_effects, by = "protein", all = TRUE, sort = FALSE)
if (exists("pca_loadings")) {
  protein_diagnostics <- merge(protein_diagnostics,
                               pca_loadings[, c("protein", "PC1", "PC2")],
                               by = "protein", all = TRUE, sort = FALSE)
}

if (interactive()) {
  View(protein_diagnostics, title = "Protein-level diagnostics")
  View(missing_by_mouse, title = "Mouse-level diagnostics")
  View(missingness_by_class, title = "Missingness by class")
}
