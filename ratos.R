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

data_path <- file.path(utils::shortPathName(getwd()), "Data_Cortex_Nuclear.xls")
if (!file.exists(data_path)) data_path <- file.choose()
data_path <- utils::shortPathName(data_path)
read_path <- file.path(tempdir(), "Data_Cortex_Nuclear_for_analysis.xls")
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
for (i in seq_len(nrow(mouse_metadata))) {
  idx <- dados$MouseID == mouse_metadata$MouseID[i]
  mouse_means[i, ] <- vapply(dados[idx, protein_names, drop = FALSE], function(z) {
    if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
  }, numeric(1))
}

class_counts <- as.data.frame(table(mouse_metadata$class), stringsAsFactors = FALSE)
names(class_counts) <- c("class", "number_of_mice")
class_measurements <- as.data.frame(table(dados$class), stringsAsFactors = FALSE)
names(class_measurements) <- c("class", "number_of_measurements")
class_counts <- merge(class_counts, class_measurements, by = "class", sort = FALSE)

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
  missing_measurements = vapply(mouse_metadata$MouseID, function(id) {
    sum(is.na(as.matrix(dados[dados$MouseID == id, protein_names, drop = FALSE])))
  }, integer(1)),
  proteins_missing_after_averaging = rowSums(is.na(mouse_means)),
  stringsAsFactors = FALSE
)

total_protein_values <- nrow(dados) * length(protein_names)
total_missing <- sum(is.na(as.matrix(dados[protein_names])))
class_colors <- setNames(grDevices::hcl.colors(nlevels(droplevels(mouse_metadata$class)), "Dark 3"),
                         levels(droplevels(mouse_metadata$class)))

cat("\nMICE PROTEIN EXPRESSION - PRELIMINARY ANALYSIS\n")
cat("Source file: Data_Cortex_Nuclear.xls\n")
cat("Measurement rows:", nrow(dados), "| Mice:", nrow(mouse_metadata),
    "| Proteins:", length(protein_names), "| Classes:", nlevels(droplevels(dados$class)), "\n")
cat("Measurements per mouse (min-max):", min(measurements_per_mouse), "-",
    max(measurements_per_mouse), "\n\n")
cat("Mouse and measurement counts by class:\n")
print(class_counts, row.names = FALSE)
cat("\nMissing protein readings:", total_missing, "of", total_protein_values,
    sprintf("(%.2f%%)\n", 100 * total_missing / total_protein_values))
cat("Most incomplete proteins:\n")
print(utils::head(missing_by_protein, 12), row.names = FALSE)
cat("\nMouse-level missingness summary:\n")
print(summary(missing_by_mouse$missing_measurements))

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

# Figure 1. Number of independent mice in each of the eight genotype-treatment-behavior classes.
counts <- table(droplevels(mouse_metadata$class))
par(mar = c(8, 5, 4, 2))
bp <- barplot(counts, col = class_colors[names(counts)], border = NA, las = 2,
              ylab = "Number of mice", main = "Mouse-level class balance")
text(bp, counts, labels = counts, pos = 3, cex = 0.9)

# Figure 2. Proteins with the greatest amount of missing measurement data.
n_show <- min(15, nrow(missing_by_protein))
miss_plot <- missing_by_protein[seq_len(n_show), ]
par(mar = c(5, 10, 4, 2))
bp <- barplot(rev(miss_plot$missing_percent), names.arg = rev(miss_plot$protein), horiz = TRUE,
              las = 1, col = "#4477AA", border = NA, xlab = "Missing measurements (%)",
              xlim = c(0, max(miss_plot$missing_percent) * 1.18),
              main = "Missingness is concentrated in a small set of proteins")
text(rev(miss_plot$missing_percent), bp,
     labels = sprintf("%.1f%%", rev(miss_plot$missing_percent)), pos = 4, cex = 0.75)

# Figure 3. Missing mouse-level protein averages across mice for the most incomplete proteins.
top_missing <- head(missing_by_protein$protein[missing_by_protein$missing_measurements > 0], 12)
if (length(top_missing)) {
  missing_pattern <- is.na(mouse_means[, top_missing, drop = FALSE]) * 1
  par(mar = c(8, 9, 4, 2))
  image(x = seq_len(ncol(missing_pattern)), y = seq_len(nrow(missing_pattern)),
        z = t(missing_pattern[nrow(missing_pattern):1, , drop = FALSE]),
        col = c("white", "#CC6677"), axes = FALSE, xlab = "Protein", ylab = "Mice",
        main = "Missing mouse-level protein averages (red = missing)")
  axis(1, at = seq_len(ncol(missing_pattern)), labels = colnames(missing_pattern), las = 2, cex.axis = 0.75)
  axis(2, at = seq_len(nrow(missing_pattern)), labels = rev(rownames(missing_pattern)), las = 2, cex.axis = 0.5)
}

# Figure 4. Mouse-level distributions of four high-variance proteins across the eight classes.
if (length(selected_proteins)) {
  oldpar <- par(mfrow = c(2, 2), mar = c(8, 4, 4, 1))
  for (p in selected_proteins) {
    boxplot(mouse_means[, p] ~ class_factor, col = class_colors[levels(class_factor)],
            las = 2, xlab = "Class", ylab = "Mouse-level mean protein expression",
            main = p, outline = TRUE)
  }
  par(oldpar)
}

if (nrow(X) >= 3 && ncol(X) >= 2) {
  pca_fit <- stats::prcomp(X, center = TRUE, scale. = TRUE)
  pca_variance <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  pca_scores <- data.frame(MouseID = mouse_metadata$MouseID, class = class_factor,
                           PC1 = pca_fit$x[, 1], PC2 = pca_fit$x[, 2])
  cat(sprintf("\nPCA after median imputation and standardization: PC1 %.1f%%; PC2 %.1f%%; combined %.1f%%.\n",
              100 * pca_variance[1], 100 * pca_variance[2], 100 * sum(pca_variance[1:2])))
  cat("PCA imputation and scaling use all mice for description only; refit them within training folds for prediction.\n")

  # Figure 5. Mouse-level PCA scores show class overlap in the first two components.
  par(mar = c(5, 5, 4, 10), xpd = TRUE)
  plot(pca_scores$PC1, pca_scores$PC2, pch = 19,
       col = class_colors[as.character(pca_scores$class)],
       xlab = sprintf("PC1 (%.1f%%)", 100 * pca_variance[1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * pca_variance[2]),
       main = "PCA of mouse-level protein averages")
  legend("topright", inset = c(-0.32, 0), legend = names(class_colors),
         col = class_colors, pch = 19, bty = "n", cex = 0.85)

  # Figure 6. Proportion of total variance explained by successive principal components.
  par(mar = c(5, 5, 4, 2), xpd = FALSE)
  n_scree <- min(20, length(pca_variance))
  plot(seq_len(n_scree), pca_variance[seq_len(n_scree)], type = "b", pch = 19,
       xlab = "Principal component", ylab = "Proportion of variance explained",
       main = "PCA scree plot")

  # Figure 7. Correlations among the most variable proteins at mouse level.
  corr_proteins <- names(sort(apply(X, 2, stats::var), decreasing = TRUE))[seq_len(min(20, ncol(X)))]
  corr <- stats::cor(X[, corr_proteins, drop = FALSE])
  cat("\nProteins included in the correlation heatmap:\n")
  print(corr_proteins)
  par(mar = c(10, 10, 4, 2))
  pal <- grDevices::colorRampPalette(c("#4477AA", "white", "#CC6677"))(101)
  image(seq_len(ncol(corr)), seq_len(nrow(corr)), t(corr[nrow(corr):1, , drop = FALSE]),
        col = pal, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "",
        main = "Correlation heatmap: most variable proteins")
  axis(1, at = seq_len(ncol(corr)), labels = colnames(corr), las = 2, cex.axis = 0.7)
  axis(2, at = seq_len(nrow(corr)), labels = rev(rownames(corr)), las = 2, cex.axis = 0.7)
}

cat("\nQuestions for the group:\n")
cat("1. Keep the eight genotype-treatment-behavior classes as the primary target, with genotype-only as a secondary comparison?\n")
cat("2. Use mouse-level means as the primary analysis unit and keep each mouse within a single validation fold?\n")
cat("3. Compare removing the six most incomplete proteins against imputation fitted within each training fold?\n")
cat("\nNext: compare supervised classifiers on original and PCA-reduced data using mouse-grouped validation.\n")
cat("With 77 proteins and only 72 mice, standard LDA/QDA covariance estimates may be singular; reduce dimension or use regularization.\n")

if (interactive()) {
  View(missing_by_protein, title = "Protein missingness")
}
