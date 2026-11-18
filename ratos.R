options(stringsAsFactors = FALSE, scipen = 999)
if (!interactive()) grDevices::pdf(NULL)

if (!file.exists("project_utils.R")) {
  stop("Run this script with the project folder as the working directory.")
}
source("project_utils.R")
am_project_packages(c("readxl", "ggplot2", "corrplot"))

plot_output_dir <- "plots_ratos"
if (!dir.exists(plot_output_dir)) dir.create(plot_output_dir, recursive = TRUE)
save_ggplot <- function(plot_object, filename, width = 10, height = 7) {
  print(plot_object)
  ggplot2::ggsave(file.path(plot_output_dir, filename), plot = plot_object,
                  width = width, height = height, units = "in", dpi = 220,
                  bg = "white")
}
save_base_plot <- function(filename, width = 10, height = 7, draw) {
  output_path <- file.path(plot_output_dir, filename)
  if (interactive()) {
    draw()
    recorded_plot <- grDevices::recordPlot()
    grDevices::png(output_path, width = width, height = height,
                   units = "in", res = 220, bg = "white")
    grDevices::replayPlot(recorded_plot)
    grDevices::dev.off()
  } else {
    grDevices::png(output_path, width = width, height = height,
                   units = "in", res = 220, bg = "white")
    draw()
    grDevices::dev.off()
  }
  invisible(output_path)
}

source_data <- am_read_mouse_data("Data_Cortex_Nuclear.xls")
dados <- source_data$data
protein_names <- source_data$protein_names
if (nrow(dados) != 1080L || length(unique(dados$MouseID)) != 72L ||
    length(protein_names) != 77L) {
  stop("The analysis expects the project workbook version with 1,080 rows, 72 mice and 77 proteins.")
}
measurement_index <- as.integer(sub("^.*_([0-9]+)$", "\\1", dados$MeasurementID))
if (anyNA(measurement_index)) stop("Measurement IDs must end with a numeric within-mouse suffix.")
protein_pairs <- combn(protein_names, 2, simplify = FALSE)
exact_duplicate_list <- Filter(function(pair) identical(dados[[pair[1]]], dados[[pair[2]]]), protein_pairs)
exact_duplicate_features <- if (length(exact_duplicate_list)) {
  as.data.frame(do.call(rbind, exact_duplicate_list), stringsAsFactors = FALSE)
} else {
  data.frame(protein_1 = character(0), protein_2 = character(0))
}
names(exact_duplicate_features) <- c("protein_1", "protein_2")

measurements_per_mouse <- table(dados$MouseID)
mouse_metadata <- am_mouse_metadata(dados)
mouse_means <- am_mouse_means(dados, protein_names, mouse_metadata$MouseID)
protein_matrix <- as.matrix(dados[, protein_names, drop = FALSE])
mouse_group <- factor(as.character(dados$MouseID), levels = mouse_metadata$MouseID)
missing_counts <- rowsum(1.0 * is.na(protein_matrix), mouse_group, reorder = FALSE)
missing_counts <- missing_counts[match(mouse_metadata$MouseID, rownames(missing_counts)), , drop = FALSE]
rownames(missing_counts) <- mouse_metadata$MouseID

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

protein_missing <- is.na(dados[protein_names])
readings_per_mouse <- as.integer(measurements_per_mouse[rownames(missing_counts)])
fully_missing_mouse_protein <- sweep(missing_counts, 1, readings_per_mouse, `==`)
partly_missing_mouse_protein <- missing_counts > 0 &
  sweep(missing_counts, 1, readings_per_mouse, `<`)
missing_by_protein <- data.frame(
  protein = protein_names,
  missing_measurements = colSums(protein_missing),
  total_measurements = nrow(dados),
  missing_percent = 100 * colMeans(protein_missing),
  mice_missing_all_readings = colSums(fully_missing_mouse_protein),
  stringsAsFactors = FALSE
)
missing_by_protein <- missing_by_protein[order(-missing_by_protein$missing_measurements,
                                               missing_by_protein$protein), ]
# Retain all proteins in the missingness summaries, but omit proteins that are
# entirely missing in at least five mice from the primary PCA.
high_missing_proteins <- missing_by_protein$protein[
  missing_by_protein$mice_missing_all_readings >= 5
]

missing_by_mouse <- data.frame(
  MouseID = mouse_metadata$MouseID,
  class = mouse_metadata$class,
  missing_measurements = rowSums(missing_counts),
  missing_protein_averages = rowSums(fully_missing_mouse_protein),
  partly_missing_proteins = rowSums(partly_missing_mouse_protein),
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
total_missing <- sum(protein_missing)
proteins_with_missing <- sum(colSums(missing_counts) > 0)
measurement_rows_with_missing <- sum(rowSums(protein_missing) > 0)
fully_missing_mouse_protein_pairs <- sum(fully_missing_mouse_protein)
partly_missing_mouse_protein_pairs <- sum(partly_missing_mouse_protein)
mouse_readings <- split(dados[protein_names], factor(dados$MouseID,
                                                     levels = unique(dados$MouseID)))
within_sd_by_protein <- vapply(protein_names, function(p) {
  mouse_sds <- vapply(mouse_readings, function(mouse_data) {
    z <- mouse_data[[p]]
    if (sum(!is.na(z)) < 2) NA_real_ else sd(z, na.rm = TRUE)
  }, numeric(1))
  if (all(is.na(mouse_sds))) NA_real_ else mean(mouse_sds, na.rm = TRUE)
}, numeric(1))
between_sd_by_protein <- apply(mouse_means, 2, sd, na.rm = TRUE)
variability_by_protein <- data.frame(
  protein = protein_names,
  mean_within_mouse_sd = within_sd_by_protein,
  between_mouse_sd_of_averages = between_sd_by_protein,
  stringsAsFactors = FALSE
)
variability_by_protein <- variability_by_protein[order(-variability_by_protein$mean_within_mouse_sd), ]

# Ahmed et al. (2015; doi:10.1371/journal.pone.0119491) describe a five-point
# dilution series printed in triplicate for each sample. Audit whether the numeric
# workbook suffix is associated with signal changes.
# This is consistent with an ordered dilution series, but cannot verify which
# suffixes map to the five dilution levels and three technical replicates.
measurement_rows_by_mouse <- split(seq_len(nrow(dados)), factor(
  dados$MouseID, levels = unique(dados$MouseID)
))
suffix_order_diagnostics <- do.call(rbind, lapply(protein_names, function(p) {
  mouse_rhos <- vapply(measurement_rows_by_mouse, function(index) {
    observed <- is.finite(dados[[p]][index]) & is.finite(measurement_index[index])
    if (sum(observed) < 5L) return(NA_real_)
    suppressWarnings(stats::cor(measurement_index[index][observed],
                                dados[[p]][index][observed], method = "spearman"))
  }, numeric(1))
  data.frame(
    protein = p,
    median_within_mouse_spearman = stats::median(mouse_rhos, na.rm = TRUE),
    mice_with_estimable_rho = sum(is.finite(mouse_rhos)),
    percent_mouse_profiles_decreasing = 100 * mean(mouse_rhos < 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

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
cat("Measurement structure: the original RPPA study describes 15 readings per sample as three technical replicates across a five-point dilution series. The workbook has numeric within-mouse suffixes but no explicit dilution or replicate columns. Arithmetic means and intervals therefore summarise assay readings across the series; they are not biological-replicate uncertainty intervals.\n")
cat("The suffix-order trend check below is exploratory only; it cannot establish the exact mapping from suffixes to dilution levels and technical replicates.\n")
cat("\nSuffix-order signal diagnostic (within-mouse Spearman correlations, descriptive):\n")
cat(sprintf("Median of protein-level median correlations: %.3f; proteins with a negative median: %d/%d; median percentage of mice with a negative profile correlation: %.1f%%.\n",
            stats::median(suffix_order_diagnostics$median_within_mouse_spearman, na.rm = TRUE),
            sum(suffix_order_diagnostics$median_within_mouse_spearman < 0, na.rm = TRUE),
            nrow(suffix_order_diagnostics),
            stats::median(suffix_order_diagnostics$percent_mouse_profiles_decreasing, na.rm = TRUE)))
cat("This positional pattern is compatible with a dilution series, but is not a substitute for dilution metadata or calibration.\n")
within_sd_lower_count <- sum(
  variability_by_protein$mean_within_mouse_sd < variability_by_protein$between_mouse_sd_of_averages,
  na.rm = TRUE
)
cat(sprintf("\nDescriptive spread check: mean within-mouse SD across the 15 assay readings is below the between-mouse SD of the 15-reading arithmetic means for %d/%d proteins. This is not a biological within-versus-between variance decomposition: the within-mouse spread includes dilution response and technical-replicate variation.\n",
            within_sd_lower_count, sum(is.finite(variability_by_protein$mean_within_mouse_sd) &
                                         is.finite(variability_by_protein$between_mouse_sd_of_averages))))
cat("Class counts (two rows: mice and measurement rows):\n")
print(class_counts)
cat("\nMouse counts by experimental factor:\n")
print(factor_balance, row.names = FALSE)
cat("\nMissing protein readings:", total_missing, "of", total_protein_values,
    sprintf("(%.2f%%)\n", 100 * total_missing / total_protein_values))
cat("Proteins with missing values:", proteins_with_missing, "of", length(protein_names), "\n")
cat("Measurement rows with at least one missing protein:", measurement_rows_with_missing,
    "of", nrow(dados), "| complete rows:", nrow(dados) - measurement_rows_with_missing, "\n")
cat("Mouse-protein pairs fully missing (all readings):", fully_missing_mouse_protein_pairs,
    "| partly missing (some readings):", partly_missing_mouse_protein_pairs, "\n")
cat("Mouse-level protein averages still missing after averaging:", sum(is.na(mouse_means)), "\n")
cat("Most incomplete proteins (top six):\n")
top_missing_console <- utils::head(missing_by_protein[, c("protein", "missing_measurements",
                                                          "missing_percent", "mice_missing_all_readings")], 6)
top_missing_console$missing_percent <- round(top_missing_console$missing_percent, 1)
print(top_missing_console, row.names = FALSE)
cat("Exact duplicate protein columns:", nrow(exact_duplicate_features), "\n")
if (nrow(exact_duplicate_features)) print(exact_duplicate_features, row.names = FALSE)
cat("Mice with at least one protein average missing:", sum(rowSums(is.na(mouse_means)) > 0),
    "of", nrow(mouse_means), "| complete mouse profiles:",
    sum(rowSums(is.na(mouse_means)) == 0), "\n")
cat("Mice with the most missing readings (raw repeated-measurement data):\n")
print(utils::head(missing_by_mouse[order(-missing_by_mouse$missing_measurements),
                                   c("MouseID", "class", "missing_measurements",
                                     "missing_protein_averages", "partly_missing_proteins")], 5),
      row.names = FALSE)
cat("\nMice with no missing protein readings:\n")
mice_without_missing <- rownames(missing_counts)[rowSums(missing_counts) == 0]
cat("Number of mice with no missing readings:", length(mice_without_missing), "of", nrow(missing_counts), "\n")
cat(if (length(mice_without_missing)) paste(mice_without_missing, collapse = ", ") else "None", "\n")
cat("\nProteins with no missing readings:\n")
proteins_without_missing <- colnames(missing_counts)[colSums(missing_counts) == 0]
cat("Number of proteins with no missing readings:", length(proteins_without_missing), "of", ncol(missing_counts), "\n")
cat(if (length(proteins_without_missing)) paste(proteins_without_missing, collapse = ", ") else "None", "\n")
cat("Number of proteins with at least one missing reading:", proteins_with_missing, "of", ncol(missing_counts), "\n")
cat("\nMice with at least one missing protein average, by class (descriptive):\n")
missingness_by_class_print <- data.frame(
  Class = missingness_by_class$class,
  Mice = missingness_by_class$mice,
  Mice_with_missing_average = missingness_by_class$mice_with_missing_averages,
  Percent_of_mice = round(missingness_by_class$percent_mice_with_missing_averages, 1)
)
print(missingness_by_class_print, row.names = FALSE)
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
cat("They use per-mouse arithmetic averages across the 15 workbook assay readings; interpretation is provisional until the dilution and replicate mapping is confirmed.\n")

cat(sprintf("\nSlide-ready summary: %d mice, %d proteins, %d classes; %d-%d measurements per mouse; %.2f%% of protein readings are missing.\n",
            nrow(mouse_metadata), length(protein_names), nlevels(droplevels(mouse_metadata$class)),
            min(measurements_per_mouse), max(measurements_per_mouse),
            100 * total_missing / total_protein_values))

class_factor <- droplevels(mouse_metadata$class)

# Higuera et al. (2015, PLoS ONE 10:e0129126, Data preprocessing) filled missing
# values with the class mean and removed a t-CS-m mouse with extensive missingness
# and atypical values from an earlier file; that animal is absent from this 72-mouse
# UCI workbook. For our PCA, avoid class-based imputation (it uses the target labels) and
# avoid imputing whole-animal gaps: exclude the six proteins entirely missing in at
# least five mice, then median-impute only the two remaining missing mouse averages.
# Tang, Jana & Fan (2024, arXiv:2408.12564, Section 6.1) also drop these six proteins;
# unlike their replicate-level complete-case analysis, we aggregate within mouse
# first and retain all 72 mice. For classification, refit this feature handling
# inside each training fold.
# Ahmed et al. (2015, PLoS ONE 10:e0119491, Statistical analysis) model replicate
# and dilution measurements within mouse using three-level mixed effects. Mouse
# means are a transparent descriptive summary here, not that inferential model.
# The course PCA notes motivate standardisation when variables have different scales.
duplicate_proteins_to_exclude_from_pca <- unique(exact_duplicate_features$protein_2)
pca_candidate_proteins <- setdiff(
  protein_names,
  c(duplicate_proteins_to_exclude_from_pca, high_missing_proteins)
)
observed_medians <- apply(mouse_means[, pca_candidate_proteins, drop = FALSE], 2,
                          stats::median, na.rm = TRUE)
usable_proteins <- names(observed_medians)[is.finite(observed_medians)]
X <- mouse_means[, usable_proteins, drop = FALSE]
n_pca_imputed <- sum(is.na(X))
for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- observed_medians[colnames(X)[j]]
variable_proteins <- apply(X, 2, stats::sd) > 0
X <- X[, variable_proteins, drop = FALSE]
cat(sprintf("\nPCA preprocessing: %d high-missingness proteins excluded (%s), %d duplicate column excluded; %d variable proteins retained; %d remaining mouse averages median-imputed.\n",
            length(high_missing_proteins), paste(high_missing_proteins, collapse = ", "),
            length(duplicate_proteins_to_exclude_from_pca), ncol(X), n_pca_imputed))

# Choose the boxplot examples by association with the eight-class label rather than
# raw variance, which mainly ranks proteins by measurement scale. This is descriptive
# feature selection; repeat it inside training folds if used for classification.
boxplot_candidates <- factor_effects[
  !factor_effects$protein %in% c(high_missing_proteins,
                                  unique(exact_duplicate_features$protein_2)),
  , drop = FALSE
]
selected_proteins <- head(
  boxplot_candidates$protein[order(boxplot_candidates$eight_class_eta2, decreasing = TRUE)],
  4
)
cat("\nProteins selected for the class boxplots (highest marginal eight-class eta-squared; descriptive):\n")
print(selected_proteins)

# Figure 1. Proteins with the greatest amount of missing measurement data.
n_show <- min(15, nrow(missing_by_protein))
miss_plot <- missing_by_protein[seq_len(n_show), ]
save_base_plot("01_missingness_by_protein.png", 11, 8, function() {
  par(mar = c(5, 10, 4, 2))
  bp <- barplot(rev(miss_plot$missing_percent), names.arg = rev(miss_plot$protein), horiz = TRUE,
                las = 1, col = "#4477AA", border = NA, xlab = "Missing measurements (%)",
                xlim = c(0, max(miss_plot$missing_percent) * 1.18),
                main = "Missingness is concentrated in a small set of proteins")
  text(rev(miss_plot$missing_percent), bp,
       labels = sprintf("%.1f%%", rev(miss_plot$missing_percent)), pos = 4, cex = 0.75)
})

# Figure 2. Missing readings across mice and proteins.
missing_plot <- missing_counts
mice_with_missing <- rownames(missing_plot)[rowSums(missing_plot) > 0]
proteins_with_missing <- colnames(missing_plot)[colSums(missing_plot) > 0]
missing_plot <- missing_plot[mice_with_missing, proteins_with_missing, drop = FALSE]
mice_selected <- mouse_metadata[match(rownames(missing_plot), mouse_metadata$MouseID), , drop = FALSE]
display_order <- order(mice_selected$class, mice_selected$MouseID)
mice_ordered <- mice_selected[display_order, , drop = FALSE]
missing_plot <- missing_plot[display_order, , drop = FALSE]
missing_plot_reversed <- missing_plot[nrow(missing_plot):1, , drop = FALSE]
missing_long <- as.data.frame(as.table(missing_plot_reversed), stringsAsFactors = FALSE)
names(missing_long) <- c("mouse", "protein", "missing_count")
missing_long$mouse <- factor(missing_long$mouse, levels = rownames(missing_plot_reversed))
missing_long$protein <- factor(missing_long$protein, levels = colnames(missing_plot_reversed))
missing_long$cell_label <- ifelse(missing_long$missing_count > 0,
                                  as.character(missing_long$missing_count), "")
missingness_map <- ggplot2::ggplot(missing_long, ggplot2::aes(x = protein, y = mouse, fill = missing_count)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
  ggplot2::geom_text(ggplot2::aes(label = cell_label, colour = missing_count >= 8), size = 2.3) +
  ggplot2::scale_fill_gradientn(
    colours = c("#FFF7F7", "#E58B9A", "#8F1230"),
    limits = c(0, max(readings_per_mouse)),
    breaks = pretty(c(0, max(readings_per_mouse)), n = 4),
    name = "Missing readings"
  ) +
  ggplot2::scale_colour_manual(values = c("#4A0913", "white"), guide = "none") +
  ggplot2::labs(
    title = "Missing readings across mice and proteins",
    subtitle = sprintf("%d mice with missing readings; %d proteins affected",
                       length(mice_with_missing), length(proteins_with_missing)),
    x = "Protein", y = "Mouse"
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5),
    axis.text.y = ggplot2::element_text(size = 4.5),
    legend.position = "top",
    legend.justification = "left",
    plot.title = ggplot2::element_text(face = "bold")
  )
class_breaks <- which(as.character(mice_ordered$class[-1]) !=
                        as.character(mice_ordered$class[-nrow(mice_ordered)]))
if (length(class_breaks)) {
  missingness_map <- missingness_map +
    ggplot2::geom_hline(
      yintercept = nrow(missing_plot_reversed) - class_breaks + 0.5,
      colour = "#555555", linewidth = 0.35
    )
}
save_ggplot(missingness_map, "02_missingness_mouse_protein_map.png", 15, 10)

# Figure 3. Experimental design as a dot matrix: one dot represents one mouse.
design_plot_data <- mouse_metadata
design_plot_data$class <- factor(
  design_plot_data$class,
  levels = c("c-CS-s", "c-CS-m", "c-SC-s", "c-SC-m",
             "t-CS-s", "t-CS-m", "t-SC-s", "t-SC-m")
)
design_plot_data$GenotypeTreatment <- factor(
  paste(design_plot_data$Genotype, design_plot_data$Treatment, sep = " | "),
  levels = c("Control | Saline", "Control | Memantine",
             "Ts65Dn | Saline", "Ts65Dn | Memantine")
)
design_plot_data$dot_index <- ave(seq_len(nrow(design_plot_data)),
                                  design_plot_data$class, FUN = seq_along)
design_plot_data$dot_x <- (design_plot_data$dot_index - 1) %% 5
design_plot_data$dot_y <- (design_plot_data$dot_index - 1) %/% 5
design_labels <- unique(design_plot_data[, c("Behavior", "GenotypeTreatment", "class")])
design_labels$n_mice <- as.integer(mice_by_class[as.character(design_labels$class)])
design_plot <- ggplot2::ggplot(design_plot_data,
                               ggplot2::aes(x = dot_x, y = dot_y, fill = class)) +
  ggplot2::geom_point(shape = 21, size = 3.6, colour = "white", stroke = 0.7) +
  ggplot2::geom_text(
    data = design_labels,
    ggplot2::aes(x = 2, y = 2.15, label = paste0(as.character(class), "\n(n = ", n_mice, ")")),
    inherit.aes = FALSE, size = 3.1
  ) +
  ggplot2::facet_grid(Behavior ~ GenotypeTreatment) +
  ggplot2::scale_fill_manual(values = class_colors, guide = "none") +
  ggplot2::scale_x_continuous(limits = c(-0.5, 4.5), breaks = NULL) +
  ggplot2::scale_y_continuous(limits = c(-0.45, 2.65), breaks = NULL) +
  ggplot2::labs(
    title = "Experimental design and mice per class",
    subtitle = "Each dot represents one mouse; groups combine genotype, behaviour and treatment",
    x = NULL, y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )
save_ggplot(design_plot, "03_experimental_design.png", 12, 8)

# Figure 4. Boxplots for the four proteins most associated with the eight-class label.
if (length(selected_proteins)) {
  save_base_plot("04_selected_protein_boxplots.png", 13, 9, function() {
  oldpar <- par(mfrow = c(2, 2), mar = c(8, 4, 4, 1))
  for (p in selected_proteins) {
    boxplot(mouse_means[, p] ~ class_factor, col = class_colors[levels(class_factor)],
            las = 2, xlab = "Class", ylab = "Per-mouse mean of observed assay readings",
            main = p, outline = TRUE,
            sub = "Descriptive arithmetic summary across the dilution series; one point per mouse")
    observed <- !is.na(mouse_means[, p])
    stripchart(split(mouse_means[observed, p], class_factor[observed]),
               vertical = TRUE, method = "jitter", jitter = 0.12, add = TRUE,
               pch = 21, bg = "white", col = "#333333", cex = 0.8)
    mtext("Points represent individual mice", side = 3, line = 0.2, cex = 0.7)
  }
  par(oldpar)
  })
}

# Figure 5. Distribution of per-protein marginal eta-squared values by factor.
save_base_plot("05_eta_squared_by_factor.png", 9, 7, function() {
  par(mar = c(8, 5, 4, 2))
  boxplot(list(Genotype = factor_effects$genotype_eta2,
               Treatment = factor_effects$treatment_eta2,
               Behavior = factor_effects$behavior_eta2,
               Eight_classes = factor_effects$eight_class_eta2),
          names = c("Genotype", "Treatment", "Behaviour", "8 classes"),
          col = c("#4477AA", "#66CCEE", "#228833", "#CC6677"),
          ylab = "Marginal eta-squared across proteins", xlab = "Factor",
          main = "Exploratory per-protein association with each factor",
          las = 1, outline = TRUE)
    mtext("Associations use per-mouse means across available assay readings",
          side = 1, line = 7, cex = 0.78)
})

# Figure 6. Proteins with the strongest marginal association with at least one experimental factor.
factor_eta_columns <- c("genotype_eta2", "treatment_eta2", "behavior_eta2")
# Exclude pS6_N from this ranking because it is an exact duplicate of ARC_N.
duplicate_redundancies <- unique(exact_duplicate_features$protein_2)
factor_effects_for_plot <- factor_effects[
  !factor_effects$protein %in% duplicate_redundancies,
  , drop = FALSE
]
factor_eta_matrix_all <- as.matrix(factor_effects_for_plot[, factor_eta_columns, drop = FALSE])
max_factor_eta2 <- apply(factor_eta_matrix_all, 1, max, na.rm = TRUE)
top_factor_indices <- order(max_factor_eta2, decreasing = TRUE)[
  seq_len(min(10, nrow(factor_effects_for_plot)))
]
top_factor_effects <- factor_effects_for_plot[top_factor_indices,
                                              c("protein", factor_eta_columns), drop = FALSE]
rownames(top_factor_effects) <- NULL
top_factor_eta_matrix <- as.matrix(top_factor_effects[, factor_eta_columns, drop = FALSE])
top_factor_eta_matrix[!is.finite(top_factor_eta_matrix)] <- NA_real_
factor_effects_long <- rbind(
  data.frame(protein = top_factor_effects$protein, Factor = "Genotype",
             eta_squared = top_factor_effects$genotype_eta2),
  data.frame(protein = top_factor_effects$protein, Factor = "Treatment",
             eta_squared = top_factor_effects$treatment_eta2),
  data.frame(protein = top_factor_effects$protein, Factor = "Behaviour",
             eta_squared = top_factor_effects$behavior_eta2)
)
factor_effects_long$protein <- factor(factor_effects_long$protein,
                                      levels = rev(top_factor_effects$protein))
factor_effects_long$Factor <- factor(factor_effects_long$Factor,
                                     levels = c("Genotype", "Treatment", "Behaviour"))
plot_max_eta2 <- max(top_factor_eta_matrix, na.rm = TRUE)
axis_max_eta2 <- ceiling(plot_max_eta2 * 10) / 10
factor_association_plot <- ggplot2::ggplot(
  factor_effects_long,
  ggplot2::aes(x = eta_squared, y = protein, colour = Factor, shape = Factor)) +
  ggplot2::geom_point(size = 2.8, position = ggplot2::position_dodge(width = 0.58)) +
  ggplot2::scale_colour_manual(values = c(Genotype = "#0072B2", Treatment = "#D55E00",
                                          Behaviour = "#009E73")) +
  ggplot2::scale_shape_manual(values = c(Genotype = 16, Treatment = 17, Behaviour = 15)) +
  ggplot2::scale_x_continuous(limits = c(0, axis_max_eta2),
                              breaks = seq(0, axis_max_eta2, by = 0.1),
                              expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::labs(title = "Top protein associations by experimental factor",
                subtitle = "Top 10 by maximum marginal eta-squared; per-mouse assay-reading means; descriptive associations",
                x = expression("Marginal " * eta^2), y = NULL,
                colour = NULL, shape = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
                  shape = "none") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                 legend.position = "top", legend.justification = "center",
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_blank())
save_ggplot(factor_association_plot, "06_top_factor_associations.png", 12, 8)
cat("\nTop proteins by marginal eta-squared for the three experimental factors:\n")
print(top_factor_effects, row.names = FALSE)

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
  # Flag extreme positions in the first two PCs for review, without excluding mice.
  pca_scores$distance_first_two_pcs <- sqrt(
    (pca_scores$PC1 / stats::sd(pca_scores$PC1))^2 +
      (pca_scores$PC2 / stats::sd(pca_scores$PC2))^2
  )
  pca_score_flags <- pca_scores[
    order(-pca_scores$distance_first_two_pcs),
    c("MouseID", "class", "PC1", "PC2", "distance_first_two_pcs"),
    drop = FALSE
  ][seq_len(min(5L, nrow(pca_scores))), , drop = FALSE]
  cat("Mice furthest from the origin in the first two standardised PC score directions (descriptive flags only):\n")
  pca_score_flags$PC1 <- round(pca_score_flags$PC1, 2)
  pca_score_flags$PC2 <- round(pca_score_flags$PC2, 2)
  pca_score_flags$distance_first_two_pcs <- round(pca_score_flags$distance_first_two_pcs, 2)
  print(pca_score_flags, row.names = FALSE)
  cat("These mice remain in the analysis; an extreme PCA score alone is not a reason for exclusion. Check assay records if a profile is unexpected.\n")
  cat(sprintf("\nPCA after median imputation and standardisation: PC1 %.1f%%; PC2 %.1f%%; combined %.1f%%.\n",
              100 * pca_variance[1], 100 * pca_variance[2], 100 * sum(pca_variance[1:2])))
  cat("Components needed for 80% / 90% cumulative variance:", n_pc_80, "/", n_pc_90, "\n")
  cat("Largest absolute loadings for PC1 and PC2 (for component interpretation):\n")
  for (component in c("PC1", "PC2")) {
    loading_order <- order(abs(pca_loadings[[component]]), decreasing = TRUE)[1:5]
    cat(component, "\n")
    print(data.frame(protein = pca_loadings$protein[loading_order],
                     loading = round(pca_loadings[[component]][loading_order], 3)), row.names = FALSE)
  }
  cat("PCA uses standardised per-mouse arithmetic means across the available workbook readings for description; this summary combines the dilution series and technical replicates. Six highly incomplete proteins and the duplicate pS6_N column are excluded, and two remaining values are median-imputed. For prediction, fit imputation, scaling and PCA on training folds only.\n")

  # Figure 7. Mouse-level protein profiles projected onto the first two PCs by class.
  save_base_plot("07_pca_by_eight_class.png", 9, 7, function() {
    par(mar = c(5, 5, 4, 2), xpd = FALSE)
    plot(pca_scores$PC1, pca_scores$PC2, pch = 19,
         col = class_colors[as.character(pca_scores$class)],
         xlab = sprintf("PC1 (%.1f%%)", 100 * pca_variance[1]),
         ylab = sprintf("PC2 (%.1f%%)", 100 * pca_variance[2]),
         main = "PCA of per-mouse assay-reading means",
         sub = "Assay-reading summaries; overlap in PC1-PC2 does not rule out later separation")
    legend("topright", legend = names(class_colors),
           col = class_colors, pch = 19, bty = "o", bg = "white",
           cex = 0.72, ncol = 2)
  })

  # Figure 8. The same PCA scores coloured by genotype.
  genotype_levels <- levels(droplevels(mouse_metadata$Genotype))
  genotype_colors <- setNames(grDevices::hcl.colors(length(genotype_levels), "Set 2"), genotype_levels)
  save_base_plot("08_pca_by_genotype.png", 9, 7, function() {
    par(mar = c(5, 5, 4, 2), xpd = FALSE)
    plot(pca_scores$PC1, pca_scores$PC2, pch = 19,
         col = genotype_colors[as.character(pca_scores$Genotype)],
         xlab = sprintf("PC1 (%.1f%%)", 100 * pca_variance[1]),
         ylab = sprintf("PC2 (%.1f%%)", 100 * pca_variance[2]),
         main = "PCA scores by genotype",
         sub = "Same two-component projection, coloured by genotype")
    legend("topright", legend = names(genotype_colors),
           col = genotype_colors, pch = 19, bty = "o", bg = "white", cex = 0.85)
  })

  # Figure 9. Individual and cumulative variance explained.
  save_base_plot("09_pca_scree_cumulative.png", 11, 8, function() {
    par(mar = c(5, 5, 5, 6), xpd = FALSE)
    n_scree <- min(20, length(pca_variance))
    scree_x <- seq_len(n_scree)
    scree_max <- max(pca_variance[seq_len(n_scree)], na.rm = TRUE) * 1.12
    red <- "#CC3333"
    plot(scree_x, pca_variance[seq_len(n_scree)], type = "b", pch = 19,
         ylim = c(0, scree_max),
         xlab = "Principal component",
         ylab = "Individual variance explained",
         main = "PCA of per-mouse assay-reading means: individual and cumulative variance",
         cex.main = 1)
    legend("top", inset = c(0, -0.13), xpd = NA,
           legend = c("Individual variance", "Cumulative variance",
                      "70%, 80%, 90% thresholds"),
           col = c("black", red, red),
           pch = c(19, 17, NA),
           lty = c(1, 1, 2),
           bty = "n", horiz = TRUE, cex = 0.78)
    abline(h = scree_max * c(0.7, 0.8, 0.9),
           lty = 2, col = red)
    lines(scree_x, pca_cumulative[seq_len(n_scree)] * scree_max,
          type = "b", pch = 17, col = red, lwd = 2)
    axis(4,
         at = scree_max * seq(0, 1, by = 0.1),
         labels = paste0(seq(0, 100, by = 10), "%"),
         las = 1, col = red, col.axis = red)
    mtext("Cumulative variance explained",
          side = 4, line = 3, col = red)
  })
  # Figure 10. Pairwise-complete correlations among the 20 most variable proteins;
  # exclude both columns in any exact-duplicate pair so the heatmap has no duplicate axes.
  protein_variances_for_corr <- apply(mouse_means, 2, stats::var, na.rm = TRUE)
  duplicate_proteins <- unique(c(exact_duplicate_features$protein_1,
                                 exact_duplicate_features$protein_2))
  corr_ranked_proteins <- names(sort(
    protein_variances_for_corr[is.finite(protein_variances_for_corr) &
                                 !names(protein_variances_for_corr) %in% duplicate_proteins],
    decreasing = TRUE
  ))
  corr_proteins <- head(corr_ranked_proteins, 20)
  corr_data <- mouse_means[, corr_proteins, drop = FALSE]
  corr <- stats::cor(corr_data, use = "pairwise.complete.obs")
  pairwise_n <- crossprod(!is.na(as.matrix(corr_data)))
  dimnames(pairwise_n) <- list(corr_proteins, corr_proteins)
  pairwise_n_table <- data.frame(protein = corr_proteins, pairwise_n,
                                 check.names = FALSE, row.names = NULL)

  all_correlations <- stats::cor(mouse_means[, corr_ranked_proteins, drop = FALSE],
                                 use = "pairwise.complete.obs")
  all_pairwise_n <- crossprod(!is.na(as.matrix(mouse_means[, corr_ranked_proteins, drop = FALSE])))
  correlation_indices <- which(upper.tri(all_correlations), arr.ind = TRUE)
  correlation_pairs <- data.frame(
    protein_1 = rownames(all_correlations)[correlation_indices[, 1]],
    protein_2 = colnames(all_correlations)[correlation_indices[, 2]],
    correlation = all_correlations[correlation_indices],
    n_complete = all_pairwise_n[correlation_indices],
    stringsAsFactors = FALSE
  )
  correlation_pairs <- correlation_pairs[is.finite(correlation_pairs$correlation), , drop = FALSE]
  correlation_pairs <- correlation_pairs[order(-abs(correlation_pairs$correlation)), ]
  correlation_pairs$correlation <- round(correlation_pairs$correlation, 3)
  cat("\nStrongest absolute protein correlations (pairwise-complete per-mouse 15-reading averages; no imputation):\n")
  print(utils::head(correlation_pairs, 8), row.names = FALSE)
  cat("\nNumber of complete mouse pairs for correlations in the heatmap (range): ",
      min(pairwise_n[upper.tri(pairwise_n)]), "-",
      max(pairwise_n[upper.tri(pairwise_n)]), "\n", sep = "")
  cat(sprintf("\nCorrelation heatmap: the %d highest-variance non-duplicate proteins; no imputation.\n",
              length(corr_proteins)))
  save_base_plot("10_correlation_heatmap_top20_variance.png", 13, 12, function() {
    corrplot::corrplot(corr, method = "color", type = "full", order = "original",
                       col = grDevices::colorRampPalette(c("#174A7E", "white", "#A50026"))(201),
                       col.lim = c(-1, 1), addCoef.col = "#222222", number.digits = 2,
                       number.cex = 0.62, tl.col = "#222222", tl.srt = 45,
                       tl.cex = 0.68, cl.pos = "r",
                title = "Correlations between per-mouse 15-reading means", mar = c(0, 0, 2, 0))
  })

  # Figure 11. Compare a robust central interval with the full observed range per mouse.
  interval_central80 <- am_mouse_intervals(
    dados, protein_names, mouse_metadata$MouseID,
    lower_probability = 0.10, upper_probability = 0.90,
    bounds = "central80"
  )
  interval_minmax <- am_mouse_intervals(dados, protein_names,
                                       mouse_metadata$MouseID, bounds = "minmax")
  interval_proteins <- head(selected_proteins, 3)
  interval_width_data <- do.call(rbind, lapply(interval_proteins, function(protein) {
    data.frame(
      MouseID = rep(mouse_metadata$MouseID, 2),
      class = rep(mouse_metadata$class, 2),
      protein = protein,
      bound_type = rep(c("Central 80%", "Full range"),
                       each = nrow(mouse_metadata)),
      width = c(2 * interval_central80$radius[, protein],
                2 * interval_minmax$radius[, protein]),
      stringsAsFactors = FALSE
    )
  }))
  interval_width_data <- interval_width_data[is.finite(interval_width_data$width), , drop = FALSE]
  interval_width_data$class <- factor(interval_width_data$class, levels = class_levels)
  interval_summary <- do.call(rbind, lapply(interval_proteins, function(protein) {
    data.frame(
      protein = protein,
      mice_with_any_readings = sum(interval_central80$observed_count[, protein] > 0),
      minimum_readings_for_one_mouse = min(interval_central80$observed_count[, protein]),
      median_central80_width = stats::median(2 * interval_central80$radius[, protein], na.rm = TRUE),
      median_full_range_width = stats::median(2 * interval_minmax$radius[, protein], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  cat("\nInterval-profile summary for the three proteins selected for descriptive boxplots:\n")
  print(interval_summary, row.names = FALSE)
  cat("Intervals summarise observed spread across available readings; they are not confidence intervals for a mouse's biological mean. Central 80% bounds reduce sensitivity to single extremes; min–max bounds are shown as a sensitivity comparison.\n")

  interval_width_plot <- ggplot2::ggplot(
    interval_width_data,
    ggplot2::aes(x = class, y = width, colour = bound_type)
  ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = interaction(class, bound_type)),
      position = ggplot2::position_dodge(width = 0.82),
      outlier.shape = NA, alpha = 0.18, linewidth = 0.4
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.82),
      size = 1.6, alpha = 0.75
    ) +
    ggplot2::facet_wrap(~ protein, scales = "free_y", ncol = 1) +
    ggplot2::scale_colour_manual(values = c(
      "Central 80%" = "#0072B2",
      "Full range" = "#D55E00"
    )) +
    ggplot2::labs(
      title = "Assay-reading spread across the 15 within-mouse measurements",
      subtitle = "Three technical replicates across five dilution levels; widths summarise assay readings, not biological uncertainty.",
      x = "Experimental class", y = "Observed interval width on the protein scale",
      colour = "Interval definition"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
  save_ggplot(interval_width_plot, "11_interval_widths_by_class.png", 13, 11)


}

# Figure 12. Within-mouse association between protein signal and workbook suffix.
suffix_order_diagnostics$trend <- ifelse(
  suffix_order_diagnostics$median_within_mouse_spearman < 0,
  "Negative median trend", "Non-negative median trend"
)
suffix_order_diagnostics$protein <- stats::reorder(
  suffix_order_diagnostics$protein,
  suffix_order_diagnostics$median_within_mouse_spearman
)
suffix_order_plot <- ggplot2::ggplot(
  suffix_order_diagnostics,
  ggplot2::aes(x = median_within_mouse_spearman, y = protein, colour = trend)
) +
  ggplot2::geom_vline(xintercept = 0, colour = "#777777", linetype = "dashed") +
  ggplot2::geom_point(size = 2.1) +
  ggplot2::scale_colour_manual(values = c(
    "Negative median trend" = "#0072B2",
    "Non-negative median trend" = "#D55E00"
  )) +
  ggplot2::labs(
    title = "Protein signals often decrease across workbook suffix order",
    subtitle = "Each dot is the median within-mouse Spearman correlation across mice; consistent with, but not proof of, dilution order.",
    x = "Median within-mouse Spearman correlation with numeric suffix (1-15)",
    y = "Protein", colour = "Trend"
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    legend.position = "top",
    plot.title = ggplot2::element_text(face = "bold")
  )
save_ggplot(suffix_order_plot, "12_measurement_suffix_trend.png", 12, 14)

cat("\nFinal-project analysis plan:\n")
cat("1. Treat eight-class prediction as primary; compare it with genotype and behaviour-protocol prediction as secondary tasks.\n")
cat("   The behaviour label records the experimental protocol, not whether a mouse successfully learnt.\n")
cat("2. Keep one mouse profile per observation and group all data by mouse in every validation fold.\n")
cat("3. Compare fold-specific missing-data strategies, direct protein profiles, PCA, and repeated-reading interval profiles.\n")
cat("4. Confirm the workbook suffix-to-dilution/replicate mapping and whether the 15 readings should be modelled by their assay structure or summarised with descriptive intervals.\n")
cat(sprintf("With %d proteins and %d mice, standard LDA/QDA covariance estimates may be singular; compare PCA and regularised classifiers instead.\n",
            length(protein_names), nrow(mouse_metadata)))
cat("\nSaved exploratory figures in:", normalizePath(plot_output_dir, winslash = "/"), "\n")

protein_diagnostics <- merge(missing_by_protein, variability_by_protein, by = "protein", all = TRUE, sort = FALSE)
protein_diagnostics <- merge(protein_diagnostics, factor_effects, by = "protein", all = TRUE, sort = FALSE)
if (exists("pca_loadings")) {
  protein_diagnostics <- merge(protein_diagnostics,
                               pca_loadings[, c("protein", "PC1", "PC2")],
                               by = "protein", all = TRUE, sort = FALSE)
}

cat("Detailed diagnostic tables are printed above.\n")
