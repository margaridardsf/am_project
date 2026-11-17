options(stringsAsFactors = FALSE, scipen = 999)
if (!interactive()) grDevices::pdf(NULL)

if (!file.exists("project_utils.R")) {
  stop("Run this script with the project folder as the working directory.")
}
source("project_utils.R")
am_project_packages(c("readxl", "ggplot2", "kernlab", "ranger"))

plot_output_dir <- "plots_classification_interval_profiles"
if (!dir.exists(plot_output_dir)) dir.create(plot_output_dir, recursive = TRUE)
save_ggplot <- function(plot_object, filename, width = 10, height = 7) {
  print(plot_object)
  ggplot2::ggsave(file.path(plot_output_dir, filename), plot = plot_object,
                  width = width, height = height, units = "in", dpi = 220,
                  bg = "white")
}

source_data <- am_read_mouse_data("Data_Cortex_Nuclear.xls")
raw_data <- source_data$data
protein_names <- source_data$protein_names
mouse_metadata <- am_mouse_metadata(raw_data)
mouse_means <- am_mouse_means(raw_data, protein_names, mouse_metadata$MouseID)
mouse_metadata <- mouse_metadata[match(rownames(mouse_means), mouse_metadata$MouseID), , drop = FALSE]

central_intervals <- am_mouse_intervals(
  raw_data, protein_names, mouse_metadata$MouseID,
  lower_probability = 0.10, upper_probability = 0.90,
  bounds = "central80"
)
full_range_intervals <- am_mouse_intervals(
  raw_data, protein_names, mouse_metadata$MouseID, bounds = "minmax"
)

make_interval_features <- function(interval_object, protein_names) {
  centre <- interval_object$centre
  radius <- interval_object$radius / sqrt(3)
  colnames(centre) <- paste0("C__", protein_names)
  colnames(radius) <- paste0("R__", protein_names)
  cbind(centre, radius)
}

feature_sets <- list(
  "Mean profile" = mouse_means,
  "Central 80% interval" = make_interval_features(central_intervals, protein_names),
  "Full min-max interval" = make_interval_features(full_range_intervals, protein_names)
)

cat("INTERVAL-PROFILE CLASSIFICATION: MOUSE-GROUPED NESTED CROSS-VALIDATION\n")
cat("Source workbook:", source_data$source_path, "\n")
cat("Each mouse contributes one feature profile; no mouse can be split across folds.\n")
cat("The 10th–90th percentile bounds describe a central 80% of available readings. Min–max bounds are a sensitivity analysis. Neither is a confidence interval for a biological mean.\n")
cat("For interval profiles, radius/sqrt(3) represents the uniform-within-interval Wasserstein coordinate; the uniform assumption is explicit and should be confirmed with the instructor.\n")
cat("The PIC iJRR scripts concern interval-valued regression, not classification. This analysis reuses centre/radius representation only; it does not transfer the regression model.\n")
cat("All feature filtering and imputation use the training fold only. Proteins with >5% missing mouse profiles are removed within that fold; remaining gaps use training-fold medians without class labels.\n")
cat("Models: RBF SVM and Random Forest on interval features, plus a distance-based kNN classifier as a simple interval-aware comparator.\n\n")

targets <- c("class", "Genotype")
classifiers <- c("SVM", "RandomForest", "kNN")
repeats <- 5L
outer_folds <- 5L
inner_folds <- 3L
analysis_results <- list()
result_index <- 0L
total_analyses <- length(targets) * length(feature_sets) * length(classifiers)

for (target_name in targets) {
  for (basis_name in names(feature_sets)) {
    for (classifier in classifiers) {
      result_index <- result_index + 1L
      cat(sprintf("[%02d/%02d] %s | %s | %s\n",
                  result_index, total_analyses, target_name, basis_name, classifier))
      result <- am_run_nested_group_cv(
        x = feature_sets[[basis_name]], metadata = mouse_metadata,
        target_name = target_name,
        missing_strategy = "remove_sparse",
        representation = "raw", method = classifier,
        repeats = repeats, outer_folds = outer_folds,
        inner_folds = inner_folds, seed = 20261117L, verbose = FALSE
      )
      result$predictions$representation <- basis_name
      result$repeat_metrics$representation <- basis_name
      result$tuning$representation <- basis_name
      analysis_results[[result_index]] <- result
    }
  }
}

model_summary <- am_summarise_cv_results(analysis_results)
model_summary_print <- model_summary
metric_columns <- grep("^(mean|sd)_", names(model_summary_print), value = TRUE)
model_summary_print[metric_columns] <- lapply(model_summary_print[metric_columns], round, 3)
cat("\nGrouped nested-CV results by feature representation:\n")
print(model_summary_print, row.names = FALSE)
cat("\nThe mean-profile baseline and interval models use the same mice and grouped folds.\n")
cat("Scores are internal estimates; differences should be read cautiously because only 72 mice are available.\n")

performance_long <- rbind(
  data.frame(model_summary[, c("target", "representation", "classifier")],
             metric = "Balanced accuracy", mean = model_summary$mean_balanced_accuracy,
             sd = model_summary$sd_balanced_accuracy),
  data.frame(model_summary[, c("target", "representation", "classifier")],
             metric = "Macro-F1", mean = model_summary$mean_macro_f1,
             sd = model_summary$sd_macro_f1)
)
performance_long$target_plot <- ifelse(
  performance_long$target == "class", "Eight-class label", "Genotype"
)
performance_plot <- ggplot2::ggplot(
  performance_long,
  ggplot2::aes(x = classifier, y = mean, colour = representation,
               group = representation, shape = representation)
) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.55), size = 2.8) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = pmax(0, mean - sd), ymax = pmin(1, mean + sd)),
    position = ggplot2::position_dodge(width = 0.55), width = 0.15
  ) +
  ggplot2::facet_grid(metric ~ target_plot) +
  ggplot2::scale_colour_manual(values = c(
    "Mean profile" = "#009E73",
    "Central 80% interval" = "#0072B2",
    "Full min-max interval" = "#D55E00"
  )) +
  ggplot2::scale_shape_manual(values = c(
    "Mean profile" = 16,
    "Central 80% interval" = 17,
    "Full min-max interval" = 15
  )) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                              expand = ggplot2::expansion(mult = c(0, 0.02))) +
  ggplot2::labs(
    title = "Does within-mouse interval information help classification?",
    subtitle = "Points show mean performance over five repeated nested grouped-CV runs; error bars show one standard deviation",
    x = NULL, y = "Performance", colour = "Mouse-level representation",
    shape = "Mouse-level representation"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")
  )
save_ggplot(performance_plot, "01_mean_vs_interval_performance.png", 14, 9)

all_predictions <- do.call(rbind, lapply(analysis_results, `[[`, "predictions"))
for (target_name in targets) {
  target_summary <- model_summary[model_summary$target == target_name, , drop = FALSE]
  selected <- target_summary[order(-target_summary$mean_macro_f1,
                                   -target_summary$mean_balanced_accuracy), ][1, ]
  selected_predictions <- all_predictions[
    all_predictions$target == target_name &
      all_predictions$representation == selected$representation &
      all_predictions$classifier == selected$classifier, , drop = FALSE
  ]
  levels_target <- levels(mouse_metadata[[target_name]])
  confusion <- table(
    Observed = factor(selected_predictions$observed, levels = levels_target),
    Predicted = factor(selected_predictions$predicted, levels = levels_target)
  )
  row_proportions <- prop.table(confusion, margin = 1)
  confusion_data <- as.data.frame(as.table(confusion), stringsAsFactors = FALSE)
  names(confusion_data) <- c("observed", "predicted", "count")
  confusion_data$proportion <- as.vector(row_proportions)
  confusion_data$label <- sprintf("%d\n%.0f%%", confusion_data$count,
                                  100 * confusion_data$proportion)
  confusion_plot <- ggplot2::ggplot(
    confusion_data,
    ggplot2::aes(x = predicted, y = observed, fill = proportion)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = label),
                       size = if (target_name == "class") 2.8 else 3.8,
                       colour = "#202020", lineheight = 0.9) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#2166AC", limits = c(0, 1),
                                 breaks = seq(0, 1, 0.25), labels = scales::percent) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste("Out-of-fold interval-profile confusion matrix:",
                    if (target_name == "class") "eight-class label" else target_name),
      subtitle = paste0("Selected after comparison: ", selected$representation,
                        " / ", selected$classifier,
                        "\nCounts and row percentages pool repeated out-of-fold predictions; descriptive only"),
      x = "Predicted class", y = "Observed class", fill = "Within-class\nproportion"
    ) +
    ggplot2::theme_minimal(base_size = if (target_name == "class") 9 else 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    )
  save_ggplot(confusion_plot, paste0("02_confusion_", target_name, ".png"),
              if (target_name == "class") 10 else 7,
              if (target_name == "class") 8 else 6)
  cat(sprintf("\nHighest mean macro-F1 configuration for %s: %s with %s; macro-F1 %.3f ± %.3f; balanced accuracy %.3f ± %.3f.\n",
              target_name, selected$representation, selected$classifier,
              selected$mean_macro_f1, selected$sd_macro_f1,
              selected$mean_balanced_accuracy, selected$sd_balanced_accuracy))
  cat("This configuration was selected after comparison, so its pooled confusion matrix is descriptive rather than an independent test estimate.\n")
}

cat("\nInterval meaning and limitations:\n")
cat("- The 15 workbook rows are treated as repeated measurements within one mouse, not as 15 independent mice.\n")
cat("- A 10th–90th percentile interval is a robust summary of the central spread among available readings.\n")
cat("- A min–max interval is sensitive to extreme assay readings and is included only as a sensitivity comparison.\n")
cat("- The centre/radius representation with radius divided by sqrt(3) corresponds to a uniform-within-interval Wasserstein distance; it is a modelling assumption, not an empirically verified distribution.\n")
cat("- The original PIC interval joint robust regression is a regression method and is not used as a classifier.\n")
cat("- Ask the instructor whether the five dilution levels and three technical replicates should be represented separately before treating the interval models as a primary result.\n")
cat("Saved interval-classification figures in:", normalizePath(plot_output_dir, winslash = "/"), "\n")
