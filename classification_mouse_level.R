options(stringsAsFactors = FALSE, scipen = 999)
if (!interactive()) grDevices::pdf(NULL)

if (!file.exists("project_utils.R")) {
  stop("Run this script with the project folder as the working directory.")
}
source("project_utils.R")
am_project_packages(c("readxl", "ggplot2", "kernlab", "ranger", "glmnet"))

plot_output_dir <- "plots_classification_mouse_level"
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
mouse_profiles <- am_mouse_means(raw_data, protein_names, mouse_metadata$MouseID)
mouse_metadata <- mouse_metadata[match(rownames(mouse_profiles), mouse_metadata$MouseID), , drop = FALSE]

if (nrow(raw_data) != 1080L || nrow(mouse_profiles) != 72L ||
    ncol(mouse_profiles) != 77L || anyDuplicated(mouse_metadata$MouseID)) {
  stop("Expected the project workbook version with 1,080 rows, 72 mice and 77 proteins.")
}
cat("MOUSE-LEVEL CLASSIFICATION: NESTED, MOUSE-GROUPED CROSS-VALIDATION\n")
cat("Source workbook:", source_data$source_path, "\n")
cat("One arithmetic mean profile of up to 15 assay readings per mouse:", nrow(mouse_profiles), "mice x",
    ncol(mouse_profiles), "proteins.\n")
cat("The source study describes three technical replicates across five dilution levels per mouse. The workbook lacks explicit dilution/replicate identifiers, so these means are pragmatic baseline features, not biological means from 15 independent replicates.\n")
cat("No class labels, behaviour labels, treatment labels or mouse IDs are predictors.\n")
cat("Each validation fold holds out whole mouse profiles; preprocessing and tuning are refitted within training folds.\n")
cat("The same class-stratified mouse folds are reused across targets so their scores are paired on the same held-out mice.\n")
cat("Comparison: all predictors with training-fold median imputation versus removal of proteins with >5% missing mouse averages; direct profiles versus PCA retaining 90% of training variance.\n")
cat("Classifiers: RBF SVM, class-weighted Random Forest and ridge/elastic-net/lasso-regularised logistic regression.\n")
cat("The eight-class logistic model is an exploratory comparator: inner training folds contain fewer than eight mice per class, so its estimates may be unstable.\n")
cat("The behaviour target is the experimental protocol label; it does not indicate whether a mouse actually learnt.\n\n")

targets <- c("class", "Genotype", "Behavior")
target_display <- c(class = "Eight-class label", Genotype = "Genotype",
                    Behavior = "Behaviour protocol")
missing_strategies <- c("retain", "remove_sparse")
representations <- c("raw", "pca90")
classifiers_by_target <- list(
  class = c("SVM", "RandomForest", "RegularisedLogistic"),
  Genotype = c("SVM", "RandomForest", "RegularisedLogistic"),
  Behavior = c("SVM", "RandomForest", "RegularisedLogistic")
)
repeats <- suppressWarnings(as.integer(Sys.getenv("AM_PROJECT_CV_REPEATS", "5")))
if (is.na(repeats) || repeats < 1L) stop("AM_PROJECT_CV_REPEATS must be a positive integer.")
outer_folds <- 5L
inner_folds <- 3L
analysis_results <- list()
result_index <- 0L
analysis_start <- proc.time()[["elapsed"]]
total_analyses <- sum(vapply(targets, function(target_name) {
  length(missing_strategies) * length(representations) *
    length(classifiers_by_target[[target_name]])
}, numeric(1)))

for (target_name in targets) {
  for (missing_strategy in missing_strategies) {
    for (representation in representations) {
      for (classifier in classifiers_by_target[[target_name]]) {
        result_index <- result_index + 1L
        cat(sprintf("[%02d/%02d] %s | %s | %s | %s\n",
                    result_index, total_analyses, target_name,
                    missing_strategy, representation, classifier))
        analysis_results[[result_index]] <- am_run_nested_group_cv(
          x = mouse_profiles, metadata = mouse_metadata,
          target_name = target_name,
          missing_strategy = missing_strategy,
          representation = representation,
          method = classifier,
          repeats = repeats, outer_folds = outer_folds,
          inner_folds = inner_folds, stratify_name = "class",
          seed = 20261117L, verbose = FALSE
        )
      }
    }
  }
}

model_summary <- am_summarise_cv_results(analysis_results)
model_summary_print <- model_summary
metric_columns <- grep("^(mean|sd)_", names(model_summary_print), value = TRUE)
model_summary_print[metric_columns] <- lapply(model_summary_print[metric_columns], round, 3)
cat("\nRepeated outer-fold performance (mean and SD across complete grouped-CV repeats):\n")
print(model_summary_print, row.names = FALSE)
cat("\nThe eight-class task is primary; genotype and behaviour-protocol tasks are secondary factor-wise comparisons.\n")
cat("Balanced accuracy and macro-F1 weight classes equally; accuracy is included for context.\n")
cat("The reported SD measures variation across repeated partitions, not a formal confidence interval.\n")
cat("Reference balanced-accuracy levels are 1/8 for the eight-class task and 1/2 for each binary task.\n")
for (target_name in targets) {
  target_counts <- table(mouse_metadata[[target_name]])
  cat(sprintf("  %s: %s mice per class; majority-class accuracy reference %.3f.\n",
              target_display[[target_name]],
              paste(paste(names(target_counts), as.integer(target_counts), sep = "="),
                    collapse = ", "),
              max(target_counts) / sum(target_counts)))
}
failed_candidate_fits <- sum(vapply(analysis_results, `[[`, integer(1), "errors"))
cat("Failed inner-fold candidate fits handled during tuning:", failed_candidate_fits, "\n")
if (failed_candidate_fits > 0L) {
  failure_messages <- unique(unlist(lapply(analysis_results, `[[`, "error_messages")))
  cat("Distinct candidate-fit messages:\n")
  print(failure_messages)
}
cat(sprintf("Nested-CV analysis elapsed time: %.1f minutes.\n",
            (proc.time()[["elapsed"]] - analysis_start) / 60))

performance_long <- rbind(
  data.frame(model_summary[, c("target", "missing_strategy", "representation", "classifier")],
             metric = "Balanced accuracy", mean = model_summary$mean_balanced_accuracy,
             sd = model_summary$sd_balanced_accuracy),
  data.frame(model_summary[, c("target", "missing_strategy", "representation", "classifier")],
             metric = "Macro-F1", mean = model_summary$mean_macro_f1,
             sd = model_summary$sd_macro_f1)
)
performance_long$target_plot <- factor(
  unname(target_display[performance_long$target]),
  levels = unname(target_display[targets])
)
performance_plot <- ggplot2::ggplot(
  performance_long,
  ggplot2::aes(x = classifier, y = mean, colour = missing_strategy,
               shape = representation, group = interaction(missing_strategy, representation))
) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.55), size = 2.7) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = pmax(0, mean - sd), ymax = pmin(1, mean + sd)),
    position = ggplot2::position_dodge(width = 0.55), width = 0.16
  ) +
  ggplot2::facet_grid(metric ~ target_plot) +
  ggplot2::scale_colour_manual(
    values = c(retain = "#0072B2", remove_sparse = "#D55E00"),
    labels = c(retain = "All predictors, median imputation",
               remove_sparse = "Remove predictors >5% missing")
  ) +
  ggplot2::scale_shape_manual(values = c(raw = 16, pca90 = 17),
                              labels = c(raw = "Protein profiles", pca90 = "PCA (90%)")) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                              expand = ggplot2::expansion(mult = c(0, 0.02))) +
  ggplot2::labs(
    title = "Mouse-level classification performance",
    subtitle = "Mouse-level means; error bars show SD across five nested grouped-CV runs",
    x = NULL, y = "Performance", colour = "Missing-data strategy", shape = "Predictor representation"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = ggplot2::element_text(face = "bold")
  )
  save_ggplot(performance_plot, "01_model_performance_comparison.png", 17, 9)

all_predictions <- do.call(rbind, lapply(analysis_results, `[[`, "predictions"))
for (target_name in targets) {
  target_summary <- model_summary[model_summary$target == target_name, , drop = FALSE]
  selected <- target_summary[order(-target_summary$mean_macro_f1,
                                   -target_summary$mean_balanced_accuracy), ][1, ]
  selected_predictions <- all_predictions[
    all_predictions$target == target_name &
      all_predictions$missing_strategy == selected$missing_strategy &
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
    ggplot2::geom_text(ggplot2::aes(label = label), size = if (target_name == "class") 2.8 else 3.8,
                       colour = "#202020", lineheight = 0.9) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#2166AC", limits = c(0, 1),
                                 breaks = seq(0, 1, 0.25), labels = scales::percent) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste("Out-of-fold confusion matrix:",
                    target_display[[target_name]]),
      subtitle = paste0("Selected after comparison: ", selected$missing_strategy,
                        " / ", selected$representation, " / ", selected$classifier,
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
  cat(sprintf("\nHighest mean macro-F1 configuration for %s (descriptive model selection): %s, %s, %s; macro-F1 %.3f ± %.3f; balanced accuracy %.3f ± %.3f.\n",
              target_display[[target_name]], selected$missing_strategy, selected$representation,
              selected$classifier, selected$mean_macro_f1, selected$sd_macro_f1,
              selected$mean_balanced_accuracy, selected$sd_balanced_accuracy))
  cat("Its confusion matrix pools one out-of-fold prediction per mouse per repeat; choosing the configuration after comparison makes this matrix descriptive, not an independent test estimate.\n")
}

selected_parameters <- do.call(rbind, lapply(analysis_results, `[[`, "tuning"))
cat("\nMost frequently selected hyperparameter settings (by outer fold):\n")
parameter_counts <- as.data.frame(table(
  selected_parameters$target, selected_parameters$missing_strategy,
  selected_parameters$representation, selected_parameters$classifier,
  selected_parameters$selected_parameters
), stringsAsFactors = FALSE)
names(parameter_counts) <- c("target", "missing_strategy", "representation",
                             "classifier", "parameters", "fold_count")
parameter_counts <- parameter_counts[parameter_counts$fold_count > 0, , drop = FALSE]
parameter_counts <- parameter_counts[order(parameter_counts$target,
                                           parameter_counts$classifier,
                                           -parameter_counts$fold_count), , drop = FALSE]
print(utils::head(parameter_counts, 20), row.names = FALSE)

cat("\nNo separate holdout cohort is available. Treat these grouped nested-CV results as internal estimates for a small dataset, and report their variability and class-wise errors.\n")
cat("Saved classification figures in:", normalizePath(plot_output_dir, winslash = "/"), "\n")
