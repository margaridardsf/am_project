# Shared data handling, grouped validation and classifier helpers.

am_project_packages <- function(packages, install_missing = TRUE) {
  user_library <- Sys.getenv("AM_PROJECT_R_LIB", unset = "")
  if (!nzchar(user_library)) {
    user_library <- file.path(
      Sys.getenv("LOCALAPPDATA", unset = path.expand("~")),
      "R", "win-library",
      paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
    )
  }
  if (!dir.exists(user_library)) dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_library, .libPaths()))

  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && install_missing) {
    install.packages(missing, lib = user_library, repos = "https://cloud.r-project.org")
  }
  still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing)) {
    stop("Required R packages are unavailable: ", paste(still_missing, collapse = ", "))
  }
  invisible(user_library)
}

am_read_mouse_data <- function(path = "Data_Cortex_Nuclear.xls") {
  if (!file.exists(path)) path <- file.choose()
  read_path <- file.path(tempdir(), "Data_Cortex_Nuclear_for_analysis.xls")
  if (!file.copy(path, read_path, overwrite = TRUE)) {
    stop("Could not copy the source workbook to a temporary path for reading.")
  }
  data <- as.data.frame(readxl::read_excel(read_path, na = c("NA", "")))
  names(data) <- trimws(names(data))

  metadata_names <- c("MouseID", "Genotype", "Treatment", "Behavior", "class")
  missing_names <- setdiff(metadata_names, names(data))
  if (length(missing_names)) {
    stop("Expected metadata columns not found: ", paste(missing_names, collapse = ", "))
  }
  data$MeasurementID <- as.character(data$MouseID)
  data$MouseID <- sub("_[0-9]+$", "", data$MeasurementID)
  for (variable in c("Genotype", "Treatment", "Behavior", "class")) {
    data[[variable]] <- factor(data[[variable]])
  }

  protein_names <- setdiff(names(data), c(metadata_names, "MeasurementID"))
  non_numeric <- protein_names[!vapply(data[protein_names], is.numeric, logical(1))]
  if (length(non_numeric)) {
    stop("Protein columns must be numeric: ", paste(non_numeric, collapse = ", "))
  }
  list(data = data, protein_names = protein_names,
       metadata_names = metadata_names, source_path = normalizePath(path, winslash = "/"))
}

am_mouse_metadata <- function(data) {
  metadata <- data[!duplicated(data$MouseID),
                   c("MouseID", "Genotype", "Treatment", "Behavior", "class"),
                   drop = FALSE]
  rownames(metadata) <- metadata$MouseID
  for (variable in c("Genotype", "Treatment", "Behavior", "class")) {
    number_of_levels <- tapply(as.character(data[[variable]]), data$MouseID,
                               function(values) length(unique(na.omit(values))))
    if (any(number_of_levels > 1)) {
      stop("Metadata variable changes within MouseID: ", variable)
    }
  }
  metadata
}

am_mouse_means <- function(data, protein_names, mouse_ids = unique(data$MouseID)) {
  result <- matrix(NA_real_, nrow = length(mouse_ids), ncol = length(protein_names),
                   dimnames = list(mouse_ids, protein_names))
  for (i in seq_along(mouse_ids)) {
    rows <- data$MouseID == mouse_ids[i]
    result[i, ] <- vapply(data[rows, protein_names, drop = FALSE], function(values) {
      if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
    }, numeric(1))
  }
  result
}

am_mouse_intervals <- function(data, protein_names, mouse_ids = unique(data$MouseID),
                               lower_probability = 0.10, upper_probability = 0.90,
                               bounds = c("central80", "minmax")) {
  bounds <- match.arg(bounds)
  centre <- radius <- lower <- upper <- matrix(
    NA_real_, nrow = length(mouse_ids), ncol = length(protein_names),
    dimnames = list(mouse_ids, protein_names)
  )
  observed_count <- matrix(0L, nrow = length(mouse_ids), ncol = length(protein_names),
                           dimnames = list(mouse_ids, protein_names))

  for (i in seq_along(mouse_ids)) {
    rows <- data$MouseID == mouse_ids[i]
    for (j in seq_along(protein_names)) {
      values <- data[[protein_names[j]]][rows]
      values <- values[is.finite(values)]
      observed_count[i, j] <- length(values)
      if (!length(values)) next
      if (bounds == "central80") {
        limits <- stats::quantile(values,
                                  probs = c(lower_probability, upper_probability),
                                  names = FALSE, type = 7)
      } else {
        limits <- range(values)
      }
      lower[i, j] <- limits[1]
      upper[i, j] <- limits[2]
      centre[i, j] <- mean(limits)
      radius[i, j] <- diff(limits) / 2
    }
  }
  list(lower = lower, upper = upper, centre = centre, radius = radius,
       observed_count = observed_count, bounds = bounds)
}

am_exact_duplicate_columns <- function(x) {
  if (ncol(x) < 2L) return(character())
  duplicates <- character()
  for (j in seq_len(ncol(x) - 1L)) {
    for (k in seq.int(j + 1L, ncol(x))) {
      if (identical(unname(x[[j]]), unname(x[[k]]))) {
        duplicates <- c(duplicates, names(x)[k])
      }
    }
  }
  unique(duplicates)
}

am_make_group_folds <- function(metadata, target, groups = metadata$MouseID,
                                folds = 5L, seed = 20261117L) {
  target <- factor(target)
  if (length(target) != nrow(metadata) || length(groups) != nrow(metadata)) {
    stop("Targets and group identifiers must have one value per metadata row.")
  }
  group_table <- data.frame(group = as.character(groups), target = as.character(target),
                            stringsAsFactors = FALSE)
  group_label_counts <- tapply(group_table$target, group_table$group,
                               function(values) length(unique(values)))
  if (any(group_label_counts != 1L)) stop("A group has more than one target label.")
  group_table <- group_table[!duplicated(group_table$group), , drop = FALSE]
  group_table$fold <- NA_integer_

  set.seed(seed)
  for (level in levels(target)) {
    indices <- which(group_table$target == level)
    if (length(indices) < folds) {
      stop("There are fewer groups than folds for target class: ", level)
    }
    assignment <- sample(rep(seq_len(folds), length.out = length(indices)))
    group_table$fold[indices] <- assignment
  }
  lapply(seq_len(folds), function(fold) {
    held_out_groups <- group_table$group[group_table$fold == fold]
    which(as.character(groups) %in% held_out_groups)
  })
}

am_preprocess_fold <- function(x_train, x_test,
                               missing_strategy = c("retain", "remove_sparse"),
                               representation = c("raw", "pca90"),
                               missing_threshold = 0.05) {
  missing_strategy <- match.arg(missing_strategy)
  representation <- match.arg(representation)
  x_train <- as.data.frame(x_train, check.names = FALSE)
  x_test <- as.data.frame(x_test, check.names = FALSE)
  if (!identical(names(x_train), names(x_test))) stop("Training and test features do not match.")

  duplicate_columns <- am_exact_duplicate_columns(x_train)
  if (length(duplicate_columns)) {
    x_train <- x_train[, !names(x_train) %in% duplicate_columns, drop = FALSE]
    x_test <- x_test[, !names(x_test) %in% duplicate_columns, drop = FALSE]
  }
  if (missing_strategy == "remove_sparse") {
    missing_fraction <- colMeans(is.na(x_train))
    remove <- names(missing_fraction)[missing_fraction > missing_threshold]
    x_train <- x_train[, !names(x_train) %in% remove, drop = FALSE]
    x_test <- x_test[, !names(x_test) %in% remove, drop = FALSE]
  } else {
    remove <- character()
  }

  all_missing <- names(x_train)[vapply(x_train, function(values) all(is.na(values)), logical(1))]
  if (length(all_missing)) {
    x_train <- x_train[, !names(x_train) %in% all_missing, drop = FALSE]
    x_test <- x_test[, !names(x_test) %in% all_missing, drop = FALSE]
  }
  medians <- vapply(x_train, stats::median, numeric(1), na.rm = TRUE)
  for (variable in names(x_train)) {
    x_train[[variable]][is.na(x_train[[variable]])] <- medians[[variable]]
    x_test[[variable]][is.na(x_test[[variable]])] <- medians[[variable]]
  }

  training_sd <- vapply(x_train, stats::sd, numeric(1))
  constant <- names(training_sd)[!is.finite(training_sd) | training_sd <= 1e-12]
  if (length(constant)) {
    x_train <- x_train[, !names(x_train) %in% constant, drop = FALSE]
    x_test <- x_test[, !names(x_test) %in% constant, drop = FALSE]
  }
  if (!ncol(x_train)) stop("No usable predictors remain after fold-specific preprocessing.")

  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  number_components <- NA_integer_
  if (representation == "pca90") {
    if (nrow(x_train) < 3L) stop("At least three training mice are needed for PCA.")
    pca_fit <- stats::prcomp(x_train, center = TRUE, scale. = TRUE)
    variance <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
    number_components <- which(cumsum(variance) >= 0.90)[1]
    if (is.na(number_components)) number_components <- length(variance)
    number_components <- max(1L, min(number_components, ncol(pca_fit$x)))
    x_train <- pca_fit$x[, seq_len(number_components), drop = FALSE]
    x_test <- stats::predict(pca_fit, newdata = x_test)[, seq_len(number_components), drop = FALSE]
    colnames(x_train) <- colnames(x_test) <- paste0("PC", seq_len(number_components))
  }

  list(x_train = x_train, x_test = x_test,
       removed_sparse = remove, removed_duplicates = duplicate_columns,
       removed_all_missing = all_missing, removed_constant = constant,
       number_components = number_components)
}

am_standardise_fold <- function(x_train, x_test) {
  centre <- colMeans(x_train)
  scale <- apply(x_train, 2, stats::sd)
  scale[!is.finite(scale) | scale <= 1e-12] <- 1
  list(x_train = sweep(sweep(x_train, 2, centre, `-`), 2, scale, `/`),
       x_test = sweep(sweep(x_test, 2, centre, `-`), 2, scale, `/`),
       centre = centre, scale = scale)
}

am_standardise_interval_fold <- function(x_train, x_test) {
  if (is.null(colnames(x_train)) ||
      !any(grepl("^C__", colnames(x_train))) ||
      !any(grepl("^R__", colnames(x_train)))) {
    return(am_standardise_fold(x_train, x_test))
  }
  used <- character()
  for (centre_name in colnames(x_train)[grepl("^C__", colnames(x_train))]) {
    protein <- sub("^C__", "", centre_name)
    radius_name <- paste0("R__", protein)
    if (!radius_name %in% colnames(x_train)) next
    centre_mean <- mean(x_train[, centre_name])
    centre_sd <- stats::sd(x_train[, centre_name])
    if (!is.finite(centre_sd) || centre_sd <= 1e-12) centre_sd <- 1
    x_train[, centre_name] <- (x_train[, centre_name] - centre_mean) / centre_sd
    x_test[, centre_name] <- (x_test[, centre_name] - centre_mean) / centre_sd
    x_train[, radius_name] <- x_train[, radius_name] / centre_sd
    x_test[, radius_name] <- x_test[, radius_name] / centre_sd
    used <- c(used, centre_name, radius_name)
  }
  remaining <- setdiff(colnames(x_train), used)
  if (length(remaining)) {
    scaled <- am_standardise_fold(x_train[, remaining, drop = FALSE],
                                  x_test[, remaining, drop = FALSE])
    x_train[, remaining] <- scaled$x_train
    x_test[, remaining] <- scaled$x_test
  }
  list(x_train = x_train, x_test = x_test)
}

am_classification_metrics <- function(observed, predicted, class_levels = levels(factor(observed))) {
  observed <- factor(as.character(observed), levels = class_levels)
  predicted <- factor(as.character(predicted), levels = class_levels)
  confusion <- table(Observed = observed, Predicted = predicted)
  recall <- diag(confusion) / pmax(rowSums(confusion), 1)
  precision <- diag(confusion) / pmax(colSums(confusion), 1)
  f1 <- ifelse(recall + precision == 0, 0, 2 * recall * precision / (recall + precision))
  c(accuracy = sum(diag(confusion)) / sum(confusion),
    balanced_accuracy = mean(recall), macro_f1 = mean(f1))
}

am_classifier_grid <- function(method) {
  switch(method,
    SVM = expand.grid(C = c(0.5, 4), sigma_multiplier = c(0.5, 2),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
    RandomForest = expand.grid(mtry_rule = c("sqrt", "third"),
                               min_node_size = 2L,
                               KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
    RegularisedLogistic = expand.grid(alpha = c(0, 1), lambda = c(0.05, 0.5),
                                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
    kNN = data.frame(k = c(3L, 5L, 7L)),
    stop("Unknown classifier: ", method)
  )
}

am_fit_predict <- function(x_train, x_test, y_train, method, parameters, seed = 1L) {
  y_train <- droplevels(factor(y_train))
  if (nlevels(y_train) < 2L) stop("A training fold contains fewer than two classes.")
  set.seed(seed)
  if (method %in% c("SVM", "RegularisedLogistic", "kNN")) {
    has_interval_features <- !is.null(colnames(x_train)) &&
      any(grepl("^C__", colnames(x_train))) && any(grepl("^R__", colnames(x_train)))
    scaled <- if (has_interval_features && method %in% c("SVM", "kNN")) {
      am_standardise_interval_fold(x_train, x_test)
    } else {
      am_standardise_fold(x_train, x_test)
    }
    x_train <- scaled$x_train
    x_test <- scaled$x_test
  }

  if (method == "SVM") {
    sigma_estimates <- kernlab::sigest(x_train, frac = 0.5)
    base_sigma <- as.numeric(sigma_estimates[2])
    if (!is.finite(base_sigma) || base_sigma <= 0) base_sigma <- 1 / ncol(x_train)
    class_counts <- table(y_train)
    class_weights <- length(y_train) / (nlevels(y_train) * class_counts)
    fit <- kernlab::ksvm(
      x = x_train, y = y_train, type = "C-svc",
      kernel = kernlab::rbfdot(sigma = base_sigma * parameters$sigma_multiplier),
      C = parameters$C, scaled = FALSE,
      class.weights = as.numeric(class_weights) |> stats::setNames(names(class_weights))
    )
    return(list(predicted = as.character(kernlab::predict(fit, x_test)),
                number_components = ncol(x_train)))
  }

  if (method == "RandomForest") {
    mtry <- if (parameters$mtry_rule == "sqrt") {
      max(1L, floor(sqrt(ncol(x_train))))
    } else {
      max(1L, floor(ncol(x_train) / 3))
    }
    class_counts <- table(y_train)
    class_weights <- length(y_train) / (nlevels(y_train) * class_counts)
    names(class_weights) <- names(class_counts)
    fit <- ranger::ranger(
      x = as.data.frame(x_train), y = y_train, num.trees = 300L,
      mtry = min(mtry, ncol(x_train)), min.node.size = parameters$min_node_size,
      classification = TRUE, class.weights = class_weights,
      importance = "none", seed = seed, num.threads = 1L
    )
    return(list(predicted = as.character(stats::predict(fit, data = as.data.frame(x_test))$predictions),
                number_components = ncol(x_train)))
  }

  if (method == "RegularisedLogistic") {
    if (nlevels(y_train) == 2L) {
      numeric_target <- as.integer(y_train == levels(y_train)[2])
      class_counts <- table(y_train)
      observation_weights <- length(y_train) / (nlevels(y_train) * class_counts[as.character(y_train)])
      fit <- glmnet::glmnet(
        x = x_train, y = numeric_target, family = "binomial", alpha = parameters$alpha,
        lambda = parameters$lambda, weights = as.numeric(observation_weights),
        standardize = FALSE
      )
      probabilities <- as.numeric(stats::predict(fit, newx = x_test, s = parameters$lambda,
                                                  type = "response"))
      predicted <- ifelse(probabilities >= 0.5, levels(y_train)[2], levels(y_train)[1])
    } else {
      class_counts <- table(y_train)
      observation_weights <- length(y_train) / (nlevels(y_train) * class_counts[as.character(y_train)])
      fit <- glmnet::glmnet(
        x = x_train, y = y_train, family = "multinomial", alpha = parameters$alpha,
        lambda = parameters$lambda, weights = as.numeric(observation_weights),
        standardize = FALSE
      )
      predicted <- as.character(stats::predict(fit, newx = x_test, s = parameters$lambda,
                                                type = "class"))
    }
    return(list(predicted = as.character(predicted), number_components = ncol(x_train)))
  }

  if (method == "kNN") {
    distances <- as.matrix(stats::dist(rbind(x_train, x_test)))
    training_rows <- seq_len(nrow(x_train))
    test_rows <- nrow(x_train) + seq_len(nrow(x_test))
    class_levels <- levels(y_train)
    predicted <- vapply(test_rows, function(row) {
      distance <- distances[row, training_rows]
      neighbours <- order(distance)[seq_len(min(parameters$k, length(distance)))]
      neighbour_weights <- 1 / pmax(distance[neighbours], 1e-8)
      votes <- vapply(class_levels, function(level) {
        sum(neighbour_weights[as.character(y_train[neighbours]) == level])
      }, numeric(1))
      tied <- class_levels[votes == max(votes)]
      if (length(tied) == 1L) tied else {
        first_tied_neighbour <- which(as.character(y_train[neighbours]) %in% tied)[1]
        as.character(y_train[neighbours[first_tied_neighbour]])
      }
    }, character(1))
    return(list(predicted = predicted, number_components = ncol(x_train)))
  }
  stop("Unknown classifier: ", method)
}

am_run_nested_group_cv <- function(x, metadata, target_name,
                                   missing_strategy = c("retain", "remove_sparse"),
                                   representation = c("raw", "pca90"),
                                   method = c("SVM", "RandomForest", "RegularisedLogistic", "kNN"),
                                   repeats = 5L, outer_folds = 5L, inner_folds = 3L,
                                   seed = 20261117L, verbose = TRUE) {
  missing_strategy <- match.arg(missing_strategy)
  representation <- match.arg(representation)
  method <- match.arg(method)
  metadata <- metadata[match(rownames(x), metadata$MouseID), , drop = FALSE]
  if (anyNA(metadata$MouseID)) stop("Mouse profiles and metadata are not aligned.")
  target <- factor(metadata[[target_name]])
  if (anyNA(target)) stop("Target labels cannot be missing.")
  if (nrow(x) != nrow(metadata) || anyDuplicated(metadata$MouseID)) {
    stop("Classification data must contain one profile per unique mouse.")
  }

  candidate_grid <- am_classifier_grid(method)
  prediction_results <- list()
  tuning_results <- list()
  error_count <- 0L
  result_index <- 0L
  tuning_index <- 0L
  for (repeat_id in seq_len(repeats)) {
    outer <- am_make_group_folds(metadata, target, metadata$MouseID, outer_folds,
                                 seed + repeat_id * 1000L)
    for (fold_id in seq_along(outer)) {
      test_index <- outer[[fold_id]]
      train_index <- setdiff(seq_len(nrow(metadata)), test_index)
      train_meta <- metadata[train_index, , drop = FALSE]
      train_target <- target[train_index]
      inner <- am_make_group_folds(train_meta, train_target, train_meta$MouseID,
                                   inner_folds,
                                   seed + repeat_id * 1000L + fold_id * 31L)
      candidate_scores <- rep(-Inf, nrow(candidate_grid))
      candidate_balanced <- rep(-Inf, nrow(candidate_grid))

      for (candidate_id in seq_len(nrow(candidate_grid))) {
        inner_observed <- character()
        inner_predicted <- character()
        candidate <- as.list(candidate_grid[candidate_id, , drop = FALSE])
        candidate_failed <- FALSE
        for (inner_fold in seq_along(inner)) {
          inner_valid <- inner[[inner_fold]]
          inner_train <- setdiff(seq_len(nrow(train_meta)), inner_valid)
          try_result <- tryCatch({
            prepared <- am_preprocess_fold(
              x[train_index[inner_train], , drop = FALSE],
              x[train_index[inner_valid], , drop = FALSE],
              missing_strategy = missing_strategy, representation = representation
            )
            fitted <- am_fit_predict(
              prepared$x_train, prepared$x_test, train_target[inner_train], method,
              candidate, seed + repeat_id * 10000L + fold_id * 100L +
                candidate_id * 10L + inner_fold
            )
            list(predicted = fitted$predicted,
                 observed = as.character(train_target[inner_valid]))
          }, error = function(error) {
            candidate_failed <<- TRUE
            error_count <<- error_count + 1L
            NULL
          })
          if (is.null(try_result)) break
          inner_observed <- c(inner_observed, try_result$observed)
          inner_predicted <- c(inner_predicted, try_result$predicted)
        }
        if (!candidate_failed && length(inner_observed) == length(train_target)) {
          scores <- am_classification_metrics(train_target, inner_predicted,
                                              class_levels = levels(target))
          candidate_scores[candidate_id] <- scores[["macro_f1"]]
          candidate_balanced[candidate_id] <- scores[["balanced_accuracy"]]
        }
      }
      best_candidate <- if (all(!is.finite(candidate_scores))) 1L else {
        order(-candidate_scores, -candidate_balanced)[1]
      }
      selected_parameters <- as.list(candidate_grid[best_candidate, , drop = FALSE])

      try_result <- tryCatch({
        prepared <- am_preprocess_fold(
          x[train_index, , drop = FALSE], x[test_index, , drop = FALSE],
          missing_strategy = missing_strategy, representation = representation
        )
        fitted <- am_fit_predict(prepared$x_train, prepared$x_test, target[train_index],
                                 method, selected_parameters,
                                 seed + repeat_id * 100000L + fold_id)
        fitted$prepared <- prepared
        fitted
      }, error = function(error) {
        error_count <<- error_count + 1L
        stop("Outer-fold model failed (", target_name, ", ", method, ", repeat ",
             repeat_id, ", fold ", fold_id, "): ", conditionMessage(error))
      })
      result_index <- result_index + 1L
      prediction_results[[result_index]] <- data.frame(
        target = target_name, missing_strategy = missing_strategy,
        representation = representation, classifier = method,
        repeat_id = repeat_id, fold = fold_id,
        MouseID = metadata$MouseID[test_index],
        observed = as.character(target[test_index]),
        predicted = fitted$predicted,
        components = fitted$prepared$number_components,
        stringsAsFactors = FALSE
      )
      tuning_index <- tuning_index + 1L
      tuning_results[[tuning_index]] <- data.frame(
        target = target_name, missing_strategy = missing_strategy,
        representation = representation, classifier = method,
        repeat_id = repeat_id, fold = fold_id,
        selected_parameters = paste(names(selected_parameters),
                                    unlist(selected_parameters), sep = "=", collapse = "; "),
        inner_macro_f1 = candidate_scores[best_candidate],
        inner_balanced_accuracy = candidate_balanced[best_candidate],
        components = fitted$prepared$number_components,
        stringsAsFactors = FALSE
      )
      if (verbose) {
        cat(sprintf("  %s | %s | %s | repeat %d/%d, fold %d/%d: macro-F1 %.3f\n",
                    target_name, missing_strategy, method, repeat_id, repeats,
                    fold_id, outer_folds,
                    am_classification_metrics(target[test_index], fitted$predicted,
                                              class_levels = levels(target))[["macro_f1"]]))
      }
    }
  }

  predictions <- do.call(rbind, prediction_results)
  repeat_results <- lapply(split(predictions, predictions$repeat_id), function(rows) {
    metrics <- am_classification_metrics(rows$observed, rows$predicted,
                                         class_levels = levels(target))
    data.frame(target = target_name, missing_strategy = missing_strategy,
               representation = representation, classifier = method,
               repeat_id = rows$repeat_id[1], t(metrics), row.names = NULL,
               check.names = FALSE)
  })
  list(predictions = predictions,
       repeat_metrics = do.call(rbind, repeat_results),
       tuning = do.call(rbind, tuning_results),
       errors = error_count)
}

am_summarise_cv_results <- function(results) {
  metrics <- do.call(rbind, lapply(results, `[[`, "repeat_metrics"))
  grouping <- interaction(metrics$target, metrics$missing_strategy,
                          metrics$representation, metrics$classifier, drop = TRUE)
  summary <- do.call(rbind, lapply(split(metrics, grouping), function(rows) {
    data.frame(target = rows$target[1], missing_strategy = rows$missing_strategy[1],
               representation = rows$representation[1], classifier = rows$classifier[1],
               mean_accuracy = mean(rows$accuracy), sd_accuracy = stats::sd(rows$accuracy),
               mean_balanced_accuracy = mean(rows$balanced_accuracy),
               sd_balanced_accuracy = stats::sd(rows$balanced_accuracy),
               mean_macro_f1 = mean(rows$macro_f1), sd_macro_f1 = stats::sd(rows$macro_f1),
               stringsAsFactors = FALSE)
  }))
  rownames(summary) <- NULL
  summary[order(summary$target, -summary$mean_macro_f1), , drop = FALSE]
}
