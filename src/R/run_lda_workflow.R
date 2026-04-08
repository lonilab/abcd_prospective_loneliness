.check_target <- function(df, target) {
    if (!target %in% names(df)) {
        stop("Target column not found in the dataframe.")
    }
}

.check_factor <- function(df, target) {
    if(!is.factor(df[[target]])) {
        stop("Target must be a factor.")
    }
}

.check_binary <- function(df, target) {
    target_vals <- df[[target]]
    unique_vals <- unique(na.omit(target_vals))
    if (length(unique_vals) != 2) {
        stop("Target must be a factor with 2 levels.")
    }
}

.perform_randomundersampling <- function(df, target, replace) {
    n_minority = min(table(df[[target]]))
    major_label <- names(which.max(table(df[[target]])))
    minor_label <- names(which.min(table(df[[target]])))

    class_major <- df |> 
        dplyr::filter(.data[[target]] == major_label) %>% 
        dplyr::slice_sample(n = n_minority, replace = replace)
    class_minor <- df |> 
        dplyr::filter(.data[[target]] == minor_label)

    return(dplyr::bind_rows(class_major, class_minor))
}

.run_bootstrap_lda <- function(df, target, predictors, n_bootstraps, replace) {
    
    lda_weights <- vector(mode = "list", length = n_bootstraps)
    lda_bal_acc <- rep(NA_real_, n_bootstraps)
    lda_auc <- rep(NA_real_, n_bootstraps)
    
    for (ith in seq_len(n_bootstraps)) {
        cat(sprintf("\rRunning Bootstrapping: %d / %d", ith, n_bootstraps))
        
        set.seed(ith)
        
        rus_df <- .perform_randomundersampling(df, target, replace)

        formula_ <- as.formula(sprintf("%s ~ .", target))
        lda_fit <- rus_df |>
            dplyr::select(
                tidyselect::all_of(target),
                tidyselect::all_of(predictors)
            ) |>
            (\(df) MASS::lda(formula = formula_, data = df))()
        
        pred <- predict(lda_fit, newdata = rus_df)
        
        eval_df <- tibble::tibble(
            truth = factor(rus_df[[target]], levels = c("No", "Yes")),
            .pred_Yes = pred$posterior[, "Yes"],
            .pred_class = factor(pred$class, levels = c("No", "Yes"))
        )
        
        model_perf <- metric_set(roc_auc, bal_accuracy)(
            eval_df,
            truth = truth,
            .pred_Yes,
            estimate = .pred_class
        )
        
        lda_weights[[ith]] <- lda_fit |>
            purrr::pluck("scaling") |>
            as.data.frame() |>
            tibble::rownames_to_column(var = "variable_name") |>
            tibble::as_tibble() |>
            dplyr::rename(!!paste0("iter", ith) := LD1)
        
        lda_bal_acc[[ith]] <- model_perf$.estimate[1]
        lda_auc[[ith]] <- model_perf$.estimate[2]
    }
    
    return(tibble(
        weights = lda_weights,
        bal_acc = lda_bal_acc,
        roc_auc = lda_auc
    ))
        
}

.summarize_results <- function(lda_weights, alpha) {
    
    bs_result <- purrr::reduce(
        lda_weights, dplyr::left_join, by = "variable_name"
    ) |>
        dplyr::rowwise() |>
        dplyr::mutate(
            mean_weight = mean(dplyr::c_across(-1)),
            ci_low = (coxed::bca(dplyr::c_across(-1), 1 - alpha/2))[1],
            ci_high = (coxed::bca(dplyr::c_across(-1), 1 - alpha/2))[2],
            is_sig = dplyr::if_else(ci_low * ci_high > 0, TRUE, FALSE)
        ) |>
        dplyr::ungroup() |>
        dplyr::select(variable_name, mean_weight, ci_low, ci_high, is_sig)
    
    return(bs_result)
}

run_lda_workflow <- function(
    df, target, predictors, n_bootstraps, alpha, replace = TRUE
) {
    
    .check_target(df, target)
    .check_factor(df, target)
    .check_binary(df, target)
    
    lda_weights <- .run_bootstrap_lda(
        df = df, 
        target = target, 
        predictors = predictors,
        n_bootstraps = n_bootstraps,
        replace = replace
    )
    
    cat("\nSummarizing results...\n")
    result = list(
        weights = .summarize_results(lda_weights$weights, alpha),
        bal_acc_mean = mean(lda_weights$bal_acc),
        bal_acc_sd = sd(lda_weights$bal_acc),
        roc_auc_mean = mean(lda_weights$roc_auc),
        roc_auc_sd = sd(lda_weights$roc_auc)
    )
    
    cat("\nDone!\n")
    
    return(result)
    
}
