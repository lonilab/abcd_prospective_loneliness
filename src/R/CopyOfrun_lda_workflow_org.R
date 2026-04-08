.perform_undersampling <- function(df, target, frac, replace) {
    
    test_indices <- sample(
        x = nrow(df), 
        size = floor(frac *  nrow(df))
    )
    test <- df[test_indices, ]
    train <- df[-test_indices, ]
    
    n_minority = min(table(train[[target]]))
    fct_levels <- levels(train[[target]])
    class0 = train |> 
        dplyr::filter(.data[[target]] == fct_levels[1]) %>% 
        dplyr::slice_sample(n = n_minority, replace = replace)
    class1 = train |> 
        dplyr::filter(.data[[target]] == fct_levels[2]) |>  
        dplyr::slice_sample(n = n_minority, replace = replace)
    train_ds <- dplyr::bind_rows(class0, class1)
    
    return(
        list(
            test = test,
            train = train,
            train_ds = train_ds
        )
    )
}

.run_lda <- function(df, target, predictors, n_bootstraps, frac, replace) {
    
    coeffs <- vector(mode = "list", length = n_bootstraps)
    train_balanced_acc <- vector(mode = "double", length = n_bootstraps)
    test_balanced_acc <- vector(mode = "double", length = n_bootstraps)
    
    for (ith in seq_len(n_bootstraps)) {
        cat(sprintf("\rRuning Bootstrapping: %d / %d", ith, n_bootstraps))
        
        set.seed(ith)
        splits <- .perform_downsample(df, target, frac, replace)
        
        train_ds <- splits$train_ds
        test <- splits$test
        
        lda_fit <- train_ds |>
            select(all_of(target), all_of(predictors)) |>
            (\(df) MASS::lda(
                formula = loneliness_fu_any ~ ., data = df
            ))()
        
        coeffs[[ith]] <- lda_fit |>
            purrr::pluck("scaling") |>
            as.data.frame() |>
            rownames_to_column(var = "predictor") |>
            tibble::as_tibble() |>
            rename(!!paste0("iter", ith) := LD1)
        
        train_balanced_acc[ith] <- caret::confusionMatrix(
            data = predict(lda_fit)$class,
            reference = train_ds[[target]]
        )$byClass["Balanced Accuracy"]
        
        test_balanced_acc[ith] <- caret::confusionMatrix(
            data = predict(lda_fit, newdata = test)$class, 
            reference = test[[target]]
        )$byClass["Balanced Accuracy"]
        
    }
    
    return(list(
        train_balanced_acc = train_balanced_acc, 
        test_balanced_acc = test_balanced_acc, 
        coeffs = coeffs
    ))
}

.summarize_results <- function(df, target, predictors, n_bootstraps, frac, alpha, replace) {
    
    lda_res <- .run_lda(
        df = df, 
        target = target, 
        predictors = predictors,
        n_bootstraps = n_bootstraps,
        frac = frac,
        replace = replace
    )
    
    coeffs <- reduce(lda_res$coeffs, left_join, by = "predictor") |>
        rowwise() |>
        mutate(
            ci_low = quantile(c_across(-1), alpha / 2),
            ci_high = quantile(c_across(-1), 1 - alpha / 2),
            mean = mean(c_across(-1)),
            is_sig = if_else(ci_low * ci_high > 0, TRUE, FALSE)
        ) |>
        ungroup() |>
        select(predictor, mean, ci_low, ci_high, is_sig)
    
    cat("\nDone!\n")

    
    return(
        list(
            train_balanced_acc = lda_res$train_balanced_acc,
            test_balanced_acc = lda_res$test_balanced_acc,
            coeffs = coeffs           
        )
    )
}

run_lda_workflow <- function(df, target, predictors, n_bootstraps, frac, alpha, replace = TRUE) {
    
    library(dplyr)
    
    if (!is.factor(df[[target]])) {
        stop("Target variable must be a factor with 2 levels.")
    }
    
    result <- .summarize_results(
        df = df,
        target = target,
        predictors = predictors,
        n_bootstraps = n_bootstraps,
        frac = frac,
        alpha = alpha,
        replace = replace
    )
    
    cat(sprintf(
        "\nAverage Balanced Accuracy in Train: %.4f\nAverage Balanced Accuracy in Test: %.4f\n",
        mean(result$train_balanced_acc, na.rm = TRUE),
        mean(result$test_balanced_acc, na.rm = TRUE)
    ))
    
    return(result)
}