run_lda_performance <- function(
    df,
    target,
    predictors,
    n_repeats = 100,
    k_folds = 5,
    seed = 1234
) {
    set.seed(seed)
    
    data <- df %>%
        dplyr::select(
            dplyr::all_of(target),
            dplyr::all_of(predictors)
        ) %>%
        tidyr::drop_na()
    
    lda_recipe <- recipe(
        as.formula(paste(target, "~ .")),
        data = data
    )
    
    lda_spec <- discrim_linear() %>%
        set_engine("MASS", prior = c(0.5, 0.5)) %>%
        set_mode("classification")
    
    lda_workflow <- workflow() %>%
        add_recipe(lda_recipe) %>%
        add_model(lda_spec)
    
    folds <- vfold_cv(
        data,
        v = k_folds,
        repeats = n_repeats,
        strata = all_of(target)
    )
    
    metric_set_used <- metric_set(roc_auc, bal_accuracy)
    
    lda_rs <- fit_resamples(
        lda_workflow,
        resamples = folds,
        metrics = metric_set_used,
        control = control_resamples(
            save_pred = TRUE,
            save_workflow = TRUE
        )
    )
    
    cv_results <- collect_metrics(lda_rs) %>%
        dplyr::select(.metric, mean, std_err)
    
    cv_predictions <- collect_predictions(lda_rs)
    
    results <- list(
        cv_summary = cv_results,
        cv_predictions = cv_predictions,
        resample_fit = lda_rs
    )
    
    return(results)
    
}
