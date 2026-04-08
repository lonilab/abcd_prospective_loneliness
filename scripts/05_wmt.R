# Environment ------------------------------------------------------------------
library(tidyverse)
library(lme4)
library(ggseg)
library(patchwork)
library(here)
library(glue)
library(ggsegJHU)

source(here("src", "R", "plot_ggseg_brain.R"))
source(here("src", "R", "run_lda_workflow.R"))


# Data I/O ---------------------------------------------------------------------
master_preprocessed_df <- here(
    "data", "processed", "master_preprocessed_df.rds"
) |> 
    read_rds()

# using data release 5.1
# please change the path to ABCD study folder
abcd_folder_path <- file.path(
    "/Users/tywong/OneDrive/opendata/abcd", 
    "abcd-data-release-5.1/core"
)

# general
mri_y_qc_clfind <- here(abcd_folder_path, "imaging", "mri_y_qc_clfind.csv") |> 
    read_csv(show_col_types = FALSE)
mri_y_qc_incl <- here(abcd_folder_path, "imaging", "mri_y_qc_incl.csv") |> 
    read_csv(show_col_types = FALSE)
mri_y_adm_info <- here(abcd_folder_path, "imaging", "mri_y_adm_info.csv") |> 
    read_csv(show_col_types = FALSE) |> 
    select(src_subject_id, eventname, mri_info_deviceserialnumber)
mri_y_qc_motion <- here(abcd_folder_path, "imaging", "mri_y_qc_motion.csv") |> 
    read_csv(show_col_types = FALSE)

# WMT variables
wmt_variables_df <- here("data", "raw", "included_variables.xlsx") |>
    readxl::read_excel(sheet = "WMT")
selected_filenames <- c(unique(pull(wmt_variables_df, table_name)))
selected_variables <- c(pull(wmt_variables_df, variable_name)) # 84

selected_files <- list()
for (ith in seq_along(selected_filenames)) {
    selected_files[[ith]] <- list.files(
        path = abcd_folder_path,
        full.names = TRUE,
        pattern = glue("{selected_filenames[ith]}.csv"),
        recursive = TRUE
    ) |>
        read_csv(show_col_types = FALSE) |>
        select(src_subject_id, eventname, any_of(selected_variables))
}

master_wmt_df <- selected_files |>
    reduce(left_join, by = c("src_subject_id", "eventname")) |> 
    left_join(mri_y_qc_clfind, by = c("src_subject_id", "eventname")) |> 
    left_join(mri_y_qc_incl, by = c("src_subject_id", "eventname")) |> 
    left_join(mri_y_adm_info, by = c("src_subject_id", "eventname")) |> 
    left_join(mri_y_qc_motion, by = c("src_subject_id", "eventname")) |> 
    filter(eventname == "baseline_year_1_arm_1") |> 
    mutate(
        across(where(is.numeric), ~na_if(., 777)),
        across(where(is.numeric), ~na_if(., 999))
    ) |> 
    filter(
        # https://docs.abcdstudy.org/latest/documentation/imaging/type_qc.html
        mrif_score %in% c(1, 2) & imgincl_t1w_include == 1,
        imgincl_dmri_include == 1
    ) |> 
    select(
        src_subject_id:dmdtifp1_137, dmri_meanmotion,
        mri_info_deviceserialnumber
    ) |> 
    right_join(master_preprocessed_df, by = "src_subject_id") |> 
    drop_na() # 7,952

wmt_labels <- master_wmt_df |> 
    select(dmdtifp1_10:dmdtifp1_137) |> 
    colnames() # 60


# Linear Discriminant Analysis -------------------------------------------------
master_wmt_resid_df <- master_wmt_df

# remove confounding effects
for (ith in seq_along(wmt_labels)) { 
    
    cat(sprintf("\rExtract region residuals: %d / %d", ith, length(wmt_labels)))
    
    region <- wmt_labels[ith]
    
    formula_ <- glue::glue(
        "{region} ~ loneliness_bl + interview_age + demo_sex_v2 + \
        race_ethnicity_3 + dmri_meanmotion + \
        (1 | mri_info_deviceserialnumber) + (1 | rel_family_id)"
    ) |> 
        as.formula()
    
    master_wmt_resid_df[[region]] <- master_wmt_df |>
        mutate(
            !!region := as.numeric(scale(.data[[region]]))
        ) |> 
        lmer(
            formula = formula_,
            control = lmerControl(
                optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)
            )
        ) |> 
        resid() |> 
        scale() |> 
        as.numeric()
}

wmt_lda_result <- run_lda_workflow(
    df = master_wmt_resid_df, 
    target = "loneliness_fu_any", 
    predictors = wmt_labels,
    n_bootstraps = 1000,
    alpha = 0.05,
    replace = TRUE
)

wmt_lda_result$bal_acc_mean
wmt_lda_result$bal_acc_sd
wmt_lda_result$roc_auc_mean
wmt_lda_result$roc_auc_sd

wmt_lda_result$weights |> 
    left_join(wmt_variables_df, by = "variable_name") |>
    select(variable_name, variable_label, hemi, region, mean_weight, ci_low, 
           ci_high, is_sig) |> 
    mutate_if(is.numeric, round, digits = 3) |> 
    rename(
        "LDA weight" = "mean_weight", 
        "95% CI (Lower)" = "ci_low",
        "95% CI (Upper)" = "ci_high",
        "Significant?" = "is_sig"
    ) |> 
    write_csv(here("outputs", "tables", "lda_results_wmm.csv"))

wmt_sig_result <- wmt_lda_result$weights |> 
    left_join(wmt_variables_df, by = "variable_name") |>
    select(variable_name, hemi, region, mean_weight, ci_low, ci_high, is_sig, variable_label) |> 
    filter(is_sig)

colfunc <- colorRampPalette(c("#395D9C", "#358CA7", "white", "#F57A17", "#CE204E"))
colors <- colfunc(100)
wmt_vis <- wmt_sig_result |>
    brain_join(tracula) |> 
    ggplot() + 
    geom_sf(aes(fill = mean_weight)) +
    scale_x_continuous(
        breaks = c(7.52, 10.7, 13.92),
        labels = c("Upper Axial", "Coronal", "Lower Axial")
    ) +
    scale_fill_gradientn(
        colors = colors,
        na.value = "grey85",
        limits = c(-3, 3),
        breaks = seq(-3, 3, 1)
    ) +
    labs(fill = "") +
    theme(
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "mm")
    )
wmt_vis
ggsave(
    plot = wmt_vis,
    filename = here("outputs", "figures", "wmt_ldaweights.pdf"),
    device = cairo_pdf,
    width = 10,
    height = 4
)


# Logistic Regression ----------------------------------------------------------
wmt_lr_res_list <- list()

for (ith in 1:(length(wmt_labels))) {
    
    cat(sprintf("\rWMT Running: %d / %d", ith, length(wmt_labels)))
    
    wmt_var <- wmt_labels[ith]
    
    f1 <- glue::glue(
        "loneliness_fu_any ~ {wmt_var} + loneliness_bl + demo_sex_v2 + \
        interview_age + race_ethnicity_3 + (1 | site_id_l) + \
        (1 | rel_family_id)"
    ) |> 
        as.formula()
    
    res <- master_wmt_df |> 
        mutate(!!wmt_var := as.numeric(scale(.data[[wmt_var]]))) |> 
        glmer(
            formula = f1, 
            family = binomial(link = "logit"),
            control = glmerControl(tolPwrss = 1e-10)
        ) |> 
        parameters::model_parameters(
            effect = "fixed", 
            verbose = FALSE
        ) |>
        as_tibble() |> 
        dplyr::slice(2) |> 
        mutate(
            cohend = effectsize::logoddsratio_to_d(Coefficient, log = TRUE),
            cohend_low = effectsize::logoddsratio_to_d(CI_low, log = TRUE),
            cohend_high = effectsize::logoddsratio_to_d(CI_high, log = TRUE),
            variable_name = wmt_var,
            .before = "p"
        )
    
    wmt_lr_res_list[[ith]] <- res
}

wmt_res_df <- tibble(
    variable_name = wmt_lda_result$weights$variable_name,
    lr_weights = wmt_lr_res_list |> reduce(bind_rows) |> pull(cohend),
    lda_weights = wmt_lda_result$weights$mean_weight
)

cor_results <- correlation::correlation(wmt_res_df)
r_val <- cor_results$r[1]
p_val <- cor_results$p[1]
label_text <- sprintf("italic(r) == %.2f*','~~italic(p) == %.2g", r_val, p_val)

wmt_comp_fig <- wmt_res_df |> 
    ggplot(aes(x = lda_weights, y = lr_weights)) +
    geom_point(size = 4, color = "gray30", alpha = 0.75, shape = 16) +
    geom_smooth(method = "lm", color = "tomato3", fill = "tomato2") +
    annotate(
        "text",
        x = -Inf, y = Inf,
        label = label_text,
        parse = TRUE,
        hjust = -0.1, vjust = 1.5,
        size = 4.5
    ) +
    labs(
        title = "White Matter Microstructure",
        x = "Mean LDA Weights (Bootstrapped, n = 1,000)",
        y = "Cohen's d from Mixed-effects Logistic Regression"
    ) +
    ggthemes::theme_pander() +
    theme(
        plot.margin = margin(5, 5, 5, 5, "mm"),
        plot.title.position = "plot"
    )
wmt_comp_fig


ggsave(
    plot = wmt_comp_fig,
    filename = here("outputs", "figures", "wmt_comp_fig.pdf"),
    width = 6,
    height = 5
)
