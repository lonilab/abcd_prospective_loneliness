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
    select(-eventname) |>
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
    drop_na()

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
    target = "loneliness_firstyear_any", 
    predictors = wmt_labels,
    n_bootstraps = 1000,
    alpha = 0.05,
    replace = TRUE
)

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
    write_csv(here("outputs", "tables", "lda_results_firstyear_wmm.csv"))

wmt_lda_orig_result <- read_csv(here("outputs", "tables", "lda_results_wmm.csv"))

corr_res <- cor.test(
    wmt_lda_orig_result$`LDA weight`,
    wmt_lda_result$weights$mean_weight
)

corr_wmm_fig <- tibble(
    x = wmt_lda_orig_result$`LDA weight`,
    y = wmt_lda_result$weights$mean_weight
) |> 
    ggplot(aes(x = x, y = y)) +
    geom_point(size = 4, alpha = 0.7, color = "grey30", shape = 16) +
    geom_smooth(method = "lm", color = "tomato3", fill = "tomato3") +
    labs(
        title = "D. Correlation of Effect Estimates Across WMM Analyses",
        x = "Primary Analysis",
        y = "Sensitivity Analysis (Restricted to First-Year)"
    ) +
    annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = paste0(
            "italic(r) == ", round(corr_res$estimate, 2),
            "*','~~italic(p) == ", signif(corr_res$p.value, 2)
        ),
        hjust = -0.15,
        vjust = 3,
        size = 5,
        parse = TRUE
    ) +
    ggthemes::theme_pander() +
    theme(
        plot.title.position = "plot",
        plot.margin = margin(5, 5, 5, 5, "mm")
    )
corr_wmm_fig
ggsave(
    plot = corr_wmm_fig,
    filename = here("outputs", "figures", "corr_wmm_firstyear_fig.pdf"),
    width = 8,
    height = 6,
    device = cairo_pdf
)


