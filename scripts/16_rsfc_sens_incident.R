# Environment ------------------------------------------------------------------
library(tidyverse)
library(lme4)
library(ggseg)
library(patchwork)
library(here)
library(glue)
library(ggraph)
library(igraph)

source(here("src", "R", "plot_ggseg_brain.R"))
source(here("src", "R", "run_lda_workflow.R"))


# Data I/O ---------------------------------------------------------------------
master_preprocessed_df <- here(
    "data", "processed", "master_preprocessed_df.rds"
) %>% 
    read_rds() |> 
    filter(loneliness_bl == "No")

# using data release 5.1
# please change the path to ABCD study folder
abcd_folder_path <- "/Users/tywong/OneDrive/opendata/abcd" %>%
    here("abcd-data-release-5.1/core")

# general
mri_y_qc_clfind <- here(abcd_folder_path, "imaging", "mri_y_qc_clfind.csv") %>% 
    read_csv(show_col_types = FALSE)
mri_y_qc_incl <- here(abcd_folder_path, "imaging", "mri_y_qc_incl.csv") %>% 
    read_csv(show_col_types = FALSE)
mri_y_adm_info <- here(abcd_folder_path, "imaging", "mri_y_adm_info.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(src_subject_id, eventname, mri_info_deviceserialnumber)
mri_y_qc_motion <- here(abcd_folder_path, "imaging", "mri_y_qc_motion.csv") %>% 
    read_csv(show_col_types = FALSE)

# RSFC variables
rsfc_variables_df <- here("data", "raw", "included_variables.xlsx") %>%
    readxl::read_excel(sheet = "RSFC")
selected_filenames <- c(unique(pull(rsfc_variables_df, table_name)))
selected_variables <- c(pull(rsfc_variables_df, variable_name)) # 416

selected_files <- list()
for (ith in seq_along(selected_filenames)) {
    selected_files[[ith]] <- list.files(
        path = abcd_folder_path,
        full.names = TRUE,
        pattern = glue("{selected_filenames[ith]}.csv"),
        recursive = TRUE
    ) %>%
        read_csv(show_col_types = FALSE) %>%
        select(src_subject_id, eventname, any_of(selected_variables))
}


master_rsfc_df <- selected_files %>%
    reduce(left_join, by = c("src_subject_id", "eventname")) %>% 
    left_join(mri_y_qc_clfind, by = c("src_subject_id", "eventname")) %>% 
    left_join(mri_y_qc_incl, by = c("src_subject_id", "eventname")) %>% 
    left_join(mri_y_adm_info, by = c("src_subject_id", "eventname")) %>% 
    left_join(mri_y_qc_motion, by = c("src_subject_id", "eventname")) %>% 
    filter(eventname == "baseline_year_1_arm_1") %>% 
    mutate(
        across(where(is.numeric), ~na_if(., 777)),
        across(where(is.numeric), ~na_if(., 999))
    ) %>% 
    filter(
        # https://docs.abcdstudy.org/latest/documentation/imaging/type_qc.html
        mrif_score %in% c(1, 2) & imgincl_t1w_include == 1,
        imgincl_rsfmri_include == 1
    ) %>% 
    select(src_subject_id:rsfmri_cor_ngd_vta_scs_vtdcrh, rsfmri_meanmotion, 
           mri_info_deviceserialnumber) %>% 
    right_join(master_preprocessed_df, by = "src_subject_id") %>% 
    drop_na()

find_identical_columns <- function(df) {
    col_names <- colnames(df)
    n <- length(col_names)
    identical_pairs <- list()
    
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            if (identical(df[[i]], df[[j]])) {
                identical_pairs[[length(identical_pairs) + 1]] <- 
                    c(col1 = col_names[i], col2 = col_names[j])
            }
        }
    }
    
    if (length(identical_pairs) == 0) {
        message("No identical columns found.")
        return(invisible(NULL))
    }
    
    do.call(rbind, identical_pairs) |> as.data.frame()
}

identical_labels <- master_rsfc_df %>% 
    find_identical_columns() %>% 
    pull(col2)

rsfc_labels <- master_rsfc_df %>% 
    select(rsfmri_c_ngd_ad_ngd_ad:rsfmri_cor_ngd_vta_scs_vtdcrh) %>% 
    select(-all_of(identical_labels)) %>% 
    colnames() # 33



master_rsfc_resid_df <- master_rsfc_df
for (ith in seq_along(rsfc_labels)) {
    
    cat(sprintf("\rExtract residuals: %d / %d", ith, length(rsfc_labels)))
    
    label <- rsfc_labels[ith]
    
    # perform Fisher's z transformation
    master_rsfc_resid_df[[label]] <- atanh(master_rsfc_resid_df[[label]])
    
    formula_ <- glue::glue(
        "{label} ~ interview_age + demo_sex_v2 + \
        race_ethnicity_3 + rsfmri_meanmotion + \
        (1 | mri_info_deviceserialnumber) + (1 | rel_family_id)"
    ) %>%
        as.formula()
    
    master_rsfc_resid_df[[label]] <- master_rsfc_df %>%
        mutate(
            !!label := as.numeric(scale(.data[[label]]))
        ) %>% 
        lmer(
            formula = formula_,
            control = lmerControl(
                optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)
            )
        ) %>% 
        resid() %>% 
        scale() %>% 
        as.numeric()
}

rsfc_lda_result <- run_lda_workflow(
    df = master_rsfc_resid_df, 
    target = "loneliness_fu_any", 
    predictors = rsfc_labels,
    n_bootstraps = 1000,
    alpha = 0.05,
    replace = TRUE
)

rsfc_lda_result$weights %>% 
    left_join(rsfc_variables_df, by = "variable_name") %>% 
    select(variable_name, from, to, mean_weight, ci_low, ci_high, is_sig) %>% 
     mutate_if(is.numeric, round, digits = 3) %>% 
     rename(
         "LDA weight" = "mean_weight", 
         "95% CI (Lower)" = "ci_low",
         "95% CI (Upper)" = "ci_high",
         "Significant?" = "is_sig"
     ) %>% 
     write_csv(here("outputs", "tables", "lda_results_incident_rsfc.csv"))

rsfc_lda_orig_result <- read_csv(
    here("outputs", "tables", "lda_results_rsfc.csv")
)

corr_res <- cor.test(
    rsfc_lda_orig_result$`LDA weight`,
    rsfc_lda_result$weights$mean_weight
)

corr_rsfc_fig <- tibble(
    x = rsfc_lda_orig_result$`LDA weight`,
    y = rsfc_lda_result$weights$mean_weight
) |> 
    ggplot(aes(x = x, y = y)) +
    geom_point(size = 4, alpha = 0.7, color = "grey30", shape = 16) +
    geom_smooth(method = "lm", color = "tomato3", fill = "tomato3") +
    labs(
        title = "E. Correlation of Effect Estimates Across RSFC Analyses",
        x = "Primary Analysis",
        y = "Sensitivity Analysis (Incident)"
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
corr_rsfc_fig
ggsave(
    plot = corr_rsfc_fig,
    filename = here("outputs", "figures", "corr_rsfc_incident_fig.pdf"),
    width = 8,
    height = 6,
    device = cairo_pdf
)
