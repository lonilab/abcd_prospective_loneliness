# Environment ------------------------------------------------------------------
library(tidyverse)
library(lme4)
library(ggseg)
library(patchwork)
library(here)
library(glue)
library(patchwork)

source(here("src", "R", "plot_ggseg_brain.R"))
source(here("src", "R", "run_lda_workflow.R"))

.winsorize_sd <- function(x, sd_cutoff = 3) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    lower <- m - sd_cutoff * s
    upper <- m + sd_cutoff * s
    x[x < lower] <- lower
    x[x > upper] <- upper
    x
}

# Data I/O ---------------------------------------------------------------------
master_preprocessed_df <- here(
    "data", "processed", "master_preprocessed_df.rds"
) %>% 
    read_rds()

# using data release 5.1
# please change the path to ABCD study folder
abcd_folder_path <- file.path(
    "/Users/tywong/OneDrive/opendata/abcd", 
    "abcd-data-release-5.1/core"
)

# general
mri_y_qc_clfind <- here(abcd_folder_path, "imaging", "mri_y_qc_clfind.csv") %>% 
    read_csv(show_col_types = FALSE)
mri_y_qc_incl <- here(abcd_folder_path, "imaging", "mri_y_qc_incl.csv") %>% 
    read_csv(show_col_types = FALSE)
mri_y_adm_info <- here(abcd_folder_path, "imaging", "mri_y_adm_info.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(src_subject_id, eventname, mri_info_deviceserialnumber)

# GMV variables
gmv_variables_df <- here("data", "raw", "included_variables.xlsx") %>%
    readxl::read_excel(sheet = "GMV")
selected_filenames <- c(unique(pull(gmv_variables_df, table_name)))
selected_variables <- c(pull(gmv_variables_df, variable_name))

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

# combine
master_gmv_df <- selected_files %>%
    reduce(left_join, by = c("src_subject_id", "eventname")) %>% 
    left_join(mri_y_qc_clfind, by = c("src_subject_id", "eventname")) %>% 
    left_join(mri_y_qc_incl, by = c("src_subject_id", "eventname")) %>% 
    left_join(mri_y_adm_info, by = c("src_subject_id", "eventname")) %>% 
    filter(eventname == "baseline_year_1_arm_1") %>% 
    select(-eventname) %>% 
    select(
        src_subject_id:smri_vol_scs_intracranialv, mrif_score, 
        imgincl_t1w_include, mri_info_deviceserialnumber
    ) %>%
    filter(
        # https://docs.abcdstudy.org/latest/documentation/imaging/type_qc.html
        mrif_score %in% c(1, 2) & imgincl_t1w_include == 1
    ) %>%
    right_join(master_preprocessed_df, by = "src_subject_id") %>% 
    drop_na() |> # 8,732 × 98
    filter(first_loneliness_period %in% c("Never", "Prepandemic"))
    


gmv_regions <- master_gmv_df %>%
    select(smri_vol_cdk_banksstslh:smri_vol_scs_tprh) %>% 
    colnames()


# Linear Discriminant Analysis -------------------------------------------------
master_gmv_resid_df <- master_gmv_df

# remove confounding of each region
for (ith in seq_along(gmv_regions)) {
    
    cat(sprintf("\rExtract region residuals: %d / %d", ith, length(gmv_regions)))
    
    region <- gmv_regions[ith]
    
    formula_ <- glue::glue(
        "{region} ~ loneliness_bl + interview_age + demo_sex_v2 + \
        race_ethnicity_3 + smri_vol_scs_intracranialv + \
        (1 | mri_info_deviceserialnumber) + (1 | rel_family_id)"
    ) %>% 
        as.formula()
    
    master_gmv_resid_df[[region]] <- master_gmv_df %>%
        mutate(
            !!region := as.numeric(scale(.data[[region]])),
            smri_vol_scs_intracranialv = as.numeric(
                scale(smri_vol_scs_intracranialv)
            ),
        ) %>% 
        lmer(
            formula = formula_,
            control = lmerControl(
                optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)
            )
        ) %>% 
        resid() %>% 
        .winsorize_sd(sd_cutoff = 3) %>% 
        scale() %>% 
        as.numeric()
}

gmv_lda_result <- run_lda_workflow(
    df = master_gmv_resid_df, 
    target = "loneliness_fu_any", 
    predictors = gmv_regions,
    n_bootstraps = 1000,
    alpha = 0.05,
    replace = TRUE
)

# save to csv
gmv_lda_result$weights %>% 
    left_join(gmv_variables_df, by = "variable_name") %>% 
    select(variable_name, variable_label, ggseg_label, mean_weight, ci_low, 
           ci_high, is_sig) %>% 
    mutate_if(is.numeric, round, digits = 3) %>% 
    rename(
        "LDA weight" = "mean_weight", 
        "95% CI (Lower)" = "ci_low",
        "95% CI (Upper)" = "ci_high",
        "Significant?" = "is_sig"
    ) %>% 
    write_csv(here("outputs", "tables", "lda_prepandemic_results_gmv.csv"))

gmv_lda_orig_result <- read_csv(here("outputs", "tables", "lda_results_gmv.csv"))

corr_res <- cor.test(
    gmv_lda_orig_result$`LDA weight`, 
    gmv_lda_result$weights$mean_weight
)

corr_gmv_fig <- tibble(
    x = gmv_lda_orig_result$`LDA weight`,
    y = gmv_lda_result$weights$mean_weight
) |> 
    ggplot(aes(x = x, y = y)) +
    geom_point(size = 4, alpha = 0.7, color = "grey30", shape = 16) +
    geom_smooth(method = "lm", color = "tomato3", fill = "tomato3") +
    labs(
        title = "C. Correlation of Effect Estimates Across GMV Analyses",
        x = "Primary Analysis",
        y = "Sensitivity Analysis (Pre-pandemic Only)"
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
corr_gmv_fig
ggsave(
    plot = corr_gmv_fig,
    filename = here("outputs", "figures", "corr_gmv_prepandemic_fig.pdf"),
    width = 8,
    height = 6,
    device = cairo_pdf
)
