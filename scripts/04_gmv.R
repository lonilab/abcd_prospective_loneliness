# Environment ------------------------------------------------------------------
library(tidyverse)
library(lme4)
library(ggseg)
library(patchwork)
library(here)
library(glue)
library(patchwork)
library(tidymodels)

source(here("src", "R", "plot_ggseg_brain.R"))
source(here("src", "R", "run_lda_performance.R"))
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
    drop_na() # 8,732 × 98
    


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

gmv_lda_result$bal_acc_mean
gmv_lda_result$bal_acc_sd
gmv_lda_result$roc_auc_mean
gmv_lda_result$roc_auc_sd


# save to csv
gmv_lda_result$lda_weights %>% 
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
    write_csv(here("outputs", "tables", "lda_results_gmv.csv"))

gmv_sig_result <- gmv_lda_result %>% 
    left_join(gmv_variables_df, by = "variable_name") %>% 
    select(ggseg_label, mean_weight, ci_low, ci_high, is_sig) %>% 
    rename(label = ggseg_label) %>% 
    filter(is_sig)
gmv_sig_result

cgmv_res <- dk %>% 
    as_tibble() %>% 
    left_join(gmv_sig_result, by = "label")

sgmv_res <- aseg %>% 
    as_tibble() %>% 
    left_join(gmv_sig_result, by = "label")


cgmv_fig <- plot_ggseg_brain(
    dat = cgmv_res,
    atlas = "dk", 
    fill = "mean_weight", 
    min = -0.4, 
    max = 0.4, 
    break_int = 0.2
)
cgmv_fig

ggsave(
    plot = cgmv_fig,
    filename = here("outputs", "figures", "cgmv_ldaweights.pdf"),
    device = cairo_pdf,
    width = 10,
    height = 4
)

# Logistic Regression ----------------------------------------------------------
lr_gmv_res_list <- list()

for (ith in 1:(length(gmv_regions))) {
    
    cat(sprintf("\rGMV Running: %d / %d", ith, length(gmv_regions)))
    
    gmv_var <- gmv_regions[ith]
    
    f1 <- glue::glue(
        "loneliness_fu_any ~ {gmv_var} + loneliness_bl + demo_sex_v2 + \
        interview_age + race_ethnicity_3 + (1 | site_id_l) + \
        (1 | rel_family_id)"
    ) |> 
        as.formula()
    
    res <- master_gmv_df |> 
        mutate(!!gmv_var := as.numeric(scale(.data[[gmv_var]]))) |> 
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
            variable_name = gmv_var,
            .before = "p"
        )
    
    lr_gmv_res_list[[ith]] <- res
}

gmv_res_df <- tibble(
    variable_name = gmv_lda_result$weights$variable_name,
    lr_weights = lr_gmv_res_list |> reduce(bind_rows) |> pull(cohend),
    lda_weights = gmv_lda_result$weights$mean_weight
)

cor_results <- correlation::correlation(gmv_res_df)
r_val <- cor_results$r[1]
p_val <- cor_results$p[1]
label_text <- sprintf("italic(r) == %.2f*','~~italic(p) == %.2g", r_val, p_val)

gmv_comp_fig <- gmv_res_df |> 
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
        title = "Gray Matter Volume",
        x = "Mean LDA Weights (Bootstrapped, n = 1,000)",
        y = "Cohen's d from Mixed-effects Logistic Regression"
    ) +
    ggthemes::theme_pander() +
    theme(
        plot.margin = margin(5, 5, 5, 5, "mm"),
        plot.title.position = "plot"
    )
gmv_comp_fig


ggsave(
    plot = gmv_comp_fig,
    filename = here("outputs", "figures", "gmv_comp_fig.pdf"),
    width = 6,
    height = 5
)
