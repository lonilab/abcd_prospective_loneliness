# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(lme4)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(gghighlight)
library(ggrepel)
library(progress)
library(glue)

source(here("src", "R", "is_binary_column.R"))

# Data I/O ---------------------------------------------------------------------
master_preprocessed_df <- here(
    "data", "processed", "master_preprocessed_df.rds"
) |> 
    read_rds() |> 
    filter(loneliness_bl == "No") # 8,444

# using data release 5.1
# please change the path to ABCD study folder
abcd_folder_path <- "/Users/tywong/OneDrive/opendata/abcd" |>
    here("abcd-data-release-5.1/core")

# exposome variables
exposome_variables_df <- here("data", "raw", "included_variables.xlsx") |>
    readxl::read_excel(sheet = "Exposome")
selected_filenames <- c(unique(pull(exposome_variables_df, table_name)))
selected_variables <- c(pull(exposome_variables_df, variable_name)) # 344

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

# combine
master_exposome_df <- selected_files |>
    reduce(left_join, by = c("src_subject_id", "eventname")) |> 
    select(-eventname) |>
    mutate(
        across(where(is.numeric), ~na_if(., 777)),
        across(where(is.numeric), ~na_if(., 999))
    )


# Preprocessing ----------------------------------------------------------------
# pull participants with more than 15% missing variables
remove_id <- master_exposome_df |>
    is.na() |>
    rowSums() |>
    as_tibble() |>
    mutate(
        src_subject_id = master_exposome_df$src_subject_id,
        prop_missing = value / dim(master_exposome_df)[2]
    ) |>
    filter(prop_missing > 0.15) |>
    pull(src_subject_id) # 170 subjects

master_exposome_cleaned_df <- master_exposome_df |>
    filter(!(src_subject_id %in% remove_id)) |>
    mutate(across(-c(src_subject_id), as.numeric)) |> 
    mutate(
        demo_comb_income_v2 = case_when(
            demo_comb_income_v2 <= 4 ~ 1,
            demo_comb_income_v2 > 4 & demo_comb_income_v2 <= 6 ~ 2,
            demo_comb_income_v2 == 7 ~ 3,
            demo_comb_income_v2 == 8 ~ 4,
            demo_comb_income_v2 == 9 ~ 5,
            demo_comb_income_v2 == 10 ~ 6,
            TRUE ~ NA
        ),
        demo_comb_income_v2 = ordered(
            demo_comb_income_v2,
            levels = 1:6,
            label = c("1. <$25K", "2. $25K-$49K", "3. $50K-$74K", 
                      "4. $75K-$99K", "5. $100K-$199K", "6. $200K+")
        ),
        famhx_ss_fath_mh_history = famhx_ss_fath_prob_alc_p + 
            famhx_ss_fath_prob_dg_p + famhx_ss_fath_prob_dprs_p + 
            famhx_ss_fath_prob_ma_p + famhx_ss_fath_prob_vs_p + 
            famhx_ss_fath_prob_trb_p + famhx_ss_fath_prob_nrv_p + 
            famhx_ss_fath_prob_prf_p + famhx_ss_fath_prob_hspd_p + 
            famhx_ss_fath_prob_scd_p,
        famhx_ss_moth_mh_history = famhx_ss_moth_prob_alc_p + 
            famhx_ss_moth_prob_dg_p + famhx_ss_moth_prob_dprs_p + 
            famhx_ss_moth_prob_ma_p + famhx_ss_moth_prob_vs_p + 
            famhx_ss_moth_prob_trb_p + famhx_ss_moth_prob_nrv_p + 
            famhx_ss_moth_prob_prf_p + famhx_ss_moth_prob_hspd_p + 
            famhx_ss_moth_prob_scd_p,
        devhx_10sum3 = devhx_10a3_p + devhx_10b3_p + devhx_10c3_p +
            devhx_10d3_p + devhx_10e3_p + devhx_10f3_p + devhx_10g3_p +
            devhx_10h3_p + devhx_10i3_p + devhx_10j3_p + devhx_10k3_p +
            devhx_10l3_p + devhx_10m3_p,
        devhx_14sum3 = devhx_14a3_p + devhx_14b3_p + devhx_14c3_p + 
            devhx_14d3_p + devhx_14e3_p + devhx_14f3_p + devhx_14g3_p +
            devhx_14h3_p,
        ksads_ptsd_sum = ksads_ptsd_raw_754_p + ksads_ptsd_raw_755_p + 
            ksads_ptsd_raw_756_p + ksads_ptsd_raw_757_p + ksads_ptsd_raw_758_p +
            ksads_ptsd_raw_759_p + ksads_ptsd_raw_760_p + ksads_ptsd_raw_761_p +
            ksads_ptsd_raw_762_p + ksads_ptsd_raw_763_p + ksads_ptsd_raw_764_p +
            ksads_ptsd_raw_765_p + ksads_ptsd_raw_766_p + ksads_ptsd_raw_767_p +
            ksads_ptsd_raw_768_p + ksads_ptsd_raw_769_p + ksads_ptsd_raw_770_p,
        devhx_16_p = if_else(devhx_16_p == 9990, NA, devhx_16_p),
        birth_weight_kg = birth_weight_lbs * 0.45359237 + birth_weight_oz * 0.0283495231,
    ) |> 
    mutate_if(is_binary_column, factor) |> 
    select(
        -c(birth_weight_lbs, birth_weight_oz)
    )

master_exposome_cleaned_df |> 
    select(-src_subject_id) |> 
    skimr::skim()

# pull variables with near zero variance 
nzv_vars <- master_exposome_cleaned_df |>
    caret::nearZeroVar(names = TRUE)

combined_exposome_df <- master_preprocessed_df |> 
    filter(!(src_subject_id %in% remove_id)) |> 
    left_join(master_exposome_cleaned_df, by = "src_subject_id") |> 
    select(-all_of(nzv_vars))


# ExWAS ------------------------------------------------------------------------
exposome_variables <- master_exposome_cleaned_df %>%
    select(-all_of(nzv_vars)) %>% 
    select(-src_subject_id) %>% 
    colnames() # 294 exposome variables

res_list <- list()
for (ith in seq_along(exposome_variables)) {
    
    cat(sprintf("\rExposome Running: %d / %d", ith, length(exposome_variables)))
    
    env <- exposome_variables[ith]
    
    if (is.numeric(combined_exposome_df[[env]])) {
        combined_exposome_df[[env]] <- combined_exposome_df[[env]] %>% 
            scale() %>% 
            as.numeric()
    }
    
    df_ <- combined_exposome_df |> 
        drop_na(all_of(env))
    
    df_$rel_family_id
    
    f1 <- glue::glue(
        "loneliness_fu_any ~ {env} + demo_sex_v2 + \
        interview_age + race_ethnicity_3 + (1 | site_id_l) + \
        (1 | rel_family_id)"
    ) %>% 
        as.formula()
    
    res_list[[ith]] <- combined_exposome_df %>% 
        glmer(
            formula = f1, 
            family = binomial(link = "logit"),
            control = glmerControl(
                optimizer = "bobyqa",
                optCtrl = list(maxfun = 2e5)
            )
        ) %>% 
        parameters::model_parameters(
            effect = "fixed", 
            verbose = FALSE
        ) %>%
        as_tibble() %>% 
        dplyr::slice(2) %>% 
        mutate(
            cohend = effectsize::logoddsratio_to_d(Coefficient, log = TRUE),
            cohend_low = effectsize::logoddsratio_to_d(CI_low, log = TRUE),
            cohend_high = effectsize::logoddsratio_to_d(CI_high, log = TRUE),
            variable_name = env,
            .before = "p"
        )
}

exposome_res <- reduce(res_list, bind_rows) %>%
    mutate(p.adj = p.adjust(p, method = "bonferroni")) %>%
    left_join(exposome_variables_df, by = "variable_name") %>%
    arrange(category) %>%
    mutate(
        x = factor(1:n()),
        category = case_when(
            variable_name == "famhx_ss_fath_mh_history" ~ "Parental Psychopathology",
            variable_name == "famhx_ss_moth_mh_history" ~ "Parental Psychopathology",
            variable_name == "devhx_10sum3" ~ "Developmental History",
            variable_name == "devhx_14sum3" ~ "Developmental History",
            variable_name == "ksads_ptsd_sum" ~ "Developmental History",
            variable_name == "birth_weight_kg" ~ "Developmental History",
            TRUE ~ category
        )
    )

write_rds(
    exposome_res,
    here("outputs", "caches", "exposome_incident_res.rds")
)
exposome_res_orig <- read_rds(
    here("outputs", "caches", "exposome_res.rds")
)

corr_res <- cor.test(exposome_res_orig$cohend, exposome_res$cohend)

corr_env_fig <- tibble(
    x = exposome_res_orig$cohend,
    y = exposome_res$cohend
) |> 
    ggplot(aes(x = x, y = y)) +
    geom_point(size = 4, alpha = 0.7, color = "grey30", shape = 16) +
    geom_smooth(method = "lm", color = "tomato3", fill = "tomato3") +
    labs(
        title = "A. Correlation of Effect Estimates Across Environemtal Analyses",
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
corr_env_fig
ggsave(
    plot = corr_env_fig,
    filename = here("outputs", "figures", "corr_env_incident_fig.pdf"),
    width = 8,
    height = 6,
    device = cairo_pdf
)

# save results to csv
exposome_res %>% 
    select(category, table_name, variable_name, variable_label, cohend, 
           cohend_low, cohend_high, p, p.adj) %>% 
    mutate(
        cohend = round(cohend, digits = 3),
        cohend_low = round(cohend_low, digits = 3),
        cohend_high = round(cohend_high, digits = 3),
        p = format(p, scientific = TRUE),
        p.adj = format(p.adj, scientific = TRUE)
    ) %>% 
    rename(
        "cohen's d" = "cohend", 
        "95% CI low cohen's d" = "cohend_low",
        "95% CI high cohen's d" = "cohend_high"
    ) %>% 
    write_csv(
        here("outputs", "tables", "lmm_results_fu_incident_environment.csv")
    )

