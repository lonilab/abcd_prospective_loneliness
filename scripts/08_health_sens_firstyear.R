# Environment ------------------------------------------------------------------
library(here)
library(tidyverse)
library(lme4)
library(ggplot2)
library(ggthemes)
library(gghighlight)
library(ggrepel)
library(progress)

source(here("src", "R", "is_binary_column.R"))

# Data I/O ---------------------------------------------------------------------
master_preprocessed_df <- here(
    "data", "processed", "master_preprocessed_df.rds"
) |> 
    read_rds()

# using data release 5.1
# please change the path to ABCD study folder
abcd_folder_path <- "/Users/tywong/OneDrive/opendata/abcd" |>
    here("abcd-data-release-5.1/core")

# health variables
health_variables_df <- here("data", "raw", "included_variables.xlsx") |>
    readxl::read_excel(sheet = "Health")
selected_filenames <- c(unique(pull(health_variables_df, table_name)))
selected_variables <- c(pull(health_variables_df, variable_name)) # 61

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
master_health_df <- selected_files |>
    reduce(left_join, by = c("src_subject_id", "eventname")) |> 
    filter(eventname == "baseline_year_1_arm_1") |> 
    select(-eventname) |>
    mutate(
        across(where(is.numeric), ~na_if(., 777)),
        across(where(is.numeric), ~na_if(., 999))
    )


# Preprocessing ----------------------------------------------------------------
# pull participants with more than 15% missing variables
remove_id <- master_health_df |>
    is.na() |>
    rowSums() |>
    as_tibble() |>
    mutate(
        src_subject_id = master_health_df$src_subject_id,
        prop_missing = value / dim(master_health_df)[2]
    ) |>
    filter(prop_missing > 0.15) |>
    pull(src_subject_id) # 157 subjects

master_health_cleaned_df <- master_health_df |>
    filter(!(src_subject_id %in% remove_id)) |>
    mutate(across(-c(src_subject_id), as.numeric)) |> 
    mutate(
        medhx_total = medhx_2a + medhx_2b + medhx_2d + medhx_2e + medhx_2f +
            medhx_2g + medhx_2h + medhx_2i + medhx_2j + medhx_2k + medhx_2l +
            medhx_2m + medhx_2n + medhx_2o + medhx_2p + medhx_2q + medhx_2r +
            medhx_2s,
        bmi = (anthroweightcalc * 703) / (anthroheightcalc^2),
        bmi = if_else(bmi > 100 | bmi < 5, NA, bmi)
    ) |> 
    mutate_if(is_binary_column, factor)


# pull variables with near zero variance 
nzv_vars <- master_health_cleaned_df |>
    caret::nearZeroVar(names = TRUE)

combined_health_df <- master_preprocessed_df |> 
    filter(!(src_subject_id %in% remove_id)) |>
    left_join(master_health_cleaned_df, by = "src_subject_id") |>
    # I can only calculate puberty stage here as require sex information
    mutate(
        pds_ss_male_category = (pds_p_ss_male_category + pds_y_ss_male_category) / 2,
        pds_ss_female_category = (pds_p_ss_female_category + pds_y_ss_female_category) / 2,
        pds_ss_category = if_else(
            demo_sex_v2 == "Male", pds_ss_male_category, pds_ss_female_category
        )
    ) |> 
    select(
        -all_of(nzv_vars),
        -c(
            pds_p_ss_male_category, pds_y_ss_male_category,
            pds_p_ss_female_category, pds_y_ss_female_category,
            pds_ss_male_category, pds_ss_female_category
        )
    )
    

# ExWAS ------------------------------------------------------------------------
health_variables <- colnames(combined_health_df)[12:60] # 49 health variables

res_list <- list()
for (ith in 1:(length(health_variables))) {
    
    cat(sprintf("\rHealth Running: %d / %d", ith, length(health_variables)))
    
    health_var <- health_variables[ith]
    
    if (is.numeric(combined_health_df[[health_var]])) {
        combined_health_df[[health_var]] <- combined_health_df[[health_var]] |> 
            scale() |> 
            as.numeric()
    }
    
    f1 <- glue::glue(
        "loneliness_firstyear_any ~ {health_var} + loneliness_bl + demo_sex_v2 + \
        interview_age + race_ethnicity_3 + (1 | site_id_l) + \
        (1 | rel_family_id)"
    ) |> 
        as.formula()
    
    res <- combined_health_df |> 
        glmer(
            formula = f1, 
            family = binomial(link = "logit"),
            control = glmerControl(
                optimizer = "bobyqa",
                optCtrl = list(maxfun = 2e5)
            )
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
            variable_name = health_var,
            .before = "p"
        )
    
    res_list[[ith]] <- res
}

health_res <- reduce(res_list, bind_rows) |>
    mutate(p.adj = p.adjust(p, method = "bonferroni")) |>
    left_join(health_variables_df, by = "variable_name") |>
    arrange(category) |>
    mutate(
        x = factor(1:n()),
        category = case_when(
            variable_name == "medhx_total" ~ "Physical Health",
            variable_name == "bmi" ~ "Physical Health",
            variable_name == "pds_ss_category" ~ "Physical Health",
            TRUE ~ category
        )
    )

write_rds(
    health_res,
    here("outputs", "caches", "health_firstyear_res.rds")
)
health_res_orig <- read_rds(
    here("outputs", "caches", "health_res.rds")
)

corr_res <- cor.test(health_res_orig$cohend, health_res$cohend)

corr_health_fig <- tibble(
    x = health_res_orig$cohend,
    y = health_res$cohend
) |> 
    ggplot(aes(x = x, y = y)) +
    geom_point(size = 4, alpha = 0.7, color = "grey30", shape = 16) +
    geom_smooth(method = "lm", color = "tomato3", fill = "tomato3") +
    labs(
        title = "B. Correlation of Effect Estimates Across Health Analyses",
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
corr_health_fig
ggsave(
    plot = corr_health_fig,
    filename = here("outputs", "figures", "corr_health_firstyear_fig.pdf"),
    width = 8,
    height = 6,
    device = cairo_pdf
)

# save results to csv
health_res |> 
    select(category, table_name, variable_name, variable_label, Coefficient,
           cohend, cohend_low, cohend_high, p, p.adj) |> 
    mutate(
        cohend = round(cohend, digits = 3),
        cohend_low = round(cohend_low, digits = 3),
        cohend_high = round(cohend_high, digits = 3),
        p = format(p, scientific = TRUE),
        p.adj = format(p.adj, scientific = TRUE)
    ) |> 
    rename(
        "beta" = "Coefficient",
        "cohen's d" = "cohend", 
        "95% CI low cohen's d" = "cohend_low",
        "95% CI high cohen's d" = "cohend_high"
    ) |> 
    write_csv(here("outputs", "tables", "lmm_results_fu_firstyear_health.csv"))
