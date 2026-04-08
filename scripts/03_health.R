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
        "loneliness_fu_any ~ {health_var} + loneliness_bl + demo_sex_v2 + \
        interview_age + race_ethnicity_3 + (1 | site_id_l) + \
        (1 | rel_family_id)"
    ) |> 
        as.formula()
    
    res <- combined_health_df |> 
        glmer(
            formula = f1, 
            family = binomial(link = "logit"),
            control = glmerControl(tolPwrss = 1e-10)
        )  |> 
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
    here("outputs", "caches", "health_res.rds")
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
    write_csv(here("outputs", "tables", "lmm_results_fu_health.csv"))


# Visualization ----------------------------------------------------------------
health_fig <- health_res |> 
    mutate(
        category = factor(category, levels = sort(unique(category))),
        x = factor(x, levels = x[order(category)])
    ) |> 
    ggplot(aes(x = x, y = cohend, color = category)) +
    geom_point(size = 4, alpha = 0.75) +
    labs(x = "", y = "Effect Size (Cohen's d)") +
    scale_x_discrete(expand = c(0.05, 0)) +
    scale_y_continuous(
        limits = c(-0.05, 0.12),
        breaks = c(-0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4)
    ) +
    theme_clean() +
    theme(
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank(),
        legend.title = element_blank(),
        legend.background = element_rect(colour = NA, fill = NA),
        legend.position = c(0.5, 0.08),
        legend.direction = "horizontal",
        legend.key.size = unit(1, "mm")
    ) + 
    gghighlight(
        p.adj < 0.05,
        unhighlighted_params = list(size = 3, colour = NULL, alpha = 0.1), 
        use_direct_label = FALSE,
        keep_scales = TRUE
    ) +
    geom_text_repel(
        data = health_res |> 
            subset(p.adj < 0.05 & category == "Physical Health") |> 
            arrange(desc(abs(cohend))) |> 
            mutate(label = c(
                "total sleep problems", "excessive somnolence", 
                "sleep-wake transition problems",
                "initiating and maintaining sleep problems",
                "any other medical conditions", "problems with vision", 
                "sleep arousal problems", "total medical conditions", 
                "sleep breathing problems", "pubertal stage", "sleep hyperhydrosis",
                "BMI"
            )),
        aes(label = label),
        nudge_x = c(
            8, 1, -8, -3, 1, -5,
            -8, -8, -8, -2, 1, -2
        ),
        nudge_y = c(
            0.05, 0.1, -0.02, 0.1, 0.06, 0, 
            0, 0.01, -0.05, -0.08, -0.06, 0
        ),
        hjust = 0.5,
        segment.size = 0.25,
        show.legend = FALSE
    ) +
    geom_text_repel(
        mapping = aes(label = label),
        data = health_res |> 
            subset(p.adj < 0.05 & category == "Mental Health") |> 
            arrange(desc(abs(cohend))) |> 
            mutate(label = c(
                "total problems", "internalizing", "externalizing", "social", 
                "anxious/depressed", "attention", "aggressive", "thought", 
                "somatic", "withdrawal depressive", "rule breaking", 
                "bipolar", "prodromal psychosis"
            )),
        nudge_x = c(
            5, -4, 4, -2, 2, 5, 5,
            0, 0, 0, 0, 0, 0
        ),
        nudge_y = c(
            0, 0.05, 0.03, 0.06, 0.05, 0.04, 0,
            -0.05, 0.04, -0.06, -0.06, -0.05, -0.05
        ),
        hjust = 0.5,
        segment.size = 0.25,
        show.legend = FALSE
    )
health_fig

ggsave(
    plot = health_fig,
    filename = here("outputs", "figures", "health_fu_loneliness.pdf"), 
    width = 10, 
    height = 6,
    device = cairo_pdf
)

