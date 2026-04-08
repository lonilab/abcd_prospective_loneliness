# Environment ------------------------------------------------------------------
library(here)
library(tidyverse)
library(gtsummary)
library(glue)
library(flextable)
library(officer)

source(here("src", "R", "is_binary_column.R"))


# Data I/O ---------------------------------------------------------------------
# ABCD data path
# abcd_path <- "/path/to/abcd/folder"
abcd_folder_path <- "/Users/tywong/OneDrive/opendata/abcd" |>
    here("abcd-data-release-5.1/core")

# basic demographics
abcd_y_lt <- here(abcd_folder_path, "abcd-general", "abcd_y_lt.csv") %>%
    read_csv(show_col_types = FALSE)
abcd_p_demo <- here(abcd_folder_path, "abcd-general", "abcd_p_demo.csv") %>%
    read_csv(show_col_types = FALSE)

interview_date_df <- abcd_y_lt |> 
    select(src_subject_id, eventname, interview_date) |> 
    mutate(interview_date = lubridate::mdy(interview_date)) |> 
    filter(
        eventname %in% c(
            "1_year_follow_up_y_arm_1",
            "2_year_follow_up_y_arm_1",
            "3_year_follow_up_y_arm_1"
        )
    ) |> 
    mutate(
        pandemic_period = case_when(
            interview_date < as.Date("2020-03-01") ~ "Prepandemic",
            interview_date >= as.Date("2020-03-01") & interview_date < as.Date("2021-07-01") ~ "Pandemic",
            interview_date >= as.Date("2021-07-01") ~ "Postpandemic",
            TRUE ~ NA_character_
        )
    )

# loneliness data
loneliness_baseline <- here(abcd_folder_path, "mental-health/mh_p_cbcl.csv") |>
    read_csv(show_col_types = FALSE) |>
    select(src_subject_id, eventname, cbcl_q12_p) |>
    drop_na(cbcl_q12_p) |>
    filter(eventname == "baseline_year_1_arm_1") |>
    mutate(loneliness_bl = if_else(cbcl_q12_p == 0, 0, 1)) |>
    select(-eventname) # baseline all prepandemic
table(loneliness_baseline$cbcl_q12_p)
# 0: 10394
# 1: 1308
# 2: 160, around 1.3%

loneliness_firstyear <- here(abcd_folder_path, "mental-health/mh_p_cbcl.csv") |>
    read_csv(show_col_types = FALSE) |>
    drop_na(cbcl_q12_p) |>
    filter(eventname == "1_year_follow_up_y_arm_1") |> 
    mutate(loneliness_firstyear_any = if_else(cbcl_q12_p > 0, 1, 0)) |> 
    select(src_subject_id, loneliness_firstyear_any) 

loneliness_followup <- here(abcd_folder_path, "mental-health/mh_p_cbcl.csv") |>
    read_csv(show_col_types = FALSE) |>
    drop_na(cbcl_q12_p) |>
    filter(
        eventname %in% c(
            "1_year_follow_up_y_arm_1",
            "2_year_follow_up_y_arm_1",
            "3_year_follow_up_y_arm_1"
        )
    ) |>
    left_join(interview_date_df, by = c("src_subject_id", "eventname")) |> 
    select(src_subject_id, eventname, cbcl_q12_p, pandemic_period) |> # only 1.3% of 2
    mutate(loneliness_bin = if_else(cbcl_q12_p > 0, 1, 0)) |>
    arrange(src_subject_id, eventname) |>
    group_by(src_subject_id) |>
    mutate(
        n = n()
    ) |> 
    filter(n == 3) |> 
    mutate(
        loneliness_fu_sum = sum(loneliness_bin),
        loneliness_fu_any = as.integer(any(loneliness_bin == 1)),
        first_loneliness_period = if_else(
            loneliness_fu_any == 1,
            pandemic_period[match(1, loneliness_bin)],
            "Never"
        )
    ) |>
    slice(1) |>
    ungroup() |>
    select(src_subject_id, loneliness_fu_sum, loneliness_fu_any, first_loneliness_period) 


# Preprocessing ----------------------------------------------------------------
master_preprocessed_df <- abcd_y_lt |>
    left_join(abcd_p_demo, by = c("src_subject_id", "eventname")) |>
    left_join(loneliness_baseline, by = "src_subject_id") |>
    left_join(loneliness_firstyear, by = "src_subject_id") |>
    left_join(loneliness_followup, by = "src_subject_id") |>
    filter(eventname == "baseline_year_1_arm_1") |> # 11,868 subjects
    filter(
        demo_sex_v2 %in% c(1, 2), # 11,865 subjects
    ) |>
    mutate(
        across(where(is.numeric), ~na_if(., 777)),
        across(where(is.numeric), ~na_if(., 999))
    ) %>%
    mutate(
        demo_sex_v2 = factor(demo_sex_v2, labels = c("Male", "Female")),
        site_id_l = factor(site_id_l),
        rel_family_id = factor(rel_family_id),
        race_white = if_else(demo_race_a_p___10 == 1, 1, 0),
        race_black = if_else(demo_race_a_p___11 == 1, 1, 0),
        race_asian = if_else(
            (demo_race_a_p___18 + demo_race_a_p___19 + demo_race_a_p___20 +
                 demo_race_a_p___21 + demo_race_a_p___22 + demo_race_a_p___23 +
                 demo_race_a_p___24) > 0, 1, 0
        ),
        race_aian = if_else((demo_race_a_p___12 + demo_race_a_p___13) > 0, 1, 0),
        race_nhpi = if_else(
            (demo_race_a_p___14 + demo_race_a_p___15 + demo_race_a_p___17 +
                 demo_race_a_p___17) > 0, 1, 0
        ),
        race_other = if_else(demo_race_a_p___25 == 1, 1, 0),
        race_mixed = if_else(
            (race_white + race_black + race_asian + race_aian + race_aian + 
                 race_other) > 1, 1, 0
        ),
        race_ethnicity = case_when(
            race_white == 1 & race_mixed == 0 & demo_ethn_v2 == 2 ~ 1,
            race_black == 1 & race_mixed == 0 ~ 2,
            race_asian == 1 & race_mixed == 0 ~ 3,
            demo_ethn_v2 == 1 & race_mixed == 0 ~ 4,
            race_mixed == 1 ~ 5,
            demo_race_a_p___77 | demo_race_a_p___99 ~ NA,
            TRUE ~ 6
        ),
        race_ethnicity_6 = factor(
            race_ethnicity, 
            levels = 1:6,
            label = c("Non-Hispanic White", "Black", "Asian", "Hispanic", 
                      "Mixed", "Other")
        ),
        race_ethnicity_3 = case_when(
            race_white == 1 & race_mixed == 0 & demo_ethn_v2 == 2 ~ 1,
            race_black == 1 & race_mixed == 0 ~ 2,
            demo_race_a_p___77 | demo_race_a_p___99 ~ NA,
            TRUE ~ 3
        ),
        race_ethnicity_3 = factor(
            race_ethnicity_3, 
            levels = 1:3,
            label = c("White", "Black", "Other")
        ),
        loneliness_bl = factor(
            loneliness_bl, levels = c(0, 1), labels = c("No", "Yes")
        ),
        loneliness_fu_any = factor(
            loneliness_fu_any, levels = c(0, 1), labels = c("No", "Yes")
        ),
        loneliness_firstyear_any = factor(
            loneliness_firstyear_any, levels = c(0, 1), labels = c("No", "Yes")
        ),
        first_loneliness_period = factor(
            first_loneliness_period, 
            levels = c("Never", "Prepandemic", "Pandemic", "Postpandemic")
        )
    ) |>
    select(
        src_subject_id, site_id_l, rel_family_id, interview_age, demo_sex_v2,
        race_ethnicity_3, race_ethnicity_6, loneliness_bl, loneliness_fu_sum,
        loneliness_fu_any, loneliness_firstyear_any, first_loneliness_period
    ) %>%
    drop_na() # 9,602 subjects

master_preprocessed_df$first_loneliness_period |> 
    table()

# Never 7098 
# Prepandemic 1662  
# Pandemic 766            
# Postpandemic 76


# Table 1 ----------------------------------------------------------------------
table01 <- master_preprocessed_df |>
    select(loneliness_bl, loneliness_fu_any, interview_age,
           demo_sex_v2, race_ethnicity_3) %>%
    tbl_summary(
        by = loneliness_fu_any,
        label = list(
            loneliness_bl ~ "Loneliness at baseline",
            interview_age ~ "Age (months)",
            demo_sex_v2 ~ "Sex at birth",
            race_ethnicity_3 ~ "Race/Ethnicity"
        ),
        missing = "no",
        statistic = list(
            all_continuous() ~ "{mean} ± {sd}"
        ),
        digits = list(
            all_continuous() ~ 2
        )
    ) |>
    add_overall() |>
    modify_header(label = "**Characteristics**") |>
    modify_spanning_header(
        c(stat_1, stat_2) ~ "**Sum of loneliness in follow-ups**"
    ) |>
    add_p()
table01

# docx page setup
sect_properties <- prop_section(
    page_size = page_size(),
    type = "continuous",
    page_margins = page_mar(
        bottom = 0.5, top = 0.5, right = 0.5, left = 0.5, gutter = 0
    )
)

table01 |>
    as_flex_table() |>
    fontsize(size = 8.5, part = "all") |>
    padding(padding.top = 1, padding.bottom = 1, part = "all") |>
    bold(part = "header") |>
    set_table_properties(width = 1, layout = "autofit") |>
    save_as_docx(
        path = here("outputs", "tables", "table01.docx"),
        pr_section = sect_properties
    )

# check differences


table_s1 <- abcd_p_demo |> 
    left_join(abcd_y_lt, by = c("src_subject_id", "eventname")) |> 
    left_join(loneliness_baseline, by = c("src_subject_id")) |> 
    mutate(included = if_else(
        src_subject_id %in% master_preprocessed_df$src_subject_id, TRUE, FALSE
    )) |> 
    filter(eventname == "baseline_year_1_arm_1") |> 
    filter(
        demo_sex_v2 %in% c(1, 2)
    ) |> 
    mutate(
        demo_sex_v2 = factor(demo_sex_v2, labels = c("Male", "Female")),
        race_white = if_else(demo_race_a_p___10 == 1, 1, 0),
        race_black = if_else(demo_race_a_p___11 == 1, 1, 0),
        race_asian = if_else(
            (demo_race_a_p___18 + demo_race_a_p___19 + demo_race_a_p___20 +
                 demo_race_a_p___21 + demo_race_a_p___22 + demo_race_a_p___23 +
                 demo_race_a_p___24) > 0, 1, 0
        ),
        race_aian = if_else((demo_race_a_p___12 + demo_race_a_p___13) > 0, 1, 0),
        race_nhpi = if_else(
            (demo_race_a_p___14 + demo_race_a_p___15 + demo_race_a_p___17 +
                 demo_race_a_p___17) > 0, 1, 0
        ),
        race_other = if_else(demo_race_a_p___25 == 1, 1, 0),
        race_mixed = if_else(
            (race_white + race_black + race_asian + race_aian + race_aian + 
                 race_other) > 1, 1, 0
        ),
        race_ethnicity = case_when(
            race_white == 1 & race_mixed == 0 & demo_ethn_v2 == 2 ~ 1,
            race_black == 1 & race_mixed == 0 ~ 2,
            race_asian == 1 & race_mixed == 0 ~ 3,
            demo_ethn_v2 == 1 & race_mixed == 0 ~ 4,
            race_mixed == 1 ~ 5,
            demo_race_a_p___77 | demo_race_a_p___99 ~ NA,
            TRUE ~ 6
        ),
        race_ethnicity_6 = factor(
            race_ethnicity, 
            levels = 1:6,
            label = c("Non-Hispanic White", "Black", "Asian", "Hispanic", 
                      "Mixed", "Other")
        ),
        race_ethnicity_3 = case_when(
            race_white == 1 & race_mixed == 0 & demo_ethn_v2 == 2 ~ 1,
            race_black == 1 & race_mixed == 0 ~ 2,
            demo_race_a_p___77 | demo_race_a_p___99 ~ NA,
            TRUE ~ 3
        ),
        race_ethnicity_3 = factor(
            race_ethnicity_3, 
            levels = 1:3,
            label = c("White", "Black", "Other")
        ),
        loneliness_bl = factor(
            loneliness_bl, levels = c(0, 1), labels = c("No", "Yes")
        )
    ) |> 
    select(included, demo_sex_v2, interview_age, race_ethnicity_3, loneliness_bl) |> 
    drop_na() |> 
    tbl_summary(
        by = included,
        label = list(
            loneliness_bl ~ "Loneliness at baseline",
            interview_age ~ "Age (months)",
            demo_sex_v2 ~ "Sex at birth",
            race_ethnicity_3 ~ "Race/Ethnicity"
        ),
        missing = "no"
    ) |> 
    add_p() |> 
    add_q() |> 
    modify_spanning_header(
        c(stat_1, stat_2) ~ "**Included?**"
    )

table_s1 |>
    as_flex_table() |>
    fontsize(size = 8.5, part = "all") |>
    padding(padding.top = 1, padding.bottom = 1, part = "all") |>
    bold(part = "header") |>
    set_table_properties(width = 1, layout = "autofit") |>
    save_as_docx(
        path = here("outputs", "tables", "table_s1.docx"),
        pr_section = sect_properties
    )





# Save -------------------------------------------------------------------------
master_preprocessed_df |>
    write_rds(file = here("data", "processed", "master_preprocessed_df.rds"))

