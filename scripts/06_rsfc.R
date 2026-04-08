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
    read_rds()

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
    drop_na() # only complete data

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
        "{label} ~ loneliness_bl + interview_age + demo_sex_v2 + \
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

rsfc_lda_result$bal_acc_mean
rsfc_lda_result$bal_acc_sd
rsfc_lda_result$roc_auc_mean
rsfc_lda_result$roc_auc_sd

rsfc_lda_result %>% 
    left_join(rsfc_variables_df, by = "variable_name") %>% 
    select(variable_name, from, to, mean_weight, ci_low, ci_high, is_sig) %>% 
     mutate_if(is.numeric, round, digits = 3) %>% 
     rename(
         "LDA weight" = "mean_weight", 
         "95% CI (Lower)" = "ci_low",
         "95% CI (Upper)" = "ci_high",
         "Significant?" = "is_sig"
     ) %>% 
     write_csv(here("outputs", "tables", "lda_results_rsfc.csv"))

rsfc_lda_sig_result <- rsfc_lda_result %>% 
    left_join(rsfc_variables_df, by = "variable_name") %>% 
    select(variable_name, from, to, mean_weight, ci_low, ci_high, is_sig) %>% 
    filter(is_sig) %>% 
    select(from, to, mean_weight) %>% 
    mutate(from = case_when(
        from == "cingulo-parietal" ~ "CPT",
        from == "fronto-parietal" ~ "FPN",
        from == "salience" ~ "SAL",
        from == "sensorimotor mouth" ~ "SML",
        from == "cingulo-opercular" ~ "COP"
    )) %>% 
    mutate(to = case_when(
        to == "visual" ~ "VIS",
        to == "none" ~ "NON",
        to == "salience" ~ "SAL",
        to == "left-pallidum" ~ "L Pallidum",
        to == "right-caudate" ~ "R Caudate"
    )) %>% 
    mutate(col = ifelse(mean_weight > 0, "#CE204E", "#395D9C"))

library(circlize)

node_colors <- c(
    CPT = "#006dfe",
    NON = "#c1c0be",
    FPN = "#f3e601",
    SAL = "#0a0a08",
    SML = "#ff7f00",
    VIS = "#3a469a",
    COP = "#810080",
    `L Pallidum` = "#666666",
    `R Caudate` = "#666666"
)

circos.clear()
pdf(here("outputs", "figures", "chord_diagram_plot.pdf"), width = 5, height = 5)
chordDiagram(
    x = rsfc_res,
    grid.col = node_colors,
    col = rsfc_res$col,
    transparency = 0.3,
    annotationTrack = "grid",
    preAllocateTracks = 1
)

# Add inner labels
circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
        circos.text(
            x = CELL_META$xcenter,
            y = CELL_META$ylim[1] - 0.25,
            labels = CELL_META$sector.index,
            facing = "bending.inside",
            col = "white",
            niceFacing = TRUE,
            adj = c(0.5, 0.5),
            cex = 0.8
        )
    },
    bg.border = NA,
    track.height = 0.1
)
dev.off()

highlight_df <- tibble::tibble(
    label = c("Left-Pallidum", "Right-Caudate"),
    highlight = TRUE
)

# Join with aseg and plot
aseg_highlighted_vis <- highlight_df %>%
    brain_join(aseg) %>%
    reposition_brain(position = "coronal") %>% 
    ggplot() +
    geom_sf(
        aes(fill = highlight), 
        color = "black"
    ) +
    scale_fill_manual(
        values = c("TRUE" = "#666666", "FALSE" = "gray85"), 
        na.value = "gray85"
    ) +
    theme_void() +
    theme(legend.position = "none")


ggsave(
    plot = aseg_highlighted_vis,
    filename = here("outputs", "figures", "aseg_highlighted.pdf"),
    device = cairo_pdf,
    width = 4,
    height = 4
)


# Logistic Regression ----------------------------------------------------------
rsfc_lr_res_list <- list()

for (ith in 1:(length(rsfc_labels))) {
    
    cat(sprintf("\rRSFC Running: %d / %d", ith, length(rsfc_labels)))
    
    rsfc_var <- rsfc_labels[ith]
    
    f1 <- glue::glue(
        "loneliness_fu_any ~ {rsfc_var} + loneliness_bl + demo_sex_v2 + \
        interview_age + race_ethnicity_3 + (1 | site_id_l) + \
        (1 | rel_family_id)"
    ) |> 
        as.formula()
    
    res <- master_rsfc_df |> 
        mutate(!!rsfc_var := as.numeric(scale(.data[[rsfc_var]]))) |> 
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
            variable_name = rsfc_var,
            .before = "p"
        )
    
    rsfc_lr_res_list[[ith]] <- res
}


rsfc_res_df <- tibble(
    variable_name = rsfc_lda_result$weights$variable_name,
    lr_weights = rsfc_lr_res_list |> reduce(bind_rows) |> pull(cohend),
    lda_weights = rsfc_lda_result$weights$mean_weight
)

cor_results <- correlation::correlation(rsfc_res_df)
r_val <- cor_results$r[1]
p_val <- cor_results$p[1]
label_text <- sprintf("italic(r) == %.2f*','~~italic(p) == %.2g", r_val, p_val)

rsfc_comp_fig <- rsfc_res_df |> 
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
        title = "Resting-state Functional Connectivity",
        x = "Mean LDA Weights (Bootstrapped, n = 1,000)",
        y = "Cohen's d from Mixed-effects Logistic Regression"
    ) +
    ggthemes::theme_pander() +
    theme(
        plot.margin = margin(5, 5, 5, 5, "mm"),
        plot.title.position = "plot"
    )
rsfc_comp_fig

ggsave(
    plot = rsfc_comp_fig,
    filename = here("outputs", "figures", "rsfc_comp_fig.pdf"),
    width = 6,
    height = 5
)
