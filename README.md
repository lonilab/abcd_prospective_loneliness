# What makes a lonely child: Environmental, health, and multimodal neuroimaging correlates of prospective loneliness in the ABCD study

## Overview
This project investigates multilevel correlates of **prospective loneliness** in late childhood using data from the Adolescent Brain Cognitive Development (ABCD) Study. The aim is to identify environmental, health, and neurobiological features associated with later loneliness, with emphasis on **pattern-level inference rather than prediction**.

Prospective loneliness is defined as the **presence of parent-reported loneliness at any follow-up (1–3 years)**, based on the CBCL item “complains of loneliness,” dichotomized and aggregated across waves.

## Study Design
- **Cohort**: ABCD Study (Release 5.1)
- **Baseline age**: 9–10 years  
- **Sample size**: ~9,600 participants (after exclusions)
- **Follow-up window**: 1–3 years  
- **Outcome**: Binary indicator of loneliness at any follow-up

## Data Domains
We examine a broad set of baseline features:

- **Environmental exposures (~300+)**  
  Family, school, neighborhood, and socioeconomic factors

- **Health indicators (~60+)**  
  Physical and mental health measures

- **Neuroimaging features (~500+)**
  - Gray matter volume (GMV)
  - White matter microstructure (WMM)
  - Resting-state functional connectivity (RSFC)

## Analytical Strategy

### 1. Mass-Univariate Screening
- Linear mixed-effects models (LMM)
- Random intercepts: site, family  
- Covariates: demographics, baseline loneliness
- Multiple testing: Bonferroni correction (stringent control)

Focus: **effect size patterns**, not isolated significance

### 2. Multivariate Neuroimaging Analysis
- Linear Discriminant Analysis (LDA)
- Applied separately to GMV, WMM, RSFC
- Handles correlated features and identifies **distributed contribution patterns**

### 3. Bootstrap Inference
- 1,000 iterations
- Stability of feature weights assessed across resamples

### 4. Sensitivity Analyses
- Early follow-up restriction (temporal stability)
- Baseline loneliness stratification (incident patterns)
- Pre-pandemic restriction (COVID robustness)

## Key Principles
- **Associational, not causal**: longitudinal design supports temporal ordering but not causal inference  
- **Inferential over predictive**: models are used to characterize structure in the data, not optimize classification  
- **Multilevel integration**: environmental, health, and brain features are examined within the same framework  

## Reproducibility
- Analyses implemented in **R**
- Fully scripted workflow for transparency
- Random seeds fixed for reproducibility where applicable

## Limitations
- Loneliness measured using a **single parent-reported CBCL item**
- Dichotomization reduces sensitivity to severity
- Aggregation across waves does not distinguish **transient vs. persistent** loneliness
- Observational design limits causal interpretation

## Citation
If you use this code or build upon this work, please cite:

> *What makes a lonely child: Environmental, health, and multimodal neuroimaging correlates of prospective loneliness in the ABCD study*

## Data Access
ABCD data are available via the **NIMH Data Archive (NDA)**.  
Access requires a data use agreement: https://nda.nih.gov/abcd

## Contact
For questions or collaboration:
- Laboratory of Neurodevelopmental Informatics (LoNI)
- Email:

