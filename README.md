# README

**Author:** Marta Mosna
**Date:** 22.04.2026

---
  
## Files Overview
  
| File | Type |
|------|------|
| `main_script.R` | R Script |
| `hyena_reproduction_dataset.csv` | Dataset |
| `maternal_lineage_with_ranks.csv` | Dataset |
| `paternal_lineage_with_ranks.csv` | Dataset |
| `res_rep_1000_app.rds` | Data saved as R Data Object |
| `LICENSE.md` | License file |

## R Scripts

### `main_script.R`

The main script to be used for running all the analyses.

  
## Data
  
### `hyena_reproduction_dataset.csv`
  
This dataset contains life-history, reproductive, social rank, and genetic contribution data for all individuals included in the study. All reproductive counts and genetic contribution metrics are computed up to January 1st, 2023.

#### Variable Description

| Variable | Description |
|----------|-------------|
| `ID` | Unique individual identifier. |
| `birth_date` | Date of birth of the individual. |
| `death_date` | Date of death of the individual (`NA` if still alive at the end of the study period). |
| `selection_first` | Date when the individual begins its reproductive career. |
| `sex` | Biological sex of the individual. |
| `tenure` | Duration (in years) during which the individual was reproductively active within the study population (from selection_first to death or end of the study period). |
| `n_off_type` | Number of offspring genetically typed for both mother and father. |
| `n_g_off` | Number of grandoffspring (offspring of genetically typed offspring). |
| `n_gg_off` | Number of great-grandoffspring (grand-offspring of genetically typed offspring). |
| `n_alive_genotyped_offspring` | Number of genetically typed offspring alive at end of the study period. |
| `rank_maternal.std` | Standardized maternal social rank at the time of the individual's birth. |
| `n_off_obs` | Total number of offspring observed, regardless of genetic typing. |
| `n_g_off_obs` | Total number of grandoffspring observed (includes those from non-genotyped offspring). |
| `gen_cont` | Long-term genetic contribution (gᵢ): expected number of gene copies present in the population at the end of the study period. |

### `maternal_lineage_with_ranks.csv` & `paternal_lineage_with_ranks.csv`

These datasets summarize the genealogical lineages and parental ranks of spotted hyenas used in the study. Each table contains information on all individuals with known ancestry in the population, including the rank of their ancestors across multiple generations.

- `maternal_lineage_with_ranks.csv` tracks maternal ancestry (mother, maternal grandmother, etc.).

- `paternal_lineage_with_ranks.csv` tracks paternal ancestry (father, paternal grandfather, etc.).


These tables form the basis for analyses of rank inheritance and sex-biased lineage effects across generations.

#### Variable Description

| Variable | Description |
|----------|-------------|
| `ID` | Unique individual identifier. |
| `n_generation` | Number of generations including the focal individual (G0) known for the ID. |
| `rank_G0` | Maternal rank of the focal individual. |
| `rank_G1`, `rank_G2`, ... | Maternal rank of the individual's ancestors in generations G1, G2, ... respectively. |

### `res_rep_1000_app.rds`
  
This file contains results from 1,000 permutation simulations, used to test the effect of disrupting maternal lineage continuity while preserving the empirical variation in reproductive output. It allows testing whether social rank inheritance significantly affects the overrepresentation of alpha-descendant lineages across generations.

#### Structure

A list of 1,000 elements, each corresponding to a single permutation simulation. Each element is itself a list with two components:
  
  - **`geometric_mean`** — Geometric mean of fold-changes in the proportion of alpha-ranked mothers across consecutive generations (same number as with the real data).

  - **`long_data_proportions`** — A tidy data frame containing:
  
| Column | Description |
|--------|-------------|
| `generation` | Generation relative to the focal individual (G0, G-1, ...). |
| `rank` | Maternal rank of each ancestor in that generation. |
| `n` | Number of individuals with that rank in the permutation. |
| `proportion` | Proportion of individuals with that rank within the generation. |
| `total_n` | Total number of individuals in the generation. |


## Software

| Package | Version |
|---------|---------|
| R | 4.5.0 (2025-04-11 ucrt) |
| ggplot2 | 4.0.3 |
| patchwork | 1.3.2 |
| ragg | 1.5.2 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| lubridate | 1.9.5 |
| tibble | 3.3.1 |
| SkewCalc | 1.0 (install via `remotes::install_github("Ctross/SkewCalc")`) |
| DescTools | 0.99.60 |
| flextable | 0.9.11 |
| glmmTMB | 1.1.14 |
| DHARMa | 0.4.7 |
| car | 3.1.5 |
| broom.mixed | 0.2.9.7 |
| diffcor | 0.8.4 |
| officer | 0.7.3 |
| ggbreak | 0.1.7 |