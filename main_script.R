
# General info ----
#
# R Script: The legacy of privilege
# 
# Purpose:
#   This script reproduces the analyses and figures reported in the
#   manuscript entitled "The Legacy of Privilege: Social Inheritance
#   Reverses Sex Differences in Reproductive Inequality in the Spotted
#   Hyena", by Mosna et al. 2026
# Author: Marta Mosna
# Date: 2026-04-22
# License: CC-BY version 4
# Contact: marta.mosna@evobio.eu 



# Load libraries ---------------------------------------------------------------

library(ggplot2)   # for plotting
library(patchwork) # for assembling plots
library(ragg)      # not called directly, but having it installed ensures consistent figure sizes
library(dplyr)     # for data wrangling
library(tidyr)     # for data wrangling
library(lubridate) # for manipulating dates
library(tibble)    # for manipulating tables
library(SkewCalc)  # for M-index computation, install using remotes::install_github("Ctross/SkewCalc")
library(DescTools) # for Gini computation
library(flextable) # for creating Word tables
library(glmmTMB)   # for statistical modeling
library(DHARMa)    # for checking modeling assumption
library(car)       # for Anova
library(broom.mixed) # for extracting model values
library(diffcor)   # to compare correlation
library(officer)   # for creating Word tables
library(ggbreak) # for fig s1-s3

# Shared theme — for ALL FIGURES adjusted to Science advances ----

theme_sadv <- function(base_size = 8, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        size  = 9,
        face  = "bold",
        hjust = 0,
        margin = margin(b = 6, unit = "pt")
      ),
      axis.title        = element_text(size = 8),
      axis.text         = element_text(size = 8),
      axis.line         = element_line(linewidth = 0.4),
      axis.ticks        = element_line(linewidth = 0.4),
      panel.border      = element_blank(),
      panel.grid        = element_blank(),
      legend.text       = element_text(size = 8),
      legend.title      = element_text(size = 8),
      plot.margin       = margin(5, 5, 5, 5, "pt"),
    )
}

sex_cols <- c(female = "#8624F5", male = "#1FC3AA")


# Create folders to store elements ---------------------------------------------
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("tables")) dir.create("tables")

# Common data loading and preparation ------------------------------------------

## Load data ----
ids_completed                <-  read.csv("data/hyena_reproduction_dataset.csv",  header = TRUE)
table_ancestor_native_r      <-  read.csv("data/maternal_lineage_with_ranks.csv", header = TRUE)
table_ancestor_native_male_r <-  read.csv("data/paternal_lineage_with_ranks.csv", header = TRUE)

data_perm <- readRDS("data/res_rep_1000_app.rds")


## Convert date columns to Date class ----
ids_completed <- ids_completed |> 
  mutate(
    birth_date      = dmy(birth_date),
    death_date      = dmy(death_date),
    selection_first = dmy(selection_first)
  )

## Define end date ----
end_genotype_date <- as.Date("2023-01-01")  

## Cohort classification ----
# Extract the year individuals reached adulthood from
# selection_first, then assign each individual to one
# of four cohort blocks: 1996-1999, 2000-2003,
# 2004-2007, and 2008-2010.

ids_completed <- ids_completed |>
  mutate(cohort = year(selection_first))
cohort_breaks <- c(1996, 2000, 2004, 2008, 2011)
cohort_labels <- c("1996–1999", "2000–2003", "2004–2007", "2008–2010")

ids_completed <- ids_completed |>
  mutate(
    cohort_c = cut(
      cohort,
      breaks = cohort_breaks,
      labels = cohort_labels,
      right  = FALSE,
      include.lowest = TRUE
    )
  )

##  Calculate completeness ----
# Calculate completeness (how complete is the count of grandoffspring)
ids_completed <- ids_completed |>
  mutate(
    completeness = dplyr::if_else(
      n_off_type > 0,
      1 - (n_alive_genotyped_offspring / n_off_type),
      1  # If no offspring, completeness is 1 — no possible grand offspring is missing
    )
  )

# Calculate average completeness per selection year
# Identify the last cohort for which great-grandoffspring counts are reliable:
# keep only cohorts with avg completeness >= 0.95, then take the most recent one.

ids_completeness <- ids_completed |> 
  group_by(cohort) |> 
  summarise(avg_completeness = mean(completeness, na.rm = TRUE)) 

last_cohort_GGO <- ids_completeness |>
  arrange(cohort) |>
  # flag where completeness drops below threshold
  mutate(meets_threshold = avg_completeness >= 0.95,
         # cumsum trick: once a FALSE appears, all subsequent rows get a higher cumsum
         consecutive = cumsum(!meets_threshold) == 0) |>
  filter(consecutive) |>
  slice_max(cohort, n = 1) |>
  pull(cohort)


# Descriptive stats ------------------------------------------------------------
total_n <- nrow(ids_completed) # Total number of individuals
total_n # 492
table_by_sex <- table(ids_completed$sex) # Number by sex
table_by_sex
# female   male 
#    287    205

## Number and percentage dead by the cutoff date ----
n_dead <- sum(ids_completed$death_date < end_genotype_date, na.rm = TRUE)
n_dead # 487
perc_dead <- round(n_dead / total_n * 100, 2)
perc_dead # 98.98


# RESULT: Reproductive success in female and male spotted hyenas ---------------

## 1. Annual offspring rate ----

# Filter individuals with tenure > 1 year
repro_data <- ids_completed |> 
  filter(tenure > 1) |> 
  mutate(offspring_rate = n_off_type / tenure)

# Summary stats by sex
summary_offspring_rate <- repro_data |> 
  group_by(sex) |> 
  summarise(
    n = n(),
    mean_rate = mean(offspring_rate),
    min_rate = min(offspring_rate),
    max_rate = max(offspring_rate),
    sd_rate = sd(offspring_rate)
  )
summary_offspring_rate
# # A tibble: 2 × 6
#   sex        n mean_rate min_rate max_rate sd_rate
#   <chr>  <int>     <dbl>    <dbl>    <dbl>   <dbl>
# 1 female   268     0.463        0     1.82   0.334
# 2 male     188     0.745        0     2.83   0.609

# Is the difference in rate production between sexes significant ? 
wilcox.test(offspring_rate ~ sex, data = repro_data)
# W = 18695, p-value = 2.622e-06


## 2. LRS number of offspring ----

summary_LRS <- ids_completed |> 
  group_by(sex) |> 
  summarise(
    n = n(),
    mean_rate = mean(n_off_type),
    min_rate = min(n_off_type),
    max_rate = max(n_off_type),
    sd_rate = sd(n_off_type)
  )
summary_LRS
# # A tibble: 2 × 6
#   sex        n mean_rate min_rate max_rate sd_rate
#   <chr>  <int>     <dbl>    <int>    <int>   <dbl>
# 1 female   287      3.80        0       19    3.71
# 2 male     205      4.95        0       26    4.73


# Is the difference in LRS between sexes significant? 
wilcox.test(n_off_type ~ sex, data = ids_completed)
# W = 25816, p-value = 0.01967

# Percentage of ids with no offspring during their lifetime
non_repro_stats <- ids_completed |> 
  group_by(sex) |> 
  summarise(
    n_total = n(),
    n_zero = sum(n_off_type == 0),
    perc_zero = round((n_zero / n_total) * 100, 1)
  )
non_repro_stats
# # A tibble: 2 × 4
#   sex    n_total n_zero perc_zero
#   <chr>    <int>  <int>     <dbl>
# 1 female     287     64      22.3
# 2 male       205     38      18.5


## 3.Grandoffspring count ----

summary_GOff <- ids_completed |> 
  group_by(sex) |> 
  summarise(
    n = n(),
    mean_rate = mean(n_g_off),
    min_rate = min(n_g_off),
    max_rate = max(n_g_off),
    sd_rate = sd(n_g_off)
  )
summary_GOff
# # A tibble: 2 × 6
#   sex        n mean_rate min_rate max_rate sd_rate
#   <chr>  <int>     <dbl>    <int>    <int>   <dbl>
# 1 female   287      8.23        0      106    15.6
# 2 male     205     11.3         0       85    15.8


# Is the difference in grandoffspring count between sexes significant? 
wilcox.test(n_g_off ~ sex, data = ids_completed)
# W = 24236, p-value = 0.000543

# Percentage of ids with no grandoffspring 
non_repro_GOff <- ids_completed |> 
  group_by(sex) |> 
  summarise(
    n_total = n(),
    n_zero = sum(n_g_off == 0),
    perc_zero = round((n_zero / n_total) * 100, 1)
  )
non_repro_GOff
# # A tibble: 2 × 4
#   sex    n_total n_zero perc_zero
#   <chr>    <int>  <int>     <dbl>
# 1 female     287    135      47  
# 2 male       205     69      33.7


## 4. Great-grandoffspring count ----

ids_completed_GGO <- ids_completed |>
  filter(cohort <= last_cohort_GGO)   

summary_by_sex_GGD <- ids_completed_GGO |>
  group_by(sex) |> 
  summarise(
    n = n(),
    mean_rate = mean(n_gg_off),
    min_rate = min(n_gg_off),
    max_rate = max(n_gg_off),
    sd_rate = sd(n_gg_off)
  )
summary_by_sex_GGD
# A tibble: 2 × 6
##sex        n mean_rate min_rate max_rate sd_rate
#  <chr>  <int>     <dbl>    <int>    <int>   <dbl>
# 1 female   168      11.7        0      152    26.0
# 2 male     123      16.5        0      155    27.8


# Is the difference in great-grandoffspring count between sexes significant? 
wilcox.test(n_gg_off ~ sex, data = ids_completed_GGO)
# W = 8690, p-value = 0.01298

# Percentage of ids with no great-grandoffspring 
non_repro_GGO <- ids_completed_GGO |> 
  group_by(sex) |> 
  summarise(
    n_total = n(),
    n_zero = sum(n_gg_off == 0),
    perc_zero = round((n_zero / n_total) * 100, 1)
  )
non_repro_GGO
# # A tibble: 2 × 4
#   sex    n_total n_zero perc_zero
#   <chr>    <int>  <int>     <dbl>
# 1 female     168    95        56.5
# 2 male       123    53        43.1


# RESULT:Reproductive inequality increases more in females across generations ----

## Helper functions ----

parse_stan_summary <- function(lines) {
  
  i_start <- grep("^\\s*M\\s", lines)[1]
  i_end   <- grep("^Samples were drawn", lines)[1] - 1
  
  block <- trimws(lines[i_start:i_end])
  block <- gsub("\\s+", " ", block)
  
  df <- read.table(text = block, header = FALSE, stringsAsFactors = FALSE)
  
  names(df) <- c(
    "param", "mean", "se_mean", "sd",
    "p2.5", "p25", "p50", "p75", "p97.5",
    "n_eff", "Rhat"
  )
  
  df
}


fit_M_index <- function(RS, Time, ..., priors = NULL) {
  
  stan_out <- capture.output(
    M_index_stan(RS, Time, ..., priors = priors)
  )
  
  tab <- parse_stan_summary(stan_out)
  M_row <- tab[tab$param == "M", ]
  
  list(
    mean = M_row$mean,
    ci   = c(M_row$`p2.5`, M_row$`p97.5`)
  )
}


## Data preparation ----

female_moy <- ids_completed |>
  filter(sex == "female",
         tenure > 0.1) # to avoid numerical issues

male_moy <- ids_completed |>
  filter(sex == "male")

female_life <- ids_completed |>
  filter(sex == "female") |>
  mutate(life = 1) # to make sure all individuals have the same `Time` in M-index

male_life <- ids_completed |>
  filter(sex == "male") |>
  mutate(life = 1)

female_life_GGO <- ids_completed_GGO |>
  filter(sex == "female") |>
  mutate(life = 1) # to make sure all individuals have the same `Time` in M-index

male_life_GGO <- ids_completed_GGO |>
  filter(sex == "male") |>
  mutate(life = 1)

## 1. Moy: annual number of offspring ----

moy_f <- fit_M_index(
  RS   = female_moy$n_off_type,
  Time = female_moy$tenure
)

moy_m <- fit_M_index(
  RS   = male_moy$n_off_type,
  Time = male_moy$tenure
)

(Moy_f_mean <- moy_f$mean) # 0.19
(Moy_f_ci   <- moy_f$ci) # 0.07 0.31 # note: results made differ slightly since seed not fixed in stan

(Moy_m_mean <- moy_m$mean) # 0.41
(Moy_m_ci   <- moy_m$ci)# 0.26 0.57 # note: results made differ slightly since seed not fixed in stan


## 2. MoL: offspring per lifetime ----

mol_f <- fit_M_index(
  RS   = female_life$n_off_type,
  Time = female_life$life
)

mol_m <- fit_M_index(
  RS   = male_life$n_off_type,
  Time = male_life$life
)

(Mol_f_mean <- mol_f$mean) # 0.76
(Mol_f_ci   <- mol_f$ci) # 0.58 0.95 # note: results made differ slightly since seed not fixed in stan

(Mol_m_mean <- mol_m$mean) # 0.78
(Mol_m_ci   <- mol_m$ci) # 0.6 1.0 # note: results made differ slightly since seed not fixed in stan 


## 3. MG: grand-offspring ----

mg_f <- fit_M_index(
  RS   = female_life$n_g_off,
  Time = female_life$life
)

mg_m <- fit_M_index(
  RS   = male_life$n_g_off,
  Time = male_life$life
)

(Mg_f_mean <- mg_f$mean) # 3.53
(Mg_f_ci   <- mg_f$ci) # 3.14 3.96

(Mg_m_mean <- mg_m$mean) # 1.91
(Mg_m_ci   <- mg_m$ci) # 1.70 2.15


## 4. MGG: great-grand-offspring ----

priors_mgg <- matrix(
  c(1, 0.1, 8, 0.2),
  nrow = 2,
  byrow = TRUE
)

mgg_f <- fit_M_index(
  RS   = female_life_GGO$n_gg_off,
  Time = female_life_GGO$life,
  adapt_delta = 0.99,
  max_treedepth = 14,
  priors = priors_mgg
)

mgg_m <- fit_M_index(
  RS   = male_life_GGO$n_gg_off,
  Time = male_life_GGO$life,
  adapt_delta = 0.99,
  max_treedepth = 14,
  priors = priors_mgg
)

(Mgg_f_mean <- mgg_f$mean) # 4.63
(Mgg_f_ci   <- mgg_f$ci) # 4.16 5.16

(Mgg_m_mean <- mgg_m$mean) # 2.62
(Mgg_m_ci   <- mgg_m$ci) # 2.34 2.92

## Preparation for plotting ----

# Final data frame for plotting

data <- data.frame(
  Label = c("MOYear", "MOTot", "MG", "MGG"),
  Male_Reproductive_Skew   = c(Moy_m_mean,  Mol_m_mean,  Mg_m_mean,  Mgg_m_mean),
  Male_Lower               = c(Moy_m_ci[1], Mol_m_ci[1], Mg_m_ci[1], Mgg_m_ci[1]),
  Male_Upper               = c(Moy_m_ci[2], Mol_m_ci[2], Mg_m_ci[2], Mgg_m_ci[2]),
  Female_Reproductive_Skew = c(Moy_f_mean,  Mol_f_mean,  Mg_f_mean,  Mgg_f_mean),
  Female_Lower             = c(Moy_f_ci[1], Mol_f_ci[1], Mg_f_ci[1], Mgg_f_ci[1]),
  Female_Upper             = c(Moy_f_ci[2], Mol_f_ci[2], Mg_f_ci[2], Mgg_f_ci[2])
)

# Axes limits (shared for 1:1 comparison)

axis_max <- max(
  data$Male_Upper,
  data$Female_Upper
  ) * 1.05

# M-index comparison plot (males vs females)

M_index_plot <- ggplot(
  data,
  aes(x = Male_Reproductive_Skew, y = Female_Reproductive_Skew)
) +
  # Vertical error bars (females)
  geom_errorbar(
    aes(ymin = Female_Lower, ymax = Female_Upper),
    width = 0.02
  ) +
  # Horizontal error bars (males)
  geom_errorbar(
    aes(xmin = Male_Lower, xmax = Male_Upper),
    width = 0.02,
    orientation = "y"
  ) +
  # Point estimates
  geom_point(
    shape = 21,
    fill = "grey",
    color = "black",
    size = 1.5
  ) +
  # Labels placed to the right of male CI
  geom_text(
    aes(x = Male_Upper, label = Label),
    hjust = -0.2,
    vjust = 0.4,
    size = 2.5
  ) +
  # 1:1 reference line
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "grey",
    linetype = "dashed"
  ) +
  # Shared axis scaling
  scale_x_continuous(
    limits = c(0, axis_max),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0, axis_max),
    expand = expansion(mult = 0)
  ) +
  # Labels
  labs(
    x = "M-index of reproductive inequality for males",
    y = "M-index of reproductive inequality\n for females"
  ) +
  # Theme
  theme_sadv()


## Compute Lorenz curves ----

compute_lorenz <- function(data,  sex = c("female", "male")) {
  data |>
    filter(sex == {{sex}}) |>
    arrange(n_off_type) |>     # Offspring
    mutate(
      L_off = 100 * cumsum(n_off_type) / sum(n_off_type),
      p_off = 100 * seq_len(n()) / n()
    ) |> 
    arrange(n_g_off) |>     # Grandoffspring
    mutate(
      L_GOff = 100 * cumsum(n_g_off) / sum(n_g_off),
      p_GOff = 100 * seq_len(n()) / n()
    )
}

compute_lorenz_ggo <- function(data,  sex = c("female", "male")) {
  data |>
    filter(sex ==  {{sex}}) |>
    arrange(n_gg_off) |>   # Great-grandoffspring
    mutate(
      L_gGOff = 100 * cumsum(n_gg_off) / sum(n_gg_off),
      p_gGOff = 100 * seq_len(n()) / n()
    ) |>
    dplyr::select(ID, sex, L_gGOff, p_gGOff)  # only keep what's needed for joining
}


lorenz_data <- bind_rows(
  compute_lorenz(ids_completed, sex = "female"),
  compute_lorenz(ids_completed, sex = "male")
) |>
  left_join(
    bind_rows(
      compute_lorenz_ggo(ids_completed_GGO, sex = "female"),
      compute_lorenz_ggo(ids_completed_GGO, sex = "male")
    ),
    by = c("ID", "sex")
  )


## Plot Lorenz curves ----

descendants_plot <- ggplot(lorenz_data) +
  geom_line(
    aes(x = p_off, y = L_off, colour = sex, linetype = "offspring"),
    linewidth = 1.5, alpha = 0.5
  ) +
  geom_line(
    aes(x = p_GOff, y = L_GOff, colour = sex, linetype = "gd offspring"),
    linewidth = 1.5, alpha = 0.5
  ) +
  geom_line(
    aes(x = p_gGOff, y = L_gGOff, colour = sex, linetype = "gd gd offspring"),
    linewidth = 1.5, alpha = 0.5
  ) +
  scale_linetype_manual(
    values = c(
      "offspring" = 1,
      "gd offspring" = 5,
      "gd gd offspring" = 3
    ),
    breaks = c("offspring", "gd offspring", "gd gd offspring"), # <- ensures correct order in legend
    labels = c(
      "offspring" = "Offspring",
      "gd offspring" = "Grandoffspring",
      "gd gd offspring" = "Great-grandoffspring"
    )
  ) +
  # remaining plot settings stay the same
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  scale_x_continuous(breaks = seq(0, 100, 10), minor_breaks = seq(0, 100, 5)) +
  scale_y_continuous(breaks = seq(0, 100, 10), minor_breaks = seq(0, 100, 5)) +
  scale_color_manual(values = c("female" = "#8624F5", "male" = "#1FC3AA")) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.7))) +
  coord_fixed(xlim = c(0, 100.1), ylim = c(0, 100), expand = FALSE) +
  labs(
    x = "% of individuals by contribution to descendants",
    y = "Cumulative share of descendants (%)",
    colour = "Sex",
    linetype = "Descendants"
  ) +
  theme_sadv() +
  theme(
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.2),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )


##  Figure 1. Comparison of reproductive inequality in female and male spotted hyenas ----

figure_1 <- M_index_plot + descendants_plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

ggsave(
  filename = "figures/figure1.pdf",
  plot     = figure_1,
  width    = 184,
  height   = 90,
  units    = "mm"
)

ggsave(
  filename = "figures/figure1.png",
  plot     = figure_1,
  width    = 184,
  height   = 90,
  units    = "mm",
  dpi = 600
)



## Compute Lorenz thresholds (50%) per sex ----

get_lorenz_threshold50 <- function(x_prop, y_cumul) {
  approx( # to interpolate
    x = unique(y_cumul),
    y = unique(x_prop)[!duplicated(y_cumul)],
    xout = 50
  )$y
}

# Split data by sex
lorenz_female <- lorenz_data[lorenz_data$sex == "female", ]
lorenz_male   <- lorenz_data[lorenz_data$sex == "male", ]

# Female top 50% contributors
top50_offspring_f  <- 100 - get_lorenz_threshold50(lorenz_female$p_off, lorenz_female$L_off)
top50_grandoff_f   <- 100 - get_lorenz_threshold50(lorenz_female$p_GOff, lorenz_female$L_GOff)
top50_greatf_f     <- 100 - get_lorenz_threshold50(
  lorenz_female$p_gGOff[!is.na(lorenz_female$p_gGOff)],
  lorenz_female$L_gGOff[!is.na(lorenz_female$L_gGOff)]
)

# Male top 50% contributors
top50_offspring_m <- 100 - get_lorenz_threshold50(lorenz_male$p_off, lorenz_male$L_off)
top50_grandoff_m  <- 100 - get_lorenz_threshold50(lorenz_male$p_GOff, lorenz_male$L_GOff)
top50_greatf_m    <- 100 - get_lorenz_threshold50(
  lorenz_male$p_gGOff[!is.na(lorenz_male$p_gGOff)],
  lorenz_male$L_gGOff[!is.na(lorenz_male$L_gGOff)]
)

lorenz_top50 <- data.frame( # Combine into a simple, readable table
  Sex = c("Female", "Male"),
  Offspring      = c(top50_offspring_f, top50_offspring_m),
  Grandoffspring = c(top50_grandoff_f, top50_grandoff_m),
  GreatGrandoff  = c(top50_greatf_f, top50_greatf_m)
)

lorenz_top50 <- lorenz_top50 |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))
lorenz_top50
#      Sex Offspring Grandoffspring GreatGrandoff
# 1 Female     19.49           8.22          6.01
# 2   Male     19.76          12.37          9.44



## Table 1. M-index and Gini coefficients for reproductive success in female and male spotted hyenas ----

## Helper function to get Gini + bootstrapped CI ----
get_gini_ci <- function(x) {
  result <- Gini(x, conf.level = 0.95, R = 1000, type = "bca") # Davison & Hinkley 197 book p 203
  print(result)
  paste0(round(result["gini"], 2),
         " (", round(result["lwr.ci"], 2), "–", round(result["upr.ci"], 2), ")")
}

## 1. Extract Gini coefficients with bootstrapped CI ----
gini_wide <- bind_rows(
  repro_data |>
    group_by(sex) |>
    summarise(
      value = get_gini_ci(offspring_rate),
      Measure = "Offspring annual rate"
    ),
  ids_completed |>
    group_by(sex) |>
    summarise(
      value = get_gini_ci(n_off_type),
      Measure = "Offspring total",
    ),
  ids_completed |>
    group_by(sex) |>
    summarise(
      value = get_gini_ci(n_g_off),
      Measure = "Grandoffspring",
    ),
  ids_completed_GGO |>
    group_by(sex) |>
    summarise(
      value = get_gini_ci(n_gg_off),
      Measure = "Great-grandoffspring",
    )
) |>
  pivot_wider(names_from = sex, values_from = value)

## 2. Format M-index data frame ----
m_index_table <- data.frame(
  Measure = c(
    "Offspring annual rate",
    "Offspring total",
    "Grandoffspring",
    "Great-grandoffspring"
  ),
  Female_M = paste0(
    round(data$Female_Reproductive_Skew, 2),
    " (", round(data$Female_Lower, 2),
    "–", round(data$Female_Upper, 2), ")"
  ),
  Male_M = paste0(
    round(data$Male_Reproductive_Skew, 2),
    " (", round(data$Male_Lower, 2),
    "–", round(data$Male_Upper, 2), ")"
  )
)

## 3. Join M-index and Gini ----
combined_table_skew <- left_join(m_index_table, gini_wide, by = "Measure") |>
  rename(
    `M-index female\n(95% CrI)*` = Female_M,
    `M-index male\n(95% CrI)*`   = Male_M,
    `Gini female\n(95% CI)†`    = female,
    `Gini male\n(95% CI)†`      = male
  )

## 4. Reformat for Word table export ----
combined_table_skew <- combined_table_skew |>
  mutate(across(2:5, ~ gsub("\\(", "\n(", .x)))

# Little preview in RStudio
ft_combined <- flextable(combined_table_skew) |>
  theme_booktabs() |>
  align(align = "center", j = 2:5, part = "body") |>
  bold(part = "header") |>
  autofit()

# A4 page with 1 inch margins = 6.27 inches content width
page_width <- 6.27

ft_combined <- flextable(combined_table_skew) |>
  theme_booktabs() |>
  align(align = "center", j = 2:5, part = "body") |>
  align(align = "left", j = 1, part = "body") |>
  bold(part = "header") |>
  font(fontname = "Times New Roman", part = "all") |>      
  fontsize(size = 12, part = "all") |>           
  width(j = 1, width = 1.2) |>                 
  width(j = 2:5, width = (page_width - 1.2) / 4) |>  
  hrule(rule = "auto") |>    
  set_table_properties(layout = "fixed")         

doc <- read_docx() |>
  body_add_flextable(ft_combined)
 
print(doc, target = "tables/Table_1.docx")



## Figure S4:M-index of reproductive inequality by cohort block ----
## M-index of reproductive inequality by cohort block
## Each cohort block contains only individuals who started the reproductive carrer
## in the same time window, ensuring ecological and demographic
## contemporaneity.

# Loop over cohort blocks and fit M-index models
# for each sex and generational level

for (coh in cohort_labels) {
  
  print(paste0("computing M-indices for the cohort ", coh, ". Be patient..."))
  
  female_moy_c  <- female_moy  |> filter(cohort_c == coh)
  male_moy_c    <- male_moy    |> filter(cohort_c == coh)
  female_life_c <- female_life |> filter(cohort_c == coh)
  male_life_c   <- male_life   |> filter(cohort_c == coh)
  female_life_GGO_c <- female_life_GGO |> filter(cohort_c == coh)
  male_life_GGO_c   <- male_life_GGO   |> filter(cohort_c == coh)
  
  # 1. MOYear
  moy_f <- fit_M_index(RS = unname(female_moy_c$n_off_type),  Time = unname(female_moy_c$tenure))
  moy_m <- fit_M_index(RS = unname(male_moy_c$n_off_type),    Time = unname(male_moy_c$tenure))
  
  # 2. MOTot
  mol_f <- fit_M_index(RS = unname(female_life_c$n_off_type), Time = unname(female_life_c$life))
  mol_m <- fit_M_index(RS = unname(male_life_c$n_off_type),   Time = unname(male_life_c$life))
  
  # 3. MG
  mg_f  <- fit_M_index(RS = unname(female_life_c$n_g_off),    Time = unname(female_life_c$life))
  mg_m  <- fit_M_index(RS = unname(male_life_c$n_g_off),      Time = unname(male_life_c$life))
  
  # 4. MGG — only fit if data is available for both sexes
  if (nrow(female_life_GGO_c) > 0 && nrow(male_life_GGO_c) > 0) {
    mgg_f <- fit_M_index(RS = unname(female_life_GGO_c$n_gg_off), Time = unname(female_life_GGO_c$life),
                         adapt_delta = 0.99, max_treedepth = 14, priors = priors_mgg)
    mgg_m <- fit_M_index(RS = unname(male_life_GGO_c$n_gg_off),   Time = unname(male_life_GGO_c$life),
                         adapt_delta = 0.99, max_treedepth = 14, priors = priors_mgg)
  } else {
    mgg_f <- list(mean = NA, ci = c(NA, NA))
    mgg_m <- list(mean = NA, ci = c(NA, NA))
  }
  
  # Build data frame
  data_cohort <- data.frame(
    Label                    = c("MOYear", "MOTot", "MG", "MGG"),
    Male_Reproductive_Skew   = c(moy_m$mean, mol_m$mean, mg_m$mean, mgg_m$mean),
    Male_Lower               = c(moy_m$ci[1], mol_m$ci[1], mg_m$ci[1], mgg_m$ci[1]),
    Male_Upper               = c(moy_m$ci[2], mol_m$ci[2], mg_m$ci[2], mgg_m$ci[2]),
    Female_Reproductive_Skew = c(moy_f$mean, mol_f$mean, mg_f$mean, mgg_f$mean),
    Female_Lower             = c(moy_f$ci[1], mol_f$ci[1], mg_f$ci[1], mgg_f$ci[1]),
    Female_Upper             = c(moy_f$ci[2], mol_f$ci[2], mg_f$ci[2], mgg_f$ci[2])
  )
  
  assign(paste0("data_", gsub("–", "_", coh)), data_cohort)
}

# Plotting function
# Replicates the main text M-index plot for a single
# cohort block. Axis limits are computed from the data
# and allowed to go below zero to accommodate CIs that
# cross zero (common in smaller cohort subsets).

plot_M_index <- function(data, title = NULL) {
  
  axis_max <- max(data$Male_Upper, data$Female_Upper, na.rm = TRUE) * 1.05
  axis_min <- min(data$Female_Lower, data$Male_Lower, na.rm = TRUE)
  axis_min <- min(axis_min, 0)  # keep 0 if all values are positive
  
  
  ggplot(data, aes(x = Male_Reproductive_Skew, y = Female_Reproductive_Skew)) +
    geom_errorbar(aes(ymin = Female_Lower, ymax = Female_Upper), width = 0) +
    geom_errorbar(aes(xmin = Male_Lower, xmax = Male_Upper), width = 0, orientation = "y") +
    geom_point(shape = 21, fill = "grey", color = "black", size = 1.5) +
    geom_text(aes(x = Male_Upper, label = Label), hjust = -0.2, vjust = 0.4, size = 2.5) +
    geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
    scale_x_continuous( 
      limits = c(axis_min, axis_max),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_y_continuous(
      limits = c(axis_min, axis_max),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    labs(
      title = title, 
      x     = "M-index of reproductive inequality for males",
      y     = "M-index of reproductive inequality\n for females"
    ) +
    theme_sadv() +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0)) 
}


# Generate one plot per cohort block and arrange
# in a 2x2 grid using patchwork

plot_1996_1999 <- plot_M_index(data_1996_1999, title = "A | 1996–1999")
plot_2000_2003 <- plot_M_index(data_2000_2003, title = "B | 2000–2003")
plot_2004_2007 <- plot_M_index(data_2004_2007, title = "C | 2004–2007")
plot_2008_2010 <- plot_M_index(data_2008_2010, title = "D | 2008–2010")

figure_s4 <- (plot_1996_1999 | plot_2000_2003) /
  (plot_2004_2007 | plot_2008_2010) 

ggsave(
  filename = "figures/figureS4.pdf",
  plot     = figure_s4,
  width    = 184,
  height   = 184,
  units    = "mm"
)

ggsave(
  filename = "figures/figureS4.png",
  plot     = figure_s4,
  width    = 184,
  height   = 184,
  units    = "mm",
  dpi = 600
)

# sample size for moy,mol,mg
ids_completed |>
  group_by(cohort_c, sex) |>
  summarise(n = n(), .groups = "drop") |>
  pivot_wider(names_from = sex, values_from = n)
# # A tibble: 4 × 3
#   cohort_c  female  male
#   <fct>      <int> <int>
# 1 1996–1999     41    36
# 2 2000–2003     79    50
# 3 2004–2007     80    64
# 4 2008–2010     87    55


# sample size for mgg
ids_completed_GGO |>
  group_by(cohort_c, sex) |>
  summarise(n = n(), .groups = "drop") |>
  pivot_wider(names_from = sex, values_from = n)
# # A tibble: 3 × 3
#   cohort_c  female  male
#   <fct>      <int> <int>
# 1 1996–1999     41    36
# 2 2000–2003     79    50
# 3 2004–2007     48    37

# RESULT:Stronger intergenerational transmission of reproductive success in females ----

##  Data preparation ----

# Build rank categories (3 bins) and set clean labels
ids_completed$rank_cat <- cut(
  ids_completed$rank_maternal.std,
  breaks = c(-1, -0.3333334, 0.3333333, 1),
  labels = c("Low", "Medium", "High"),
  right = FALSE,
  include.lowest = TRUE
)

ids_completed$rank_cat <- as.character(ids_completed$rank_cat)
ids_completed$rank_cat[is.na(ids_completed$rank_cat)] <- "NA"
ids_completed$rank_cat <- factor(
  ids_completed$rank_cat,
  levels = c("High", "Medium", "Low", "NA")
)

# Exclude individuals with zero offspring
ids_GOmodel_no0 <- ids_completed |>
  filter(n_off_type != 0)

ids_GOmodel_no0 |>
  count(sex, name = "n")
#      sex   n
# 1 female 223
# 2   male 167

# Transform variables
ids_GOmodel_no0 <- ids_GOmodel_no0 |>
  mutate(
    log_offspring  = log(n_off_type),
    completeness_z = as.numeric(scale(completeness)),
    sex            = factor(sex, levels = c("female", "male"))
  )


## Candidate models (as described in Methods) and model selection ----

#  Poisson 
model_pois <- glmmTMB(
  n_g_off ~  sex *(log_offspring + completeness_z),
  family = poisson,
  data = ids_GOmodel_no0,
  na.action = na.omit
)

#  Zero-inflated Poisson (ZIP)
model_zip <- glmmTMB(
  n_g_off ~  sex *(log_offspring + completeness_z),
  ziformula = ~1,  
  family = poisson,
  data = ids_GOmodel_no0,
  na.action = na.omit
)

# Negative binomial (NB2)
model_nb2 <- glmmTMB(
  n_g_off ~  sex *(log_offspring + completeness_z),
  family = nbinom2,
  data = ids_GOmodel_no0,
  na.action = na.omit
)

# Zero-inflated NB2 (intercept only)
model_zinb0 <- glmmTMB(
  n_g_off ~  sex *(log_offspring + completeness_z),
  ziformula = ~1,
  family = nbinom2, 
  data = ids_GOmodel_no0,
  na.action = na.omit
)

# Zero-inflated NB2 (intercept + completeness)
model_zinb_full <- glmmTMB(
  n_g_off ~  sex *(log_offspring + completeness_z),
  ziformula =  ~1 + completeness_z,
  family = nbinom2, 
  data = ids_GOmodel_no0,
  na.action = na.omit
)

## Model diagnostics (DHARMa)
sim_pois <- simulateResiduals(model_pois, n = 1000)
sim_zip <- simulateResiduals(model_zip, n = 1000)
sim_nb2 <- simulateResiduals(model_nb2, n = 1000)
sim_zinb0 <- simulateResiduals(model_zinb0, n = 1000)
sim_zinb_full <- simulateResiduals(model_zinb_full, n = 1000)

testDispersion(sim_pois) # yes
testDispersion(sim_zip) # yes
testDispersion(sim_nb2)  # yes
testDispersion(sim_zinb0)  # yes 0.002
testDispersion(sim_zinb_full) # no 0.07

testZeroInflation(sim_pois) #yes
testZeroInflation(sim_zip) #no
testZeroInflation(sim_nb2) #yes
testZeroInflation(sim_zinb0) #yes
testZeroInflation(sim_zinb_full) #no 0.052


anova(model_zinb0, model_zinb_full) ## Compare zero-inflation formulations


##  Final model inference (Type II LRT)
summary(model_zinb_full)
LRS_model <- Anova(model_zinb_full, type = "II") 


## Table S1. Parameter estimates from the zero-inflated negative binomial (ZINB) model of grandoffspring count ----

# Extract and clean fixed-effect coefficients from the ZINB model

model_coef <- tidy(model_zinb_full, effects = "fixed") |> 
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 2)
  ) |> 
  dplyr::select(-effect, -statistic, -p.value)

# Prepare likelihood ratio test (LRT) results from Type II Anova

anova_tbl <- LRS_model |> 
  as.data.frame() |>             
  rownames_to_column("term") |>       
  mutate(
    # Standardize term names to match coefficient table
    term = sub("^sex$", "sexmale", term),
    term = sub("^sex:", "sexmale:", term),
    Chisq = round(Chisq, 3)                    # Round Chi-squared statistic to 3 decimals
  ) |> 
  mutate(
    # Format p-values for readability
    `Pr(>Chisq)` = case_when(
      is.na(`Pr(>Chisq)`) ~ NA_character_,       # keep NA for missing
      `Pr(>Chisq)` < 0.001 ~ "< 0.001",
      TRUE                 ~ as.character(round(`Pr(>Chisq)`, 3)) # round to 3 digits
    )
  )

# Join coefficient estimates with LRT statistics
model_coef2 <- model_coef |> 
  left_join(anova_tbl, by = "term") |> 
  mutate(
    chisq  = ifelse(component == "cond", Chisq,  NA_real_),
    df     = ifelse(component == "cond", Df,     NA_real_),
    `pr(>Chisq)` = ifelse(component == "cond", `Pr(>Chisq)`, NA_real_),
    component = replace_values(component, "cond" ~ "count")) |> 
  dplyr::select(-Chisq, -Df, -`Pr(>Chisq)`)    # Remove intermediate columns

# Rename terms and columns
term_labels <- c(
  "(Intercept)" = "Intercept",
  "log_offspring" = "Log number of offspring",
  "sexmale" = "Sex (male)",
  "completeness_z" = "Completeness", # or   "completeness" = "Completeness"
  "sexmale:log_offspring" = "Sex (male) × Log number of offspring",
  "sexmale:completeness_z" = "Sex (male) × Completeness"
)

model_coef2$term <- term_labels[model_coef2$term]
colnames(model_coef2) <- c( "Component", "Term", "Estimate", "SE", "chisq","df", "p-value")

ft_coef <- flextable(model_coef2) |> 
  autofit() |> 
  set_table_properties(width = 0.9, layout = "autofit")

ft_coef <- flextable(model_coef2) |>
  theme_booktabs() |>
  bold(part = "header") |>
  font(fontname = "Times New Roman", part = "all") |>      
  fontsize(size = 12, part = "all") |>           
  autofit() |> 
  hrule(rule = "auto") |>                        
  set_table_properties(layout = "fixed")         

doc <- read_docx() |>
  body_add_flextable(ft_coef)

print(doc, target = "tables/Table_S1.docx")


## Correlation between offspring and grandoffspring for high-fertility individuals ----

# Filter dataset for individuals with >= 10 offspring
female_cor <- ids_GOmodel_no0 |> 
  filter(sex == "female", n_off_type >= 10)  # females
male_cor <- ids_GOmodel_no0 |> 
  filter(sex == "male", n_off_type >= 10)    # males

# Sample sizes
(n_female <- nrow(female_cor)) # 21
(n_male   <- nrow(male_cor))   # 34

# Spearman correlations between offspring and grandoffspring
(rho_fem <- with(female_cor, cor(x = n_off_type, y = n_g_off, method = "spearman"))) # 0.6791629
(rho_mal <- with(male_cor, cor(x = n_off_type, y = n_g_off, method = "spearman"))) # 0.3128143

# Test whether the correlation differs significantly between sexes
# (Fisher z-test using diffcor package)
diff_test <- diffcor.two(rho_fem, rho_mal,
                         n1 = nrow(female_cor),
                         n2 = nrow(male_cor), digit = 5)
diff_test[c("z", "p")]
#         z       p
# 1 1.70044 0.04452

# Proportion of high-ranking mothers among individuals with >=10 offspring
(prop_hr_female <- mean(female_cor$rank_cat == "High")) # 0.7142857
(prop_hr_male   <- mean(male_cor$rank_cat == "High"))   # 0.4117647


## Figure 2 – Relationship between number of offspring and grandoffspring ----

# Compute buffered axis limits for visualization
y_min    <- min(ids_GOmodel_no0$n_g_off)
y_max    <- max(ids_GOmodel_no0$n_g_off)
y_span   <- y_max - y_min
buffer   <- y_span * 0.05
new_min  <- y_min - buffer
new_max  <- y_max + buffer
x_limits <- range(ids_GOmodel_no0$n_off_type) + c(-1, 1)

size_vals <- c("High" = 4, "Medium" = 2.5, "Low" = 1.5, "NA" = 0.5)

# Female panel

p_f_response <- ggplot(
  subset(ids_GOmodel_no0, sex == "female" & n_off_type >= 1),
  aes(x = n_off_type, y = n_g_off, size = rank_cat)
) +
  geom_point(alpha = 0.4, color = sex_cols["female"]) +
  scale_size_manual(values = size_vals,
                    breaks = c("High", "Medium", "Low", "NA"),
                    name   = "Maternal rank") +
  scale_y_continuous(limits = c(new_min, new_max), expand = c(0, 0)) +
  scale_x_continuous(limits = x_limits, expand = expansion(add = c(1, 1))) +
  labs(y = "Number of grandoffspring") +
  theme_sadv() +
  theme(
    axis.title.x   = element_blank(),
    legend.position = "none"        # hide legend in female panel
  )

# Male panel

p_m_response <- ggplot(
  subset(ids_GOmodel_no0, sex == "male" & n_off_type >= 1),
  aes(x = n_off_type, y = n_g_off, size = rank_cat)
) +
  geom_point(alpha = 0.4, color = sex_cols["male"]) +
  scale_size_manual(values = size_vals,
                    breaks = c("High", "Medium", "Low", "NA"),
                    name   = "Maternal rank") +
  guides(
    size  = guide_legend(override.aes = list(color = "black")),
    color = "none"
  ) +
  scale_y_continuous(limits = c(new_min, new_max), expand = c(0, 0)) +
  scale_x_continuous(limits = x_limits, expand = expansion(add = c(1, 1))) +
  theme_sadv() +
  theme(
    axis.title.x         = element_blank(),
    axis.title.y         = element_blank(),
    legend.position      = c(1, 0.95),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(
      fill      = "white",
      colour    = "black",
      linewidth = 0.5
    ),
    legend.key.size  = unit(0.3, "cm")
  )

# Combine the two panels
final_response <- (p_f_response | p_m_response) +
  plot_annotation() +
  plot_layout() &
  theme(plot.margin = margin(b = 0))

final_response <- wrap_elements(final_response) +
  labs(tag = "Number of offspring") +
  theme(
    plot.tag          = element_text(size = 8, face = "plain", hjust = 0.5),
    plot.tag.position = "bottom"
  )

## Model predictions ----
# Compute predicted mean numbers of grandoffspring and 95% CIs
# from the zero-inflated NB2 model (ZINB2) at completeness = 1
# (z-transformed "mean-centered reference value").


# Compute z-transformed completeness corresponding to raw completeness = 1
mean_comp <- mean(ids_GOmodel_no0$completeness, na.rm = TRUE)
sd_comp   <- sd(ids_GOmodel_no0$completeness, na.rm = TRUE)
max_completeness_z <- (1 - mean_comp) / sd_comp

# Build prediction grid across sexes and offspring numbers
offs_seq <- seq(1, 20, length.out = 200)
pred_grid <- expand.grid(
  sex        = c("female", "male"),
  n_off_type = offs_seq
)
pred_grid$log_offspring  <- log(pred_grid$n_off_type)
pred_grid$completeness_z <- max_completeness_z

# Predict on link scale with standard errors
pred_mat <- predict(
  model_zinb_full,
  newdata = pred_grid,
  type   = "link",
  se.fit = TRUE
)

# Organize predictions and transform to response scale
pred_df <- cbind(
  pred_grid,
  fit = as.numeric(pred_mat$fit),
  se  = as.numeric(pred_mat$se.fit)
)
pred_df$lwr <- pred_df$fit + qnorm(0.025) * pred_df$se
pred_df$upr <- pred_df$fit + qnorm(0.975) * pred_df$se

pred_df$fit_response <- model_zinb_full$modelInfo$family$linkinv(pred_df$fit)
pred_df$lwr_response <- model_zinb_full$modelInfo$family$linkinv(pred_df$lwr)
pred_df$upr_response <- model_zinb_full$modelInfo$family$linkinv(pred_df$upr)

# Split prediction data by sex and define colors
pred_sex_female <- subset(pred_df, sex == "female")
pred_sex_male   <- subset(pred_df, sex == "male")


##  Figure 2 with prediction----
# Add predicted mean lines and 95% CI ribbons to each panel
p_f_response + 
  geom_ribbon(
    data = pred_sex_female,
    aes(x = n_off_type, ymin = lwr_response, ymax = upr_response),
    inherit.aes = FALSE,
    alpha = 0.20, fill = sex_cols["female"]
  ) +
  geom_line(
    data = pred_sex_female,
    aes(x = n_off_type, y = fit_response),
    inherit.aes = FALSE,
    linewidth = 1, color = sex_cols["female"]
  )  -> female_response

p_m_response + 
  geom_ribbon(
    data = pred_sex_male,
    aes(x = n_off_type, ymin = lwr_response, ymax = upr_response),
    inherit.aes = FALSE,
    alpha = 0.20, fill = sex_cols["male"]
  ) +
  geom_line(
    data = pred_sex_male,
    aes(x = n_off_type, y = fit_response),
    inherit.aes = FALSE,
    linewidth = 1, color = sex_cols["male"]
  )  -> male_response


final_response <- (female_response | male_response) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(size = 9, face = "bold"))

final_response <- wrap_elements(final_response) +
  labs(tag = "Number of offspring") +
  theme(
    plot.tag          = element_text(size = 9, face = "plain", hjust = 0.5),
    plot.tag.position = "bottom"
  )

ggsave(
  filename = "figures/figure2.pdf",
  plot     = final_response,
  width    = 184,
  height   = 90,
  units    = "mm"
)

ggsave(
  filename = "figures/figure2.png",
  plot     = final_response,
  width    = 184,
  height   = 90,
  units    = "mm",
  dpi = 600
)

# RESULT: Numbers of descendants of high-ranking mothers increase more along maternal than paternal lineages across generation----

## Calculate alpha-rank proportions and fold-changes per generation-----
# Combine maternal and paternal ancestor tables
table_ancestor_native_r$Lineage_Type <- "Maternal"
table_ancestor_native_male_r$Lineage_Type <- "Paternal"
# Ensure identical columns for merging
missing_cols <- setdiff(names(table_ancestor_native_r), names(table_ancestor_native_male_r))

for (col in missing_cols) {
  table_ancestor_native_male_r[[col]] <- NA
}
# Order males columns as for female
table_ancestor_native_male_r <- table_ancestor_native_male_r[, names(table_ancestor_native_r)]

# Merge tables
combined_table <- rbind(table_ancestor_native_r, table_ancestor_native_male_r)

# Pivot to long format
long_lineage <- combined_table |> 
  pivot_longer(
    cols = starts_with("rank_G"),
    names_to = c(".value", "Generation"),
    names_pattern = "(rank)_G(\\d+)"
  ) |> 
  mutate(
    Generation = -as.integer(Generation),
    Lineage_Type = as.factor(Lineage_Type)
  ) |> 
  rename(
    MaternalRank_at_Birth = rank
  ) |> 
  filter(!is.na(MaternalRank_at_Birth))

# G0 summary for maternal lineage (used to correct paternal G0)
p <- table_ancestor_native_r |> 
  summarise(
    n_total = sum(!is.na(rank_G0)),
    n_alpha = sum(rank_G0 == 1, na.rm = TRUE),
    proportion_alpha = n_alpha / n_total
  )  

# Proportion of alpha-ranked individuals by generation
alpha_proportions_by_generation <- long_lineage |> 
  group_by(Generation, Lineage_Type) |> 
  summarise(
    n_total = sum(!is.na(MaternalRank_at_Birth)),
    n_alpha = sum(MaternalRank_at_Birth == 1, na.rm = TRUE),
    proportion_alpha = n_alpha / n_total,
    .groups = "drop"
  ) |> 
  arrange(Lineage_Type, Generation)

# Correct paternal G0 using female summary since no need to restrict 
# that specific generation to individuals with known fathers
alpha_proportions_by_generation <- alpha_proportions_by_generation |>
  mutate(
    n_total = ifelse(Generation == 0 & Lineage_Type == "Paternal",
                     p$n_total, n_total),
    n_alpha = ifelse(Generation == 0 & Lineage_Type == "Paternal",
                     p$n_alpha, n_alpha),
    proportion_alpha = ifelse(Generation == 0 & Lineage_Type == "Paternal",
                              p$proportion_alpha, proportion_alpha)
  ) |> 
  filter(n_total > 20)

# Calculate fold-changes in alpha proportions across generations
maternal_props <- alpha_proportions_by_generation |>
  filter(Lineage_Type == "Maternal") |>
  arrange(desc(Generation)) |>  # G0 to G-4
  pull(proportion_alpha)

paternal_props <- alpha_proportions_by_generation |>
  filter(Lineage_Type == "Paternal") |>
  arrange(desc(Generation)) |>  # G0 to G-3
  pull(proportion_alpha)

# Calculate fold-changes between generations (n < 20 individuals)
maternal_fold_changes <- maternal_props[-1] / maternal_props[-length(maternal_props)]
maternal_avg_fold_change_sub <- exp(mean(log(maternal_fold_changes), na.rm = TRUE))

paternal_fold_changes <- paternal_props[-1] / paternal_props[-length(paternal_props)]
paternal_avg_fold_change_sub <- exp(mean(log(paternal_fold_changes), na.rm = TRUE))


## Figure 3. Effect of social rank inheritance: distribution of social ranks among ancestors in maternal and paternal lineages-----

# Maternal lineage: 
# Reshape ancestor rank data to long format and compute
# proportional rank distributions per generation


# Identify rank columns (rank_G0, rank_G1, ...)
rank_columns <- names(table_ancestor_native_r)[grepl("^rank_", names(table_ancestor_native_r))]

# Define plotting order (oldest generation first)
names_sequence <-  rev(gsub("rank_", "", rank_columns)) 

# Convert wide ancestor table to long format:
# one row per ancestor per generation
long_data <- table_ancestor_native_r |> 
  dplyr::select(all_of(rank_columns)) |> 
  pivot_longer(cols = everything(), names_to = "generation", values_to = "rank") |>
  mutate(
    generation = gsub("rank_", "", generation),
    generation = factor(generation, levels = names_sequence),
    rank = factor(rank, levels = sort(unique(rank), decreasing = TRUE))
  )
# Remove rows without rank information
long_data_filtered <- long_data |> filter(!is.na(rank))

# Sample size per generation
sample_sizes_female <- long_data_filtered |> 
  group_by(generation) |> 
  summarise(total_n = n(), .groups = "drop")

# Compute proportional distribution of ranks within each generation
long_data_proportions <- long_data_filtered |> 
  group_by(generation, rank) |> 
  summarise(n = n(), .groups = "drop") |> 
  group_by(generation) |> 
  mutate(proportion = n / sum(n)) |> 
  ungroup()

# Paternal lineage
# Apply identical transformation to paternal ancestors.

rank_columns_m <- names(table_ancestor_native_male_r)[grepl("^rank_", names(table_ancestor_native_male_r))]
names_sequence_m <- rev(gsub("rank_", "", rank_columns_m))

long_data_male <- table_ancestor_native_male_r |> 
  dplyr::select(all_of(rank_columns_m)) |> 
  pivot_longer(cols = everything(), names_to = "generation", values_to = "rank") |>
  mutate(
    generation = gsub("rank_", "", generation),
    generation = factor(generation, levels = names_sequence_m),
    rank = factor(rank, levels = sort(unique(rank), decreasing = TRUE))
  )

long_data_filtered_male <- long_data_male |> filter(!is.na(rank))

long_data_proportions_male <- long_data_filtered_male |> 
  group_by(generation, rank) |> 
  summarise(n = n(), .groups = "drop") |> 
  group_by(generation) |> 
  mutate(proportion = n / sum(n)) |> 
  ungroup()


# Standardize G0 across lineages
# Generation 0 corresponds to focal individuals and must
# be identical for maternal and paternal panels.

# Remove male G0 (only IDs with known father, but we should include also the one with unknown father)
long_data_proportions_male <- long_data_proportions_male |> filter(generation != "G0")

# Extract female G0 to add to male dataset
long_data_proportions_G0 <- long_data_proportions |> filter(generation == "G0")

long_data_proportions_male <- bind_rows(long_data_proportions_male, long_data_proportions_G0)

# Update male sample sizes
sample_sizes_male <- long_data_proportions_male |> 
  group_by(generation) |> 
  summarise(total_n = sum(n), .groups = "drop")

# Collapse rare high ranks (≥ 39) and prepare combined dataset
# for visualization

# Combine female and male datasets
full_generations <- sort(unique(long_data_proportions$generation))

# Group ranks ≥ 39 
long_data_proportions1 <- long_data_proportions |> 
  mutate(rank = ifelse(as.numeric(as.character(rank)) >= 39, 39, as.numeric(as.character(rank))),
         generation = factor(generation, levels = full_generations))

long_data_proportions_male1 <- long_data_proportions_male |> 
  mutate(rank = ifelse(as.numeric(as.character(rank)) >= 39, 39, as.numeric(as.character(rank))),
         generation = factor(generation, levels = full_generations))

# Make sure rank is factor (discrete)
long_data_proportions1$rank <- factor(long_data_proportions1$rank, levels = 39:1)
long_data_proportions_male1$rank <- factor(long_data_proportions_male1$rank, levels = 39:1)

# Add 'sex' columns
long_data_proportions1$sex <- "Female"
long_data_proportions_male1$sex <- "Male"

# Combine datasets
long_data_combined <- bind_rows(long_data_proportions1, long_data_proportions_male1)
long_data_combined$sex <- factor(long_data_combined$sex, levels = c("Female", "Male"))

# find the top rank value 
top_rank_val <- min(as.numeric(as.character(long_data_combined$rank)), na.rm = TRUE)
long_data_combined <- long_data_combined |>
  mutate(
    opacity = ifelse(as.numeric(as.character(rank)) == top_rank_val, 1, 0.6)
  )

sample_sizes_combined <- bind_rows(
  sample_sizes_male   |> mutate(sex = "Male"),
  sample_sizes_female |> mutate(sex = "Female")
)
sample_sizes_combined$sex <- factor(sample_sizes_combined$sex,
                                    levels = levels(long_data_combined$sex))

# Select representative simulated lineage:
# Choose the permutation whose geometric mean is closest
# to the overall mean of the null distribution.

# Extract geometric means from permutation results (to be uploaded)
geometric_means_random <- vapply(data_perm, `[[`, numeric(1), "geometric_mean")
geometric_means_random <- geometric_means_random[is.finite(geometric_means_random)]
mean_random <- mean(geometric_means_random, na.rm = TRUE)  # Mean of permutation distribution

# Find index of the permutation closest to this mean
closest_idx <- which.min(abs(geometric_means_random - mean_random))
# Extract that permutation
closest_simulation <- data_perm[[closest_idx]]

closest_simulation$geometric_mean# see its geometric mean
sim <- closest_simulation$long_data_proportions
sim$generation <- gsub("-", "", sim$generation)


# Plot stacked proportional rank distributions
# Each bar represents the distribution of maternal social ranks
# among ancestors in a given generation.

# Ensure the simulated data matches the generation order in the empirical plot
gen_levels <- if (is.factor(long_data_combined$generation)) {
  levels(long_data_combined$generation)
} else {
  sort(unique(long_data_combined$generation))
}

# Ensure rank levels match empirical data; fallback to 39:1 for consistency
rank_levels <- if (is.factor(long_data_combined$rank)) {
  levels(long_data_combined$rank)
} else {
  as.character(39:1)
}

# Format simulated lineage data to match empirical dataset
sim_for_plot <- sim |> 
  mutate(
    rank       = as.numeric(as.character(rank)),
    rank       = ifelse(rank >= 39, 39L, rank),
    rank       = factor(as.character(rank), levels = rank_levels),
    generation = factor(generation, levels = gen_levels),
    sex        = "Simulated"
  )

# Highlight top-ranked ancestors (rank 1) with full opacity
sim_for_plot <- sim_for_plot |>
  mutate(opacity = ifelse(as.numeric(as.character(rank)) == top_rank_val, 1, 0.6))

# Combine empirical and simulated datasets for plotting
long_data_combined3 <- bind_rows(long_data_combined, sim_for_plot)

# Ensure facets are ordered: Female (maternal), Simulated, Male (paternal)
long_data_combined3$sex <- factor(long_data_combined3$sex,
                                  levels = c("Female", "Simulated", "Male"))

# Compute sample sizes for simulated lineage per generation
sample_sizes_simu <- sim_for_plot |>
  group_by(generation) |>
  summarise(total_n = sum(n), .groups = "drop") |>
  mutate(sex = factor("Simulated", levels = levels(long_data_combined3$sex)))

sample_sizes_combined3 <- bind_rows(
  sample_sizes_combined |>
    mutate(sex = factor(as.character(sex), levels = levels(long_data_combined3$sex))),
  sample_sizes_simu
)

# Define facet labels for clarity in the figure
nice_labels <- c(
  "Female"    = "A | Maternal lineage",
  "Simulated" = "B | Simulated lineage with no social inheritance of rank",
  "Male"      = "C | Paternal lineage"
)

# Generate stacked bar plot
long_data_combined3 <- long_data_combined3 |> 
  mutate(generation = factor(generation, levels = rev(levels(long_data_combined3$generation))))

# Same for sample sizes
sample_sizes_combined3 <- sample_sizes_combined3 |> 
  mutate(generation = factor(generation, levels = rev(levels(long_data_combined3$generation))))

plot_combined <- ggplot(long_data_combined3,
                        aes(x = generation, y = proportion, fill = rank, alpha = opacity)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  geom_text(data = sample_sizes_combined3,
            aes(x = generation, y = 1.05, label = total_n),
            inherit.aes = FALSE, size = 9/.pt) +
  scale_fill_viridis_d(
    option = "turbo",
    labels = function(x) { x[x == "39"] <- "39+"; x }
  ) +
  scale_alpha_identity(guide = "none") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, by = 0.1), limits = c(0, 1.15)) +
  labs(x = "Generation", y = "Proportion", fill = "Rank") +
  facet_wrap(~ sex, ncol = 1, labeller = as_labeller(nice_labels)) +
  theme_sadv() +                                    
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold", hjust = 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.spacing.x = unit(0.1, "cm")
  ) +
  guides(fill = guide_legend(ncol = 20, reverse = TRUE))

ggsave(
  filename = "figures/figure3.pdf",
  plot     = plot_combined,
  width    = 184,
  height   = 140,
  units    = "mm"
)

ggsave(
  filename = "figures/figure3.png",
  plot     = plot_combined,
  width    = 184,
  height   = 140,
  units    = "mm",
  dpi = 600
)

# RESULT: Social inheritance increases the prevalence of top-ranking females among individuals' ancestors -----
## Figure S5. Permutation distribution of geometric means for maternal-line rank transmission -----

# Extract geometric means from permutation results (to be uploaded)
geometric_means_random <- vapply(data_perm, `[[`, numeric(1), "geometric_mean")
geometric_means_random <- geometric_means_random[is.finite(geometric_means_random)]

# Observed geometric mean from real maternal-line data
geometric_mean_truth <- maternal_avg_fold_change_sub

# Prepare data for plotting and compute summary statistics
df <- tibble(geometric_mean = geometric_means_random)   # Data frame for ggplot
mean_random <- mean(geometric_means_random, na.rm = TRUE)  # Mean of permutation distribution
# Compute one-sided permutation p-value
p_value <- (sum(geometric_means_random >= geometric_mean_truth, na.rm = TRUE) + 1) /
  (length(geometric_means_random) + 1)

# Plot permutation distribution
plot_perm <- ggplot(df, aes(x = geometric_mean)) +
  geom_histogram(fill = "lightgray", color = "white", bins = 30) +
  geom_vline(xintercept = geometric_mean_truth, color = "red",   linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = mean_random,          color = "black", linetype = "solid",  linewidth = 0.5) +
  annotate("text",
           x = geometric_mean_truth, y = Inf,
           label = paste0("p = ", round(p_value, 4)),
           hjust = 1.1, vjust = 2, size = 9/.pt, color = "red") +
  annotate("text",
           x = mean_random, y = Inf,
           label = "Mean of permutations",
           hjust = 1.1, vjust = 2, size = 9/.pt, color = "black") +
  labs(x = "Geometric Mean", y = "Frequency") +
  theme_sadv()  

# Save figure

ggsave(
  filename = "figures/figureS5.pdf",
  plot     = plot_perm,
  width    = 90,
  height   = 90,
  units    = "mm"
)

ggsave(
  filename = "figures/figureS5.png",
  plot     = plot_perm,
  width    = 90,
  height   = 90,
  units    = "mm",
  dpi = 600
)


#RESULT: Grandoffspring number better predicts long-term genetic contribution than offspring number ----

##  Data preparation ----
ids_completed |> filter(n_off_obs > 0 ) -> data_genet_complete

##  Correlation ----
# Overall correlations
rho_overall_off  <- cor(data_genet_complete$n_off_obs, data_genet_complete$gen_cont, method = "spearman") 
rho_overall_gdoff <- cor(data_genet_complete$n_g_off_obs, data_genet_complete$gen_cont, method = "spearman") 

# Female-specific correlations
rho_f_off  <- cor(data_genet_complete$n_off_obs[data_genet_complete$sex == "female"], data_genet_complete$gen_cont[data_genet_complete$sex == "female"],
    method = "spearman") 
rho_f_gdoff  <- cor(data_genet_complete$n_g_off_obs[data_genet_complete$sex == "female"], data_genet_complete$gen_cont[data_genet_complete$sex == "female"],
    method = "spearman") 

# Male-specific correlations
rho_m_off <- cor(data_genet_complete$n_off_obs[data_genet_complete$sex == "male"], data_genet_complete$gen_cont[data_genet_complete$sex == "male"],
    method = "spearman") 
rho_m_gdoff <- cor(data_genet_complete$n_g_off_obs[data_genet_complete$sex == "male"], data_genet_complete$gen_cont[data_genet_complete$sex == "male"],
    method = "spearman")

# Comparing correlations
n_females <- sum(data_genet_complete$sex == "female") # 239
n_males <- sum(data_genet_complete$sex == "male") # 167

cor_f_off_gdoff <- diffcor.two(rho_f_off, rho_f_gdoff,
            n1 = n_females,
            n2 = n_females, digit = 3)
cor_f_off_gdoff[c("z", "p")]
#       z p
# 1 -4.72 0

cor_m_off_gdoff <- diffcor.two(rho_m_off, rho_m_gdoff,
            n1 = n_males,
            n2 = n_males, digit = 5)
cor_m_off_gdoff[c("z", "p")]
#          z       p
# 1 -2.96235 0.00153

# Combine results in a table
cor_table <- data.frame(
  Sex = c("Female", "Male"),
  Offspring = c(rho_f_off, rho_m_off),
  Grandoffspring = c(rho_f_gdoff, rho_m_gdoff)
) |> left_join(rbind(cbind(cor_f_off_gdoff[c("LL1", "LL2", "UL1", "UL2")], Sex = "Female"),
                     cbind(cor_m_off_gdoff[c("LL1", "LL2", "UL1", "UL2")], Sex = "Male"))) |> 
  rename(lwr_off = LL1, lwr_gdoff = LL2, upr_off = UL1, upr_gdoff = UL2)
  

cor_table
#      Sex Offspring Grandoffspring lwr_off lwr_gdoff upr_off upr_gdoff
# 1 Female 0.6639132      0.8438082 0.58600   0.80300 0.72900   0.87700
# 2   Male 0.6655421      0.8109627 0.57143   0.75167 0.74237   0.85725

##  Figure 4 ----
cor_plot_data <- cor_table |>
  pivot_longer(
    cols      = c(Offspring, Grandoffspring),
    names_to  = "Generation",
    values_to = "Correlation"
  ) |>
  mutate(Generation = factor(Generation, 
                             levels = c("Offspring", "Grandoffspring")),
         lwr = if_else(Generation == "Offspring", lwr_off, lwr_gdoff),
         upr = if_else(Generation == "Offspring", upr_off, upr_gdoff)) |> 
  dplyr::select(-lwr_off, -lwr_gdoff, -upr_off, -upr_gdoff)

figure4 <- ggplot(cor_plot_data, aes(x = Generation, y = Correlation, 
                          color = Sex, group = Sex, shape = Sex, label = round(Correlation, 3))) +
  geom_point(size = 9/.pt, position = position_dodge(width = .2)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), position = position_dodge(width = .2), width = 0.1) +
  geom_line(linewidth = 0.8, position = position_dodge(width = .2)) +
  scale_color_manual(values = c("Female" = "#8624F5", "Male" = "#1FC3AA")) +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(
    x     = "Measure of reproductive success",
    y     = "Correlation with expected genetic contribution",
    color = "Sex"
  ) +
  theme_sadv()+
  theme(
    legend.position      = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background    = element_rect(
      fill      = "white",
      colour    = "black",
      linewidth = 0.5
    )
  )+
  scale_x_discrete(labels = c("Offspring (n)", "Grandoffspring (n)")) +
  geom_text(
    hjust = 1.4,
    vjust = -1,
    size  = 3,
    show.legend = FALSE,
    data = cor_plot_data |> filter(Sex == "Female")
  ) +
  geom_text(
    hjust = -0.7,
    vjust = 1,
    size  = 3,
    show.legend = FALSE,
    data = cor_plot_data |> filter(Sex == "Male")
  ) 

figure4

# Save figure

ggsave(
  filename = "figures/figure4.pdf",
  plot     = figure4,
  width    = 90,
  height   = 90,
  units    = "mm"
)

ggsave(
  filename = "figures/figure4.png",
  plot     = figure4,
  width    = 90,
  height   = 90,
  units    = "mm",
  dpi = 600
)

## ---- Table S2: Descriptive statistics of long-term genetic contribution g_i  ----

# Helper function: compute key descriptive stats for a numeric vector
.desc_stats <- function(x) {
  tibble(
    N          = length(x),
    Mean       = mean(x),
    Median     = median(x),
    SD         = sd(x),
    Min        = min(x),
    Max        = max(x),
    Q1 = quantile(x, 0.25, names = FALSE),
    Q3 = quantile(x, 0.75, names = FALSE)
  )
}

# Prepare data: keep g_i values; standardize sex labels; remove missing g_i
df_gi <- data_genet_complete |> 
  mutate(
    sex = case_when(
      tolower(sex) %in% c("f","female") ~ "Female",
      tolower(sex) %in% c("m","male")   ~ "Male",
      TRUE                              ~ NA_character_
    ),
    gi  = as.numeric(gen_cont),
    .keep = "none"
  ) 

# Compute descriptive stats by sex
tab_by_sex <- df_gi |> 
  filter(!is.na(sex)) |> 
  group_by(Sex = sex) |> 
  summarise(.desc_stats(gi), .groups = "drop")

# Compute pooled "All" row 
tab_all <- df_gi |> 
  summarise(Sex = "All", .desc_stats(gi))

# Combine tables, set factor levels, arrange, and round numeric values
table_S2 <- bind_rows(tab_by_sex, tab_all) |> 
  mutate(Sex = factor(Sex, levels = c("Female", "Male", "All"))) |> 
  arrange(Sex) |> 
  mutate(across(where(is.numeric), ~round(.x, 2)))

# Format for publication: flextable
num_cols <- names(table_S2)[sapply(table_S2, is.numeric)]

ft_S2 <- flextable(table_S2) |>
  theme_booktabs() |>
  bold(part = "header") |>
  align(j = num_cols, align = "center", part = "all") |>
  font(fontname = "Times New Roman", part = "all") |>
  fontsize(size = 12, part = "all") |>
  colformat_num(j = num_cols, digits = 2) |>
  set_table_properties(layout = "fixed") |>
  fit_to_width(max_width = 6.27)  # fit to A4 content width

# Export to Word
doc <- read_docx() |>
  body_add_par(
    "Table S2. Descriptive statistics of expected genetic contribution for all individuals and separated by sex. ",
    style = "Normal"
  ) |>
  body_add_flextable(ft_S2) 

print(doc, target = "tables/Table_S2.docx")

# SUPPLEMENTARY ANALYSIS ----
##Completeness of grandoffspring counts across sexes----

# --- Summary 1: Overall zero-alive-offspring stats --- #
zero_alive_summary <- ids_completed |>
  summarise(
    total = n(),
    n_zero_alive = sum(n_alive_genotyped_offspring == 0),
    prop_zero_alive = mean(n_alive_genotyped_offspring == 0)
  )
zero_alive_summary
#   total n_zero_alive prop_zero_alive
# 1   492          381       0.7743902

# --- Summary 2: Proportion with 0 alive offspring by sex --- #
zero_by_sex <- ids_completed |>
  group_by(sex) |>
  summarise(
    n = n(),
    n_zero_alive = sum(n_alive_genotyped_offspring == 0),
    prop_zero_alive = mean(n_alive_genotyped_offspring == 0)
  )
zero_by_sex
# # A tibble: 2 × 4
#   sex        n n_zero_alive prop_zero_alive
#   <chr>  <int>        <int>           <dbl>
# 1 female   287          229           0.798
# 2 male     205          152           0.741

# --- Chi-squared test: 2x2 table of sex vs. zero-alive-offspring --- #
table_zero <- ids_completed |>
  mutate(has_zero_alive = n_alive_genotyped_offspring == 0) |>
  count(sex, has_zero_alive) |>
  pivot_wider(names_from = has_zero_alive, values_from = n, values_fill = 0) |>
  column_to_rownames("sex") |>
  as.matrix()

fisher_result <- fisher.test(table_zero)
fisher_result
# Fisher's Exact Test for Count Data
# 
# data:  table_zero
# p-value = 0.1553
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.464573 1.138192
# sample estimates:
# odds ratio 
#  0.7268686 


# --- Summary 3: Among individuals with >0 alive offspring, get stats by sex --- #
mean_alive_among_nonzero <- ids_completed |>
  dplyr::filter(n_alive_genotyped_offspring > 0) |>
  dplyr::group_by(sex) |>
  dplyr::summarise(
    mean_n_alive = mean(n_alive_genotyped_offspring),
    median_n_alive = median(n_alive_genotyped_offspring),
    sd_n_alive = sd(n_alive_genotyped_offspring),
    max_n_alive = max(n_alive_genotyped_offspring),
    .groups = "drop"
  )

# --- Wilcoxon test for difference in alive offspring counts by sex --- #
df_nonzero <- ids_completed |>
  dplyr::filter(n_alive_genotyped_offspring > 0)

wilcox_result <- wilcox.test(n_alive_genotyped_offspring ~ sex, data = df_nonzero)
wilcox_result
# W = 1363.5, p-value = 0.2579


# --- Summary 4:completeness stats  --- #

completeness_summary <- ids_completed |>
  summarise(
    mean_completeness = mean(completeness, na.rm = TRUE),
    median_completeness = median(completeness, na.rm = TRUE),
    sd_completeness = sd(completeness, na.rm = TRUE),
    min_completeness = min(completeness, na.rm = TRUE),
    max_completeness = max(completeness, na.rm = TRUE),
    n_missing = sum(is.na(completeness)),
    n_total = n()
  )
completeness_summary
#   mean_completeness median_completeness sd_completeness min_completeness max_completeness n_missing n_total
# 1         0.9316829                   1       0.1617974                0                1         0     492

completeness_by_sex <- ids_completed |>
  group_by(sex) |>
  summarise(
    mean_completeness   = mean(completeness, na.rm = TRUE),
    median_completeness = median(completeness, na.rm = TRUE),
    sd_completeness     = sd(completeness, na.rm = TRUE),
    min_completeness    = min(completeness, na.rm = TRUE),
    max_completeness    = max(completeness, na.rm = TRUE),
    n_missing           = sum(is.na(completeness)),
    n_total             = n(),
    .groups = "drop"
  )
completeness_by_sex
# # A tibble: 2 × 8
#   sex    mean_completeness median_completeness sd_completeness min_completeness max_completeness n_missing n_total
#   <chr>              <dbl>               <dbl>           <dbl>            <dbl>            <dbl>     <int>   <int>
# 1 female             0.944                   1           0.135                0                1         0     287
# 2 male               0.914                   1           0.192                0                1         0     205


wilcox.test(completeness ~ sex, data = ids_completed)
# W = 31220, p-value = 0.1133


## Mean reproductive tenure in females and males ----

# Mean tenure 
tenure_summary <- ids_completed |> 
  group_by(sex) |> 
  summarise(
    n = n(),
    mean_tenure = mean(tenure, na.rm = TRUE),
    sd_tenure = sd(tenure, na.rm = TRUE),
    min_tenure = min(tenure, na.rm = TRUE),
    max_tenure = max(tenure, na.rm = TRUE)
  )
tenure_summary
# # A tibble: 2 × 6
#   sex        n mean_tenure sd_tenure min_tenure max_tenure
#   <chr>  <int>       <dbl>     <dbl>      <dbl>      <dbl>
# 1 female   287        7.32      4.35    0.00274       16.5
# 2 male     205        6.46      3.67    0.110         15.3


# Test the difference between sex with parametric test
wilcox.test(tenure ~ sex, data = ids_completed)
#W = 32684, p-value = 0.03567

# Mean tenure among the individuals with tenure > 1 year:
tenure_1y_summary <- repro_data |> 
  group_by(sex) |> 
  summarise(
    n = n(),
    mean_tenure = mean(tenure, na.rm = TRUE),
    sd_tenure = sd(tenure, na.rm = TRUE),
    min_tenure = min(tenure, na.rm = TRUE),
    max_tenure = max(tenure, na.rm = TRUE)
  )
tenure_1y_summary
# # A tibble: 2 × 6
#   sex        n mean_tenure sd_tenure min_tenure max_tenure
#   <chr>  <int>       <dbl>     <dbl>      <dbl>      <dbl>
# 1 female   268        7.80      4.08       1.07       16.5
# 2 male     188        6.99      3.36       1.03       15.3

# Test the difference between sex with parametric test
wilcox.test(tenure ~ sex, data = repro_data)
#W = 27992, p-value = 0.04324

## Fig S1-S3: data preparation ----

ids_barplot <- ids_completed |>
  mutate(
    sex        = factor(sex, levels = c("female", "male")),
    off_fac    = factor(n_off_type, levels = 0:max(n_off_type, na.rm = TRUE)),# create factor versions for plotting (so all levels are represented in bars)
    go_fac     = factor(n_g_off,    levels = 0:max(n_g_off,    na.rm = TRUE))
  )

ids_barplot_GGO <- ids_completed_GGO |>
  mutate(
    sex        = factor(sex, levels = c("female", "male")),
    ggo_fac    = factor(n_gg_off,   levels = 0:max(n_gg_off,   na.rm = TRUE))
  )

cols <- c("female" = "#8624F5", "male" = "#1FC3AA")

## Figure S1: Number of Offspring----
# Female panel
p_female <- ggplot(subset(ids_barplot, sex == "female"), aes(x = off_fac)) +
  geom_bar(fill = cols["female"], colour = "black", linewidth = 0.3, width = 0.7) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65), breaks = seq(0, 65, by = 5)) +
  labs(x = NULL, y = "Count") +
  theme_sadv() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )
# Male panel
p_male <- ggplot(subset(ids_barplot, sex == "male"), aes(x = off_fac)) +
  geom_bar(fill = cols["male"], colour = "black", linewidth = 0.3, width = 0.7) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 65), breaks = seq(0, 65, by = 5)) +
  labs(x = "Number of Offspring", y = "Count") +
  theme_sadv() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

# Combine panels vertically
combined_plot <- p_female / p_male +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

ggsave(
  filename = "figures/figureS1.pdf",
  plot     = combined_plot,
  width    = 184,
  height   = 184,
  units    = "mm"
)

ggsave(
  filename = "figures/figureS1.png",
  plot     = combined_plot,
  width    = 184,
  height   = 184,
  units    = "mm",
  dpi = 600
)

## Figure S2: Number of Grandoffspring----

# Female panel
p_female_go <- ggplot(
  subset(ids_barplot, sex == "female"),
  aes(x = go_fac)
) +
  geom_bar(fill = cols["female"], colour = "black", linewidth = 0.3, width = 0.7) +
  scale_x_discrete(
    drop = FALSE,
    breaks = as.character(seq(0, 106, by = 5))  # Tick labels every 5
  ) +
  scale_y_continuous(
    breaks = c(0, 5, 10, 15, 134, 135),
    expand = c(0, 0)
  ) +
  scale_y_break(c(15, 134), scales = 0.3) +
  labs(x = NULL, y = "Count") +
  theme_sadv() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Male panel
p_male_go <- ggplot(
  subset(ids_barplot, sex == "male"),
  aes(x = go_fac)
) +
  geom_bar(fill = cols["male"], colour = "black", linewidth = 0.3, width = 0.7) +
  scale_x_discrete(
    drop = FALSE,
    breaks = as.character(seq(0, 106, by = 5))
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = c(0, 5, 10, 15, 68, 69)
  ) +
  scale_y_break(c(16, 68), scales = 0.3) +
  labs(x = "Number of Grandoffspring", y = "Count") +
  theme_sadv() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )


# Combine panels vertically
combined_plot_go <- p_female_go / p_male_go +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

ggsave(
  filename = "figures/figureS2.pdf",
  plot     = combined_plot_go,
  width    = 184,
  height   = 184,
  units    = "mm"
)

ggsave(
  filename = "figures/figureS2.png",
  plot     = combined_plot_go,
  width    = 184,
  height   = 184,
  units    = "mm",
  dpi = 600
)

## Figure S3: Number of Great-grandoffspring ----
# Female panel

p_female_great <- ggplot(
  subset(ids_barplot_GGO, sex == "female"),
  aes(x = ggo_fac)
) +
  geom_bar(fill = cols["female"], colour = "black", linewidth = 0.3, width = 0.7,   alpha = 0.5) +
  scale_x_discrete(
    drop = FALSE,
    breaks = as.character(seq(0, 156, by = 5))  # Tick labels every 5
  ) +
  scale_y_continuous(
    breaks = c(0, 5, 10, 94, 95), # last one is the max number I have of counts
    expand = c(0, 0)
  ) +
  scale_y_break(c(10, 94), scales = 0.3) +
  labs(x = NULL, y = "Count") +
  theme_sadv() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Male panel
p_male_great <- ggplot(
  subset(ids_barplot_GGO, sex == "male"),
  aes(x = ggo_fac)
) +
  geom_bar(fill = cols["male"], colour = "black", linewidth = 0.3, width = 0.7, alpha = 0.5) +
  scale_x_discrete(
    drop = FALSE,
    breaks = as.character(seq(0, 156, by = 5))  # Tick labels every 5
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = c(0, 5, 10, 52, 53)
  ) +
  scale_y_break(c(10, 52), scales = 0.3) +
  labs(x = "Number of Great-grandoffspring", y = "Count") +
  theme_sadv() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Combine panels vertically
combined_plot_great <- p_female_great / p_male_great +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

ggsave(
  filename = "figures/figureS3.pdf",
  plot     = combined_plot_great,
  width    = 184,
  height   = 184,
  units    = "mm"
)

ggsave(
  filename = "figures/figureS3.png",
  plot     = combined_plot_great,
  width    = 184,
  height   = 184,
  units    = "mm",
  dpi = 600
)
