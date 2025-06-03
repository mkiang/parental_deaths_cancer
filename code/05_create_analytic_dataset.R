## 03_create_analytic_dataset.R ----
## 
## Create analytic dataset.

## Imports ----
library(tidyverse)
library(here)
library(fs)
library(foreach)
library(doParallel)
library(narcan)  # remotes::install_github("mkiang/narcan")
library(future)
library(furrr)
source(here::here("code", "utils.R"))

## Bug workaround
## See: https://github.com/rstudio/rstudio/issues/6692
## Revert to 'sequential' setup of PSOCK cluster in RStudio Console on macOS and R 4.0.0
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
    parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}

## Constants ----
DATA_FOLDER <- here::here("data")

## Data ----
all_pops_long <- readRDS(here::here(
    DATA_FOLDER,
    "pop_est_national_single_year_age_1990to2020_long.RDS"
))
summarized_deaths <- readRDS(here::here(DATA_FOLDER, "summarized_deaths_1990_2020.RDS"))

## Create the actual analytic dataset ----
pop_skeleton <- all_pops_long |> 
    dplyr::mutate(race_eth = dplyr::case_when(
        hispanic == 1 ~ "hispanic",
        TRUE ~ race
    ),
    sex = dplyr::case_when(female == 1 ~ "female",
                    female == 0 ~ "male",
                    TRUE ~ NA_character_)) |> 
    dplyr::rename(age_years = age) |> 
    dplyr::group_by(year, race_eth, sex, age_years) |> 
    dplyr::summarize(pop = sum(pop_est)) |> 
    dplyr::ungroup()

analytic_df <- pop_skeleton |> 
    dplyr::left_join(summarized_deaths) |> 
    dplyr::mutate(across(starts_with("n_"), .fns = ~replace_na(., 0))) |>
    dplyr::mutate(across(n_cancer:n_multi, .fns = ~ifelse(year < 1999, NA, .)))

## Save ----
saveRDS(
    analytic_df, 
    here::here(
        DATA_FOLDER,
        "national_year_age-sex-race_cancer-mortality.RDS"
    ),
    compress = "xz"
)

## Create life tables ----
df.dth <- bind_rows(
    analytic_df |>
        rename("age" = age_years),
    analytic_df |>
        rename("age" = age_years) |>
        group_by(year, sex, age) |>
        summarize(across(pop:n_deaths, .fns = ~ sum(.))) |>
        ungroup() |>
        mutate(race_eth = "total")
)

## Extract constants ----
years <- unique(df.dth$year)
n.years <- length(years)
ages <- unique(df.dth$age)
n.ages <- length(ages)
races <- unique(df.dth$race_eth)
n.races <- length(races)
radix <- 100000

## Make a life table
lt.US <- df.dth |>
    # Assume same ax for all races/ethnicities
    group_by(year, race_eth, sex) |> 
    mutate(
        mx = n_deaths/pop,
        ax = case_when(age == 0 ~ 0.07 + 1.7*mx[1],
                       age == max(ages) ~ Inf,
                       TRUE ~ 0.5),
        n = case_when(age == max(ages) ~ Inf,
                      TRUE ~ 1),
        qx = (n * mx)/(1 + (n - ax) * mx),
        qx = ifelse(age == max(ages), 1, qx),
        px = 1 - qx,
        lx = cumprod(c(1, head(px, -1)))*radix,
        dx = lx - lead(lx, 1),
        dx = ifelse(age == max(ages), lx, dx),
        Lx = lead(lx)*n + ax*dx,
        # Last age gp is NA due to lead()
        Lx = ifelse(age == max(ages), dx/mx, Lx),
        ex = rev(cumsum(rev(Lx)))/lx
    ) |>
    ungroup()

## Save life table
saveRDS(
    lt.US,
    here::here(
        DATA_FOLDER,
        "national_lifetable.RDS"),
    compress = "xz"
)
