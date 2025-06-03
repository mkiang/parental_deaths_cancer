## Imports ----
library(tidyverse)
library(here)
library(fs)
library(foreach)
library(haven)
library(doParallel)
library(narcan)  # remotes::install_github("mkiang/narcan")
library(future)
library(furrr)
source(here::here("code", "utils.R"))

## Bug workaround
## See: https://github.com/rstudio/rstudio/issues/6692
## Revert to 'sequential' setup of PSOCK cluster in RStudio Console on macOS and R 4.0.0
if (Sys.getenv("RSTUDIO") == "1" &&
    !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
    parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}

## Constants ----
START_YEAR <- 1999
END_YEAR <-  2020
NCORES <-  4
RAW_FOLDER <- here::here("data_raw")
DATA_FOLDER <- here::here("data_processed")
FORCE_REFRESH <- TRUE
KEEP_ZIPS <- TRUE
years <- START_YEAR:END_YEAR

## Demographics of cancer deaths
if (!file_exists(here::here("data", "death_demographics.RDS"))) {
    future::plan(future::multisession(workers = NCORES))
    summarized_deaths <- furrr::future_map_dfr(
        .x = years,
        .f = ~ {
            ### Read in processed death file ----
            temp_x <- readRDS(here::here(DATA_FOLDER,
                sprintf("processed_mcod_%i.RDS", .x)))

            ### Recode race to match population data ----
            temp_x <- temp_x |>
                dplyr::mutate(race_eth = dplyr::case_when(hspanicr %in% c(1:5, 9) ~ "hispanic",
                    TRUE ~ race_bridged))

            ### Flag death types ----
            ### Malignant neoplasms C00-C97
            temp_x <- temp_x |>
                dplyr::mutate(all_cancer_death = 0 + grepl(paste0(
                    "\\<C[0-8]{1}[0-9]{1}|", "\\<C9[0-7]{1}"
                ), ucod),
                all_deaths = 1)

            ### Summarize ----
            temp_x |>
                dplyr::group_by(year, age_years, sex, race_eth) |>
                dplyr::summarize(
                    n_cancer = sum(all_cancer_death),
                    n_deaths = sum(all_deaths)
                ) |>
                dplyr::ungroup() |>
                dplyr::arrange(year, race_eth, sex, age_years)
        })

    ### Close out ----
    doParallel::stopImplicitCluster()
    closeAllConnections()

    ### Save ----
    saveRDS(
        summarized_deaths,
        here::here("data", "death_demographics.RDS"),
        compress = "xz"
    )
}
