## 01_process_fertility_data.R ----
##
## Takes the raw fertility data and returns a cleaned up data set. Note that
## the raw data come from:
##
##      1. An excel file (brth1980_2019.xlsx) from CDC for 1990-2019 (female)
##      2. A PDF file (nvsr72-01.pdf) from NVSS for 2020 (female)
##      3. A text delimited file from FertilityData.org for 1990-2020(male)
##
## NOTE: Code comes from Ben's github repo based on our JAMA paper.
## See https://github.com/benjisamschlu/parental_deaths for original.

## Imports ----
library(tidyverse)
library(here)
library(openxlsx)
library(pdftools)

## 1. Female fertilty data from CDC for 1990-2019 ----
## https://www.cdc.gov/nchs/hus/data-finder.htm?year=2020-2021&table=Table%20Brth
## Downloaded on the 24/05/2023 (original repo)
## NOTE: last age interval is 10 year wide (45-54)

### Excel file constants ----
columns.xlsx    <- c(1, 4, 5, 8:13)
age.gp.xlsx     <- seq(10, 45, 5)
columns.names   <- c("year", as.character(age.gp.xlsx))
rows.total.xlsx <- 21:50   # Total
rows.nhw.xlsx   <- 109:142 # White, not Hispanic or Latina
rows.nhb.xlsx   <- 190:223 # Black or African American, not Hispanic or Latina
rows.aian.xlsx  <- 262:285 # American Indian or Alaska Native
rows.api.xlsx   <- 320:339 # Asian or Pacific Islander, not Hispanic or Latina
rows.h.xlsx     <- 352:381 # Hispanic or Latina
ethnic.gp       <- c("total", "white", "black", "aian", "api", "hispanic")
years.xlsx      <- 1990:2019
n.years.xlsx    <- length(years.xlsx)

### Read in ----
fx <- openxlsx::read.xlsx(
    here::here("data_raw", "brth1980_2019.xlsx"),
    rows = c(
        rows.total.xlsx,
        rows.nhw.xlsx,
        rows.nhb.xlsx,
        rows.aian.xlsx,
        rows.api.xlsx,
        rows.h.xlsx
    ),
    cols = columns.xlsx,
    colNames = FALSE
) |>
    tibble::as_tibble()

### More descriptive (temporary) column names
fx <- fx |>
    dplyr::rename_at(paste0("X", 1:9), ~ columns.names)

### Remove years with "single race" estimates
fx <- fx |>
    dplyr::filter(!grepl("(single race)", year))

### Change suppressed birth rates to 0
fx <- fx |>
    dplyr::mutate(dplyr::across(dplyr::all_of(columns.names), 
                                ~ sub("*", 0, .x, fixed = TRUE)))

### Add ethnic group variable according to number of years
fx <- fx |>
    dplyr::mutate(race_eth = rep(ethnic.gp, c(
        rep(n.years.xlsx, 3), rep(n.years.xlsx - 10, 2), n.years.xlsx
    )))

### Reshape to long format
fx <- fx |>
    tidyr::pivot_longer(as.character(age.gp.xlsx),
                        values_to = "fx",
                        names_to = "age") %>%
    dplyr::mutate(dplyr::across(c(year, age, fx), ~ as.numeric(.)))

## 2. Fertility data from CDC for 2020 (pdf) ----
## Get fertility rates for 2020 from latest NVSS report (Jan 31 2023)
## Source: https://www.cdc.gov/nchs/data/nvsr/nvsr72/nvsr72-01.pdf
## Consulted on the 24/05/2023

### PDF file constants ----
age.fx <- unique(fx$age)
n.age.fx <- length(age.fx)
page_num <- 13:14 # Pages of Table 2
cols_table <- c(3:4, 7:12) # Cols of table with fx values
id_race <- c(16, 30, 37, 44, 16) # Vector items corresponding to racial/ethnic groups (2020)

### Import pdf ----
pdf_doc <- pdftools::pdf_text(here::here("data_raw", "nvsr72-01.pdf"))

### Loop and extract ----
### NOTE: API not available
fx.latest <- tibble::tibble()
for (r in ethnic.gp[ethnic.gp != "api"]) {
    if (r == "hispanic") {
        p <- 2
    } else {
        p <- 1
    }
    
    # Extract table
    table_pdf <- pdf_doc[page_num[p]]
    
    # Create vector with elements before \n
    table_pdf <- strsplit(table_pdf, "\n")[[1]]
    
    # Extract fx for year 2020
    fx_2020 <- table_pdf[id_race[which(ethnic.gp[ethnic.gp != "api"] == r)]]
    
    # Cleaning
    fx_2020 <- strsplit(fx_2020, " ")[[1]]
    fx_2020 <- fx_2020[!fx_2020 %in% c("", ".")]
    fx_2020_out <- tibble::tibble(
        year = as.integer(fx_2020[1]),
        age = age.fx,
        fx = as.numeric(fx_2020[cols_table]),
        race_eth = r
    )
    fx.latest <- bind_rows(fx.latest, fx_2020_out)
}

## 3. Male fertility data from human fertility collection ----
## Source: https://www.fertilitydata.org/Country/Country?code=USA
fx.male <- utils::read.table(here::here("data_raw", "m_ASFR_USA.txt"),
                             header = TRUE) |>
    dplyr::rename_with(tolower) |>
    dplyr::filter(year >= 1990) |>
    dplyr::mutate(sex = "male", race_eth = "total") |> 
    tibble::as_tibble()

## Save ----
saveRDS(
    bind_rows(fx, fx.latest) %>%
        # fx expressed per 1000 women
        dplyr::mutate(fx = fx / 1000),
    here::here("data", "fertility_female_processed.RDS")
)

saveRDS(fx.male, here::here("data", "fertility_male_processed.RDS"))
