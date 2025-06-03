## 06_estimate_overdispersion.R ----
## 
## Instead of assuming a Poisson distribution for births and deaths, we want
## to assume a negative binomial distribution since we know that births
## and deaths are (usually) overdispersed and a Poisson distribution will 
## give us smaller confidence intervals than is appropriate. This file 
## estimates the overdispersion parameter from the data and saves it for the
## simulations. 
## 
## Code comes from Ben's github repo based on our JAMA paper:
##      https://github.com/benjisamschlu/parental_deaths

## Imports ----
library(tidyverse)
library(here)

## Data ----
lt.US <- readRDS(here::here("data", "national_lifetable.RDS"))
df.fx <- readRDS(here::here("data", "fx_US.rds"))

## Overdispersion in the mortality data ----
theta.dth <- lt.US |>
    filter(year >= 1999,
           race_eth %in% c("total", "black", "white", "hispanic")) |>
    dplyr::select(year, race_eth, sex, age, n_deaths) |>
    summarize(
        .by = c(race_eth, sex, age),
        mu = mean(n_deaths),
        var = var(n_deaths)
    ) |>
    ungroup() |>
    ## There are 4 ages with underdispersion — likely due to small sample
    ## size. Assume a normal poisson variance here. 
    mutate(theta = (mu ^ 2) / (var - mu),
           theta = ifelse(theta < 0, 1e6, theta))

## Overdispersion in birth data ----
theta.bth <- df.fx |> 
    # Add population counts
    left_join(
        lt.US |> 
            dplyr::select(
                year, race_eth, sex, age, pop
            ),
        by = c("year", "race_eth", "sex", "age")
    ) |>  
    mutate(
        birth = fx * pop
    ) |> 
    summarise(
        .by = c(race_eth, sex, age),
        
        mu = mean(birth),
        var = var(birth)
    ) |> 
    ungroup() |> 
    ## Assume no reproduction at very young ages (so theta does not matter) and
    ## assume no underdispersion (instead use Poisson). 
    mutate(
        theta = (mu^2) / (var - mu),
        theta = case_when(
            is.na(theta) ~ 1e6, 
            theta < 0 ~ 1e6, 
            TRUE ~ theta
        )
    )

## Save ----
saveRDS(theta.dth,
        here("data", "overdispersion_deaths.RDS"),
        compress = "xz")

# Overdispersion births
saveRDS(theta.bth,
        here("data", "overdispersion_births.RDS"),
        compress = "xz")
