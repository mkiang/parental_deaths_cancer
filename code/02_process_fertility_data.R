## 02_process_fertility_data.R ----
##
## Takes the extracted fertility data and (1) interpolates the annual
## fertility rates and (2) generates the modeled male fertility rates.
##
## NOTE: Code comes from Ben's github repo based on our JAMA paper.
## See https://github.com/benjisamschlu/parental_deaths for original.

## Imports ----
library(tidyverse)
library(here)
library(pracma)
source(here::here("code", "utils.R"))

## Data ----
fx <- readRDS(here::here("data", "fertility_female_processed.RDS"))
fx.male <- readRDS(here::here("data", "fertility_male_processed.RDS"))

## 1. Interpolate female fertility rates ----
fx.female <- fx |>
    # Mid age class
    dplyr::mutate(age.m = dplyr::case_when(
        age == 45 ~ age + 5,
        # Last age group=45-54 yo
        TRUE ~ age + 2.5
    )) %>%
    dplyr::arrange(race_eth, year, age) %>%
    dplyr::group_by(race_eth, year) %>%
    dplyr::group_modify(~ intra_fx(.x, sex = "female")) |>
    dplyr::ungroup() |>
    dplyr::filter(race_eth %in% c("black", "hispanic", "white", "total")) |>
    dplyr::select(race_eth , year, age, asfr = fx, sex)

## 2. Model male fertility rates ----
### Compute TFR for male/female model ----
tfr <- dplyr::bind_rows(
    fx.female |>
        dplyr::filter(race_eth == "total") |>
        dplyr::summarize(.by = c(year), tfr = sum(asfr)) |>
        dplyr::mutate(sex = "female"),
    fx.male |>
        dplyr::ungroup() |>
        dplyr::filter(race_eth == "total") |>
        dplyr::summarize(.by = c(year), tfr = sum(asfr)) |>
        dplyr::mutate(sex = "male")
)

### Compute TFR Ratio ----
df.y.pred <- tibble::tibble(year = 2016:2020)

## In line with Schoumaker (2019), the ratio is often < 1.
tfr.ratio <- tfr |>
    tidyr::pivot_wider(names_from = sex, values_from = tfr) |>
    dplyr::mutate(ratio = male / female)

### Linear model on TFR ratio 2010-2015
ratio.lm <- stats::lm(ratio ~ 1 + year, data = tfr.ratio |>
                          dplyr::filter(year >= 2010, !is.na(ratio)))

### Get predicted tfr ratio for 2016-2020
ratio.tfr.2016_20 <- tibble::tibble(year = 2016:2020, ratio = predict(ratio.lm, df.y.pred))

### Combine observed data with prediction
tfr.ratio.1990_2020 <- dplyr::bind_rows(
    tfr.ratio |>
        dplyr::filter(year < 2016) |>
        dplyr::select(year, ratio),
    ratio.tfr.2016_20
)

### Extract tfr female 2016-2020 to rescale forecast of normalized male fx
tfr.f.2016_20 <- fx.female |>
    dplyr::filter(race_eth == "total", year > 2015) |>
    dplyr::summarize(.by = c(year), tfr = sum(asfr)) |>
    # Add linearly forecasted tfr ratio 2016-2020
    dplyr::mutate(ratio.tfr = ratio.tfr.2016_20$ratio)

### Forecast male fertility rates
fx.male.frcst <- fx.male |>
    dplyr::filter(year >= 2010) |>
    # Forecast fx using 2010-2015 period, on each age
    dplyr::group_by(age) |>
    tidyr::nest() |>
    dplyr::mutate(
        # Model on log scale to avoid negative value
        model = purrr::map(data, function(df)
            stats::lm(log(asfr) ~ 1 + year, data = df)),
        # Get prediction
        pred = purrr::map(model, function(.lm)
            stats::predict(.lm, df.y.pred))
    ) |>
    dplyr::select(-c(model, data)) |>
    tidyr::unnest(cols = c(pred)) |>
    dplyr::mutate(
        asfr = exp(pred),
        year = df.y.pred$year,
        sex = "male",
        race_eth = "total"
    ) |>
    dplyr::ungroup() |>
    # Normalise forcasted male fx
    dplyr::mutate(.by = (year), asfr.std = asfr / sum(asfr)) |>
    # Add female TFR and 2015 tfr ratio to scale forecast for male
    dplyr::left_join(tfr.f.2016_20, by = c("year")) |>
    dplyr::mutate(
        # Rescale normalized male fx
        asfr = asfr.std * tfr * ratio.tfr
    )

### Combine data with forecast
fx.male.1990_2020 <-
    dplyr::bind_rows(fx.male.frcst, fx.male)

### Estimate linear model on ratio between female fx race and female
### fx population
df.ratio.modeling <- fx.female |>
    dplyr::filter(race_eth != "total") |>
    dplyr::select(-sex) |>
    dplyr::rename("asfr.race" = asfr) |>
    # Add female fx at population level
    dplyr::left_join(
        fx.female |>
            dplyr::filter(race_eth == "total") |>
            dplyr::select(-c(race_eth, sex)),
        by = c("year", "age")
    ) |>
    dplyr::mutate(
        .by = c(race_eth, year),
        
        asfr.race.std = asfr.race / sum(asfr.race),
        asfr.std = asfr / sum(asfr)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(ratio = asfr.race.std / asfr.std) |>
    # Cases where 0/0 or 1/0
    dplyr::filter(!is.na(ratio), !is.infinite(ratio)) |>
    dplyr::group_by(race_eth, age) |>
    tidyr::nest() |>
    # Linear model of the ratio between fx race and fx population
    dplyr::mutate(model = purrr::map(data, function(df)
        stats::lm(ratio ~ 1 + year, data = df)))

### Extract SE for model uncertainty propagation of race-specific male fx
df.ratio.se <- df.ratio.modeling |>
    dplyr::mutate(
        pred = purrr::map2(model, data, function(.lm, .data)
            stats::predict(.lm, .data)),
        se = purrr::map2(model, data, function(.lm, .data)
            stats::predict(.lm, .data, se.fit = TRUE)$se.fit)
    ) |>
    dplyr::select(-model) |>
    tidyr::unnest(cols = c(data, pred, se))

### Extract betas from model
df.ratio.pars <-
    df.ratio.modeling |>
    dplyr::mutate(lm_tidy = purrr::map(model, broom::tidy)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(data, model)) |>
    tidyr::unnest(cols = c(lm_tidy))

### Now apply modeled estimates on male fertility
### CONSTANTS
years <- 1990:2020
n.years <- length(years)
ages <- unique(fx.male$age)

# Container
df.fx.male.fit <- tibble::tibble()
for (r in c("black", "hispanic", "white")) {
    for (x in ages) {
        # No ratio model above 55 yo since last age available for women
        # and 55 yo is characterize by high noise.
        # Use 54 yo estimates for all ages >= 55
        if (x < 55) {
            beta <- df.ratio.pars |>
                dplyr::filter(race_eth == r, age == x) |>
                dplyr::pull(estimate)
            
            # Get se for model uncertainty propagation
            prediction.se <-
                df.ratio.se |>
                dplyr::filter(race_eth == r, age == x) |>
                dplyr::arrange(year) |>
                dplyr::pull(se)
            
        } else {
            beta <- df.ratio.pars |>
                dplyr::filter(race_eth == r, age == 54) |>
                dplyr::pull(estimate)
            
            prediction.se <-
                df.ratio.se |>
                dplyr::filter(race_eth == r, age == 54) |>
                dplyr::arrange(year) |>
                dplyr::pull(se)
        }
        
        # Design matrix
        X <- matrix(cbind(rep(1, n.years), years), nrow = n.years, ncol = 2)
        
        # Apply female ratio estimates to male
        ratio <- X %*% beta
        
        # Store as df
        df.out <-
            tibble::tibble(
                race_eth = rep(r, n.years),
                age = rep(x, n.years),
                year = years,
                fit.r = as.numeric(ratio),
                se = prediction.se
            )
        
        df.fx.male.fit <- rbind(df.fx.male.fit, df.out)
    }
}

# Add known asfr from male at population level
df.fx.male.fit <-
    df.fx.male.fit |>
    dplyr::left_join(
        fx.male.1990_2020 |>
            dplyr::select(year, age, asfr.total = asfr, sex),
        by = c("year", "age")
    ) |>
    # standardize fx male population level
    dplyr::mutate(.by = c(race_eth, year),
                  
                  
                  asfr.std = asfr.total / sum(asfr.total)) |>
    dplyr::ungroup() |>
    # Add race-specific TFR from female
    dplyr::left_join(
        fx.female |>
            dplyr::filter(race_eth != "total") |>
            dplyr::summarise(.by = c(race_eth, year), tfr.female = sum(asfr)),
        by = c("race_eth", "year")
    ) |>
    # Add ratio between male and female tfr at pop level
    dplyr::left_join(tfr.ratio.1990_2020, by = c("year")) |>
    # Compute fx for male of each race/ethnicity
    # 95% CI bounds for uncertainty propagation in MC
    dplyr::mutate(
        asfr.male.fit = fit.r * asfr.std * tfr.female * ratio,
        
        asfr.male.fit.l95 = (fit.r - 1.96 * se) * asfr.std * tfr.female * ratio,
        asfr.male.fit.u95 = (fit.r + 1.96 * se) * asfr.std * tfr.female * ratio
    )

# Bind female and male fx
df.fx <- dplyr::bind_rows(
    # Female section
    fx.female |>
        dplyr::rename("fx" = asfr),
    # Male section
    dplyr::bind_rows(
        df.fx.male.fit |>
            dplyr::select(
                race_eth,
                sex,
                age,
                year,
                fx = asfr.male.fit,
                fx.l95 = asfr.male.fit.l95,
                fx.u95 = asfr.male.fit.u95
            ),
        # Add population level male fx
        fx.male.1990_2020 |>
            dplyr::select(race_eth, sex, age, year, fx = asfr)
    )
)

## Save ----
saveRDS(df.fx, here::here("data", "fx_US.rds"), compress = "xz")
