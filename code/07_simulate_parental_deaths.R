## Imports ----
library(MASS)
library(tidyverse)
library(here)
library(fs)
library(doParallel)
library(future)
library(future.apply)
source(here("code", "utils.R"))

## Data ----
lt.US <- readRDS(here("data", "national_lifetable.RDS"))

### Overdispersion parameters for NBin in MC
theta.bth <- readRDS(here("data", "overdispersion_births.RDS"))

### Fertility data (new input for male fertility rates)
fx.US <- readRDS(here("data", "fx_US.rds"))

## CONSTANTS ----
N_SIM <- 2000
N_CORES <- 16
FORCE_REFRESH <- FALSE # TRUE

## Clean up life table 
lt.US <- lt.US %>%
    filter(race_eth %in% c("black", "white", "hispanic", "total")) |>
    mutate(n_other = n_deaths - n_cancer) |> 
    mutate(
        mx_cancer = n_cancer / pop,
        mx_bone = n_bone / pop,
        mx_breast = n_breast / pop,
        mx_digestive = n_digestive / pop,
        mx_eye = n_eye / pop,
        mx_female = n_female / pop,
        mx_lip = n_lip / pop,
        mx_lymphoid = n_lymphoid / pop,
        mx_male = n_male / pop,
        mx_meso = n_meso / pop,
        mx_resp = n_resp / pop,
        mx_skin = n_skin / pop,
        mx_thyroid = n_thyroid / pop,
        mx_urinary = n_urinary / pop,
        mx_other_cancer = n_other_cancer / pop,
        mx_other_cancer_grouped = n_other_cancer_grouped / pop,
        mx_multi = n_multi / pop,
        mx_other = n_other / pop
    )

## Setting up the simulation ----
### Need a scaled version of fertility df with overdispersion ----
scaled_fx <- create_scaled_fx_df(fx.US, lt.US) |>
    left_join(theta.bth |>
                  dplyr::select(-c(mu, var)), by = c("race_eth", "sex", "age"))

### Create PHI matrices for male and female parents ----
phi_matrices <- list("female" = get_phi(unique(lt.US$age), "female"),
                     "male" = get_phi(unique(lt.US$age), "male"))

## Set up the parallization ----
### Create grid to iterate over ----
param_grid <- expand_grid(
    fertility_multipliers = c(1, seq(.75, .95, .05)),
    cause_of_death = sort(unique(scaled_fx$cause)),
    race_eth = c("total", "black", "hispanic", "white")
) |>
    mutate(random_seed = round(runif(n(), 0, 1) * 1000000)) |>
    arrange(fertility_multipliers, cause_of_death, race_eth)

### Make file paths now
param_grid <- param_grid |>
    mutate(f_path = case_when(
        fertility_multipliers == 1 ~ here(
            "output",
            "simulations",
            race_eth,
            sprintf(
                "%s_%s_scenario%i_sims%05d.RDS",
                race_eth,
                cause_of_death,
                round(fertility_multipliers * 100),
                N_SIM
            )
        ),
        fertility_multipliers != 1 ~ here(
            "output",
            "simulations",
            race_eth,
            "sensitivity",
            sprintf(
                "%s_%s_scenario%i_sims%05d.RDS",
                race_eth,
                cause_of_death,
                round(fertility_multipliers * 100),
                N_SIM
            )
        )
    ))

param_grid <- param_grid |>
    filter(fertility_multipliers %in% c(1))

for (i in 1:NROW(param_grid)) {
    race_eth <- param_grid$race_eth[i]
    cod <- param_grid$cause_of_death[i]
    sc.fx <- param_grid$fertility_multipliers[i]
    f_path <- param_grid$f_path[i]
    random_seed <- param_grid$random_seed[i]
    
    if (!file_exists(f_path) | FORCE_REFRESH) {
        fs::dir_create(dirname(f_path))
        
        set.seed(random_seed)
        
        future::plan(future::multisession(workers = N_CORES))
        
        timer <- tictoc::tic()
        
        ## NOTE: calculate_parental_deaths() will subset to just our years
        ## of interest *but* will keep all child ages so we will need to
        ## subset ages later. 
        female_results <- future.apply::future_replicate(N_SIM, {
            calculate_parental_deaths(
                lifetable_df = lt.US,
                scaled_fertility_df = scaled_fx,
                phi_matrix = phi_matrices,
                years_of_interest = 1999:2020,
                parent_sex = "female",
                race_x = race_eth,
                cod_x = cod,
                fx_multiplier = sc.fx
            )
        })
        
        male_results <- future.apply::future_replicate(N_SIM, {
            calculate_parental_deaths(
                lifetable_df = lt.US,
                scaled_fertility_df = scaled_fx,
                phi_matrix = phi_matrices,
                years_of_interest = 1999:2020,
                parent_sex = "male",
                race_x = race_eth,
                cod_x = cod,
                fx_multiplier = sc.fx
            )
        })
        
        tictoc::toc(quiet = TRUE)
        
        ## Close the parallel backend
        doParallel::stopImplicitCluster()
        closeAllConnections()
        
        metainfo <- tibble(
            f_name = basename(f_path),
            random_seed = random_seed,
            datetime = Sys.time(),
            race_eth = race_eth,
            cause_of_death = cod,
            fertility_scenario = sc.fx,
            n_boots = N_SIM,
            n_cores = N_CORES,
            clocktime = as.numeric(timer)
        )
        
        results <- list(metainfo = metainfo,
                        female = female_results,
                        male = male_results)
        
        ## Parallelized saveRDS
        if (.Platform$OS.type == "unix") {
            saveRDS_xz(results, f_path, threads = N_CORES)
        } else {
            saveRDS(results, f_path, compress = "xz")
        }
        
        ## Clean up
        rm(
            results,
            male_results,
            female_results,
            metainfo,
            race_eth,
            cod,
            sc.fx,
            f_path,
            random_seed
        )
        gc()
    }
}
