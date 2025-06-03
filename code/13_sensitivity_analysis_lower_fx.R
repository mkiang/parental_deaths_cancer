## Imports ----
library(MASS)
library(tidyverse)
library(here)
library(fs)
library(abind)
library(doParallel)
library(future)
library(future.apply)
source(here("code", "utils.R"))

## CONSTANTS ----
N_SIM <- 2000
N_CORES <- 16
FORCE_REFRESH <- FALSE 
AGE_IX <- c(1:18, 87:104) # 1:18 females. 87:104 males

## Data ----
lt.US <- readRDS(here("data", "national_lifetable.RDS"))

### Overdispersion parameters for NBin in MC
theta.bth <- readRDS(here("data", "overdispersion_births.RDS"))
theta.dth <- readRDS(here("data", "overdispersion_deaths.RDS")) |>
    dplyr::select(-c(mu, var))

### Fertility data (new input for male fertility rates)
fx.US <- readRDS(here("data", "fx_US.rds"))

### Reshape life tables ----
lt.US.long <- lt.US |>
    filter(year >= 1999,
           race_eth %in% c("black", "white", "hispanic", "total")) |>
    mutate(n_other = n_deaths - (n_cancer)) |>
    dplyr::select(year, race_eth, sex, age, all = n_deaths, n_other, n_cancer:n_multi) |>
    pivot_longer(all:n_multi, names_to = "cause", values_to = "dth") |>
    mutate(cause = gsub("n_", "", cause)) |>
    left_join(theta.dth, by = c("race_eth", "sex", "age"))

### Clean up original (wide) life table 
lt.US <- lt.US |> 
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

## Part 1. Run the simulation ----
### Setting up the simulation ----
#### Need a scaled version of fertility df with overdispersion ----
scaled_fx <- create_scaled_fx_df(fx.US, lt.US) |>
    left_join(theta.bth |>
                  dplyr::select(-c(mu, var)), by = c("race_eth", "sex", "age"))

#### Create PHI matrices for male and female parents ----
phi_matrices <- list("female" = get_phi(unique(lt.US$age), "female"),
                     "male" = get_phi(unique(lt.US$age), "male"))

#### Set up the parallelization ----
##### Create grid to iterate over ----
param_grid <- expand_grid(
    fertility_multipliers = seq(.8, .95, .05),
    cause_of_death = sort(unique(scaled_fx$cause)),
    race_eth = c("total", "black", "hispanic", "white")
) |>
    mutate(random_seed = round(runif(n(), 0, 1) * 1000000)) |>
    arrange(fertility_multipliers, cause_of_death, race_eth)

#### Make file paths now
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
    filter(cause_of_death %in% c("cancer"))

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

## Part 2. Extract parameters of interest from simulations ----
### Create grid to loop over ----
f_list <- dir_ls(
    here("output", "simulations"),
    recurse = TRUE,
    type = "file",
    glob = "*.RDS"
)

param_grid <- tibble(
    race_eth = stringr::str_split_i(basename(f_list), "_", 1),
    cause = case_when(
        stringr::str_split_i(basename(f_list), "_", 4) == "grouped" ~ "other_cancer_grouped",
        stringr::str_split_i(basename(f_list), "_", 3) == "cancer" ~ "other_cancer",
        TRUE ~ stringr::str_split_i(basename(f_list), "_", 2)
    ),
    f_path = f_list
) |>
    mutate(sensitivity = grepl("sensitivity", f_path, fixed = TRUE))

all_grid <- param_grid |> 
    filter(cause == "all")

param_grid <- param_grid |>
    filter(sensitivity)

doParallel::registerDoParallel(cores = N_CORES)
foreach::foreach(i = 1:NROW(param_grid), .inorder = FALSE) %dopar% {
    COD <- param_grid$cause[i]
    R <- param_grid$race_eth[i]
    F_PATH <- param_grid$f_path[i]
    F_ALL <- all_grid |> filter(race_eth == R) |> pull(f_path)
    SC_FX <- gsub("(.*scenario)([0-9]*)(.*)", "\\2", basename(F_PATH))
    SC_FX <- as.numeric(SC_FX) / 100
    
    f_loss_by_cause <- here(
        "output",
        "summaries_sensitivity",
        sprintf("df_n_loss_by_cause_%s_%s_scenario%i.RDS", COD, R, SC_FX *
                    100)
    )
    
    f_incidence_by_cause <- here(
        "output",
        "summaries_sensitivity",
        sprintf(
            "df_incidence_parental_loss_by_cause_%s_%s_scenario%i.RDS",
            COD,
            R,
            SC_FX * 100
        )
    )
    
    f_incidence_by_sex <- here(
        "output",
        "summaries_sensitivity",
        sprintf(
            "df_incidence_parental_loss_by_cause_sex_%s_%s_scenario%i.RDS",
            COD,
            R,
            SC_FX * 100
        )
    )
    
    f_loss_by_sex <- here(
        "output",
        "summaries_sensitivity",
        sprintf(
            "df_n_loss_by_cause_sex_%s_%s_scenario%i.RDS",
            COD,
            R,
            SC_FX * 100
        )
    )
    
    if (!all(file_exists(c(f_loss_by_cause,
                           f_incidence_by_cause, 
                           f_incidence_by_sex,
                           f_loss_by_sex)))) {
        
        ## Reshape the "all cause" array
        temp_x <- readRDS(F_ALL)
        temp_x$female <- temp_x$female[AGE_IX, , , ]
        temp_x$male <- temp_x$male[AGE_IX, , , ]
        
        ## Combined female/male arrays to match Ben's code
        binded_array_all <- abind::abind(temp_x$female, 
                                         temp_x$male,
                                         along = 2.5)
        dimnames(binded_array_all) <- list(
            "child age" = dimnames(binded_array_all)[[1]],
            "focal age" = dimnames(binded_array_all)[[2]],
            "focal sex" = c("female", "male"),
            year = dimnames(binded_array_all)[[4]],
            simulation = 1:NROW(binded_array_all[1, 1, 1, 1, ])
        )
        meta_df <- temp_x$metainfo
        
        rm(temp_x)
        gc()
        
        ## Reshape the COD array
        temp_x <- readRDS(F_PATH)
        temp_x$female <- temp_x$female[AGE_IX, , , ]
        temp_x$male <- temp_x$male[AGE_IX, , , ]
        
        ## Combined female/male arrays to match Ben's code
        binded_array_cod <- abind::abind(temp_x$female,
                                         temp_x$male, 
                                         along = 2.5)
        dimnames(binded_array_cod) <- list(
            "child age" = dimnames(binded_array_cod)[[1]],
            "focal age" = dimnames(binded_array_cod)[[2]],
            "focal sex" = c("female", "male"),
            year = dimnames(binded_array_cod)[[4]],
            simulation = 1:NROW(binded_array_cod[1, 1, 1, 1, ])
        )
        
        COD <- temp_x$metainfo$cause_of_death
        
        rm(temp_x)
        gc()
        
        #### Extract constants we will need later ----
        N_BOOTS <- meta_df$n_boots[1]
        RACE_X <- meta_df$race_eth[1]
        N_YEARS <- NROW(dimnames(binded_array_cod)[["year"]])
        YEARS <- as.numeric(unique(dimnames(binded_array_cod)[["year"]]))
        N_SEX <- NROW(dimnames(binded_array_cod)[["focal sex"]])
        SEXES <- unique(dimnames(binded_array_cod)[["focal sex"]])
        
        ## Get the number of living children by age of parent and year
        a_child_t_all <- apply(binded_array_all, 2:5, sum)
        a_child_t_cod <- apply(binded_array_cod, 2:5, sum)
        
        rm(binded_array_all, binded_array_cod)
        gc()
        
        ### Make empty containers for quantities of interest ----
        ### Container for the number of maternal/paternal orphans
        n_orphans <- array(
            NA,
            dim = c(N_SEX, N_YEARS, 1, N_BOOTS),
            dimnames = list(
                "sex" = SEXES,
                "year" = YEARS,
                "cause" = COD,
                "sim" = 1:N_BOOTS
            )
        )
        
        ### Container for the total number of children
        n_children <- array(
            NA,
            dim = c(N_SEX, N_YEARS, 1, N_BOOTS),
            dimnames = list(
                "sex" = SEXES,
                "year" = YEARS,
                "cause" = COD,
                "sim" = 1:N_BOOTS
            )
        )
        
        # Loop over simulations
        for (j in 1:N_BOOTS) {
            ## Simulate death counts
            lt.US.i <- lt.US.long |>
                filter(race_eth == RACE_X) |>
                # Generate deaths at each MC iteration
                mutate(D.sim = map2_int(dth, theta, ~ suppressWarnings(rnegbin(1, .x, .y))))
            
            ## Iterate through the array
            for (s in SEXES) {
                for (y in 1:N_YEARS) {
                    # Av. nber of child per individual age x for whole pop
                    child_per_ind <- a_child_t_all[, s, y, j]
                    
                    # Population counts
                    N_t <- lt.US |>
                        filter(sex == s, year == 1998 + y, race_eth == RACE_X) |>
                        pull(pop)
                    
                    
                    # Average number of child per individual age x if parent die from cause i
                    child_per_ind_cod <- a_child_t_cod[, s, y, j]
                    
                    # Death counts by cause of death
                    D_t <- lt.US.i |>
                        filter(sex == s,
                               year == 1998 + y,
                               race_eth == RACE_X,
                               cause == COD) |>
                        pull(D.sim)
                    
                    n_orphans[s, y, COD, j] <- child_per_ind_cod %*% D_t
                    n_children[s, y, COD, j] <- child_per_ind %*% N_t
                    
                }
            }
        }
        
        ### Clean up
        rm(a_child_t_all, a_child_t_cod)
        gc()
        
        # Number of orphans by parent sex and cause
        df.n.orphans <- as.data.frame.table(n_orphans) |>
            rename("n" = Freq) |>
            as_tibble() |> 
            mutate(
                # Make categories coherent with code for figures
                year = as.character(year) |> as.numeric(),
                sex = ifelse(sex == "female", "mother", "father"),
                race = RACE_X
            ) |>
            group_by(sex, race, year, cause) |>
            group_modify(get_quantiles, var = "n") |>
            ungroup() |> 
            mutate(scenario = SC_FX)
        
        # Orphanhood incidence by parent sex
        incidence <- n_orphans / n_children
        
        df.inc.orphans <- as.data.frame.table(incidence) |>
            rename("p" = Freq) |>
            as_tibble() |>
            # Make categories coherent with code for figures
            mutate(
                year = as.character(year) |> as.numeric(),
                sex = ifelse(sex == "female", "mother", "father"),
                race = RACE_X
            ) |>
            group_by(sex, race, year, cause) |>
            group_modify(get_quantiles, var = "p") |>
            ungroup() |> 
            mutate(scenario = SC_FX)
        
        # Orphanhood incidence of at least one parent
        df.inc.orphans.at.least <- as.data.frame.table(incidence) |>
            rename("p" = Freq) |>
            as_tibble() |>
            mutate(year = as.character(year) |> as.numeric(),
                   race = RACE_X) |>
            pivot_wider(names_from = sex, values_from = p) |>
            mutate(p = 1 - ((1 - female) * (1 - male))) |>
            group_by(race, year, cause) |>
            group_modify(get_quantiles, var = "p") |>
            ungroup() |> 
            mutate(scenario = SC_FX)
        
        # Number of orphans (using incidence at least from above)
        df.n.orphans.at.least <- as.data.frame.table(incidence) |>
            rename("p" = Freq) |>
            as_tibble() |>
            mutate(year = as.character(year) |> as.numeric(),
                   race = RACE_X,
                   scenario = SC_FX) |>
            pivot_wider(names_from = sex, values_from = p) |>
            # Prob losing at least one parent
            mutate(p = 1 - ((1 - female) * (1 - male))) |>
            left_join(
                # Sum over sex and ages
                lt.US |>
                    rename("race" = race_eth) |>
                    filter(age < 18) |>
                    summarise(.by = c(year, race),
                              pop = sum(pop)),
                by = c("year", "race")
            ) |>
            # Number of children that lose at least a parent
            mutate(n = p * pop) |>
            group_by(race, year, cause) |>
            group_modify(get_quantiles, var = "n") |>
            ungroup() |> 
            mutate(scenario = SC_FX)
        
        # Save ----
        dir_create(here("output", "summaries_sensitivity"))
        saveRDS(df.n.orphans.at.least, f_loss_by_cause, compress = "xz")
        saveRDS(df.inc.orphans.at.least, f_incidence_by_cause, compress = "xz")
        saveRDS(df.inc.orphans, f_incidence_by_sex, compress = "xz")
        saveRDS(df.n.orphans, f_loss_by_sex, compress = "xz")
        
        ## Clean up
        rm(df.n.orphans.at.least,
           df.inc.orphans.at.least,
           df.inc.orphans,
           df.n.orphans)
        gc()
    }
}

## Close the parallel backend
doParallel::stopImplicitCluster()
closeAllConnections()

## Part 3. Gather summaries together ----
### Get file lists ----
f_incidence <- dir_ls(here("output", "summaries_sensitivity"),
                      recurse = FALSE,
                      regexp = "\\<df_incidence")
f_number <- dir_ls(here("output", "summaries_sensitivity"),
                   recurse = FALSE,
                   regexp = "\\<df_n_loss")

### Read in each independently ----
df_incidence <- map_dfr(.x = f_incidence, .f = ~ {
    temp_x <- readRDS(.x) |>
        mutate(metric = "percent")
    if (!has_name(temp_x, "sex")) {
        temp_x |>
            mutate(sex = "both")
    } else {
        temp_x
    }
})

df_number <- map_dfr(.x = f_number, .f = ~ {
    temp_x <- readRDS(.x) |>
        mutate(metric = "number")
    if (!has_name(temp_x, "sex")) {
        temp_x |>
            mutate(sex = "both")
    } else {
        temp_x
    }
})

### Save ----
saveRDS(bind_rows(df_incidence, df_number),
        here("data", "summarized_results_lowerfx.RDS"),
        compress = "xz")
