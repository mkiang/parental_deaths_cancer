## Imports ----
library(MASS)
library(tidyverse)
library(here)
library(abind)
library(fs)
library(foreach)
library(doParallel)
source(here("code", "utils.R"))

## CONSTANTS ----
N_CORES <- 8
AGE_IX <- c(1:18, 87:104) # 1:18 females. 87:104 males

## Data ----
### Mortality data
lt.US <- readRDS(here("data", "national_lifetable.RDS"))

### Overdispersion parameters for NBin in MC
theta.dth <- readRDS(here("data", "overdispersion_deaths.RDS")) |>
    dplyr::select(-c(mu, var))

### Fertility data (new input for male fertility)
fx.US <- readRDS(here("data", "fx_US.rds"))

## Reshape life tables ----
lt.US.long <- lt.US |>
    filter(year >= 1999,
           race_eth %in% c("black", "white", "hispanic", "total")) |>
    mutate(n_other = n_deaths - (n_cancer)) |>
    dplyr::select(year, race_eth, sex, age, all = n_deaths, n_other, n_cancer:n_multi) |>
    pivot_longer(all:n_multi, names_to = "cause", values_to = "dth") |>
    mutate(cause = gsub("n_", "", cause)) |>
    left_join(theta.dth, by = c("race_eth", "sex", "age"))

## Loop through the race/ethnicity and cause-specific results ----
## NOTE: We are going to do each cause separately and this appears to be 
## less computationally efficient because we need to read in the "all cause"
## file every time, but because this is much more memory efficient, we are
## able to run it across more cores simultaneously. I think overall this
## results in a faster approach but we can revisit this at some point.

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
    mutate(sensitivity = grepl("sensitivity", f_path, fixed = TRUE)) |>
    filter(!sensitivity)

all_grid <- param_grid |> 
    filter(cause == "all")

### Do the actual extractions for each cause, noting that we will
### save "all causes" multiple times. 
doParallel::registerDoParallel(cores = N_CORES)
foreach::foreach(i = 1:NROW(param_grid), .inorder = FALSE) %dopar% {
    COD <- param_grid$cause[i]
    R <- param_grid$race_eth[i]
    F_PATH <- param_grid$f_path[i]
    F_ALL <- all_grid |> filter(race_eth == R) |> pull(f_path)
    
    f_loss_by_cause <- here("output",
                            "summaries",
                            sprintf("df_n_loss_by_cause_%s_%s.RDS", COD, R))
    
    f_incidence_by_cause <- here(
        "output",
        "summaries",
        sprintf("df_incidence_parental_loss_by_cause_%s_%s.RDS", COD, R)
    )
    
    f_incidence_by_sex <- here(
        "output",
        "summaries",
        sprintf("df_incidence_parental_loss_by_cause_sex_%s_%s.RDS", COD, R)
    )
    
    f_loss_by_sex <- here("output",
                          "summaries",
                          sprintf("df_n_loss_by_cause_sex_%s_%s.RDS", COD, R))
    
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
        
        ## Extract constants we will need later ----
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
            ungroup()
        
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
            ungroup()
        
        # Orphanhood incidence of at least one parent
        df.inc.orphans.at.least <- as.data.frame.table(incidence) |>
            rename("p" = Freq) |>
            as_tibble() |>
            mutate(year = as.character(year) |> as.numeric(), race = RACE_X) |>
            pivot_wider(names_from = sex, values_from = p) |>
            mutate(p = 1 - ((1 - female) * (1 - male))) |>
            group_by(race, year, cause) |>
            group_modify(get_quantiles, var = "p") |>
            ungroup()
        
        # Number of orphans (using incidence at least from above)
        df.n.orphans.at.least <- as.data.frame.table(incidence) |>
            rename("p" = Freq) |>
            as_tibble() |>
            mutate(year = as.character(year) |> as.numeric(), race = RACE_X) |>
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
            ungroup()
        
        # Save ----
        dir_create(here("output", "summaries"))
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

## Gather summaries together ----
### Get file lists ----
f_incidence <- dir_ls(here("output", "summaries"),
                      recurse = FALSE,
                      regexp = "\\<df_incidence")
f_number <- dir_ls(here("output", "summaries"),
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
        here("data", "summarized_results.RDS"),
        compress = "xz")
