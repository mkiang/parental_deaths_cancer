## Imports ----
library(tidyverse)
library(parallel)
library(pracma)
library(Matrix)

## Death certificate helper functions ----
create_big_race <- function(df) {
    df %>% 
        mutate(race_bridged = case_when(
            race %in% c(1) ~ "white",
            race %in% c(2) ~ "black",
            race %in% c(3) ~ "aian",
            race %in% c(4:7, 18, 28, 38, 48, 58, 68, 78) ~ "api",
            TRUE ~ NA_character_
        ))
}

recode_age <- function(df) {
    df %>% 
        mutate(age_years = case_when(
            age %in% c(9999, 1999, 2999, 4999, 5999, 6999) ~ NA_real_, 
            age > 1999 ~ 0,
            age < 1999 ~ age - 1000,
            TRUE ~ NA_real_
        ))
}

categorize_race <- function(df) {
    df |>
        mutate(race_cat = factor(
            race,
            levels = c("hispanic", "black", "white", "total"),
            labels = c(
                "Hispanic",
                "Non-Hispanic Black",
                "Non-Hispanic White",
                "Total"
            ),
            ordered = TRUE
        )) |>
        mutate(race_cat_rev = factor(
            race,
            levels = rev(c("hispanic", "black", "white", "total")),
            labels = rev(c(
                "Hispanic",
                "Non-Hispanic Black",
                "Non-Hispanic White",
                "Total"
            )),
            ordered = TRUE
        )) |>
        arrange(race_cat)
}

categorize_deaths <- function(df) {
    df |> 
        mutate(death_cat = factor(
            cause, 
            levels = c(
                "all",
                "cancer",
                "breast",
                "digestive",
                "eye",
                "female",
                "lymphoid",
                "resp",
                "urinary",
                "other_cancer_grouped",
                "other", 
                "bone",
                "lip",
                "male",
                "meso",
                "thyroid",
                "multi",
                "other_cancer",
                "skin"
            ),
            labels = c(
                "All causes",
                "All cancers",
                "Breast cancer",
                "Digestive cancer",
                "Eye/Brain cancer",
                "Female organs cancer",
                "Lymphoid cancer",
                "Respiratory cancer",
                "Urinary cancer",
                "All other cancers",
                "All other causes",
                "Bone cancer",
                "Lip cancer",
                "Male organs cancer",
                "Mesothelioma",
                "Thyroid cancer",
                "Multisite cancer",
                "Other cancer",
                "Skin cancer"
            ),
            ordered = TRUE
        )) |> 
        mutate(death_cat_rev = factor(
            cause, 
            levels = rev(c(
                "all",
                "cancer",
                "breast",
                "digestive",
                "eye",
                "female",
                "lymphoid",
                "resp",
                "urinary",
                "other_cancer_grouped",
                "other", 
                "bone",
                "lip",
                "male",
                "meso",
                "thyroid",
                "multi",
                "other_cancer",
                "skin"
            )),
            labels = rev(c(
                "All causes",
                "All cancers",
                "Breast cancer",
                "Digestive cancer",
                "Eye/Brain cancer",
                "Female organs cancer",
                "Lymphoid cancer",
                "Respiratory cancer",
                "Urinary cancer",
                "All other cancers",
                "All other causes",
                "Bone cancer",
                "Lip cancer",
                "Male organs cancer",
                "Mesothelioma",
                "Thyroid cancer",
                "Multisite cancer",
                "Other cancer",
                "Skin cancer"
            )),
            ordered = TRUE
        )) |> 
        mutate(death_cat_alpha_rev = factor(
            cause, 
            levels = rev(c(
                "all",
                "cancer",
                "other", 
                "breast",
                "bone",
                "digestive",
                "eye",
                "female",
                "lymphoid",
                "resp",
                "urinary",
                "lip",
                "male",
                "meso",
                "skin",
                "thyroid",
                "multi",
                # "other_cancer",
                "other_cancer_grouped"
            )),
            labels = rev(c(
                "All causes",
                "All cancers",
                "All other causes",
                "Breast",
                "Bone",
                "Digestive organs",
                "Eye, Brain, and CNS",
                "Female genital organs",
                "Lymphoid",
                "Respiratory",
                "Urinary",
                "Lip, Oral cavity, and Pharynx",
                "Male genital organs",
                "Mesothelial",
                "Melanoma", 
                "Thyroid",
                "Multisite",
                # "Other cancer",
                "All other cancers"
            )),
            ordered = TRUE
        ))
}

## Life table helper functions ----
get_U <- function(s, y, r, lifeTable) {
    ## Inputs:
    ##  s: sex (male or female)
    ##  y: year
    ##  r: race_eth
    ##  lifeTable: a lifeTable dataframe
    ##
    ## Returns: a matrix U where the subdiagonal contains the survival
    ## probabilities for this sex-year-race combination
    
    ## Extract survival probs
    px <- lifeTable %>%
        filter(sex == s, year == y, race_eth == r) %>%
        pull(px)
    
    ## Dim of A
    omega <- NROW(px)
    
    ## Creation of U
    U <- matrix(0, nrow = omega, ncol = omega)
    
    ## Put the survival prob (except the last age) on subdiagonal
    # ASSUMPTION: Everybody die at last age
    U[row(U) - col(U) == 1] <- head(px, -1)
    
    return(U)
}

get_U_2sex <- function(lifeTable, y, r) {
    ## Create a block diagonal matrix of U for both sexes
    
    ## U for both sexes
    U <- lapply(c("female", "male"), get_U, y, r, lifeTable)
    
    U_2sex <- bdiag(U)
    
    return(as.matrix(U_2sex))
}

get_F <- function(s, y, r, ages, asfr) {
    ## Inputs:
    ##  s: sex (male or female)
    ##  y: year
    ##  r: race_eth
    ##  asfr: fertility dataframe
    ## 
    ## Make a fertility matrix F from fx dataframe
    
    omega <- length(ages)
    
    ## Extract asfr
    fx <- asfr %>%
        filter(sex == s,
               year == y,
               race_eth == r) %>%
        pull(fx)
    
    ## Reproduction ages present in asfr
    ages.repro <- asfr %>%
        filter(sex == s,
               year == y,
               race_eth == r) %>%
        pull(age)
    
    ## Creation of F
    Fdf <- matrix(0,
                 nrow = omega,
                 ncol = omega)
    ## ASFR on 1st row
    Fdf[1, (ages.repro + 1)] <- fx
    
    return(Fdf)
}

## Construct F for both sexes
get_F_2sex <- function(birth_female = 1 / 2.04, ages, asfr, y, r) {
    ## U for both sexes
    Fdf <- sapply(c("female", "male"), function(s) {
        Fdf <- get_F(s, y, r, ages, asfr)
        
        return(Fdf)
    },
    simplify = FALSE,
    USE.NAMES = TRUE)
    
    # Create block matrix when focal is female or male
    block_f <- cbind(birth_female * Fdf[["female"]], 
                     birth_female * Fdf[["male"]])
    block_m <- cbind((1 - birth_female) * Fdf[["female"]],
                     (1 - birth_female) * Fdf[["male"]])
    
    F_2sex <- rbind(block_f, block_m)
    
    return(as.matrix(F_2sex))
}

## Construct phi depending on Focal sex (ie mother or father)
get_phi <- function(ages, sex = "female") {
        
        omega <- length(ages)
        
        phi <- diag(rep(1, omega))
        zero <- diag(rep(0, omega))
        
        if (sex == "female") {
                PHI <- rbind(phi, zero)
        } else {
                PHI <- rbind(zero, phi)
        }
        return(PHI)
}

intra_fx <- function(df, sex) {
    ## From https://github.com/benjisamschlu/parental_deaths
    ## Ben's function to interpolate age-specific fertility rates using a
    ## Hermitte cubic. Takes as input a dataframe of year, race_eth, age, and
    ## fertility rate (raw data from NCHS) and returns as output 
    
    if (sex == "female") {
        age.bd = c(9, 55)
    } else {
        age.bd = c(9, 65)
    }
    
    fx.int <- interp1(c(age.bd[1], df$age.m, age.bd[2]),
                      # Set fx=0 at ages 9 and 55 yo.
                      c(0, df$fx, 0),
                      age.bd[1]:age.bd[2],
                      method = "cubic")
    
    fx.int <- tibble(age = age.bd[1]:age.bd[2], fx = fx.int) %>%
        # Make sure interpolation are >0
        mutate(fx = ifelse(fx < 0, 0, fx), sex = sex)
    return(fx.int)
}

create_exposure_df <- function(lifetable_df) {
    ## The Monte Carlo simulations need scaled population counts by cause
    ## of death. This just takes the life table and returns the proportion
    ## of deaths for each cause by year, race_eth, sex, and age.
    ##
    ## NOTE: This is not generalized intentionally. For every new cause, you
    ## must change this manually. However, you should only change the deaths
    ## in between the columns "other" and "all" if you want the column 
    ## selection to work in downstream functions. 
    ## 
    ## Note that we also add a scalar of 1 for total population. 
    
    lifetable_df |>
        # No cause-specific death counts before 1999 because of ICD change
        dplyr::filter(year >= 1999) |>
        dplyr::group_by(year, race_eth, sex) |>
        dplyr::mutate(
            cancer = sum(n_cancer) / sum(pop),
            bone = sum(n_bone) / sum(pop),
            breast = sum(n_breast) / sum(pop),
            digestive = sum(n_digestive) / sum(pop),
            eye = sum(n_eye) / sum(pop),
            female = sum(n_female) / sum(pop),
            lip = sum(n_lip) / sum(pop),
            lymphoid = sum(n_lymphoid) / sum(pop),
            male = sum(n_male) / sum(pop),
            meso = sum(n_meso) / sum(pop),
            resp = sum(n_resp) / sum(pop),
            skin = sum(n_skin) / sum(pop),
            thyroid = sum(n_thyroid) / sum(pop),
            urinary = sum(n_urinary) / sum(pop),
            other_cancer = sum(n_other_cancer) / sum(pop),
            other_cancer_grouped = sum(n_other_cancer_grouped) / sum(pop),
            multi = sum(n_multi) / sum(pop),
            other = sum(n_other) / sum(pop), 
            all = 1
        ) |>
        dplyr::ungroup() |>
        dplyr::select(year, race_eth, sex, age, cancer:all)
}

create_scaled_fx_df <- function(fx_df, lifetable_df) {
    ## Get scalars
    scalar_df <- create_exposure_df(lifetable_df)
    
    ## Need to differentiate pre- and post-1999
    fx.US.post.1999 <- fx_df |>
        filter(year >= 1999) |>
        # Add scaling factors back in
        left_join(scalar_df, by = c("year", "race_eth", "sex", "age"))
    
    # Assume scaling factors of 1999 for years < 1999
    fx.US.pre.1999 <- fx_df |>
        filter(year < 1999) |>
        left_join(
            scalar_df |>
                filter(year == 1999) |>
                dplyr::select(-year),
            by = c("race_eth", "sex", "age")
        )
    
    # Merge both
    fx.US <- fx.US.pre.1999 |>
        bind_rows(fx.US.post.1999) |>
        pivot_longer(cancer:all, names_to = "cause", values_to = "scalar.exp")
    
    
    # Add population counts to fx
    fx.US <- fx.US |>
        left_join(
            lifetable_df |>
                dplyr::select(year, race_eth, sex, age, pop),
            by = c("year", "race_eth", "sex", "age")
        )
    
    # Get exposures by cause of death by applying the scalars
    # mu is mean parameter of negbin in simulations
    # bounds are only meaningful when fertility is modeled,
    # else assume accurately captured.
    fx.US <- fx.US |>
        mutate(
            exp = scalar.exp * pop,
            mu_nb = fx * exp,
            mu_nb.l95 = ifelse(is.na(fx.l95), mu_nb, fx.l95 * exp),
            mu_nb.u95 = ifelse(is.na(fx.u95), mu_nb, fx.u95 * exp)
        )
    
    fx.US
}

calculate_parental_deaths <- function(lifetable_df,
                                      scaled_fertility_df,
                                      phi_matrix, 
                                      years_of_interest = 1999:2020,
                                      parent_sex = "female", 
                                      race_x = "total",
                                      cod_x = "all", 
                                      fx_multiplier = 1) {
    ## CONSTANTS ----
    AGES <- unique(lifetable_df$age)
    N_AGES <- length(AGES)
    ALL_YEARS <- unique(lifetable_df$year)
    N_YEARS <- length(ALL_YEARS)
    
    ## Create an array ----
    # Container for the average number of living children over 
    # age, year, sex, and race/ethnicity of their parent (=Focal)
    a_x_t <- array(
        NA,
        dim = c(
            2 * N_AGES, 
            N_AGES, 
            N_YEARS + 1),
        dimnames = list(
            "children" = 1:(2 * N_AGES),
            "focal age" = AGES,
            "year" = c("stable", ALL_YEARS)
        )
    )
    
    # Boundary conditions
    # No child at birth in every year
    a_x_t[, 1, ] <- 0
    
    ## Create simulated fertility 
    fx.US.i <- scaled_fertility_df |>
        filter(cause == cod_x, 
               race_eth == race_x) |>
        dplyr::select(
            race_eth,
            year, 
            age,
            sex, 
            theta,
            exp, 
            mu_nb.l95, 
            mu_nb.u95, 
            fx) |>
        mutate(
            mu_nb.sim = map2_dbl(mu_nb.l95, mu_nb.u95, ~ qrunif(.x, .y)),
            B.sim = map2_int(mu_nb.sim, theta, ~ qrnegbin(sc.fx * .x, .y)),
            fx.sim = ifelse(exp == 0, 0, B.sim / exp)
        ) |>
        dplyr::select(race_eth, year, age, sex, fx = fx.sim)
    
    # Get a(x,0) for all x assuming stable pop in the year 
    # 1990 (Focal could be any age in the year 1990)
    U_2sex <- get_U_2sex(lifetable_df, min(ALL_YEARS), race_x)
    F_2sex <- get_F_2sex(birth_female = 1 / 2.04, 
                         AGES, 
                         fx.US.i,
                         min(ALL_YEARS), 
                         race_x)
    
    ## Loop over ages 
    for (i in 2:N_AGES) {
        a_x_t[, i, 1] <- U_2sex %*% a_x_t[, i - 1, 1] +
            F_2sex %*% phi_matrix[[parent_sex]][, i - 1]
    }
    
    # Start loop over the years
    for (i in 1:N_YEARS) {
        # Projection of children
        U_2sex <- get_U_2sex(lifetable_df, ALL_YEARS[i], race_x)
        F_2sex <- get_F_2sex(birth_female = 1 / 2.04, 
                             AGES, 
                             fx.US.i, 
                             ALL_YEARS[i], 
                             race_x)
        
        a_x_t[, 2:N_AGES, (i + 1)] <-
            U_2sex %*% a_x_t[, 1:(N_AGES - 1), i] +
            F_2sex %*% phi_matrix[[parent_sex]][, 1:(N_AGES - 1)]
    }
    
    # Do not consider years before 1999 (no death by cause)
    # and focus on children (both sexes) < 18 yo
    # a_x_t[c(1:18, 87:104), , which(ALL_YEARS %in% years_of_interest) + 1]
    ## Create year indices: 
    YEAR_IX <- which(dimnames(a_x_t)$year %in% years_of_interest)
    a_x_t[, , YEAR_IX]
}


## Random number generators ----
qrunif <- function(.x, .y) {
    suppressWarnings(runif(1, .x, .y))
}

qrnegbin <- function(.x, .y) {
    suppressWarnings(rnegbin(1, .x, .y))
}

# Function to obtain quantiles of quantities of interest
get_quantiles <- function(x, var, ...) {
    
    qs <- quantile(x[[var]],
                   probs = c(0.025, 0.5, 0.975))
    out <- tibble(
        median = qs[2],
        lower95 = qs[1],
        upper95 = qs[3]
    )
}

## RDS helper functions ----
saveRDS_xz <- function(object, file, threads = parallel::detectCores() - 1) {
    ## https://gist.github.com/ShanSabri/b1bdf0951efa0dfee0edeb5509f87e88
    ## If xz version 5 or higher is installed, parallelize
    if (any(grepl("(XZ Utils) 5.", system("xz -V", intern = TRUE), fixed = TRUE))) {
        con <- pipe(paste0("xz -T", threads, " > ", file), "wb")
        saveRDS(object, file = con)
        close(con)
    } else {
        saveRDS(object, file = file, compress = "xz")
    }
}

