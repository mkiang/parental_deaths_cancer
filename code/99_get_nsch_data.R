## 99_get_nsch_data.R ----
##
## NOTE: Code comes from Ben's github repo based on our JAMA paper.
## See https://github.com/benjisamschlu/parental_deaths for original.

## Imports ----
library(tidyverse)
library(here)
library(fs)
library(rvest)

### Scrape NSCH data if we don't already have it ----
if (!fs::file_exists(here::here("data", "nsch_estimates_allcause_parental_death.RDS"))) {
    ## NSCH % about "Parent or guardian died"
    # Scrape html table
    # Url link
    urls_nat <- list(
        "2016" = "https://nschdata.org/browse/survey/results?q=4786&r=1",
        "2017" = "https://nschdata.org/browse/survey/results?q=6763&r=1",
        "2018" = "https://nschdata.org/browse/survey/results?q=7445&r=1",
        "2019" = "https://nschdata.org/browse/survey/results?q=8334&r=1",
        "2020" = "https://nschdata.org/browse/survey/results?q=9135&r=1"
    )
    urls_race <- list(
        "2016" = "https://nschdata.org/browse/survey/results?q=4786&r=1&g=606",
        "2017" = "https://nschdata.org/browse/survey/results?q=6763&r=1&g=675",
        "2018" = "https://nschdata.org/browse/survey/results?q=7445&r=1&g=757",
        "2019" = "https://nschdata.org/browse/survey/results?q=8334&r=1&g=832",
        "2020" = "https://nschdata.org/browse/survey/results?q=9135&r=1&g=935"
        
    )
    # Xpath associated to url link
    xpath_table <-
        '//*[@id="PageContent_PageContent_C001_tblResults"]'
    
    df.nsch <- tibble::tibble()
    for (y in as.character(2016:2020)) {
        # National level
        scraped_table <- urls_nat[[y]] |>
            rvest::read_html() |>
            rvest::html_nodes(xpath = xpath_table) |>
            rvest::html_table()
        scraped_table <- scraped_table[[1]]
        p_nat <- scraped_table[1, 2]
        lower95_nat <- substr(scraped_table[2, 2], 1, 3)
        upper95_nat <- substr(scraped_table[2, 2], 7, 9)
        
        # Race/ethnic level
        scraped_table_race <- urls_race[[y]] |>
            rvest::read_html() |>
            rvest::html_nodes(xpath = xpath_table) |>
            rvest::html_table()
        scraped_table_race <- scraped_table_race[[1]]
        names(scraped_table_race)[1:3] <-
            c("race", "metric", "Qtrue")
        p_race <-
            scraped_table_race |>
            dplyr::filter(race %in% c(
                "Hispanic",
                paste0("White, ", c("n", "N"), "on-Hispanic"),
                paste0("Black, ", c("n", "N"), "on-Hispanic")
            ), metric == "%") |>
            dplyr::pull(Qtrue)
        CI <-
            scraped_table_race |>
            dplyr::filter(race %in% c(
                "Hispanic",
                paste0("White, ", c("n", "N"), "on-Hispanic"),
                paste0("Black, ", c("n", "N"), "on-Hispanic")
            ),
            metric == "C.I.") |>
            dplyr::pull(Qtrue)
        lower95_race <- substr(CI, 1, 3)
        upper95_race <- substr(CI, 7, 9)
        
        
        df.out <- tibble::tibble(
            year = as.numeric(y),
            median = c(p_nat, p_race),
            # not sure it is the median but makes compatible with data of estimates
            lower95 = c(lower95_nat, lower95_race),
            upper95 = c(upper95_nat, upper95_race),
            race = c("Total", "Hispanic", "NH White", "NH Black"),
            source = "NSCH"
        ) |>
            dplyr::mutate(dplyr::across(median:upper95, ~ as.numeric(.x)))
        df.nsch <- rbind(df.nsch, df.out)
        
    }
    
    # Save
    saveRDS(
        df.nsch,
        here::here("data", "nsch_estimates_allcause_parental_death.RDS"),
        compress = "xz"
    )
}
