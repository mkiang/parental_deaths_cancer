## Imports ----
library(tidyverse)
library(here)
library(ggsci)
source(here::here("code", "utils.R"))

## Data ----
results_df <- 
    bind_rows(
        readRDS(here("data", "summarized_results.RDS")) |>
            mutate(scenario = 1), 
        readRDS(here("data", "summarized_results_lowerfx.RDS"))
    ) |> 
    filter(cause == "cancer", sex == "both", race == "total", metric == "number") |>
    categorize_race()

results_df <- results_df |>
    group_by(race, cause, scenario) |>
    mutate(
        cume_median = cumsum(median),
        cume_lower95 = cumsum(lower95),
        cume_upper95 = cumsum(upper95)
    ) |>
    ungroup() |>
    mutate(scenario_cat = factor(
        scenario,
        levels = rev(seq(.8, 1, .05)),
        labels = c(
            "Baseline model",
            "5% lower fertility",
            "10% lower fertility",
            "15% lower fertility",
            "20% lower ferility"
        ),
        ordered = TRUE
    ))

p1 <- ggplot(
    results_df |> filter(metric == "number"),
    aes(
        x = year,
        y = cume_median,
        ymin = cume_lower95,
        ymax = cume_upper95,
        color = scenario_cat,
        fill = scenario_cat,
        group = scenario_cat
    )
) +
    geom_ribbon(color = NA, alpha = .1) +
    geom_line() + 
    scale_color_jama(name = "Fertility scenario") + 
    scale_fill_jama(name = "Fertility scenario") + 
    scale_x_continuous(NULL, expand = c(0, 0)) + 
    scale_y_continuous("Cumulative number of youth (in millions)",
                       expand = c(0, 0),
                       labels = function(x) round(x / 1000000, 1)) + 
    theme_bw()

## Save ----
ggplot2::ggsave(
    here::here("figures", "figS3_lower_fertility.jpg"),
    p1,
    width = 7,
    height = 4,
    scale = 1,
    dpi = 600,
    units = "in"
)
ggplot2::ggsave(
    here::here("figures", "figS3_lower_fertility.pdf"),
    p1,
    width = 7,
    height = 4,
    scale = 1,
    device = grDevices::cairo_pdf,
    units = "in"
)

readr::write_csv(
    results_df |>
        filter(metric == "number") |> 
        select(race_cat, cause, year, starts_with("cume_")),
    here::here("output", "figS3_lower_fertility.csv")
)
