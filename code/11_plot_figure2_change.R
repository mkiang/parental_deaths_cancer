## Imports ----
library(tidyverse)
library(here)
library(patchwork)
library(ggsci)
source(here::here("code", "utils.R"))

## Data ----
results_df <- readRDS(here("data", "summarized_results.RDS")) |>
    categorize_deaths() |>
    categorize_race()

## Data cleaning ----
t1 <- results_df |>
    filter(year %in% c(1999, 2020), cause == "cancer", metric == "percent") |>
    group_by(race, sex) |>
    mutate(
        rel_median = (median - lag(median)) / lag(median),
        rel_upper = (lower95 - lag(lower95)) / lag(lower95),
        rel_lower = (upper95 - lag(upper95)) / lag(upper95)
    ) |>
    filter(!is.na(rel_median)) |>
    mutate(sex_cat = factor(
        sex,
        levels = c("both", "father", "mother"),
        labels = c("Either parent", "Father", "Mother"),
        ordered = TRUE
    )) |>
    mutate(y_pos = as.numeric(race_cat_rev)) |>
    mutate(y_pos = case_when(sex == "both" ~ y_pos + 0.2, sex == "mother" ~ y_pos - 0.2, TRUE ~ y_pos)) |>
    ungroup() |>
    mutate(y_lab = sprintf(
        "%0.1f (%0.1f to %0.1f)",
        round(rel_median * 100, 1),
        round(rel_lower * 100, 1),
        round(rel_upper * 100, 1)
    )) 

p1 <- ggplot(t1,
             aes(
                 y = y_pos,
                 x = rel_median,
                 xmin = rel_lower,
                 xmax = rel_upper,
                 color = sex_cat
             )) +
    geom_vline(xintercept = 0, alpha = .5) +
    geom_errorbarh(height = 0) +
    geom_point() +
    theme_bw() +
   theme(legend.position = "bottom",
          axis.ticks.y.right = element_blank()) +
    scale_color_jama(name = NULL) +
    scale_x_continuous(NULL,
                       labels = scales::percent,
                       limits = c(-.45, .30125)) +
    scale_y_continuous(
        NULL,
        breaks = 1:4,
        labels = levels(t1$race_cat_rev)
    )

p1_num <- ggplot(t1, aes(y = y_pos, x = 1, label = y_lab)) +
    geom_text(size = 3) + 
    theme_void()

t2 <- results_df |>
    filter(
        year %in% c(1999, 2020),
        sex == "both",
        cause != "cancer",
        cause != "all",
        cause != "other_cancer",
        cause != "other",
        metric == "percent",
        race == "total"
    ) |>
    group_by(race, cause) |>
    mutate(
        rel_median = (median - lag(median)) / lag(median),
        rel_upper = (lower95 - lag(lower95)) / lag(lower95),
        rel_lower = (upper95 - lag(upper95)) / lag(upper95)
    ) |>
    filter(!is.na(rel_median)) |>
    ungroup() |>
    mutate(y_lab = sprintf(
        "%0.1f (%0.1f to %0.1f)",
        round(rel_median * 100, 1),
        round(rel_lower * 100, 1),
        round(rel_upper * 100, 1)
    )) 

p2 <- ggplot(t2,
             aes(
                 y = death_cat_alpha_rev,
                 x = rel_median,
                 xmin = rel_lower,
                 xmax = rel_upper
             )) +
    geom_vline(xintercept = 0, alpha = .5) +
    geom_errorbarh(height = 0, color = ggsci::pal_jama()(3)[1]) +
    geom_point(color = ggsci::pal_jama()(3)[1]) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.ticks.y.right = element_blank()) +
    scale_x_continuous(
        "Relative change in parental cancer deaths per 1,000 youth, %",
        labels = scales::percent,
        limits = c(-.45, .30125)
    ) +
    scale_y_discrete(NULL)

p2_num <- ggplot(t2, aes(y = death_cat_alpha_rev, x = 1, label = y_lab)) +
    geom_text(size = 3) +
    theme_void()

design <- 
"ABBB
ABBB
CDDD
CDDD
CDDD
CDDD"


p1_full <- p1_num + p1 + p2_num + p2 +
    # patchwork::plot_layout(ncol = 1, heights = c(2, 5), widths = c(1, 2))
    patchwork::plot_layout(
        design = design, 
        guides = "collect"
    ) & theme(legend.position = 'bottom')
p1_full

## Save ----
ggplot2::ggsave(
    here::here("figures", "fig2_change.jpg"),
    p1_full,
    width = 8,
    height = 7,
    scale = 1,
    dpi = 600,
    units = "in"
)
ggplot2::ggsave(
    here::here("figures", "fig2_change.pdf"),
    p1_full,
    width = 8,
    height = 7,
    scale = 1,
    device = grDevices::cairo_pdf,
    units = "in"
)

readr::write_csv(
    bind_rows(t1, t2) |> select(death_cat, sex, race_cat, starts_with("rel_")),
    here::here("output", "fig2_change.csv")
)
