library(tidyverse)
library(here)
source(here("code", "mk_nytimes.R"))

death_df <- readRDS(here("data", "summarized_deaths_1990_2020.RDS"))
fx_df <- readRDS(here("data", "fx_US.rds"))

sub_fx <- fx_df |>
    filter(race_eth == "total", year >= 1999) |>
    rename(age_years = age) |>
    group_by(age_years, sex) |>
    summarize(fx = mean(fx)) |>
    group_by(sex) |> 
    mutate(cume_total_fx = cumsum(fx) / sum(fx)) |> 
    mutate(remaining_fx = 1 - cume_total_fx) |>
    ungroup() |>
    arrange(sex, age_years)

sub_df <- death_df |>
    filter(year >= 1999) |>
    group_by(age_years, sex) |>
    summarize(n_cancer = sum(n_cancer),
              n_deaths = sum(n_deaths)) |>
    mutate(n_noncancer = n_deaths - n_cancer) |>
    pivot_longer(starts_with("n_"), names_to = "death_type", values_to = "n") |>
    mutate(death_type = gsub("n_", "", death_type)) |>
    arrange(death_type, sex, age_years) |>
    group_by(death_type, sex) |>
    mutate(prop = n / sum(n)) |>
    mutate(cume_prop = cumsum(prop)) |>
    ungroup()

sub_df <- sub_df |>
    left_join(sub_fx) |>
    filter(!is.na(fx)) |>
    mutate(
        sex_cat = factor(
            sex,
            levels = c("female", "male"),
            labels = c("Female", "Male"),
            ordered = TRUE
        ),
        death_cat = factor(
            death_type,
            levels = c("cancer", "noncancer", "deaths"),
            labels = c("Cancer", "Non-cancer", "All deaths"),
            ordered = TRUE
        )
    )

p1 <- ggplot() +
    geom_col(
        data = sub_df |>
            filter(death_type == "cancer"),
        aes(x = age_years, y = fx),
        fill = "black",
        color = NA,
        alpha = .25
    ) +
    geom_line(data = sub_df, aes(x = age_years, y = cume_prop, color = death_cat)) +
    
    facet_grid( ~ sex_cat) +
    coord_cartesian(xlim = c(10, 50), ylim = c(0, .2)) +
    mk_nytimes(legend.position = "bottom") +
    scale_color_brewer("Death type", palette = "Set1") +
    scale_x_continuous("Age (years)") +
    scale_y_continuous("Cumulative proportion of deaths (lines)\nAge-specific fertility rate (bars)")

p2 <- ggplot(data = sub_df |> filter(death_type == "deaths"),
             aes(x = age_years, y = remaining_fx)) +
    geom_line() +
    facet_wrap(~ sex_cat) +
    coord_cartesian(xlim = c(10, 50), ylim = c(0, 1)) +
    mk_nytimes(legend.position = "bottom") +
    scale_color_brewer("Death type", palette = "Set1") +
    scale_x_continuous("Age (years)", expand = c(0, .01)) +
    scale_y_continuous("Proportionate remaining fertility", expand = c(0, .01))

ggsave(
    here("figures", "figS1_age_distribution.pdf"),
    p1,
    width = 7,
    height = 4,
    scale = 1,
    device = cairo_pdf
)

ggsave(
    here("figures", "figS1_age_distribution.jpg"),
    p1,
    width = 7,
    height = 4,
    scale = 1,
    dpi = 600
)

ggsave(
    here("figures", "figS2_remaining_fertility.pdf"),
    p2,
    width = 7,
    height = 4,
    scale = 1,
    device = cairo_pdf
)

ggsave(
    here("figures", "figS2_remaining_fertility.jpg"),
    p2,
    width = 7,
    height = 4,
    scale = 1,
    dpi = 600
)

write_csv(
    sub_df |>
        filter(between(age_years, 10, 50)) |>
        select(death_cat, sex_cat, age_years, fx, n, prop, cume_prop),
    here("output", "figS1_age_distribution.csv")
)

write_csv(
    sub_df |> 
        filter(death_type == "deaths") |> 
        select(age_years, remaining_fx, fx, sex),
    here("output", "figS2_remaining_fertility.csv")
)
