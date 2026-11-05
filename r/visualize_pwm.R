library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

#Load the Python-generated PWM

pwm<- read_csv(
    "results/pwm_log_odds.csv",
    show_col_types = FALSE
)

#Rename the first column to the nucleotide
names(pwm)[1] <- "nucleotide"

print(pwm)

#Convert the matrix into long format
pwm_long <- pwm %>%
    pivot_longer(
        cols = -nucleotide,
        names_to = "position",
        values_to = "score"
    ) %>%
    mutate(
        position = as.integer(position)
    )

print(head(pwm_long))

#Creat the PWM heatmap

pwm_plot <-ggplot(
    pwm_long,
    aes(
        x = position,
        y = nucleotide,
        fill = score
    )
) +
    geom_text(
        aes(label = round(score,2)),
        size = 4
    ) +
    scale_fill_gradient2(
        low = "steelblue",
        mid = "white",
        high = "firebrick",
        midpoint = 0
    ) +
    scale_x_continuous(
        breaks = sort(unique(pwm_long$position))
    ) +
    labs(
        title = "DNA Position Weight Matrix",
        subtitle = "Nucleotide log-odds scores",
        x = "Motif Position",
        y = "Nucleotide",
        fill = "Log2 Odds"
    ) +
    theme_minimal()

print(pwm_plot)

ggsave(
    "figures/pwm_heatmap.png",
    plot = pwm_plot,
    width = 8,
    height = 5,
    dpi = 300
)