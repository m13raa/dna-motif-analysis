library (readr)
library (dplyr)
library (ggplot2)

#Load the results of the motif analysis
motif_data <- read_csv(
    "results/motif_summary.csv",
    show_col_types = FALSE
)

#Inspect dataset
print(head(motif_data))
print(names(motif_data))

#Create a directory for the output

dir.create(
    "figures",
    showWarnings = FALSE,
    recursive = TRUE
)

#Summarize the motif frequencies

motif_totals <-motif_data %>%
    group_by(motif_name) %>%
    summarize(
        total_count = sum(count),
        mean_count = mean(count),
        .groups = "drop"
   ) %>%
   arrange(desc(total_count))
print(motif_totals)

#Plot the total motif ocunts
p1 <- ggplot(
    motif_totals,
    aes(
        x = reorder(motif_name, total_count),
        y = total_count,
    fill = motif_name
    )
) +
    geom_col() +
    coord_flip() + 
    labs(
        title = "Total DNA Motif Frequencies",
        x = "DNA Motif",
        y = "Total Occurences"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

print(p1)

ggsave(
    "figures/motif_total_frequencies.png",
    plot = p1,
    width = 9,
    height = 6,
    dpi = 300
)

#Compare the motif densities

p2 <- ggplot(
    motif_data,
    aes(
        x = motif_name,
        y = density_per_kb,
        fill = motif_name
    )
) +
    geom_boxplot(outlier.alpha = 0.3) +
    coord_flip() +
    labs(
        title = "DNA Motif Density Distribution",
        x = "DNA Motif",
        y = "Occurrences per kb"
    ) +
    theme_minimal() + 
    theme(legend.position = "none")
print(p2)

ggsave(
    "figures/motif_density_boxplot.png",
    plot = p2,
    width = 9,
    height = 6, 
    dpi = 300
)

motif_positions <-read_csv(
    "results/motif_positions.csv",
    show_col_types = FALSE
)

print(names(motif_positions))
print(head(motif_positions))

p3 <- ggplot(
    motif_positions,
    aes(
        x = `relative position`)
) +
    geom_histogram(
        bins = 30,
        fill = "steelblue",
        color = "white"
    ) +
    facet_wrap(
        ~ motif_name,
        scales = "free_y"
    ) +
    labs(
        title = "Distribution of DNA Motif Positions",
        x = "Relative Position Along Sequence",
        y = "Number of Occurrences"
    ) +
    theme_minimal() 

print(p3)

ggsave(
    "figures/motif_position_distribution.png",
    plot = p3,
    width = 11,
    height = 7,
    dpi = 300
)

#Join the sequence lengths to the motif positional data

position_data <- motif_positions %>%
    left_join(
        motif_data %>%
            distinct(
                sequence_id,
                sequence_length
            ),
        by = "sequence_id"
    ) %>%
    mutate(
        relative_position = start_0based / sequence_length
    )

#Plot the normalized motif positions
p4 <- ggplot(
    position_data,
    aes(
        x = relative_position,
        fill = motif_name
    )
) +
    geom_histogram(
        bins = 30,
        color = "white"
    ) +
    facet_wrap(
        ~ motif_name,
        scales = "free_y"
    ) +
    labs(
        title = "Relative DNA Motif Positions",
        x = "Relative Position Along Sequence",
        y = "Number of Occurrences"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

print(p4)

ggsave(
    "figures/motif_relative_positions.png",
    plot = p4,
    width = 11,
    height = 7,
    dpi = 300
)

cat("\nAll four motif visualizations are saved successfully!\n")

#Add scatter plot for GC content relationships

p5 <- ggplot(
    motif_data,
    aes(
        x = gc_content,
        y = density_per_kb,
        color = motif_name
    )
) +
    geom_point(
        alpha = 0.25,
        size = 1
    ) +
    geom_smooth(
        method = "lm",
        se = TRUE
    ) +
    facet_wrap(
        ~ motif_name,
        scales = "free_y"
    ) +
    labs(
        title = "GC Content vs DNA Motif Density",
        x = "GC Conent (%)",
        y = "Motif Occurrences per kb"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

print(p5)

ggsave(
    "figures/gc_concent_vs_motif_density.png",
    plot = p5,
    width = 12,
    height = 7,
    dpi = 300
)

#Calculate the correlations 

gc_correlations <- motif_data %>%
    group_by(motif_name) %>%
    summarize(
        correlation = cor(
            gc_content,
            density_per_kb,
            method = "spearman",
            use = "complete.obs"
        ),
        .groups = "drop"
    )
print(gc_correlations)

write_csv(
    gc_correlations,
    "results/gc_motif_correlations.csv"
)

#Calculate the average motif density per class

class_summary <- motif_data %>%
    group_by(class, motif_name) %>%
    summarize(
        mean_density = mean(density_per_kb),
        median_density = median(density_per_kb),
        mean_count = mean(count),
        n_sequences = n(),
        .groups = "drop"
    )

print(class_summary)

write_csv(
    class_summary,
    "results/motif_class_summary.csv"
)

#Create a heatmap

p6 <- ggplot(
    class_summary,
    aes(
        x = motif_name,
        y = factor(class),
        fill = mean_density
    )
) +
    geom_tile(color = "white") +
    scale_fill_gradient(
        low = "white",
        high = "steelblue"
    ) +
    labs(
        title = "Mean Motif Density Across Gene Families",
        x = "DNA Motif",
        y = "Gene-Family Class",
        fill = "Motifs/kb"
    ) +
    theme_minimal()

print(p6)

ggsave(
    "figures/motif_gene_family_heatmap.png",
    plot = p6,
    width = 10,
    height = 6,
    dpi = 300
)

#Grouped boxplots

p7 <- ggplot(
    motif_data,
    aes(
        x = factor(class),
        y = density_per_kb,
        fill = factor(class)
    )
) + 
    geom_boxplot(outlier.alpha = 0.2) +
    facet_wrap(
        ~ motif_name,
        scales = "free_y"
    ) +
    labs(
        title = "Motif Density by Gene-Family Class",
        x = "Gene- Gamily Class",
        y = "Motifs per kb"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

print(p7)

ggsave(
    "figures/motif_density_by_gene_family.png",
    plot =p7,
    width = 12,
    height = 7,
    dpi = 300
)
