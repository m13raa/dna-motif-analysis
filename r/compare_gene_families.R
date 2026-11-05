library(readr)
library(dplyr)
library(ggplot2)

motif_data <-read_csv(
    "results/motif_summary.csv",
    show_col_types = FALSE
)

motif_data <- motif_data %>%
    mutate(class = factor(class))

print(table(motif_data$class))

#Kruskal-Wallis test

kruskal_results <- motif_data %>%
    group_by(motif_name) %>%
    group_modify(~ {
        test <- kruskal.test(
            density_per_kb ~ class,
            data = .x
        )
        data.frame(
            statistic = unname(test$statistic),
            p_value = test$p.value
        )
    }) %>%
    ungroup()

print(kruskal_results)

#Correct for multiple testing

kruskal_results <- kruskal_results %>%
    mutate(
        adjusted_p = p.adjust(
            p_value,
            method = "BH"
        ), 
        significant = adjusted_p < 0.05
    )

print(kruskal_results)

write_csv(
    kruskal_results,
    "results/kruskal_wallis_results.csv"
)

#Pairwise comparisons

significant_motifs <- kruskal_results %>%
    filter(significant) %>%
    pull(motif_name)

pairwise_results <-list()

for (motif in significant_motifs) {

    subset_data <- motif_data %>%
    filter(motif_name == motif)

    test <-pairwise.wilcox.test(
        x = subset_data$density_per_kb,
        g = subset_data$class,
        p.adjust.method = "BH",
        exact = FALSE
    )

    print(paste("Pairwise comparisons for:", motif))
    print(test)

    pairwise_results[[motif]] <- test$p.value
}

saveRDS(
    pairwise_results,
    "results/pairwise_gene_family_tests.rds"
)