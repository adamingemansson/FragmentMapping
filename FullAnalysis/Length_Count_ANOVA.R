group_means <- fragment_length_by_enrichment %>%
  group_by(name) %>%
  summarise(mean_value = mean(values))

sample_sizes <- fragment_length_by_enrichment %>%
  group_by(name) %>%
  summarise(n = n())

whisker_max <- fragment_length_by_enrichment %>%
  group_by(name) %>%
  summarise(
    Q1 = quantile(`values`, 0.25),
    Q3 = quantile(`values`, 0.75),
    IQR = IQR(`values`),
    upper_whisker = max(`values`[`values` <= (Q3 + 1.5 * IQR)])
  )

y_limit <- max(whisker_max$upper_whisker)

# Skapa ggplot
ggplot(fragment_length_by_enrichment, aes(x = name, y = values, fill = name)) +
  geom_boxplot(width = 0.6, staplewidth = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = c("#6e1b1b", "#a65e5e")) +
  geom_point(data = group_means,
             aes(x = name, y = mean_value),
             color = "black", size = 2) +
  geom_text(data = group_means,
            aes(x = name, y = mean_value, label = paste0("Mean: ", round(mean_value, 2))),
            color = "black", vjust = -0.8, size = 3) +
  geom_text(data = sample_sizes,
            aes(x = name, y = 0, label = paste0("n = ", n)),
            vjust = 1.5, size = 3, color = "black") +
  coord_cartesian(ylim = c(0, y_limit * 1.025))+
  labs(
    title = "Fragment Length by Tissue Enrichment",
    y = NULL,
    x = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

wilcox.test(values ~ name,
            data = fragment_length_by_enrichment,
            alternative = "greater")

group_means <- fragment_count_by_enrichment %>%
  group_by(name) %>%
  summarise(mean_value = mean(values))

sample_sizes <- fragment_count_by_enrichment %>%
  group_by(name) %>%
  summarise(n = n())

whisker_max <- fragment_count_by_enrichment %>%
  group_by(name) %>%
  summarise(
    Q1 = quantile(`values`, 0.25),
    Q3 = quantile(`values`, 0.75),
    IQR = IQR(`values`),
    upper_whisker = max(`values`[`values` <= (Q3 + 1.5 * IQR)])
  )

y_limit <- max(whisker_max$upper_whisker)

# Skapa ggplot
ggplot(fragment_count_by_enrichment, aes(x = name, y = values, fill = name)) +
  geom_boxplot(width = 0.6, staplewidth = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = c("#6e1b1b", "#a65e5e")) +
  geom_point(data = group_means,
             aes(x = name, y = mean_value),
             color = "black", size = 2) +
  geom_text(data = group_means,
            aes(x = name, y = mean_value, label = paste0("Mean: ", round(mean_value, 2))),
            color = "black", vjust = 0.5, hjust = -0.2, size = 3) +
  geom_text(data = sample_sizes,
            aes(x = name, y = 0, label = paste0("n = ", n)),
            vjust = 1.5, size = 3, color = "black") +
  coord_cartesian(ylim = c(0, y_limit * 1.025))+
  labs(
    title = "Number of Fragments by Tissue Enrichment",
    y = NULL,
    x = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

wilcox.test(values ~ name,
            data = fragment_count_by_enrichment,
            alternative = "greater")

