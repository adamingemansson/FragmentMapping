# Histogram över fragmentens längd för alla proteiner i complete_hansson
ggplot(hansson_non_brain_enriched, aes(x = `Fragment Sequence Length`)) +
  geom_histogram(binwidth = 5, fill = "orange", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Hansson Full",
       x = "Fragment Sequence Length",
       y = "Fragment Count")

# Histogram över fragmentens längd för proteiner i common_proteins
ggplot(hansson_brain_enriched, aes(x = `Fragment Sequence Length`)) +
  geom_histogram(binwidth = 5, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Brain-Enriched",
       x = "Fragment Sequence Length",
       y = "Fragment Count")