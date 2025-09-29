frag_bar <- ggplot(cp20, aes(x = factor(cp20$"Gene", levels = cp20$"Gene"), y = cp20$"Hansson Count", fill = cp20$"Hansson Count")) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = paste0("n=", cp20$"Hansson Count")), hjust = -0.1, size = 3, angle = 90) +
  scale_fill_gradient(low="skyblue", high ="darkblue", guide = "none") +
  theme_minimal() +
  labs(
    title ="Fragment Distribution of Brain-Enriched Proteins", 
    x = "Protein",
    y = "Number of Fragments",
  ) +
  coord_cartesian(ylim = c(0, max(cp20$"Hansson Count")*1.05)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
frag_bar

hansson_bar <- ggplot(hansson_full_50, aes(x = factor(hansson_full_50$"Gene", levels = hansson_full_50$"Gene"), y = hansson_full_50$"Hansson Count", fill = hansson_full_50$"Hansson Count")) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = paste0("n=", hansson_full_50$"Hansson Count")), hjust = -0.1, size = 3, angle = 90) +
  scale_fill_gradient(low="orange", high ="brown", guide = "none") +
  theme_minimal() +
  labs(
    title ="Fragment Distribution of All Proteins", 
    x = "Protein",
    y = "Number of Fragments",
  ) +
  coord_cartesian(ylim = c(0, max(hansson_full_50$"Hansson Count")*1.05)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
hansson_bar

