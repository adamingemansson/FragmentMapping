# hansson_counts <- complete_hansson %>% count(Gene, Uniprot, name = "Hansson Count")
exp_df <- complete_hansson %>% select("Gene", "Uniprot", "Fragment Sequence Length", "Protein class")
exp_df_choice <- exp_df[grepl("intracellular|membrane|secreted", exp_df$"Protein class", ignore.case = TRUE), ]

extract_single_class <- function(x) {
  target_classes <- c("Predicted intracellular proteins", "Predicted membrane proteins", "Predicted secreted proteins")
  matches <- grep(paste(target_classes, collapse = "|"), x, value = TRUE)
  found <- target_classes[sapply(target_classes, function(t) grepl(t, x, fixed = TRUE))]
  if (length(found) == 1) {
    return(found)
  } else {
    return(NA)
  }
}

# Skapa ny kolumn med bara den träffade klassen
exp_df_choice$Filtered_class <- sapply(exp_df_choice$`Protein class`, extract_single_class)
exp_df_choice <- exp_df_choice[!is.na(exp_df_choice$Filtered_class), ]
exp_df_choice$`Protein class` <- NULL
exp_df_choice$`Uniprot` <- NULL

exp_count <- exp_df_choice %>% count(Gene, name = "Hansson Count")
exp_count <- inner_join(exp_count, exp_df_choice, by = "Gene")
exp_count <- exp_count %>% select(-"Fragment Sequence Length")
exp_count <- distinct(exp_count)

kruskal.test(`Fragment Sequence Length` ~ Filtered_class, data = exp_df_choice)

kruskal.test(`Hansson Count` ~ Filtered_class, data = exp_count)

anova_result <- aov(`Fragment Sequence Length` ~ Filtered_class, exp_df_choice)
summary(anova_result)

anova_result <- aov(`Hansson Count` ~ Filtered_class, exp_count)
summary(anova_result)

length_mean <- aggregate(`Fragment Sequence Length` ~ `Filtered_class`, data = exp_df_choice , FUN = mean)
count_mean <- aggregate(`Hansson Count` ~ `Filtered_class`, data = exp_count , FUN = mean)

# Sample size per kategori
sample_sizes <- exp_df_choice %>%
  group_by(Filtered_class) %>%
  summarise(n = n())

# Övre whisker per kategori (max inom 1.5 * IQR)
whisker_max <- exp_df_choice %>%
  group_by(Filtered_class) %>%
  summarise(
    Q1 = quantile(`Fragment Sequence Length`, 0.25),
    Q3 = quantile(`Fragment Sequence Length`, 0.75),
    IQR = IQR(`Fragment Sequence Length`),
    upper_whisker = max(`Fragment Sequence Length`[`Fragment Sequence Length` <= (Q3 + 1.5 * IQR)])
  )

y_limit <- max(whisker_max$upper_whisker)

# Plotten
ggplot(exp_df_choice, aes(x = Filtered_class, y = `Fragment Sequence Length`, fill = `Filtered_class`)) +
  geom_boxplot(width = 0.5, staplewidth = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = c(
    "Predicted intracellular proteins" = "#FBD7BB",
    "Predicted membrane proteins" = "#daffe7",
    "Predicted secreted proteins" = "#fadaeb"
  )) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "#D37676") +
  stat_summary(fun = mean, geom = "text", aes(label = paste("Mean:", round(..y.., 2))), 
               vjust = -0.8, hjust = 0.5, color = "black", size = 3) +
  geom_text(data = sample_sizes, aes(x = Filtered_class, y = 0, label = paste0("n = ", n)), 
            vjust = 1.5, size = 3, color = "black") +
  coord_cartesian(ylim = c(0, y_limit * 1.02)) +
  theme_minimal() +
  labs(x = "", y = "", title = "Fragment Length by Subcellular Localisation") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )


# Beräkna sample size (n) per grupp
sample_sizes <- exp_count %>%
  group_by(Filtered_class) %>%
  summarise(n = n())

# Beräkna övre whiskers för att sätta ylim
whisker_max <- exp_count %>%
  group_by(Filtered_class) %>%
  summarise(
    Q1 = quantile(`Hansson Count`, 0.25),
    Q3 = quantile(`Hansson Count`, 0.75),
    IQR = IQR(`Hansson Count`),
    upper_whisker = max(`Hansson Count`[`Hansson Count` <= (Q3 + 1.5 * IQR)])
  )

y_limit <- max(whisker_max$upper_whisker)

# Rita plotten
ggplot(exp_count, aes(x = Filtered_class, y = `Hansson Count`, fill = `Filtered_class`)) +
  geom_boxplot(width = 0.5, staplewidth = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = c(
    "Predicted intracellular proteins" = "#FBD7BB",
    "Predicted membrane proteins" = "#daffe7",
    "Predicted secreted proteins" = "#fadaeb"
  )) +
  theme_minimal() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "#D37676") +
  stat_summary(fun = mean, geom = "text", aes(label = paste("Mean:", round(..y.., 2))), 
               vjust = 0.4, hjust = -0.2, color = "black", size = 3) +
  geom_text(data = sample_sizes, aes(x = Filtered_class, y = 0, label = paste0("n = ", n)), 
            vjust = 1.5, size = 3, color = "black") +
  labs(x = "", y = "", title = "Number of Fragments by Subcellular Localisation") +
  coord_cartesian(ylim = c(0, y_limit * 1.025)) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

