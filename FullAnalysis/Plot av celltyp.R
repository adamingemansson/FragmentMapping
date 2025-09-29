cell_names <- sub(":.*", "", common_proteins_sorted$"RNA Single Cell nTPM")
cell_counts <- sort(table(cell_names), decreasing = TRUE)
labels <- paste(names(cell_counts), "(", cell_counts, ")", sep = "")

# Palett: 7 första färgade, resten ljusgrå
main_colors <- brewer.pal(7, "Set2")
extra_colors <- rep("gray90", length(cell_counts) - 7)
bar_colors <- c(main_colors, extra_colors)

# Barplot
barplot(cell_counts,
        las = 2,
        col = bar_colors,
        main = "Förekomst av celltyper",
        ylab = "Antal",
        cex.names = 0.8)

# Piechart
pie(cell_counts,
    labels = labels,
    main = "Förekomst av celltyper (antal, sorterat)",
    col = bar_colors,
    cex = 0.7,
    clockwise = TRUE)

df_filtered <- common_proteins_sorted
df_filtered$cell_type_clean <- str_extract(df_filtered$first_cell_type, "^[^:]+")
df_filtered <- df_filtered %>%
  filter(!is.na(`Hansson Count`),
         !is.na(cell_type_clean),
         `Hansson Count` <= 300) %>%
  group_by(cell_type_clean) %>%
  filter(n() > 4)
df_filtered <- df_filtered[order(df_filtered$"cell_type_clean"), ]


# Sample size
sample_sizes <- df_filtered %>%
  group_by(cell_type_clean) %>%
  summarise(n = n())

# Övre whisker
whisker_max <- df_filtered %>%
  group_by(cell_type_clean) %>%
  summarise(
    Q1 = quantile(`Hansson Count`, 0.25),
    Q3 = quantile(`Hansson Count`, 0.75),
    IQR = IQR(`Hansson Count`),
    upper_whisker = max(`Hansson Count`[`Hansson Count` <= (Q3 + 1.5 * IQR)])
  )

y_limit <- max(whisker_max$upper_whisker)

# Färgpalett
cell_types <- unique(df_filtered$cell_type_clean)
main_colors <- brewer.pal(7, "Pastel2")
fill_colors <- setNames(c(main_colors, rep("gray90", length(cell_types) - 7)), cell_types)

# Boxplot
ggplot(df_filtered, aes(x = cell_type_clean, y = `Hansson Count`, fill = cell_type_clean)) +
  geom_boxplot(width = 0.8, alpha = 1, staplewidth = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = fill_colors) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "#D37676") +
  stat_summary(fun = mean, geom = "text", aes(label = paste("Mean:", round(..y.., 2))),
               vjust = -1, hjust = 0.58, color = "black", size = 3) +
  geom_text(data = sample_sizes, aes(x = cell_type_clean, y = 0, label = paste0("n = ", n)),
            vjust = 1.5, size = 3, color = "black") +
  coord_cartesian(ylim = c(0, y_limit * 1.05)) +
  theme_minimal() +
  labs(x = "", y = "", title = "Number of Fragments by Cell Type") +
  theme(axis.text.x = element_text(angle = 15, hjust = 0.5),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
        )

kruskal.test(`Hansson Count` ~ cell_type_clean, data = df_filtered)

df_cell_length <- inner_join(complete_hansson, brain_enriched_desc_uni, by = "Uniprot")
df_cell_length$cell_type_clean <- str_extract(df_cell_length$first_cell_type, "^[^:]+")
df_cell_length <- df_cell_length %>%
  filter(!is.na(`Fragment Sequence Length`),
         !is.na(cell_type_clean),
         `Fragment Sequence Length` <= 250) %>%
  group_by(cell_type_clean) %>%
  filter(n() > 60)
df_cell_length <- df_cell_length[order(df_cell_length$"cell_type_clean"), ]

# Sample size
sample_sizes_length <- df_cell_length %>%
  group_by(cell_type_clean) %>%
  summarise(n = n())

# Övre whisker
whisker_max_length <- df_cell_length %>%
  group_by(cell_type_clean) %>%
  summarise(
    Q1 = quantile(`Fragment Sequence Length`, 0.25),
    Q3 = quantile(`Fragment Sequence Length`, 0.75),
    IQR = IQR(`Fragment Sequence Length`),
    upper_whisker = max(`Fragment Sequence Length`[`Fragment Sequence Length` <= (Q3 + 1.5 * IQR)])
  )

y_limit_length <- max(whisker_max_length$upper_whisker)

# Färgpalett
cell_types_length <- unique(df_cell_length$cell_type_clean)
main_colors <- brewer.pal(7, "Pastel2")
fill_colors_length <- setNames(c(main_colors, rep("gray90", length(cell_types_length) - 7)), cell_types_length)

# Boxplot
ggplot(df_cell_length, aes(x = cell_type_clean, y = `Fragment Sequence Length`, fill = cell_type_clean)) +
  geom_boxplot(width = 0.8, alpha = 1, staplewidth = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = fill_colors_length) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "#D37676") +
  stat_summary(fun = mean, geom = "text", aes(label = paste("Mean:", round(..y.., 2))),
               vjust = -1, hjust = 0.5, color = "black", size = 3) +
  geom_text(data = sample_sizes_length, aes(x = cell_type_clean, y = 0, label = paste0("n = ", n)),
            vjust = 1.5, size = 3, color = "black") +
  coord_cartesian(ylim = c(0, y_limit_length * 1.01)) +
  theme_minimal() +
  labs(x = "", y = "", title ="Fragment Length by Cell Type") +
  theme(axis.text.x = element_text(angle = 15, hjust = 0.5),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

kruskal.test(`Fragment Sequence Length` ~ cell_type_clean, data = df_cell_length)
