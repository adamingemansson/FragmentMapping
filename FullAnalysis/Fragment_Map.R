stor_sekvens <- protein_sep$"Full Sequence"[1]
mindre_sekvenser <- antigen_sep$"Antigen Sequence"

match_resultat <- lapply(mindre_sekvenser, function(sekvens) {
  sekvens_längd <- nchar(sekvens)
  bästa_distans <- Inf
  bästa_position <- NA
  
  # Loopa genom stor_sekvens och hitta delsekvens med minst avstånd
  for (i in 1:(nchar(stor_sekvens) - sekvens_längd + 1)) {
    delsekvens <- substr(stor_sekvens, i, i + sekvens_längd - 1)
    distans <- stringdist(delsekvens, sekvens, method = "lv")
    
    if (distans < bästa_distans) {
      bästa_distans <- distans
      bästa_position <- i
    }
  }
  
  # Tillåt upp till 90% fel (dvs krävs minst 10% rätta tecken)
  max_fel <- floor(sekvens_längd * 0.9)
  
  if (!is.na(bästa_position) && bästa_distans <= max_fel) {
    data.frame(Sekvens = sekvens,
               Start_Position = bästa_position,
               End_Position = bästa_position + sekvens_längd - 1,
               Perfekt_Match = (bästa_distans == 0))
  } else {
    data.frame(Sekvens = sekvens,
               Start_Position = NA,
               End_Position = NA,
               Perfekt_Match = FALSE)
  }
})

match_resultat_df <- do.call(rbind, match_resultat)

print(match_resultat_df)

data <- protein_sep

data <- data %>% mutate(Start = as.numeric(gsub('\\[|\\]|-.*', '', Location)),
                        End = as.numeric(gsub('\\[.*-|\\]', '', Location)),
                        Fragment_Length = End - Start + 1)
data <- data %>% arrange(Start)
#data <- data %>% arrange(desc(Fragment_Length), Start)

data$Row <- seq_len(nrow(data))

full_length <- data$"Full Length"[1]

overlap_counts <- data.frame(Position = 1:full_length, Overlap = 0)
for (i in 1:nrow(data)) {
  overlap_counts$Overlap[data$Start[i]:data$End[i]] <- overlap_counts$Overlap[data$Start[i]:data$End[i]] + 1
}

antal_fragment <- nrow(data)
linje_tjocklek <- antal_fragment * 0.04
offset_y <- antal_fragment * 0.025

gen_namn <- protein_sep$Gene[1]

ggplot() +
  geom_tile(data = overlap_counts, aes(x = Position, y = -offset_y, fill = Overlap), height = linje_tjocklek) +
  geom_segment(aes(x = 1, xend = full_length, y = -offset_y, yend = -offset_y), color = 'black', size = 3) +
  geom_segment(data = data, aes(x = Start, xend = End, y = Row, yend = Row, color = Fragment_Length), size = 0.5) +
  geom_segment(data = data, aes(x = Start, xend = End, y = Row, yend = Row), color = 'black', linewidth = 0.08) +
  geom_segment(data = match_resultat_df, aes(x = Start_Position, xend = End_Position, y = -offset_y, yend = -offset_y), color = '#008000', size = 2) +
  geom_text(aes(x = 1, y = -offset_y, label = "1"), color = "black", size =4, hjust = 1.7) +
  geom_text(aes(x = full_length, y = -offset_y, label = full_length), color = "black", size = 4, hjust = -0.3) +
  geom_text(data = match_resultat_df %>% distinct(Sekvens, .keep_all = TRUE), 
            aes(x = Start_Position, y = offset_y-antal_fragment*0.05, label = Start_Position), 
            color = "black", size = 4, vjust = 2, hjust = 0.5) +
  geom_text(data = match_resultat_df %>% distinct(Sekvens, .keep_all = TRUE), 
            aes(x = End_Position, y = offset_y-antal_fragment*0.05, label = End_Position), 
            color = "black", size = 4, vjust = -1, hjust = 0.5) +
  scale_fill_gradient(low = 'white', high = 'red', name = "Fragment Overlap") +
  scale_color_gradient(low = "lightblue", high = "darkblue", name = "Fragment Length") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = paste0(' Fragment Mapping for ', gen_namn, ' (Uniprot ID: ', uni_choice, ')\n Total Fragments: ', antal_fragment))




