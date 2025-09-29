library("ggplot2")
library("dplyr")
library("readxl")
library("dplyr")
library("tidyr")
library("readr")
library("XML")
library("xml2")
library("tidyverse")
library("stringdist")
library("hrbrthemes")
library("viridis")
library("plotrix")
library("stringr")
library("RColorBrewer")
library("rcompanion")

uni_choice <- "Q9NYQ8"

hansson_data_fragment <- read_excel("hansson.xlsx", sheet = 3)
hansson_data_fragment <- hansson_data_fragment %>% rename("Fragment Sequence" = "Sequence")
hansson_data_fragment <- hansson_data_fragment %>% rename("Fragment Sequence Length" = "Sequence Length")
hansson_data_fragment <- hansson_data_fragment %>% mutate("Positions in Master Proteins" = gsub("\\]; \\[", "], [", hansson_data_fragment$"Positions in Master Proteins"))
hansson_data_fragment <- hansson_data_fragment %>% separate_rows("Positions in Master Proteins", sep = "; ")
hansson_data_fragment <- hansson_data_fragment %>% separate("Positions in Master Proteins", into = c("Uniprot", "Location"), sep = "\\ ")
hansson_data_fragment <- hansson_data_fragment %>% filter("Location" != "" & !is.na("Location"))
hansson_data_fragment$"Uniprot" <- sub("\\-.*", "", hansson_data_fragment$"Uniprot")
hansson_data_fragment <- hansson_data_fragment %>% distinct(hansson_data_fragment$"Fragment Sequence", hansson_data_fragment$"Uniprot", .keep_all = TRUE)

uniprot_data <- read_tsv("uniprot_data.tsv")
colnames(uniprot_data) <- c("Uniprot", "Full Length", "Full Sequence", "Gene")
ptm_data <- read_tsv("ptm.tsv")
colnames(ptm_data)[1] <- "Uniprot"
uniprot_data <- inner_join(uniprot_data, ptm_data, by = "Uniprot")
uniprot_data$"Gene" <- sub("\\ . *", "", uniprot_data$"Gene")

antibodies <- read.csv("HPA_antibodies.csv", sep=";")
antibodies <- antibodies %>% separate_rows("Uniprot", sep = ";")

antigen <- read.csv("antigen.csv", sep=";")
antigen <- antigen %>% separate_rows("Uniprot", sep = ";")
antigen <- antigen %>% rename("Antigen Sequence" = "PrEST.seq..aa.")
antigen <- antigen %>% select("Antigen Sequence", "Uniprot", "Single.cell.type.Enriched.Enhanced.Tissue")


diff <- data.frame(setdiff(hansson_data_fragment$"Uniprot", uniprot_data$"Uniprot"))

hpa <- read_tsv("proteinatlas.tsv")
hpa <- hpa %>% select("Gene", "Uniprot", "Protein class")

complete_hansson <- inner_join(uniprot_data, hansson_data_fragment, by = "Uniprot")
complete_hansson <- complete_hansson %>% select(-"Gene")
complete_hansson <- inner_join(hpa,complete_hansson, by = "Uniprot")

complete_hansson <- complete_hansson %>% select("Gene", "Uniprot", "Full Length", "Full Sequence", "Fragment Sequence", "Fragment Sequence Length", "Location", "Protein class")

hansson_counts <- complete_hansson %>% count(Gene, Uniprot, name = "Hansson Count")

brain_enriched <- read.csv("brain_enriched.tsv", sep = "\t")
brain_enriched_desc_uni <- select(brain_enriched, "Uniprot", "RNA.tissue.specific.nTPM", "RNA.single.cell.type.specific.nTPM", "Protein.class")
brain_enriched_desc_uni <- brain_enriched_desc_uni %>% rename("RNA Tissue nTPM" = "RNA.tissue.specific.nTPM")
brain_enriched_desc_uni <- brain_enriched_desc_uni %>% rename("RNA Single Cell nTPM" = "RNA.single.cell.type.specific.nTPM")

brain_enriched_desc_uni$"RNA Single Cell nTPM" <- sapply(brain_enriched_desc_uni$"RNA Single Cell nTPM", function(row) {
  if (row == "" | row == ":") {
    return("")  # Om tom, returnera tom sträng
  }
  pairs <- unlist(strsplit(row, ";"))
  split_pairs <- strsplit(pairs, ": ")
  names <- sapply(split_pairs, `[`, 1)
  values <- as.numeric(sapply(split_pairs, `[`, 2))
  if (all(is.na(values))) {
    return("")
  }
  sorted_indices <- order(values, decreasing = TRUE)
  paste0(names[sorted_indices], ": ", values[sorted_indices], collapse = ";")
})

brain_enriched_desc_uni <- brain_enriched_desc_uni %>%
  mutate(
    # Om cell_types är tomt, sätt både first_cell_type och other_cell_types till tomma strängar
    first_cell_type = ifelse(brain_enriched_desc_uni$"RNA Single Cell nTPM" == "", "", sub("(^[^;]+);.*", "\\1", brain_enriched_desc_uni$"RNA Single Cell nTPM")),  # Extrahera första celltypen och dess uttryck
    other_cell_types = ifelse(brain_enriched_desc_uni$"RNA Single Cell nTPM" == "", "", sub("^[^;]+;(.*)", "\\1", brain_enriched_desc_uni$"RNA Single Cell nTPM"))  # Extrahera de andra celltyperna
  ) %>%
  # Om det inte finns några andra celltyper (bara en celltyp), sätt other_cell_types till en tom sträng
  mutate(
    other_cell_types = ifelse(first_cell_type != "" & !grepl(";", brain_enriched_desc_uni$"RNA Single Cell nTPM"), "", other_cell_types)
  )

common_proteins <- inner_join(hansson_counts, brain_enriched_desc_uni, by = "Uniprot")
common_proteins_sorted <- common_proteins[order(-common_proteins$"Hansson Count"), ]

protein_sep <- complete_hansson %>% filter(grepl(uni_choice, complete_hansson$"Uniprot"))
antigen_sep <- antigen %>% filter(grepl(uni_choice, antigen$"Uniprot"))

cp20 <- common_proteins_sorted[c(1:50),]

hansson_full <- hansson_counts[order(-hansson_counts$"Hansson Count"), ]
hansson_counts <- hansson_counts[order(-hansson_counts$"Hansson Count"), ]
hansson_full_50 <- hansson_counts[c(1:50),]

hansson_unique <- complete_hansson %>% distinct(Uniprot)
common_unique <- common_proteins_sorted %>% distinct(Uniprot)
hansson_normal <- setdiff(hansson_unique,common_unique)

hansson_normal_count <- nrow(hansson_normal)
common_unique_count <- nrow(common_unique)


hansson_brain_enriched <- complete_hansson %>% filter(Uniprot %in% common_proteins$Uniprot)
hansson_non_brain_enriched <- complete_hansson %>% filter(!Uniprot %in% common_proteins$Uniprot)


non_brain_enriched_fragments <- hansson_non_brain_enriched %>% select("Gene","Uniprot", "Fragment Sequence Length")
brain_enriched_fragments <- hansson_brain_enriched %>% select("Gene","Uniprot", "Fragment Sequence Length")

non_brain_enriched_fragment_length <- data_frame(name="Non-Brain-Enriched", values = non_brain_enriched_fragments$"Fragment Sequence Length")
brain_enriched_fragment_length <- data_frame(name="Brain-Enriched", values = brain_enriched_fragments$"Fragment Sequence Length")

fragment_length_by_enrichment <- rbind(non_brain_enriched_fragment_length, brain_enriched_fragment_length)

non_brain_enriched_fragment_count <- non_brain_enriched_fragments %>% count(Gene, Uniprot, name = "Hansson Count")
brain_enriched_fragment_count <- brain_enriched_fragments %>% count(Gene, Uniprot, name = "Hansson Count")

non_brain_enriched_fragment_count <- data_frame(name="Non-Brain-Enriched", values = non_brain_enriched_fragment_count$"Hansson Count")
brain_enriched_fragment_count <- data_frame(name="Brain-Enriched", values = brain_enriched_fragment_count$"Hansson Count")

fragment_count_by_enrichment <- rbind(non_brain_enriched_fragment_count, brain_enriched_fragment_count)

a <- read_excel("hansson.xlsx", sheet = 3)
a <- a %>% mutate("Ma" = gsub("\\]; \\[", "], [", a$"Master Protein Accessions"))
a <- a %>% separate_rows("Master Protein Accessions", sep = "; ")
a <- a %>% distinct(a$"Master Protein Accessions", .keep_all = TRUE)

print(mean(hansson_counts$`Hansson Count`))

