df <- uniprot_data
df <- df %>% select(-`Gene`,-`Full Sequence`, -`Post-translational modification`)
gene_name <- hpa %>% select(Gene, Uniprot)
df <- inner_join(gene_name, df, by = "Uniprot")
prot <- df %>% filter(Uniprot == uni_choice)
prot <- prot %>%
  mutate(`Initiator methionine` = ifelse(is.na(`Initiator methionine`), NA, "Removed"))

extract_site_and_note_filtered <- function(vec, keyword = "SITE") {
  sapply(vec, function(x) {
    if (is.na(x)) return(NA)
    x <- str_remove_all(x, "/evidence=\"[^\"]*\";?")
    str_split(x, paste0("(?=", keyword, " )"))[[1]] |>
      lapply(function(block) {
        if (!str_starts(str_trim(block), keyword)) return(NULL)
        site <- str_extract(block, "\\d+(\\.\\.\\d+)?")
        note <- str_extract(block, '/note="[^"]+"') |>
          str_remove('^/note=') |>
          str_remove_all('"') |>
          str_replace_all(";", ":") |>
          str_remove_all(" ?\\(.*?\\)") |>
          str_trim()
        if (!is.na(note)) paste(site, note) else site
      }) |>
      unlist() |>
      {\(res) if (length(res) == 0) NA else paste(res, collapse = "; ")}()
  })
}

prot$`Disulfide bond` <- extract_site_and_note_filtered(prot$`Disulfide bond`, "DISULFID")
prot$`Modified residue` <- extract_site_and_note_filtered(prot$`Modified residue`, "MOD_RES")
prot$`Glycosylation` <- extract_site_and_note_filtered(prot$`Glycosylation`, "CARBOHYD")
prot$`Lipidation` <- extract_site_and_note_filtered(prot$`Lipidation`, "LIPID")
prot$`Transit peptide` <- extract_site_and_note_filtered(prot$`Transit peptide`, "TRANSIT")
prot$`Peptide` <- extract_site_and_note_filtered(prot$`Peptide`, "PEPTIDE")
prot$`Propeptide` <- extract_site_and_note_filtered(prot$`Propeptide`, "PROPEP")
prot$`Signal peptide` <- extract_site_and_note_filtered(prot$`Signal peptide`, "SIGNAL")
prot$`Cross-link` <- extract_site_and_note_filtered(prot$`Cross-link`, "CROSSLNK")
prot$`Site` <- extract_site_and_note_filtered(prot$`Site`, "SITE")
prot$`Chain` <- extract_site_and_note_filtered(prot$`Chain`, "CHAIN")
