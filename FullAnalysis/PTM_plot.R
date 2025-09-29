# Filtrera protein
seq_length <- prot$`Full Length`

extract_positions <- function(raw, key = NULL) {
  if (is.na(raw) || raw == "") return(tibble(type = character(), pos = integer(), start = integer(), end = integer()))
  
  entries <- str_split(raw, ";\\s*")[[1]] %>% str_trim()
  if (!is.null(key)) entries <- entries[str_detect(entries, regex(key, ignore_case = TRUE))]
  if (length(entries) == 0) return(tibble(type = character(), pos = integer(), start = integer(), end = integer()))
  
  map_dfr(entries, ~{
    m_pair <- str_match(.x, "(?<![a-zA-Z])(\\d+)\\.\\.(\\d+)(?!\\d)")
    m_single <- str_match(.x, "(?<![a-zA-Z])(\\d+)(?!\\d)")
    
    if (!is.na(m_pair[2])) tibble(type = "pair", pos = NA_integer_, start = as.integer(m_pair[2]), end = as.integer(m_pair[3]))
    else if (!is.na(m_single[2])) tibble(type = "single", pos = as.integer(m_single[2]), start = NA_integer_, end = NA_integer_)
    else tibble(type = character(), pos = integer(), start = integer(), end = integer())
  }) %>% filter(!is.na(type) & type != "")
}

# Extrahera modifieringar
cross <- extract_positions(prot$`Cross-link`)
disulf <- extract_positions(prot$`Disulfide bond`)

gly_N <- extract_positions(prot$Glycosylation, "N-linked")
gly_O <- extract_positions(prot$Glycosylation, "O-linked")
gly_C <- extract_positions(prot$Glycosylation, "C-linked")
gly_S <- extract_positions(prot$Glycosylation, "S-linked")

lip_N <- extract_positions(prot$Lipidation, "N-")
lip_O <- extract_positions(prot$Lipidation, "O-")
lip_S <- extract_positions(prot$Lipidation, "S-")
lip_GPI <- extract_positions(prot$Lipidation, "GPI-anchor")

mod_phos  <- extract_positions(prot$`Modified residue`, "Phospho")
mod_acet  <- extract_positions(prot$`Modified residue`, "acetyl")
mod_methyl <- extract_positions(prot$`Modified residue`, "methyl")
mod_sulfo <- extract_positions(prot$`Modified residue`, "Sulfo")
mod_citru <- extract_positions(prot$`Modified residue`, "Citru")

site_gamma  <- extract_positions(prot$Site, "gamma-secretase")
site_theta  <- extract_positions(prot$Site, "theta-secretase")
site_beta   <- extract_positions(prot$Site, "beta-secretase")
site_alpha  <- extract_positions(prot$Site, "alpha-secretase")

site_trypsin    <- extract_positions(prot$Site, "trypsin")
site_cathepsin  <- extract_positions(prot$Site, "cathepsin")
site_thrombin   <- extract_positions(prot$Site, "thrombin")

site_ace    <- extract_positions(prot$Site, "ACE")
site_adam   <- extract_positions(prot$Site, "ADAM")
site_casp   <- extract_positions(prot$Site, "CASP|caspase")  # kombinerat regex
site_capn   <- extract_positions(prot$Site, "CAPN|calpain")  # kombinerat regex
site_mmp    <- extract_positions(prot$Site, "MMP")

AB40 <- extract_positions(prot$"Chain", "Amyloid-beta protein 40")
AB42 <- extract_positions(prot$"Chain", "Amyloid-beta protein 42")

all_chains <- bind_rows(
  AB42 %>% filter(type == "pair") %>% mutate(type = "Aβ42 (672-713)"),
  AB40 %>% filter(type == "pair") %>% mutate(type = "Aβ40 (672-711)")
)

all_mods <- bind_rows(
  cross %>% filter(type == "single") %>% mutate(type = "Cross-link"),
  disulf %>% filter(type == "single") %>% mutate(type = "Disulfide Bond"),
  gly_N %>% filter(type == "single") %>% mutate(type = "N-Glycosylation"),
  gly_O %>% filter(type == "single") %>% mutate(type = "O-Glycosylation"),
  gly_C %>% filter(type == "single") %>% mutate(type = "C-Glycosylation"),
  gly_S %>% filter(type == "single") %>% mutate(type = "S-Glycosylation"),
  lip_N %>% filter(type == "single") %>% mutate(type = "N-Lipidation"),
  lip_O %>% filter(type == "single") %>% mutate(type = "O-Lipidation"),
  lip_S %>% filter(type == "single") %>% mutate(type = "S-Lipidation"),
  lip_GPI %>% filter(type == "single") %>% mutate(type = "GPI-anchor"),
  mod_phos %>% filter(type == "single") %>% mutate(type = "Phosphorylation"),
  mod_acet %>% filter(type == "single") %>% mutate(type = "Acetylation"),
  mod_methyl %>% filter(type == "single") %>% mutate(type = "Methylation"),
  mod_sulfo %>% filter(type == "single") %>% mutate(type = "Sulfation"),
  mod_citru %>% filter(type == "single") %>% mutate(type = "Citrulline"),
)

cleavage <- bind_rows(
  site_gamma %>% filter(type == "pair") %>% mutate(type = "γ-Secretase Cleavage"),
  site_beta %>% filter(type == "pair") %>% mutate(type = "β-Secretase Cleavage"),
  site_theta %>% filter(type == "pair") %>% mutate(type = "θ-Secretase Cleavage"),
  site_alpha %>% filter(type == "pair") %>% mutate(type = "α-Secretase Cleavage"),
  site_trypsin %>% filter(type == "pair") %>% mutate(type = "Trypsin Cleavage"),
  site_cathepsin %>% filter(type == "pair") %>% mutate(type = "Cathepsin Cleavage"),
  site_thrombin %>% filter(type == "pair") %>% mutate(type = "Thrombin Cleavage"),
  site_ace %>% filter(type == "pair") %>% mutate(type = "ACE Cleavage"),
  site_adam %>% filter(type == "pair") %>% mutate(type = "ADAM Cleavage"),
  site_casp %>% filter(type == "pair") %>% mutate(type = "Caspase Cleavage"),
  site_capn %>% filter(type == "pair") %>% mutate(type = "Calpain Cleavage"),
  site_mmp %>% filter(type == "pair") %>% mutate(type = "MMP Cleavage"),
)

paired_mods <- bind_rows(
  cross %>% filter(type == "pair") %>% mutate(type = "Cross-link"),
  disulf %>% filter(type == "pair") %>% mutate(type = "Disulfide Bond")
)

# Plotta allt med legend
ggplot() +
  geom_rect(data = all_chains, aes(xmin = start, xmax = end, ymin = 0, ymax = 0.4, fill = type), alpha = 0.5) +
  
  # Parvis (cross + disulfid)
  geom_curve(data = paired_mods, aes(x = start, xend = end, y = 0, yend = 0, color = type),
             curvature = -1.5, linewidth = 0.8)+
  
  # Modifieringar med färglegend
  geom_segment(data = all_mods, aes(x = pos, y = 0.4, yend = 0, color = type), linewidth = 0.7) +
  
  # Cleavage
  geom_segment(data = cleavage, aes(x = start, y = 0.4, yend = 0, color = type),
               linewidth = 0.7, linetype = "dashed") +
  
  scale_color_manual(values = c(
    "Cross-link" = "#f7870f",        # stark röd
    "Disulfide Bond" = "#ed1a07",    # neutral grå
    "N-Glycosylation" = "#12e319",   # grön
    "O-Glycosylation" = "#336e2a",   # lila
    "C-Glycosylation" = "#eda407",   # orange
    "S-Glycosylation" = "#0cb05b",   # rosa
    "N-Lipidation" = "#a65628",      # brun
    "O-Lipidation" = "#a47ad6",      # cyan
    "S-Lipidation" = "#eda407",      # blå
    "GPI-anchor" = "#000000",        # svart
    "Phosphorylation" = "#08519c",   # djup blå (klar kontrast mot grön/lila/röd)
    "Acetylation" = "#d95f02",       # olivgrön
    "Methylation" = "#f781bf",       # lavendel
    "Sulfation" = "#b3a740",         # orange-brun
    "γ-Secretase Cleavage" = "#e41a1c",   # röd
    "β-Secretase Cleavage" = "#377eb8",   # blå
    "θ-Secretase Cleavage" = "#4daf4a",   # grön
    "α-Secretase Cleavage" = "#984ea3",   # lila
    "Trypsin Cleavage"     = "#ff7f00",   # orange
    "Cathepsin Cleavage"   = "#ffff33",   # gul
    "Thrombin Cleavage"    = "#a65628",   # brun
    "ACE Cleavage"         = "#f781bf",   # rosa
    "ADAM Cleavage"        = "#999999",   # grå
    "Caspase Cleavage"     = "#999999",   # turkos
    "Calpain Cleavage"     = "#8dd3c7",   # ljus turkos
    "MMP Cleavage"         = "#b3de69",    # limegrön
    "Citrulline" = "#ffd700"
  )) +
  
  scale_fill_manual(values = c(
    "Aβ40 (672-711)" = "pink",
    "Aβ42 (672-713)" = "lightblue"
  )) +
  
  labs(title = paste0("Post-Translational Modifications for ", prot$`Gene`, " (Uniprot ID: ", uni_choice, ")"), x = "Amino Acid Position", y = NULL, color = "", fill="") +
  scale_x_continuous(limits = c(0, seq_length), breaks = sort(unique(c(seq(0, seq_length, 50), seq_length))), expand = c(0, 0),sec.axis = dup_axis()) +
  scale_y_continuous(limits = c(0, 0.4), expand = c(0, 0) ,sec.axis = dup_axis()) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x.top = element_blank(),
    axis.title.x.top = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.line.y.right = element_line(color = "black", linewidth = 0.5),
    axis.line.x.top = element_line(color = "black", linewidth = 0.5),
    legend.position = "right",
    legend.margin = margin(t = -12.5, r = 0, b = -12.5, l = 0)
  )
