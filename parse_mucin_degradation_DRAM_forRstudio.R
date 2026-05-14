# RStudio-friendly version.
# Edit these values if needed, then run the whole script.

suppressPackageStartupMessages(library(tidyverse))

INPUT <- "data/annotations.tsv"
OUTDIR <- "results"

message("🔎 Reading DRAM annotations from: ", INPUT)

if (!file.exists(INPUT)) {
  stop("Input file not found: ", INPUT)
}

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

annotations <- readr::read_tsv(INPUT, show_col_types = FALSE, progress = FALSE, name_repair = "minimal")

first_column_name <- names(annotations)[1]
if (is.na(first_column_name) || first_column_name == "" || first_column_name == "...1") {
  names(annotations)[1] <- "GeneID"
} else if (first_column_name != "GeneID") {
  names(annotations)[1] <- "GeneID"
}

if (!"cazy_hits" %in% names(annotations)) {
  annotations$cazy_hits <- NA_character_
}
if (!"cazy_ids" %in% names(annotations)) {
  annotations$cazy_ids <- NA_character_
}

annotation_columns <- c(
  "gene", "product", "description", "kegg_hit", "pfam_hits", "ec",
  "ko_id", "peptidase_id", "peptidase_family", "peptidase_hit",
  "cazy_best_hit", "cazy_subfam_ec"
)
annotation_columns <- intersect(annotation_columns, names(annotations))

if (length(annotation_columns) == 0) {
  annotations$annotation_text <- NA_character_
} else {
  annotations <- annotations %>%
    mutate(
      annotation_text = do.call(
        paste,
        c(across(all_of(annotation_columns), ~ replace_na(as.character(.x), "")), sep = " | ")
      ),
      annotation_text = str_squish(annotation_text),
      annotation_text = na_if(annotation_text, "")
    )
}

annotations <- annotations %>%
  mutate(
    GeneID = as.character(GeneID),
    cazy_hits = as.character(cazy_hits),
    cazy_ids = as.character(cazy_ids),
    cazy_text = str_squish(paste(replace_na(cazy_hits, ""), replace_na(cazy_ids, ""), sep = " | "))
  )

cazy_markers <- tribble(
  ~matched_marker, ~functional_category, ~confidence,
  "GH33",  "sialic_acid_removal",          "high",
  "GH29",  "fucose_removal",               "high",
  "GH95",  "fucose_removal",               "high",
  "GH101", "O_glycan_core_degradation",    "high",
  "GH129", "O_glycan_core_degradation",    "high",
  "GH2",   "terminal_glycan_trimming",     "medium",
  "GH20",  "terminal_glycan_trimming",     "medium",
  "GH35",  "terminal_glycan_trimming",     "medium",
  "GH42",  "terminal_glycan_trimming",     "medium",
  "GH84",  "terminal_glycan_trimming",     "medium",
  "GH85",  "O_glycan_core_degradation",    "medium",
  "GH89",  "terminal_glycan_trimming",     "medium",
  "GH16",  "LacNAc_polyLacNAc_degradation", "medium",
  "GH98",  "LacNAc_polyLacNAc_degradation", "medium"
)

keyword_markers <- tribble(
  ~matched_marker, ~keyword_pattern, ~functional_category, ~confidence,
  "mucin",                    "\\bmucin\\b",                       "other_mucin_related",          "high",
  "sialidase",                "\\bsialidase\\b",                   "sialic_acid_removal",          "high",
  "neuraminidase",            "\\bneuraminidase\\b",               "sialic_acid_removal",          "high",
  "fucosidase",               "\\bfucosidase\\b",                  "fucose_removal",               "high",
  "sulfatase",                "\\bsulfatase\\b",                   "sulfate_removal",              "high",
  "arylsulfatase",            "\\barylsulfatase\\b",               "sulfate_removal",              "high",
  "beta-galactosidase",       "\\bbeta[- ]galactosidase\\b",       "terminal_glycan_trimming",     "medium",
  "alpha-galactosidase",      "\\balpha[- ]galactosidase\\b",      "terminal_glycan_trimming",     "medium",
  "hexosaminidase",           "\\bhexosaminidase\\b",              "terminal_glycan_trimming",     "medium",
  "N-acetylglucosaminidase",  "\\bn[- ]acetylglucosaminidase\\b",  "terminal_glycan_trimming",     "medium",
  "N-acetylgalactosaminidase","\\bn[- ]acetylgalactosaminidase\\b","O_glycan_core_degradation",    "medium",
  "glycopeptidase",           "\\bglycopeptidase\\b",              "glycopeptide_mucinase",        "high",
  "O-glycopeptidase",         "\\bo[- ]glycopeptidase\\b",         "glycopeptide_mucinase",        "high",
  "mucinase",                 "\\bmucinase\\b",                    "glycopeptide_mucinase",        "high",
  "nanH",                     "\\bnanH\\b",                        "sialic_acid_removal",          "high",
  "fuc",                      "\\bfuc[A-Za-z0-9_-]*\\b",           "sugar_catabolism_transport",   "low",
  "nag",                      "\\bnag[A-Za-z0-9_-]*\\b",           "sugar_catabolism_transport",   "low",
  "aga",                      "\\baga[A-Za-z0-9_-]*\\b",           "sugar_catabolism_transport",   "low",
  "gal",                      "\\bgal[A-Za-z0-9_-]*\\b",           "sugar_catabolism_transport",   "low"
)

message("🧬 Searching CAZy marker families")

cazy_hits_long <- cazy_markers %>%
  mutate(marker_regex = paste0("(^|[^A-Za-z0-9])", matched_marker, "([^0-9]|$)")) %>%
  pmap_dfr(function(matched_marker, functional_category, confidence, marker_regex) {
    annotations %>%
      filter(str_detect(cazy_text, regex(marker_regex, ignore_case = TRUE))) %>%
      transmute(
        GeneID,
        matched_marker,
        marker_type = "CAZy_family",
        functional_category,
        confidence,
        cazy_hits,
        cazy_ids,
        annotation_text
      )
  })

message("📝 Searching annotation keywords")

keyword_hits_long <- keyword_markers %>%
  pmap_dfr(function(matched_marker, keyword_pattern, functional_category, confidence) {
    annotations %>%
      filter(str_detect(replace_na(annotation_text, ""), regex(keyword_pattern, ignore_case = TRUE))) %>%
      transmute(
        GeneID,
        matched_marker,
        marker_type = "keyword",
        functional_category,
        confidence,
        cazy_hits,
        cazy_ids,
        annotation_text
      )
  })

gene_hits <- bind_rows(cazy_hits_long, keyword_hits_long) %>%
  distinct() %>%
  group_by(GeneID) %>%
  filter(!(confidence == "low" & any(confidence %in% c("high", "medium")))) %>%
  ungroup() %>%
  arrange(GeneID, desc(confidence), marker_type, matched_marker)

confidence_rank <- c(high = 3, medium = 2, low = 1)

gene_hits_wide <- gene_hits %>%
  mutate(confidence_score = confidence_rank[confidence]) %>%
  group_by(GeneID) %>%
  summarise(
    matched_markers = paste(sort(unique(matched_marker)), collapse = ";"),
    marker_types = paste(sort(unique(marker_type)), collapse = ";"),
    functional_categories = paste(sort(unique(functional_category)), collapse = ";"),
    confidence = names(which.max(tapply(confidence_score, confidence, max))),
    cazy_hits = first(na.omit(cazy_hits), default = NA_character_),
    cazy_ids = first(na.omit(cazy_ids), default = NA_character_),
    annotation_text = first(na.omit(annotation_text), default = NA_character_),
    .groups = "drop"
  ) %>%
  arrange(GeneID)

summary_by_marker <- gene_hits %>%
  count(matched_marker, marker_type, confidence, sort = TRUE, name = "n_hits") %>%
  arrange(desc(n_hits), matched_marker)

summary_by_category <- gene_hits %>%
  count(functional_category, sort = TRUE, name = "n_hits")

summary_by_confidence <- gene_hits %>%
  count(confidence, sort = TRUE, name = "n_hits") %>%
  mutate(confidence = factor(confidence, levels = c("high", "medium", "low"))) %>%
  arrange(confidence)

readme_text <- c(
  "Mucin degradation marker screen from DRAM annotations",
  "",
  "Searched high-confidence CAZy families: GH33, GH29, GH95, GH101, GH129.",
  "Searched supporting CAZy families: GH2, GH20, GH35, GH42, GH84, GH85, GH89, GH16, GH98.",
  "",
  "Keyword searches included mucin, sialidase, neuraminidase, fucosidase, sulfatase, arylsulfatase, beta-galactosidase, alpha-galactosidase, hexosaminidase, N-acetylglucosaminidase, N-acetylgalactosaminidase, glycopeptidase, O-glycopeptidase, mucinase, nanH, fuc, nag, aga, and gal.",
  "",
  "Confidence assignment:",
  "- high: GH33, GH29, GH95, GH101, GH129, or clear mucin/sialidase/fucosidase/sulfatase/glycopeptidase/O-glycopeptidase/mucinase keywords.",
  "- medium: supporting CAZy families and specific terminal glycan trimming keywords.",
  "- low: broad sugar metabolism or transport terms only, such as fuc, nag, aga, or gal without stronger mucin evidence.",
  "",
  "Interpretation note: this is a screening result, not proof of mucin degradation.",
  "Stronger evidence requires gene context, secretion/localization evidence, transporters, and growth or experimental validation."
)

message("💾 Writing result tables to: ", OUTDIR)

readr::write_tsv(gene_hits, file.path(OUTDIR, "mucin_degradation_gene_hits.tsv"))
readr::write_tsv(gene_hits_wide, file.path(OUTDIR, "mucin_degradation_gene_hits_wide.tsv"))
readr::write_tsv(summary_by_marker, file.path(OUTDIR, "mucin_degradation_summary_by_marker.tsv"))
readr::write_tsv(summary_by_category, file.path(OUTDIR, "mucin_degradation_summary_by_category.tsv"))
readr::write_tsv(summary_by_confidence, file.path(OUTDIR, "mucin_degradation_summary_by_confidence.tsv"))
writeLines(readme_text, file.path(OUTDIR, "mucin_degradation_README.txt"))

message("✅ Done")
message("Genes with at least one mucin-related marker: ", n_distinct(gene_hits$GeneID))
message("Total marker hits: ", nrow(gene_hits))
message("High-confidence hits: ", sum(gene_hits$confidence == "high"))
message("Medium-confidence hits: ", sum(gene_hits$confidence == "medium"))
message("Low-confidence hits: ", sum(gene_hits$confidence == "low"))
