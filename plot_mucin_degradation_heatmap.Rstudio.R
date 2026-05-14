# RStudio-friendly version.
# Edit these values if needed, then run the whole script.

required_packages <- c("tidyverse", "pheatmap", "RColorBrewer")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them in the R environment, for example:\n",
    "install.packages(c('", paste(missing_packages, collapse = "', '"), "'))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

HITS_FILE <- "results/mucin_degradation_gene_hits.tsv"
HITS_WIDE_FILE <- "results/mucin_degradation_gene_hits_wide.tsv"
FRESHWATER_GTDB_FILE <- "data/freshwater.gtdbtk.bac120.summary.tsv"
MUCUS_GTDB_FILE <- "data/mucus.gtdbtk.bac120.summary.tsv"
SALTWATER_GTDB_FILE <- "data/saltwater.gtdbtk.bac120.summary.tsv"
OUTDIR <- "results"
USE_BINARY <- FALSE

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
font_cache_dir <- file.path(tempdir(), "fontconfig-cache")
dir.create(font_cache_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(XDG_CACHE_HOME = font_cache_dir)

input_files <- c(
  HITS_FILE,
  HITS_WIDE_FILE,
  FRESHWATER_GTDB_FILE,
  MUCUS_GTDB_FILE,
  SALTWATER_GTDB_FILE
)
missing_input_files <- input_files[!file.exists(input_files)]
if (length(missing_input_files) > 0) {
  stop("Missing input file(s): ", paste(missing_input_files, collapse = ", "), call. = FALSE)
}

message("📥 Reading mucin hit tables...")
mucin_hits <- readr::read_tsv(HITS_FILE, show_col_types = FALSE)
invisible(readr::read_tsv(HITS_WIDE_FILE, show_col_types = FALSE, progress = FALSE))

if (nrow(mucin_hits) == 0) {
  stop("No mucin hits found in ", HITS_FILE, call. = FALSE)
}

required_hit_columns <- c("GeneID", "matched_marker", "functional_category", "confidence")
missing_hit_columns <- setdiff(required_hit_columns, names(mucin_hits))
if (length(missing_hit_columns) > 0) {
  stop("Mucin hit table is missing column(s): ", paste(missing_hit_columns, collapse = ", "), call. = FALSE)
}

message("🧬 Parsing GTDB-Tk taxonomy...")

parse_gtdb <- function(file, treatment_name) {
  read_tsv(file, show_col_types = FALSE) %>%
    select(user_genome, classification) %>%
    separate(
      classification,
      into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = ";",
      fill = "right"
    ) %>%
    mutate(across(Domain:Species, ~ str_remove(.x, "^[a-z]__"))) %>%
    mutate(Treatment = treatment_name) %>%
    rename(MAG = user_genome)
}

clean_taxonomy_value <- function(x) {
  x <- as.character(x)
  x <- if_else(is.na(x) | str_trim(x) == "" | x == "NA" | x == "N/A", "Unclassified", x)
  x
}

GTDB <- bind_rows(
  parse_gtdb(FRESHWATER_GTDB_FILE, "freshwater"),
  parse_gtdb(MUCUS_GTDB_FILE, "mucus"),
  parse_gtdb(SALTWATER_GTDB_FILE, "saltwater")
) %>%
  mutate(across(c(Domain, Phylum, Class, Order, Family, Genus, Species), clean_taxonomy_value)) %>%
  distinct(MAG, .keep_all = TRUE) %>%
  arrange(Treatment, Phylum, Class, Order, Family, Genus, Species, MAG)

GTDB <- GTDB %>%
  mutate(
    Label = case_when(
      is.na(Species) | str_trim(Species) == "" | Species == "Unclassified" ~ Genus,
      TRUE ~ Species
    )
  ) %>%
  mutate(
    Label_unique = ave(
      Label, Label,
      FUN = function(x) if (length(x) > 1) paste(x, seq_along(x)) else x
    )
  )

message("🔗 Assigning genes to MAGs...")

assign_mag_from_gene <- function(gene_id, mag_names) {
  matches <- mag_names[str_starts(gene_id, fixed(mag_names))]
  if (length(matches) == 0) {
    return(NA_character_)
  }
  matches[which.max(nchar(matches))]
}

mag_names <- GTDB$MAG

mucin_hits_with_mag <- mucin_hits %>%
  mutate(MAG = map_chr(GeneID, assign_mag_from_gene, mag_names = mag_names))

unassigned_genes <- mucin_hits_with_mag %>%
  filter(is.na(MAG)) %>%
  distinct(GeneID, matched_marker, marker_type, functional_category, confidence)

readr::write_tsv(unassigned_genes, file.path(OUTDIR, "mucin_degradation_unassigned_genes.tsv"))
if (nrow(unassigned_genes) > 0) {
  message("⚠️  Unassigned genes written to results/mucin_degradation_unassigned_genes.tsv: ", n_distinct(unassigned_genes$GeneID))
}

assigned_hits <- mucin_hits_with_mag %>%
  filter(!is.na(MAG))

if (nrow(assigned_hits) == 0) {
  stop("No mucin hits could be assigned to MAGs after GeneID parsing.", call. = FALSE)
}

marker_order <- c(
  "GH33", "GH29", "GH95", "GH101", "GH129",
  "GH2", "GH20", "GH35", "GH42", "GH84", "GH85", "GH89", "GH16", "GH98",
  "sulfatase", "arylsulfatase", "mucin", "sialidase", "neuraminidase", "fucosidase",
  "beta-galactosidase", "alpha-galactosidase", "hexosaminidase",
  "N-acetylglucosaminidase", "N-acetylgalactosaminidase",
  "glycopeptidase", "O-glycopeptidase", "mucinase", "nanH",
  "fuc", "nag", "aga", "gal"
)

present_markers <- unique(assigned_hits$matched_marker)
ordered_markers <- c(intersect(marker_order, present_markers), setdiff(sort(present_markers), marker_order))

message("🔥 Building mucin degradation heatmap...")

marker_counts <- assigned_hits %>%
  group_by(MAG, matched_marker) %>%
  summarise(gene_count = n_distinct(GeneID), .groups = "drop") %>%
  complete(
    MAG = GTDB$MAG,
    matched_marker = ordered_markers,
    fill = list(gene_count = 0)
  )

heatmap_matrix_table <- marker_counts %>%
  mutate(matched_marker = factor(matched_marker, levels = ordered_markers)) %>%
  arrange(match(MAG, GTDB$MAG), matched_marker) %>%
  pivot_wider(names_from = matched_marker, values_from = gene_count, values_fill = 0) %>%
  arrange(match(MAG, GTDB$MAG))

mat_counts <- heatmap_matrix_table %>%
  column_to_rownames("MAG") %>%
  as.matrix()

if (USE_BINARY) {
  mat_plot <- ifelse(mat_counts > 0, 1, 0)
  legend_title <- "Presence"
} else {
  mat_plot <- log10(mat_counts + 1)
  legend_title <- "log10(count + 1)"
}

rownames(mat_plot) <- GTDB$Label_unique[match(rownames(mat_plot), GTDB$MAG)]

marker_categories <- assigned_hits %>%
  count(matched_marker, functional_category, confidence, name = "n_hits") %>%
  mutate(confidence_rank = recode(confidence, high = 3L, medium = 2L, low = 1L, .default = 0L)) %>%
  arrange(matched_marker, desc(n_hits), desc(confidence_rank), functional_category) %>%
  group_by(matched_marker) %>%
  summarise(
    functional_category = first(functional_category),
    confidence = confidence[which.max(confidence_rank)],
    .groups = "drop"
  ) %>%
  mutate(matched_marker = factor(matched_marker, levels = ordered_markers)) %>%
  arrange(matched_marker) %>%
  mutate(matched_marker = as.character(matched_marker))

readr::write_tsv(marker_categories, file.path(OUTDIR, "mucin_degradation_heatmap_marker_categories.tsv"))

row_annotation <- GTDB %>%
  select(MAG, Treatment, Phylum, Class, Order) %>%
  arrange(match(MAG, rownames(mat_plot)))

functional_category_order <- marker_categories %>%
  pull(functional_category) %>%
  unique()

column_annotation <- marker_categories %>%
  select(matched_marker, functional_category) %>%
  mutate(functional_category = factor(functional_category, levels = functional_category_order)) %>%
  column_to_rownames("matched_marker")

readr::write_tsv(heatmap_matrix_table, file.path(OUTDIR, "mucin_degradation_heatmap_matrix.tsv"))
readr::write_tsv(
  row_annotation %>% select(MAG, Treatment, Order),
  file.path(OUTDIR, "mucin_degradation_heatmap_annotation.tsv")
)

row_annotation_for_plot <- row_annotation %>%
  select(MAG, Treatment, Order) %>%
  mutate(MAG = GTDB$Label_unique[match(MAG, GTDB$MAG)]) %>%
  column_to_rownames("MAG")

make_palette <- function(values, palette_name, sort_values = TRUE) {
  values <- unique(as.character(values))
  values <- values[!is.na(values)]
  if (sort_values) {
    values <- sort(values)
  }
  max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
  if (length(values) <= max_colors) {
    colors <- RColorBrewer::brewer.pal(max(3, length(values)), palette_name)[seq_along(values)]
  } else {
    colors <- colorRampPalette(RColorBrewer::brewer.pal(max_colors, palette_name))(length(values))
  }
  names(colors) <- values
  colors
}

treatment_colors <- c(
  mucus = "#a2ecc2ff",
  freshwater = "#66b0e4ff",
  saltwater = "#f3716fff"
)
treatment_colors <- treatment_colors[intersect(names(treatment_colors), unique(row_annotation_for_plot$Treatment))]

annotation_colors <- list(
  Treatment = treatment_colors,
  functional_category = make_palette(functional_category_order, "Dark2", sort_values = FALSE)
)

unique_orders <- unique(row_annotation_for_plot$Order)
paired_max_colors <- RColorBrewer::brewer.pal.info["Paired", "maxcolors"]
if (length(unique_orders) <= paired_max_colors) {
  annotation_colors$Order <- make_palette(row_annotation_for_plot$Order, "Paired")
} else {
  message("⚠️  More Orders than the Paired palette supports; Order will use automatic pheatmap colors.")
}

heatmap_colors <- rev(hcl.colors(500, "Purples 3"))
plot_width <- max(10, 0.35 * ncol(mat_plot) + 5)
plot_height <- max(8, 0.08 * nrow(mat_plot) + 4)

heatmap_plot <- pheatmap(
  mat_plot,
  color = heatmap_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annotation_for_plot,
  annotation_col = column_annotation,
  annotation_colors = annotation_colors,
  show_rownames = F,
  angle_col = 45,
  border_color = NA,
  main = paste("Mucin-degradation marker counts (", legend_title, ")", sep = ""),
  silent = TRUE
)

save_pheatmap <- function(pheatmap_result, filename, canvas_width = 16, canvas_height = 10, dpi = 300) {
  if (file.exists(filename)) {
    file.remove(filename)
  }
  heatmap_canvas <- grid::grobTree(
    pheatmap_result$gtable,
    vp = grid::viewport(
      x = grid::unit(0.02, "npc"),
      y = grid::unit(0.5, "npc"),
      width = grid::unit(plot_width, "in"),
      height = grid::unit(plot_height, "in"),
      just = c("left", "center")
    )
  )
  ggplot2::ggsave(
    filename = filename,
    plot = heatmap_canvas,
    width = canvas_width,
    height = canvas_height,
    units = "in",
    dpi = dpi,
    limitsize = FALSE,
    bg = "white"
  )
}

save_pheatmap(
  heatmap_plot,
  file.path(OUTDIR, "mucin_degradation_heatmap.pdf"),
  canvas_width = 16,
  canvas_height = 10
)

save_pheatmap(
  heatmap_plot,
  file.path(OUTDIR, "mucin_degradation_heatmap.png"),
  canvas_width = 16,
  canvas_height = 10
)

message("✅ Done! Heatmap written to results/")
message("MAGs in heatmap: ", nrow(mat_plot))
message("Markers in heatmap: ", ncol(mat_plot))
message("Assigned mucin marker hits: ", nrow(assigned_hits))
