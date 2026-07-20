##############################
## FCER1G gnomAD exon plot  ##
##############################

## Load required packages
library(readxl)    # read_xlsx
library(readr)     # read_csv
library(dplyr)     # data manipulation
library(ggplot2)   # plotting
library(scales)    # For trans_breaks, trans_format (log plot bg)
library(biomaRt)   # Ensembl exon coordinates
library(plotly)    # interactive plot
library(htmlwidgets) # saveWidget
library(ggrepel)   # for non-overlapping text labels
library(stringr)   # regex string manipulations for labels

## -------- Parameters --------

# Input gnomAD file (edit path if needed)
gnomad_file <- "gnomad.FCER1G_variants.xlsx"
gnomad_file <- "gnomAD_v4.1.1_ENSG00000158869_2026_07_01_20_18_59.genome.csv"
gnomad_file <- "gnomAD_v4.1.1_ENSG00000158869_2026_07_01_21_27_57.exome-genome.csv"


# Output files
pdf_out    <- "figure/FCER1G_Gnomad_variants.debug_txn.pdf"
png_out    <- "figure/FCER1G_Gnomad_variants.debug_txn.png"
pdf_noex_out    <- "figure/FCER1G_Gnomad_variants.no_txn.pdf"
png_noex_out    <- "figure/FCER1G_Gnomad_variants.no_txn.png"

html_out   <- "FCER1G_Gnomad_variants_plotly.html"

# Transcripts of interest (without version numbers for Ensembl)
tx_ids <- c("ENST00000289902", "ENST00000367992")

# Ensembl version
ensembl_ver = 115

# list latest Ensembl verison
bm_info = listEnsembl(GRCh = 38)
print(paste("Latest Ensembl Release:",listEnsembl(GRCh = 38)[bm_info$biomart == "genes", "version"]))
print(paste("Using Ensembl Release: ",listEnsembl(version=ensembl_ver)[bm_info$biomart == "genes", "version"]))

## -------- Get exon/transcript layout from Ensembl --------


ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = ensembl_ver
  #GRCh    = 38,   # change to 37 if needed, else defaults to 38
)

exons <- getBM(
  attributes = c(
    "ensembl_transcript_id",
    "chromosome_name",
    "exon_chrom_start",
    "exon_chrom_end",
    "rank"
  ),
  filters  = "ensembl_transcript_id",
  values   = tx_ids,
  mart     = ensembl
)

tx_info <- getBM(
  attributes = c(
    "ensembl_transcript_id",
    "ensembl_transcript_id_version",
    "ensembl_gene_id",
    "hgnc_symbol"
  ),
  filters  = "ensembl_transcript_id",
  values   = tx_ids,
  mart     = ensembl
)

exons <- exons %>%
  left_join(tx_info, by = "ensembl_transcript_id")

# 1) Pull unique transcript/version pairs
tx_versions <- exons %>%
  distinct(ensembl_transcript_id, ensembl_transcript_id_version)

# 2) Define human-readable names
tx_labels <- tibble::tibble(
  ensembl_transcript_id = c("ENST00000289902", "ENST00000367992"),
  base_label            = c("WT-Gamma",        "AS-Gamma")
)

# 3) Join and build final label: "WT-Gamma (ENST00000289902.2)"
tx_labels <- tx_labels %>%
  left_join(tx_versions, by = "ensembl_transcript_id") %>%
  mutate(isoform = paste0(base_label, " (", ensembl_transcript_id_version, ")"))

# 4) Attach back to exon data
exons <- exons %>%
  left_join(
    tx_labels %>% select(ensembl_transcript_id, isoform),
    by = "ensembl_transcript_id"
  )

# Create y-bands for the two isoform exon tracks BELOW the variant points
# WT-Gamma closer to zero, AS-Gamma a bit lower
exons <- exons %>%
  mutate(
    #ymin = ifelse(grepl("WT-Gamma", isoform), -0.04, -0.09),
    ymin = ifelse(grepl("WT-Gamma", isoform), 0.00001,0.0001),
    ymax = ymin *5
  )

## Determine region spanned by these transcripts
chr_fcer1g   <- unique(exons$chromosome_name)
region_start <- min(exons$exon_chrom_start)
region_end   <- max(exons$exon_chrom_end)

## -------- Read gnomAD variant data from earlier XLSX dump  --------
##
## would be nicer to read this directly from Gnomad
## Mais, j'ai des autres chats a fouetter.
##
ext <- tolower(tools::file_ext(gnomad_file))

gnomad_raw <- switch(
  ext,
  csv  = read_csv(gnomad_file, show_col_types = FALSE),
  xlsx = read_xlsx(gnomad_file),
  stop("Unsupported file extension: ", ext)
)

src_cols = c(
  'chr'        = 'Chromosome',
  'pos'       = 'Position',
  'af'         = 'Allele Frequency',
  'annotation' = 'VEP Annotation',
  'protein_cons' = 'Protein Consequence'
)
for(key in src_cols) {
  print(paste(key, ":", grep(pattern=key, names(gnomad_raw))))
}
gnomad <- gnomad_raw %>%
  rename(
    chr        = `Chromosome`,
    pos        = `Position`,
    af         = `Allele Frequency`,
    annotation = `VEP Annotation`,
    protein_cons = `Protein Consequence`
  ) %>%
  mutate(
    pos        = as.numeric(pos),
    af         = as.numeric(af),
    annotation = as.factor(annotation),
    protein_cons = as.character(protein_cons)
  )

# Restrict to FCER1G region & chromosome if needed
gnomad_region <- gnomad %>%
  filter(
    chr %in% chr_fcer1g,
    pos >= region_start,
    pos <= region_end,
  )

# "flip" any AF > 0.5 to be 1-AF (fixes our 99% "variant")
gnomad <- gnomad %>%
  mutate(
    af_orig = af,
    af = if_else(af > 0.3, 1.0 - af, af)
  )

# Filter out very rare variants
gnomad_region <- gnomad %>%
  filter(
    af > 5.0e-04,
    !(annotation %in% c("intron_variant", 
                        "non_coding_transcript_exon_variant",
                        "3_prime_UTR_variant"))
  )

## Safety check in case filtering nukes everything
if (nrow(gnomad_region) == 0) {
  warning("No gnomAD variants found in the FCER1G transcript region with this filter. Using full file instead.")
  gnomad_region <- gnomad
}

## ==============================================================
## ==                     Build ggplots                       ==
## ==============================================================

# ---------------------------------------------------------------
# internal usage plot to show exons and variants on genomic scale
# ---------------------------------------------------------------

min_y_exon <- min(exons$ymin)
max_y_var  <- max(gnomad_region$af, na.rm = TRUE)
y_upper    <- 1.0 # max_y_var * 1.1

p <- ggplot() +
  ## Exon rectangles
  geom_rect(
    data = exons,
    aes(
      xmin = exon_chrom_start,
      xmax = exon_chrom_end,
      ymin = ymin,
      ymax = ymax,
      fill = isoform
    ),
    color = "black",
    alpha = 0.8
  ) +
  ## Variant points
  geom_point(
    data = gnomad_region,
    aes(
      x     = pos,
      y     = af,
      color = annotation,
      annotation = annotation,
      protein_cons = protein_cons
    ),
    size  = 2,
    alpha = 0.8
  ) +
  geom_text_repel(
#    data = subset(gnomad_region, af > 0.01),
    data = gnomad_region,
    aes(
      x     = pos,
      y     = af,
      label = protein_cons,
      color = annotation
    ),
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
#  scale_y_continuous(
#    name   = "Allele frequency (gnomAD)",
#    limits = c(min_y_exon, y_upper),
#    expand = c(0, 0)
#  ) +
  scale_y_log10(
    name = " gnomAD Allele frequency (log10)",
    limits = c(0.00001, 1.0),
  ) +
  labs(
    x     = "Genomic position (bp)",
    title = "Gnomad variants FCER1G"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

## -------- Save PDF --------
ggsave(pdf_out, plot = p, width = 10, height = 5)
ggsave(png_out, plot = p, width = 10, height = 5)

## -------- Make Plotly interactive plot --------
p_plotly <- ggplotly(p, tooltip = c("x", "y", "annotation","protein_cons"))

saveWidget(
  widget = p_plotly,
  file   = html_out,
  selfcontained = TRUE
)

p_plotly
p

# ---------------------------------------------------------------
# ACTUAL FIGURE plot just variants on, still on a genomic scale
# 
# this will be labeled and combined with isoform diagram in pptx.
# ---------------------------------------------------------------

# format labels - bold/underline the AA that is the most common
# (some of these SNPs have AF > 0.5, and we plot them as (1.0-AF))
hgvs_parse <- str_match(gnomad_region$`HGVS Consequence`,
               "^p\\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3})$")

gnomad_region$protein_cons <- ifelse(
  gnomad_region$af_orig <= 0.5,
  sprintf('paste("p.",underline(bold("%s")),"%s%s")', hgvs_parse[, 2], hgvs_parse[, 3], hgvs_parse[, 4]),
  sprintf('paste("p.%s%s",underline(bold("%s")))', hgvs_parse[, 2], hgvs_parse[, 3], hgvs_parse[, 4])
)

p2 <-
  ggplot() +
  # ## Exon rectangles
  # geom_rect(
  #   data = exons,
  #   aes(
  #     xmin = exon_chrom_start,
  #     xmax = exon_chrom_end,
  #     ymin = ymin,
  #     ymax = ymax,
  #     fill = isoform
  #   ),
  #   color = "black",
  #   alpha = 0.8
  # ) +
  ## Variant points
  geom_point(
    data = gnomad_region,
    aes(
      x     = pos,
      y     = af,
      #color = "black",
      annotation = annotation,
      protein_cons = protein_cons
    ),
    size  = 5,
    alpha = 0.8
  ) +
  geom_text_repel(
#    data = subset(gnomad_region, af > 0.01),
    data = gnomad_region,
    aes(
      x     = pos,
      y     = af,
      label = protein_cons
      #color = annotation
    ),
    size = 5,
    max.overlaps = Inf,
    min.segment.length = 0,
    nudge_y      = 0.3,
    box.padding = 0.5,
    point.padding = 0.4,
    force = 5,
    parse = TRUE
  ) +
#  scale_y_continuous(
#    name   = "Allele frequency (gnomAD)",
#    limits = c(min_y_exon, y_upper),
#    expand = c(0, 0)
#  ) +
#  scale_y_log10(
#    name = " gnomAD Allele frequency (log10)",
#    limits = c(0.001, 1.0),
#  ) +
  scale_y_continuous(trans='log10',
                     breaks = c(0.001, 0.01, 0.1, 0.5, 1),
                     labels = c("0.001", "0.01", "0.1", "0.5", "1.0"),
                     minor_breaks = c(
                       0.0025, 0.005, 0.0075,
                       0.025, 0.05, 0.075,
                       0.2, 0.3, 0.4,
                       0.6, 0.7, 0.8, 0.9
                     ),
                     name = " gnomAD Allele frequency (log10)",
                     #minor_breaks = NULL, # remove minor
                     limits = c(0.001,1.0)
                      ) + # Format labels
  theme(panel.grid.major.y = element_line(color="black"), # Log-spaced major grid
        panel.grid.minor.y = element_line(color="gray95", linetype="dotted") # Minor log grid
  ) +
  labs(
    x     = "Genomic position (bp)",
    title = "Gnomad variants FCER1G"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5),
    legend.title = element_blank(),
    axis.title   = element_text(size = 24),
    axis.title.y = element_text(size = 24, margin = margin(r = 25)),
    axis.text    = element_text(size = 24),
    panel.grid.major.y = element_line(
      color = "grey80",
      linewidth = 0.8
    ),
    panel.grid.minor.y = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
   ) + 
  # hide some things
  theme(legend.position = "none") +
  theme(
    plot.title   = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank()
  ) + 
  theme(
    axis.ticks.y.minor = element_blank(),
    axis.ticks.x       = element_blank()
  ) + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
p2

## -------- Save PDF --------
ggsave(pdf_noex_out, plot = p2, width = 10, height = 5)
ggsave(png_noex_out, plot = p2, width = 10, height = 5)

## -------- Make Plotly interactive plot --------

p_plotly <- ggplotly(p2, tooltip = c("x", "y", "annotation","protein_cons"))
saveWidget(
  widget = p_plotly,
  file   = html_out,
  selfcontained = TRUE
)

