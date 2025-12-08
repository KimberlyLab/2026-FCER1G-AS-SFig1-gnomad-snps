##############################
## FCER1G gnomAD exon plot  ##
##############################

## Load required packages
library(readxl)    # read_xlsx
library(dplyr)     # data manipulation
library(ggplot2)   # plotting
library(biomaRt)   # Ensembl exon coordinates
library(plotly)    # interactive plot
library(htmlwidgets) # saveWidget
library(ggrepel)   # for non-overlapping text labels

## -------- Parameters --------

# Input gnomAD file (edit path if needed)
gnomad_file <- "gnomad.FCER1G_canonical_1perc_snps.xlsx"

# Output files
pdf_out    <- "FCER1G_Gnomad_variants.pdf"
png_out    <- "FCER1G_Gnomad_variants.png"
html_out   <- "FCER1G_Gnomad_variants_plotly.html"

# Transcripts of interest (without version numbers for Ensembl)
tx_ids <- c("ENST00000289902", "ENST00000367992")

## -------- Get exon layout from Ensembl --------
## If your gnomAD data is GRCh37, set GRCh = 37 instead.
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  GRCh    = 38   # change to 37 if needed
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

# Map transcripts to human-readable isoform names
exons <- exons %>%
  mutate(
    isoform = case_when(
      ensembl_transcript_id == "ENST00000289902" ~ "WT-Gamma (ENST00000289902.2)",
      ensembl_transcript_id == "ENST00000367992" ~ "AS-Gamma (ENST00000367992.7)",
      TRUE ~ ensembl_transcript_id
    )
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

## -------- Read gnomAD variant data --------
gnomad_raw <- read_xlsx(gnomad_file)

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
    af = if_else(af > 0.5, 1.0 - af, af)
  )

# Filter out very rare variants
gnomad_region <- gnomad %>%
  filter(
    af > 0.005,
    !(annotation %in% c("intron_variant", "non_coding_transcript_exon_variant"))
  )

## Safety check in case filtering nukes everything
if (nrow(gnomad_region) == 0) {
  warning("No gnomAD variants found in the FCER1G transcript region with this filter. Using full file instead.")
  gnomad_region <- gnomad
}

## -------- Build ggplot --------

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
    data = subset(gnomad_region, af > 0.01),
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