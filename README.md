# fcer1g_gnomad_snps

Supplemental Figure visualizing the frequency of coding SNPs in gene FCER1G's WT (ENST00000289902.2) and AS (ENST00000367992.7) transcripts with an alternate allele frequence (AF) > 0.001, based on AFs data exported from the Genome Aggregation Database (gnomAD) v4.1 (Broad Institute of MIT and Harvard, accessed April 28, 2025, [https://gnomad.broadinstitute.org](https://gnomad.broadinstitute.org), RRID: [SCR_014964]()https://scicrunch.org/resolver/SCR_014964). For methodological details on the gnomAD resource, see [Karczewski et al. (2020)](https://pubmed.ncbi.nlm.nih.gov/32461654/) 

Protocol: 
  1. Manual export from gnomAD was saved to the file [gnomad.FCER1G_canonical_1perc_snps.xlsx](gnomad.FCER1G_canonical_1perc_snps.xlsx).
  2. xlsx was processed and visualization creatd with R script [FCER1G_snps_exome_plot.R](FCER1G_snps_exome_plot.R) producing [figure/FCER1G_Gnomad_variants.no_txn.png](figure/FCER1G_Gnomad_variants.no_txn.png)
  3. Resulting figure was annotated in Microsoft PowerPoint [FCER1G_Gnomad_variants.paste_up.pptx](FCER1G_Gnomad_variants.paste_up.pptx) slide 4
