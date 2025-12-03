library(Hypatia)
library(tidyverse)

load("data/gbm/processing/isoform/gbm_input.rda")
load("data/rcc/processing/isoform/rcc_input.rda")
load("data/heart/processing/isoform/heart_input.rda")

# CreateSCE
gbm <- CreateSCE(
  countData = gbm_countData,
  colData = gbm_colData,
  rowData = gbm_rowData
)

rcc <- CreateSCE(
  countData = rcc_countData,
  colData = rcc_colData,
  rowData = rcc_rowData
)

heart <- CreateSCE(
  countData = heart_countData,
  colData = heart_colData,
  rowData = heart_rowData
)

rm(gbm_countData, gbm_colData, gbm_rowData, gbm_gtf)
rm(rcc_countData, rcc_colData, rcc_rowData, rcc_gtf)
rm(heart_countData, heart_colData, heart_rowData, heart_gtf)

# Set defaults
gbm <- SetGenes(gbm, "gene_name")
gbm <- SetTranscripts(gbm, "transcript_name")
gbm <- SetGroups(gbm, "cell_type")

rcc <- SetGenes(rcc, "gene_name")
rcc <- SetTranscripts(rcc, "transcript_name")
rcc <- SetGroups(rcc, "cell_type")

heart <- SetGenes(heart, "gene_name")
heart <- SetTranscripts(heart, "transcript_name")
heart <- SetGroups(heart, "cell_type")

# Filter transcripts
gbm <- SubsetTranscripts(gbm, subset = nCell >= 3)
gbm <- SubsetTranscripts(gbm, subset = ( (structural_category == "full-splice_match" & !RTS_stage) | (within_CAGE_peak & within_polyA_site & !RTS_stage) ))

rcc <- SubsetTranscripts(rcc, subset = nCell >= 3)
rcc <- SubsetTranscripts(rcc, subset = ( (structural_category == "full-splice_match" & !RTS_stage) | (within_CAGE_peak & within_polyA_site & !RTS_stage) ))

heart <- SubsetTranscripts(heart, subset = nCell >= 3)
heart <- SubsetTranscripts(heart, subset = ( (structural_category == "full-splice_match" & !RTS_stage) | (within_CAGE_peak & within_polyA_site & !RTS_stage) ))

# Filter cells
gbm <- SubsetCells(object = gbm, subset = nTranscript >= 50)
rcc <- SubsetCells(object = rcc, subset = nTranscript >= 50)
heart <- SubsetCells(object = heart, subset = nTranscript >= 50)

# Normalize
gbm <- NormalizeCounts(gbm, scale.factor = 10000)
rcc <- NormalizeCounts(rcc, scale.factor = 10000)
heart <- NormalizeCounts(heart, scale.factor = 10000)

# Save
saveRDS(gbm, file = "data/gbm/processing/isoform/gbm.rds")
saveRDS(rcc, file = "data/rcc/processing/isoform/rcc.rds")
saveRDS(heart, file = "data/heart/processing/isoform/heart.rds")
