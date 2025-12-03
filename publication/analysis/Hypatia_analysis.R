library(Hypatia)
library(tidyverse)

gbm <- readRDS(file = "data/gbm/processing/isoform/gbm.rds")
rcc <- readRDS(file = "data/rcc/processing/isoform/rcc.rds")
heart <- readRDS(file = "data/heart/processing/isoform/heart.rds")

## DEI
gbm_dei <- RunDEI(gbm, only.pos = TRUE)
write_tsv(gbm_dei, file = "data/gbm/processing/isoform/deis.tsv")

rcc_dei <- RunDEI(rcc, only.pos = TRUE)
write_tsv(rcc_dei, file = "data/rcc/processing/isoform/deis.tsv")

heart_dei <- RunDEI(heart, only.pos = TRUE)
write_tsv(heart_dei, file = "data/heart/processing/isoform/deis.tsv")

## DIU
gbm_diu <- RunDIU(gbm, min.gene.pct = 0.05)
write_tsv(gbm_diu$stats, file = "data/gbm/processing/isoform/dius.tsv")

rcc_diu <- RunDIU(rcc, min.gene.pct = 0.05)
write_tsv(rcc_diu$stats, file = "data/rcc/processing/isoform/dius.tsv")

heart_diu <- RunDIU(heart, min.gene.pct = 0.05)
write_tsv(heart_diu$stats, file = "data/heart/processing/isoform/dius.tsv")

## DIV
gbm_div <- RunDIV(gbm, min.gene.pct = 0.05)
write_tsv(gbm_div$data, file = "data/gbm/processing/isoform/divs.tsv")

rcc_div <- RunDIV(rcc, min.gene.pct = 0.05)
write_tsv(rcc_div$data, file = "data/rcc/processing/isoform/divs.tsv")

heart_div <- RunDIV(heart, min.gene.pct = 0.05)
write_tsv(heart_div$data, file = "data/heart/processing/isoform/divs.tsv")
