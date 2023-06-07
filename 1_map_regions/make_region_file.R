# this script initializes a region file that is read and updated downstream
# initial alignable region file is based on Chaoran's input
library(tidyverse)
load("./processed_data/indexed_dhs_hg38_mm10_mapped/regions/hg38_mm10_regions.rda")
hg38_regions$identifier = regions$id
mm10_regions$identifier = regions$id
save(hg38_regions, mm10_regions, file = "./metadata_processed/indexed_aligned_regions.rda")
# TODO: what is the last column here?
hg38_mm10_regions = hg38_mm10_regions %>% rename(index = X1, id = X2, hg38_region = X3, mm10_region = X4, axnet_id = X5, gap = X6)
# adding phastCons and percent similarity from Chaoran
phast_scores = readr::read_tsv("./processed_data/indexed_dhs_hg38_mm10_mapped/regions/hg38_phastcon&phylop.txt", col_types = c(id = "c", hg38_region = "c"))
hg38_mm10_regions = left_join(hg38_mm10_regions, select(phast_scores, -hg38_region), by = "id")
saveRDS(hg38_mm10_regions, file = "./output/indexed_dhs_mapped_regions.rds")
