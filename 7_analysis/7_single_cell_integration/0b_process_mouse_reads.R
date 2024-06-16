library(liftOver)
library(GenomicAlignments)
chrlen_df = readr::read_tsv("../shared/chrom.sizes/processed/mm10.sizes", col_names = F)
# load("../shared/encode_compiled/Sep20/hg38_mm10_regions.rda")

# Convert mm10 regions into mm9. The single cell data is in mm9
mm10_regions@elementMetadata$id = hg38_mm10_regions$X2
ch = import.chain("mm10ToMm9.over.chain")
seqlevelsStyle(mm10_regions) = "UCSC"
mm9_regions_dup = unlist(liftOver(mm10_regions, ch))
mm9_regions = mm9_regions_dup[width(mm9_regions_dup) == 200]
saveRDS(mm9_regions, file = "mm9_regions.rds")

# Each bam file is proecssed in the same way as below
mm9_regions = readRDS("mm9_regions.rds")
cell_meta = readr::read_tsv("/dcl02/hongkai/data/wfang/mouse_atac_atlas/cell_metadata.txt")
reads = readGAlignmentPairs("/dcl02/hongkai/data/wfang/mouse_atac_atlas/BoneMarrow_62016.bam", use.names = T)
str_m = stringr::str_match(names(reads), "^(.*):(\\d+)#(.*)$")

cell_meta_subset = cell_meta[cell_meta$tissue.replicate == "BoneMarrow_62016", ]
cell_id = str_m[, 2]
cts_list = lapply(cell_meta_subset$cell, function(x) {
	reads_x = reads[cell_id %in% x]
	cts = ghelper::align2region(reads_x, mm9_regions, chrlen_df)
	as(cbind(cts), "dgCMatrix")
	})
saveRDS(cts_list, file = "temp_cts_list_bm_62016.rds")
