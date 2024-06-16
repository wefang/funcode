# this script save different versions of the raw DHS files
library(GenomicRanges)
library(dplyr)

dhs = readr::read_tsv("./metadata/ENCODE_DHS/ENCFF503GCK.tsv", col_types = "???c??????")
dhs_gr = makeGRangesFromDataFrame(dplyr::select(dhs, seqname, start, end))
dhs_gr$identifier = dhs$identifier
dhs_gr = dhs_gr[seqnames(dhs_gr) != "chrY"]

dhs_summit = makeGRangesFromDataFrame(transmute(dhs, seqname, start = summit, end = summit))
dhs_summit$identifier = dhs$identifier
dhs_summit = dhs_summit[!seqnames(dhs_summit) %in% c("chrY")]
dhs_summit_extended = resize(dhs_summit, width = 200, fix = 'center')

# dhs_summit = select(dhs, identifier, seqname, summit)
# readr::write_tsv(dhs_summit, "../metadata/indexed_dhs_summit.tsv", col_names = F)

# lifting over to hg19
# BiocManager::install("liftOver")
library(liftOver)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
seqlevelsStyle(dhs_gr) = "UCSC"
seqlevelsStyle(dhs_summit_extended) = "UCSC"

dhs_lifover = unlist(liftOver(dhs_gr, ch))
dhs_id_mapped = dhs_lifover$identifier
dhs_id_mapped_tb = table(dhs_id_mapped)
dhs_id_sel = names(dhs_id_mapped_tb)[dhs_id_mapped_tb == 1]
dhs_lifover = dhs_lifover[match(dhs_id_sel, dhs_lifover$identifier)]
dhs_hg19_lifover = dhs_lifover[width(dhs_lifover) == width(dhs_gr[match(dhs_id_sel, dhs_gr$identifier)])]

dhs_summit_lifover = unlist(liftOver(dhs_summit_extended, ch))
dhs_id_mapped = dhs_summit_lifover$identifier
dhs_id_mapped_tb = table(dhs_id_mapped)
dhs_id_sel = names(dhs_id_mapped_tb)[dhs_id_mapped_tb == 1]
dhs_summit_lifover = dhs_summit_lifover[match(dhs_id_sel, dhs_summit_lifover$identifier)]
dhs_summit_extended_hg19_lifover = dhs_summit_lifover[width(dhs_summit_lifover) == width(dhs_summit_extended[match(dhs_id_sel, dhs_summit_extended$identifier)])]

save(dhs_gr, dhs_summit_extended,
     dhs_hg19_lifover, dhs_summit_extended_hg19_lifover,
     file = "./metadata_processed/human_dhs_noY.rda")

# mouse DHS below
dhs = readr::read_tsv("./metadata/ENCODE_DHS//ENCFF910SRW.tsv", col_types = "???c??????")
dhs_gr = makeGRangesFromDataFrame(dplyr::select(dhs, chr, start, end))
dhs_gr$identifier = dhs$index_num
dhs_gr = dhs_gr[seqnames(dhs_gr) != "chrY"]
dhs_center_extended = resize(dhs_gr, width = 200, fix = "center")
save(dhs_gr, dhs_center_extended, file = "./metadata_processed//mouse_dhs_noY.rda")

# saving sets of region files
load("./metadata_processed/human_dhs_noY.rda")
human_regions = rlang::duplicate(dhs_gr)
load("./metadata_processed//mouse_dhs_noY.rda")
mouse_regions = rlang::duplicate(dhs_gr)
save(human_regions, mouse_regions, file = "./metadata_processed/human_mouse_dhs_raw.rda")

load("./metadata_processed/human_dhs_noY.rda")
human_regions = rlang::duplicate(dhs_summit_extended)
load("./metadata_processed//mouse_dhs_noY.rda")
mouse_regions = rlang::duplicate(dhs_center_extended)
save(human_regions, mouse_regions, file = "./metadata_processed/human_mouse_dhs_200bp.rda")

