library(seqinr)
library(tidyverse)
library(GenomicRanges)
compute_gc_content <- function(gr, genome_seq) {
        j = 0
        gr_width = width(gr)
        regions_gc = map_dbl(1:length(gr), function(i) {
                j <<- j + 1
                if (j %% 1e5 == 0) {
                        message(j)
                }
                x = gr[i]
                sum(genome_seq[[as.character(seqnames(x))]][start(x):end(x)] %in% c("g", "c")) / gr_width[i]        
        })
        tibble(id = gr$identifier,
               gc_content = regions_gc)
}

genome_seq = read.fasta("./metadata/ENCODEv4_Reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz")
load("./metadata_processed/indexed_aligned_regions.rda")
gc_tb = compute_gc_content(hg38_regions, genome_seq)
saveRDS(gc_tb, file = "metadata_processed/indexed_DHS_GRCh38_gc_content.rds")

load("./metadata_processed/DHS/human_mouse_dhs_200bp.rda")
gc_tb = compute_gc_content(human_regions, genome_seq)
saveRDS(gc_tb, file = "metadata_processed/indexed_All_DHS_200bp_GRCh38_gc_content.rds")


genome_seq = read.fasta("./metadata/ENCODEv4_Reference/mm10_no_alt_analysis_set_ENCODE.fasta.gz")
load("./metadata_processed/indexed_aligned_regions.rda")
gc_tb = compute_gc_content(mm10_regions, genome_seq)
saveRDS(gc_tb, file = "metadata_processed/indexed_DHS_mm10_gc_content.rds")

load("./metadata_processed/DHS/human_mouse_dhs_200bp.rda")
gc_tb = compute_gc_content(mouse_regions, genome_seq)
saveRDS(gc_tb, file = "metadata_processed/indexed_All_DHS_200bp_mm10_gc_content.rds")

# adding gc content to region file
gc_tb_hg38 = readRDS("metadata_processed/indexed_DHS_GRCh38_gc_content.rds")
gc_tb_mm10 = readRDS("metadata_processed/indexed_DHS_mm10_gc_content.rds")
gc_tb = full_join(gc_tb_hg38, gc_tb_mm10, by = "id", suffix = c("_hg38", "_mmm10"))
source("R_production/config.R")
regions = readRDS(region_file)
regions = left_join(regions, gc_tb)
saveRDS(regions, file = region_file)

