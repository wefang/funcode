# annotated DHS by gene reference
#### Process gene annotation file ####
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("./metadata/ENCODEv4_Reference/ENCFF159KBI.gtf.gz", format = 'gtf')

all.introns <- tidyIntrons(txdb)
all.exons <- tidyExons(txdb)
all.promoters <- promoters(txdb, upstream=1000, downstream=0)
all.tss <- promoters(txdb, upstream=0, downstream=1)
all.transcripts <- tidyTranscripts(txdb)
all.genes <- genes(txdb)
cds <- cdsBy(txdb, by='tx',use.names=T)
all.cds <- unlist(cds)

utr5_list = fiveUTRsByTranscript(txdb, use.names=T)
utr3_list = threeUTRsByTranscript(txdb, use.names=T)

all.utr5 = unlist(fiveUTRsByTranscript(txdb, use.names=T))
all.utr3 = unlist(threeUTRsByTranscript(txdb, use.names=T))

save(all.cds, all.utr5, all.utr3, all.promoters, all.tss, all.introns, all.genes,
     file = "./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_gene_hg38_gene_features.rda")
#### End Process ####

#### Repeat for mouse regions ####
txdb <- makeTxDbFromGFF("./metadata/ENCODEv4_Reference/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz", format = 'gtf')
all.introns <- tidyIntrons(txdb)
all.exons <- tidyExons(txdb)
all.promoters <- promoters(txdb, upstream=1000, downstream=0)
all.tss <- promoters(txdb, upstream=0, downstream=1)
all.transcripts <- tidyTranscripts(txdb)
all.genes <- genes(txdb)
cds <- cdsBy(txdb, by='tx',use.names=T)
all.cds <- unlist(cds)

utr5_list = fiveUTRsByTranscript(txdb, use.names=T)
utr3_list = threeUTRsByTranscript(txdb, use.names=T)

all.utr5 = unlist(fiveUTRsByTranscript(txdb, use.names=T))
all.utr3 = unlist(threeUTRsByTranscript(txdb, use.names=T))

save(all.cds, all.utr5, all.utr3, all.promoters, all.tss, all.introns, all.genes,
     file = "./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_gene_mm10_gene_features.rda")
# End Process

#### Annotating DHS ####
annotate_regions <- function(regions) {
        checklist = c('CDS','5UTR','3UTR','Promoter','Intron')
        gencode_list = GRangesList(all.cds, all.utr5, all.utr3, all.promoters, all.introns)
        region_tb = tibble(id = regions$identifier)
        for(i in 1:5){
                region_tb[paste0('is.',checklist[i])] = FALSE
                gencode = gencode_list[[i]]
                hits <- findOverlaps(regions, gencode)
                overlaps <- pintersect(regions[queryHits(hits)], gencode[subjectHits(hits)])
                percentOverlap <- width(overlaps) / width(regions[queryHits(hits)])
                hits <- hits[percentOverlap >= 0.5]
                ids = mcols(regions[queryHits(hits)])$identifier
                region_tb[region_tb$id %in% ids,][paste0('is.',checklist[i])] = TRUE
        }
        region_tb['is.Intergenic'] = TRUE
        for(i in 1:5){
                region_tb[['is.Intergenic']][region_tb[[paste0('is.',checklist[i])]]] = F
        }
        # remove any overlaps with all all genes
        hits <- findOverlaps(regions, all.genes)
        ids = mcols(regions[queryHits(hits)])$identifier
        region_tb[region_tb$id %in% ids,]['is.Intergenic'] = FALSE
        
        dist_tss = distanceToNearest(regions, all.tss, ignore.strand = T)
        region_tss_tb = tibble(id = regions$identifier[queryHits(dist_tss)],
                               nearest_tx = all.tss$tx_name[subjectHits(dist_tss)],
                               dist_tss = mcols(dist_tss)$distance)
        region_tb = left_join(region_tb, region_tss_tb)
        region_tb
}

load("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_gene_hg38_gene_features.rda")
load("./metadata_processed/indexed_aligned_regions.rda")
region_tb = annotate_regions(hg38_regions)
saveRDS(region_tb, "./metadata_processed/indexed_DHS_hg38_gencode_v29_gene_context.rds")

load("./metadata_processed/DHS/human_mouse_dhs_200bp.rda")
region_tb = annotate_regions(human_regions)
saveRDS(region_tb, "./metadata_processed/indexed_All_DHS_200bp_hg38_gencode_v29_gene_context.rds")
#### End Annotation ####

# annotated mouse regions
load("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_gene_mm10_gene_features.rda")
load("./metadata_processed/indexed_aligned_regions.rda")
region_tb = annotate_regions(mm10_regions)
saveRDS(region_tb, "./metadata_processed/indexed_DHS_mm10_gencode_vM21_gene_context.rds")

load("./metadata_processed/DHS/human_mouse_dhs_200bp.rda")
region_tb = annotate_regions(mouse_regions)
saveRDS(region_tb, "./metadata_processed/indexed_All_DHS_200bp_mm10_gencode_vM21_gene_context.rds")

#### Adding annotations to master region file ####
source("R_production/config.R")
regions = readRDS(region_file)
gene_feature_tb = full_join(readRDS("./metadata_processed/indexed_DHS_hg38_gencode_v29_gene_context.rds"),
                            readRDS("./metadata_processed/indexed_DHS_mm10_gencode_vM21_gene_context.rds"),
                            by = "id",
                            suffix = c("_hg38", "_mm10"))
regions = left_join(regions, gene_feature_tb)
saveRDS(regions, file = region_file)
#### End adding to mater region file
