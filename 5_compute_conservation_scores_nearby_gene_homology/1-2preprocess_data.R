## preprocess data: retain matched samples with anchor score > 0.5
library(dplyr)
setwd('/data/hji7/ENCODE4_unaligned/')

for(dataset in c('H3K27ac','H3K4me1','H3K4me3', 'chromatin')){ # chromatin
    print(dataset)

    #################
    ## Prepare data
    #################
    anchors_sel <- readRDS(paste0('./raw_data/anchors_',dataset,'_v1.rds'))
    print(dim(anchors_sel))
    # anchors_sel <- filter(anchors, score > 0.5)
    # saveRDS(anchors_sel,paste0('./proc_data/anchors_',dataset,'.rds'))

    # read in normalized matrix
    if(dataset == 'chromatin'){
        h <- readRDS(paste0('./raw_data/encodev4_human_chromatin_all_dhs_correct.rds')) # encodev4_human_chromatin_all_dhs_correct.rds
        m <- readRDS(paste0('./raw_data/encodev4_mouse_chromatin_all_dhs_correct.rds')) # encodev4_mouse_chromatin_all_dhs_correct.rds
    }else{
        h <- readRDS(paste0('./raw_data/Homo_sapiens_',dataset,'_norm_correct.rds'))
        m <- readRDS(paste0('./raw_data/Mus_musculus_',dataset,'_norm_correct.rds'))
    }

    print('==full matrix:')
    print(dim(h))
    print(dim(m))
    # matched anchors
    human_gene_mat <- h[, anchors_sel$accession1]
    mouse_gene_mat <- m[, anchors_sel$accession2]
    print('==anchor')
    print(dim(human_gene_mat))
    print(identical(ncol(human_gene_mat),nrow(anchors_sel)))
    saveRDS(human_gene_mat,paste0('./proc_data/human_',dataset,'.rds'))
    saveRDS(mouse_gene_mat,paste0('./proc_data/mouse_',dataset,'.rds'))

}

## H3K27ac (raw: human - 371; mouse - 97; processed: human/mouse - 208)
## H3K4me1 (raw: human - 323; mouse - 104; processed: human/mouse - 160)
## H3K4me3 (raw: human - 432; mouse - 103; processed: human/mouse - 168)
## chromatin (raw: human - 1226; mouse - 206; processed: human/mouse - 718)
