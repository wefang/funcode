## Calculate conservation scores
suppressMessages(library(dplyr))
suppressMessages(library(matrixStats))
setwd('/data/hji7/ENCODE4_unaligned/')
source('cnsv_score.R')

dataset <- 'chromatin' # can be replaced by other data modalities: H3K27ac, H3K4me1, H3K4me3
if(!dir.exists(paste0('./result/',dataset,'/'))) dir.create(paste0('./result/',dataset,'/'))
output_dir <- paste0('./result/',dataset,'/score/')
if(!dir.exists(output_dir)) dir.create(output_dir)

human_gene_mat <- readRDS(paste0('./proc_data/human_',dataset,'.rds'))
mouse_gene_mat <- readRDS(paste0('./proc_data/mouse_',dataset,'.rds'))
anchors_sel <- readRDS(paste0('./raw_data/anchors_',dataset,'_v1.rds'))

# h-value
prob_vec = seq(0.01, 1.0, by = 0.01)
set.seed(123)
human_quantile_vec = quantile(human_gene_mat[sample(nrow(human_gene_mat), 1e5), anchors_sel$accession1], prob_vec)
set.seed(123)
mouse_quantile_vec = quantile(mouse_gene_mat[sample(nrow(mouse_gene_mat), 1e5), anchors_sel$accession2], prob_vec)


i <- commandArgs(trailingOnly=TRUE)[1] %>% as.integer
if(!file.exists(paste0(output_dir,i,'.rds'))){
    print('===')
    print(i)
    annot <- readRDS(paste0('./result/annot_split/',i,'.rds')) %>%
        select(hid,human_id,mouse_id,human_symbol,mouse_symbol,human_dist, mouse_dist, group)

    ws.score <- sapply(1:nrow(annot), function(x){
        if(x %% 500 == 0) print(x)
        weightedSpearman(annot$human_id[x], annot$mouse_id[x])
    })

    h_value <- compute_hval(human_gene_mat[annot$human_id, anchors_sel$accession1],
                                mouse_gene_mat[annot$mouse_id, anchors_sel$accession2],
                                human_quantile_vec,
                                mouse_quantile_vec,
                                prob_vec)
    if(length(ws.score) == length(h_value) & length(h_value) == nrow(annot)){
        annot$weighted_spearman <- ws.score
        annot$h_value <- h_value
        saveRDS(as_tibble(annot),paste0(output_dir,i,'.rds'))
    }
}

