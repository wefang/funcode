## Create set of nulls for each category
## human_dist_group: 3
## mouse_dist_group: 3
## total: 9 groups (11,12,13,21,22,23,31,32,33)

library(dplyr)
library(parallel)
library(matrixStats)
setwd('/data/hji7/ENCODE4_unaligned/')
source('cnsv_score.R')

dataset <- 'chromatin' # can be replaced by other data modalities: H3K27ac, H3K4me1, H3K4me3
output_dir <- paste0('./result/',dataset,'/null9/')
if(!dir.exists(output_dir)) dir.create(output_dir)

human_gene_mat <- readRDS(paste0('./proc_data/human_',dataset,'.rds'))
mouse_gene_mat <- readRDS(paste0('./proc_data/mouse_',dataset,'.rds'))
anchors_sel <- readRDS(paste0('./raw_data/anchors_',dataset,'_v1.rds'))

## compute NULL
raw.human <- readRDS('./result/annot/comb/human_dhs_info.rds')
raw.mouse <- readRDS('./result/annot/comb/mouse_dhs_info.rds')
proc.human <- split(raw.human,factor(raw.human$human_dist_group))
proc.mouse <- split(raw.mouse,factor(raw.mouse$mouse_dist_group))


categories <- expand.grid(human_dist_group=1:3,mouse_dist_group=1:3) %>% mutate(group = paste0(human_dist_group,mouse_dist_group))
samplesize <- 1e5

prob_vec = seq(0.01, 1.0, by = 0.01)
set.seed(123)
human_quantile_vec = quantile(human_gene_mat[sample(nrow(human_gene_mat), 1e5), anchors_sel$accession1], prob_vec)
set.seed(123)
mouse_quantile_vec = quantile(mouse_gene_mat[sample(nrow(mouse_gene_mat), 1e5), anchors_sel$accession2], prob_vec)


k <- commandArgs(trailingOnly=TRUE)[1] %>% as.integer # 1:9
print(k)
sub.human <- proc.human[[categories$human_dist_group[k]]]
sub.mouse <- proc.mouse[[categories$mouse_dist_group[k]]]
set.seed(k)
human.select <- sample(sub.human$human_id,samplesize,replace = T)
set.seed(k+1234)
mouse.select <- sample(sub.mouse$mouse_id,samplesize,replace = T)
ws.dt <- sapply(1:samplesize,function(x){
    if(x %% 10000 == 0) print(x)
    weightedSpearman(human.select[x], mouse.select[x])
})
saveRDS(ws.dt, paste0(output_dir,'/spearman_',categories$group[k],'.rds'))

hval.dt <- compute_hval(human_gene_mat[human.select, anchors_sel$accession1],
                        mouse_gene_mat[mouse.select, anchors_sel$accession2],
                        human_quantile_vec,
                        mouse_quantile_vec,
                        prob_vec)
saveRDS(hval.dt, paste0(output_dir,'/hvalue_',categories$group[k],'.rds'))


