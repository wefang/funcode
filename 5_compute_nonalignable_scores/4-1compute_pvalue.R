## compute p-value
library(data.table)
library(dplyr)
library(qvalue)
setwd('/data/hji7/ENCODE4_unaligned/')

dataset <- 'chromatin' # can be replaced by other data modalities: H3K27ac, H3K4me1, H3K4me3
output_dir <- paste0('./result/',dataset,'/pvalue/')
if(!dir.exists(output_dir)) dir.create(output_dir)

alignable <- fread('./raw_data/aligned_dhs_pairs.tsv') %>%
    mutate(alignable = 1, mouse_dhs = as.character(mouse_dhs), human_dhs = as.character(human_dhs)) %>%
    dplyr::rename(mouse_id = mouse_dhs, human_id = human_dhs)

run.dt <- data.frame(id = 1:131872)
n <- nrow(run.dt)
chunk <- 5e3
run.dt$run  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]

run <- commandArgs(trailingOnly=TRUE)[1] %>% as.numeric
print(paste0('=== run: ', run))
sub <- run.dt$id[run.dt$run == run]
print(length(sub))
print(range(sub))

comb.ls <- lapply(sub,function(f){
    print(f)
    readRDS(paste0('./result/chromatin/score/',f,'.rds'))
})

sapply(comb.ls,nrow) %>% sum %>% print # 50000000
comb <- do.call(rbind,comb.ls)
dim(comb)
# saveRDS(comb,paste0(output_dir,'comb_sub.rds'))

# split by category
cat.ls <- split(comb,factor(comb$group))
# saveRDS(cat.ls,paste0(output_dir,'comb_sub_bycategory.rds'))
sapply(cat.ls,nrow) %>% print

# p-value
new.ls <- lapply(1:9,function(k){
    print(paste0('===',k))

    category <- c('11','12','13','21','22','23','31','32','33')[k]
    print(category)
    spearman.null <- readRDS(paste0('./result/chromatin/null9/spearman_',category,'.rds'))
    hvalue.null <- readRDS(paste0('./result/chromatin/null9//hvalue_',category,'.rds'))
    samplesize <- length(spearman.null)

    dt <- cat.ls[[category]]
    print(paste0('total:', nrow(dt)))
    dt$weighted_spearman_p <- sapply(dt$weighted_spearman, function(x) sum(spearman.null >= x))/samplesize
    dt$h_value_p <- sapply(dt$h_value, function(x) sum(hvalue.null >= x))/samplesize

    dt <- suppressMessages(left_join(dt,alignable) %>%
                               mutate(alignable = ifelse(is.na(alignable),0,alignable)))

    saveRDS(dt,paste0(output_dir,category,'_',run,'.rds'))
})

