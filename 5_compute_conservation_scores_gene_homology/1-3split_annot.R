## split annotation files into smaller chunks
library(dplyr)
library(tidyr)
library(parallel)
setwd('/data/hji7/ENCODE4_unaligned/')

input_dir <- './result/annot/'
output_dir <- './result/annot_split/'
if(!dir.exists(output_dir)) dir.create(output_dir)

## combine all annot files
all.ls <- lapply(1:78021,function(i){
    message(i)
    readRDS(paste0(input_dir,i,'.rds')) %>% mutate(group = paste0(human_dist_group,mouse_dist_group))
})
if(!dir.exists(paste0(output_dir,'/comb/'))) dir.create(paste0(output_dir,'/comb/'))
saveRDS(all.ls, paste0(output_dir,'/comb/all_annotation_ls.rds'))

all.nrow <- sapply(all.ls, nrow)
saveRDS(all.nrow, paste0(output_dir,'/comb/all_annotation_nrow.rds'))
print(sum(all.nrow))

## split by 1e4 per file
all <- do.call(rbind,all.ls)
chunk <- 1e4
n <- nrow(all)
all$file  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
print(dim(all))
print(paste0('=== max: ', max(unique(all$file))))
saveRDS(all, paste0(output_dir,'/comb/all_annotation.rds'))

## save in annot_split
split.ls <- split(all, factor(all$file))
split.nrow <- sapply(split.ls,nrow)
idx <- which(split.nrow != chunk)
print('=== smaller chunks:')
print(idx)
w.ls <- lapply(1:max(unique(all$file)), function(i){
    saveRDS(split.ls[[i]],paste0(output_dir,i,'.rds'))
})
