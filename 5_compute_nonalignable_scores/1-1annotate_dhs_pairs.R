## compile and annotate DHS pairs
library(dplyr)
library(tidyr)
library(purrr)
setwd('/data/hji7/ENCODE4_unaligned/')

output_dir <- './result/annot/'
if(!dir.exists(output_dir)) dir.create(output_dir)

# === read in
# (1) DHS pairs and related info
load('./raw_data/all_dhs_pairs.rda')
hidset_pairs # human_hid_set, mouse_hid_set (78021 pairs)

## summary
human_hid_set$human_nrow <- sapply(1:nrow(human_hid_set),function(x) nrow(human_hid_set$human_dhs[x][[1]]))
mouse_hid_set$mouse_nrow <- sapply(1:nrow(mouse_hid_set),function(x) nrow(mouse_hid_set$mouse_dhs[x][[1]]))
all.summary <- inner_join(hidset_pairs,human_hid_set %>% select(hid_set,human_nrow) %>% rename(human_hid_set = hid_set)) %>%
    inner_join(mouse_hid_set %>% select(hid_set,mouse_nrow) %>% rename(mouse_hid_set = hid_set)) %>%
    mutate(total_nrow = human_nrow * mouse_nrow)
saveRDS(all.summary,'./proc_data/all_dhs_pairs_summary.rds')


# (2) GC content
human.gc <- readRDS('./raw_data/DHS_annotations/indexed_All_DHS_200bp_GRCh38_gc_content.rds') %>%
    mutate(gc_group = cut(gc_content,breaks = quantile(gc_content), include.lowest = T, labels = 1:4) %>% as.character %>% as.integer)
colnames(human.gc) <- paste0('human_',colnames(human.gc))
mouse.gc <- readRDS('./raw_data/DHS_annotations/indexed_All_DHS_200bp_mm10_gc_content.rds') %>%
    mutate(gc_group = cut(gc_content,breaks = quantile(gc_content), include.lowest = T, labels = 1:4) %>% as.character %>% as.integer)
colnames(mouse.gc) <- paste0('mouse_',colnames(mouse.gc))

# === process data (individual human-mouse pairs per hid)
run_hidset_pair <- function(i) {
    message(i)
    human_ind = which(human_hid_set$hid_set == hidset_pairs$human_hid_set[i])
    mouse_ind = which(mouse_hid_set$hid_set == hidset_pairs$mouse_hid_set[i])
    gene_int = intersect(human_hid_set$hid_list[[human_ind]],
                         mouse_hid_set$hid_list[[mouse_ind]])
    human_dhs_tb = human_hid_set$human_dhs[[human_ind]]
    mouse_dhs_tb = mouse_hid_set$mouse_dhs[[mouse_ind]]

    dhs_pairs = suppressMessages(as_tibble(expand.grid(human_id = human_dhs_tb$id,
                                                       mouse_id = mouse_dhs_tb$id)) %>%
                                     left_join(human_dhs_tb, by = c("human_id" = "id")) %>%
                                     left_join(mouse_dhs_tb, by = c("mouse_id" = "id")) %>%
                                     mutate(gene_tb = purrr::map2(human_gene, mouse_gene, function(x, y) {
                                         out = inner_join(dplyr::rename(x,
                                                                        human_symbol = gene,
                                                                        human_dist = dist_to_tss),
                                                          dplyr::rename(y,
                                                                        mouse_symbol = gene,
                                                                        mouse_dist = dist_to_tss),
                                                          by = "hid")
                                         out = out[which.min(out$human_dist + out$mouse_dist), ][1, ]
                                         out
                                     })) %>% dplyr::select(human_id, mouse_id, gene_tb) %>%
                                     tidyr::unnest(cols = "gene_tb") %>%
                                     mutate(human_dist_group = ifelse(human_dist <= 200,1,ifelse(human_dist <= 2000,2,3)),
                                            mouse_dist_group = ifelse(mouse_dist <= 200,1,ifelse(mouse_dist <= 2000,2,3))) %>%
                                     left_join(human.gc) %>%
                                     left_join(mouse.gc) %>%
                                     mutate(group = paste0(human_dist_group,human_gc_group,
                                                           mouse_dist_group,mouse_gc_group)))
    return(dhs_pairs)
}

# 78021 hidset_pairs
dat.ls <- mclapply(1:nrow(hidset_pairs),function(i){
    if(!file.exists(paste0(output_dir,i,'.rds'))){
        res <- run_hidset_pair(i)
        saveRDS(res,paste0(output_dir,i,'.rds'))
    }
},mc.cores = 20)

# compile human_dhs_info; mouse_dhs_info
new <- all.summary %>% mutate(ID = 1:nrow(all.summary))
new.human <- new %>% select(ID,human_hid_set,human_nrow)
human.dat <- new.human[!duplicated(new.human[,-1]),] # 26496
sum(human.dat$human_nrow) # 3395914

new.mouse <- new %>% select(ID,mouse_hid_set,mouse_nrow)
mouse.dat <- new.mouse[!duplicated(new.mouse[,-1]),] # 26478
sum(mouse.dat$mouse_nrow) # 1692289

human.ls <- lapply(human.dat$ID,function(i){
    dt <- readRDS(paste0(output_dir,i,'.rds')) %>%
        select(human_id, hid, human_symbol,
               human_dist, human_dist_group,
               human_gc_content, human_gc_group) %>%
        mutate(human_group = paste0(human_dist_group,human_gc_group)) %>%
        unique
    dt
})

comb_dir <- paste0(output_dir,'comb/')
if(!dir.exists(comb_dir)) dir.create(comb_dir)
human.info <- do.call(rbind,human.ls)
saveRDS(human.info, paste0(comb_dir,'human_dhs_info.rds'))

mouse.ls <- lapply(mouse.dat$ID,function(i){
    dt <- readRDS(paste0(output_dir,i,'.rds')) %>%
        select(mouse_id, hid, mouse_symbol,
               mouse_dist, mouse_dist_group,
               mouse_gc_content, mouse_gc_group) %>%
        mutate(mouse_group = paste0(mouse_dist_group,mouse_gc_group)) %>%
        unique
    dt
})

mouse.info <- do.call(rbind,mouse.ls)
saveRDS(mouse.info, paste0(comb_dir,'mouse_dhs_info.rds'))

