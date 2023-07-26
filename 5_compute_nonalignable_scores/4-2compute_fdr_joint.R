## compute fdr (combine all categories and compute fdr)
library(data.table)
library(dplyr)
library(purrr)
library(readr)

setwd('/data/hji7/ENCODE4_unaligned/')
options(stringsAsFactors = F)

dataset <- 'chromatin' # can be replaced by other data modalities: H3K27ac, H3K4me1, H3K4me3
output_dir <- paste0('./result/',dataset,'/pvalue/')
new_dir <- paste0('./result/',dataset,'/fdr/')

# combine all 9 categories
dt.ls <- lapply(c('11','12','13','21','22','23','31','32','33'), function(category){
    print(category)
    comb.ls <- lapply(1:27,function(x){
        print(x)
        readRDS(paste0(output_dir,category,'_',x,'.rds'))
    })
    comb <- do.call(rbind,comb.ls)
})
dt <- do.call(rbind, dt.ls)

# compute FDR after pooling all categories
dt$weighted_spearman_fdr <- p.adjust(dt$weighted_spearman_p, method = "fdr")
dt$h_value_fdr <- p.adjust(dt$h_value_p, method = "fdr")

fdr_cutoff <- 0.25 # previous cutoff: 0.2
dt$weighted_spearman_sign <- 0
dt$weighted_spearman_sign[which(dt$weighted_spearman_fdr <= fdr_cutoff)] <- 1
dt$h_value_sign <- 0
dt$h_value_sign[which(dt$h_value_fdr <= fdr_cutoff)] <- 1
dt$weighted_spearman_cnsv <- paste0(dt$weighted_spearman_sign, dt$alignable)
dt$h_value_cnsv <- paste0(dt$h_value_sign, dt$alignable)

# combine with coordinates file
human_coord <- readRDS('./proc_data/coordinate_human.rds')
mouse_coord <- readRDS('./proc_data/coordinate_mouse.rds')

# extract significant
if(sum(dt$weighted_spearman_sign)>0){
    ws.sign <- dt[which(dt$weighted_spearman_sign == 1),]
    nrow(ws.sign)
    ws.sign2 <- suppressMessages(inner_join(ws.sign,human_coord) %>% inner_join(mouse_coord))
    if(nrow(ws.sign2) == nrow(ws.sign)){
        saveRDS(ws.sign2, paste0(new_dir, 'total_ws_sign.rds'))
        write_tsv(ws.sign2, paste0(new_dir, 'total_ws_sign.tsv'))
    }else{
        stop('ws nrow incomp')
    }
}

if(sum(dt$h_value_sign)>0){
    h.sign <- dt[which(dt$h_value_sign == 1),]
    nrow(h.sign)
    h.sign2 <- suppressMessages(inner_join(h.sign,human_coord) %>% inner_join(mouse_coord))
    if(nrow(h.sign2) == nrow(h.sign)){
        saveRDS(h.sign2, paste0(new_dir, 'total_h_sign.rds'))
        write_tsv(h.sign2, paste0(new_dir, 'total_h_sign.tsv'))
    }
}

# summarize by category
final.ls <- split(dt, factor(dt$group))
rm(dt)
l <- lapply(final.ls, function(x){
    group <- x$group[1]
    print('=====')
    print(group)

    summary <- data.frame(group = group,
                          total = nrow(x),
                          total_align = sum(x$alignable, na.rm = T),
                          total_nonalign = sum(x$alignable == 0, na.rm = T),

                          ws_sign = sum(x$weighted_spearman_sign, na.rm = T),
                          ws_sign_align = sum(x$weighted_spearman_cnsv == '11', na.rm = T),
                          ws_sign_nonalign = sum(x$weighted_spearman_cnsv == '10', na.rm = T),
                          ws_ns_align = sum(x$weighted_spearman_cnsv == '01', na.rm = T),
                          ws_ns_nonalign = sum(x$weighted_spearman_cnsv == '00', na.rm = T),

                          h_sign = sum(x$h_value_sign, na.rm = T),
                          h_sign_align = sum(x$h_value_cnsv == '11', na.rm = T),
                          h_sign_nonalign = sum(x$h_value_cnsv == '10', na.rm = T),
                          h_ns_align = sum(x$h_value_cnsv == '01', na.rm = T),
                          h_ns_nonalign = sum(x$h_value_cnsv == '00', na.rm = T))
    if(group == '11'){
        write.table(summary, paste0(new_dir,'summary.csv'), sep = ",",
                    append = T, quote = F, row.names = F)
    }else{
        write.table(summary, paste0(new_dir,'summary.csv'), sep = ",",
                    append = T, quote = F, row.names = F, col.names = F)
    }

    # add coordinates and save results
    if(group != '33'){
        new <- suppressMessages(inner_join(x,human_coord) %>% inner_join(mouse_coord))
        if(nrow(new) == nrow(x)){
            saveRDS(new, paste0(new_dir, group,'.rds'))
        }
    }else{
        # for group 33, further split into smaller chuncks
        n <- nrow(x)
        chunk <- 1e7
        start <- seq(1,n,by = chunk)
        end <- c(seq(chunk,n,by = chunk),n)

        for(i in 1:length(start)){
            print(paste0('==33 -',i))
            sub <- x[(start[i]:end[i]),]
            new <- suppressMessages(inner_join(sub,human_coord) %>% inner_join(mouse_coord))
            if(nrow(sub) == nrow(new)){
                print(nrow(new))
                saveRDS(new, paste0(new_dir, group,'_',i,'.rds'))
            }
        }

    }
})
