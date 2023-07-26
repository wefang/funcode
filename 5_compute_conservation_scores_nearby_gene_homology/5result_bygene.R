# split results by gene
library(data.table)
library(dplyr)
library(purrr)
library(readr)

setwd('/data/hji7/ENCODE4_unaligned/')
options(stringsAsFactors = F)

dataset <- 'chromatin' # can be replaced by other data modalities: H3K27ac, H3K4me1, H3K4me3
new_dir <- paste0('./result/',dataset,'/fdr/')
output_dir <- paste0('./result/',dataset,'/by_gene/')
if(!dir.exists(output_dir)) dir.create(output_dir)

files <- c('11', '22', '12', '13', '23', '21', '31', '32',
           paste0('33_',1:125))

dt.ls <- list()
sum_nrow <- 0
for(i in 1:length(files)){
    f <- files[i]
    print(f)
    dt.ls[[i]] <- readRDS(paste0(new_dir,f,'.rds'))
    # count row
    sum_nrow <- sum_nrow + nrow(dt.ls[[i]])
}
dt <- do.call(rbind,dt.ls)
rm(dt.ls)
print(sum_nrow)
print(nrow(dt))
identical(as.integer(sum_nrow),as.integer(nrow(dt)))
identical(1318719814,as.integer(sum_nrow))

# dt <- dt[,c(1,2,4,6,20:23,
#             3,5,7,24:27,
#             8,9,11,14,15,13,16)]
dt$split <- paste0(dt$hid,'_',dt$human_symbol,'_',dt$mouse_symbol)

comb.ls <- split(dt,factor(dt$split))
saveRDS(comb.ls, paste0(new_dir,'comb_ws_bygene_ls.rds'))
max_col <- ncol(comb.ls[[1]])
print(length(comb.ls))
rm(dt)

new_nrow <- 0
for(i in 1:length(comb.ls)){
    print(i)
    final <- comb.ls[[i]]
    data.table::fwrite(final[,-max_col], paste0(output_dir,names(comb.ls)[i],'.tsv'))
    
    summary <- data.frame(file = names(comb.ls)[i],
                          n_human = length(unique(final$human_id)),
                          n_mouse = length(unique(final$mouse_id)),
                          n_total = nrow(final))
    new_nrow <- new_nrow + summary$n_total
    
    if(i == 1){
        write.table(summary, paste0(output_dir,'summary.csv'), sep = ",",
                    append = T, quote = F, row.names = F)
    }else{
        write.table(summary, paste0(output_dir,'summary.csv'), sep = ",",
                    append = T, quote = F, row.names = F, col.names = F)
    }
    rm(final)
}
print(new_nrow)
identical(sum_nrow,new_nrow)
