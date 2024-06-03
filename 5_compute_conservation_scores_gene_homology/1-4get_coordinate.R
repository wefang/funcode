# prepare coordinates
library(dplyr)
setwd('/data/hji7/ENCODE4_unaligned/')

load('./raw_data/ENCODE_human_mouse_DHS_noY.rda')
human <- human_regions %>% data.frame %>% select(-c(strand)) %>% as_tibble()
colnames(human)
colnames(human) <- c('human_chromosome','human_start','human_end','human_width','human_id')
saveRDS(human, './proc_data/coordinate_human.rds')


mouse <- mouse_regions %>% data.frame %>% select(-c(strand)) %>% as_tibble()
colnames(mouse)
colnames(mouse) <- c('mouse_chromosome','mouse_start','mouse_end','mouse_width','mouse_id')
saveRDS(mouse, './proc_data/coordinate_mouse.rds')


