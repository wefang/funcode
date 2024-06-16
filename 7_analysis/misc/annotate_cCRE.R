# this script wraps annotation of cCRE for genomic regions, based on current ENCODE cCRE files
library(tidyverse)
library(GenomicRanges)
# input DHS
# INPUT cCRE file
meta_dir = "./metadata/ENCODE_cCRE_annotation/"
# annotation file metadata
cre_labels = c("dELS", "DNase-H3K4me3", "pELS", "CTCF-only", "PLS")
# collapse 'CTCF bound' labels
cre_label_mapper = c(cre_labels, cre_labels)
names(cre_label_mapper) = c(cre_labels, paste0(cre_labels, ",CTCF-bound"))
# assign priority to different cCRE labels
cre_label_priority = 1:5
names(cre_label_priority) = c(cre_labels[c(4, 2, 1, 3, 5)])

assign_human_cre <- function(human_gr) {
        hg38_el_files = readr::read_tsv(paste0(meta_dir, "/hg38_metadata.tsv"))
        human_dhs_cre_map = list()
        human_dhs_cre_class_map = list()
        human_cre_df = list()
        for (j in 1:5) { # loop over cre labels
                # read annotation and make ranges
                x_df = readr::read_tsv(paste0(meta_dir, hg38_el_files$`File accession`, ".bed")[j], col_names = F)
                x_gr = x_df %>% transmute(seqnames = X1, start = X2, end = X3) %>% makeGRangesFromDataFrame()
                
                # finding overlaps
                x_ol = findOverlaps(human_gr, x_gr)
                # cre maps to CRE ID
                # cre class maps to type of CRE
                cre_map = split(x_df$X4[subjectHits(x_ol)], human_gr$identifier[queryHits(x_ol)])
                cre_class_map = map(split(x_df$X10[subjectHits(x_ol)], human_gr$identifier[queryHits(x_ol)]), unique)
                human_dhs_cre_map[[cre_labels[j]]] = cre_map
                human_dhs_cre_class_map[[cre_labels[j]]] = cre_class_map
                
                dhs_map = split(human_gr$identifier[queryHits(x_ol)], x_df$X4[subjectHits(x_ol)])
                dhs_map_df = tibble(X4 = names(dhs_map), dhs = dhs_map)
                x_df = left_join(x_df, dhs_map_df)
                human_cre_df[[j]] = x_df
        }
        human_dhs_cre_df = tibble(id = sort(unique(do.call(c, map(human_dhs_cre_map, names)))))
        human_dhs_cre_df$human_cre = map(1:nrow(human_dhs_cre_df), function(i) character())
        human_dhs_cre_df$human_cre_class = map(1:nrow(human_dhs_cre_df), function(i) character())
        
        # make columns for each category
        for (j in 1:5) {
                temp_tb = tibble(id = names(human_dhs_cre_map[[j]]))
                temp_tb = add_column(temp_tb, !!(cre_labels[j]) := human_dhs_cre_map[[j]])
                human_dhs_cre_df = left_join(human_dhs_cre_df, temp_tb)
                human_dhs_cre_df = mutate(human_dhs_cre_df, human_cre = map2(human_cre, !! rlang::sym(cre_labels[j]),
                                                                             function(a, b) {
                                                                                     c(a, b)
                                                                             }))
                temp_tb1 = tibble(id = names(human_dhs_cre_class_map[[j]]))
                temp_tb1 = add_column(temp_tb1, !!(paste0(cre_labels[j], "_class")) := human_dhs_cre_class_map[[j]])
                human_dhs_cre_df = left_join(human_dhs_cre_df, temp_tb1)
                human_dhs_cre_df = mutate(human_dhs_cre_df, human_cre_class = map2(human_cre_class, !! rlang::sym(paste0(cre_labels[j], "_class")),
                                                                                   function(a, b) {
                                                                                           c(a, b)
                                                                                   }))
        }
        human_dhs_cre_df = mutate(human_dhs_cre_df, human_cre_class_collapse = map(human_cre_class, function(x) {
                cre_label_mapper[x]
        }))
        human_dhs_cre_df = human_dhs_cre_df %>% mutate(human_cre_class_prio = map_chr(human_cre_class_collapse, function(x) {
                if (length(x) == 0) {
                        return("non-cCRE")
                }
                names(cre_label_priority)[max(cre_label_priority[x])]
        }))
        human_dhs_cre_df$human_cre_class_prio = factor(human_dhs_cre_df$human_cre_class_prio, levels = c("non-cCRE", names(cre_label_priority)))
        human_dhs_cre_df
}

assign_mouse_cre <- function(mouse_gr) {
        mm10_el_files = readr::read_tsv(paste0(meta_dir, "/mm10_metadata.tsv"))
        mouse_dhs_cre_map = list()
        mouse_dhs_cre_class_map = list()
        mouse_cre_df = list()
        for (j in 1:5) {
                x_df = readr::read_tsv(paste0(meta_dir, mm10_el_files$`File accession`, ".bed")[j], col_names = F)
                x_gr = x_df %>% transmute(seqnames = X1, start = X2, end = X3) %>% makeGRangesFromDataFrame()
                x_ol = findOverlaps(mouse_gr, x_gr)
                
                cre_map = split(x_df$X4[subjectHits(x_ol)], mouse_gr$identifier[queryHits(x_ol)])
                cre_class_map = map(split(x_df$X10[subjectHits(x_ol)], mouse_gr$identifier[queryHits(x_ol)]), unique)
                mouse_dhs_cre_map[[cre_labels[j]]] = cre_map
                mouse_dhs_cre_class_map[[cre_labels[j]]] = cre_class_map
                
                dhs_map = split(mouse_gr$identifier[queryHits(x_ol)], x_df$X4[subjectHits(x_ol)])
                dhs_map_df = tibble(X4 = names(dhs_map), dhs = dhs_map)
                x_df = left_join(x_df, dhs_map_df)
                mouse_cre_df[[j]] = x_df
        }
        mouse_dhs_cre_df = tibble(id = sort(unique(do.call(c, map(mouse_dhs_cre_map, names)))))
        mouse_dhs_cre_df$mouse_cre = lapply(1:nrow(mouse_dhs_cre_df), function(i) character())
        mouse_dhs_cre_df$mouse_cre_class = lapply(1:nrow(mouse_dhs_cre_df), function(i) character())
        for (j in 1:5) {
                temp_tb = tibble(id = names(mouse_dhs_cre_map[[j]]))
                temp_tb = add_column(temp_tb, !!(cre_labels[j]) := mouse_dhs_cre_map[[j]])
                mouse_dhs_cre_df = left_join(mouse_dhs_cre_df, temp_tb)
                mouse_dhs_cre_df = mutate(mouse_dhs_cre_df, mouse_cre = map2(mouse_cre, !! rlang::sym(cre_labels[j]),
                                                                             function(a, b) {
                                                                                     c(a, b)
                                                                             }))
                temp_tb1 = tibble(id = names(mouse_dhs_cre_class_map[[j]]))
                temp_tb1 = add_column(temp_tb1, !!(paste0(cre_labels[j], "_class")) := mouse_dhs_cre_class_map[[j]])
                mouse_dhs_cre_df = left_join(mouse_dhs_cre_df, temp_tb1)
                mouse_dhs_cre_df = mutate(mouse_dhs_cre_df, mouse_cre_class = map2(mouse_cre_class, !! rlang::sym(paste0(cre_labels[j], "_class")),
                                                                                   function(a, b) {
                                                                                           c(a, b)
                                                                                   }))
        }
        mouse_dhs_cre_df = mutate(mouse_dhs_cre_df, mouse_cre_class_collapse = map(mouse_cre_class, function(x) {
                cre_label_mapper[x]
        }))
        mouse_dhs_cre_df = mouse_dhs_cre_df %>% mutate(mouse_cre_class_prio = map_chr(mouse_cre_class_collapse, function(x) {
                if (length(x) == 0) {
                        return("non-cCRE")
                }
                names(cre_label_priority)[max(cre_label_priority[x])]
        }))
        mouse_dhs_cre_df$mouse_cre_class_prio = factor(mouse_dhs_cre_df$mouse_cre_class_prio, levels = c("non-cCRE", names(cre_label_priority)))
        mouse_dhs_cre_df
}

# All DHS
load("./metadata_processed/DHS/human_mouse_dhs_raw.rda")
human_cre_tb = assign_human_cre(human_regions)
mouse_cre_tb = assign_mouse_cre(mouse_regions)
save(human_cre_tb, mouse_cre_tb, file = "./metadata_processed/DHS_cCRE/dhs_raw_cCRE.rda")

load("./metadata_processed/DHS/human_mouse_dhs_200bp.rda")
human_cre_tb = assign_human_cre(human_regions)
mouse_cre_tb = assign_mouse_cre(mouse_regions)
save(human_cre_tb, mouse_cre_tb, file = "./metadata_processed/DHS_cCRE/dhs_200bp_cCRE.rda")

load("./processed_data/indexed_dhs_hg38_mm10_mapped/regions/hg38_mm10_regions.rda")
hg38_regions$identifier = regions$id
mm10_regions$identifier = regions$id
human_cre_tb = assign_human_cre(hg38_regions)
mouse_cre_tb = assign_mouse_cre(mm10_regions)
save(human_cre_tb, mouse_cre_tb, file = "./metadata_processed/DHS_cCRE/aligned_cCRE.rda")

# updating region file
region_file = "metadata_processed/indexed_aligned_combined.rda"
regions = readRDS(region_file)
regions = left_join(regions, dplyr::select(human_cre_tb, id, human_cre_class_prio))
regions$human_cre_class_prio[is.na(regions$human_cre_class_prio)] = "non-cCRE"
regions = left_join(regions, dplyr::select(mouse_cre_tb, id, mouse_cre_class_prio))
regions$mouse_cre_class_prio[is.na(regions$mouse_cre_class_prio)] = "non-cCRE"
saveRDS(regions, file = region_file)



