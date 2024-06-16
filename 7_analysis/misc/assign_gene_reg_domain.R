# assign genes to DHS
library(tidyverse)
library(GenomicRanges)
# Input below for aligned regions
# load("./metadata_processed/indexed_aligned_regions.rda")
# human_regions = hg38_regions
# mouse_regions = mm10_regions
# Input below for all DHS
print(load("./metadata_processed/DHS/human_mouse_dhs_200bp.rda"))
# END input

# begin gene annotation
human_genes_domain = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_protein_coding_gene_domains.rds")
mouse_genes_domain = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_protein_coding_gene_domains.rds")
human_genes_gr = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_protein_coding_genes.rds")
mouse_genes_gr = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_protein_coding_genes.rds")
# end gene annotation

assign_gene <- function(regions, genes_domain, genes_gr) {
        gene_promoters = promoters(genes_gr, 0, 1)
        domain_ol = findOverlaps(regions, genes_domain)
        tb = data.frame(id = regions$identifier[queryHits(domain_ol)],
                        gene = as.character(genes_domain$gene_name[subjectHits(domain_ol)]),
                        dist_to_tss = distance(regions[queryHits(domain_ol)],
                                               gene_promoters[subjectHits(domain_ol)],
                                               ignore.strand = T),
                        stringsAsFactors = F)
        tb
}
human_tb = assign_gene(human_regions, human_genes_domain, human_genes_gr)
mouse_tb = assign_gene(mouse_regions, mouse_genes_domain, mouse_genes_gr)

hist(log2(human_tb$dist_to_tss))

# homologous genes
source("R_misc/load_homolog.R")
common_indices = which((human_symbols %in% human_genes_gr$gene_name) &
                               (mouse_symbols %in% mouse_genes_gr$gene_name))
homolog_id = homolog_id[common_indices]
human_symbols = human_symbols[common_indices]
mouse_symbols = mouse_symbols[common_indices]

# below if only unique gene names/id
# human_gene_name_count = table(human_genes_gr$gene_name)
# human_gene_name_unqiue = names(human_gene_name_count)[human_gene_name_count == 1]
# mouse_gene_name_count = table(mouse_genes_gr$gene_name)
# mouse_gene_name_unqiue = names(mouse_gene_name_count)[mouse_gene_name_count == 1]
# end homologous genes

#### This section specific to aligned DHS ####
# adding hid to genes
human_tb$hid = homolog_id[match(human_tb$gene, human_symbols)]
mouse_tb$hid = homolog_id[match(mouse_tb$gene, mouse_symbols)]
human_tb_nest = nest(human_tb, human_gene_tb = -id)
mouse_tb_nest = nest(mouse_tb, mouse_gene_tb = -id)

source("R_production/config.R")
regions = readRDS(region_file)
regions = left_join(regions, human_tb_nest)
regions = left_join(regions, mouse_tb_nest)

regions$hid_count = map2_dbl(regions$human_gene_tb, regions$mouse_gene_tb, function(h_tb, m_tb) {
        if (!is.null(h_tb) & !is.null(m_tb)) {
                return(sum(h_tb$hid %in% m_tb$hid))
        } else {
                return(0)
        }
})
j = 0
regions$homolog_gene_tb = map2(regions$human_gene_tb, regions$mouse_gene_tb, function(h_tb, m_tb) {
        j <<- j + 1
        if (j %% 10000 == 0) {
                message(j)
        }
        if (!is.null(h_tb) & !is.null(m_tb)) {
                if (sum(h_tb$hid %in% m_tb$hid) > 0) {
                        h_tb = filter(h_tb, !is.na(hid))
                        b_tb = filter(m_tb, !is.na(hid))
                        hid_tb = inner_join(h_tb, m_tb, by = "hid", suffix = c("_human", "_mouse"))
                        return(hid_tb)
                } else {
                        return(NULL)
                }
        } else {
                return(NULL)
        }
})
saveRDS(regions, file = region_file)
#### end section aligned DHS ####

#### SECTION below specific to ALL (unaligned) DHS ####
# removing some duplicated gene annotations
human_tb_sel = dplyr::filter(human_tb, gene %in% human_symbols)
human_tb_sel = nest(human_tb_sel, data = -c(id, hid))
human_tb_sel$data = map(human_tb_sel$data, function(x) {
        if(nrow(x) > 1) {
                x[which.min(x$dist_to_tss), ][1, ]
        } else {
                return(x)
        }
})
human_tb_sel = human_tb_sel %>% unnest(cols = data)

mouse_tb_sel = dplyr::filter(mouse_tb, gene %in% mouse_symbols)
mouse_tb_sel = nest(mouse_tb_sel, data = -c(id, hid))
mouse_tb_sel$data = map(mouse_tb_sel$data, function(x) {
        if(nrow(x) > 1) {
                x[which.min(x$dist_to_tss), ][1, ]
        } else {
                return(x)
        }
})
mouse_tb_sel = mouse_tb_sel %>% unnest(cols = data)
human_tb_sel$hid = homolog_id[match(human_tb_sel$gene, human_symbols)]
mouse_tb_sel$hid = homolog_id[match(mouse_tb_sel$gene, mouse_symbols)]

# put each dhs into a category of hid set
human_hid_set = human_tb_sel %>% nest(human_gene = -id) %>%
        mutate(hid_set = map_chr(human_gene, function(x) paste0(sort(x$hid), collapse = "_"))) %>%
        nest(human_dhs = -hid_set)
mouse_hid_set = mouse_tb_sel %>% nest(mouse_gene = -id) %>%
        mutate(hid_set = map_chr(mouse_gene, function(x) paste0(sort(x$hid), collapse = "_"))) %>%
        nest(mouse_dhs = -hid_set)

# first each human and mouse dhs classified into hid sets
# generate all hid set combinations
# criteria for hid set pair: overlap of hids
human_hid_set$hid_list = strsplit(human_hid_set$hid_set, "_")
mouse_hid_set$hid_list = strsplit(mouse_hid_set$hid_set, "_")

library(furrr)
plan(multisession, workers = 12)
hidset_pairs = bind_rows(future_map(1:nrow(human_hid_set), function(i) {
        message(i)
        x = human_hid_set$hid_list[[i]]
        tibble(human_hid_set = human_hid_set$hid_set[i],
               mouse_hid_set = mouse_hid_set$hid_set[map_lgl(mouse_hid_set$hid_list, function(y) all(y %in% x) | all(x %in% y))])
}, .progress = T))
save(human_hid_set, mouse_hid_set, hidset_pairs, file = "./output/all_dhs_pairs.rda")

load("./output/all_dhs_pairs.rda")
# human_hid_set$human_dhs_close = map(human_hid_set$human_dhs, function(x) {
#         x$human_gene = map(x$human_gene, function(y) {
#                 y = y[which.min(y$dist_to_tss), ][1, ]
#                 unnest(y, cols = human_gene)
#         })
# })
run_hidset_pair <- function(i) {
        human_ind = which(human_hid_set$hid_set == hidset_pairs$human_hid_set[i])
        mouse_ind = which(mouse_hid_set$hid_set == hidset_pairs$mouse_hid_set[i])
        gene_int = intersect(human_hid_set$hid_list[[human_ind]],
                             mouse_hid_set$hid_list[[mouse_ind]])
        human_dhs_tb = human_hid_set$human_dhs[[human_ind]]
        mouse_dhs_tb = mouse_hid_set$mouse_dhs[[mouse_ind]]
        
        dhs_pairs = as_tibble(expand.grid(human_id = human_dhs_tb$id,
                              mouse_id = mouse_dhs_tb$id)) %>%
                left_join(human_dhs_tb, by = c("human_id" = "id")) %>%
                left_join(mouse_dhs_tb, by = c("mouse_id" = "id")) %>%
                mutate(gene_tb = map2(human_gene, mouse_gene, function(x, y) {
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
                unnest(cols = "gene_tb")
        return(dhs_pairs)
}



matched_dhs = full_join(nest(human_tb_sel, human_dhs = -hid),
                        nest(mouse_tb_sel, mouse_dhs = -hid))
matched_dhs$n_human = map_dbl(matched_dhs$human_dhs, function(x) {
        if (!is.null(x)) {
                nrow(x)
        } else {
                return(NA)
        }
})
matched_dhs$n_mouse = map_dbl(matched_dhs$mouse_dhs, function(x) {
        if (!is.null(x)) {
                nrow(x)
        } else {
                return(NA)
        }
})
# filter(matched_dhs, n_mouse < 3000) %>% ggscatter(x = "n_human", y = "n_mouse") + geom_smooth()
median(matched_dhs$n_human, na.rm = T)
median(matched_dhs$n_mouse, na.rm = T)
sum(matched_dhs$n_human * matched_dhs$n_mouse, na.rm = T)/1e8
# human_regions_tb = as_tibble(human_regions) %>% select(c(1, 2, 3, 6)) %>% rename(hg38_seqnames = seqnames, hg38_start = start, hg38_end = end)
# mouse_regions_tb = as_tibble(mouse_regions) %>% select(c(1, 2, 3, 6)) %>% rename(mm10_seqnames = seqnames, mm10_start = start, mm10_end = end)

# for each human DHS with multiple genes:
# list all mouse DHS that spans the same set of genes
# Alternatively, index by the unique hid combinations
mouse_id_gene_dup$hid_set = map_chr(mouse_id_gene_dup$mouse_gene, function(x) paste0(sort(x$hid), collapse = "_"))
length(unique(mouse_id_gene_dup$hid_set))

# all (non-aligned DHS)
matched_dhs$human_symbols = human_symbols[matched_dhs$hid]
matched_dhs$mouse_symbols = mouse_symbols[matched_dhs$hid]
saveRDS(matched_dhs, file = "./metadata_processed/human_mouse_all_DHS_matched_protein_coding.rds")
#### END section ALL DHS ####

#### Section repeat for the list of human and mouse mapped DHS ####
library(tidyverse)
library(GenomicRanges)
# Input below for aligned regions
print(load("metadata_processed/indexed_aligned_combined.rda"))
# END input

# begin gene annotation
human_genes_domain = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_protein_coding_gene_domains.rds")
mouse_genes_domain = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_protein_coding_gene_domains.rds")
human_genes_gr = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_v29_protein_coding_genes.rds")
mouse_genes_gr = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_protein_coding_genes.rds")
# end gene annotation

assign_gene <- function(regions, genes_domain, genes_gr) {
        gene_promoters = promoters(genes_gr, 0, 1)
        domain_ol = findOverlaps(regions, genes_domain)
        tb = data.frame(id = regions$identifier[queryHits(domain_ol)],
                        gene = as.character(genes_domain$gene_name[subjectHits(domain_ol)]),
                        dist_to_tss = distance(regions[queryHits(domain_ol)],
                                               gene_promoters[subjectHits(domain_ol)],
                                               ignore.strand = T),
                        stringsAsFactors = F)
        tb
}
human_tb = assign_gene(hg38_regions, human_genes_domain, human_genes_gr)
mouse_tb = assign_gene(mm10_regions, mouse_genes_domain, mouse_genes_gr)

# homologous genes
source("7_analysis/2_gene_expression_conservation/load_homolog.R")
common_indices = which((human_symbols %in% human_genes_gr$gene_name) &
                               (mouse_symbols %in% mouse_genes_gr$gene_name))
homolog_id = homolog_id[common_indices]
human_symbols = human_symbols[common_indices]
mouse_symbols = mouse_symbols[common_indices]

# below if only unique gene names/id
# human_gene_name_count = table(human_genes_gr$gene_name)
# human_gene_name_unqiue = names(human_gene_name_count)[human_gene_name_count == 1]
# mouse_gene_name_count = table(mouse_genes_gr$gene_name)
# mouse_gene_name_unqiue = names(mouse_gene_name_count)[mouse_gene_name_count == 1]
# end homologous genes

#### This section specific to aligned DHS ####
# adding hid to genes
human_tb$hid = homolog_id[match(human_tb$gene, human_symbols)]
human_tb_nest = nest(human_tb, human_gene_tb = -id)
mouse_tb$hid = homolog_id[match(mouse_tb$gene, mouse_symbols)]
mouse_tb_nest = nest(mouse_tb, mouse_gene_tb = -id)

regions = readRDS("./output/indexed_dhs_mapped_regions_combined.rds")
regions$id = hg38_regions$identifier
regions = left_join(regions, human_tb_nest)
regions = left_join(regions, mouse_tb_nest)

regions$hid_count = map2_dbl(regions$human_gene_tb, regions$mouse_gene_tb, function(h_tb, m_tb) {
        if (!is.null(h_tb) & !is.null(m_tb)) {
                return(sum(h_tb$hid %in% m_tb$hid))
        } else {
                return(0)
        }
})
j = 0
regions$homolog_gene_tb = map2(regions$human_gene_tb, regions$mouse_gene_tb, function(h_tb, m_tb) {
        j <<- j + 1
        if (j %% 100000 == 0) {
                message(j)
        }
        if (!is.null(h_tb) & !is.null(m_tb)) {
                if (sum(h_tb$hid %in% m_tb$hid) > 0) {
                        h_tb = filter(h_tb, !is.na(hid))
                        b_tb = filter(m_tb, !is.na(hid))
                        hid_tb = inner_join(h_tb, m_tb, by = "hid", suffix = c("_human", "_mouse"))
                        return(hid_tb)
                } else {
                        return(NULL)
                }
        } else {
                return(NULL)
        }
})
regions_gene = regions[which(regions$hid_count > 0), ]
saveRDS(regions_gene, file = "./intermediate_data/regions_combined_gene.rds")
#### End Section ####
