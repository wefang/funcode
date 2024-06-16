library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(patchwork)
source("./helper/helper.R")
# load mouse gene regulatory domains
mouse_genes_domain = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_protein_coding_gene_domains.rds")
mouse_genes_gr = readRDS("./metadata_processed/gene_annotations/ENCODE_Ref_gencode_vM21_protein_coding_genes.rds")

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
print(load("metadata_processed/indexed_aligned_combined.rda"))
mouse_tb = assign_gene(mm10_regions, mouse_genes_domain, mouse_genes_gr)
# load ABC interactions for human
print(load("./processed_data/abc/131_ABCscore_DHS_overlap_info.rda"))
human_tb = human.to.mouse %>%
        transmute(id = DHS.id, gene = TargetGene, dist_to_tss = distance) %>%
        group_by(id ,gene) %>%
        summarise(dist_to_tss = mean(dist_to_tss))

# homologous genes
source("./load_homolog.R")
common_indices = which((human_symbols %in% human_tb$gene) &
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

human_tb$hid = homolog_id[match(human_tb$gene, human_symbols)]
mouse_tb$hid = homolog_id[match(mouse_tb$gene, mouse_symbols)]
human_tb_nest = nest(human_tb, human_gene_tb = -id)
mouse_tb_nest = nest(mouse_tb, mouse_gene_tb = -id)

regions = readRDS("./output/indexed_dhs_mapped_regions_v1_manual.rds")
# need to remove older results
regions = regions[!names(regions) %in% c("human_gene_tb", "mouse_gene_tb")]
regions = left_join(regions, human_tb_nest)
regions = left_join(regions, mouse_tb_nest)

# mouse id need to be matched manually
regions$n_mouse = map_dbl(regions$mouse_gene_tb, function(x) {
        if (is.null(x)) {
                return(0)
        } else {
                return(nrow(x))
        }
})
regions$n_human = map_dbl(regions$human_gene_tb, function(x) {
        if (is.null(x)) {
                return(0)
        } else {
                return(nrow(x))
        }
})
# table(regions$n_human > 0, regions$n_mouse > 0)
regions$hid_count = map2_dbl(regions$human_gene_tb, regions$mouse_gene_tb, function(h_tb, m_tb) {
        if (!is.null(h_tb) & !is.null(m_tb)) {
                return(sum(h_tb$hid %in% m_tb$hid))
        } else {
                return(0)
        }
})
regions$consist_target = regions$hid_count > 0
g1 = regions %>% ggviolin(x = "consist_target", y = "cov_chromatin_accessibility",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-V-CA")
        # stat_compare_means()
g2 = regions %>% ggviolin(x = "consist_target", y = "cob_chromatin_accessibility",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-B-CA")
        # stat_compare_means()
(g1 + g2) %>% push_pdf(file_name = "abc_dhs_co_ca", width = 4., h = 3.5)

g1 = regions %>% ggviolin(x = "consist_target", y = "cov_H3K4me3",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-V-H3K4me3")
        # stat_compare_means()
g2 = regions %>% ggviolin(x = "consist_target", y = "cob_H3K4me3",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-B-H3K4me3")
        # stat_compare_means()
(g1 + g2) %>% push_pdf(file_name = "abc_dhs_co_H3K4me3", width = 4., h = 3.5)

g1 = regions %>% ggviolin(x = "consist_target", y = "cov_H3K27ac",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-V-H3K27ac")
# stat_compare_means()
g2 = regions %>% ggviolin(x = "consist_target", y = "cob_H3K27ac",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-B-H3K27ac")
# stat_compare_means()
(g1 + g2) %>% push_pdf(file_name = "abc_dhs_co_H3K27ac", width = 4., h = 3.5)

g1 = regions %>% ggviolin(x = "consist_target", y = "cov_H3K4me1",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-V-H3K4me1")
# stat_compare_means()
g2 = regions %>% ggviolin(x = "consist_target", y = "cob_H3K4me1",
                          size = 1, fill = "consist_target",
                          draw_quantiles = 0.5) +
        scale_fill_manual(values = c("#4453a4", "#f15a29")) +
        ylab("CO-B-H3K4me1")
# stat_compare_means()
(g1 + g2) %>% push_pdf(file_name = "abc_dhs_co_H3K4me1", width = 4., h = 3.5)


regions_gene = regions[which(regions$hid_count >=1 & regions$hid_count <= 2), ]
j = 0
regions_gene$homolog_gene_tb = map2(regions_gene$human_gene_tb, regions_gene$mouse_gene_tb, function(h_tb, m_tb) {
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
saveRDS(regions_gene, file = "./intermediate_data/regions_gene_v1.rds")


regions_gene = readRDS("./intermediate_data/regions_gene_v1.rds")
# Method 2: select overall closest gene
j = 0
regions_gene$homolog_gene_filter_tb = map(regions_gene$homolog_gene_tb, function(tb) {
        j <<- j + 1
        if (j %% 1e4 == 0) {
                message(j)
        }
        tb = mutate(tb,
                    dist_to_tss_human_rank = rank(dist_to_tss_human),
                    dist_to_tss_mouse_rank = rank(dist_to_tss_mouse),
                    dist_to_tss_rank = dist_to_tss_human_rank + dist_to_tss_mouse_rank)
        tb = filter(tb, dist_to_tss_rank == min(dist_to_tss_rank))
        tb  
})
regions_gene_ext = select(regions_gene, id,
                          human_cre_class_prio_v4,
                          mouse_cre_class_prio_v4,
                          cov_mean,
                          cob_mean,
                          lecif,
                          cov_chromatin_accessibility,
                          cov_H3K4me1,
                          cov_H3K4me3,
                          cov_H3K27ac,
                          cob_chromatin_accessibility,
                          cob_H3K4me1,
                          cob_H3K4me3,
                          cob_H3K27ac,
                          percentage, phastCons4way, phyloP4way,
                          cov_manual_chromatin_accessiblity,
                          cov_manual_H3K27ac,
                          cov_manual_H3K4me1,
                          cov_manual_H3K4me3,
                          homolog_gene_filter_tb) %>% 
        unnest(cols = "homolog_gene_filter_tb")

gene_tb = readRDS("./output/rna_conservation/rna_conservation_scores.rds")
regions_gene_ext = left_join(regions_gene_ext, gene_tb)
# regions_gene_ext = filter(regions_gene_ext, !is.na(gene_spearman) & !is.na(spearman_global))

# regions_ind = regions_gene_ext$cre_pair %in% pair_plot
produce_gene_cor <- function(score) {
        out = regions_gene_ext[regions_ind & !is.na(regions_gene_ext[, score]), ] %>%
                nest(data = -cre_pair) %>%
                transmute(cre_pair = cre_pair,
                          dhs_gene_cor = map_dbl(data, function(x) {
                                  indices = which(!is.na(x$gene_spearman))
                                  cor(x[indices, score], x$gene_spearman[indices], method = "spearman")
                          }),
                          Method = score)
        out
}

regions_gene_ext = left_join(regions_gene_ext,
                             select(regions, id, human_cre_class_prio_v4, mouse_cre_class_prio_v4, high_spearman, lecif, phyloP4way),
                             by = "id")

regions_gene_ext$is_human_cre = factor(regions_gene_ext$human_cre_class_prio_v4 != "non-cCRE", labels = c("non-cCRE", "cCRE"))
regions_gene_ext$is_mouse_cre = factor(regions_gene_ext$mouse_cre_class_prio_v4 != "non-cCRE", labels = c("non-cCRE", "cCRE"))

g_theme = theme(text = element_text(size = 8),
                legend.position = "none",
                strip.background = element_blank(),
                axis.text.x = element_text(angle = 45))
regions_gene_ext$high_spearman = factor(regions_gene_ext$high_spearman, labels = c("Other", "High CACO-V"))
g1 = regions_gene_ext %>%
        ggboxplot(x = "high_spearman",
                  y = "gene_spearman",
                  fill = "high_spearman",
                  facet.by = "is_human_cre",
                  size = 0.25,
                  outlier.size = 0.25) +
        stat_compare_means() +
        ylab("GECO-V") + xlab("") +
        g_theme
g2 = regions_gene_ext %>%
        ggboxplot(x = "high_spearman",
                  y = "gene_spearman",
                  fill = "high_spearman",
                  facet.by = "is_mouse_cre",
                  size = 0.25,
                  outlier.size = 0.25) +
        stat_compare_means() +
        ylab("") +
        xlab("") +
        g_theme
# (g1 + g2 + plot_layout(guide = "collect")) %>% push_pdf("geco-v_eval", w = 3.2, h = 3.)
# 
# # group_by(human_cre_class_prio_v4, high_spearman) %>%
# # summarise(mean = mean(gene_spearman, na.rm = T),
# #           se = sd(gene_spearman, na.rm = T)) %>%
# # ggbarplot(, y = "mean", position = position_dodge())

produce_gene_cor_cov <- function(score) {
        out = regions_gene_ext[!is.na(regions_gene_ext[, score]), ] %>%
                filter(mouse_cre_class_prio_v4 == human_cre_class_prio_v4) %>%
                nest(data = -mouse_cre_class_prio_v4) %>%
                transmute(mouse_cre_class_prio_v4 = mouse_cre_class_prio_v4,
                          dhs_gene_cor = map_dbl(data, function(x) {
                                  indices = which(!is.na(x$gene_spearman))
                                  cor(x[indices, score], x$gene_spearman[indices], method = "spearman")
                          }),
                          Method = score)
        out
}
produce_gene_cor_cob <- function(score) {
        out = regions_gene_ext[!is.na(regions_gene_ext[, score]), ] %>%
                filter(mouse_cre_class_prio_v4 == human_cre_class_prio_v4) %>%
                nest(data = -mouse_cre_class_prio_v4) %>%
                transmute(#human_cre_class_prio_v4 = human_cre_class_prio_v4,
                        mouse_cre_class_prio_v4 = mouse_cre_class_prio_v4,
                        dhs_gene_cor = map_dbl(data, function(x) {
                                indices = which(!is.na(x$gene_hval))
                                cor(x[indices, score], x$gene_hval[indices], method = "spearman")
                        }),
                        Method = score)
        out
}

method_col = c("cov_manual_chromatin_accessiblity" = "#0c2c84",
               "cov_manual_H3K27ac" = "#225ea8",
               "cov_manual_H3K4me3" = "#1d91c0",
               "cov_manual_H3K4me1" = "#41b6c4",
               "cov_chromatin_accessibility" = "#d73027",
               "cov_H3K27ac" = "#f46d43",
               "cov_H3K4me3" = "#fdae61",
               "cov_H3K4me1" = "#fee08b",
               "cob_mean" = "#006837",
               "cob_chromatin_accessibility" = "#1a9850",
               "cob_H3K4me3" = "#66bd63",
               "cob_H3K27ac" = "#a6d96a",
               "cob_H3K4me1" = "#d9ef8b",
               "percentage"  = "#737373",
               "phastCons4way" = "#bdbdbd",
               "phylop4way" = "#989898",
               "lecif" = "#d9d9d9")

g1 = bind_rows(list(
        produce_gene_cor_cob("cov_chromatin_accessibility"),
        produce_gene_cor_cob("cov_H3K4me1"),
        produce_gene_cor_cob("cov_H3K4me3"),
        produce_gene_cor_cob("cov_H3K27ac"),
        produce_gene_cor_cob("cob_chromatin_accessibility"),
        produce_gene_cor_cob("cob_H3K4me1"),
        produce_gene_cor_cob("cob_H3K4me3"),
        produce_gene_cor_cob("percentage"),
        produce_gene_cor_cob("phastCons4way"),
        produce_gene_cor_cob("phyloP4way"),
        produce_gene_cor_cob("lecif")
)) %>%
        filter(!mouse_cre_class_prio_v4 %in% c("CA", "TF", "CA-TF", "CA-CTCF", "CA-H3K4me3")) %>%
        # filter(mouse_cre_class_prio_v4 != "non-cCRE") %>%
        mutate(Method = factor(Method, levels = names(method_col))) %>%
        ggbarplot(x = "mouse_cre_class_prio_v4",
                  y = "dhs_gene_cor",
                  fill = "Method",
                  color = NA,
                  position = position_dodge()) +
        scale_fill_manual(values = method_col) + 
        theme(legend.position = "right",
              axis.text.x = element_text(angle = 30),
              text = element_text(size = 10),
              legend.key.size = unit(0.25, 'cm'))

g2 = bind_rows(list(#produce_gene_cor_cov("cov_mean"),
        produce_gene_cor_cov("cov_chromatin_accessibility"),
        produce_gene_cor_cov("cov_H3K4me1"),
        produce_gene_cor_cov("cov_H3K4me3"),
        produce_gene_cor_cov("cov_H3K27ac"),
        produce_gene_cor_cov("cob_chromatin_accessibility"),
        #produce_gene_cor_cov("cob_mean"),
        produce_gene_cor_cov("cob_H3K4me1"),
        produce_gene_cor_cov("cob_H3K4me3"),
        produce_gene_cor_cov("cob_H3K27ac"),
        produce_gene_cor_cov("percentage"),
        produce_gene_cor_cov("phastCons4way"),
        produce_gene_cor_cov("phyloP4way"),
        produce_gene_cor_cov("lecif")
        )) %>%
        filter(!mouse_cre_class_prio_v4 %in% c("CA", "TF", "CA-TF", "CA-CTCF", "CA-H3K4me3")) %>%
        mutate(Method = factor(Method, levels = names(method_col))) %>%
        # filter(mouse_cre_class_prio_v4 != "non-cCRE") %>%
        ggbarplot(x = "mouse_cre_class_prio_v4",
                  y = "dhs_gene_cor",
                  color = NA,
                  fill = "Method",
                  position = position_dodge()) +
        scale_fill_manual(values = method_col) +
        theme(legend.position = "right",
              axis.text.x = element_text(angle = 30),
              text = element_text(size = 10),
              legend.key.size = unit(0.25, 'cm'))

(g1 + g2 + plot_layout(ncol = 1, guides = "collect")) %>%
        push_pdf(file_name = "abc_domain_dhs_gene_cor_mouse2", width = 5.5, h = 5.)

produce_gene_cor_cov <- function(score) {
        out = regions_gene_ext[!is.na(regions_gene_ext[, score]), ] %>%
                # filter(mouse_cre_class_prio_v4 == human_cre_class_prio_v4) %>%
                nest(data = -human_cre_class_prio_v4) %>%
                transmute(human_cre_class_prio_v4 = human_cre_class_prio_v4,
                          dhs_gene_cor = map_dbl(data, function(x) {
                                  indices = which(!is.na(x$gene_spearman))
                                  cor(x[indices, score], x$gene_spearman[indices], method = "spearman")
                          }),
                          Method = score)
        out
}
produce_gene_cor_cob <- function(score) {
        out = regions_gene_ext[!is.na(regions_gene_ext[, score]), ] %>%
                # filter(mouse_cre_class_prio_v4 == human_cre_class_prio_v4) %>%
                nest(data = -human_cre_class_prio_v4) %>%
                transmute(#human_cre_class_prio_v4 = human_cre_class_prio_v4,
                        human_cre_class_prio_v4 = human_cre_class_prio_v4,
                        dhs_gene_cor = map_dbl(data, function(x) {
                                indices = which(!is.na(x$gene_hval))
                                cor(x[indices, score], x$gene_hval[indices], method = "spearman")
                        }),
                        Method = score)
        out
}

regions_gene_ext %>% count(mouse_cre_class_prio_v4)
regions_gene_ext %>% count(human_cre_class_prio_v4)

g1 = bind_rows(list(produce_gene_cor_cob("cov_chromatin_accessibility"),
        produce_gene_cor_cob("cov_H3K4me1"),
        produce_gene_cor_cob("cov_H3K4me3"),
        produce_gene_cor_cob("cov_H3K27ac"),
        produce_gene_cor_cob("cob_chromatin_accessibility"),
        produce_gene_cor_cob("cob_H3K4me1"),
        produce_gene_cor_cob("cob_H3K4me3"),
        produce_gene_cor_cob("cob_H3K27ac"),
        produce_gene_cor_cob("percentage"),
        produce_gene_cor_cob("phastCons4way"),
        produce_gene_cor_cob("phyloP4way"),
        produce_gene_cor_cob("lecif")
)) %>%
        filter(!human_cre_class_prio_v4 %in% c("CA", "TF", "CA-TF", "CA-CTCF", "CA-H3K4me3")) %>%
        # filter(human_cre_class_prio_v4 != "non-cCRE") %>%
        ggbarplot(x = "human_cre_class_prio_v4",
                  y = "dhs_gene_cor",
                  color = NA,
                  fill = "Method",
                  position = position_dodge()) +
        scale_fill_manual(values = method_col) + 
        theme(legend.position = "right")

g2 = bind_rows(list(#produce_gene_cor_cov("cov_mean"),
        produce_gene_cor_cov("cov_chromatin_accessibility"),
        produce_gene_cor_cov("cov_H3K4me1"),
        produce_gene_cor_cov("cov_H3K4me3"),
        produce_gene_cor_cov("cov_H3K27ac"),
        produce_gene_cor_cov("cob_chromatin_accessibility"),
        produce_gene_cor_cov("cob_H3K4me1"),
        produce_gene_cor_cov("cob_H3K4me3"),
        produce_gene_cor_cov("cob_H3K27ac"),
        produce_gene_cor_cov("percentage"),
        produce_gene_cor_cov("phastCons4way"),
        produce_gene_cor_cov("phyloP4way"),
        produce_gene_cor_cov("lecif")
)) %>%
        filter(!human_cre_class_prio_v4 %in% c("CA", "TF", "CA-TF", "CA-CTCF", "CA-H3K4me3")) %>%
        # filter(human_cre_class_prio_v4 != "non-cCRE") %>%
        ggbarplot(x = "human_cre_class_prio_v4",
                  y = "dhs_gene_cor",
                  color = NA,
                  fill = "Method",
                  position = position_dodge()) +
        scale_fill_manual(values = method_col) + 
        theme(legend.position = "right")
# end repeating grouped with human cCRE
(g1 + g2 + plot_layout(ncol = 1, guides = "collect")) %>%
        push_pdf(file_name = "abc_domain_dhs_gene_cor_human2", width = 6., h = 5.)


