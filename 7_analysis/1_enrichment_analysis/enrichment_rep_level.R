library(tidyverse)
library(ggpubr)
library(patchwork)
source("./helper/helper.R")

# color and naming definitions
all_assays = c("dnase_atac", "H3K27ac", "H3K4me1","H3K4me3")
assay_name_mapping = c("dnase_atac" = "chromatin_accessibility",
                       "H3K27ac" ="H3K27ac",
                       "H3K4me1" = "H3K4me1",
                       "H3K4me3" = "H3K4me3")
co_col = c("High CO-V" = "#9B1C4F",
           "High CO-B" = "#FFD56A",
           "Other" = "#868686")

# Comibine forward and reverse mapped regions for core set
regions =readRDS("./output/indexed_dhs_mapped_regions_v1.rds")
regions_mm10 =readRDS("./output/indexed_mouse_dhs_mapped_regions.rds")
print(load("./intermediate_data/mouse_DHS_alignment_overlap.rda"))
indices_add = which(filtering_df$sum[match(regions_mm10$id, paste0("MouseDHS-", filtering_df$identifier))] < 2)
regions_mm10_fil = regions_mm10[indices_add, ]
common_names = intersect(names(regions), names(regions_mm10_fil))
regions_all = bind_rows(dplyr::rename(regions[c(common_names, "mouse_id")], human_id = id),
                        dplyr::rename(regions_mm10_fil[c(common_names, "human_id")], mouse_id = id))
saveRDS(regions_all, file = "./output/indexed_dhs_mapped_regions_combined.rds")

# make corresponding GRanges object
print(load("./metadata_processed//mm10_hg38_regions_cleaned.rda"))
hg38_regions_add = hg38_regions[indices_add]
mm10_regions_add = mm10_regions[indices_add]
print(load("metadata_processed/indexed_aligned_regions.rda"))
hg38_regions = c(hg38_regions, hg38_regions_add)
mm10_regions = c(mm10_regions, mm10_regions_add)
save(hg38_regions, mm10_regions, file = "metadata_processed/indexed_aligned_combined.rda")
# combined region
regions_all = readRDS("./output/indexed_dhs_mapped_regions_combined.rds")
regions_all$cov_any3 = (regions_all$cov_H3K27ac_sign + regions_all$cov_H3K4me3_sign + regions_all$cov_chromatin_accessibility_sign) > 0


# Figure S6: UpSet plot for combinations of conservation modalities
library(UpSetR)
cov_sign_list = list(cov_ca = which(regions_all$cov_chromatin_accessibility_sign == 1),
                     cov_H3K4me3 = which(regions_all$cov_H3K4me3_sign == 1),
                     cov_H3K4me1 = which(regions_all$cov_H3K4me1_sign == 1),
                     cov_H3K27ac = which(regions_all$cov_H3K27ac_sign == 1))
cob_sign_list = list(cob_ca = which(regions_all$cob_chromatin_accessibility_sign == 1),
                     cob_H3K4me3 = which(regions_all$cob_H3K4me3_sign == 1),
                     cob_H3K27ac = which(regions_all$cob_H3K27ac_sign == 1),
                     cob_H3K4me1 = which(regions_all$cob_H3K4me1_sign == 1))
u1 = upset(fromList(cov_sign_list), order.by = "freq")
u2 = upset(fromList(cob_sign_list), order.by = "freq")
script_plot_dir = "./plots_v1/"
push_pdf(u1, file_name = "upset_cov_wide", w = 7, h = 3)
push_pdf(u2, file_name = "upset_cob_wide", w = 7, h = 3)

# Figure 2f: CO-V overlap between modalities
g1 = regions_all %>% count(cov_chromatin_accessibility_sign, cov_H3K27ac_sign) %>%
        group_by(cov_H3K27ac_sign) %>%
        mutate(cov_chromatin_accessibility_sign = factor(cov_chromatin_accessibility_sign),
               frac = n / sum(n)) %>%
        ggbarplot(x = "cov_H3K27ac_sign", fill = "cov_chromatin_accessibility_sign", y = "frac", color = NA)
g2 = regions_all %>% count(cov_chromatin_accessibility_sign, cov_H3K4me1_sign) %>%
        group_by(cov_H3K4me1_sign) %>%
        mutate(cov_chromatin_accessibility_sign = factor(cov_chromatin_accessibility_sign),
               frac = n / sum(n)) %>%
        ggbarplot(x = "cov_H3K4me1_sign", fill = "cov_chromatin_accessibility_sign", y = "frac", color = NA)
g3 = regions_all %>% count(cov_chromatin_accessibility_sign, cov_H3K4me3_sign) %>%
        group_by(cov_H3K4me3_sign) %>%
        mutate(cov_chromatin_accessibility_sign = factor(cov_chromatin_accessibility_sign),
               frac = n / sum(n)) %>%
        ggbarplot(x = "cov_H3K4me3_sign", fill = "cov_chromatin_accessibility_sign", y = "frac", color = NA)
(g1 + g2 + g3) %>% push_pdf("multimodal_enrich_cov", w=  5, h = 3)

# CO-B synergy between modalities (Not included)
g1 = regions_all %>% count(cob_chromatin_accessibility_sign, cob_H3K27ac_sign) %>%
        group_by(cob_H3K27ac_sign) %>%
        mutate(cob_chromatin_accessibility_sign = factor(cob_chromatin_accessibility_sign),
               frac = n / sum(n)) %>%
        ggbarplot(x = "cob_H3K27ac_sign", fill = "cob_chromatin_accessibility_sign", y = "frac", color = NA)
g2 = regions_all %>% count(cob_chromatin_accessibility_sign, cob_H3K4me1_sign) %>%
        group_by(cob_H3K4me1_sign) %>%
        mutate(cob_chromatin_accessibility_sign = factor(cob_chromatin_accessibility_sign),
               frac = n / sum(n)) %>%
        ggbarplot(x = "cob_H3K4me1_sign", fill = "cob_chromatin_accessibility_sign", y = "frac", color = NA)
g3 = regions_all %>% count(cob_chromatin_accessibility_sign, cob_H3K4me3_sign) %>%
        group_by(cob_H3K4me3_sign) %>%
        mutate(cob_chromatin_accessibility_sign = factor(cob_chromatin_accessibility_sign),
               frac = n / sum(n)) %>%
        ggbarplot(x = "cob_H3K4me3_sign", fill = "cob_chromatin_accessibility_sign", y = "frac", color = NA)
(g1 + g2 + g3) %>% push_pdf("multimodal_enrich_cob", w=  5, h = 3)

# Raw number of elements called, not included
cov_sum = regions_all %>% filter(cov_any) %>%
        mutate(multi = (cov_chromatin_accessibility_sign + cov_H3K27ac_sign + cov_H3K4me3_sign + cov_H3K4me1_sign >= 2)) %>%
        count(multi) %>%
        mutate(type = "cov")
cob_sum = regions_all %>% filter(cob_any) %>%
        mutate(multi = (cob_chromatin_accessibility_sign + cob_H3K27ac_sign + cob_H3K4me3_sign + cob_H3K4me1_sign >= 2)) %>%
        count(multi) %>%
        mutate(type = "cob")

bind_rows(tibble(modality = all_assays,
                 count = c(sum(regions_all$cov_chromatin_accessibility_sign),
                           sum(regions_all$cov_H3K27ac_sign),
                           sum(regions_all$cov_H3K4me1_sign),
                           sum(regions_all$cov_H3K4me3_sign)),
                 type = "cov"),
          tibble(modality = all_assays,
                 count = c(sum(regions_all$cob_chromatin_accessibility_sign),
                           sum(regions_all$cob_H3K27ac_sign),
                           sum(regions_all$cob_H3K4me1_sign),
                           sum(regions_all$cob_H3K4me3_sign)),
                 type = "cob")) %>%
        ggbarplot(x = "modality",
                  y = "count",
                  fill = "modality",
                  facet.by = "type") +
        theme(legend.position = "right", axis.text.x = element_text(angle = 90))

# Figure S4: Venn diagram showing overlaps of CO-V and CO-B
library(ggVennDiagram)
library(RColorBrewer)
regions_all$human_cre = factor(regions_all$human_cre_class_prio_v4 != "non-cCRE", labels = c("non-cCRE", "cCRE"))
regions_all$mouse_cre = factor(regions_all$mouse_cre_class_prio_v4 != "non-cCRE", labels = c("non-cCRE", "cCRE"))
regions_all$high_spearman = factor(regions_all$cov_any == 1, labels = c("Other", "High CO-V"))
regions_all$high_housekeeping = factor(regions_all$cob_any == 1, labels = c("Other", "High CO-B"))

table(regions_all$human_cre, regions_all$high_spearman == "High CO-V" & regions_all$high_housekeeping == "High CO-B")
table(regions_all$mouse_cre, regions_all$high_spearman != "High CO-V" & regions_all$high_housekeeping == "High CO-B")

g_theme = theme(text = element_text(size = 8),
                legend.position = "none",
                strip.background = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(),
                axis.text.x = element_text(angle = 60, margin = margin(t = 20)))

v1 = ggVennDiagram(list("CO-V" = which(regions_all$cov_chromatin_accessibility_sign == 1),
                        "CO-B" = which(regions_all$cob_chromatin_accessibility_sign == 1))) +
        scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9, "Blues"))(9), 
                             name = "Values", 
                             na.value = "white") +
        scale_color_manual(values = c("#9f1d52", "#fecf0f"))
v2 = ggVennDiagram(list("CO-V" = which(regions_all$cov_H3K4me3_sign == 1),
                        "CO-B" = which(regions_all$cob_H3K4me3_sign == 1))) +
        scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9, "Blues"))(9), 
                             name = "Values", 
                             na.value = "white") +
        scale_color_manual(values = c("#9f1d52", "#fecf0f"))
v3 = ggVennDiagram(list("CO-V" = which(regions_all$cov_H3K27ac_sign == 1),
                        "CO-B" = which(regions_all$cob_H3K27ac_sign == 1))) +
        scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9, "Blues"))(9), 
                             name = "Values", 
                             na.value = "white") +
        scale_color_manual(values = c("#9f1d52", "#fecf0f"))
v4 = ggVennDiagram(list("CO-V" = which(regions_all$cov_H3K4me1_sign == 1),
                        "CO-B" = which(regions_all$cob_H3K4me1_sign == 1))) +
        scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9, "Blues"))(9), 
                             name = "Values", 
                             na.value = "white") +
        scale_color_manual(values = c("#9f1d52", "#fecf0f"))
(v1 / v2 / v3 / v4) %>% push_pdf("cov_cob_overlap_venn", w=  4, h = 8)

# Evaluate calling de novo functional elements
regions_all$cov_any4 = rowSums(as.matrix(regions_all[, paste0("cov_", assay_name_mapping[c(1, 2, 3, 4)], "_sign")])) > 0
regions_all$cob_any4 = rowSums(as.matrix(regions_all[, paste0("cob_", assay_name_mapping[c(1, 2, 3, 4)], "_sign")])) > 0
regions_all$high_spearman = factor(regions_all$cov_any4 == 1, labels = c("Other", "High CO-V"))
regions_all$high_housekeeping = factor(regions_all$cob_any4 == 1, labels = c("Other", "High CO-B"))
regions_all %>%
        filter(!((high_spearman == "High CO-V") & (high_housekeeping != "High CO-B"))) %>%
        group_by(human_cre, high_housekeeping) %>%
        summarise(mean_mouse_cre = mean(mouse_cre == "cCRE", na.rm = T))
regions_all %>%
        filter(!((high_spearman == "High CO-V") & (high_housekeeping != "High CO-B"))) %>%
        group_by(human_cre, high_spearman) %>%
        summarise(mean_mouse_cre = mean(mouse_cre == "cCRE", na.rm = T),
                  count = n())
regions_all %>%
        filter(!((high_spearman == "High CO-V") & (high_housekeeping != "High CO-B"))) %>%
        group_by(mouse_cre, high_spearman) %>%
        summarise(mean_human_cre = mean(human_cre == "cCRE", na.rm = T),
                  count = n())
# de novo elements: subset to conserved elements
regions_all$high_any = ifelse(regions_all$high_spearman == "High CO-V" | regions_all$high_housekeeping == "High CO-B",
                              "Conserved", "non-Conserved")
regions_hico = regions_all[regions_all$high_spearman == "High CO-V" | regions_all$high_housekeeping == "High CO-B", ]
regions_all$human_dhs = ifelse(!is.na(regions_all$human_id), "Human DHS", "Human non-DHS")
regions_all$mouse_dhs = ifelse(!is.na(regions_all$mouse_id), "Mouse DHS", "Mouse non-DHS")

regions_all$mouse_cre_dhs = factor(paste0(regions_all$mouse_cre, "_", regions_all$mouse_dhs),
                                   labels = c("DHS & cCRE", "Only cCRE", "Only DHS", "non-DHS, non-cCRE"))
regions_all$human_cre_dhs = factor(paste0(regions_all$human_cre, "_", regions_all$human_dhs),
                                   labels = c("DHS & cCRE", "Only cCRE", "Only DHS", "non-DHS, non-cCRE"))

# Figure 5a: Pie chart breaking down current annotations
g1 = regions_all %>% count(human_cre_dhs) %>%
        # mutate(label = paste0(human_cre_dhs, ": ", n)) %>%
        ggpie(x = "n", label = "n", fill = "human_cre_dhs") +
        scale_fill_brewer(type = "qual") +
        theme(legend.position = "right")

g2 = regions_all %>% count(human_cre_dhs, high_any) %>%
        filter(high_any == "Conserved") %>%
        # mutate(label = paste0(human_cre_dhs, ": ", n)) %>%
        ggpie(x = "n", label = "n", fill = "human_cre_dhs") +
        scale_fill_brewer(type = "qual") +
        theme(legend.position = "right")

g3 = regions_all %>% count(mouse_cre_dhs) %>%
        ggpie(x = "n", label = "n", fill = "mouse_cre_dhs") +
        scale_fill_brewer(type = "qual") +
        theme(legend.position = "right")

g4 = regions_all %>% count(mouse_cre_dhs, high_any) %>%
        filter(high_any == "Conserved") %>%
        ggpie(x = "n", label = "n", fill = "mouse_cre_dhs") +
        scale_fill_brewer(type = "qual") +
        theme(legend.position = "right")
((g1 + g2 + g3 + g4) + plot_layout(guide = "collect", nrow = 1)) %>%
        push_pdf("dhs_cre_pie", w = 8.5, h = 2.5)

# Figure 5b: alignment rates for de novo elements
g1 = regions_all %>%
        group_by(human_cre_dhs, high_housekeeping) %>% summarise(mean_mouse_cre = mean(mouse_cre == "cCRE", na.rm = T)) %>%
        ggbarplot(x = "human_cre_dhs",
                  y = "mean_mouse_cre",
                  fill = "high_housekeeping",
                  position = position_dodge(),
                  col = NA,
                  xlab = "",
                  ylab = "Fraction aligned to mouse cCRE",
                  ylim = c(0, 1)) +
        # facet(facet.by = "human_cre_dhs", nrow = 1) + 
        scale_fill_manual(values = co_col) +
        g_theme

g2 = regions_all %>%
        group_by(mouse_cre_dhs, high_housekeeping) %>% summarise(mean_human_cre = mean(human_cre == "cCRE", na.rm = T)) %>%
        ggbarplot(x = "mouse_cre_dhs",
                  y = "mean_human_cre",
                  fill = "high_housekeeping",
                  position = position_dodge(),
                  col = NA,
                  ylim = c(0, 1),
                  xlab = "",
                  ylab = "Fraction aligned to human cCRE") +
        # facet(facet.by = "mouse_cre_dhs", nrow = 1) +
        scale_fill_manual(values = co_col) +
        g_theme

g3 = regions_all %>%
        group_by(human_cre_dhs, high_spearman) %>% summarise(mean_mouse_cre = mean(mouse_cre == "cCRE", na.rm = T)) %>%
        ggbarplot(x = "human_cre_dhs",
                  y = "mean_mouse_cre",
                  fill = "high_spearman",
                  position = position_dodge(),
                  col = NA,
                  xlab = "",
                  ylab = "Fraction aligned to mouse cCRE",
                  ylim = c(0, 1)) +
        # facet(facet.by = "human_cre_dhs", nrow = 1) +
        scale_fill_manual(values = co_col) +
        g_theme

g4 = regions_all %>%
        group_by(mouse_cre_dhs, high_spearman) %>% summarise(mean_human_cre = mean(human_cre == "cCRE", na.rm = T)) %>%
        ggbarplot(x = "mouse_cre_dhs",
                  y = "mean_human_cre",
                  fill = "high_spearman",
                  position = position_dodge(),
                  col = NA,
                  xlab = "",
                  ylab = "Fraction aligned to human cCRE",
                  ylim = c(0, 1)) +
        # facet(facet.by = "mouse_cre_dhs", nrow = 1) +
        scale_fill_manual(values = co_col) + 
        g_theme
((g3 + g4 + g1 + g2) + patchwork::plot_layout(nrow = 1)) %>% 
        push_pdf("denove_cre_aligned_cre_final", w = 8.5, h = 3.)

# Figure 5c: H3K4me1 evaluations
regions_all$cov_any3 = rowSums(as.matrix(regions_all[, paste0("cov_", assay_name_mapping[c(1, 2, 4)], "_sign")])) > 0
regions_all$cob_any3 = rowSums(as.matrix(regions_all[, paste0("cob_", assay_name_mapping[c(1, 2, 4)], "_sign")])) > 0
regions_all$high_spearman = factor(regions_all$cov_any3 == 1, labels = c("Other", "High CO-V"))
regions_all$high_housekeeping = factor(regions_all$cob_any3 == 1, labels = c("Other", "High CO-B"))

assay = "H3K4me1"
g1 = ggboxplot(regions_all, x = "human_cre_dhs", fill = "high_spearman", y = "cov_H3K4me1", ylim = c(-0.5, 1),
               size = 0.25, outlier.size = 0.01,
               ylab = paste0("CO-V (", assay, ")")) + scale_fill_manual(values = co_col) + g_theme
g2 = ggboxplot(regions_all, x = "mouse_cre_dhs", fill = "high_spearman", y = "cov_H3K4me1", ylim = c(-0.5, 1),
               size = 0.25, outlier.size = 0.01,
               ylab = paste0("CO-V (", assay, ")")) + scale_fill_manual(values = co_col) +g_theme
g3 = ggboxplot(regions_all, x = "human_cre_dhs", fill = "high_housekeeping", y = "cob_H3K4me1", ylim = c(0, 1),
               size = 0.25, outlier.size = 0.01,
               ylab = paste0("CO-B (", assay, ")")) + scale_fill_manual(values = co_col) + g_theme
g4 = ggboxplot(regions_all, x = "mouse_cre_dhs", fill = "high_housekeeping", y = "cob_H3K4me1", ylim = c(0, 1),
               size = 0.25, outlier.size = 0.01,
               ylab = paste0("CO-B (", assay, ")")) + scale_fill_manual(values = co_col) +g_theme
((g1 + g2 + g3 + g4) + plot_layout(guide = "collect", nrow = 1)) %>% 
        push_png(file_name = "h3k4me1_eval_final_v1", width = 8.5, height = 3.)

# Add GWAS disease variants
snp_all = readRDS("./metadata_processed/gwas/pruned_gwas_snps_noss.rds")
snp_all_unique = unique(snp_all)
regions = regions_all 
regions$human_cre = factor(regions$human_cre_class_prio_v4 != "non-cCRE", labels = c("non-cCRE", "cCRE"))
regions$mouse_cre = factor(regions$mouse_cre_class_prio_v4 != "non-cCRE", labels = c("non-cCRE", "cCRE"))
regions$human_dhs = ifelse(!is.na(regions$human_id), "Human DHS", "Human non-DHS")
regions$mouse_dhs = ifelse(!is.na(regions$mouse_id), "Mouse DHS", "Mouse non-DHS")
regions$mouse_cre_dhs = factor(paste0(regions$mouse_cre, "_", regions$mouse_dhs),
                               labels = c("DHS & cCRE", "Only cCRE", "Only DHS", "non-DHS, non-cCRE"))
regions$human_cre_dhs = factor(paste0(regions$human_cre, "_", regions$human_dhs),
                               labels = c("DHS & cCRE", "Only cCRE", "Only DHS", "non-DHS, non-cCRE"))
regions$cov_any4 = rowSums(as.matrix(regions[, paste0("cov_", assay_name_mapping[c(1, 2, 3, 4)], "_sign")])) > 0
regions$cob_any4 = rowSums(as.matrix(regions[, paste0("cob_", assay_name_mapping[c(1, 2, 3, 4)], "_sign")])) > 0
regions$high_spearman = factor(regions$cov_any4 == 1, labels = c("Other", "High CO-V"))
regions$high_housekeeping = factor(regions$cob_any4 == 1, labels = c("Other", "High CO-B"))
regions$high_any = ifelse(regions$high_spearman == "High CO-V" | regions$high_housekeeping == "High CO-B",
                          "Conserved", "non-Conserved")
regions$gwas_snp = ifelse(1:nrow(regions) %in% subjectHits(findOverlaps(snp_all_unique, hg38_regions)), 
                          1, 0)

regions_gwas = filter(regions, gwas_snp == 1)

regions %>% count(gwas_snp, high_any) %>%
         group_by(gwas_snp) %>%
        mutate(frac = n / sum(n)) %>% 
        filter(high_any  == "Conserved") %>%
        ggbarplot(x = "gwas_snp", y = "frac")


# Table S4: GREAT analysis
# GREAT analysis
library(rGREAT)
set.seed(73)
global_bg_indices = which(!is.na(regions_all$cov_any))
bg_indices = global_bg_indices
top_indices = which(regions_all$cov_any)
bg_indices = sample(bg_indices, size = 6e5)
top_indices = top_indices[top_indices %in% bg_indices]

great_job = submitGreatJob(gr = hg38_regions[top_indices],
                           bg = hg38_regions[bg_indices],
                           species = "hg38",
                           request_interval = 60)
tb_top = getEnrichmentTables(great_job)

tb_plot_top = tb_top$`GO Biological Process` %>% mutate(val = -log10(Hyper_Adjp_BH))
g_top = ggbarplot(tb_plot_top[25:1, ], x = "name", y = "val", orientation = "horiz",
                  fill = "Hyper_Fold_Enrichment",
                  xscale = "none",
                  xlab = "",
                  ylab = "-log10(Adj Hyper P-value)") +
        scale_fill_distiller(palette = "Oranges", direction = 1, limits = c(1, 2)) + theme(text = element_text(size = 7))
# push_pdf(g_top, file_name = "high_spearman_go_bp", ps = 7, w = 6, h = 4)
write_tsv(tb_top$`GO Molecular Function`, file = "./output/go_final/cov_mf.tsv")
write_tsv(tb_top$`GO Biological Process`, file = "./output/go_final/cov_bp.tsv")
write_tsv(tb_top$`GO Cellular Component`, file = "./output/go_final/cov_cc.tsv")
saveRDS(tb_top, file = "./output/go_final/cov.rds")
# GREAT CO-B:
set.seed(73)
global_bg_indices = which(!is.na(regions_all$cob_any))
bg_indices = global_bg_indices
top_indices = which(regions_all$cob_any)
bg_indices = sample(bg_indices, size = 6e5)
top_indices = top_indices[top_indices %in% bg_indices]

great_job = submitGreatJob(gr = hg38_regions[top_indices],
                           bg = hg38_regions[bg_indices],
                           species = "hg38",
                           request_interval = 60)
tb_top = getEnrichmentTables(great_job)

tb_plot_top = tb_top$`GO Biological Process` %>% mutate(val = -log10(Hyper_Adjp_BH))
g_top = ggbarplot(tb_plot_top[25:1, ], x = "name", y = "val", orientation = "horiz",
                  fill = "Hyper_Fold_Enrichment",
                  xscale = "none",
                  xlab = "",
                  ylab = "-log10(Adj Hyper P-value)") +
        scale_fill_distiller(palette = "Oranges", direction = 1, limits = c(1, 2)) + theme(text = element_text(size = 7))
# push_pdf(g_top, file_name = "high_spearman_go_bp", ps = 7, w = 6, h = 4)
write_tsv(tb_top$`GO Molecular Function`, file = "./output/go_final/cob_mf.tsv")
write_tsv(tb_top$`GO Biological Process`, file = "./output/go_final/cob_bp.tsv")
write_tsv(tb_top$`GO Cellular Component`, file = "./output/go_final/cob_cc.tsv")
saveRDS(tb_top, file = "./output/go_final/cob.rds")
