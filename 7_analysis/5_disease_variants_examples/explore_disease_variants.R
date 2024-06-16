# Visualisze examples of disease variants
library(GenomicRanges)
library(ComplexHeatmap)
library(tidyverse)
library(ggpubr)
source("./helper/helper.R")

human_correct = readRDS("./intermediate_data/human_chromatin_correct.rds")
mouse_correct = readRDS("./intermediate_data/mouse_chromatin_correct.rds")

anchors_all = read_tsv(paste0("./output/chromatin_conservation/", "anchors.tsv"))
# anchors <- readRDS(paste0('./non_alignable_env/raw_data/anchors_chromatin_v1.rds'))
print(load("metadata_processed/indexed_aligned_combined.rda"))
regions_all = readRDS("./output/indexed_dhs_mapped_regions_v1_manual.rds")
regions = regions_all

# Flagship Figure: ATOH example
region_index = subjectHits(findOverlaps(ghelper::str_to_gr("chr10:68251200-68253000"),
                                        hg38_regions))
region_index = region_index[region_index < nrow(human_correct)]
# region_index = region_index[regions_all$cov_chromatin_accessibility[region_index] > 0.25]
# regions_all$cov_chromatin_accessibility_fdr[region_index]
h_mat = t(scale(t(log2(1 + human_correct[region_index, anchors_all$accession1]))))
m_mat = t(scale(t(log2(1 + mouse_correct[region_index, anchors_all$accession2]))))

Heatmap(h_mat) + 
        Heatmap(m_mat)

anchors_all$`Biosample term name_1`[order(colSums(h_mat), decreasing = T)][1:50]
human_correct_ext = readRDS("./intermediate_data/encodev4_human_chromatin_mouse_mapped_dhs_correct.rds")
mouse_correct_ext = readRDS("./intermediate_data/encodev4_mouse_chromatin_mouse_mapped_dhs_correct.rds")
human_correct_ext = human_correct_ext[regions_all$mouse_id[which(grepl("MouseDHS-", regions_all$mouse_id))], ]
mouse_correct_ext = mouse_correct_ext[regions_all$mouse_id[which(grepl("MouseDHS-", regions_all$mouse_id))], ]

x_val = colMeans(human_correct[region_index, anchors_all$accession1])
y_val = colMeans(mouse_correct[region_index, anchors_all$accession2])
region_index
library(plotly)
signal_tb = tibble(x = log2(1+x_val),
                   y = log2(1+y_val),
                   biosample1 = paste0(anchors_all$`Life stage_1`, "_", anchors_all$`Biosample term name_1`),
                   biosample2 = paste0(anchors_all$`Life stage_2`, "_", anchors_all$`Biosample term name_2`)) %>%
        mutate(biosample_pair = paste0(biosample1, "-", biosample2))
relevant_bio = c("embryonic_retina", "adult_retina", "embryonic_eye", "adult_eye", "postnatal_retina",
                 "child_WERI-Rb-1", "embryonic_forebrain", "child_BE2C", "embryonic_neural tube")
signal_tb$biosample_pair[!(signal_tb$biosample1 %in%  relevant_bio &
                                 signal_tb$biosample2 %in% relevant_bio)] = "zother"
(signal_tb %>%
        ggscatter(x = "x", y=  "y", color = "biosample_pair") + theme(legend.position = "none")) %>%
        ggplotly()
c("#FFFFCC", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D",
           "#238443", "#006837", "#004529", "#0C2C84", "#081D58")
(signal_tb %>%
                ggscatter(x = "x", y=  "y", color = "biosample_pair") +
                scale_color_manual(values = rev(c("#9e9e9e", "#FFFFCC", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D",
                                                           "#238443", "#006837", "#004529", "#0C2C84", "#081D58"))
                ) +
                geom_smooth() + theme(legend.position = "right")) %>%
        push_pdf("anchor_signals_retina_Atoh7", w = 6, h = 3)


# Figure 4f: SCN10A Example
region_index = subjectHits(findOverlaps(ghelper::str_to_gr("chr3:38726000-38726001"),
                                        hg38_regions))
region_index_rev = region_index[region_index > nrow(human_correct)]
region_index_rev = region_index_rev - nrow(human_correct)

anchors_rev <- readRDS(paste0('./non_alignable_env/raw_data/anchors_chromatin_v1.rds'))

h_mat = t(scale(t(log2(1 + human_correct_ext[region_index_rev, anchors_rev$accession1]))))
m_mat = t(scale(t(log2(1 + mouse_correct_ext[region_index_rev, anchors_rev$accession2]))))

h_mat = log2(1 + human_correct_ext[region_index_rev, anchors_rev$accession1])
m_mat = log2(1 + mouse_correct_ext[region_index_rev, anchors_rev$accession2])

Heatmap(h_mat) + 
        Heatmap(m_mat)
anchors_rev$Biosample.term.name_1[order(h_mat, decreasing = T)][1:50]
anchors_rev$Biosample.term.name_2[order(m_mat, decreasing = T)][1:50]

# human_correct_ext = readRDS("./intermediate_data/encodev4_human_chromatin_mouse_mapped_dhs_correct.rds")
# mouse_correct_ext = readRDS("./intermediate_data/encodev4_mouse_chromatin_mouse_mapped_dhs_correct.rds")
# human_correct_ext = human_correct_ext[regions_all$mouse_id[which(grepl("MouseDHS-", regions_all$mouse_id))], ]
# mouse_correct_ext = mouse_correct_ext[regions_all$mouse_id[which(grepl("MouseDHS-", regions_all$mouse_id))], ]
x_val = colMeans(human_correct_ext[region_index_rev, anchors_rev$accession1, drop = F])
y_val = colMeans(mouse_correct_ext[region_index_rev, anchors_rev$accession2, drop = F])
library(plotly)
signal_tb = tibble(x = log2(1+x_val),
                   y = log2(1+y_val),
                   biosample1 = paste0(anchors_rev$Life.stage_1, "_", anchors_rev$Biosample.term.name_1),
                   biosample2 = paste0(anchors_rev$Life.stage_2, "_", anchors_rev$Biosample.term.name_2)) %>%
        mutate(biosample_pair = paste0(biosample1, "-", biosample2))

relevant_bio = c("embryonic_heart", "adult_heart", "unknown_heart", "adult_heart left ventricle",
                 "adult_heart right ventricle",
                 "embryonic_cardiac muscle cell", "adult_Right ventricle myocardium inferior", "child heart",
                 "postnatal_heart")
signal_tb$biosample_pair[!(signal_tb$biosample1 %in%  relevant_bio &
                                   signal_tb$biosample2 %in% relevant_bio)] = "zother"
(signal_tb %>%
                ggscatter(x = "x", y=  "y", color = "biosample_pair") + geom_smooth() + theme(legend.position = "none")) %>%
        ggplotly()
(signal_tb %>%
                ggscatter(x = "x", y=  "y", color = "biosample_pair") +
                scale_color_manual(values = rev(c("#9e9e9e", "#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1",
                                                           "#DD3497", "#AE017E", "#7A0177", "#49006A", "#34003B", "#220020"))
                ) +
                geom_smooth() + theme(legend.position = "right")) %>%
        push_pdf("anchor_signals_heart_SCN5A", w = 7, h = 3)

# PRS example
# below is converted to hg38 from the original publication: Gordon, C. T. et al., 2014
test_gr = ghelper::str_to_gr(c("chr17:70661789-70662446",
                               "chr17:70674519-70675103",
                               "chr17:70738632-70739724",
                               "chr17:70761402-70762170",
                               "chr17:70776217-70777095",
                               "chr17:70934150-70934797",
                               "chr17:71181119-71181867",
                               "chr17:71709401-71713190",
                               "chr17:72119227-72119680",
                               "chr17:72375216-72377287"))
test_ind = c(F, T, F, F, F, F, F, T, T, T)
region_index1 = subjectHits(findOverlaps(test_gr[test_ind],
                                         hg38_regions))
region_index1 = region_index1[region_index1 < nrow(human_correct)]
region_index2 = subjectHits(findOverlaps(test_gr[!test_ind],
                                         hg38_regions))
region_index2 = region_index2[region_index2 < nrow(human_correct)]

summary_tb = tibble(regions_index = c(region_index1, region_index2),
                    cov_ca = c(regions$cov_chromatin_accessibility[region_index1],
                            regions$cov_chromatin_accessibility[region_index2]),
                    cov_ca_sign = c(regions$cov_chromatin_accessibility_sign[region_index1],
                                 regions$cov_chromatin_accessibility_sign[region_index2]),
                    cov_H3K27ac = c(regions$cov_H3K27ac[region_index1],
                               regions$cov_H3K27ac[region_index2]),
                    cov_H3K27ac_sign = c(regions$cov_H3K27ac_sign[region_index1],
                                    regions$cov_H3K27ac_sign[region_index2]),
                    phastCons4way = c(regions$phastCons4way[region_index1],
                                      regions$phastCons4way[region_index2]),
                    phyloP4way = c(regions$phyloP4way[region_index1],
                                   regions$phyloP4way[region_index2]),
                    ind = c(rep(T, length(region_index1)),
                            rep(F, length(region_index2))))


summary_tb %>% group_by(cov_ca_sign | cov_H3K27ac_sign) %>%
        summarize(mean(ind))
summary_tb %>% group_by(cov_H3K27ac_sign) %>%
        summarize(mean(ind))
cov_ca_cutoff = min(regions_all$cov_chromatin_accessibility[regions_all$cov_chromatin_accessibility_sign > 0])
cov_H3K27ac_cutoff = min(regions_all$cov_H3K27ac[regions_all$cov_H3K27ac_sign > 0])

summary_tb %>% arrange(desc(cov_ca)) %>% mutate(rank = 1:nrow(summary_tb)) %>%
        ggline(x = "rank", y = "cov_ca", point.color = "ind") +
        geom_hline(yintercept = cov_ca_cutoff)

summary_tb %>% arrange(desc(lecif)) %>% mutate(rank = 1:nrow(summary_tb)) %>%
        ggline(x = "rank", y = "lecif", point.color = "ind") #+
        # geom_hline(yintercept = cov_ca_cutoff)
(summary_tb %>% arrange(ind) %>%
        ggscatterhist(x = "cov_ca", y = "cov_H3K27ac", color = "ind", margin.plot = "boxplot",
                      margin.params = list(fill = "ind", color = "black", size = 0.2)) +
                geom_hline(yintercept = cov_H3K27ac_cutoff) +
                geom_vline(xintercept = cov_ca_cutoff)) %>%
        push_pdf("PRS_scatter", width = 3.5, height = 3.5)
g1 = (summary_tb %>%
              ggboxplot(x = "ind", y = "cov_ca", add = "jitter", add.params = list(color = "ind"))+
              stat_compare_means())

g2 = (summary_tb %>%
              ggboxplot(x = "ind", y = "cov_H3K27ac", add = "jitter", add.params = list(color = "ind")) +
              stat_compare_means())

((g1 + g2)& theme(legend.position = "top")) + plot_layout(guide = "collect")
# end PRS example
