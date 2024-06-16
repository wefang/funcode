library(tidyverse)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
source("./helper/helper.R")
region_file = "./output/indexed_dhs_mapped_regions.rds"
# Chromatin accessibility
regions = readRDS(region_file)
exp_matched_biosample_new = readRDS(paste0("./custom_metadata/manual_matched_ca.rds"))
human_matched_acc = purrr::reduce(map(exp_matched_biosample_new$exp_tb.x, "Accession"), c)
human_matched_biosample = rep(paste0(exp_matched_biosample_new$`Life stage`, "_",
                                     exp_matched_biosample_new$biosample_mapped), map_dbl(exp_matched_biosample_new$exp_tb.x, nrow))
mouse_matched_acc = reduce(map(exp_matched_biosample_new$exp_tb.y, "Accession"), c)
mouse_matched_biosample = rep(paste0(exp_matched_biosample_new$`Life stage`, "_",
                                     exp_matched_biosample_new$biosample_mapped), map_dbl(exp_matched_biosample_new$exp_tb.y, nrow))
human_correct_manual = t(ghelper::aveMatFac(t(human_correct[, human_matched_acc]),
                                            human_matched_biosample))
mouse_correct_manual = t(ghelper::aveMatFac(t(mouse_correct[, mouse_matched_acc]),
                                            mouse_matched_biosample))
matched_biosample = intersect(colnames(human_correct_manual), colnames(mouse_correct_manual))
loci_cor_manual = map_dbl(1:nrow(human_correct), function(j) {
        if (j %% 1e5 == 0) message(j)
        x = human_correct_manual[j, matched_biosample]
        y = mouse_correct_manual[j, matched_biosample]
        if (var(x) > 0 & var(y) > 0) {
                cor(x, y, method = "spearman")
        } else {
                return(NA)
        }
})
regions$cov_manual_chromatin_accessiblity = loci_cor_manual
saveRDS(regions, file = "./output/indexed_dhs_mapped_regions_v1_manual.rds")
# H3K4me3
regions = readRDS("./output/indexed_dhs_mapped_regions_v1_manual.rds")
assay = "H3K27ac"
assay = "H3K4me3"
assay = "H3K4me1"
exp_matched_biosample_new = readRDS(paste0("./custom_metadata/manual_matched_", assay, ".rds"))
human_histone = readRDS(paste0("./processed_data/Histone_norm_correct/Homo_sapiens_", assay, "_norm_correct.rds"))
mouse_histone = readRDS(paste0("./processed_data/Histone_norm_correct/Mus_musculus_", assay, "_norm_correct.rds"))

human_matched_acc = purrr::reduce(map(exp_matched_biosample_new$exp_tb.x, "Accession"), c)
human_matched_biosample = rep(paste0(exp_matched_biosample_new$`Life stage`, "_",
                                     exp_matched_biosample_new$biosample_mapped), map_dbl(exp_matched_biosample_new$exp_tb.x, nrow))
mouse_matched_acc = reduce(map(exp_matched_biosample_new$exp_tb.y, "Accession"), c)
mouse_matched_biosample = rep(paste0(exp_matched_biosample_new$`Life stage`, "_",
                                     exp_matched_biosample_new$biosample_mapped), map_dbl(exp_matched_biosample_new$exp_tb.y, nrow))


human_correct_manual = t(ghelper::aveMatFac(t(human_histone[, human_matched_acc]),
                                            human_matched_biosample))
mouse_correct_manual = t(ghelper::aveMatFac(t(mouse_histone[, mouse_matched_acc]),
                                            mouse_matched_biosample))
matched_biosample = intersect(colnames(human_correct_manual), colnames(mouse_correct_manual))

loci_cor_manual = map_dbl(1:nrow(human_histone), function(j) {
        if (j %% 1e5 == 0) message(j)
        x = human_correct_manual[j, matched_biosample]
        y = mouse_correct_manual[j, matched_biosample]
        if (var(x) > 0 & var(y) > 0) {
                cor(x, y, method = "spearman")
        } else {
                return(NA)
        }
})
regions[[paste0("cov_manual_", assay)]] = loci_cor_manual
saveRDS(regions, file = "./output/indexed_dhs_mapped_regions_v1_manual.rds")
