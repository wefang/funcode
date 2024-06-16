library(tidyverse)
library(Seurat)
library(ggpubr)
library(ComplexHeatmap)
source("./helper/helper.R")
source("./helper/def_color.R")
region_file = "./output/indexed_dhs_mapped_regions.rds"
regions = readRDS(region_file)
script_output_dir = "./output/ENCODE4_conservation/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
script_plot_dir = "./plots/ENCODE4_conservation/"
if (!dir.exists(script_plot_dir)) {
        dir.create(script_plot_dir)
}
# batch corrected data
human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_correct.rds")
mouse_correct = readRDS("./intermediate_data/encodev4_mouse_chromatin_correct.rds")
# human_correct_ext = readRDS("./intermediate_data/encodev4_human_chromatin_mouse_mapped_dhs_correct.rds")
# mouse_correct_ext = readRDS("./intermediate_data/encodev4_mouse_chromatin_mouse_mapped_dhs_correct.rds")
rownames(human_correct) = rownames(mouse_correct) = regions$id

meta = read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
meta$`Tissue/cell types` = NA
region_filter = matrixStats::rowMaxs(human_correct) > 20 &
        matrixStats::rowMaxs(mouse_correct) > 20
human_var_std = calc_var_std_ver1(human_correct)
mouse_var_std = calc_var_std_ver1(mouse_correct)
human_mean = rowMeans(human_correct)
mouse_mean = rowMeans(mouse_correct)
data_cross_var_feature = regions$id[which(rank(human_var_std) > 1.3e6 &
                                                  rank(mouse_var_std) > 1.3e6 &
                                                  region_filter)]
saveRDS(data_cross_var_feature, file = "./intermediate_data/ENCODE4_human_mouse_chromatin_var_features.rds")

data_cross_var_feature = readRDS("./intermediate_data/ENCODE4_human_mouse_chromatin_var_features.rds")
human_correct_obj = CreateSeuratObject(counts = human_correct[data_cross_var_feature, ],
                                       assay = "OpenChromatin")
human_correct_obj = add_metadata(human_correct_obj, meta)
mouse_correct_obj = CreateSeuratObject(counts = mouse_correct[data_cross_var_feature, ],
                                       assay = "OpenChromatin")
mouse_correct_obj = add_metadata(mouse_correct_obj, meta)
human_correct_obj = FindVariableFeatures(human_correct_obj, nfeatures = 40000)
mouse_correct_obj = FindVariableFeatures(mouse_correct_obj, nfeatures = 40000)
anchors <-
        FindIntegrationAnchors(
                list(human = human_correct_obj,
                     mouse = mouse_correct_obj),
                anchor.features = 10000,
                k.anchor = 5,
                k.filter = 20,
                k.score = 30,
                dims = 1:30)
integrated = IntegrateData(anchors,
                           dims = 1:15,
                           k.weight = 20)
integrated$organism = ifelse(integrated$Accession %in% human_correct_obj$Accession, "human", "mouse")
integrated$source = paste0(integrated$organism, "_", integrated$`Assay title`)
integrated = quick_process(integrated, anchors@anchor.features, dims = 1:15)
# integrated <- FindNeighbors(integrated, dims = 1:15)
# integrated <- FindClusters(integrated)
# split(paste0(integrated$`Biosample term name`, "...",integrated$source, integrated$`Life stage`), integrated$seurat_clusters)
# table(integrated$source, integrated$seurat_clusters)

(UMAPPlot(integrated, pt.size = 0.5, group.by = "organism") +
                # scale_color_manual(values = chromatin_col) +
                theme(text = element_text(size = 10), legend.position = "top", legend.text = element_text(size = 8))) %>%
        push_pdf(., file_name = paste0("umap_integrated_split_ca"), w = 3, h = 3.5, ps = 10)
}

(UMAPPlot(integrated, pt.size = 0.5, group.by = "source") +
                # scale_color_manual(values = chromatin_col) +
                theme(text = element_text(size = 10), legend.position = "top", legend.text = element_text(size = 8))) # %>%
        push_pdf(., file_name = "umap_integrated", w = 3.2, h = 3.5, ps = 10)

HoverLocator(plot = UMAPPlot(integrated, pt.size = 3, group.by = "source"), # + scale_color_manual(values = chromatin_col),
             information = FetchData(integrated,
                                     vars = c("source", "Biosample term name", "Life stage"))) # %>% 
        push_widget(file_name = "integrated")
# defining in silico matched pairs
anchor_df = make_anchor_tb(anchors, human_correct_obj, mouse_correct_obj)
write_tsv(anchor_df, file = paste0(script_output_dir, "anchors.tsv"))
anchor_df = read_tsv(paste0(script_output_dir, "anchors.tsv"))

# Computing initial correlation scores for visualization
loci_cor = sapply(1:nrow(human_correct), function(j) {
        if (j %% 1e5 == 0) message(j)
        x = human_correct[j, anchor_df$accession1]
        y = mouse_correct[j, anchor_df$accession2]
        if (var(x) > 0 & var(y) > 0) {
                wCorr::weightedCorr(x,
                                    y,
                                    weights = anchor_df$score,
                                    method = "spearman")
        } else {
                return(NA)
        }
})
# Visualize scatter plot
set.seed(73)
val_max = 13
g_theme = theme(text = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 30))
plot_index = sample(which(abs(regions$spearman_global) > 0.8 & regions$hk_pct < 0.5), 1)
print(regions$spearman_global[plot_index])
print(regions$hk_pct[plot_index])
g1 = (tibble(human = log2(1+human_correct[plot_index, anchor_df$accession1]),
             mouse = log2(1+mouse_correct[plot_index, anchor_df$accession2])) %>%
              ggscatter(x = "human", y = "mouse", xlim = c(0, val_max), ylim = c(0, val_max),
                        xlab = "Human", ylab = "Mouse", size = 0.25,
              )  + geom_abline(a = 1, b = 0) +
              geom_smooth(method = "lm", se = F) + g_theme)
plot_index = sample(which(abs(regions$spearman_global) < 0.05 & (regions$hk_pct < 0.3 & regions$hk_pct > 0.1)), 1)
print(regions$spearman_global[plot_index])
print(regions$hk_pct[plot_index])
g2 = (tibble(human = log2(1+human_correct[plot_index, anchor_df$accession1]),
             mouse = log2(1+mouse_correct[plot_index, anchor_df$accession2])) %>%
              ggscatter(x = "human", y = "mouse",
                        xlab = "Human", ylab = "", size = 0.25,
                        xlim = c(0, val_max), ylim = c(0, val_max)
              )   + geom_abline(a = 1, b = 0) + geom_smooth(method = "lm", se = F) + g_theme)

plot_index = which.max(regions$hk_pct)
print(regions$spearman_global[plot_index])
print(regions$hk_pct[plot_index])
g3 = (tibble(human = log2(1+human_correct[plot_index, anchor_df$accession1]),
             mouse = log2(1+mouse_correct[plot_index, anchor_df$accession2])) %>%
              ggscatter(x = "human", y = "mouse", xlim = c(0, val_max), ylim = c(0, val_max),
                        xlab = "Human", ylab = "", size = 0.25
              )   + geom_abline(a = 1, b = 0) + geom_smooth(method = "lm", se = F) + g_theme)
(g1 + g3 + g2) %>% push_pdf(file_name = "caco_scatter", w = 4.86, h = 1.6)

# Visualize heatmaps
anchor_tb = read_tsv(paste0(script_output_dir, "anchors.tsv"))
cbind(anchor_tb$`Biosample term name_1`, anchor_tb$`Biosample term name_2`)

library(ComplexHeatmap)
set.seed(73)
val_col = circlize::colorRamp2(c(-2, 0, 2),
                               colors = c("blue", "white", "red"))
anchor_tb_sel = filter(anchor_tb, score > 0.25)
sample_indices0 = sample(which(regions$high_spearman &
                                       regions$human_chromatin_var_std_global_anchors > 1.2 &
                                       regions$mouse_chromatin_var_std_global_anchors > 1.2),
                         size = 1000)
sample_indices1 = sample(which(regions$spearman_global < 0.2 &
                                       regions$human_chromatin_var_std_global_anchors > 1.2 &
                                       regions$mouse_chromatin_var_std_global_anchors > 1.2), size = 250)
heatmap_wrap(sample_indices0) %>% push_png(file_name = "heatmap_high_spearman", w = 5., h = 4., res = 300)
heatmap_wrap(sample_indices1) %>% push_png(file_name = "heatmap_low_spearman", w = 5., h = 1.5, res = 300)
