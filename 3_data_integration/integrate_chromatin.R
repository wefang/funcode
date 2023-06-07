devtools::load_all()
library(tidyverse)
library(Seurat)
library(ggpubr)
library(ComplexHeatmap)
source("R_production/config.R")
source("R_production/def_color.R")

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

# region_filter = matrixStats::rowMaxs(human_correct) > 20 &
#         matrixStats::rowMaxs(mouse_correct) > 20
# human_var_std = calc_var_std_ver1(human_correct)
# mouse_var_std = calc_var_std_ver1(mouse_correct)
# human_mean = rowMeans(human_correct)
# mouse_mean = rowMeans(mouse_correct)
# data_cross_var_feature = regions$id[which(rank(human_var_std) > 1.3e6 &
#                                                   rank(mouse_var_std) > 1.3e6 &
#                                                   region_filter)]
# saveRDS(data_cross_var_feature, file = "./intermediate_data/ENCODE4_human_mouse_chromatin_var_features.rds")

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

# loci correlation
anchor_df = make_anchor_tb(anchors, human_correct_obj, mouse_correct_obj)
write_tsv(anchor_df, file = paste0(script_output_dir, "anchors.tsv"))

anchor_df = read_tsv(paste0(script_output_dir, "anchors.tsv"))
# computing weighted spearman correlation
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

# calling conserved elements
# computing empirical null distribution
null_size = 5e5
sample_indices1 = sample(nrow(human_correct), null_size)
sample_indices2 = sample(nrow(human_correct), null_size)
human_null = human_correct[sample_indices1, ]
mouse_null = mouse_correct[sample_indices2, ]

loci_cor_null = map_dbl(1:null_size, function(j) {
        if (j %% 1e5 == 0) message(j)
        x = human_null[j, anchor_df$accession1]
        y = mouse_null[j, anchor_df$accession2]
        if (var(x) > 0 & var(y) > 0) {
                wCorr::weightedCorr(x,
                                    y,
                                    weights = anchor_df$score,
                                    method = "spearman")
        } else {
                return(NA)
        }
})
# saveRDS(loci_cor_null, file = output_rds_path("spearman_null"))
spearman_pval = qvalue::empPvals(loci_cor, loci_cor_null)
spearman_pval_adj = p.adjust(spearman_pval, method = "BH")
sum(spearman_pval_adj < 0.1)
qval = qvalue::qvalue(spearman_pval)$qvalue
sum(qval < 0.1)


# computing stats with anchors
human_var_std = calc_var_std_ver1(human_correct[, anchor_df$accession1])
mouse_var_std = calc_var_std_ver1(mouse_correct[, anchor_df$accession2])
human_mean = rowMeans(human_correct[, anchor_df$accession1])
mouse_mean = rowMeans(mouse_correct[, anchor_df$accession2])

# save(loci_cor, spearman_pval_adj, region_filter,
#      human_mean, mouse_mean,
#      human_var_std, mouse_var_std,
#      file = paste0(script_output_dir, "encodev4_caco_scores.rda"))
# cor(regions$spearman_global, loci_cor, use = 'complete.obs')


# setting columns for region file
# regions$spearman_global = loci_cor
# regions$spearman_global_fdr = spearman_pval_adj
# regions$filter_global = region_filter
# regions$high_spearman = regions$spearman_global > min(regions$spearman_global[regions$spearman_global_fdr < 0.05])
# regions$human_chromatin_mean_global_anchors = human_mean
# regions$mouse_chromatin_mean_global_anchors = mouse_mean
# regions$human_chromatin_var_std_global_anchors = human_var_std
# regions$mouse_chromatin_var_std_global_anchors = mouse_var_std
# saveRDS(regions, file = region_file)
# end setting columns for region file

# max(regions$spearman_global_fdr[regions$spearman_global > 0.5], na.rm =  T)
sum(regions$high_spearman, na.rm = T)

(tibble(score = c(regions$spearman_global,
                  loci_cor_null),
        source = c(rep("ZAlt", nrow(regions)),
                   rep("Null", length(loci_cor_null)))) %>%
                ggdensity(x = "score", color = "source", fill = "source", xlab = "Variable Conservation Score") +
                geom_vline(xintercept = min(regions$spearman_global[which(regions$high_spearman)]), col = "#005792") +
                scale_color_manual(values = c("#868686FF", "#99154e"))+
                scale_fill_manual(values = c("#868686FF", "#99154e")) + theme(text = element_text(size = 10), legend.position = "none")) %>%
        push_pdf(file_name = "hist_spearman_0d1", w = 2.8, h = 2.0)

hk_cutoff = 0.9
# hk_cutoff = 0.8
# conserved housekeeping
human_cutoff = quantile(regions$human_chromatin_mean_global_anchors, hk_cutoff)
mouse_cutoff = quantile(regions$mouse_chromatin_mean_global_anchors, hk_cutoff)

hk_pct_human = matrixStats::rowWeightedMeans(human_correct[, anchor_df$accession1] > human_cutoff, anchor_df$score)
hk_pct_mouse = matrixStats::rowWeightedMeans(mouse_correct[, anchor_df$accession2] > mouse_cutoff, anchor_df$score)
hk_pct = pmin(hk_pct_human, hk_pct_mouse)
# constructing null
null_size = 5e5
sample_indices1 = sample(nrow(regions), null_size)
sample_indices2 = sample(nrow(regions), null_size)
house_keeping_score_null =
        pmin(matrixStats::rowWeightedMeans(human_correct[sample_indices1, anchor_df$accession1] > human_cutoff, anchor_df$score),
             matrixStats::rowWeightedMeans(mouse_correct[sample_indices2, anchor_df$accession2] > mouse_cutoff, anchor_df$score))
hk_pval = qvalue::empPvals(hk_pct, house_keeping_score_null)
hk_pval_adj = p.adjust(hk_pval, method = "BH")

# updateing region file housekeeping
regions$hk_pct = hk_pct
regions$hk_pct_fdr = hk_pval_adj
regions$high_housekeeping = regions$hk_pct > min(regions$hk_pct[which(regions$hk_pct_fdr < 0.1)])
# saveRDS(regions, file = region_file)
# end updateing region file

pval = spearman_pval
pval_adj = spearman_pval_adj
save(loci_cor, pval, pval_adj, hk_pct, hk_pval, hk_pval_adj,
     file = paste0("./output/histone_conservation/", "temp_original_", "dnase_atac", "_scores.rds"))

sum(regions$high_housekeeping, na.rm = T)

(tibble(score = c(regions$hk_pct,
                  house_keeping_score_null),
        source = c(rep("ZAlt", nrow(regions)),
                   rep("Null", length(house_keeping_score_null)))) %>%
                ggdensity(x = "score", color = "source", fill = "source", xlab = "Conserved Housekeeping Score") + xlim(c(0., 1)) + coord_cartesian(ylim = c(0, 2.)) +
                geom_vline(xintercept = min(regions$hk_pct[which(regions$high_housekeeping)]), col = "#005792") +
                scale_color_manual(values = c("#868686FF", "#ffd56b"))+
                scale_fill_manual(values = c("#868686FF", "#ffd56b")) + theme(text = element_text(size = 10), legend.position = "none")) %>%
        push_pdf(file_name = "hist_hk_0d1", w = 2.8, h = 2.0)

# TODO: Heatmap below
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

heatmap_wrap <- function(sample_indices) {
        human_mat = t(scale(t(log2(1+human_correct[sample_indices, anchor_tb_sel$accession1]))))
        mouse_mat = t(scale(t(log2(1+mouse_correct[sample_indices, anchor_tb_sel$accession2]))))
        col_ord = hclust(dist(t(rbind(human_mat, mouse_mat))))
        
        Heatmap(human_mat,
                cluster_columns = col_ord,
                show_column_names = F,
                show_row_names = F,
                col = val_col,
                name = "Scaled Count",
                show_row_dend = F,
                show_column_dend = F,
                bottom_annotation = HeatmapAnnotation(Assay = anchor_tb_sel$`Assay title_1`,
                                                      col = list(Assay = assay_col_hu),
                                                      annotation_label = "Human Assay",
                                                      show_annotation_name = F)
        ) + 
                Heatmap(mouse_mat,
                        cluster_columns = col_ord,
                        show_column_names = F,
                        show_row_names = F,
                        col = val_col,
                        name = "Scaled Count",
                        show_row_dend = F,
                        show_column_dend = F,
                        bottom_annotation = HeatmapAnnotation(Assay =  anchor_tb_sel$`Assay title_2`,
                                                              col = list(Assay = assay_col_mm),
                                                              annotation_label = "Mouse_Assay",
                                                              show_annotation_name = F))#,
        # right_annotation = rowAnnotation(Component = regions$component[sample_indices0],
        #                                  col = list(Component = comp_col)))
}
# TODO: NEED to append some old annotation data for mouse
heatmap_wrap(sample_indices0) %>% push_png(file_name = "heatmap_high_spearman", w = 5., h = 4., res = 300)
heatmap_wrap(sample_indices1) %>% push_png(file_name = "heatmap_low_spearman", w = 5., h = 1.5, res = 300)


example_loci = which(regions$spearman_global > 0.5 &
                             regions$hk_pct < 0.15)
sort(table(reduce(map(regions$homolog_gene_tb[example_loci],
                      function(x) x$gene_human), c)), decreasing = T)[1:20]

example_gene = "PBX1"
hg38_regions[example_loci[
        map_lgl(regions[example_loci, ]$homolog_gene_tb, function(x) example_gene %in% x$gene_human)]]
mm10_regions[example_loci[
        map_lgl(regions[example_loci, ]$homolog_gene_tb, function(x) example_gene %in% x$gene_human)]]

ol = findOverlaps(mm10_regions,
                  GRanges(seqnames = "chr12", IRanges(start = 12810000,
                                                      end = 12815000)))
mm10_regions[queryHits(ol)]
hg38_regions[queryHits(ol)]

