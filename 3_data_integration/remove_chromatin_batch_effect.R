library(Seurat)
library(ComplexHeatmap)
library(tidyverse)
source("./R_production/config.R")
source("./R_production/def_color.R")
source("./R/helper.R")

script_plot_dir = "./plots/batch_correct_encodev4/"
if (!dir.exists(script_plot_dir)) {
        dir.create(script_plot_dir)
}

regions = readRDS(region_file)
# human_dnase = readRDS("E:\\GlobusDownload/Histone_norm_correct/Homo_sapiens_DNase_norm_correct.rds")
# human_atac = readRDS("E:\\GlobusDownload/Histone_norm_correct/Homo_sapiens_ATAC_norm_correct.rds")
# mouse_dnase = readRDS("E:\\GlobusDownload/Histone_norm_correct/Mus_musculus_DNase_norm_correct.rds")
# mouse_atac = readRDS("E:\\GlobusDownload/Histone_norm_correct/Mus_musculus_ATAC_norm_correct.rds")
# rownames(human_dnase) = rownames(human_atac) = regions$id
# rownames(mouse_dnase) = rownames(mouse_atac) = regions$id

# human_dnase = readRDS("E:\\GlobusDownload/Histone_norm_all_dhs/Homo_sapiens_DNase_norm_correct.rds")
# human_atac = readRDS("E:\\GlobusDownload/Histone_norm_all_dhs/Homo_sapiens_ATAC_norm_correct.rds")
mouse_dnase = readRDS("E:\\GlobusDownload/Histone_norm_all_dhs/Mus_musculus_DNase_norm_correct.rds")
mouse_atac = readRDS("E:\\GlobusDownload/Histone_norm_all_dhs/Mus_musculus_ATAC_norm_correct.rds")

rm(regions)
gc()
# metadata
meta = read_tsv("E:\\GlobusDownload/Histone_norm_correct/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
meta$`Tissue/cell types` = NA

region_ind_human_dnase = logical()
num_batches = 20
batch_size = ceiling(nrow(human_dnase) / num_batches)
for (j in 1:num_batches) {
        message(j)
        indices = (1:nrow(human_dnase))
        indices_batch = which(ceiling(indices / batch_size) == j)
        region_ind_human_dnase = c(region_ind_human_dnase,
                                   matrixStats::rowMaxs(human_dnase[indices_batch, ]) > 20)
}
region_ind_human_atac = matrixStats::rowMaxs(human_atac) > 20
region_filter_human = which(region_ind_human_dnase | region_ind_human_atac)
region_filter_mouse = which(matrixStats::rowMaxs(mouse_dnase) > 20 |
                                    matrixStats::rowMaxs(mouse_atac) > 20)

calc_var_std_ver2 <- function(mat, f = 1/10, plot_sample_size = 5000) {
        num_batches = 50
        batch_size = ceiling(nrow(mat)/num_batches)
        means = map(1:num_batches, function(j) {
                message(j)
                indices = (1:nrow(mat))
                indices_batch = which(ceiling(indices / batch_size) == j)
                out = rowMeans(mat[indices_batch, ])
        }) %>% purrr::reduce(c)
        vars = map(1:num_batches, function(j) {
                message(j)
                indices = (1:length(means))
                indices_batch = which(ceiling(indices / batch_size) == j)
                out = matrixStats::rowVars(mat[indices_batch, ])
        }) %>% purrr::reduce(c)
        sel = means !=0 & vars != 0
        
        means = log10(means[sel])
        vars = log10(vars[sel])
        lowess_fit = lowess(means, vars, f = f)
        vars_exp = approxfun(lowess_fit)(means)
        # negative vars_exp may be produced
        # mat_std = (mat[sel, ] - 10^means)/sqrt(10^vars_exp)
        # mat_std = pmin(mat_std, sqrt(ncol(mat_std)))
        # vars_std1 = matrixStats::rowVars(mat_std)
        vars_std = map(1:num_batches, function(j) {
                message(j)
                indices = (1:length(means))
                indices_batch = which(ceiling(indices / batch_size) == j)
                mat_std = (mat[which(sel)[indices_batch], ] - 10^means[indices_batch])/sqrt(10^vars_exp[indices_batch])
                mat_std = pmin(mat_std, sqrt(ncol(mat_std)))
                vars_std = matrixStats::rowVars(mat_std)
        }) %>% purrr::reduce(c)
        # sample_idx = sample(nrow(mat), size = plot_sample_size)
        # smoothScatter(means[sample_idx], vars[sample_idx])
        # lines(lowess_fit)
        vars_std_all = numeric(nrow(mat))
        vars_std_all[sel] = vars_std
        vars_std_all[!sel] = NA
        vars_std_all
}
        
data_vars_std = list(calc_var_std_ver2(human_dnase),
                     calc_var_std_ver2(human_atac),
                     calc_var_std_ver2(mouse_dnase),
                     calc_var_std_ver2(mouse_atac))
data_vars_rank = data_vars_std %>% map(rank)

data_vars_std_human_dnase = calc_var_std_ver2(human_dnase)
data_vars_std_human_atac = calc_var_std_ver2(human_atac)
data_human_var_feature = rownames(human_dnase)[order(data_vars_rank[[1]] + data_vars_rank[[2]], decreasing = T)[1:25000]]
data_human_var_feature = data_human_var_feature[data_human_var_feature %in% rownames(human_dnase)[region_filter_human]]

data_vars_std_mouse_dnase = calc_var_std_ver1(mouse_dnase)
data_vars_std_mouse_atac = calc_var_std_ver1(mouse_atac)
data_mouse_var_feature = rownames(mouse_dnase)[order(rank(data_vars_std_mouse_dnase) + rank(data_vars_std_mouse_atac), decreasing = T)[1:25000]]
data_mouse_var_feature = data_mouse_var_feature[data_mouse_var_feature %in% rownames(mouse_dnase)[region_filter_mouse]]

saveRDS(data_human_var_feature, file = "./intermediate_data/temp_encode4_human_chromatin_all_dhs_batch_correct_anchor_features.rds")
saveRDS(data_mouse_var_feature, file = "./intermediate_data/temp_encode4_mouse_chromatin_all_dhs_batch_correct_anchor_features.rds")

data_human_var_feature = readRDS("./intermediate_data/temp_encode4_human_chromatin_all_dhs_batch_correct_anchor_features.rds")
data_mouse_var_feature = readRDS("./intermediate_data/temp_encode4_mouse_chromatin_all_dhs_batch_correct_anchor_features.rds")

all(colnames(human_dnase) %in% meta$Accession)
all(colnames(human_atac) %in% meta$Accession)
all(colnames(mouse_dnase) %in% meta$Accession)
all(colnames(mouse_atac) %in% meta$Accession)

data_seurat_human = map2(list(human_dnase, human_atac), c("Human Dnase", "Human ATAC"), function(x, source) {
        y = CreateSeuratObject(counts = x[data_human_var_feature, ],
                               assay = "OpenChromatin")
        y$source = source
        y
})
data_seurat_human = map(data_seurat_human, function(x) {
        add_metadata(x, meta)
})
names(data_seurat_human) = c("Human Dnase", "Human ATAC")

# DEFINE common tissues to show #
# TODO: to be updated with metadta input, depends on later code now
# ct_assay_tissue = table(human_obj$tissue, human_obj$source)
# ct_assay_tissue = ct_assay_tissue[order(ct_assay_tissue[, 1]), ]
# common_tissue_human = rownames(ct_assay_tissue)[ct_assay_tissue[, 1] > 0]
# common_tissue_human = common_tissue_human[common_tissue_human != "NA"]
# 
# ct_assay_tissue = table(mouse_obj$tissue, mouse_obj$source)
# ct_assay_tissue = ct_assay_tissue[order(ct_assay_tissue[, 1]), ]
# common_tissue_mouse = rownames(ct_assay_tissue)[ct_assay_tissue[, 1] > 0]
# common_tissue_mouse = common_tissue_mouse[common_tissue_mouse != "NA"]
# 
# common_tissue_names = sort(union(common_tissue_mouse, common_tissue_human))
# define consistent color palette
# common_tissue_col = c25[1:23]
# names(common_tissue_col) = common_tissue_names

# plot UMAP and heatmap before batch correction:
human_obj = merge(data_seurat_human[[1]], data_seurat_human[[2]])
human_obj = FindVariableFeatures(human_obj, nfeatures = 20000)
human_obj = ScaleData(human_obj, features = data_human_var_feature)
human_obj = RunPCA(human_obj, features = data_human_var_feature)
human_obj = RunUMAP(human_obj, dims = 1:30)

anchors_human = FindIntegrationAnchors(data_seurat_human,
                                       anchor.features = data_human_var_feature,
                                       reduction = "cca",
                                       k.anchor = 8,
                                       k.filter = 20,
                                       dims = 1:30)

# Integrated DATA
human_integrated = IntegrateData(anchors_human, features.to.integrate = anchors_human@anchor.features,
                                 dims = 1:15,
                                 k.weight = 10)
human_integrated = quick_process(human_integrated, var_feature = data_human_var_feature, dims = 1:30)

# Figure BC_UHA
(UMAPPlot(human_integrated, pt.size = 0.5, group.by = "source") + theme(legend.position = "top",
                                                                        text = element_text(size = 10))) %>%
        push_pdf(., file_name = "human_oc_after_correct", w = 3.2, h = 3.5)
HoverLocator(plot = UMAPPlot(human_integrated, pt.size = 3, group.by = "source"),
             information = FetchData(human_integrated,
                                     vars = c("source", "Biosample term name", "Life stage"))) %>%
        push_widget(file_name = "human_oc_after_correct")

# human_obj$tissue_stage = str_replace_all(paste0(human_obj$`Tissue/cell types`, " ", ifelse(human_obj$`Life stage` != "adult", human_obj$`Life stage`, "")), " ", "")
# ct_assay_tissue = table(human_obj$tissue_stage, human_obj$source)
# ct_assay_tissue = ct_assay_tissue[order(ct_assay_tissue[, 1]), ]
# common_tissue = rownames(ct_assay_tissue)[ct_assay_tissue[, 1] > 0]
# human_obj$tissue_stage[!human_obj$tissue_stage %in% common_tissue] = NA
# 
# human_obj$tissue = human_obj$`Tissue/cell types`
# human_obj$tissue[!human_obj$tissue %in% common_tissue_human] = NA
# 
# human_integrated$tissue = human_integrated$`Tissue/cell types`
# human_integrated$tissue[!human_integrated$tissue %in% common_tissue_human] = NA       
# Figure BC_UHB
(UMAPPlot(human_obj, pt.size = 0.5, group.by = "Assay title") + 
                scale_color_manual(values = assay_col) + theme(legend.position = "top",
                                                               text = element_text(size = 10))) %>%
        push_pdf(., file_name = "human_oc_before_correct", w = 3.2, h = 3.5)
HoverLocator(plot = UMAPPlot(human_obj, pt.size = 1, group.by = "source"),
             information = FetchData(human_obj,
                                     vars = c("source", "Biosample term name", "Life stage"))) %>%
        push_widget(file_name = "human_oc_before_correct")

# (UMAPPlot(human_obj, pt.size = 0.5, group.by = "tissue") + 
#                 scale_color_manual(values = common_tissue_col,
#                                    na.value = "565656") +
#                 theme(legend.position = "top",
#                       text = element_text(size = 10),
#                       legend.text = element_text(size = 6))) %>%
#         push_pdf(., file_name = "human_oc_before_correct_tissue", w = 4, h = 5)
# 
# HoverLocator(plot = UMAPPlot(human_obj, pt.size = 1, group.by = "tissue"),
#              information = FetchData(human_obj,
#                                      vars = c("source", "Biosample term name", "Life stage"))) %>%
#         push_widget(file_name = "human_oc_before_correct_tissue")

# (UMAPPlot(human_integrated, pt.size = 0.5, group.by = "tissue") +
#                 scale_color_manual(values = common_tissue_col,                                    na.value = "565656")) + theme(legend.position = "top",
#                                                                                                                                 text = element_text(size = 10),
#                                                                                                                                 legend.text = element_text(size = 6))) %>%
#         push_pdf(., file_name = "human_oc_after_correct_tissue", w = 4, h = 5)
# HoverLocator(plot = UMAPPlot(human_integrated, pt.size = 3, group.by = "tissue_stage"),
#              information = FetchData(human_integrated,
#                                      vars = c("source", "Biosample term name", "Life stage"))) %>%
#         push_widget(file_name = "human_oc_after_correct_tissue")
# 
# (UMAPPlot(human_integrated, pt.size = 0.5, group.by = "Assay title") + 
#                 scale_color_manual(values = assay_col) +theme(legend.position = "top",
#                                                               text = element_text(size = 10))) %>%
#         push_pdf(., file_name = "human_oc_after_correct", w = 3.2, h = 3.5)

# anchor_tb = make_anchor_tb(anchors_human, data_seurat_human[[1]], data_seurat_human[[2]])
# write_tsv(anchor_tb, "./output/encode4_histone_enrichment/human_chromatin_all_dhs_batch_correct_anchors.tsv")


# make dummy seurat objects and anchor for integration
num_batches = 40
batch_size = ceiling(nrow(human_dnase) / num_batches)
for (j in 2:num_batches) {
        message(j)
        indices = (1:nrow(human_dnase))
        indices_batch = which(ceiling(indices / batch_size) == j)
        dhs_id_batch = rownames(human_dnase)[indices_batch]
        
        dummy_seurat = map2(list(human_dnase, human_atac), c("Human Dnase", "Human ATAC"), function(x, source) {
                y = CreateSeuratObject(counts = x[union(data_human_var_feature, dhs_id_batch), ],
                                       assay = "OpenChromatin")
                y$source = source
                y$Accession = colnames(x)
                y
        })
        dummy_anchors = FindIntegrationAnchors(dummy_seurat,
                                               anchor.features = data_human_var_feature,
                                               reduction = "cca",
                                               k.anchor = 8,
                                               k.filter = 20,
                                               dims = 1:30)
        dummy_integrated = IntegrateData(dummy_anchors,
                                         features.to.integrate = dhs_id_batch,
                                         dims = 1:15,
                                         k.weight = 10)
        out_mat = as.matrix(dummy_integrated@assays$integrated@data)
        out_mat[out_mat < 0] = 0
        saveRDS(out_mat, paste0("./intermediate_data/temp_encodev4_human_chromatin_all_dhs_correct_part_", j, ".rds"))
        gc()
}

j = 1
human_correct = readRDS(paste0("./intermediate_data/temp_encodev4_human_chromatin_all_dhs_correct_part_", j, ".rds"))
for (j in 2:40) {
        message(j)
        x = readRDS(paste0("./intermediate_data/temp_encodev4_human_chromatin_all_dhs_correct_part_", j, ".rds"))
        human_correct = rbind(human_correct, x)
        gc()
}
saveRDS(human_correct, file = "./intermediate_data/encodev4_human_chromatin_all_dhs_correct.rds")
# saveRDS(human_correct, file = "../human_correct_ver1.rds") # ver1 used low var dhs for dnase/atac integration

# Heatmap Comparison
# human_correct = readRDS("./intermediate_data/human_chromatin_correct.rds")
sample_indices = rownames(human_dnase)[order(data_vars_rank[[1]] + data_vars_rank[[2]], decreasing = F)[1:1000]]
# sample_indices = rownames(human_dnase)[order(data_vars_rank[[1]], decreasing = T)[1:1000]]
# sample_indices = sample(which(regions$spearman_global < 0.2 & regions$human_var_std_global > 1.2 & regions$mouse_var_std_global > 1.2), size = 250)

val_col = circlize::colorRamp2(c(-2, 0, 2), colors = c("blue", "white", "red"))
human_norm = cbind(human_dnase[sample_indices, ],
                   human_atac[sample_indices, ])
human_correct_scale = t(scale(t(log2(1+human_correct[sample_indices, ]))))
human_norm_scale = t(scale(t(log2(1+human_norm))))

all(human_obj$Accession == colnames(human_integrated))
# Figure
# plot 1000 variable features
(Heatmap(human_norm_scale, show_column_names = F, show_row_names = F,
         column_title = "Before",
         col = val_col,
         name = "Scaled Count",
         show_row_dend = F,
         show_column_dend = F,
         bottom_annotation = HeatmapAnnotation(Assay = human_obj$`Assay title`,
                                               Tissue = human_obj$tissue,
                                               col = list(Assay = assay_col,
                                                          Tissue = common_tissue_col),
                                               annotation_label = " ",
                                               show_legend = F)) +
                Heatmap(human_correct_scale, show_column_names = F, show_row_names = F,
                        column_title = "After",
                        col = val_col,
                        name = "Scaled Count",
                        show_row_dend = F,
                        show_column_dend = F,
                        bottom_annotation = HeatmapAnnotation(Assay = human_integrated$`Assay title`,
                                                              Tissue = human_integrated$tissue,
                                                              col = list(Assay = assay_col,
                                                                         Tissue = common_tissue_col),
                                                              annotation_label = " ",
                                                              show_legend = F)
                )) %>%
        push_png(file_name = "heatmap_hg38_chromatin_before_after_correct", w = 7.2, h = 3.6, res = 600)

# this verification is deprecatd
# getting matched samples for DNase and ATAC for verification
exp_biosample = purrr::map(exp_df_list, function(x) {
        x = dplyr::select(x, `Biosample type`, `Biosample term name`, `Biosample treatment`, `Life stage`, `Accession`)
        x
})
exp_biosample_known = inner_join(exp_biosample[[1]] %>% 
                                         nest(exp_tb = -c(`Biosample type`, `Biosample term name`, `Biosample treatment`, `Life stage`)),
                                 exp_biosample[[2]] %>% 
                                         nest(exp_tb = -c(`Biosample type`, `Biosample term name`, `Biosample treatment`, `Life stage`)),
                                 by = c("Biosample type", "Biosample term name", "Biosample treatment", "Life stage"))
exp_biosample_known = exp_biosample_known %>% mutate(dnase_acc = map_chr(exp_tb.x, function(x) x$Accession[1]),
                                                     atac_acc = map_chr(exp_tb.y, function(x) x$Accession[1]))
# DNase ATAC known matches

# verifying using dendrogram
exp_indices = c(exp_biosample_known$dnase_acc,
                exp_biosample_known$atac_acc)
cl_biosample = rep(paste0(exp_biosample_known$`Biosample term name`, "_", exp_biosample_known$`Life stage`, "_", exp_biosample_known$`Biosample treatment`), 2)
library(dendextend)
par(oma=c(3,3,3,10), mar = rep(0, 4))

dend_correct = as.dendrogram(hclust(dist(human_integrated@reductions$umap@cell.embeddings[exp_indices, ])))
dend_correct = dend_correct %>% set("labels_col", value = ifelse(labels(dend_correct) %in% exp_biosample_known$dnase_acc, assay_col[1], assay_col[2]))
labels(dend_correct) = cl_biosample[match(labels(dend_correct), exp_indices)]
dend_correct %>% set("labels_cex", 0.5) %>% plot(horiz = T) %>% push_pdf(file_name = "human_dendro_after", w = 6, h = 6, ps = 8)

dend_norm = as.dendrogram(hclust(dist(human_obj@reductions$umap@cell.embeddings[exp_indices, ])))
dend_norm = dend_norm %>% set("labels_col", value = ifelse(labels(dend_norm) %in% exp_biosample_known$dnase_acc, assay_col[1], assay_col[2]))
labels(dend_norm) = cl_biosample[match(labels(dend_norm), exp_indices)]
dend_norm %>% set("labels_cex", 0.5) %>% plot(horiz = T) %>% push_pdf(file_name = "human_dendro_before", w = 6, h = 6, ps = 8)

# cl_correct = kmeans(human_integrated@reductions$umap@cell.embeddings[exp_indices, ], 30)$cluster
# cl_norm = kmeans(human_obj@reductions$umap@cell.embeddings[exp_indices, ], 30)$cluster
# mclust::adjustedRandIndex(cl_correct, cl_biosample)
# end human Heatmap

# Mouse correction new below
data_seurat_mouse = map2(list(mouse_dnase, mouse_atac), c("Mouse Dnase", "Mouse ATAC"), function(x, source) {
        y = CreateSeuratObject(counts = x[data_mouse_var_feature, ],
                               assay = "OpenChromatin")
        y$source = source
        y
})
data_seurat_mouse = map(data_seurat_mouse, function(x) {
        add_metadata(x, meta)
})
names(data_seurat_mouse) = c("Mouse Dnase", "Mouse ATAC")

# plot UMAP and heatmap before batch correction:
# mouse_obj = merge(data_seurat_mouse[[1]], data_seurat_mouse[[2]])
# mouse_obj = FindVariableFeatures(mouse_obj, nfeatures = 20000)
# mouse_obj = ScaleData(mouse_obj, features = data_mouse_var_feature)
# mouse_obj = RunPCA(mouse_obj, features = data_mouse_var_feature)
# mouse_obj = RunUMAP(mouse_obj, dims = 1:30)

anchors_mouse = FindIntegrationAnchors(data_seurat_mouse,
                                       anchor.features = data_mouse_var_feature,
                                       reduction = "cca",
                                       k.anchor = 3,
                                       k.filter = 10,
                                       dims = 1:30)
anchor_tb = make_anchor_tb(anchors_mouse, data_seurat_mouse[[1]], data_seurat_mouse[[2]])
write_tsv(anchor_tb, "./output/encode4_histone_enrichment/mouse_chromatin_all_dhs_batch_correct_anchors.tsv")

# producing results for batches
num_batches = 40
batch_size = ceiling(nrow(mouse_dnase)/num_batches)
for (j in 1:num_batches) {
        message(j)
        indices = (1:nrow(mouse_dnase))
        indices_batch = which(ceiling(indices / batch_size) == j)
        dhs_id_batch = rownames(mouse_dnase)[indices_batch]
        
        dummy_seurat = map2(list(mouse_dnase, mouse_atac), c("Mouse Dnase", "Mouse ATAC"), function(x, source) {
                y = CreateSeuratObject(counts = x[union(data_mouse_var_feature, dhs_id_batch), ],
                                       assay = "OpenChromatin")
                y$source = source
                y$Accession = colnames(x)
                y
        })
        dummy_anchors = FindIntegrationAnchors(dummy_seurat,
                                               anchor.features = data_mouse_var_feature,
                                               reduction = "cca",
                                               k.anchor = 3,
                                               k.filter = 10,
                                               dims = 1:30)
        dummy_integrated = IntegrateData(dummy_anchors, features.to.integrate = dhs_id_batch,
                                         dims = 1:15,
                                         k.weight = 20)
        out_mat = as.matrix(dummy_integrated@assays$integrated@data)
        out_mat[out_mat < 0] = 0
        saveRDS(out_mat, paste0("./intermediate_data/temp_encodev4_all_dhs_mouse_chromatin_correct_part_", j, ".rds"))
}
j = 1
mouse_correct = readRDS(paste0("./intermediate_data/temp_encodev4_mouse_chromatin_correct_part_", j, ".rds"))
for (j in 2:40) {
        message(j)
        x = readRDS(paste0("./intermediate_data/temp_encodev4_mouse_chromatin_correct_part_", j, ".rds"))
        mouse_correct = rbind(mouse_correct, x)
        gc()
}
saveRDS(mouse_correct, file = "./intermediate_data/encodev4_all_dhs_mouse_chromatin_correct.rds")
# end Mosue correction new

# MOUSE corretion below ------------------------------------------------------------------------
data_seurat_mouse = map2(data_norm[3:4], names(data_norm)[3:4], function(x, source) {
        y = CreateSeuratObject(counts = x[data_mouse_var_feature, ],
                               assay = "OpenChromatin")
        y$source = source
        y$Accession = colnames(x)
        y
})
data_seurat_mouse = map2(data_seurat_mouse, exp_df_list[3:4], function(x, exp_df) {
        add_metadata(x, exp_df)
})
names(data_seurat_mouse) = names(data_norm)[3:4]

mouse_obj = merge(data_seurat_mouse[[1]], data_seurat_mouse[[2]])
mouse_obj = FindVariableFeatures(mouse_obj, nfeatures = 20000)
mouse_obj = ScaleData(mouse_obj, features = data_human_var_feature)
mouse_obj = RunPCA(mouse_obj, features = data_human_var_feature)
mouse_obj = RunUMAP(mouse_obj, dims = 1:30)

mouse_obj$tissue = mouse_obj$`Tissue/cell types`
mouse_obj$tissue[!mouse_obj$tissue %in% common_tissue_human] = NA

mouse_integrated$tissue = mouse_integrated$`Tissue/cell types`
mouse_integrated$tissue[!mouse_integrated$tissue %in% common_tissue_human] = NA       

# Figure BC_UMB
(UMAPPlot(mouse_obj, pt.size = 0.5, group.by = "Assay title") +
                scale_color_manual(values = assay_col) +
                theme(legend.position = "top",
                      text = element_text(size = 10))) %>%
        push_pdf(., file_name = "mouse_oc_before_correct", w = 3.2, h = 3.5)
HoverLocator(plot = UMAPPlot(mouse_obj, pt.size = 1, group.by = "source"),
             information = FetchData(mouse_obj,
                                     vars = c("source", "Biosample term name", "Life stage"))) %>%
        push_widget(file_name = "mouse_oc_before_correct")

# Integrated
anchors_mouse = FindIntegrationAnchors(data_seurat_mouse,
                                       anchor.features = data_mouse_var_feature,
                                       reduction = "cca",
                                       k.anchor = 3,
                                       k.filter = 10,
                                       dims = 1:30)
# repeat mouse new pipeline for integrating across assay
mouse_integrated = IntegrateData(anchors_mouse, features.to.integrate = data_mouse_var_feature,
                                 dims = 1:15,
                                 k.weight = 20)
mouse_integrated = quick_process(mouse_integrated, var_feature = data_mouse_var_feature, dims = 1:15)

mouse_obj = fix_metadata(mouse_obj)
mouse_integrated = fix_metadata(mouse_integrated)

(UMAPPlot(mouse_integrated, pt.size = 0.5, group.by = "Assay title") +
                scale_color_manual(values = assay_col) +
                theme(legend.position = "top",
                      text = element_text(size = 10))) %>%
        push_pdf(., file_name = "mouse_oc_after_correct", w = 3.2, h = 3.5)
HoverLocator(plot = UMAPPlot(mouse_integrated, pt.size = 3, group.by = "source"),
             information = FetchData(mouse_integrated,
                                     vars = c("source", "Biosample term name", "Life stage"))) %>% 
        push_widget(file_name = "mouse_oc_after_correct")

(UMAPPlot(mouse_obj, pt.size = 0.5, group.by = "tissue") + 
                scale_color_manual(values = common_tissue_col,                                    na.value = "565656")) +
        theme(legend.position = "top",
              text = element_text(size = 10),
              legend.text = element_text(size = 6))) %>%
        push_pdf(., file_name = "mouse_oc_before_correct_tissue", w = 4, h = 5)

HoverLocator(plot = UMAPPlot(mouse_obj, pt.size = 1, group.by = "tissue"),
             information = FetchData(mouse_obj,
                                     vars = c("source", "Biosample term name", "Life stage"))) %>%
        push_widget(file_name = "mouse_oc_before_correct_tissue")

(UMAPPlot(mouse_integrated, pt.size = 0.5, group.by = "tissue") + 
                scale_color_manual(values = common_tissue_col,                                    na.value = "565656")) + theme(legend.position = "top",
                                                                                                                                text = element_text(size = 10),
                                                                                                                                legend.text = element_text(size = 6))) %>%
        push_pdf(., file_name = "mouse_oc_after_correct_tissue", w = 4, h = 5)
HoverLocator(plot = UMAPPlot(mouse_integrated, pt.size = 3, group.by = "tissue"),
             information = FetchData(mouse_integrated,
                                     vars = c("source", "Biosample term name", "Life stage"))) %>%
        push_widget(file_name = "mouse_oc_after_correct_tissue")

# make dummy seurat objects and anchor for integration
num_batches = 20
batch_size = ceiling(nrow(human_dnase)/num_batches)
int_mat_list = map(1:num_batches, function(j) {
        message(j)
        indices = (1:nrow(human_dnase))
        indices_batch = which(ceiling(indices / batch_size) == j)
        dhs_id_batch = rownames(human_dnase)[indices_batch]
        
        dummy_seurat = map2(data_norm[3:4], names(data_norm)[3:4], function(x, source) {
                y = CreateSeuratObject(counts = x[union(anchors_human@anchor.features, dhs_id_batch), ],
                                       assay = "OpenChromatin")
                y$source = source
                y$Accession = colnames(x)
                y
        })
        dummy_anchors = FindIntegrationAnchors(dummy_seurat,
                                               anchor.features = data_mouse_var_feature,
                                               reduction = "cca",
                                               k.anchor = 3,
                                               k.filter = 10,
                                               dims = 1:30)
        dummy_integrated = IntegrateData(dummy_anchors, features.to.integrate = dhs_id_batch,
                                         dims = 1:15,
                                         k.weight = 20)
        out_mat = as.matrix(dummy_integrated@assays$integrated@data)
        out_mat[out_mat < 0] = 0
        out_mat
})
mouse_correct = reduce(int_mat_list, rbind)
saveRDS(mouse_correct, file = "./intermediate_data/mouse_chromatin_correct.rds")
rm(int_mat_list); gc()

sample_indices = rownames(human_dnase)[order(data_vars_rank[[3]] + data_vars_rank[[4]], decreasing = F)[1:1000]]
val_col = circlize::colorRamp2(c(-1.5, 0, 1.5), colors = c("blue", "white", "red"))
mouse_norm = cbind(mouse_dnase[sample_indices, ],
                   mouse_atac[sample_indices, ])
mouse_correct_scale = t(scale(t(log2(1+mouse_correct[sample_indices, ]))))
mouse_norm_scale = t(scale(t(log2(1+mouse_norm))))

# Figure
# plot first 1000 features

(Heatmap(mouse_norm_scale, show_column_names = F, show_row_names = F,
         column_title = "Before",
         col = val_col,
         name = "Scaled Count",
         show_row_dend = F,
         show_column_dend = F,
         bottom_annotation = HeatmapAnnotation(Assay = mouse_obj$`Assay title`,
                                               Tissue = mouse_obj$tissue,
                                               col = list(Assay = assay_col,
                                                          Tissue = common_tissue_col),
                                               annotation_label = " ",
                                               show_legend = F)) +
                (Heatmap(mouse_correct_scale, show_column_names = F, show_row_names = F,
                         column_title = "After",
                         col = val_col,
                         name = "Scaled Count",
                         show_row_dend = F,
                         show_column_dend = F,
                         bottom_annotation = HeatmapAnnotation(Assay = mouse_integrated$`Assay title`,
                                                               Tissue = mouse_integrated$tissue,
                                                               col = list(Assay = assay_col,
                                                                          Tissue = common_tissue_col),
                                                               annotation_label = " ",
                                                               show_legend = F)
                ))) %>%
        push_png(file_name = "heatmap_mm10_chromatin_before_after_correct", w = 7.2, h = 3.6, res = 600)

# getting matched samples for DNase and ATAC for verification
# exp_df_list$mouse_dnase = filter(exp_df_list$mouse_dnase, Accession %in% colnames(data_seurat_mouse$mouse_dnase))
# exp_df_list$mouse_atac = filter(exp_df_list$mouse_atac, Accession %in% colnames(data_seurat_mouse$mouse_atac))
# 
# exp_biosample = purrr::map(exp_df_list, function(x) {
#         x = dplyr::select(x, `Biosample type`, `Biosample term name`, `Biosample treatment`, `Life stage`, `Accession`)
#         x
# })
# exp_biosample_known = inner_join(exp_biosample[[3]] %>% 
#                                          nest(exp_tb = -c(`Biosample type`, `Biosample term name`, `Biosample treatment`, `Life stage`)),
#                                  exp_biosample[[4]] %>% 
#                                          nest(exp_tb = -c(`Biosample type`, `Biosample term name`, `Biosample treatment`, `Life stage`)),
#                                  by = c("Biosample type", "Biosample term name", "Biosample treatment", "Life stage"))
# exp_biosample_known = exp_biosample_known %>% mutate(dnase_acc = map_chr(exp_tb.x, function(x) x$Accession[1]),
#                                                      atac_acc = map_chr(exp_tb.y, function(x) x$Accession[1]))
# DNase ATAC known matches

# verifying using dendrogram
# exp_indices = c(exp_biosample_known$dnase_acc,
#                 exp_biosample_known$atac_acc)
# cl_biosample = rep(paste0(exp_biosample_known$`Biosample term name`, "_", exp_biosample_known$`Life stage`, "_", exp_biosample_known$`Biosample treatment`), 2)
# library(dendextend)
# par(oma=c(3,3,3,10), mar = rep(0, 4))
# 
# dend_correct = as.dendrogram(hclust(dist(mouse_integrated@reductions$umap@cell.embeddings[exp_indices, ])))
# dend_correct = dend_correct %>% set("labels_col", value = ifelse(labels(dend_correct) %in% exp_biosample_known$dnase_acc, assay_col[1], assay_col[2]))
# labels(dend_correct) = cl_biosample[match(labels(dend_correct), exp_indices)]
# dend_correct %>% set("labels_cex", 0.5) %>% plot(horiz = T) %>% push_pdf(file_name = "mouse_dendro_after", w = 6, h = 6, ps = 8)
# 
# dend_norm = as.dendrogram(hclust(dist(mouse_obj@reductions$umap@cell.embeddings[exp_indices, ])))
# dend_norm = dend_norm %>% set("labels_col", value = ifelse(labels(dend_norm) %in% exp_biosample_known$dnase_acc, assay_col[1], assay_col[2]))
# labels(dend_norm) = cl_biosample[match(labels(dend_norm), exp_indices)]
# dend_norm %>% set("labels_cex", 0.5) %>% plot(horiz = T) %>% push_pdf(file_name = "mouse_dendro_before", w = 6, h = 6, ps = 8)
# 


