library(Seurat)
library(tidyverse)
source("./helper/helper.R")
script_output_dir = "./output/histone_conservation_v1/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
script_plot_dir = "./plots/histone_chromatin_v1/"
if (!dir.exists(script_plot_dir)) {
        dir.create(script_plot_dir)
}

human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_correct.rds")
mouse_correct = readRDS("./intermediate_data/encodev4_mouse_chromatin_correct.rds")
data_cross_var_feature = readRDS("./intermediate_data/ENCODE4_human_mouse_chromatin_var_features.rds")
# human_mat_all = do.call(cbind, c(list(human_correct[data_cross_var_feature, ]), map(c("H3K4me1", "H3K4me3", "H3K27ac"), function(assay) {
#         readRDS(paste0("./intermediate_data/histone_rev_transfer/", "human_rev_transfer_", assay, ".rds"))
# })))
# mouse_mat_all = do.call(cbind, c(list(mouse_correct[data_cross_var_feature, ]), map(c("H3K4me1", "H3K4me3", "H3K27ac"), function(assay) {
#         readRDS(paste0("./intermediate_data/histone_rev_transfer/", "mouse_rev_transfer_", assay, ".rds"))
# })))
# Set assay name
# assay = "H3K4me1"
# assay = "H3K27ac"
# assay = "H3K4me3"
# assay = "methyl"
human_im_mat = readRDS(paste0("./intermediate_data/histone_rev_transfer_v1/", "human_rev_transfer_", assay, ".rds"))
mouse_im_mat = readRDS(paste0("./intermediate_data/histone_rev_transfer_v1/", "mouse_rev_transfer_", assay, ".rds"))
human_mat_all = cbind(human_correct[data_cross_var_feature, ], human_im_mat)
mouse_mat_all = cbind(mouse_correct[data_cross_var_feature, ], mouse_im_mat)
exp_meta = readr::read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
exp_meta = exp_meta[exp_meta$`Target of assay` %in% c("H3K4me1", "H3K4me3", "H3K27ac") &
                            exp_meta$`Assay name` != "Mint-ChIP-seq", ]
exp_meta$`Tissue/cell types` = NA
human_correct_obj = CreateSeuratObject(counts = human_mat_all,
                                       assay = "OpenChromatin")
mouse_correct_obj = CreateSeuratObject(counts = mouse_mat_all,
                                       assay = "OpenChromatin")
human_correct_obj = add_metadata(human_correct_obj, exp_meta)
mouse_correct_obj = add_metadata(mouse_correct_obj, exp_meta)
# removing over-represented sampels for human H3K27ac
if (assay == "H3K27ac") {
        human_correct_obj = human_correct_obj[, !human_correct_obj$`Biosample term name` %in% c("skin epidermis", "middle frontal area 46")]
}
human_correct_obj = FindVariableFeatures(human_correct_obj, nfeatures = 100000)
mouse_correct_obj = FindVariableFeatures(mouse_correct_obj, nfeatures = 100000)
human_correct_obj$organism = "human"
mouse_correct_obj$organism = "mouse"
anchors <-
        FindIntegrationAnchors(
                list(human = human_correct_obj,
                     mouse = mouse_correct_obj),
                anchor.features = 10000,
                k.anchor = 5,
                k.filter = 20,
                k.score = 30,
                dims = 1:30)
gc()
integrated = IntegrateData(anchors,
                           dims = 1:15,
                           k.weight = 20)
integrated$source = paste0(integrated$organism, "_", integrated$`Assay title`)
integrated = quick_process(integrated, anchors@anchor.features, dims = 1:15)
umap_integrated = integrated@reductions$umap@cell.embeddings
saveRDS(umap_integrated,
        file = "./intermediate_data/CA_H3K4me1_UMAP_embed.rds")

(UMAPPlot(integrated, pt.size = 0.5, group.by = "organism", split.by = "Assay title") +
                # scale_color_manual(values = chromatin_col) +
                theme(text = element_text(size = 10), legend.position = "top", legend.text = element_text(size = 8))) %>%
        push_pdf(., file_name = paste0("umap_integrated_split", assay), w = 8.5, h = 3.5, ps = 10)
}

(UMAPPlot(integrated, pt.size = 0.5, group.by = "source") +
                # scale_color_manual(values = chromatin_col) +
                theme(text = element_text(size = 10), legend.position = "top", legend.text = element_text(size = 8))) %>%
        push_pdf(., file_name = paste0("umap_integrated_", assay), w = 3.2, h = 3.5, ps = 10)

(UMAPPlot(integrated[, integrated$`Assay title` == "WGBS"], pt.size = 0.5, group.by = "source") +
                # scale_color_manual(values = chromatin_col) +
                theme(text = element_text(size = 10), legend.position = "top", legend.text = element_text(size = 8))) %>%
        push_pdf(., file_name = paste0("umap_integrated_", assay, "_only"), w = 3.2, h = 3.5, ps = 10)


HoverLocator(plot = UMAPPlot(integrated, pt.size = 3, group.by = "source"), # + 
             # scale_color_manual(values = chromatin_col),
             information = FetchData(integrated,
                                     vars = c("source", "Biosample term name", "Life stage"))) %>% 
        push_widget(file_name = paste0("integrated_", assay))

integrated_list = SplitObject(integrated[, integrated$`Target of assay` == assay], split.by = "organism")
assay_anchors = FindIntegrationAnchors(integrated_list,
                                       assay = c("integrated", "integrated"),
                                       anchor.features = anchors@anchor.features,
                                       reduction = "rpca",
                                       k.anchor = 3,
                                       k.filter = 25,
                                       k.score = 30,
                                       dims = 1:15)
assay_anchor_tb = make_anchor_tb(assay_anchors, data_obj1 = integrated_list[[1]], 
                                 data_obj2 = integrated_list[[2]])
assay_anchor_tb_sel = assay_anchor_tb
b = cbind(assay_anchor_tb_sel$`Biosample term name_1`, assay_anchor_tb_sel$`Biosample term name_2`)
b[order(b[, 2]), ]
saveRDS(assay_anchor_tb, paste0(script_output_dir, "anchors_", assay, "_v1.rds"))
write_tsv(assay_anchor_tb, file = paste0(script_output_dir, "anchors_", assay, "_v1.tsv"))