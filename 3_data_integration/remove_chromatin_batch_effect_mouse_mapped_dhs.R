library(tidyverse)
source("./R_production/config.R")
source("./R_production/def_color.R")
source("./R/helper.R")

script_plot_dir = "./plots/batch_correct_encodev4/"
if (!dir.exists(script_plot_dir)) {
        dir.create(script_plot_dir)
}

# human_dnase_ext = readRDS("E:\\GlobusDownload/mouse_human_mapped_sep12/Homo_sapiens_DNase_norm_correct.rds")
# human_atac_ext = readRDS("E:\\GlobusDownload/mouse_human_mapped_sep12/Homo_sapiens_ATAC_norm_correct.rds")
mouse_dnase_ext = readRDS("E:\\GlobusDownload/mouse_human_mapped_sep12/Mus_musculus_DNase_norm_correct.rds")
mouse_atac_ext = readRDS("E:\\GlobusDownload/mouse_human_mapped_sep12/Mus_musculus_ATAC_norm_correct.rds")

# human_dnase = readRDS("E:\\GlobusDownload/Histone_norm_all_dhs/Homo_sapiens_DNase_norm_correct.rds")
# human_atac = readRDS("E:\\GlobusDownload/Histone_norm_all_dhs/Homo_sapiens_ATAC_norm_correct.rds")
mouse_dnase = readRDS("E:\\GlobusDownload/Histone_norm_all_dhs/Mus_musculus_DNase_norm_correct.rds")
mouse_atac = readRDS("E:\\GlobusDownload/Histone_norm_all_dhs/Mus_musculus_ATAC_norm_correct.rds")

# all(colnames(human_dnase) == colnames(human_dnase_ext))
# all(colnames(human_atac) == colnames(human_atac_ext))

all(colnames(mouse_dnase) == colnames(mouse_dnase_ext))
all(colnames(mouse_atac) == colnames(mouse_atac_ext))

print(load("./intermediate_data/mm10_hg38_regions.rda"))
# rownames(human_dnase_ext) = rownames(human_atac_ext) = paste0("MouseDHS-", hg38_regions$identifier)
rownames(mouse_dnase_ext) = rownames(mouse_atac_ext) = paste0("MouseDHS-", mm10_regions$identifier)

# data_human_var_feature = readRDS("./intermediate_data/temp_encode4_human_chromatin_all_dhs_batch_correct_anchor_features.rds")
data_mouse_var_feature = readRDS("./intermediate_data/temp_encode4_mouse_chromatin_all_dhs_batch_correct_anchor_features.rds")

library(Seurat)
num_batches = 20
batch_size = ceiling(nrow(human_dnase_ext) / num_batches)
for (j in 7:num_batches) {
        message(j)
        indices = (1:nrow(human_dnase_ext))
        indices_batch = which(ceiling(indices / batch_size) == j)
        dhs_id_batch = rownames(human_dnase_ext)[indices_batch]
        
        dummy_seurat = list()
        data_mat = rbind(human_dnase[data_human_var_feature, ],
                         human_dnase_ext[dhs_id_batch, ])
        x = CreateSeuratObject(counts = data_mat,
                               assay = "OpenChromatin")
        x$source = "Human DNase"
        x$Accession = colnames(data_mat)
        dummy_seurat[[1]] = x
        
        data_mat = rbind(human_atac[data_human_var_feature, ],
                         human_atac_ext[dhs_id_batch, ])
        y = CreateSeuratObject(counts = data_mat,
                               assay = "OpenChromatin")
        y$source = "Human ATAC"
        y$Accession = colnames(data_mat)
        dummy_seurat[[2]] = y
                        
        dummy_anchors = FindIntegrationAnchors(dummy_seurat,
                                               anchor.features = data_human_var_feature,
                                               reduction = "cca",
                                               k.anchor = 8,
                                               k.filter = 20,
                                               dims = 1:30)
        dummy_anchors@anchor.features = dhs_id_batch
        dummy_integrated = CustomIntegrateReference(dummy_anchors,
                                                    features = data_human_var_feature,
                                                    features.to.integrate = dhs_id_batch,
                                                    dims = 1:15,
                                                    k.weight = 10)
        out_mat = as.matrix(dummy_integrated)
        out_mat[out_mat < 0] = 0
        saveRDS(out_mat, paste0("./intermediate_data/temp_encodev4_human_chromatin_mouse_mapped_dhs_correct_part_", j, ".rds"))
        gc()
}

j = 1
human_correct = readRDS(paste0("./intermediate_data/temp_encodev4_human_chromatin_mouse_mapped_dhs_correct_part_", j, ".rds"))
for (j in 2:20) {
        message(j)
        x = readRDS(paste0("./intermediate_data/temp_encodev4_human_chromatin_mouse_mapped_dhs_correct_part_", j, ".rds"))
        human_correct = rbind(human_correct, x)
        gc()
}
saveRDS(human_correct, file = "./intermediate_data/encodev4_human_chromatin_mouse_mapped_dhs_correct.rds")

library(Seurat)
num_batches = 20
batch_size = ceiling(nrow(mouse_dnase_ext) / num_batches)
for (j in 1:num_batches) {
        message(j)
        indices = (1:nrow(mouse_dnase_ext))
        indices_batch = which(ceiling(indices / batch_size) == j)
        dhs_id_batch = rownames(mouse_dnase_ext)[indices_batch]
        
        dummy_seurat = list()
        data_mat = rbind(mouse_dnase[data_mouse_var_feature, ],
                         mouse_dnase_ext[dhs_id_batch, ])
        x = CreateSeuratObject(counts = data_mat,
                               assay = "OpenChromatin")
        x$source = "Mouse DNase"
        x$Accession = colnames(data_mat)
        dummy_seurat[[1]] = x
        
        data_mat = rbind(mouse_atac[data_mouse_var_feature, ],
                         mouse_atac_ext[dhs_id_batch, ])
        y = CreateSeuratObject(counts = data_mat,
                               assay = "OpenChromatin")
        y$source = "Mouse ATAC"
        y$Accession = colnames(data_mat)
        dummy_seurat[[2]] = y
        
        dummy_anchors = FindIntegrationAnchors(dummy_seurat,
                                               anchor.features = data_mouse_var_feature,
                                               reduction = "cca",
                                               k.anchor = 8,
                                               k.filter = 20,
                                               dims = 1:30)
        dummy_anchors@anchor.features = dhs_id_batch
        dummy_integrated = CustomIntegrateReference(dummy_anchors,
                                                    features.to.integrate = dhs_id_batch,
                                                    dims = 1:15,
                                                    k.weight = 10)
        out_mat = as.matrix(dummy_integrated)
        out_mat[out_mat < 0] = 0
        saveRDS(out_mat, paste0("./intermediate_data/temp_encodev4_mouse_chromatin_mouse_mapped_dhs_correct_part_", j, ".rds"))
        gc()
}

j = 1
mouse_correct = readRDS(paste0("./intermediate_data/temp_encodev4_mouse_chromatin_mouse_mapped_dhs_correct_part_", j, ".rds"))
for (j in 2:20) {
        message(j)
        x = readRDS(paste0("./intermediate_data/temp_encodev4_mouse_chromatin_mouse_mapped_dhs_correct_part_", j, ".rds"))
        mouse_correct = rbind(mouse_correct, x)
        gc()
}
saveRDS(mouse_correct, file = "./intermediate_data/encodev4_mouse_chromatin_mouse_mapped_dhs_correct.rds")

