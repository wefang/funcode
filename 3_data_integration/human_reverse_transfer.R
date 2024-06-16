library(Seurat)
library(tidyverse)
source("./helper/helper.R")
file_meta = readr::read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_file_metadata.tsv")
exp_meta = readr::read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
exp_meta = exp_meta[exp_meta$`Target of assay` %in% c("H3K4me1", "H3K4me3", "H3K27ac") &
                            exp_meta$`Assay name` != "Mint-ChIP-seq", ]
exp_meta$`Tissue/cell types` = NA
human_correct = readRDS("./intermediate_data/encodev4_human_chromatin_correct.rds")
human_dnase_var_std = calc_var_std_ver1(human_correct)
all(exp_meta$Accession %in% exp_meta_wc$Accession)
# set assay name
# assay = "H3K4me3"
# assay = "H3K27ac"
# assay = "H3K4me1"
# Data processed from '2_process_fragments'
histone_mat = readRDS(paste0("E:\\GlobusDownload/histone_norm_correct/Homo_sapiens_", assay, "_norm_correct.rds"))
human_var_std = calc_var_std_ver1(histone_mat)
# control_mat_all = readRDS("E:\\GlobusDownload/histone_norm_correct/Homo_sapiens_Control_norm_correct.rds")
exp_meta_raw = readr::read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
exp_meta_raw$`Tissue/cell types` = NA

histone_mat = histone_mat[, exp_meta_raw$`Assay title`[match(colnames(histone_mat), exp_meta_raw$Accession)] == "Histone ChIP-seq"]
all(colnames(histone_mat) %in% exp_meta$Accession)
exp_meta_sel = exp_meta[match(colnames(histone_mat), exp_meta$Accession), ]
all(colnames(human_correct) %in% exp_meta_raw$Accession)
region_filter = matrixStats::rowMaxs(histone_mat) > 5 &
        matrixStats::rowMaxs(human_correct) > 10
human_mod_var_feature = rownames(human_correct)[which(rank(human_var_std) > 1.6e6 & rank(human_dnase_var_std) > 1.6e6 & region_filter)]
rownames(histone_mat) = rownames(human_correct)

human_correct_obj = CreateSeuratObject(counts = histone_mat[human_mod_var_feature, ],
                                       assay = "OpenChromatin")
human_correct_obj$organism = "human"
human_correct_obj$modality = assay
human_correct_obj = add_metadata(human_correct_obj, exp_meta_sel)
human_dnase_correct_obj = CreateSeuratObject(counts = human_correct[human_mod_var_feature, ],
                                             assay = "OpenChromatin")
human_dnase_correct_obj = add_metadata(human_dnase_correct_obj, exp_meta_raw)
human_dnase_correct_obj$modality = "DNase/ATAC"
human_dnase_correct_obj$organism = "human"
# end object definition -------------------------------------------------------

# reverse transfer here to get imputed data
anchors <- FindTransferAnchors(reference = human_dnase_correct_obj,
                               query = human_correct_obj,
                               features = human_mod_var_feature,
                               reduction = "cca",
                               k.anchor = 8,
                               k.filter = 20,
                               k.score = min(40, dim(human_correct_obj)[2]-1),
                               dims = 1:30)
anchor_tb = make_anchor_tb(anchors, human_dnase_correct_obj, human_correct_obj)
cbind(anchor_tb$`Biosample term name_1`, anchor_tb$`Biosample term name_2`)
write_csv(anchor_tb, file = paste0("./intermediate_data/histone_fc_v1/", assay, "_chromatin_integration_anchors.csv"))

data_cross_var_feature = readRDS("./intermediate_data/ENCODE4_human_mouse_chromatin_var_features.rds")
# Below are imputed data
dnase_im = TransferData(anchors,
                        refdata = human_correct[data_cross_var_feature, ],
                        weight.reduction = "cca",
                        slot = "counts",
                        k.weight = 10)
dnase_im_mat = as.matrix(dnase_im@counts)
saveRDS(dnase_im_mat, paste0("./intermediate_data/histone_rev_transfer_v1/", "human_rev_transfer_", assay, ".rds"))


