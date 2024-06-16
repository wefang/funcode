library(tidyverse)
library(Seurat)
source("./helper/helper.R")
file_meta = readr::read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_file_metadata.tsv")
exp_meta = readr::read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
exp_meta = exp_meta[exp_meta$`Target of assay` %in% c("H3K4me1", "H3K4me3", "H3K27ac") &
                            exp_meta$`Assay name` != "Mint-ChIP-seq", ]
mouse_correct = readRDS("./intermediate_data/encodev4_mouse_chromatin_correct.rds")
mouse_dnase_var_std = calc_var_std_ver1(mouse_correct)
# set assay
assay = "H3K4me1"
# assay = "H3K27ac"
# assay = "H3K4me3"
histone_mat = readRDS(paste0("E:\\GlobusDownload/histone_norm_correct/Mus_musculus_", assay, "_norm_correct.rds"))
mouse_var_std = calc_var_std_ver1(histone_mat)
exp_meta_raw = readr::read_tsv("C:\\Projects/scdl_v0/metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
exp_meta_raw$`Tissue/cell types` = NA
histone_mat = histone_mat[, exp_meta_raw$`Assay title`[match(colnames(histone_mat), exp_meta_raw$Accession)] == "Histone ChIP-seq"]
all(colnames(histone_mat) %in% exp_meta$Accession)
exp_meta_sel = exp_meta[match(colnames(histone_mat), exp_meta$Accession), ]

region_filter = matrixStats::rowMaxs(histone_mat) > 5 &
        matrixStats::rowMaxs(mouse_correct) > 10
mouse_mod_var_feature = rownames(mouse_correct)[which(rank(mouse_var_std) > 1.6e6 & rank(mouse_dnase_var_std) > 1.6e6 & region_filter)]
length(mouse_mod_var_feature)

rownames(histone_mat) = rownames(mouse_correct)
mouse_correct_obj = CreateSeuratObject(counts = histone_mat[mouse_mod_var_feature, ],
                                       assay = "OpenChromatin")
mouse_correct_obj$organism = "mouse"
mouse_correct_obj$modality = assay
exp_meta_sel$`Tissue/cell types` = NA
mouse_correct_obj = add_metadata(mouse_correct_obj, exp_meta_sel)

mouse_dnase_correct_obj = CreateSeuratObject(counts = mouse_correct[mouse_mod_var_feature, ],
                                             assay = "OpenChromatin")
mouse_dnase_correct_obj = add_metadata(mouse_dnase_correct_obj, exp_meta_raw)
mouse_dnase_correct_obj$modality = "DNase/ATAC"
mouse_dnase_correct_obj$organism = "mouse"
# end object definition -------------------------------------------------------

# reverse transfer here to get imputed data
anchors <- FindTransferAnchors(reference = mouse_dnase_correct_obj,
                               query = mouse_correct_obj,
                               features = mouse_mod_var_feature,
                               reduction = "cca",
                               k.anchor = 5,
                               k.filter = 20,
                               k.score = min(40, dim(mouse_correct_obj)[2]-1),
                               dims = 1:30)
anchor_tb = make_anchor_tb(anchors, mouse_dnase_correct_obj, mouse_correct_obj)
cbind(anchor_tb$`Biosample term name_1`, anchor_tb$`Biosample term name_2`)
write_csv(anchor_tb, file = paste0("./intermediate_data/histone_fc_v1/mouse_", assay, "_chromatin_integration_anchors.csv"))

data_cross_var_feature = readRDS("./intermediate_data/ENCODE4_human_mouse_chromatin_var_features.rds")
# Below are imputed data
dnase_im = TransferData(anchors,
                        refdata = mouse_correct[data_cross_var_feature, ],
                        weight.reduction = "cca",
                        slot = "counts",
                        k.weight = 10)
dnase_im_mat = as.matrix(dnase_im@counts)
saveRDS(dnase_im_mat, paste0("./intermediate_data/histone_rev_transfer_v1/", "mouse_rev_transfer_", assay, ".rds"))


