library(tidyverse)
library(ComplexHeatmap)
library(GenomicRanges)
source("./helper/helper.R")

# NMF applied to TFBS motif mapping
library(furrr)
plan(multisession, workers = 6)
set.seed(73)
# read raw motif count results
print(load("./metadata_processed/unalignable_motif/motifs.RData"))
# dhs_sample_indices = sample(3584137, size = 1e5)
dhs_sample_indices = 1:3584137
human_motif_ind = do.call(cbind, future_map(1:736, function(k) {
        load(paste0("./metadata_processed/unalignable_motif/target_overlap_",k,"_h_.RData"))
        target_overlap[dhs_sample_indices]
}, .progress = T))
colnames(human_motif_ind) = motifs
saveRDS(human_motif_ind, file = "./metadata_processed/unalignable_motif/human_motif_ind.rds")
mouse_motif_ind = do.call(cbind, future_map(1:736, function(k) {
        load(paste0("./metadata_processed/unalignable_motif/target_overlap_",k,"_m_.RData"))
        target_overlap
}, .progress = T))
colnames(mouse_motif_ind) = motifs
saveRDS(mouse_motif_ind, file = "./metadata_processed/unalignable_motif/mouse_motif_ind.rds")

# remove regions with no motifs
human_motif_ind = human_motif_ind[rowSums(human_motif_ind) > 1, ]
human_motif_ind = human_motif_ind[, colSums(human_motif_ind) > 10]
# Run NMF and save results
library(NMF)
human_motif_subset = human_motif_ind[1:10000, ]
motif_sel = colSums(human_motif_subset) >= 2
human_motif_subset = human_motif_subset[, motif_sel]
res = nmf(t(human_motif_subset), rank = 25)
saveRDS(res, file = "./intermediate_data/motif/motif_nmf_rev.rds")
rownames(mouse_motif_red) = paste0("MouseDHS-", rownames(mouse_motif_red))
# NMF identifies V^T = H^T %*% W^T.
# Because V is sparse, we have W = DW' where D is diagonal and W' is approximately orthogoal
# We apply the below transformation instead: H^T = V^T %*% W to approximate H
# Process in batch
human_motif_ind = readRDS("./metadata_processed/unalignable_motif/human_motif_ind.rds")
mouse_motif_ind = readRDS("./metadata_processed/unalignable_motif/mouse_motif_ind.rds")
batch_size <- 100000
num_rows <- nrow(human_motif_ind)
temp_list <- list()
for (i in seq(1, num_rows, by = batch_size)) {
  end <- min(i + batch_size - 1, num_rows)
  batch <- human_motif_ind[i:end, ]
  batch_result <- batch[, rownames(res@fit@W)] %*% res@fit@W
  temp_list[[length(temp_list) + 1]] <- batch_result
}
human_motif_red <- do.call(rbind, temp_list)
saveRDS(human_motif_red, "metadata_processed/unalignable_motif/human_reduced.rds")
batch_size <- 100000 
num_rows <- nrow(mouse_motif_ind)
temp_list <- list()
for (i in seq(1, num_rows, by = batch_size)) {
  end <- min(i + batch_size - 1, num_rows)
  batch <- mouse_motif_ind[i:end, ]
  batch_result <-  batch[, rownames(res@fit@W)] %*% res@fit@W
  temp_list[[length(temp_list) + 1]] <- batch_result
}
mouse_motif_red <- do.call(rbind, temp_list)
saveRDS(mouse_motif_red, "metadata_processed/unalignable_motif/mouse_reduced.rds")

# Figure S5e, f: ploting motif reduction
script_plot_dir = "./plots/motif_red/"
set.seed(37)
human_motif_red_plot = human_motif_red[dhs_sample_indices[1:500], ]
human_motif_red_scaled = t(scale(t(human_motif_red_plot)))
dhs_subset = which(!apply(human_motif_red_scaled, 1, anyNA))
dhs_hcl = hclust(dist(human_motif_red_plot[dhs_subset, ]))
h_h = Heatmap(human_motif_red_plot[dhs_subset, ],
              show_row_dend = F, show_column_dend = F, show_row_names = F,
              cluster_rows = dhs_hcl,
              cluster_columns = F,
              heatmap_legend_param = list(
                      title = "Coeficient", 
                      title_position = "topcenter",
                      direction = "horizontal",
                      legend_position = "top"
              ))
red_hcl = hclust(dist(res@fit@W))
h_w = Heatmap(t(res@fit@W),
              cluster_columns = red_hcl,
              show_column_names = F,
              show_row_dend = F, show_column_dend = F,
              heatmap_legend_param = list(
                      title = "Basis", 
                      title_position = "topcenter",
                      direction = "horizontal",
                      legend_position = "top"))
h_v = Heatmap(1 * (human_motif_ind[(1:500)[dhs_subset], motif_sel]> 0),
              show_row_names = F, show_column_names = F,
              show_row_dend = F,
              cluster_rows = dhs_hcl,
              cluster_columns = red_hcl,
              show_column_dend = F,
              heatmap_legend_param = list(
                      title = "Motif Indicator", 
                      title_position = "topcenter",
                      direction = "horizontal",
                      legend_position = "top"
              ))
v = (h %*% t(w))

h_v1 = Heatmap(1 * (v > 5e5),
              show_row_names = F, show_column_names = F,
              show_row_dend = F,
              cluster_rows = dhs_hcl,
              cluster_columns = red_hcl,
              show_column_dend = F,
              heatmap_legend_param = list(
                      title = "Motif Indicator", 
                      title_position = "topcenter",
                      direction = "horizontal",
                      legend_position = "top"
              ))


h_h = draw(h_h, heatmap_legend_side = "top")
h_w = draw(h_w, heatmap_legend_side = "top")
h_v = draw(h_v, heatmap_legend_side = "top")
push_png(h_h, file_name = "motif_h", height = 4., width = 3.0)
push_png(h_w, file_name = "motif_w", height = 3.0, width = 5.)
push_png(h_v, file_name = "motif_v", height = 4.0, width = 5.)

h_h_bin = Heatmap(1 * (human_motif_red_scaled[dhs_subset, ] > 2.),
                  show_row_dend = F, show_column_dend = F, show_row_names = F,
                  cluster_rows = dhs_hcl,
                  cluster_columns = F,
                  heatmap_legend_param = list(
                          title = "Binarized coeficient", 
                          title_position = "topcenter",
                          direction = "horizontal",
                          legend_position = "top"
                          ))
h_h_bin = draw(h_h_bin, heatmap_legend_side = "top")
push_png(h_h_bin, file_name = "motif_h_bin", height = 4, width = 3.0)
# end ploting motif reduction

