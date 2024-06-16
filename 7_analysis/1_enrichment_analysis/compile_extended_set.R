library(tidyverse)
library(ComplexHeatmap)
library(GenomicRanges)
source("./helper/helper.R")

sign_tb_list = c("./non_alignable_env/chromatin_total_ws_sign_fdr0.25.tsv",
                 "./non_alignable_env/H3K27ac_total_ws_sign_fdr0.25.tsv",
                 "./non_alignable_env/H3K4me3_total_ws_sign_fdr0.25.tsv",
                 "./non_alignable_env/H3K4me1_total_ws_sign_fdr0.25.tsv") %>%
        map(., readr::read_tsv, col_types = c("human_id" = "c", "mouse_id" = "c"))
sign_tb = bind_rows(list(mutate(sign_tb_list[[1]], modality = "chromatin"),
                         mutate(sign_tb_list[[2]], modality = "H3K27ac"),
                         mutate(sign_tb_list[[3]], modality = "H3K4me3"),
                         mutate(sign_tb_list[[4]], modality = "H3K4me1")))
table(is.na(sign_tb$weighted_spearman))

sign_tb$mouse_id = paste0("MouseDHS-", sign_tb$mouse_id)
all(sign_tb$mouse_id %in% rownames(mouse_motif_red))
human_indices = match(sign_tb$human_id, rownames(human_motif_red))
mouse_indices = match(sign_tb$mouse_id, rownames(mouse_motif_red))
sign_tb$human_motif = purrr::map(human_indices, function(i) {
        if (any(human_motif_red[i, ] > 0)) {
                x = human_motif_red[i, ]
                which(scale(x) > 2.)
        } else {
                return(NA)
        }
})
sign_tb$mouse_motif = purrr::map(mouse_indices, function(i) {
        if (any(mouse_motif_red[i, ] > 0)) {
                x = mouse_motif_red[i, ]
                which(scale(x) > 2.)
        } else {
                return(numeric())
        }
})
sign_tb$motif_shared = map2(sign_tb$human_motif, sign_tb$mouse_motif, function(x, y) {
        intersect(x, y)
})
sign_tb$n_motif_shared = map_dbl(sign_tb$motif_shared, length)
table(sign_tb$n_motif_shared)
saveRDS(sign_tb, file = "./output/nonalign_sign_combined.rds")

# Figure 2d: plot fraction with motif analysis
source("./helper/def_color.R")
sign_tb = readRDS("./output/nonalign_sign_combined.rds")
sign_tb$motif_sign = sign_tb$n_motif_shared >= 1.
sign_tb$weighted_spearman_bin = cut(sign_tb$weighted_spearman, breaks = c(-Inf, seq(0, 0.75, by = 0.05), 1))
g_frac_motif = sign_tb %>% group_by(weighted_spearman_bin, modality) %>% filter(!is.na(weighted_spearman_bin)) %>%
        summarise(mean_motif = mean(motif_sign),
                  se_motif = sqrt(mean_motif * (1 - mean_motif) / n()),
                  n = n()) %>%
        ggplot(data = .) + geom_line(aes(x = as.numeric(weighted_spearman_bin), y = mean_motif, color = modality)) +
        geom_errorbar(aes(x = as.numeric(weighted_spearman_bin),
                          ymin = mean_motif - se_motif, ymax = mean_motif + se_motif,
                          color = modality)) +
        scale_color_manual(values = data_mod_col1) +
        scale_x_continuous(breaks = unique(as.numeric(sign_tb$weighted_spearman_bin)), labels = unique(sign_tb$weighted_spearman_bin)) +
        labs(x = "GH-CO-V", y = "Fraction with shared TFBS motif")+
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 30, vjust = 0.7),
              legend.position = "right")
push_pdf(g_frac_motif, file = "frac_shared_motif", w = 4, h = 2.5)
