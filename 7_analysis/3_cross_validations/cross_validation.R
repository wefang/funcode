library(tidyverse)
library(ggpubr)
source("./helper/helper.R")
source("./helper/def_color.R")
regions = readRDS("./output/indexed_dhs_mapped_regions_v1_manual.rds")
exp_meta = readr::read_tsv("./metadata/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)

#### begin chromatin accessibility ####
human_correct = readRDS("./intermediate_data/human_chromatin_correct.rds")
mouse_correct = readRDS("./intermediate_data/mouse_chromatin_correct.rds")
human_histone = human_correct
mouse_histone = mouse_correct
# load anchors and matched biosamples
anchors_all = read_tsv(paste0("./output/chromatin_conservation/", "anchors.tsv"))
exp_matched_biosample_new = readRDS(paste0("./custom_metadata/manual_matched_ca.rds"))

# CO-V validation
valid_tb = readr::read_csv("custom_metadata/validation_set.csv")
valid_tb$index1 = match(paste0(valid_tb$biosample1, "_", valid_tb$lifestage1),
                         paste0(exp_matched_biosample_new$biosample_mapped, "_", exp_matched_biosample_new$`Life stage`))
valid_tb$index2 = match(paste0(valid_tb$biosample2, "_", valid_tb$lifestage2),
                         paste0(exp_matched_biosample_new$biosample_mapped, "_", exp_matched_biosample_new$`Life stage`))
valid_tb = valid_tb[!is.na(valid_tb$index1) & !is.na(valid_tb$index2), ]

for (i_valid in 1:nrow(valid_tb)) {
        message(paste0('running valid ', i_valid))
        matched_index1 = valid_tb$index1[i_valid]
        matched_index2 = valid_tb$index2[i_valid]
        exclude_indices1 = which(anchors_all$accession1 %in% exp_matched_biosample_new$exp_tb.x[[matched_index1]]$Accession |
                                         anchors_all$accession2 %in% exp_matched_biosample_new$exp_tb.y[[matched_index1]]$Accession)
        exclude_indices2 = which(anchors_all$accession1 %in% exp_matched_biosample_new$exp_tb.x[[matched_index2]]$Accession |
                                         anchors_all$accession2 %in% exp_matched_biosample_new$exp_tb.y[[matched_index2]]$Accession)
        if (length(c(exclude_indices1, exclude_indices2)) == 0) {
                anchors = anchors_all
        } else {
                anchors = anchors_all[-union(exclude_indices1, exclude_indices2), ]
        }
        set.seed(37)
        null_size = 5e5
        loci_cor = map_dbl(1:nrow(human_histone), function(j) {
                if (j %% 1e5 == 0) message(j)
                x = human_histone[j, anchors$accession1]
                y = mouse_histone[j, anchors$accession2]
                if (var(x) > 0 & var(y) > 0) {
                        wCorr::weightedCorr(x,
                                            y,
                                            weights = anchors$score,
                                            method = "spearman")
                } else {
                        return(NA)
                }
        })
        human_var_std = calc_var_std_ver1(human_histone[, anchors$accession1])
        mouse_var_std = calc_var_std_ver1(mouse_histone[, anchors$accession2])

        matched_sample_excluded = exp_matched_biosample_new[-c(matched_index1, matched_index2), ]
        human_matched_acc = purrr::reduce(map(matched_sample_excluded$exp_tb.x, "Accession"), c)
        human_matched_biosample = rep(paste0(matched_sample_excluded$`Life stage`, "_",
                                             matched_sample_excluded$biosample_mapped), map_dbl(matched_sample_excluded$exp_tb.x, nrow))
        mouse_matched_acc = reduce(map(matched_sample_excluded$exp_tb.y, "Accession"), c)
        mouse_matched_biosample = rep(paste0(matched_sample_excluded$`Life stage`, "_",
                                             matched_sample_excluded$biosample_mapped), map_dbl(matched_sample_excluded$exp_tb.y, nrow))
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
        human_var_std_manual = calc_var_std_ver1(human_correct_manual[, matched_biosample])
        mouse_var_std_manual = calc_var_std_ver1(mouse_correct_manual[, matched_biosample])

        output_tb = tibble(id = rownames(human_histone),
                           cov = loci_cor,
                           cov_manual = loci_cor_manual,
                           human_var_std = human_var_std,
                           mouse_var_std = mouse_var_std,
                           human_var_std_manual = human_var_std_manual,
                           mouse_var_std_manual = mouse_var_std_manual)
        saveRDS(output_tb, output_rds_path(paste0(i_valid, "_score")))

        output_tb = readRDS(output_rds_path(paste0(i_valid, "_score")))
        num_sim = 20
        loci_sample_size = 1000
        human_diff = rowMeans(log2(1+human_histone[, exp_matched_biosample_new$exp_tb.x[[matched_index1]]$Accession, drop = F])) -
                rowMeans(log2(1+human_histone[, exp_matched_biosample_new$exp_tb.x[[matched_index2]]$Accession, drop = F]))
        mouse_diff = rowMeans(log2(1+mouse_histone[, exp_matched_biosample_new$exp_tb.y[[matched_index1]]$Accession, drop = F])) -
                rowMeans(log2(1+mouse_histone[, exp_matched_biosample_new$exp_tb.y[[matched_index2]]$Accession, drop = F]))
        
        bg_set_list = list("CO-V" = order(rank(output_tb$human_var_std) +
                                                         rank(output_tb$mouse_var_std), decreasing = T)[1:5e5],
                           "CO-V-Manual" = order(rank(output_tb$human_var_std_manual) +
                                                                rank(output_tb$mouse_var_std_manual), decreasing = T)[1:5e5])
        sim_df = as_tibble(expand.grid(
                Method = c("CO-V",
                           "CO-V-Manual",
                           "phastCons4way",
                           "phastCons7way",
                           "phastCons17way",
                           "phastCons20way",
                           "phastCons30way",
                           "phastCons100way",
                           "phyloP4way",
                           "phyloP7way",
                           "phyloP17way",
                           "phyloP20way",
                           "phyloP30way",
                           "phyloP100way",
                           "percentage",
                           "lecif"),
                min_rank = c(5e3, 10e3, 1.5e4, 2e4, 2.5e4, 5e4, 7.5e4, 10e4),
                sim_id = 1:num_sim))
        sim_df = arrange(sim_df, Method)
        sim_df$data = purrr::map(1:nrow(sim_df), function(i) {
                if (i %% 100 == 0) message(i)
                method_name = as.character(sim_df$Method[i])
                if (method_name %in% names(bg_set_list)) {
                        bg_set = bg_set_list[[method_name]]
                } else {
                        bg_set = 1:nrow(regions)
                }
                if (method_name == "CO-V") {
                        bg_score = output_tb$cov[bg_set]
                }
                if (method_name == "CO-V-Manual") {
                        bg_score = output_tb$cov_manual[bg_set]
                }
                if (!method_name %in% c("CO-V", "CO-V-Manual")) {
                        bg_score = regions[[method_name]][bg_set]
                }
                min_rank = sim_df$min_rank[i]
                fg_set = sample(bg_set[order(bg_score, decreasing = T)][1:min_rank], loci_sample_size)
                # sample loci from bg that matches mouse fg logFC
                fg_strat_count = table(as_strata(mouse_diff[fg_set]))
                bg_strat = as_strata(mouse_diff[bg_set])
                bg_set_matched = levels(bg_strat) %>% purrr::map(function(x) {
                        # make sure there are enough loci in the bg strata
                        sample(bg_set[which(bg_strat == x)], size = fg_strat_count[x])
                }) %>% purrr::reduce(c)
                tibble(Conserved = cor(human_diff[fg_set], mouse_diff[fg_set], method = "spearman"),
                       Control = cor(human_diff[bg_set_matched], mouse_diff[bg_set_matched], method = "spearman"))
        })
        sim_df = sim_df %>% unnest(cols = c(data))
        sim_df_long = gather(sim_df, key = type, value = cor, c(Conserved, Control))
        saveRDS(sim_df_long, file = output_rds_path(paste0("res_complete_", i_valid)))
}

valid_res_tb = bind_rows(map(1:nrow(valid_tb), function(i_valid) {
        mutate(readRDS(output_rds_path(paste0("res_complete_", i_valid))), index = i_valid)
}))

rank_lab_braks = c(5, 10, 50, 100, 150)
(valid_res_tb %>%
        group_by(Method, type, min_rank) %>%
        summarize(mean_cor = mean(cor, na.rm =T),
                  se_cor = sd(cor, na.rm =T)/sqrt(n())) %>%
        filter(type == "Conserved") %>%
        ggline(x = "min_rank", y = "mean_cor", color = "Method", add = "mean_se", numeric.x.axis = T) +
        scale_color_manual(values = score_col) + 
        scale_x_continuous(labels = paste0(rank_lab_braks, "k"),
                             breaks = 10^3 * rank_lab_braks) +
                theme(legend.position = "right",
                        legend.text = element_text(size = 10),
                        text = element_text(size = 12))) %>%
        push_pdf("eval_fc_ave_ca_complete", width = 5., height = 3.5)

valid_res_tb %>%
        group_by(Method, type, min_rank) %>%
        summarize(mean_cor = mean(cor, na.rm =T),
                  se_cor = sd(cor, na.rm =T)/sqrt(n())) %>%
        filter(type == "Conserved") %>%
        ggbarplot(x = "min_rank", y = "mean_cor", fill = "Method", position = position_dodge())

# CO-B validation
library(furrr)
plan(multisession, workers = 12)
options(future.globals.maxSize= 8*1024^3)
for (i_valid in 19:nrow(exp_matched_biosample_new)) {
        message(paste0('running valid ', i_valid))
        matched_index1 = i_valid
        exclude_indices1 = which(anchors_all$accession1 %in% exp_matched_biosample_new$exp_tb.x[[matched_index1]]$Accession |
                                         anchors_all$accession2 %in% exp_matched_biosample_new$exp_tb.y[[matched_index1]]$Accession)
        if (length(exclude_indices1) == 0) {
                anchors = anchors_all
        } else {
                anchors = anchors_all[-exclude_indices1, ]
        }
        set.seed(37)
        prob_vec = seq(0.1, 1.0, by = 0.01)
        human_quantile_vec = quantile(human_histone[sample(nrow(human_histone), 1e5), anchors$accession1], prob_vec)
        mouse_quantile_vec = quantile(mouse_histone[sample(nrow(mouse_histone), 1e5), anchors$accession2], prob_vec)
        compute_hval <- function(human_mat, mouse_mat, human_quantile_vec, mouse_quantile_vec, prob_vec) {
                human_hval_ind = do.call(cbind, map(1:length(prob_vec), function(i) {
                        rowMeans(human_mat > human_quantile_vec[i]) > prob_vec[i]
                }))
                human_hval = apply(human_hval_ind, 1, function(x) {
                        if (!any(x)) {
                                return(0)
                        }
                        prob_vec[tail(which(x), 1)]
                })
                mouse_hval_ind = do.call(cbind, map(1:length(prob_vec), function(i) {
                        rowMeans(mouse_mat > mouse_quantile_vec[i]) > prob_vec[i]
                }))
                mouse_hval = apply(mouse_hval_ind, 1, function(x) {
                        if (!any(x)) {
                                return(0)
                        }
                        prob_vec[tail(which(x), 1)]
                })
                hval = pmin(human_hval, mouse_hval)
                hval
        }
        loci_hval = compute_hval(human_histone[, anchors$accession1],
                                 mouse_histone[, anchors$accession2],
                                 human_quantile_vec, mouse_quantile_vec, prob_vec)

        human_mean_train = rowMeans(human_histone[, anchors$accession1])
        mouse_mean_train = rowMeans(mouse_histone[, anchors$accession2])

        output_tb = tibble(id = rownames(human_histone),
                           cob = loci_hval,
                           human_mean = human_mean_train,
                           mouse_mean = mouse_mean_train)
        saveRDS(output_tb, output_rds_path(paste0(i_valid, "_cob_score")))

        output_tb = readRDS(output_rds_path(paste0(i_valid, "_cob_score")))
        human_mean = rowMeans(log2(1+human_histone[, exp_matched_biosample_new$exp_tb.x[[matched_index1]]$Accession, drop = F]))
        mouse_mean = rowMeans(log2(1+mouse_histone[, exp_matched_biosample_new$exp_tb.y[[matched_index1]]$Accession, drop = F]))
        names(human_mean) = regions$id
        names(mouse_mean) = regions$id
        
        human_cutoff = quantile(human_mean, 0.95)
        mouse_cutoff = quantile(mouse_mean, 0.95)
        
        num_sim = 20
        loci_sample_size = 1000
        sim_df = as_tibble(expand.grid(
                Method = c("CO-B",
                           # "CO-B Human",
                           # "CO-B Mouse",
                           "phastCons4way",
                           "phastCons7way",
                           "phastCons17way",
                           "phastCons20way",
                           "phastCons30way",
                           "phastCons100way",
                           "phyloP4way",
                           "phyloP7way",
                           "phyloP17way",
                           "phyloP20way",
                           "phyloP30way",
                           "phyloP100way",
                           "percentage",
                           "lecif"),
                min_rank = c(1e4, 2.5e4, 5e4, 7.5e4, 10e4, 20e4),
                sim_id = 1:num_sim))
        sim_df = arrange(sim_df, Method)
        sim_df$data = future_map(1:nrow(sim_df), function(i) {
                if (i %% 100 == 0) message(i)
                bg_set = regions$id[!is.na(regions$lecif)]
                bg_indices = which(!is.na(regions$lecif))
                method_name = as.character(sim_df$Method[i])
                if (method_name == "CO-B") {
                        bg_score = output_tb$cob[bg_indices]
                }
                if (!method_name %in% "CO-B") {
                        bg_score = regions[[method_name]][bg_indices]
                }
                min_rank = sim_df$min_rank[i]
                fg_set = sample(bg_set[order(bg_score, decreasing = T)][1:min_rank], loci_sample_size)
                fg_strat_count = table(as_strata(mouse_mean[fg_set],
                                                 c(seq(0, 10, by = 1), Inf)))
                bg_strat = as_strata(mouse_mean[bg_set],
                                     c(seq(0, 10, by = 1), Inf))
                bg_set_matched = levels(bg_strat) %>% purrr::map(function(x) {
                        # make sure there are enough loci in the bg strata
                        sample(bg_set[which(bg_strat == x)], size = fg_strat_count[x])
                }) %>% purrr::reduce(c)
                tibble(Conserved = sum(human_mean[fg_set] > human_cutoff & mouse_mean[fg_set] > mouse_cutoff) /
                               sum(human_mean[fg_set] > human_cutoff | mouse_mean[fg_set] > mouse_cutoff),
                       Control = sum(human_mean[bg_set_matched] > human_cutoff & mouse_mean[bg_set_matched] > mouse_cutoff) /
                               sum(human_mean[bg_set_matched] > human_cutoff | mouse_mean[bg_set_matched] > mouse_cutoff))
        }, .progress = T)
        sim_df = sim_df %>% unnest(cols = c(data))
        sim_df_long = gather(sim_df, key = type, value = accuracy, c(Conserved, Control))
        saveRDS(sim_df_long, file = output_rds_path(paste0("res_complete_cob_", i_valid)))
}
valid_res_tb = bind_rows(map(1:nrow(exp_matched_biosample_new), function(i_valid) {
        mutate(readRDS(output_rds_path(paste0("res_complete_cob_", i_valid))), index = i_valid)
}))

valid_set_proc_ave_sum = valid_res_tb %>% group_by(Method, type, min_rank) %>%
        summarize(mean_cor = mean(accuracy, na.rm =T),
                  se_cor = sd(accuracy, na.rm =T)/sqrt(n()))
# valid_set_proc_ave_sum$Method = factor(valid_set_proc_ave_sum$Method,
#                                        levels = rev(levels(valid_set_proc_ave_sum$Method)))
valid_set_proc_ave_sum$type = factor(valid_set_proc_ave_sum$type, levels = c("Control", "Conserved"))
rank_lab_braks = c(5, 10, 50, 100)

(ggplot(valid_set_proc_ave_sum[valid_set_proc_ave_sum$type == "Conserved" & valid_set_proc_ave_sum$min_rank <= 1e5, ],
        aes(x = min_rank,
            y = mean_cor,
            # alpha = type,
            color = Method,
            # linetype = type
        )) +
                geom_line(size = 0.5) +
                geom_point(size = 0.5) +
                geom_errorbar(aes(ymin = mean_cor - se_cor,
                                  ymax = mean_cor + se_cor,
                                  color = Method), size = 1, width = 0.) +
                scale_color_manual(values = score_col) +
                ylab("Jaccard index") +
                xlab("Element Rank") +
                theme_pubr() + theme(legend.position = "right",
                                     legend.text = element_text(size = 10),
                                     text = element_text(size = 12)) + 
                scale_x_continuous(labels = paste0(rank_lab_braks, "k"),
                                   breaks = 10^3 * rank_lab_braks)) %>%
        push_pdf(paste0("cross_valid_ca_cob_ave_complete"), width = 5.0, height = 3.5)
#### end chromatin accessibility ####

#### begin histone modification
assay = "H3K4me1"
# assay = "H3K27ac"
# assay = "H3K4me3"
human_histone = readRDS(paste0("./processed_data/Histone_norm_correct/Homo_sapiens_", assay, "_norm_correct.rds"))
mouse_histone = readRDS(paste0("./processed_data/Histone_norm_correct/Mus_musculus_", assay, "_norm_correct.rds"))

human_exp_tb = filter(exp_meta, Organism == "Homo sapiens" & exp_meta$`Target of assay` == assay)
mouse_exp_tb = filter(exp_meta, Organism == "Mus musculus" & exp_meta$`Target of assay` == assay)

exp_matched_biosample_new = readRDS(paste0("./custom_metadata/manual_matched_", assay, ".rds"))
anchors_all = readRDS(paste0("./output/histone_conservation_v1/", "anchors_", assay, "_v1.rds"))
valid_tb = readRDS(paste0("custom_metadata/validation_set_", assay, ".rds"))
script_output_dir = paste0("./output/cross_valid/", assay, "/")
# dir.create(script_output_dir, recursive = T)

library(furrr)
plan(multisession, workers = 12)
options(future.globals.maxSize= 8*1024^3)
for (i_valid in 1:nrow(valid_tb)) {
        message(paste0('running valid ', i_valid))
        matched_index1 = valid_tb$index1[i_valid]
        matched_index2 = valid_tb$index2[i_valid]
        exclude_indices1 = which(anchors_all$accession1 %in% exp_matched_biosample_new$exp_tb.x[[matched_index1]]$Accession |
                                         anchors_all$accession2 %in% exp_matched_biosample_new$exp_tb.y[[matched_index1]]$Accession)
        exclude_indices2 = which(anchors_all$accession1 %in% exp_matched_biosample_new$exp_tb.x[[matched_index2]]$Accession |
                                         anchors_all$accession2 %in% exp_matched_biosample_new$exp_tb.y[[matched_index2]]$Accession)
        if (length(c(exclude_indices1, exclude_indices2)) == 0) {
                anchors = anchors_all
        } else {
                anchors = anchors_all[-union(exclude_indices1, exclude_indices2), ]
        }
        set.seed(37)
        null_size = 5e5
        loci_cor = map_dbl(1:nrow(human_histone), function(j) {
                if (j %% 1e5 == 0) message(j)
                x = human_histone[j, anchors$accession1]
                y = mouse_histone[j, anchors$accession2]
                if (var(x) > 0 & var(y) > 0) {
                        wCorr::weightedCorr(x,
                                            y,
                                            weights = anchors$score,
                                            method = "spearman")
                } else {
                        return(NA)
                }
        })
        human_var_std = calc_var_std_ver1(human_histone[, anchors$accession1])
        mouse_var_std = calc_var_std_ver1(mouse_histone[, anchors$accession2])

        matched_sample_excluded = exp_matched_biosample_new[-c(matched_index1, matched_index2), ]
        human_matched_acc = purrr::reduce(map(matched_sample_excluded$exp_tb.x, "Accession"), c)
        human_matched_biosample = rep(paste0(matched_sample_excluded$`Life stage`, "_",
                                             matched_sample_excluded$biosample_mapped), map_dbl(matched_sample_excluded$exp_tb.x, nrow))
        mouse_matched_acc = purrr::reduce(map(matched_sample_excluded$exp_tb.y, "Accession"), c)
        mouse_matched_biosample = rep(paste0(matched_sample_excluded$`Life stage`, "_",
                                             matched_sample_excluded$biosample_mapped), map_dbl(matched_sample_excluded$exp_tb.y, nrow))
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
        human_var_std_manual = calc_var_std_ver1(human_correct_manual[, matched_biosample])
        mouse_var_std_manual = calc_var_std_ver1(mouse_correct_manual[, matched_biosample])

        output_tb = tibble(id = rownames(human_histone),
                           cov = loci_cor,
                           cov_manual = loci_cor_manual,
                           human_var_std = human_var_std,
                           mouse_var_std = mouse_var_std,
                           human_var_std_manual = human_var_std_manual,
                           mouse_var_std_manual = mouse_var_std_manual)
        saveRDS(output_tb, output_rds_path(paste0(i_valid, "_score")))

        output_tb = readRDS(output_rds_path(paste0(i_valid, "_score")))
        num_sim = 20
        loci_sample_size = 1000
        human_diff = rowMeans(log2(1+human_histone[, exp_matched_biosample_new$exp_tb.x[[matched_index1]]$Accession, drop = F])) -
                rowMeans(log2(1+human_histone[, exp_matched_biosample_new$exp_tb.x[[matched_index2]]$Accession, drop = F]))
        mouse_diff = rowMeans(log2(1+mouse_histone[, exp_matched_biosample_new$exp_tb.y[[matched_index1]]$Accession, drop = F])) -
                rowMeans(log2(1+mouse_histone[, exp_matched_biosample_new$exp_tb.y[[matched_index2]]$Accession, drop = F]))
        
        bg_set_list = list("CO-V" = order(rank(output_tb$human_var_std) +
                                                  rank(output_tb$mouse_var_std), decreasing = T)[1:5e5],
                           "CO-V-Manual" = order(rank(output_tb$human_var_std_manual) +
                                                         rank(output_tb$mouse_var_std_manual), decreasing = T)[1:5e5])
        sim_df = as_tibble(expand.grid(
                Method = c("CO-V",
                           "CO-V-Manual",
                           "phastCons4way",
                           "phastCons7way",
                           "phastCons17way",
                           "phastCons20way",
                           "phastCons30way",
                           "phastCons100way",
                           "phyloP4way",
                           "phyloP7way",
                           "phyloP17way",
                           "phyloP20way",
                           "phyloP30way",
                           "phyloP100way",
                           "percentage",
                           "lecif"),
                min_rank = c(5e3, 10e3, 1.5e4, 2e4, 2.5e4, 5e4, 7.5e4, 10e4),
                sim_id = 1:num_sim))
        
        sim_df = arrange(sim_df, Method)
        sim_df$data = future_map(1:nrow(sim_df), function(i) {
                if (i %% 100 == 0) message(i)
                method_name = as.character(sim_df$Method[i])
                if (method_name %in% names(bg_set_list)) {
                        bg_set = bg_set_list[[method_name]]
                } else {
                        bg_set = 1:nrow(regions)
                }
                if (method_name == "CO-V") {
                        bg_score = output_tb$cov[bg_set]
                }
                if (method_name == "CO-V-Manual") {
                        bg_score = output_tb$cov_manual[bg_set]
                }
                if (!method_name %in% c("CO-V", "CO-V-Manual")) {
                        bg_score = regions[[method_name]][bg_set]
                }
                min_rank = sim_df$min_rank[i]
                fg_set = sample(bg_set[order(bg_score, decreasing = T)][1:min_rank], loci_sample_size)
                
                # sample loci from bg that matches mouse fg logFC
                # fg_strat_count = table(as_strata(mouse_diff[fg_set]))
                # bg_strat = as_strata(mouse_diff[bg_set])
                # bg_set_matched = levels(bg_strat) %>% purrr::map(function(x) {
                #         # make sure there are enough loci in the bg strata
                #         sample(bg_set[which(bg_strat == x)], size = fg_strat_count[x])
                # }) %>% purrr::reduce(c)
                tibble(Conserved = cor(human_diff[fg_set], mouse_diff[fg_set], method = "spearman"))#,
                # Control = cor(human_diff[bg_set_matched], mouse_diff[bg_set_matched], method = "spearman"))
        }, .progress = T)
        sim_df = sim_df %>% unnest(cols = c(data))
        sim_df_long = gather(sim_df, key = type, value = cor, c(Conserved))
        saveRDS(sim_df_long, file = output_rds_path(paste0("res_complete_", i_valid))) # v1 includes PhastCons100way
}

# reading output and produce results
valid_res_tb = bind_rows(map(1:nrow(valid_tb), function(i_valid) {
        mutate(readRDS(output_rds_path(paste0("res_complete_", i_valid))), index = i_valid)
}))
valid_res_tb %>% group_by(index, Method, type, min_rank) %>%
        summarize(mean_cor = mean(cor, na.rm =T),
                  se_cor = sd(cor, na.rm =T)/sqrt(n())) %>%
        ggline(x = "min_rank", y = "mean_cor", color = "Method", add = "mean_se",
               facet.by = "index")

valid_set_proc_ave_sum = valid_res_tb %>% group_by(Method, type, min_rank) %>%
        summarize(mean_cor = mean(cor, na.rm =T),
                  se_cor = sd(cor, na.rm =T)/sqrt(n()))
valid_set_proc_ave_sum$type = factor(valid_set_proc_ave_sum$type, levels = c("Conserved"))
valid_set_proc_ave_sum$Method = factor(valid_set_proc_ave_sum$Method, levels = names(score_col))
rank_lab_braks = c(5, 10, 50, 100, 150)
(ggplot(valid_set_proc_ave_sum[valid_set_proc_ave_sum$min_rank <= 10e4 & valid_set_proc_ave_sum$type == "Conserved", ],
        aes(x = min_rank, y = mean_cor,
            # alpha = type, linetype = type,
            color = Method)) +
                geom_line(linewidth = 0.5) +
                geom_point(size = 0.5) +
                # geom_errorbar(aes(ymin = mean_cor - se_cor,
                #                   ymax = mean_cor + se_cor,
                #                   color = Method), size = 0.25, width = 0.0) +
                scale_color_manual(values = score_col) +
                # scale_alpha_manual(values = c("Control" = 0.5, "Conserved" = 1.0)) +
                ylab("Correlation of Fold Changes") +
                # scale_linetype_manual(values = c("Control" = 2, "Conserved" = 1)) +
                # ylab("Sensitivity") +
                xlab("DHS Rank") +
                theme_pubr() + theme(legend.position = "right",
                                     legend.text = element_text(size = 10),
                                     text = element_text(size = 12)) +
                scale_x_continuous(labels = paste0(rank_lab_braks, "k"),
                                   breaks = 10^3 * rank_lab_braks)) %>% 
        push_pdf(paste0("cov_eval_", assay, "_complete"), width = 4.5, height = 3.)
# end CO-V analysis

# begin CO-B analysis
for (i_valid in 1:nrow(exp_matched_biosample_new)) {
        message(paste0('running valid ', i_valid))
        matched_index1 = i_valid
        exclude_indices1 = which(anchors_all$accession1 %in% exp_matched_biosample_new$exp_tb.x[[matched_index1]]$Accession |
                                         anchors_all$accession2 %in% exp_matched_biosample_new$exp_tb.y[[matched_index1]]$Accession)
        if (length(exclude_indices1) == 0) {
                anchors = anchors_all
        } else {
                anchors = anchors_all[-exclude_indices1, ]
        }
        set.seed(37)
        prob_vec = seq(0.1, 1.0, by = 0.01)
        human_quantile_vec = quantile(human_histone[sample(nrow(human_histone), 1e5), anchors$accession1], prob_vec)
        mouse_quantile_vec = quantile(mouse_histone[sample(nrow(mouse_histone), 1e5), anchors$accession2], prob_vec)
        compute_hval <- function(human_mat, mouse_mat, human_quantile_vec, mouse_quantile_vec, prob_vec) {
                human_hval_ind = do.call(cbind, map(1:length(prob_vec), function(i) {
                        rowMeans(human_mat > human_quantile_vec[i]) > prob_vec[i]
                }))
                human_hval = apply(human_hval_ind, 1, function(x) {
                        if (!any(x)) {
                                return(0)
                        }
                        prob_vec[tail(which(x), 1)]
                })
                mouse_hval_ind = do.call(cbind, map(1:length(prob_vec), function(i) {
                        rowMeans(mouse_mat > mouse_quantile_vec[i]) > prob_vec[i]
                }))
                mouse_hval = apply(mouse_hval_ind, 1, function(x) {
                        if (!any(x)) {
                                return(0)
                        }
                        prob_vec[tail(which(x), 1)]
                })
                hval = pmin(human_hval, mouse_hval)
                hval
        }
        loci_hval = compute_hval(human_histone[, anchors$accession1],
                                 mouse_histone[, anchors$accession2],
                                 human_quantile_vec, mouse_quantile_vec, prob_vec)

        human_mean_train = rowMeans(human_histone[, anchors$accession1])
        mouse_mean_train = rowMeans(mouse_histone[, anchors$accession2])

        output_tb = tibble(id = rownames(human_histone),
                           cob = loci_hval,
                           human_mean = human_mean_train,
                           mouse_mean = mouse_mean_train)
        saveRDS(output_tb, output_rds_path(paste0(i_valid, "_cob_score")))
        
        output_tb = readRDS(output_rds_path(paste0(i_valid, "_cob_score")))
        human_mean = rowMeans(log2(1+human_histone[, exp_matched_biosample_new$exp_tb.x[[matched_index1]]$Accession, drop = F]))
        mouse_mean = rowMeans(log2(1+mouse_histone[, exp_matched_biosample_new$exp_tb.y[[matched_index1]]$Accession, drop = F]))
        names(human_mean) = regions$id
        names(mouse_mean) = regions$id
        
        human_cutoff = quantile(human_mean, 0.95)
        mouse_cutoff = quantile(mouse_mean, 0.95)
        
        num_sim = 20
        loci_sample_size = 1000
        sim_df = as_tibble(expand.grid(
                Method = c("CO-B",
                           # "CO-B Human",
                           # "CO-B Mouse",
                           "phastCons4way",
                           "phastCons7way",
                           "phastCons17way",
                           "phastCons20way",
                           "phastCons30way",
                           "phastCons100way",
                           "phyloP4way",
                           "phyloP7way",
                           "phyloP17way",
                           "phyloP20way",
                           "phyloP30way",
                           "phyloP100way",
                           "percentage",
                           "lecif"),
                min_rank = c(1e4, 2.5e4, 5e4, 7.5e4, 10e4, 20e4),
                sim_id = 1:num_sim))
        sim_df = arrange(sim_df, Method)
        sim_df$data = future_map(1:nrow(sim_df), function(i) {
                if (i %% 10 == 0) message(i)
                bg_set = regions$id[!is.na(regions$lecif)]
                bg_indices = which(!is.na(regions$lecif))
                method_name = as.character(sim_df$Method[i])
                if (method_name == "CO-B") {
                        bg_score = output_tb$cob[bg_indices]
                }
                if (!method_name %in% "CO-B") {
                        bg_score = regions[[method_name]][bg_indices]
                }
                min_rank = sim_df$min_rank[i]
                fg_set = sample(bg_set[order(bg_score, decreasing = T)][1:min_rank], loci_sample_size)
                fg_strat_count = table(as_strata(mouse_mean[fg_set],
                                                 c(seq(0, 10, by = 1), Inf)))
                bg_strat = as_strata(mouse_mean[bg_set],
                                     c(seq(0, 10, by = 1), Inf))
                bg_set_matched = levels(bg_strat) %>% purrr::map(function(x) {
                        # make sure there are enough loci in the bg strata
                        sample(bg_set[which(bg_strat == x)], size = fg_strat_count[x])
                }) %>% purrr::reduce(c)
                tibble(Conserved = sum(human_mean[fg_set] > human_cutoff & mouse_mean[fg_set] > mouse_cutoff) /
                               sum(human_mean[fg_set] > human_cutoff | mouse_mean[fg_set] > mouse_cutoff),
                       Control = sum(human_mean[bg_set_matched] > human_cutoff & mouse_mean[bg_set_matched] > mouse_cutoff) /
                               sum(human_mean[bg_set_matched] > human_cutoff | mouse_mean[bg_set_matched] > mouse_cutoff))
        }, .progress = T)
        sim_df = sim_df %>% unnest(cols = c(data))
        sim_df_long = gather(sim_df, key = type, value = accuracy, c(Conserved, Control))
        saveRDS(sim_df_long, file = output_rds_path(paste0("res_complete_cob_", i_valid)))
}

valid_res_tb = bind_rows(map(1:nrow(exp_matched_biosample_new), function(i_valid) {
        mutate(readRDS(output_rds_path(paste0("res_complete_cob_", i_valid))), index = i_valid)
}))

valid_set_proc_ave_sum = valid_res_tb %>% group_by(Method, type, min_rank) %>%
        summarize(mean_cor = mean(accuracy, na.rm =T),
                  se_cor = sd(accuracy, na.rm =T)/sqrt(n()))
valid_set_proc_ave_sum$Method = factor(valid_set_proc_ave_sum$Method,
                                       levels = rev(levels(valid_set_proc_ave_sum$Method)))
valid_set_proc_ave_sum$type = factor(valid_set_proc_ave_sum$type, levels = c("Control", "Conserved"))
rank_lab_braks = c(5, 10, 50, 100)

(ggplot(valid_set_proc_ave_sum[valid_set_proc_ave_sum$type == "Conserved" & valid_set_proc_ave_sum$min_rank <= 1e5, ],
        aes(x = min_rank,
            y = mean_cor,
            # alpha = type,
            # linetype = type,
            color = Method
        )) +
                geom_line(size = 0.5) +
                geom_point(size = 0.5) +
                # geom_errorbar(aes(ymin = mean_cor - se_cor,
                #                   ymax = mean_cor + se_cor,
                #                   color = Method), size = 1, width = 0.) +
                scale_color_manual(values = score_col) +
                # scale_alpha_manual(values = c("Control" = 0.5, "Conserved" = 1.0)) +
                ylab("Jaccard index") +
                # scale_linetype_manual(values = c("Control" = 2, "Conserved" = 1)) +
                # ylab("Sensitivity") +
                xlab("DNA Element Rank") +
                theme_pubr() + theme(legend.position = "right",
                                     legend.text = element_text(size = 10),
                                     text = element_text(size = 12)) + 
                scale_x_continuous(labels = paste0(rank_lab_braks, "k"),
                                   breaks = 10^3 * rank_lab_braks)) %>%
        push_pdf(paste0("cross_valid_", assay, "_cob_ave_complete"), width = 4.5, height = 3.)
#### end histone modification ####
