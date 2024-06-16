library(rtracklayer)
library(Seurat)
library(ComplexHeatmap)
# ENCODE specific gene annotations
GRCh38_gencode_ref = "./metadata/ENCODEv4_Reference/ENCFF159KBI.gtf.gz"
mm10_gencode_ref = "./metadata/ENCODEv4_Reference/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz"
source("./load_homolog.R")
source("./helper/helper.R")

script_output_dir = "./output/rna_conservation/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}

script_plot_dir = "./plots/rna_conservation/"
if (!dir.exists(script_plot_dir)) {
        dir.create(script_plot_dir)
}
hg38_tpm = readRDS("./processed_data/ENCODE_RNAseq/hg38_tpm_exp.rds")
mm10_tpm = readRDS("./processed_data/ENCODE_RNAseq/mm10_tpm_exp.rds")
hg38_exp = readr::read_tsv("./processed_data/ENCODE_RNAseq/hg38_experiments.tsv", skip = 1)
mm10_exp = readr::read_tsv("./processed_data/ENCODE_RNAseq/mm10_experiments.tsv", skip = 1)
hg38_gene_ids = rownames(hg38_tpm)
mm10_gene_ids = rownames(mm10_tpm)

hg38_gtf = make_genes_gr(import(GRCh38_gencode_ref))
mm10_gtf = make_genes_gr(import(mm10_gencode_ref))

hg38_gene_names = gene_id_to_names(hg38_gene_ids, hg38_gtf)
mm10_gene_names = gene_id_to_names(mm10_gene_ids, mm10_gtf)

cor_mat = cor(t(hg38_tpm[1:1000, ]))
cor_mat[is.na(cor_mat)] = -1
Heatmap(cor_mat)

common_indices = which((human_symbols %in% hg38_gene_names) & (mouse_symbols %in% mm10_gene_names))
homolog_id = homolog_id[common_indices]
human_symbols = human_symbols[common_indices]
mouse_symbols = mouse_symbols[common_indices]

human_match_indices = match(human_symbols, hg38_gene_names)
mouse_match_indices = match(mouse_symbols, mm10_gene_names)
hg38_rna_matched = hg38_tpm[human_match_indices, ]
mm10_rna_matched = mm10_tpm[mouse_match_indices, ]
rownames(hg38_rna_matched) = homolog_id
rownames(mm10_rna_matched) = homolog_id

# Seurat objects
hg38_rna_obj = CreateSeuratObject(counts = hg38_rna_matched)
mm10_rna_obj = CreateSeuratObject(counts = mm10_rna_matched)

hg38_rna_obj = FindVariableFeatures(hg38_rna_obj, nfeatures = 2000)
mm10_rna_obj = FindVariableFeatures(mm10_rna_obj, nfeatures = 2000)
hg38_rna_obj = NormalizeData(hg38_rna_obj)
mm10_rna_obj = NormalizeData(mm10_rna_obj)
hg38_rna_obj$source = "Human RNA-seq"
mm10_rna_obj$source = "Mouse RNA-seq"

hg38_rna_obj = add_metadata(hg38_rna_obj, hg38_exp)
mm10_rna_obj = add_metadata(mm10_rna_obj, mm10_exp)

hg38_rna_obj = quick_process(hg38_rna_obj, var_feature = VariableFeatures(hg38_rna_obj))
mm10_rna_obj = quick_process(mm10_rna_obj, var_feature = VariableFeatures(mm10_rna_obj))
# UMAPPlot(hg38_rna_obj, group.by = "Lab")
# UMAPPlot(hg38_rna_obj, group.by = "Assay title")

# begin integration
anchors_rna <-
        FindIntegrationAnchors(
                object.list = list(
                        hg38_rna = hg38_rna_obj,
                        mm10_rna = mm10_rna_obj
                ),
                anchor.features = 2000,
                k.filter = 40,
                dims = 1:30
        )
integrated_rna <- IntegrateData(anchorset = anchors_rna,
                                dims = 1:30)
integrated_rna  = quick_process(integrated_rna, var_feature = VariableFeatures(integrated_rna))

(UMAPPlot(integrated_rna, group.by = "source", pt.size = 0.5) +
        scale_color_manual(values = c("#b2df8a", "#1f78b4")) + 
        theme(legend.position = "top", text = element_text(size = 8))) %>%
        push_pdf("umap_rna_integrated", width = 2.7, height = 2.8)
HoverLocator(plot = UMAPPlot(integrated_rna, group.by = "source", pt.size = 3) +
                     scale_color_manual(values = c("#b2df8a", "#1f78b4")),
             information = FetchData(integrated_rna, vars = c("source", "Biosample term name", "Life stage", "Assay title", "Lab"))) %>%
        push_widget(file_name = "umap_rna_integrated")

anchor_tb = make_anchor_tb(anchors_rna, hg38_rna_obj, mm10_rna_obj)
push_tsv(anchor_tb, "anchors")

# begin computing spearman
# libsize normalizing samples
hg38_rna_matched = t(t(hg38_rna_matched) / colSums(hg38_rna_matched) * 1e6)
mm10_rna_matched = t(t(mm10_rna_matched) / colSums(mm10_rna_matched) * 1e6)
gene_cor = sapply(1:nrow(hg38_rna_matched), function(j) {
        wCorr::weightedCorr(hg38_rna_matched[j, anchor_tb$accession1],
                            mm10_rna_matched[j, anchor_tb$accession2],
                            weights = anchor_tb$score,
                            method = "spearman")
})
# constructing null
set.seed(73)
null_size = 1e5
sample_indices1 = sample(nrow(hg38_rna_matched), null_size, replace = T)
sample_indices2 = sample(nrow(hg38_rna_matched), null_size, replace = T)
human_null = hg38_rna_matched[sample_indices1, anchor_tb$accession1]
mouse_null = mm10_rna_matched[sample_indices2, anchor_tb$accession2]
gene_cor_null = map_dbl(1:nrow(gene_tb), function(j) {
        x = human_null[j, ]
        y = mouse_null[j, ]
        if (var(x) > 0 & var(y) > 0) {
                wCorr::weightedCorr(x,
                                    y,
                                    weights = anchor_tb$score,
                                    method = "spearman")
        } else {
                return(NA)
        }
})
spearman_pval = qvalue::empPvals(gene_cor, gene_cor_null)
spearman_pval_adj = p.adjust(spearman_pval, method = "BH")
# end computing spearman

# Updated version of computing hval (CO-B) scores
anchor_tb = read_tsv("./output/rna_conservation/anchors.tsv")
prob_vec = seq(0.1, 1.0, by = 0.01)
human_quantile_vec = quantile(hg38_rna_matched[, anchor_tb$accession1], prob_vec)
mouse_quantile_vec = quantile(mm10_rna_matched[, anchor_tb$accession2], prob_vec)
loci_hval = compute_hval(hg38_rna_matched[, anchor_tb$accession1],
                         mm10_rna_matched[, anchor_tb$accession2],
                         human_quantile_vec,
                         mouse_quantile_vec,
                         prob_vec)
null_size = 5e5
sample_indices1 = sample(nrow(hg38_rna_matched), null_size, replace = T)
sample_indices2 = sample(nrow(hg38_rna_matched), null_size, replace = T)
loci_hval_null = compute_hval(hg38_rna_matched[sample_indices1, anchor_tb$accession1],
                              mm10_rna_matched[sample_indices2, anchor_tb$accession2],
                              human_quantile_vec,
                              mouse_quantile_vec,
                              prob_vec)
hval_pval = qvalue::empPvals(loci_hval, loci_hval_null)
hval_pval_adj = p.adjust(hval_pval, method = "BH")

gene_tb = readRDS("./output/rna_conservation/rna_conservation_scores.rds")
gene_tb$gene_hval = loci_hval
gene_tb$gene_hval_adj = hval_pval_adj
gene_tb$high_gene_hval = hval_pval_adj < 0.1
(tibble(score = c(gene_tb$gene_hval,
                  loci_hval_null),
        source = c(rep("ZAlt", nrow(gene_tb)),
                   rep("Null", length(loci_hval_null)))) %>%
                ggdensity(x = "score", color = "source", fill = "source", xlab = "CO-B-GE", xlim = c(0.0, 1)) + 
                geom_vline(xintercept = min(gene_tb$gene_hval[gene_tb$high_gene_hval]), col = "#005792", size = 1) +
                scale_color_manual(values = c("#868686FF", "#ffd56b"))+
                scale_fill_manual(values = c("#868686FF", "#ffd56b")) + theme(text = element_text(size = 12), legend.position = "none")) %>%
        push_pdf("hist_gene_hval", w = 3, h = 2.5)
saveRDS(gene_tb, "./output/rna_conservation/rna_conservation_scores.rds")
push_tsv(gene_tb, file = "genes_summary")
# end update version hval

gene_tb = tibble(hid = rownames(hg38_rna_obj),
                 gene_spearman = gene_cor,
                 human_gene_symbols = human_symbols,
                 mouse_gene_symbols = mouse_symbols,
                 gene_spearman_pval_adj = spearman_pval_adj,
                 high_gene_spearman = gene_spearman_pval_adj < 0.1,
                 gene_hk_pct = hk_pct,
                 gene_hk_pval_adj = hk_pval_adj,
                 high_gene_housekeeping = hk_pval_adj < 0.1,
                 human_gene_mean_anchors = human_means,
                 mouse_gene_mean_anchors = mouse_means,
                 human_gene_var_std_anchors = human_var_std,
                 mouse_gene_var_std_anchors = mouse_var_std)
sum(gene_tb$high_gene_spearman)
sum(gene_tb$high_gene_housekeeping)

# plotting CO-V histogram
(tibble(score = c(gene_tb$gene_spearman,
                 gene_cor_null),
       source = c(rep("ZAlt", nrow(gene_tb)),
                  rep("Null", length(gene_cor_null)))) %>%
        ggdensity(x = "score", color = "source", fill = "source", xlab = "Gene Expr Conserved Var Score") + 
        geom_vline(xintercept = min(gene_tb$gene_spearman[gene_tb$high_gene_spearman]), col = "#005792", size = 1) +
        scale_color_manual(values = c("#868686FF", "#99154e"))+
        scale_fill_manual(values = c("#868686FF", "#99154e")) + theme(text = element_text(size = 12), legend.position = "none")) %>%
        push_pdf("hist_gene_spearman", w = 3, h = 2.5)

saveRDS(gene_tb, file = output_rds_path("rna_conservation_scores"))
push_tsv(gene_tb, file = "genes_summary")

# Figure 3b: Conserved vs non Conserved Heatmap
library(ComplexHeatmap)
set.seed(73)
val_col = circlize::colorRamp2(c(-2, 0, 2),
                               colors = c("blue", "white", "red"))
anchor_tb_sel = filter(anchor_tb, score > 0.5)
sample_indices0 = sample(which(gene_tb$high_gene_spearman &
                                      gene_tb$human_gene_var_std_anchors > 1.2 &
                                      gene_tb$mouse_gene_var_std_anchors > 1.2), size = 1000)
sample_indices1 = sample(which(gene_tb$gene_spearman < 0.2 &
                                      gene_tb$human_gene_var_std_anchors > 1.2 &
                                      gene_tb$mouse_gene_var_std_anchors > 1.2), size = 250)
heatmap_wrap <- function(sample_indices) {
        Heatmap(t(scale(t(log2(1+hg38_rna_matched[sample_indices, anchor_tb_sel$accession1])))),
                show_column_names = F,
                show_row_names = F,
                col = val_col,
                name = "Scaled Count",
                show_row_dend = F,
                show_column_dend = F) + 
                Heatmap(t(scale(t(log2(1+mm10_rna_matched[sample_indices, anchor_tb_sel$accession2])))),
                        show_column_names = F,
                        show_row_names = F,
                        col = val_col,
                        name = "Scaled Count",
                        show_row_dend = F,
                        show_column_dend = F)
}
heatmap_wrap(sample_indices0) %>% push_png(file_name = "heatmap_high_spearman", w = 4., h = 4., res = 1200)
heatmap_wrap(sample_indices1) %>% push_png(file_name = "heatmap_low_spearman", w = 4., h = 1., res = 1200)

# GO Analysis on conserved genes
library(enirchR)
select()
dbs <- c("GO_Biological_Process_2015")
enriched <- enrichr(gene_tb$human_gene_symbols[gene_tb$high_gene_spearman],
                    dbs)
plotEnrich(enriched[[1]],
           showTerms = 50,
           numChar = 40,
           y = "Count",
           orderBy = "P.value") %>%
        push_pdf(file_name = "geco_v_enirch", w = 6, h = 4.5)

enriched <- enrichr(gene_tb$human_gene_symbols[gene_tb$high_gene_housekeeping],
                    dbs)
plotEnrich(enriched[[1]],
           showTerms = 50,
           numChar = 40,
           y = "Count",
           orderBy = "P.value") %>%
        push_pdf(file_name = "geco_h_enirch", w = 6, h = 4.5)












