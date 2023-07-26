# cnsv score functions: weighted spearman and hval
library(wCorr)
library(matrixStats)
library(purrr)

# weighted spearman
weightedSpearman <- function(human_dhs_id,mouse_dhs_id){
    x = human_gene_mat[human_dhs_id, anchors_sel$accession1]
    y = mouse_gene_mat[mouse_dhs_id, anchors_sel$accession2]
    if (var(x) > 0 & var(y) > 0) {
        wCorr::weightedCorr(x,
                            y,
                            weights = anchors_sel$score,
                            method = "spearman")
    } else {
        return(NA)
    }

}

# weighted housekeep
WeightedHouseKeep <- function(human_dhs_id,mouse_dhs_id){
    hk_pct_human = matrixStats::rowWeightedMeans(1*(human_gene_mat[human_dhs_id, anchors_sel$accession1, drop = F] > human_cutoff),
                                                 anchors_sel$score)
    hk_pct_mouse = matrixStats::rowWeightedMeans(1*(mouse_gene_mat[mouse_dhs_id, anchors_sel$accession2, drop = F] > mouse_cutoff),
                                                 anchors_sel$score)
    hk_pct = pmin(hk_pct_human, hk_pct_mouse)
    return(hk_pct)

}
# new: hval
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
