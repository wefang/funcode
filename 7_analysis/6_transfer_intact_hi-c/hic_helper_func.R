map_pair_mouse_index <- function(human_dhs_id, mapping_cat) {
        align_ind = mapping_cat %in% c("align", "align_sign", "align_sign_cov")
        align_index = paste0("align_", match(human_dhs_id[align_ind], regions_all$human_id))
        nonalign_motif_ind = mapping_cat == "noalign_sign_motif"
        if (sum(nonalign_motif_ind) == 0) {
                nonalign_motif_index = NULL
        } else {
                regions_noalign_sample = regions_noalign %>% 
                        filter(human_id %in% human_dhs_id[nonalign_motif_ind]) %>% 
                        filter(motif_ind) %>% 
                        group_by(human_id) %>% 
                        sample_n(size = 1) %>%
                        ungroup() %>% 
                        transmute(mouse_id, human_id)
                nonalign_motif_index = regions_noalign_sample$mouse_id[match(human_dhs_id[nonalign_motif_ind], regions_noalign_sample$human_id)]
        }
        nonalign_ind = mapping_cat == "noalign_sign"
        if (sum(nonalign_ind) == 0) {
                nonalign_index = NULL
        } else {
                regions_noalign_sample = regions_noalign %>% 
                        filter(human_id %in% human_dhs_id[nonalign_ind]) %>% 
                        group_by(human_id) %>% 
                        sample_n(size = 1) %>%
                        ungroup() %>% 
                        transmute(mouse_id, human_id)
                nonalign_index = regions_noalign_sample$mouse_id[match(human_dhs_id[nonalign_ind], regions_noalign_sample$human_id)]
        }
        out_index = character(length(human_dhs_id))
        out_index[align_ind] = align_index
        out_index[nonalign_motif_ind] = nonalign_motif_index
        out_index[nonalign_ind] = nonalign_index
        out_index
}
decide_dhs_category <- function(human_dhs_id) {
        out = rep("other", length(human_dhs_id))
        out[human_dhs_id %in% regions_noalign$human_id] = "noalign_sign"
        out[human_dhs_id %in% regions_noalign$human_id[regions_noalign$n_motif_shared >= 1]] = "noalign_sign_motif"
        out[human_dhs_id %in% regions_all$human_id] = "align"
        out[human_dhs_id %in% regions_all$human_id[regions_all$cov_any | regions_all$cob_any]] = "align_sign"
        out[human_dhs_id %in% regions_all$human_id[regions_all$cov_any]] = "align_sign_cov"
        out
}
check_loop <- function(input_region1, input_region2, hic_loop_first, hic_loop_second) {
        first_ol = findOverlaps(input_region1,
                                hic_loop_first)
        second_ol = findOverlaps(input_region2,
                                 hic_loop_second)
        first_split = split(subjectHits(first_ol), queryHits(first_ol))
        second_split = split(subjectHits(second_ol), queryHits(second_ol))
        
        pair_names = intersect(names(first_split), names(second_split))
        pair_ind = map_lgl(pair_names, function(x) {
                length(intersect(first_split[[x]], second_split[[x]])) > 0
        })
        
        first_ol1 = findOverlaps(input_region2,
                                 hic_loop_first)
        second_ol1 = findOverlaps(input_region1,
                                  hic_loop_second)
        first_split1 = split(subjectHits(first_ol1), queryHits(first_ol1))
        second_split1 = split(subjectHits(second_ol1), queryHits(second_ol1))
        
        pair_names1 = intersect(names(first_split1), names(second_split1))
        pair_ind1 = map_lgl(pair_names1, function(x) {
                length(intersect(first_split1[[x]], second_split1[[x]])) > 0
        })
        out_ind = logical(length(input_region1))
        out_ind[as.numeric(pair_names)[pair_ind]] = T
        out_ind[as.numeric(pair_names1)[pair_ind1]] = T
        out_ind
}
calculateDistances <- function(gr1, gr2) {
        # Ensure the two GRanges objects are of the same length
        if (length(gr1) != length(gr2)) {
                stop("gr1 and gr2 must be of the same length")
        }
        
        # Calculate distances
        startDists <- abs(start(gr1) - start(gr2))
        endDists <- abs(end(gr1) - end(gr2))
        closestStart = pmin(startDists, abs(start(gr1) - end(gr2)))
        closestEnd = pmin(endDists, abs(end(gr1) - start(gr2)))
        closestDist = pmin(closestStart, closestEnd)
        
        return(closestDist)
}