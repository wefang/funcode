#' produce score track
#' # requires latest GenomeInfoDB
produce_score_track <- function(val, genome = c("hg38", "mm10")) {
        fcon_track = switch(genome,
                            hg38 = hg38_regions,
                            mm10 = mm10_regions)
        fcon_track$score = val
        fcon_track = fcon_track[!is.na(fcon_track$score)]
        fcon_disjoin = disjoin(fcon_track)
        fcon_ol <- findOverlaps(fcon_track, fcon_disjoin)
        score(fcon_disjoin) <- aggregate(score(fcon_track[queryHits(fcon_ol)]), list(subjectHits(fcon_ol)), max)$x
        genome(fcon_disjoin) = genome
        seq_info = GenomeInfoDb::getChromInfoFromUCSC(genome)
        seqlengths(fcon_disjoin) = seq_info$size[match(names(seqlengths(fcon_disjoin)), seq_info$chrom)]
        fcon_disjoin
}
cat_track_str <- function(assay, type, genome, url) {
        score_label = paste0(type, "_", assay, "_", genome)
        cat(paste0("track ", tolower(score_label))); cat("\n")
        cat("type bigWig"); cat("\n")
        cat(paste0("shortLabel ", score_label, " score")); cat("\n")
        cat(paste0("longLabel ",  genome, " ", assay_name_mapping[assay], " conservation score")); cat("\n")
        cat("visibility full"); cat("\n")
        if (type == "cov") {
                cat("color 230,97,1"); cat("\n")
        } else {
                cat("color 94,60,153"); cat("\n")
        }
        cat("windowingFunction mean "); cat("\n")
        cat("autoScale off"); cat("\n")
        if (type == "cov") {
                cat("viewLimits -0.2:0.6"); cat("\n")
        } else {
                cat("viewLimits 0:0.9"); cat("\n")
        }
        cat("graphType bar"); cat("\n")
        cat("maxHeightPixels 50"); cat("\n")
        cat(paste0("bigDataUrl ", url)); cat("\n")
        cat("\n")
}