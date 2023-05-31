library(GenomicRanges)

# region set 1
# load("/dcl01/hongkai/data/wfang//shared/encode_compiled/Sep20/hg38_mm10_regions.rda")
# region_set = "human_mouse_mapped"
# input_regions_hg38 = hg38_regions
# input_regions_mm10 = mm10_regions
# end region set 1
# regions set 2
load("../input/ENCODE_human_mouse_DHS_noY.rda")
region_set = "human_mouse_DHS" # name of region set
input_regions_hg38 = human_regions # region
input_regions_mm10 = mouse_regions # region
# end region set 2

output_dir = paste0("../output/", region_set, "_matrices/")

# all_assays = c("DNase", "ATAC", "H3K27ac", "H3K27me3",
# 	"H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9ac", "H3K9me3", "Control")
all_assays = c("DNase", "ATAC", "H3K27ac", "H3K4me1", "H3K4me3", "Control")

# for (org in c("Homo sapiens", "Mus musculus")) {
for (org in c("Homo sapiens")) {
	message(org)
	if (org == "Homo sapiens") {
		region_obj = input_regions_hg38	
	}
	if (org == "Mus musculus") {
		region_obj = input_regions_mm10
	}
	for (assay in all_assays) {
		# assay = "Control"
		message(assay)
		data_name = paste0(stringr::str_replace_all(org, " ", "_"), "_", assay)
		data_counts = readRDS(paste0(output_dir, data_name, "_counts.rds"))
		data_libsize = readRDS(paste0(output_dir, data_name, "_libsize_correct.rds"))
		if (ncol(data_counts) == length(data_libsize)) {
							data_norm = t(t(data_counts) / data_libsize * 1e8)
		} else {
			stopifnot(nrow(data_counts) == length(data_libsize))
			data_norm = t(data_counts / data_libsize * 1e8)
		}
		assertthat::assert_that(nrow(data_norm) == length(region_obj))
		rownames(data_norm) = region_obj$identifier
		saveRDS(data_norm, paste0(output_dir, data_name, "_norm_correct.rds"))
		gc()
	}
}

