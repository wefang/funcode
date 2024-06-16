# Summarize raw human fragments into read counts
library(rtracklayer)
library(tidyverse)
job_id = as.numeric(Sys.getenv("SGE_TASK_ID"))
# load regulatory element regions
load("../../shared/encode_compiled/Sep20/hg38_mm10_regions.rda")
# load("../../gencode_v24_promoters.rda")
target_regions = hg38_regions
# output_dir = "../gene_promoters/"
output_dir = "../mapped_regions/"
# files downloaded from CATLAS
cell_meta = readr::read_tsv("../meta.tsv")
sample_files = list.files("../fragments/", full = T)
print(length(sample_files))
i = job_id
sample_name = str_match(basename(sample_files[i]), "(.*)_rep1_fragments\\.bed\\.gz")[, 2]
output_file = paste0(output_dir, sample_name, ".rds")
if (!file.exists(output_file)) {
	frag = import.bed(sample_files[i])
	cellid = paste0(sample_name, "+", frag$name)
	cellid_sel = unique(cellid[cellid %in% cell_meta$cellID])

	cts_cells = map(cellid_sel, function(x) {
		message(x)
		frag_sel = frag[cellid %in% x]
		frag_sel_ctr = resize(frag_sel, fix = 'center', width = 1)
		cts = GenomicRanges::countOverlaps(target_regions, frag_sel_ctr)
		as(cbind(cts), "dgCMatrix")
	})
	cts_mat = do.call(cbind, cts_cells)
	colnames(cts_mat) = cellid_sel
	saveRDS(cts_mat, file = output_file)
}