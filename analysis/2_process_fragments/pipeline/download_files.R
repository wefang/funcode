options(timeout = max(1000, getOption("timeout")))
library(Rsamtools)
library(GenomicAlignments)
job_id = as.numeric(Sys.getenv("SGE_TASK_ID"))

# this script takes file metadata as input
# download all se and pe samples
# save fragments as .rds formats
# processed file metadata:
data_files = readr::read_tsv("../input/11Aug_hg38_mm10_DNase_ATAC_Histone_file_metadata.tsv")
data_files = data_files[data_files$`File assembly` %in% c("GRCh38", "mm10") & data_files$`Output type` == "alignments" & data_files$`File Status` == "released", ]
# parameters
job_batch_num = 20
data_dir = "../data/"
bam_dir = paste0(data_dir, "bam/") # temporary dir for bam files
target_se_dir = paste0(data_dir, "se_frag/")
target_pe_dir = paste0(data_dir, "pe_frag/")
# end parameters
job_indices = 1:nrow(data_files)
job_assign = rep(1:job_batch_num, 10000)[1:length(job_indices)]

dir.create(bam_dir, recursive = T)
dir.create(target_se_dir, recursive = T)
dir.create(target_pe_dir, recursive = T)
for (file_url in data_files$`File download URL`[job_indices[job_assign == job_id]]) {
	message(file_url)
	file_path = file.path(bam_dir, basename(file_url))
	output_files = c(se = file.path(target_se_dir, paste0(strsplit(basename(file_url), "\\.")[[1]][1], ".rds")),
			pe = file.path(target_pe_dir, paste0(strsplit(basename(file_url), "\\.")[[1]][1], ".rds")))
	if (any(file.exists(output_files))) {
		message('output exists')
		next
	}
	download.file(file_url, file_path)
	indexBam(file_path)
	is_pe = testPairedEndBam(file_path)
	if (!is_pe) {
		output_file = file.path(target_se_dir, paste0(strsplit(basename(file_url), "\\.")[[1]][1], ".rds"))
		if (file.exists(output_file)) {
			next
		} else {
			align = readGAlignments(file_path)
		}
	} else {
		output_file = file.path(target_pe_dir, paste0(strsplit(basename(file_url), "\\.")[[1]][1], ".rds"))
		align = readGAlignmentPairs(file_path)
	}
	saveRDS(align, file = output_file)
	file.remove(file_path)
}
