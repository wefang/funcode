library(GenomicRanges)
library(GenomicAlignments)
library(ghelper)
library(tidyverse)
job_id = as.numeric(Sys.getenv("SGE_TASK_ID")) # total 20 batches
message(job_id)
dyn.load("~/R/x86_64-pc-linux-gnu-library/3.6/ghelper/libs/ghelper.so")

# input: regions to summarize counts over
load("/dcl01/hongkai/data/wfang//shared/encode_compiled/Sep20/hg38_mm10_regions.rda")
region_set = "human_mouse_mapped"
input_regions_hg38 = hg38_regions
input_regions_mm10 = mm10_regions
# end input
# input: regions to summarize counts over
# load("../input/ENCODE_human_mouse_DHS_noY.rda")
# region_set = "human_mouse_DHS" # name of region set
# input_regions_hg38 = human_regions # region
# input_regions_mm10 = mouse_regions # region
# end input
# input:
# load("../human_dhs_noY.rda")
# input_regions_hg38 = dhs_summit_extended
# end input
chrlen_df_hg38 = readr::read_tsv("/dcs04/hongkai/data/wfang/ENCODEv4_reference/GRCh38.sizes", col_names = F)
chrlen_df_hg38 = filter(chrlen_df_hg38, X1 %in% paste0("chr", c(1:22, "X")))
chrlen_df_mm10 = readr::read_tsv("/dcs04/hongkai/data/wfang/ENCODEv4_reference/mm10.sizes", col_names = F)
chrlen_df_mm10 = filter(chrlen_df_mm10, X1 %in% paste0("chr", c(1:20, "X")))

rds_files = c(list.files("../data/pe_frag/", full = T),
	list.files("../data/se_frag/", full = T))
file_accs = tools::file_path_sans_ext(basename(rds_files))

# all files
# data_files = readRDS("../output/all_alignment_files_metadata.rds")
# control files
data_files = readRDS("../input/26March22_control_ChIPseq_files_filtered.rds")

data_files_rds = rds_files[match(data_files$`File accession`, file_accs)]
ref_genome = data_files$`File assembly`

out_dir = paste0("../data/", region_set, "/bin/")
libsize_dir = "../data/libsize_correct/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
if (!dir.exists(libsize_dir)) dir.create(libsize_dir, recursive = T)

#### for resubmitting job ####
# out_files = paste0(out_dir, data_files$`File accession`, ".bin")
# libsize_files = paste0(libsize_dir, data_files$`File accession`, ".rds")
# out_files_exists = file.exists(out_files)
# libsize_files_exists = file.exists(libsize_files)
# run_indices = which(!out_files_exists | !libsize_files_exists)
# freeze run_indices 
# saveRDS(run_indices, file = "temp_process_run_indices.rds")
# run_indices = readRDS("temp_process_run_indices.rds")
#### end resubmitting ####

# original submission
run_indices = 1:length(data_files_rds)
# end original submission

set.seed(37)
job_indices = sample(1:20, length(run_indices), replace = T)
for (i in run_indices[which(job_indices == job_id)]) {
	message(data_files$`File accession`[i])
	if (ref_genome[i] == "GRCh38") {
		input_regions = input_regions_hg38
		chrlen_df = chrlen_df_hg38
	}
	if (ref_genome[i] == "mm10") {
		input_regions = input_regions_mm10
		chrlen_df = chrlen_df_mm10
	}
	if (is.na(data_files_rds[i])) {
		message('file not downloaded..')
		next
	}
	out_file = paste0(out_dir, data_files$`File accession`[i], ".bin")
	libsize_file = paste0(libsize_dir, data_files$`File accession`[i], ".rds")
	if (!file.exists(libsize_file) & file.exists(out_file)) {
		message('output libsize only')
		align = tryCatch({
			readRDS(data_files_rds[i])			
			}, error = function(condition) {
				message('file read failed..')
				})
		if(inherits(align, "error")) next
		libsize = sum(seqnames(align) %in% chrlen_df$X1)
		saveRDS(libsize, file = libsize_file)
	}
	if (!file.exists(out_file)) {
		message('processing')
		align = tryCatch({
			readRDS(data_files_rds[i])			
			}, error = function(condition) {
				message('file read failed..')
				})
		if(inherits(align, "error")) next
		if(is.null(align)) next
		if (!file.exists(libsize_file)) {
			libsize = sum(seqnames(align) %in% chrlen_df$X1)
			saveRDS(libsize, file = libsize_file)	
		}
		cts = align2region(align, input_regions, chrlen_df)
		writeBin(as.integer(cts), out_file, size = 4)		
	}
}

# data_files$`File accession` got named to file_acc instead
# trying fixing the bug by renaming files
# out_dir = paste0("../data/", region_set, "/bin/")
# libsize_dir = "../data/libsize/"
# out_dir_correct = paste0("../data/", region_set, "/bin_correct/")
# libsize_dir_correct ="../data/libsize_correct/"

# out_files_wrong =  paste0(out_dir, file_accs, ".bin")
# out_files_correct =  paste0(out_dir_correct, data_files$`File accession`, ".bin")
# libsize_file_wrong = paste0(libsize_dir, file_accs, ".rds")
# libsize_file_correct = paste0(libsize_dir_correct, data_files$`File accession`, ".rds")
# for (i in 1:length(out_files_wrong)) {
# 	file.rename(out_files_wrong[i], out_files_correct[i])
# 	file.rename(libsize_file_wrong[i], libsize_file_correct[i])
# }
