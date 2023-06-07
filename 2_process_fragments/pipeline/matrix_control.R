library(tidyverse)
read_bin_file <- function(bin_file) {
	con = file(bin_file, open = "rb")
	cts = readBin(con, integer(), 1e8, size = 4)
	close(con)
	cts
}

# input
# data_files = readRDS("../output/all_alignment_files_metadata.rds")
data_files = readRDS("../input/26March22_control_ChIPseq_files_filtered.rds")
region_set = "human_mouse_mapped"
# region_set = "human_mouse_DHS"
# end input
libsize_dir = "../data/libsize_correct/"
# data_files$file_exists = file.exists(paste0("../data/", region_set, "/bin_correct/", data_files$`File accession`, ".bin"))
data_files$file_exists = file.exists(paste0("../data/", region_set, "/bin/", data_files$`File accession`, ".bin"))
table(data_files$file_exists)

output_dir = paste0("../output/", region_set, "_matrices/")
assay = "Control"
if (!dir.exists(output_dir)) dir.create(output_dir)
for (org in c("Mus musculus", "Homo sapiens")) {
    message(org)
    message(assay)
    data_files_sel = data_files[data_files$`Biosample organism` == org, ]
    bin_files = paste0("../data/", region_set, "/bin/", data_files_sel$`File accession`, ".bin")
	assertthat::assert_that(all(file.exists(bin_files)))
	print(length(bin_files))
	cts_mat = do.call(cbind, lapply(1:length(bin_files), function(i) {
        x = bin_files[i]
    	message(paste0(i, ": ", x))
        	read_bin_file(x)
            	}))
    cts_exp = ghelper::sumMatFac(t(cts_mat), data_files_sel$`Experiment accession`)
    saveRDS(cts_exp, file = paste0(output_dir, str_replace_all(org, " ", "_"), "_", assay, "_counts.rds"))

    libsize_files = paste0(libsize_dir, data_files_sel$`File accession`, ".rds")
    stopifnot(all(file.exists(libsize_files)))
    libsize = sapply(libsize_files, readRDS)
    libsize_exp = map_dbl(split(libsize, data_files_sel$`Experiment accession`), function(x) {
    	sum(x)
        	})
    all(names(libsize_exp) == colnames(cts_exp))
    saveRDS(libsize_exp, file = paste0(output_dir, str_replace_all(org, " ", "_"), "_", assay, "_libsize_correct.rds"))
}
