library(tidyverse)
read_bin_file <- function(bin_file) {
	con = file(bin_file, open = "rb")
	cts = readBin(con, integer(), 1e8, size = 4)
	close(con)
	cts
}

# input
data_files = readRDS("../output/all_alignment_files_metadata.rds")
# data_files = readRDS("../input/26March22_control_ChIPseq_files_filtered.rds")
region_set = "human_mouse_mapped"
# region_set = "human_mouse_DHS"
# end input
libsize_dir = "../data/libsize_correct/"
data_files$file_exists = file.exists(paste0("../data/", region_set, "/bin/", data_files$`File accession`, ".bin"))
# data_files$file_exists = file.exists(paste0("../data/", region_set, "/bin/", data_files$`File accession`, ".bin"))
table(data_files$file_exists)
# try correcting rds file read errors by downloading again. Delete those now to be processed again.
# data_dir = "../data/"
# target_se_dir = paste0(data_dir, "se_frag/")
# target_pe_dir = paste0(data_dir, "pe_frag/")
# walk(data_files$`File download URL`[!data_files$file_exists], function(file_url) {
# 	output_files = c(se = file.path(target_se_dir, paste0(strsplit(basename(file_url), "\\.")[[1]][1], ".rds")),
# 				pe = file.path(target_pe_dir, paste0(strsplit(basename(file_url), "\\.")[[1]][1], ".rds")))	
# 	print(file.exists(output_files))
# 	file.remove(output_files[file.exists(output_files)])
# 	})
# end correct failed downloads

data_exp = readr::read_tsv("../input/11Aug_hg38_mm10_DNase_ATAC_Histone_exp_metadata.tsv", skip = 1)
assertthat::assert_that(all(data_files$`Experiment accession` %in% data_exp$Accession))
data_files = filter(data_files, file_exists)
data_exp = filter(data_exp, Accession %in% data_files$`Experiment accession`)
data_exp = mutate(data_exp, Assay = paste0(data_exp$`Assay title`, "_", data_exp$`Target of assay`))

# criteria (histone target # of exp > 100, human or mouse)
histone_target_total = table(data_exp$`Target of assay`)
histone_targets = names(histone_target_total)[histone_target_total >= 100]

# table(data_exp$`Assay title`, data_exp$`Target of assay`)[, histone_targets]
# table(data_exp$Organism, data_exp$`Target of assay`)[, histone_targets]

histone_assays = map(histone_targets, function(x) paste0(c("Histone ChIP-seq", "Mint-ChIP-seq"), "_", x))
names(histone_assays) = histone_targets
all_assays = c(list("DNase" = paste0("DNase-seq", "_", NA),
	"ATAC" = paste0("ATAC-seq", "_", NA)),
histone_assays)
assertthat::assert_that(all(map_lgl(all_assays, function(x) any(x %in% data_exp$Assay))))

output_dir = paste0("../output/", region_set, "_matrices/")
if (!dir.exists(output_dir)) dir.create(output_dir)
for (org in c("Homo sapiens", "Mus musculus")) {
    message(org)
	for (assay in names(all_assays)) {
        message(assay)
		data_exp_sel = filter(data_exp, Organism == org & Assay %in% all_assays[[assay]])
		data_files_sel = filter(data_files, `Experiment accession` %in% data_exp_sel$Accession)
		# bin_files = paste0("../data/", region_set, "/bin_correct/", data_files_sel$`File accession`, ".bin")
  #       bin_files = paste0("../data/", region_set, "/bin/", data_files_sel$`File accession`, ".bin")
		# assertthat::assert_that(all(file.exists(bin_files)))
		# print(length(bin_files))
		cts_mat = do.call(cbind, lapply(bin_files, function(x) {
        	message(x)
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
    # all(names(libsize_exp) == colnames(cts_exp))
    saveRDS(libsize_exp, file = paste0(output_dir, str_replace_all(org, " ", "_"), "_", assay, "_libsize_correct.rds"))
    	}
    }

# org = "Homo sapiens"
# assay = "DNase"
# libsize_old = readRDS(paste0(output_dir, str_replace_all(org, " ", "_"), "_", assay, "_libsize.rds"))
# libsize_new = readRDS(paste0(output_dir, str_replace_all(org, " ", "_"), "_", assay, "_libsize_correct.rds"))