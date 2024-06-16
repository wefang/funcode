# write all bigwig to single folder and make a hub dir automatically
# load UCSC browser automatically
# wirte bed files for ENCODE DCC
library(GenomicRanges)
library(dplyr)
library(purrr)
library(rtracklayer)
library(stringr)
library(readr)
source("./helper/helper.R")
source("./6_produce_score_tracks/score_track_helper.R")
script_output_dir = "./output/ENCODE4_score_tracks_final/"
if (!dir.exists(script_output_dir)) {
        dir.create(script_output_dir)
}
regions = readRDS("./output/indexed_dhs_mapped_regions.rds")
regions_mm10 =readRDS("./output/indexed_mouse_dhs_mapped_regions.rds")
load("./intermediate_data/mouse_DHS_alignment_overlap.rda")
indices_add = which(filtering_df$sum[match(regions_mm10$id, paste0("MouseDHS-", filtering_df$identifier))] < 2)
regions_mm10_fil = regions_mm10[indices_add, ]
common_names = intersect(names(regions), names(regions_mm10_fil))
regions_all = bind_rows(dplyr::rename(regions[common_names], human_id = id),
                        dplyr::rename(regions_mm10_fil[common_names], mouse_id = id))
print(load("./metadata_processed//mm10_hg38_regions_cleaned.rda"))
hg38_regions_add = hg38_regions[indices_add]
mm10_regions_add = mm10_regions[indices_add]
print(load("metadata_processed/indexed_aligned_regions.rda"))
hg38_regions = c(hg38_regions, hg38_regions_add)
mm10_regions = c(mm10_regions, mm10_regions_add)
save(hg38_regions, mm10_regions, file = "metadata_processed/indexed_aligned_combined.rda")
regions_all = readRDS("./output/indexed_dhs_mapped_regions_combined.rds")
load("metadata_processed/indexed_aligned_combined.rda")
# output all tracks
score_data_dir = "./ucsc_tracks/data_v1/"
all_score_files = list()
for (type in c("cov")) {
# for (type in c("cob")) {
        message(type)
        # all_score_files[[type]] = list()
        for (assay in all_assays) {
                message(assay)
                # all_score_files[[type]][[assay]] = list()
                for (genome in c("hg38", "mm10")) {
                        message(genome)
                        bw_file = paste0(score_data_dir, paste0(type, "_", assay_name_mapping[assay], "_", genome, ".bw"))
                        all_score_files[[type]][[assay]][[genome]] = bw_file
                        if (type == "cob") {
                                track_obj = produce_score_track((regions_all[[paste0(type, "_", assay_name_mapping[assay])]])^5, genome = genome)
                        }
                        if (type == "cov") {
                                data_vec = regions_all[[paste0(type, "_", assay_name_mapping[assay])]]
                                track_obj = produce_score_track( sign(data_vec)*(data_vec)^2, genome = genome)
                        }
                        rtracklayer::export(track_obj, bw_file)
                }
        }
}
# Saving bed files
bed_hg38_tb = as_tibble(str_match(regions_all$hg38_region, "(chr.+):(\\d+)-(\\d+)")[, 2:4])
bed_hg38_tb$id = regions_all$human_id
bed_mm10_tb = as_tibble(str_match(regions_all$mm10_region, "(chr.+):(\\d+)-(\\d+)")[, 2:4])
bed_mm10_tb$id = regions_all$mouse_id
bed_data_dir = "./ucsc_tracks/bed_files/"
all_score_files = list()
for (type in c("cov", "cob")) {
        message(type)
        all_score_files[[type]] = list()
        for (assay in all_assays) {
                message(assay)
                all_score_files[[type]][[assay]] = list()
                for (genome in c("hg38", "mm10")) {
                        message(genome)

                        bed_file = paste0(bed_data_dir, paste0(type, "_", assay_name_mapping[assay], "_", genome, ".bed"))
                        all_score_files[[type]][[assay]][[genome]] = bed_file
                        
                        if (genome == "hg38") {
                                bed_tb_score = mutate(bed_hg38_tb, V2 = as.numeric(V2), V3 = as.numeric(V3))
                        }
                        if (genome == "mm10") {
                                bed_tb_score = mutate(bed_mm10_tb, V2 = as.numeric(V2), V3 = as.numeric(V3))
                        }
                        bed_tb_score$score = regions_all[[paste0(type, "_", assay_name_mapping[assay])]]
                        bed_tb_score$pval_adj = regions_all[[paste0(type, "_", assay_name_mapping[assay], "_fdr")]]
                        bed_tb_score = bed_tb_score[-1, ]
                        
                        write_tsv(bed_tb_score, bed_file, col_names = F)
                        R.utils::gzip(bed_file, overwrite = T)
                }
        }
}
# Madking a data hub for UCSC
base_url = "https://github.com/wefang/funcode/raw/main/track_hubs/data/"
base_dir = "./ucsc_tracks/"
hub_name = "co_mm10_hub"
genome = "mm10"
default_pos = "chr20:43879226-44106558"
hub_dir = paste0(base_dir, hub_name, "/")
score_data_dir = paste0(base_dir, "data/")
genome_dir = paste0(base_dir, hub_name, "/", genome, "/")
dir.create(hub_dir)
dir.create(genome_dir)
sink(paste0(hub_dir, "genomes.txt"))
cat(paste0("genome ", genome, "\n",
           "trackDb ", genome, "/trackDb.txt\n",
           "defaultPos ", default_pos, "\n"))
sink()
sink(paste0(hub_dir, "hub.txt"))
cat(paste0("hub ", genome, " FUNCODE scores
shortLabel ", genome, " FUNCODE scores
longLabel ", genome, " FUNCODE scores
genomesFile genomes.txt
email wfang9@jh.edu"))
sink()
sink(paste0(genome_dir, "trackDb.txt"))
for (type in c("cov", "cob")) {
        for (assay in all_assays) {
                cat_track_str(assay, type, genome,
                              url = paste0(base_url, basename(all_score_files[[type]][[assay]][[genome]])))
        }
}
sink()
hub_url = paste0("https://raw.githubusercontent.com/wefang/funcode/main/track_hubs/",
                 hub_name, "/hub.txt")
print(hub_url)