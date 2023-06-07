region_file = "./output/indexed_dhs_mapped_regions_v1.rds"
# region_file = "./output/indexed_dhs_mapped_regions.rds"
GRCh38_gencode_ref = "./metadata/ENCODEv4_Reference/ENCFF159KBI.gtf.gz"
mm10_gencode_ref = "./metadata/ENCODEv4_Reference/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz"

all_assays = c("dnase_atac", "H3K27ac", "H3K4me1","H3K4me3")
assay_name_mapping = c("dnase_atac" = "chromatin_accessibility",
                       "H3K27ac" ="H3K27ac",
                       "H3K4me1" = "H3K4me1",
                       "H3K4me3" = "H3K4me3")