# TODO: check if color is actually used
assay_col = c(DNase = "#778868",
              ATAC = "#de8a5a",
              "DNase-seq" = "#778868",
              "ATAC-seq" = "#de8a5a")
assay_col_hu = c("DNase-seq" = "#70a494",
                 "ATAC-seq" = "#008080")
assay_col_mm = c("DNase-seq" = "#de8a5a",
                 "ATAC-seq" = "#ca562c")

# DHS component colors
comp_col_df = readxl::read_excel("./custom_metadata/dhs_component_color.xlsx", col_names = F)
comp_col = comp_col_df$...2
names(comp_col) = comp_col_df$...1

library(RColorBrewer)
cCRE_col = c("non-cCRE" = "#808080",
             "CTCF-only" = "#00b0f0",
             "DNase-H3K4me3" = "#ffa8a8",
             "dELS" = "#ffcc00",
             "pELS" = "#ffa600",
             "PLS" = "#ff0000"
             )

data_mod_col = c("dnase_atac" = "#06DA93",
                 "cov_chromatin_accessibility" = "#06DA93",
                 "cov_manual_chromatin_accessiblity" = "#7AD5C7",
                 "H3K4me3" = "#ff0000",
                 "cov_H3K4me3" = "#ff0000",
                 "H3K4me2" = "#FF7272",
                 "H3K9ac" = "#FFB5B5",
                 "H3K27ac" = "#ffcc00",
                 "cov_H3K27ac" = "#ffcc00",
                 "H3K4me1" = "#ffa600",
                 "cov_H3K4me1" = "#ffa600",
                 "IdentBasePct" = "#66806A",
                 "percentage" = "#66806A",
                 "phastCons4way" = "#B4C6A6",
                 "Average" = "#1C6DD0",
                 "cov_mean3" = "#1C6DD0",
                 "lecif" = "#d9d9d9",
                 "phastCons100way" = "#bdbdbd",
                 "phyloP4way" = "#868686",
                 "phyloP100way" = "#989898")
# method_col = c("CACO-V" = "#4daf4a",
#                "CACO-V_Manual" = "#ffcc29",
#                "OC Cnsv Var" = "#4daf4a",
#                "OC Cnsv Var Manual" = "#ffcc29",
#                "OC House Keeping Cnsv" = "#e41a1c",
#                "OC House Keeping Human" = "#b4aee8",
#                "OC House Keeping Mouse" = "#ffc478",
#                "LECIF" = "#086E7D",
#                "PhastCons4Way" = "#377eb8",
#                "PhyloP4Way" = "#EDCDBB",
#                "BasePercent" = "#ead3cb",
#                "Shared" = "#7897AB")
method_col = c("CACO-V" = "#9a1c4f",
               "OC Cnsv Var" = "#9a1c4f",
               "CACO-V_Manual" = "#d6604d",
               "OC Cnsv Var Manual" = "#d6604d",
               "OC House Keeping Cnsv" = "#ffd56b",
               "CACO-B" = "#ffd56b",
               "OC House Keeping Human" = "#b4aee8",
               "OC House Keeping Mouse" = "#ffc478",
               "LECIF" = "#984ea3",
               "PhastCons4Way" = "#377eb8",
               "PhyloP4Way" = "#a6cee3",
               "BasePercent" = "#139487",
               "Shared" = "#868686")
method_col_v1 = c("CACO-V" = "#d73027",
                  "OC Cnsv Var" = "#d73027",
                  "CACO-V_Manual" = "#f0634c",
                  "OC Cnsv Var Manual" = "#f0634c",
                  "OC House Keeping Cnsv" = "#1a9850",
                  "CACO-B" = "#1a9850",
                  "OC House Keeping Human" = "#33c24d",
                  "OC House Keeping Mouse" = "#5eda67",
                  "LECIF" = "#d9d9d9",
                  "PhastCons4Way" = "#bdbdbd",
                  "PhyloP4Way" = "#989898",
                  "BasePercent" = "#737373",
                  "Shared" = "#868686")

method_mapper = c("phastCons4way" = "PhastCons4Way",
                  "percentage" = "BasePercent",
                  "lecif" = "LECIF",
                  "phyloP4way" = "PhyloP4Way",
                  "spearman_global" = "CACO-V",
                  "spearman_manual" = "CACO-V_Manual",
                  "hk_pct" = "CACO-B")


cnsv_lvl_col = c("non-alignable" = "#C8C2BC",
                 "non-cCRE" = "#BDC7C9",
                 "GRCh38-mm10 Alignable" = "#BDC7C9",
                 "GRCh38 only cCRE" = "#3C8DAD",
                 "mm10 only cCRE" = "#5E8B7E",
                 "GRCh38 & mm10 cCRE" = "#926E6F"
)
chromatin_col = c("human_ATAC" = "#008080",
                  "human_DNase" = "#70a494",
                  "mouse_ATAC" = "#ca562c",
                  "mouse_DNase" = "#de8a5a")
human_mouse_col = c("Human" = "#008080", "Mouse" = "#ca562c")

cnsv_col = c("Repurposed" = brewer.pal(9, "Blues")[5],
             "HighCnsv" = brewer.pal(9, "Oranges")[5], "Other" = "#595959")
phast_cnsv_col = c("LowCnsv" = brewer.pal(9, "Purples")[5], "HighCnsv" = brewer.pal(9, "RdPu")[5], "Other" = "#595959")

caco_col = c("High CACO-V" = "#9B1C4F",
             "High CACO-B" = "#FFD56A",
             "Other" = "#868686")

c25 <- c(
        "dodgerblue2", "#E31A1C", # red
        "green4",
        "#6A3D9A", # purple
        "#FF7F00", # orange
        "black", "gold1",
        "skyblue2", "#FB9A99", # lt pink
        "palegreen2",
        "#CAB2D6", # lt purple
        "#FDBF6F", # lt orange
        "gray70", "khaki2",
        "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
        "darkturquoise", "green1", "yellow4", "yellow3",
        "darkorange4", "brown"
)
