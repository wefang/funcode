library(RColorBrewer)
# color for cCRE
cCRE_v4_col = c("non-cCRE" = "#808080",
                "CA" = "#b3cde3",
                "TF" = "#ffffcc",
                "CA-TF" = "#fed9a6",
                "CA-CTCF" = "#00b0f0",
                "CA-H3K4me3" = "#ffa8a8",
                "dELS" = "#ffcc00",
                "pELS" = "#ffa600",
                "PLS" = "#ff0000")
# color for conservation types
sign_type_col = c("non-alignable non-conserved" = "#a3a3a3",
                  "Alignable non-conserved" = "#868586",
                  "non-alignable CO-V" = "#a8cee1",
                  "non-alignable CO-V motif" = "#2178b4",
                  "Alignable CO-B" = "#ffce0b",
                  "Alignable CO-V" = "#9d1d50")

assay_col = c(DNase = "#778868",
              ATAC = "#de8a5a",
              "DNase-seq" = "#778868",
              "ATAC-seq" = "#de8a5a")
assay_col_hu = c("DNase-seq" = "#70a494",
                 "ATAC-seq" = "#008080")
assay_col_mm = c("DNase-seq" = "#de8a5a",
                 "ATAC-seq" = "#ca562c")

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
data_mod_col1 = c(chromatin = "#008080",
                  H3K4me3 = "#FF7272",
                  H3K27ac = "#D4C304",
                  H3K4me1 = "#4169E1")

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
co_col = c("High CO-V" = "#9B1C4F",
           "High CO-B" = "#FFD56A",
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

method_col_v2 = c("cov_mean" = "#a50026",
               "cov_manual_chromatin_accessiblity" = "#0c2c84",
               "cov_manual_H3K27ac" = "#225ea8",
               "cov_manual_H3K4me3" = "#1d91c0",
               "cov_manual_H3K4me1" = "#41b6c4",
               "cov_chromatin_accessibility" = "#d73027",
               "cov_H3K27ac" = "#f46d43",
               "cov_H3K4me3" = "#fdae61",
               "cov_H3K4me1" = "#fee08b",
               "cob_mean" = "#006837",
               "cob_chromatin_accessibility" = "#1a9850",
               "cob_H3K4me3" = "#66bd63",
               "cob_H3K27ac" = "#a6d96a",
               "cob_H3K4me1" = "#d9ef8b",
               "percentage"  = "#737373",
               "phastCons4way" = "#bdbdbd",
               "phyloP4way" = "#bdbdbd",
               "lecif" = "#d9d9d9")
