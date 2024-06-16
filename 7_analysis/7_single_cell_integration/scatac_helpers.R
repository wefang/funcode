matched_label = list("Cardiomyocytes" = c("Cam", "Cardiomyocytes"),
                     "Hepatocytes" = c("Hpc", "Hepatocytes"),
                             "Pneumocytes" = c("Pal", "Type I pneumocytes", "Type II pneumocytes"),
                     "Enterocytes" = c("Enc.1", "Enc.2", "Enc.3", "Enterocytes"),
                     "Endothelial cells" = c("End.1", "End.2", c(paste0("Fib.", 1:6)), "Endothelial II cells", "Endothelial I cells"),
                     # "Smm" = "Smm",
                     # "Gbl" = "Gbl.1",
                     "Macphages" = c("Macrophages", "Alveolar macrophages", "Mac.1", "Mac.2", "Mac.3"),
                     "T cells" = c("CD4", "CD8", "T cells", "Regulatory T cells", "Tly.1", "Tly.2"),
                     "B cells" = c("Bly", "B cells", "Activated B cells", "Immature B cells"))
matched_col = RColorBrewer::brewer.pal(12, name = "Paired")
# matched_col = circlize::rand_color(8, luminosity = "birght")
names(matched_col) = names(matched_label)

matched_label_col = purrr::reduce(map(1:length(matched_label), function(i) {
        out = rep(matched_col[i], length(matched_label[[i]]))
        names(out) = matched_label[[i]]
        out
}), c)
matched_label_mapper = purrr::reduce(map(names(matched_label), function(x) {
        query = matched_label[[x]]
        out = rep(x, length(query))
        names(out) = query
        out
}), c)

list_ren_samples <- function() {
        print(map(split(meta$sample, meta$tissue), unique))
}
load_ren_sample <- function(sample_name,
                            min.cells.pct = 0.001,
                            min.features = 500) {
        cts_mat = readRDS(paste0(path_data,"mapped_regions/", sample_name, "_mapped_regions.rds"))
        print(sort(table(meta$`cell type`[match(colnames(cts_mat), meta$cellID)])))
        # all(colnames(cts_mat) %in% meta$cellID)
        rownames(cts_mat) = regions$id
        cts_mat@x = rep(1, length(cts_mat@x))
        human_obj = CreateSeuratObject(cts_mat,
                                       min.cells = min.cells.pct * ncol(cts_mat),
                                       min.features = min.features,
                                       assay = "DHS")
        human_obj$cell_type = meta$`cell type`[match(colnames(human_obj), meta$cellID)]
        human_obj$sample = sample_name
        human_obj
}
# this version works with sparse matrix
aveMatFac <- function(mat, fac) {
        if (class(fac) != "factor") 
                fac <- factor(fac)
        rown <- length(levels(fac))
        coln <- dim(mat)[2]
        out <- matrix(, rown, coln)
        ind <- as.numeric(fac)
        for (i in 1:rown) {
                out[i, ] <- Matrix::colMeans(mat[ind == i, , drop = F], na.rm = T)
        }
        rownames(out) <- levels(fac)
        return(out)
}
sumMatFac <- function(mat, fac) {
        if (class(fac) != "factor") 
                fac <- factor(fac)
        rown <- length(levels(fac))
        coln <- dim(mat)[2]
        out <- matrix(, rown, coln)
        ind <- as.numeric(fac)
        for (i in 1:rown) {
                out[i, ] <- Matrix::colSums(mat[ind == i, , drop = F], na.rm = T)
        }
        rownames(out) <- levels(fac)
        return(out)
}
quick_process <- function(integrated, dims = 1:15, var_feature = NULL) {
        if (is.null(var_feature)) {
                var_feature = rownames(integrated)
        }
        VariableFeatures(integrated) = var_feature
        integrated <-
                ScaleData(integrated,
                        features = var_feature,
                        verbose = FALSE)
        integrated <-
                RunPCA(integrated, features = var_feature, npcs = 30)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = dims)
}
integrate_harmony <- function(dhs, nclust = 20, max.iter.harmony = 5, max.iter.cluster = 10) {
        human_mouse_obj = process_scatac(human_mouse_obj,
                                         features = dhs,
                                         normalize = T,
                                         run_umap = F)
        human_mouse_obj = L2Dim(human_mouse_obj, reduction = "lsi")        
        human_mouse_harmony = RunHarmony(human_mouse_obj,
                                         group.by.vars = c("species"), reduction = "lsi.l2",
                                         assay.use = "DHS",
                                         nclust = nclust,
                                         max.iter.harmony = max.iter.harmony,
                                         max.iter.cluster = max.iter.cluster,
                                         epsilon.cluster = -Inf,
                                         epsilon.harmony = -Inf,
                                         project.dim = FALSE)
        human_mouse_harmony
}
mm9_regions = readRDS(paste0("/dcs04/hongkai/data/rzhao/sciATAC/sciATAC_human_mouse_immune/mm9_mapped_regions.rds"))
mouse_cell_meta = readr::read_tsv(paste0("/dcs04/hongkai/data/rzhao/sciATAC/sciATAC_human_mouse_immune/mouse_data/cell_metadata.txt"))
list_mouse_samples <- function(){
        print(list.files(paste0("/dcs04/hongkai/data/rzhao/sciATAC/sciATAC_human_mouse_immune/mouse_data/", pattern = ".*_new\\.rds")))
}
load_mouse_sample <- function(sample_name, min.cells.pct = 0.001) {
        mouse_counts = readRDS(paste0("/dcs04/hongkai/data/rzhao/sciATAC/sciATAC_human_mouse_immune/mouse_data/", sample_name, "_new.rds"))
        mouse_cell_meta_sub = mouse_cell_meta[mouse_cell_meta$tissue.replicate == sample_name, ]
        rownames(mouse_counts) = mm9_regions$id
        colnames(mouse_counts) = mouse_cell_meta_sub$cell
        mouse_obj = CreateSeuratObject(mouse_counts,
                                       min.cells = ceiling(ncol(mouse_counts) * min.cells.pct),
                                       assay = "DHS")
        mouse_obj$cell_type = mouse_cell_meta_sub$cell_label[match(colnames(mouse_obj), mouse_cell_meta_sub$cell)]
        # mouse_obj = mouse_obj[, mouse_obj$cell_type %in% names(sort(table(mouse_obj$cell_type), decreasing = T))[1:5]]
        mouse_obj$sample = sample_name
        # mouse_obj = RunTFIDF(mouse_obj)
        mouse_obj = mouse_obj[, mouse_obj$cell_type != "Collisions"]
        mouse_obj
}
process_scatac <- function(data_obj, features = NULL, dims = 2:30, normalize = T, run_umap = T) {
        # data_obj = NormalizeData(data_obj)
        if (normalize) {
                data_obj = RunTFIDF(data_obj)                
        }
        data_obj = FindTopFeatures(data_obj, min.cutoff = 'q0')
        data_obj = RunSVD(data_obj, features = features)
        if (run_umap) {
                data_obj = RunUMAP(object = data_obj, reduction = 'lsi', dims = dims)
        }
        data_obj
}