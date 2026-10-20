# Read the raw annotated Xenium object, prepare metadata, then recluster EC/PC and reconstruct vessels. 
# This script will take the seurat object processed through the Xenium Large Tissue Analysis scripts, subcluster EC/PC into Ec and PC, use DbScan taking Endothelial cells as anchor to detect Vessels assigning Vessel Ids to each individual vessel clusters.

library(Seurat)
library(SeuratObject)
library(harmony)
library(UCell)
library(dbscan)
library(dplyr)
library(ggplot2)
library(patchwork)

setwd("/Users/apirani/Documents/Data_projects/10X_Xenium_Analysis/NVU_Manuscript/NVU_Brain_Metastases/NVU_BBB_Manuscript/")

input_file <- "data/All.Big.Tissue.Xenium_Module_Scores.rds"
out <- "results"

# If you dont want to save each individual checkpoint rds files - make this FALSE
save_prepared_checkpoint <- TRUE

dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20260721)

fov_samples <- c(Lung_fov = "Lung", Breast_fov = "Breast", Normal_fov = "Normal",
                 Thyroid_fov = "Thyroid", Melanoma_fov = "Melanoma", Kidney_fov = "Kidney",
                 fov = "Gr_3", fov.3 = "GBM")
cohort_map <- c(Normal = "Normal", Gr_3 = "Primary", GBM = "Primary",
                Breast = "Mets", Kidney = "Mets", Lung = "Mets",
                Melanoma = "Mets", Thyroid = "Mets")
grade_map <- c(Normal = 0L, Gr_3 = 3L, GBM = 4L)


markers <- list(
  EC = c("PECAM1","CDH5","CLDN5","VWF","FLT1","KDR","TEK","ESAM",
         "PLVAP","ENG","NOS3","GJA4","NRP1","NOTCH1","ADGRL4","PODXL"),
  Mural = c("PDGFRB","RGS5","NOTCH3","KCNJ8"),
  Fibro = c("DCN","COL12A1","LAMA2","FBLN1","CSPG4","ABCC9")
)



obj <- readRDS(input_file)

Seurat::DefaultAssay(obj) <- "Xenium"

# Exactly the FOV membership lookup used in the original preparation step.
fov_rows <- lapply(names(fov_samples), function(fov) {
  if (!fov %in% names(obj@images)) {
    return(data.frame(fov = fov, sample = fov_samples[[fov]],
                      cell = character(), present = FALSE))
  }
  data.frame(fov = fov, sample = fov_samples[[fov]],
             cell = SeuratObject::Cells(obj[[fov]]), present = TRUE)
})
fov_cells <- dplyr::bind_rows(fov_rows)
missing_fovs <- setdiff(names(fov_samples), names(obj@images))
if (length(missing_fovs) > 0)
  warning("Configured FOVs missing from object: ", paste(missing_fovs, collapse = ", "))

# Fill ONLY missing sample labels, in the original metadata-row order.
cell_sample_map <- stats::setNames(fov_cells$sample, fov_cells$cell)

missing_sample <- is.na(obj@meta.data[["sample"]]) |
  !nzchar(as.character(obj@meta.data[["sample"]]))
obj@meta.data[["sample"]][missing_sample] <-
  unname(cell_sample_map[rownames(obj@meta.data)[missing_sample]])
obj@meta.data[["sample"]] <- as.character(obj@meta.data[["sample"]])

obj@meta.data$cohort <- unname(cohort_map[obj@meta.data$sample])
obj@meta.data$cohort[is.na(obj@meta.data$cohort)] <- "Unknown"
obj@meta.data$grade <- unname(grade_map[obj@meta.data$sample])


qc <- list(
  object = data.frame(n_cells = ncol(obj), n_genes = nrow(obj),
                      default_assay = Seurat::DefaultAssay(obj)),
  sample_counts = obj@meta.data |>
    dplyr::count(.data[["sample"]], name = "n_cells"),
  celltype_counts = obj@meta.data |>
    dplyr::count(.data[["sample"]], .data[["clusters"]], name = "n_cells"),
  fov_counts = fov_cells |>
    dplyr::count(fov, sample, present, name = "n_cells")
)

write.csv(qc$sample_counts, file.path(out, "01_sample_counts.csv"), row.names = FALSE)
write.csv(qc$celltype_counts, file.path(out, "01_celltype_counts.csv"), row.names = FALSE)
write.csv(qc$fov_counts, file.path(out, "01_fov_counts.csv"), row.names = FALSE)

if (save_prepared_checkpoint) {
  saveRDS(list(seurat_obj = obj, fov_cells = fov_cells, qc = qc),
          file.path(out, "01_prepared_input.rds"))
}


parent <- rownames(obj@meta.data)[as.character(obj$clusters) == "EC/PC"]
if (length(parent) < 20L) stop("Too few EC/PC cells")

ecpc <- obj
ecpc@graphs <- list()
ecpc@neighbors <- list()
ecpc@reductions <- list()
ecpc <- tryCatch(ecpc[, parent], error = function(e) subset(ecpc, cells = parent))
DefaultAssay(ecpc) <- "Xenium"
ecpc[["Xenium"]] <- tryCatch(JoinLayers(ecpc[["Xenium"]]),
                              error = function(e) ecpc[["Xenium"]])
ecpc <- ScaleData(ecpc, features = rownames(ecpc), verbose = FALSE)
ecpc <- RunPCA(ecpc, features = rownames(ecpc), npcs = 30, verbose = FALSE)
ecpc <- harmony::RunHarmony(ecpc, group.by.vars = "sample", verbose = FALSE)
ecpc <- FindNeighbors(ecpc, reduction = "harmony", dims = 1:20, verbose = FALSE)
ecpc <- FindClusters(ecpc, resolution = 0.30, verbose = FALSE)
ecpc <- RunUMAP(ecpc, reduction = "harmony", dims = 1:20, verbose = FALSE)

# Restore the original full marker panel and the original pre-mapping QC figures.
marker_sets <- lapply(markers, intersect, y = rownames(ecpc))
ecpc <- UCell::AddModuleScore_UCell(ecpc, features = marker_sets)
cluster_summary <- ecpc@meta.data %>% group_by(seurat_clusters) %>%
  summarise(n = n(), across(ends_with("_UCell"), mean, na.rm = TRUE), .groups = "drop")
write.csv(cluster_summary, file.path(out, "02_ECPC_cluster_marker_summary.csv"), row.names = FALSE)
write.csv(as.data.frame(table(ecpc$sample, ecpc$seurat_clusters)),
          file.path(out, "02_ECPC_sample_cluster_counts.csv"), row.names = FALSE)
p_umap <- DimPlot(ecpc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("EC/PC reclustering")
p_dot <- DotPlot(ecpc, features = unique(unlist(marker_sets))) + RotatedAxis()
ggsave(file.path(out, "02_ECPC_cluster_UMAP.pdf"), p_umap, width = 14, height = 12)
ggsave(file.path(out, "02_ECPC_marker_DotPlot.pdf"), p_dot, width = 14, height = 7)

cluster_map <- c(`0` = "Pericyte", `1` = "EC", `2` = "EC", `3` = "Pericyte",
                 `4` = "Pericyte", `5` = "LowQual_drop", `6` = "VLMC")


ecpc$vascular_subtype <- unname(cluster_map[as.character(ecpc$seurat_clusters)])
ecpc$EC_subtype <- ifelse(as.character(ecpc$seurat_clusters) == "2", "EC_brain",
                          ifelse(as.character(ecpc$seurat_clusters) == "1", "EC_brain_mets", NA_character_))
obj$vascular_subtype <- NA_character_
obj$EC_subtype <- NA_character_
obj@meta.data[colnames(ecpc), "vascular_subtype"] <- as.character(ecpc$vascular_subtype)
obj@meta.data[colnames(ecpc), "EC_subtype"] <- as.character(ecpc$EC_subtype)
obj$celltype_refined <- ifelse(as.character(obj$clusters) == "EC/PC",
                              as.character(obj$vascular_subtype), as.character(obj$clusters))


p1 <- DimPlot(ecpc, group.by = "seurat_clusters", label = TRUE)
p2 <- DimPlot(ecpc, group.by = "vascular_subtype", label = TRUE)
p3 <- as.data.frame(table(ecpc$vascular_subtype)) %>%
  ggplot(aes(x = Var1, y = Freq, fill = Var1)) +
  geom_col() + geom_text(aes(label = Freq), vjust = -0.4, size = 3.5) +
  labs(x = "Vascular Subtype", y = "Cell Count", fill = "Subtype") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
Idents(ecpc) <- "vascular_subtype"
vm <- FindAllMarkers(ecpc, only.pos = TRUE, min.pct = 0.25,
                     logfc.threshold = 0.25, test.use = "wilcox")
if (nrow(vm)) {
  top5 <- vm %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 5) %>%
    pull(gene) %>% unique()
  p7 <- DotPlot(ecpc, features = top5, group.by = "vascular_subtype") +
    RotatedAxis() + scale_color_gradient2(low = "blue", mid = "lightgrey", high = "red")
  p8 <- DoHeatmap(ecpc, features = top5, group.by = "vascular_subtype", angle = 45) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red")
  p10 <- DimPlot(obj, group.by = "celltype_refined", label = TRUE)
  ggsave(file.path(out, "03_ECPC_cluster_UMAP.pdf"),
         p1 + p2 + p3 + p7 + p8 + p10, width = 20, height = 12)
} else {
  stop("No vascular subtype markers detected")
}


obj$x_um <- NA_real_
obj$y_um <- NA_real_
for (fov in names(fov_samples)) {
  if (!fov %in% names(obj@images)) next
  cells <- intersect(Cells(obj[[fov]]), colnames(obj))
  obj@meta.data[cells, "fov"] <- fov
  xy <- as.data.frame(GetTissueCoordinates(obj, image = fov))
  if (!"cell" %in% colnames(xy)) xy$cell <- rownames(xy)
  rownames(xy) <- as.character(xy$cell)
  xcol <- intersect(c("x", "imagecol", "col", "center_x", "centroid_x"), colnames(xy))[1]
  ycol <- intersect(c("y", "imagerow", "row", "center_y", "centroid_y"), colnames(xy))[1]
  if (is.na(xcol) || is.na(ycol)) stop("Missing spatial coordinates for ", fov)
  xy <- xy[is.finite(xy[[xcol]]) & is.finite(xy[[ycol]]), , drop = FALSE]
  hit <- intersect(cells, rownames(xy))
  obj@meta.data[hit, "x_um"] <- as.numeric(xy[hit, xcol])
  obj@meta.data[hit, "y_um"] <- as.numeric(xy[hit, ycol])
}


obj$vessel_id <- NA_character_
summary <- list()
areas <- list()
for (fov in names(fov_samples)) {
  if (!fov %in% names(obj@images)) next

  xy <- as.data.frame(GetTissueCoordinates(obj, image = fov))
  if (!"cell" %in% colnames(xy)) xy$cell <- rownames(xy)
  rownames(xy) <- as.character(xy$cell)
  xcol <- intersect(c("x", "imagecol", "col", "center_x", "centroid_x"), colnames(xy))[1]
  ycol <- intersect(c("y", "imagerow", "row", "center_y", "centroid_y"), colnames(xy))[1]
  if (is.na(xcol) || is.na(ycol)) stop("Missing spatial coordinates for ", fov)
  xy <- xy[is.finite(xy[[xcol]]) & is.finite(xy[[ycol]]), , drop = FALSE]

  # Area uses ALL tissue cell coordinates, not only vessels.
  if (nrow(xy) >= 3L) {
    hull <- grDevices::chull(xy[[xcol]], xy[[ycol]])
    poly <- xy[c(hull, hull[1]), c(xcol, ycol), drop = FALSE]
    area_mm2 <- abs(sum(poly[-nrow(poly), 1] * poly[-1, 2] -
                        poly[-1, 1] * poly[-nrow(poly), 2])) / 2e6
  } else area_mm2 <- NA_real_
  areas[[fov]] <- data.frame(sample = fov_samples[[fov]], fov = fov,
                             area_mm2 = area_mm2)
  
  fov_cells <- intersect(Cells(obj[[fov]]), rownames(obj@meta.data))
  core <- fov_cells[obj@meta.data[fov_cells, "vascular_subtype"] %in% c("EC", "Pericyte")]
  core <- intersect(core, rownames(xy))
  if (length(core) < 5L) next  # exact vessel helper gate
  ec_core <- as.character(obj@meta.data[core, "vascular_subtype"]) == "EC"
  if (sum(ec_core) < 2L) next
  coords <- as.matrix(xy[core, c(xcol, ycol), drop = FALSE])
  ec_nn <- dbscan::kNN(coords[ec_core, , drop = FALSE], k = 1)$dist
  median_ec_nn_um <- stats::median(ec_nn, na.rm = TRUE)
  eps <- if (fov == "Normal_fov") 2 * median_ec_nn_um else 23.8
  if (!is.finite(eps) || eps <= 0) stop("Invalid DBSCAN eps for ", fov)

  clusters <- dbscan::dbscan(coords, eps = eps, minPts = 3L)$cluster
  # Retain only DBSCAN clusters containing at least one endothelial cell;
  # leave noise and pericyte-only clusters unassigned.
  ec_ok <- names(which(tapply(ec_core, clusters, any)))
  clusters[!(as.character(clusters) %in% ec_ok)] <- 0L
  ids <- ifelse(clusters == 0L, NA_character_,
                paste0(fov_samples[[fov]], "_V", clusters))
  obj@meta.data[core, "vessel_id"] <- ids

  summary[[fov]] <- data.frame(sample = fov_samples[[fov]], fov = fov,
     n_core = length(core), median_ec_nn_um = median_ec_nn_um, eps_um = eps,
     n_vessels = length(unique(stats::na.omit(ids))),
     pct_core_noise = 100 * mean(is.na(ids)))
}
vessel_summary <- dplyr::bind_rows(summary) %>%
  left_join(dplyr::bind_rows(areas), by = c("sample", "fov")) %>%
  mutate(vessels_per_mm2 = n_vessels / area_mm2)
write.csv(vessel_summary, file.path(out, "02_vessel_detection_summary.csv"), row.names = FALSE)
md <- obj@meta.data
assigned <- md %>% filter(vascular_subtype == "EC") %>% group_by(sample) %>%
  summarise(n_EC = n(), pct_EC_vessel_associated = 100 * mean(!is.na(vessel_id)),
            .groups = "drop")
write.csv(assigned, file.path(out, "02_EC_vessel_association.csv"), row.names = FALSE)
write.csv(md %>% filter(!is.na(vessel_id)) %>% group_by(sample, vessel_id) %>%
  summarise(n_EC = sum(vascular_subtype == "EC"),
            n_PC = sum(vascular_subtype == "Pericyte"), .groups = "drop"),
          file.path(out, "02_vessel_cell_counts.csv"), row.names = FALSE)


saveRDS(list(seurat_obj = obj), file.path(out, "01_vessels.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "01_R_sessionInfo.txt"))
message("Vessel-detection completed")
