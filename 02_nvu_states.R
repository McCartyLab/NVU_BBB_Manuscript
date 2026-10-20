# 02_nvu_states.R


library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(FNN)

out <- "results"
obj <- readRDS(file.path(out, "01_vessels.rds"))$seurat_obj
DefaultAssay(obj) <- "Xenium"

# Test BBB Stability Score by helding out "MFSD2A","SLCO1C1","LSR" genes. It will be BBB Indep Score
sets <- list(
  BBB_holdout = c("MFSD2A","SLCO1C1","LSR"),
  Tight_Junctions = c("CLDN5","OCLN","TJP1","TJP2","TJP3","F11R","JAM2","JAM3","LSR"),
  Tight_Junctions_indep = c("CLDN5","OCLN","TJP1","TJP2","TJP3","F11R","JAM2","JAM3"),
  BBB_Transport = c("ABCB1","ABCG2","SLC2A1","SLCO1A2","SLCO1C1","LRP1","MFSD2A"),
  BBB_Transport_indep = c("ABCB1","ABCG2","SLC2A1","SLCO1A2","LRP1"),
  Angiogenesis = c("ANGPT2","VEGFA","DLL4","NDP","KDR","TEK"),
  Vascular_Inflam = c("ICAM1","VCAM1","IL1B","IL6","CCL2","S1PR3"),
  ECM_Remodeling = c("COL4A1","LAMB1","FN1","VCL","ITGB1","MMP2","MMP9","VIM"),
  Stem_like = c("SOX2","SOX4","NES","HES1","NNAT"),
  Hypoxia = c("HILPDA","VEGFA","SLC2A1"),
  Mesenchymal = c("LOX","CAV1","ANXA1","FN1"),
  Invasive = c("BCAN","TTYH1","MMP2","MMP9"),
  Proliferative = c("TOP2A","MKI67","PCNA","EGFR"),
  Prolif_EC_QC = c("MKI67","TOP2A","CDK1","PCNA","CENPF","CCNB2","CCNA1"),
  AngioEC_QC = c("DLL4","KDR","FLT1","NRP1","ANGPT2","PLVAP","NOTCH1"),
  Ast_NonMalignant_QC = c("ALDH1L1","GFAP","GJA1","FGFR3","MLC1","SLC1A2","SLC1A3","SLC4A4","SLC38A3"),
  Ast_Malignant_QC = c("EGFR","BCAN","TTYH1","SOX2","SOX4","SOX11","NES","HES1","IGFBP3","IGFBP5",
                       "MKI67","TOP2A","PCNA","CDK1","CCNB2","CENPF")
)

# Just for safety drop signatures with fewer than two available genes - Nothing is dropped.
# delete old UCell score columns before calculating the new scores.
sets <- lapply(sets, intersect, y = rownames(obj))
found <- lengths(sets)
if (any(found < 2)) warning("Insufficient panel genes: ", paste(names(found)[found < 2], collapse = ", "))
sets <- sets[found >= 2]
write.csv(data.frame(signature = names(found), n_genes = unname(found)),
          file.path(out, "03_signature_gene_counts.csv"), row.names = FALSE)

# Score cells with AddModuleScore_UCell
stale <- intersect(paste0(names(sets), "_UCell"), colnames(obj@meta.data))
if (length(stale)) obj@meta.data[, stale] <- NULL
obj <- UCell::AddModuleScore_UCell(obj, features = sets)
obj$BBB_Stability_full <- with(obj@meta.data,
  Tight_Junctions_UCell + BBB_Transport_UCell -
    (Angiogenesis_UCell + Vascular_Inflam_UCell + ECM_Remodeling_UCell))
obj$BBB_Stability_indep <- with(obj@meta.data,
  Tight_Junctions_indep_UCell + BBB_Transport_indep_UCell -
    (Angiogenesis_UCell + Vascular_Inflam_UCell + ECM_Remodeling_UCell))

# Safety check: composite values outside resolved ECs must be NA (NOT an all-cell score).
non_ec <- is.na(obj$vascular_subtype) | obj$vascular_subtype != "EC"
obj@meta.data$BBB_Stability_full[non_ec] <- NA_real_
obj@meta.data$BBB_Stability_indep[non_ec] <- NA_real_

# # QC-only audit of broad Ast-like cells; NO exclusion from NVU coverage.
# md <- obj@meta.data
# thresholds <- list()
# for (s in unique(na.omit(md$sample))) {
#   same <- md$sample == s & !is.na(md$sample)
#   ast <- same & as.character(md$clusters) == "Ast-like"
#   tumor <- same & as.character(md$clusters) == "Tumor"
#   nonast <- same & as.character(md$clusters) != "Ast-like"
#   nm <- md$Ast_NonMalignant_QC_UCell
#   mg <- md$Ast_Malignant_QC_UCell
#   nm_cut <- if (sum(is.finite(nm[ast]), na.rm = TRUE) >= 10) quantile(nm[ast], .25, na.rm = TRUE) else
#     if (sum(is.finite(nm[same]), na.rm = TRUE) >= 10) quantile(nm[same], .75, na.rm = TRUE) else .02
#   mg_cut <- if (sum(is.finite(mg[tumor]), na.rm = TRUE) >= 10) quantile(mg[tumor], .75, na.rm = TRUE) else
#     if (sum(is.finite(mg[nonast]), na.rm = TRUE) >= 10) quantile(mg[nonast], .95, na.rm = TRUE) else .02
#   thresholds[[s]] <- data.frame(sample = s, ast_nonmal_cut = max(.02, nm_cut),
#                                 ast_malignant_cut = max(.02, mg_cut), n_ast = sum(ast, na.rm = TRUE), n_tumor = sum(tumor, na.rm = TRUE))
# }
# thresholds <- bind_rows(thresholds)
# write.csv(thresholds, file.path(out, "03_Ast_audit_thresholds.csv"), row.names = FALSE)
# md <- md %>% left_join(thresholds, by = "sample")
# obj$ast_like_audit_class <- NA_character_
# ast <- as.character(md$clusters) == "Ast-like"
# nonmal <- md$Ast_NonMalignant_QC_UCell >= md$ast_nonmal_cut &
#   md$Ast_NonMalignant_QC_UCell >= md$Ast_Malignant_QC_UCell - .02
# mal <- md$Ast_Malignant_QC_UCell >= md$ast_malignant_cut &
#   md$Ast_Malignant_QC_UCell > md$Ast_NonMalignant_QC_UCell + .02
# obj@meta.data$ast_like_audit_class[which(ast %in% TRUE)] <- "Uncertain_Ast_like"
# obj@meta.data$ast_like_audit_class[which(ast & mal %in% TRUE)] <- "Malignant_like_Ast_like"
# obj@meta.data$ast_like_audit_class[which(ast & nonmal %in% TRUE)] <- "Non_malignant_Ast_like"
# write.csv(as.data.frame(table(obj$sample, obj$ast_like_audit_class)),
#           file.path(out, "03_Ast_audit_summary.csv"), row.names = FALSE)

# Nearest vessel-ASSOCIATED EC
obj$dist_to_vessel <- NA_real_
obj$nearest_vessel <- NA_character_
fov_samples <- c(Lung_fov="Lung", Breast_fov="Breast", Normal_fov="Normal",
                 Thyroid_fov="Thyroid", Melanoma_fov="Melanoma", Kidney_fov="Kidney",
                 fov="Gr_3", fov.3="GBM")
for (fov in names(fov_samples)) {
  if (!fov %in% names(obj@images)) next
  xy <- as.data.frame(GetTissueCoordinates(obj, image = fov))
  if (!"cell" %in% colnames(xy)) xy$cell <- rownames(xy)
  rownames(xy) <- as.character(xy$cell)
  xcol <- intersect(c("x", "imagecol", "col", "center_x", "centroid_x"), colnames(xy))[1]
  ycol <- intersect(c("y", "imagerow", "row", "center_y", "centroid_y"), colnames(xy))[1]
  if (is.na(xcol) || is.na(ycol)) stop("Missing spatial coordinates for ", fov)
  cells_all <- intersect(Cells(obj[[fov]]), rownames(obj@meta.data))
  xy <- xy[intersect(cells_all, rownames(xy)), , drop = FALSE]
  xy <- xy[is.finite(xy[[xcol]]) & is.finite(xy[[ycol]]), , drop = FALSE]
  cells <- intersect(rownames(xy), rownames(obj@meta.data))
  if (!length(cells)) next
  md_fov <- obj@meta.data[cells, , drop = FALSE]
  ref <- cells[md_fov$vascular_subtype %in% "EC" & !is.na(md_fov$vessel_id)]
  if (!length(ref)) next
  qmat <- as.matrix(xy[cells, c(xcol, ycol), drop = FALSE])
  rmat <- as.matrix(xy[ref, c(xcol, ycol), drop = FALSE])
  nn <- FNN::get.knnx(rmat, qmat, k = 1)
  obj@meta.data[cells, "dist_to_vessel"] <- nn$nn.dist[, 1]
  obj@meta.data[cells, "nearest_vessel"] <-
    as.character(obj@meta.data[ref[nn$nn.index[, 1]], "vessel_id"])
}

# sample-specific radius QC: pericyte Q90, reference median,
# selected samples, 20-30 um perivascular radius and separate 15-um tumor-touch annotation.
md <- obj@meta.data
radius_raw <- md %>% filter(vascular_subtype == "Pericyte", is.finite(dist_to_vessel)) %>%
  group_by(sample) %>% summarise(n_pericyte = n(), pericyte_median_um = median(dist_to_vessel),
                                 raw_radius_um = as.numeric(quantile(dist_to_vessel, .90)), .groups = "drop")
parenchyma <- md %>% filter(as.character(clusters) %in% c("Neu-like", "OPC/OL-like"),
                           is.finite(dist_to_vessel)) %>% group_by(sample) %>%
  summarise(n_parenchyma = n(), parenchyma_median_um = median(dist_to_vessel), .groups = "drop")
reference <- c("GBM", "Gr_3", "Thyroid")
derived <- c("GBM", "Gr_3", "Thyroid", "Breast", "Melanoma")
consensus <- median(radius_raw$raw_radius_um[radius_raw$sample %in% reference], na.rm = TRUE)
if (!is.finite(consensus)) consensus <- median(radius_raw$raw_radius_um, na.rm = TRUE)
if (!is.finite(consensus)) consensus <- 25
sample_order <- c("Normal", "Gr_3", "GBM", "Breast", "Kidney", "Lung", "Melanoma", "Thyroid")
radius <- data.frame(sample = sample_order) %>% left_join(radius_raw, by = "sample") %>%
  left_join(parenchyma, by = "sample") %>%
  mutate(source = ifelse(sample %in% derived & is.finite(raw_radius_um), "data_derived", "consensus_fallback"),
         selected_radius_um = ifelse(source == "data_derived", raw_radius_um, consensus),
         periv_radius_um = pmin(pmax(selected_radius_um, 20), 30),
         separation_ratio = parenchyma_median_um / periv_radius_um)
write.csv(radius, file.path(out, "04_perivascular_radius_QC.csv"), row.names = FALSE)
obj$periv_radius_um <- radius$periv_radius_um[match(obj$sample, radius$sample)]

md <- obj@meta.data
core_ec <- !is.na(md$vessel_id) & md$vascular_subtype %in% "EC"
core_pc <- !is.na(md$vessel_id) & md$vascular_subtype %in% "Pericyte"
shell <- is.finite(md$dist_to_vessel) & is.finite(md$periv_radius_um) &
         md$dist_to_vessel > 0 & md$dist_to_vessel <= md$periv_radius_um
obj$nvu_role <- NA_character_
obj@meta.data$nvu_role[core_ec] <- "endothelial"
obj@meta.data$nvu_role[core_pc] <- "pericyte"
obj@meta.data$nvu_role[shell & !(core_ec | core_pc) & md$clusters == "Tumor" & md$sample != "Normal"] <- "perivascular-tumor"
obj@meta.data$nvu_role[shell & !(core_ec | core_pc) & md$clusters == "Ast-like"] <- "perivascular-astrocyte"
obj$nvu_unit_id <- NA_character_
obj@meta.data$nvu_unit_id[core_ec | core_pc] <- md$vessel_id[core_ec | core_pc]
recruited <- !is.na(obj$nvu_role) & is.na(md$vessel_id) & !is.na(md$nearest_vessel)
obj@meta.data$nvu_unit_id[recruited] <- md$nearest_vessel[recruited]
obj$in_NVU <- !is.na(obj$nvu_unit_id) & !is.na(obj$nvu_role)

unit <- obj@meta.data %>% filter(!is.na(nvu_unit_id)) %>% group_by(sample, nvu_unit_id) %>%
  summarise(n_EC = sum(nvu_role == "endothelial", na.rm = TRUE),
            n_PC = sum(nvu_role == "pericyte", na.rm = TRUE),
            n_astro = sum(nvu_role == "perivascular-astrocyte", na.rm = TRUE),
            n_tumor = sum(nvu_role == "perivascular-tumor", na.rm = TRUE),
            n_tumor_touching = sum(nvu_role == "perivascular-tumor" & dist_to_vessel <= 15, na.rm = TRUE),
            .groups = "drop") %>% filter(n_EC >= 1) %>%
  mutate(tumor_contact = n_tumor >= 1, tumor_touching = n_tumor_touching >= 1,
         pericyte_present = n_PC >= 1, astro_present = n_astro >= 1, pericyte_ratio = n_PC / n_EC,
         nvu8_id = case_when(
           pericyte_present & astro_present & !tumor_contact ~ "canonical_NVU_like",
           !pericyte_present & astro_present & !tumor_contact ~ "canonical_astro_no_pericyte",
           pericyte_present & astro_present & tumor_contact ~ "tumor_associated_intact_NVU",
           !pericyte_present & astro_present & tumor_contact ~ "tumor_contact_pericyte_dropout",
           pericyte_present & !astro_present & tumor_contact ~ "tumor_vessel_contact_no_astro",
           !pericyte_present & !astro_present & tumor_contact ~ "tumor_contact_bare_vessel",
           pericyte_present & !astro_present & !tumor_contact ~ "EC_PC_vessel_only",
           TRUE ~ "bare_EC_vessel"),
         nvu4_id = case_when(
           nvu8_id %in% c("canonical_NVU_like", "canonical_astro_no_pericyte") ~ "canonical_or_astro_NVU",
           nvu8_id == "tumor_associated_intact_NVU" ~ "tumor_intact_NVU",
           nvu8_id %in% c("tumor_contact_pericyte_dropout", "tumor_vessel_contact_no_astro",
                          "tumor_contact_bare_vessel") ~ "disrupted_tumor_vessel",
           TRUE ~ "vessel_only"))

# Final manuscript NVU names. Keep nvu8_id unchanged for the downstream analyses.
# These describe spatial association, not confirmed physical cellular coverage.
nvu_label <- c(
  canonical_NVU_like             = "Non-Tumor: Intact",
  EC_PC_vessel_only              = "Non-Tumor: Astrocyte negative",
  canonical_astro_no_pericyte    = "Non-Tumor: Pericyte negative",
  bare_EC_vessel                 = "Non-Tumor: Endothelial only",
  tumor_associated_intact_NVU    = "Tumor-associated: Intact",
  tumor_vessel_contact_no_astro  = "Tumor-associated: Astrocyte negative",
  tumor_contact_pericyte_dropout = "Tumor-associated: Pericyte negative",
  tumor_contact_bare_vessel      = "Tumor-associated: Endothelial only"
)

nvu_name <- c(canonical_NVU_like = "cNVU", canonical_astro_no_pericyte = "rNVU-1",
  tumor_associated_intact_NVU = "rNVU-4a", tumor_contact_pericyte_dropout = "rNVU-2",
  tumor_vessel_contact_no_astro = "rNVU-3", tumor_contact_bare_vessel = "rNVU-4b",
  EC_PC_vessel_only = "rNVU-5", bare_EC_vessel = "rNVU-6")

# One readable NVU-state field for the RDS and CSV outputs.
unit$nvu_state <- unname(nvu_label[unit$nvu8_id])
if (anyNA(unit$nvu_state)) stop("Missing manuscript label for an NVU state")
unit$long_label <- unit$nvu_state
unit$short_label <- unit$nvu_state
# Legacy IDs are retained for compatibility, not manuscript display.
unit$rnvu_name <- unname(nvu_name[unit$nvu8_id])
unit$nvu4_label <- unname(c(
 canonical_or_astro_NVU="Canonical/astro-supported NVU",
 tumor_intact_NVU="Tumor-contacted intact NVU",
 disrupted_tumor_vessel="Disrupted tumor-contacted vessel",
 vessel_only="Vessel only")[unit$nvu4_id])
composition_map <- c(
  canonical_NVU_like = "EC + PC + Ast, no tumor",
  canonical_astro_no_pericyte = "EC + Ast, no PC, no tumor",
  tumor_associated_intact_NVU = "EC + PC + Ast + Tumor",
  tumor_contact_pericyte_dropout = "EC + Ast + Tumor, no PC",
  tumor_vessel_contact_no_astro = "EC + PC + Tumor, no Ast",
  tumor_contact_bare_vessel = "EC + Tumor, no PC or Ast",
  EC_PC_vessel_only = "EC + PC, no Ast or tumor",
  bare_EC_vessel = "EC only, no PC, Ast, or tumor")
unit$composition <- unname(composition_map[unit$nvu8_id])

# Map both sample and unit ID; 
key <- paste(obj$sample, obj$nvu_unit_id, sep = "|||")
reference_key <- paste(unit$sample, unit$nvu_unit_id, sep = "|||")
for (field in c("nvu8_id", "nvu4_id", "tumor_contact", "tumor_touching",
                "pericyte_present", "astro_present", "nvu_state",
                "short_label", "long_label", "rnvu_name", "composition")) {
  obj@meta.data[[field]] <- unit[[field]][match(key, reference_key)]
}
write.csv(unit, file.path(out, "04_NVU_unit_master.csv"), row.names = FALSE)
write.csv(unit %>% count(sample, nvu8_id, nvu_state, name = "n_vessels"),
          file.path(out, "04_NVU8_counts_by_sample.csv"), row.names = FALSE)



saveRDS(list(seurat_obj = obj, nvu_unit_master = unit, perivascular_radius_qc = radius),
        file.path(out, "02_nvu_states.rds"))
message("NVU state assignment complete")
