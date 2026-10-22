# 03_bbb_analysis.R
library(Seurat)
library(dplyr)
library(tidyr)
library(FNN)

out <- "results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
run <- readRDS(file.path(out, "02_nvu_states.rds"))
obj <- run$seurat_obj
unit <- run$nvu_unit_master
md <- obj@meta.data


needed <- c("sample", "clusters", "vascular_subtype", "nvu_unit_id", "nvu_role",
            "nvu8_id", "nvu4_id", "nvu_state", "BBB_Stability", "BBB_Zone",
            "BBB_Stability_full", "BBB_Stability_indep", "BBB_holdout_UCell",
            "pericyte_present", "astro_present", "vessel_id", "nearest_vessel",
            "fov", "x_um", "y_um")

missing <- setdiff(needed, names(md))
if (length(missing)) stop("Missing metadata: ", paste(missing, collapse = ", "))
if (!"nvu_state" %in% names(unit)) stop("Run the updated 02_nvu_states.R first")
legacy_scored <- is.finite(as.numeric(md$BBB_Stability))
if (!any(legacy_scored & md$vascular_subtype == "EC" & !is.na(md$nvu_unit_id), na.rm = TRUE))
  stop("Historical BBB Stability is missing for vessel-associated ECs")


obj$BBB_Stability_legacy_raw <- as.numeric(md$BBB_Stability)
obj$BBB_Zone_legacy <- as.character(md$BBB_Zone)
obj$BBB_legacy_scoring_universe <- ifelse(legacy_scored, "Legacy EC/PC", NA_character_)


ec <- obj@meta.data %>% tibble::rownames_to_column("cell") %>%
  filter(vascular_subtype == "EC", is.finite(BBB_Stability_indep)) %>%
  group_by(sample) %>%
  mutate(
    BBB_EC_percentile_within_sample = if (n() == 1L) 50 else
      100 * percent_rank(BBB_Stability_indep),
    low_cut = if (n() >= 20) as.numeric(quantile(BBB_Stability_indep, .25)) else NA_real_,
    high_cut = if (n() >= 20) as.numeric(quantile(BBB_Stability_indep, .75)) else NA_real_,
    BBB_EC_status_within_sample = case_when(
      n() < 20 ~ NA_character_,
      BBB_Stability_indep >= high_cut ~ "High",  # High wins if cutoffs coincide
      BBB_Stability_indep <= low_cut ~ "Low",
      TRUE ~ "Mid"
    )
  ) %>% ungroup() %>% select(-low_cut, -high_cut)

obj$BBB_EC_status_within_sample <- NA_character_
obj$BBB_EC_percentile_within_sample <- NA_real_
obj@meta.data[ec$cell, "BBB_EC_status_within_sample"] <- ec$BBB_EC_status_within_sample
obj@meta.data[ec$cell, "BBB_EC_percentile_within_sample"] <- ec$BBB_EC_percentile_within_sample

md <- obj@meta.data
vessels <- md %>%
  filter(vascular_subtype == "EC", !is.na(nvu_unit_id), is.finite(BBB_Stability_indep)) %>%
  group_by(sample, nvu_unit_id) %>%
  summarise(
    n_EC_scored = n(),
    BBB_vessel_score = mean(BBB_Stability_indep),
    BBB_vessel_median = median(BBB_Stability_indep),
    BBB_full_vessel_score = if (all(!is.finite(BBB_Stability_full))) NA_real_ else
      mean(BBB_Stability_full, na.rm = TRUE),
    BBB_holdout_vessel_score = if (all(!is.finite(BBB_holdout_UCell))) NA_real_ else
      mean(BBB_holdout_UCell, na.rm = TRUE),
    .groups = "drop"
  ) %>% group_by(sample) %>%
  mutate(
    BBB_vessel_percentile_within_sample = if (n() == 1L) 50 else
      100 * percent_rank(BBB_vessel_score),
    low_cut = if (n() >= 20) as.numeric(quantile(BBB_vessel_score, .25)) else NA_real_,
    high_cut = if (n() >= 20) as.numeric(quantile(BBB_vessel_score, .75)) else NA_real_,
    BBB_vessel_status_within_sample = case_when(
      n() < 20 ~ NA_character_,
      BBB_vessel_score >= high_cut ~ "High",
      BBB_vessel_score <= low_cut ~ "Low",
      TRUE ~ "Mid"
    )
  ) %>% ungroup() %>% select(-low_cut, -high_cut)


legacy_vessels <- md %>%
  filter(vascular_subtype == "EC", !is.na(nvu_unit_id),
         is.finite(BBB_Stability_legacy_raw)) %>%
  group_by(sample, nvu_unit_id) %>%
  summarise(
    BBB_Stability_legacy_vessel_n_EC = n(),
    BBB_Stability_legacy_vessel_mean = mean(BBB_Stability_legacy_raw),
    BBB_Stability_legacy_vessel_median = median(BBB_Stability_legacy_raw),
    .groups = "drop"
  )

vm <- unit %>% left_join(vessels, by = c("sample", "nvu_unit_id")) %>%
  left_join(legacy_vessels, by = c("sample", "nvu_unit_id"))

# Compare original and independent scores for reviewer validation
score_comparison <- vm %>% group_by(sample) %>% summarise(
  n_original_vessels = sum(is.finite(BBB_Stability_legacy_vessel_mean)),
  n_independent_vessels = sum(is.finite(BBB_vessel_score)),
  rho_original_vs_independent = {
    x <- BBB_Stability_legacy_vessel_mean
    y <- BBB_vessel_score
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) NA_real_
    else cor(x[ok], y[ok], method = "spearman")
  }, .groups = "drop")
write.csv(score_comparison, file.path(out, "05_original_vs_independent_BBB_QC.csv"), row.names = FALSE)


centroids <- obj@meta.data %>%
  filter(vascular_subtype == "EC", !is.na(nvu_unit_id),
         !is.na(fov), is.finite(x_um), is.finite(y_um)) %>%
  group_by(sample, fov, nvu_unit_id) %>%
  summarise(x = mean(x_um), y = mean(y_um), .groups = "drop")
vm <- left_join(vm, centroids, by = c("sample", "nvu_unit_id"))

key <- paste(obj$sample, obj$nvu_unit_id, sep = "|||")
vm_key <- paste(vm$sample, vm$nvu_unit_id, sep = "|||")
for (field in c("BBB_vessel_score", "BBB_vessel_status_within_sample",
                "BBB_vessel_percentile_within_sample",
                "BBB_Stability_legacy_vessel_n_EC",
                "BBB_Stability_legacy_vessel_mean",
                "BBB_Stability_legacy_vessel_median")) {
  obj@meta.data[[field]] <- vm[[field]][match(key, vm_key)]
}
obj$EC_BBB_vessel_status_within_sample <- ifelse(
  obj$vascular_subtype %in% "EC", obj$BBB_vessel_status_within_sample, NA_character_)


validation <- list()
for (nm in c("BBB_full_vessel_score", "BBB_holdout_vessel_score")) {
  keep <- is.finite(vm$BBB_vessel_score) & is.finite(vm[[nm]])
  if (sum(keep) < 5) next
  z <- suppressWarnings(cor.test(vm$BBB_vessel_score[keep], vm[[nm]][keep],
                                 method = "spearman", exact = FALSE))
  validation[[nm]] <- data.frame(comparison = nm, rho = unname(z$estimate), p = z$p.value)
}
write.csv(bind_rows(validation), file.path(out, "05_BBB_validation_correlations.csv"), row.names = FALSE)
cell_bbb <- obj@meta.data %>% tibble::rownames_to_column("cell") %>%
  filter(vascular_subtype == "EC") %>%
  select(cell, sample, nvu_unit_id, nvu8_id, nvu_state,
         BBB_Stability, BBB_Zone, BBB_Stability_legacy_raw, BBB_Zone_legacy,
         BBB_legacy_scoring_universe, BBB_Stability_legacy_vessel_mean,
         BBB_Stability_legacy_vessel_median, BBB_Stability_full,
         BBB_Stability_indep, BBB_EC_percentile_within_sample,
         BBB_EC_status_within_sample, BBB_vessel_score,
         BBB_vessel_percentile_within_sample, BBB_vessel_status_within_sample)
write.csv(cell_bbb, file.path(out, "05_EC_BBB_states.csv"), row.names = FALSE)
write.csv(vm, file.path(out, "05_vessel_BBB_master.csv"), row.names = FALSE)
write.csv(obj@meta.data %>% filter(vascular_subtype == "EC") %>%
            count(sample, BBB_EC_status_within_sample, name = "n_EC"),
          file.path(out, "05_EC_status_counts_by_sample.csv"), row.names = FALSE)
write.csv(vm %>% count(sample, BBB_vessel_status_within_sample, name = "n_vessels"),
          file.path(out, "05_vessel_status_counts_by_sample.csv"), row.names = FALSE)



flip_p <- function(x) {
  x <- x[is.finite(x) & x != 0]
  if (length(x) < 2) return(NA_real_)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(x)), KEEP.OUT.ATTRS = FALSE))
  mean(abs(rowMeans(sweep(signs, 2, x, `*`))) >= abs(mean(x)) - 1e-12)
}

summarise_effects <- function(df, by, effect) {
  if (!nrow(df)) return(tibble::tibble())
  ans <- df %>% group_by(across(all_of(by))) %>%
    summarise(
      n_samples = n_distinct(sample),
      median_effect = median(.data[[effect]], na.rm = TRUE),
      mean_effect = mean(.data[[effect]], na.rm = TRUE),
      fraction_positive = mean(.data[[effect]] > 0, na.rm = TRUE),
      fraction_negative = mean(.data[[effect]] < 0, na.rm = TRUE),
      signflip_p = if (n_distinct(sample) >= 3) flip_p(.data[[effect]]) else NA_real_,
      ci = list({
        x <- .data[[effect]][is.finite(.data[[effect]])]
        if (length(x) < 2) c(NA_real_, NA_real_) else {
          set.seed(20260721)
          boot <- replicate(5000, median(base::sample(x, length(x), replace = TRUE)))
          as.numeric(quantile(boot, c(.025, .975)))
        }
      }),
      .groups = "drop"
    )
  ans$median_ci_low <- vapply(ans$ci, `[`, numeric(1), 1)
  ans$median_ci_high <- vapply(ans$ci, `[`, numeric(1), 2)
  ans$ci <- NULL
  ans$FDR <- p.adjust(ans$signflip_p, "BH")
  ans
}


# Tumor programs relative to the nearest EC
fov_sample <- c(Lung_fov = "Lung", Breast_fov = "Breast", Normal_fov = "Normal",
                Thyroid_fov = "Thyroid", Melanoma_fov = "Melanoma", Kidney_fov = "Kidney",
                fov = "Gr_3", fov.3 = "GBM")
programs <- c("Hypoxia", "Mesenchymal", "Stem_like", "Invasive", "Proliferative")
program_cols <- paste0(programs, "_UCell")
if (!all(program_cols %in% names(obj@meta.data))) stop("Tumor UCell program columns missing")
md <- obj@meta.data
maps <- list()
for (fov in names(fov_sample)) {
  if (!fov %in% names(obj@images)) next
  xy <- as.data.frame(GetTissueCoordinates(obj, image = fov))
  if (!"cell" %in% colnames(xy)) xy$cell <- rownames(xy)
  rownames(xy) <- as.character(xy$cell)
  xc <- intersect(c("x","imagecol","col","center_x","centroid_x"), names(xy))[1]
  yc <- intersect(c("y","imagerow","row","center_y","centroid_y"), names(xy))[1]
  if (is.na(xc) || is.na(yc)) stop("Coordinates missing in ", fov)
  cells <- intersect(rownames(xy), rownames(md))
  cells <- cells[is.finite(xy[cells, xc]) & is.finite(xy[cells, yc])]
  m <- md[cells, , drop = FALSE]
  tumor <- cells[which(as.character(m$clusters) == "Tumor" &
                         m$nvu_role %in% "perivascular-tumor")]
  ec_cells <- cells[which(m$vascular_subtype %in% "EC" &
                            is.finite(m$BBB_Stability_indep))]
  if (!length(tumor) || !length(ec_cells)) next
  k <- min(5L, length(ec_cells))
  nn <- FNN::get.knnx(as.matrix(xy[ec_cells, c(xc, yc), drop = FALSE]),
                      as.matrix(xy[tumor, c(xc, yc), drop = FALSE]), k = k)
  nearest <- ec_cells[nn$nn.index[, 1]]
  scores <- matrix(as.numeric(md[ec_cells, "BBB_Stability_indep"])[nn$nn.index],
                   nrow = length(tumor), ncol = k)
  scores[nn$nn.dist > 30] <- NA_real_
  same_unit <- as.character(md[tumor, "nvu_unit_id"]) ==
    as.character(md[nearest, "nvu_unit_id"])
  same_unit[is.na(md[tumor, "nvu_unit_id"]) |
              is.na(md[nearest, "nvu_unit_id"])] <- NA
  maps[[fov]] <- data.frame(
    cell = tumor, sample = as.character(md[tumor, "sample"]), fov = fov,
    tumor_nearest_vessel = as.character(md[tumor, "nearest_vessel"]),
    tumor_nvu_unit_id = as.character(md[tumor, "nvu_unit_id"]),
    tumor_nvu8_id = as.character(md[tumor, "nvu8_id"]),
    tumor_nvu4_id = as.character(md[tumor, "nvu4_id"]),
    nearest_EC_cell = nearest, nearest_EC_distance_um = nn$nn.dist[, 1],
    nearest_EC_vessel_id = as.character(md[nearest, "vessel_id"]),
    nearest_EC_nvu_unit_id = as.character(md[nearest, "nvu_unit_id"]),
    nearest_EC_nvu8_id = as.character(md[nearest, "nvu8_id"]),
    nearest_EC_nvu4_id = as.character(md[nearest, "nvu4_id"]),
    nearest_EC_nvu_state = as.character(md[nearest, "nvu_state"]),
    nearest_EC_pericyte_present = as.logical(md[nearest, "pericyte_present"]),
    nearest_EC_astro_present = as.logical(md[nearest, "astro_present"]),
    nearest_EC_has_vessel = !is.na(md[nearest, "vessel_id"]),
    same_nvu_unit = same_unit,
    nearest_EC_BBB_score = as.numeric(md[nearest, "BBB_Stability_indep"]),
    nearest_EC_BBB_percentile = as.numeric(md[nearest, "BBB_EC_percentile_within_sample"]),
    nearest_EC_BBB_status = as.character(md[nearest, "BBB_EC_status_within_sample"]),
    mean_kNN_EC_BBB_score = rowMeans(scores, na.rm = TRUE),
    n_EC_within_radius = rowSums(is.finite(scores))
  )
}

mapping <- bind_rows(maps)
if (!nrow(mapping)) stop("No eligible perivascular tumor-to-EC mappings")
mapping$mean_kNN_EC_BBB_score[!is.finite(mapping$mean_kNN_EC_BBB_score)] <- NA_real_
mapping <- mapping %>% filter(is.finite(nearest_EC_distance_um), nearest_EC_distance_um <= 30)
if (!nrow(mapping)) stop("No tumor cells mapped within 30 um of a scored EC")
mapping$nearest_NVU_pericyte_status <- ifelse(is.na(mapping$nearest_EC_pericyte_present),
                                              NA_character_, ifelse(mapping$nearest_EC_pericyte_present, "Retained", "Lost"))
mapping$nearest_NVU_astro_status <- ifelse(is.na(mapping$nearest_EC_astro_present),
                                           NA_character_, ifelse(mapping$nearest_EC_astro_present, "Retained", "Lost"))
mapping <- cbind(mapping, md[mapping$cell, program_cols, drop = FALSE])
write.csv(mapping, file.path(out, "06_tumor_to_local_EC_BBB_mapping.csv"), row.names = FALSE)

mapping_qc <- mapping %>% group_by(sample) %>% summarise(
  n_tumor_mapped = n(),
  n_nearest_EC_with_vessel = sum(nearest_EC_has_vessel),
  n_nearest_EC_without_vessel = sum(!nearest_EC_has_vessel),
  pct_nearest_EC_with_vessel = 100 * mean(nearest_EC_has_vessel),
  n_same_nvu_evaluable = sum(!is.na(same_nvu_unit)),
  n_same_nvu = sum(same_nvu_unit %in% TRUE),
  n_different_nvu = sum(same_nvu_unit %in% FALSE),
  pct_same_nvu_among_evaluable = ifelse(n_same_nvu_evaluable > 0,
                                       100 * n_same_nvu / n_same_nvu_evaluable, NA_real_),
  median_nearest_EC_distance_um = median(nearest_EC_distance_um),
  p90_nearest_EC_distance_um = as.numeric(quantile(nearest_EC_distance_um, .9)),
  .groups = "drop")
write.csv(mapping_qc, file.path(out, "06_tumor_EC_NVU_mapping_QC.csv"), row.names = FALSE)

long <- mapping %>% pivot_longer(all_of(program_cols), names_to = "program_col",
                                values_to = "tumor_program_score") %>%
  mutate(program = sub("_UCell$", "", program_col),
         program_label = recode(program, Stem_like = "Stem-like", .default = program))
write.csv(long, file.path(out, "06_tumor_programs_local_EC_BBB_long.csv"), row.names = FALSE)

high_low <- function(df, groups = c("sample", "program", "program_label")) {
  df %>% filter(nearest_EC_BBB_status %in% c("High", "Low"), is.finite(tumor_program_score)) %>%
    group_by(across(all_of(groups))) %>%
    summarise(
      n_high = sum(nearest_EC_BBB_status == "High"),
      n_low = sum(nearest_EC_BBB_status == "Low"),
      mean_high = mean(tumor_program_score[nearest_EC_BBB_status == "High"]),
      mean_low = mean(tumor_program_score[nearest_EC_BBB_status == "Low"]),
      effect_high_minus_low = mean_high - mean_low,
      .groups = "drop"
    ) %>% filter(n_high >= 10, n_low >= 10, is.finite(effect_high_minus_low))
}
primary <- high_low(long)
primary_summary <- summarise_effects(primary, c("program", "program_label"), "effect_high_minus_low")
same_nvu <- high_low(long %>% filter(same_nvu_unit %in% TRUE))
same_nvu_summary <- summarise_effects(same_nvu, c("program", "program_label"), "effect_high_minus_low")

write.csv(primary, file.path(out, "06_tumor_program_BBB_effects_by_sample.csv"), row.names = FALSE)
write.csv(primary_summary, file.path(out, "06_tumor_program_BBB_effect_summary.csv"), row.names = FALSE)
write.csv(same_nvu, file.path(out, "06_tumor_program_BBB_effects_same_NVU_by_sample.csv"), row.names = FALSE)
write.csv(same_nvu_summary, file.path(out, "06_tumor_program_BBB_effect_summary_same_NVU.csv"), row.names = FALSE)

# Per-sample Spearman correlations: nearest EC BBB score and within-sample EC percentile.
for (measure in c("nearest_EC_BBB_score", "nearest_EC_BBB_percentile")) {
  corr <- long %>% group_by(sample, program, program_label) %>%
    group_modify(~{
      ok <- is.finite(.x[[measure]]) & is.finite(.x$tumor_program_score)
      z <- if (sum(ok) >= 20) try(suppressWarnings(cor.test(
        .x[[measure]][ok], .x$tumor_program_score[ok], method = "spearman", exact = FALSE)), silent = TRUE) else NULL
      data.frame(BBB_predictor = if (measure == "nearest_EC_BBB_score") "BBB_Stability_indep" else "BBB_EC_percentile_within_sample",
                 rho = if (inherits(z, "htest")) unname(z$estimate) else NA_real_,
                 p = if (inherits(z, "htest")) z$p.value else NA_real_, n_cells = sum(ok))
    }) %>% ungroup()
  corr_summary <- summarise_effects(corr %>% rename(effect = rho),
                                    c("program", "program_label", "BBB_predictor"), "effect")
  if (measure == "nearest_EC_BBB_score") {
    write.csv(corr, file.path(out, "06_tumor_program_continuous_correlations_by_sample.csv"), row.names = FALSE)
    write.csv(corr_summary, file.path(out, "06_tumor_program_continuous_correlation_summary.csv"), row.names = FALSE)
  } else {
    write.csv(corr, file.path(out, "06_tumor_program_EC_BBB_percentile_correlations_by_sample.csv"), row.names = FALSE)
    write.csv(corr_summary, file.path(out, "06_tumor_program_EC_BBB_percentile_correlation_summary.csv"), row.names = FALSE)
  }
}


support <- bind_rows(
  long %>% mutate(support_axis = "Pericyte", support_status = nearest_NVU_pericyte_status),
  long %>% mutate(support_axis = "Astrocyte", support_status = nearest_NVU_astro_status)
) %>% filter(same_nvu_unit %in% TRUE, support_status %in% c("Retained", "Lost"),
             is.finite(tumor_program_score))
support_effects <- support %>% group_by(sample, program, program_label, support_axis) %>%
  summarise(n_retained = sum(support_status == "Retained"),
            n_lost = sum(support_status == "Lost"),
            mean_retained = mean(tumor_program_score[support_status == "Retained"]),
            mean_lost = mean(tumor_program_score[support_status == "Lost"]),
            effect_retained_minus_lost = mean_retained - mean_lost, .groups = "drop") %>%
  filter(n_retained >= 10, n_lost >= 10, is.finite(effect_retained_minus_lost))
support_summary <- summarise_effects(support_effects,
                                     c("program", "program_label", "support_axis"), "effect_retained_minus_lost")
stratified <- high_low(support, c("sample", "program", "program_label", "support_axis", "support_status"))
stratified_summary <- summarise_effects(stratified,
                                        c("program", "program_label", "support_axis", "support_status"),
                                        "effect_high_minus_low")
write.csv(support_effects, file.path(out, "06_tumor_program_by_NVU_support_axis_by_sample.csv"), row.names = FALSE)
write.csv(support_summary, file.path(out, "06_tumor_program_by_NVU_support_axis_summary.csv"), row.names = FALSE)
write.csv(stratified, file.path(out, "06_tumor_program_BBB_effects_stratified_by_NVU_support.csv"), row.names = FALSE)
write.csv(stratified_summary, file.path(out, "06_tumor_program_BBB_effect_summary_stratified_by_NVU_support.csv"), row.names = FALSE)



class_summary <- vm %>% group_by(sample, nvu8_id, nvu_state) %>%
  summarise(n_vessels = n(), median_BBB_score = median(BBB_vessel_score, na.rm = TRUE),
            median_BBB_percentile = median(BBB_vessel_percentile_within_sample, na.rm = TRUE),
            .groups = "drop")
status_composition <- vm %>% count(sample, nvu8_id, nvu_state,
                                   BBB_vessel_status_within_sample, name = "n_vessels") %>%
  group_by(sample, nvu8_id, nvu_state) %>%
  mutate(percentage = 100 * n_vessels / sum(n_vessels)) %>% ungroup()

axes <- vm %>% select(sample, nvu_unit_id, BBB_vessel_percentile_within_sample,
                      pericyte_present, astro_present, tumor_contact) %>%
  pivot_longer(c(pericyte_present, astro_present, tumor_contact),
               names_to = "axis", values_to = "present") %>%
  group_by(sample, axis) %>% summarise(
    n_present = sum(present), n_absent = sum(!present),
    mean_present = mean(BBB_vessel_percentile_within_sample[present], na.rm = TRUE),
    mean_absent = mean(BBB_vessel_percentile_within_sample[!present], na.rm = TRUE),
    effect_present_minus_absent = mean_present - mean_absent,
    .groups = "drop") %>% filter(n_present >= 3, n_absent >= 3)
axis_summary <- summarise_effects(axes, "axis", "effect_present_minus_absent")


state_order <- c("canonical_NVU_like", "canonical_astro_no_pericyte",
                 "tumor_associated_intact_NVU", "tumor_contact_pericyte_dropout",
                 "tumor_vessel_contact_no_astro", "tumor_contact_bare_vessel",
                 "EC_PC_vessel_only", "bare_EC_vessel")
state_map <- unit[!duplicated(unit$nvu8_id), c("nvu8_id", "nvu_state")]
state_labels <- setNames(state_map$nvu_state, state_map$nvu8_id)
pair_data <- bind_rows(lapply(combn(state_order, 2, simplify = FALSE), function(p) {
  a <- class_summary %>% filter(nvu8_id == p[1], n_vessels >= 3) %>%
    select(sample, a = median_BBB_percentile, n_a = n_vessels)
  b <- class_summary %>% filter(nvu8_id == p[2], n_vessels >= 3) %>%
    select(sample, b = median_BBB_percentile, n_b = n_vessels)
  inner_join(a, b, by = "sample") %>%
    mutate(class_a = p[1], class_b = p[2], effect_a_minus_b = a - b)
}))
pair_summary <- summarise_effects(pair_data, c("class_a", "class_b"), "effect_a_minus_b")
if (nrow(pair_data)) {
  pair_data$state_a <- unname(state_labels[pair_data$class_a])
  pair_data$state_b <- unname(state_labels[pair_data$class_b])
}
if (nrow(pair_summary)) {
  pair_summary$state_a <- unname(state_labels[pair_summary$class_a])
  pair_summary$state_b <- unname(state_labels[pair_summary$class_b])
}

hl <- vm %>% filter(BBB_vessel_status_within_sample %in% c("High", "Low")) %>%
  mutate(BBB_high = as.integer(BBB_vessel_status_within_sample == "High"),
         sample_factor = factor(sample))
logistic <- try({
  fit <- glm(BBB_high ~ astro_present + pericyte_present + tumor_contact + sample_factor,
             data = hl, family = binomial())
  co <- summary(fit)$coefficients
  data.frame(term = rownames(co), estimate = co[, 1], std_error = co[, 2],
             z = co[, 3], p = co[, 4]) %>% mutate(FDR = p.adjust(p, "BH"))
}, silent = TRUE)
if (inherits(logistic, "try-error")) logistic <- data.frame(term = "MODEL_FAILED")

write.csv(class_summary, file.path(out, "07_NVU8_BBB_summary_by_sample.csv"), row.names = FALSE)
write.csv(status_composition, file.path(out, "07_NVU8_BBB_status_composition.csv"), row.names = FALSE)
write.csv(axes, file.path(out, "07_structural_axis_effects_by_sample.csv"), row.names = FALSE)
write.csv(axis_summary, file.path(out, "07_structural_axis_effect_summary.csv"), row.names = FALSE)
write.csv(pair_data, file.path(out, "07_NVU8_pairwise_effects_by_sample.csv"), row.names = FALSE)
write.csv(pair_summary, file.path(out, "07_NVU8_pairwise_effect_summary.csv"), row.names = FALSE)
write.csv(logistic, file.path(out, "07_BBB_high_adjusted_logistic_model.csv"), row.names = FALSE)


# # 6. Optional numerical check against the original checkpoint
# reference <- "../../NVU_Brain_Metastases/2026_07_21_NVU_analysis_rerun/checkpoints//05_BBB_states.rds"
# if (file.exists(reference)) {
#   old_vm <- readRDS(reference)$vessel_master
#   old_key <- paste(old_vm$sample, old_vm$nvu_unit_id, sep = "|||")
#   if (!setequal(old_key, vm_key)) stop("Vessel IDs differ from the original BBB checkpoint")
#   x <- old_vm$BBB_vessel_score[match(vm_key, old_key)]
#   if (!isTRUE(all.equal(vm$BBB_vessel_score, x, tolerance = 1e-8)))
#     stop("Independent-score vessel means differ from original clean rerun")
#   message("PASS: independent-score vessel means match original clean rerun")
#   if ("BBB_Stability_legacy_vessel_mean" %in% names(old_vm)) {
#     y <- old_vm$BBB_Stability_legacy_vessel_mean[match(vm_key, old_key)]
#     if (!isTRUE(all.equal(vm$BBB_Stability_legacy_vessel_mean, y,
#                           tolerance = 1e-8)))
#       stop("Original raw BBB vessel means differ from reference checkpoint")
#     message("PASS: original raw BBB vessel means match reference checkpoint")
#   }
# }


stopifnot(
  identical(obj$BBB_Stability, run$seurat_obj$BBB_Stability),
  identical(as.character(obj$BBB_Zone), as.character(run$seurat_obj$BBB_Zone))
)

saveRDS(list(seurat_obj = obj, nvu_unit_master = unit, vessel_master = vm,
             mapping = mapping, mapping_qc = mapping_qc,
             tumor_program_results = primary, NVU_BBB_results = axis_summary),
        file.path(out, "03_bbb_analysis.rds"))
message("Step 03 complete: original manuscript BBB retained; clean-rerun independent analyses retained")
