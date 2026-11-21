# ---- Day 12 Avastin Multiplex: IgG QC, BBB classification, and area normalization ----

library(tidyverse)
library(readxl)
library(writexl)

# ---- 1. Define input and output folders ----

day12_input_dir <- "Arpan Multiplex/Day12"
day12_area_norm_dir <- "Arpan Multiplex/Day12_Area_Normalized"

dir.create(
  day12_area_norm_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

day12_excel_files <- list.files(
  path = day12_input_dir,
  pattern = "Vessel_measurements and Perivascular IgG\\.xlsx$",
  full.names = TRUE
)

if (length(day12_excel_files) != 6) {
  stop(
    "Expected 6 Day 12 Excel files, but found ",
    length(day12_excel_files),
    ". Check the input folder and filenames."
  )
}

print(basename(day12_excel_files))

# ---- 2. Read and combine Perivascular IgG tabs ----

igg_raw_Day12 <- map_dfr(day12_excel_files, function(file) {
  
  filename <- basename(file)
  
  treatment <- str_extract(filename, "^(Avastin|Saline)")
  tumor_number <- str_extract(filename, "(?<=Tumor_)\\d+")
  roi_number <- str_extract(filename, "(?<=ROI_)\\d+")
  
  read_excel(
    file,
    sheet = "Perivascular IgG"
  ) %>%
    mutate(
      File = filename,
      Treatment = treatment,
      Tumor = paste0("Tumor-", tumor_number),
      ROI = paste0("ROI-", roi_number),
      Sample_ID = paste(Treatment, Tumor, ROI, sep = "-")
    )
})

# ---- 3. Check QC columns ----

cat("\nUnique vessel centroid status values:\n")
print(unique(igg_raw_Day12$Vessel_centroid_inside_tumor_ROI))

cat("\nUnique PV_QC_status values:\n")
print(unique(igg_raw_Day12$PV_QC_status))

# ---- 4. Apply required IgG QC filters ----

igg_data_Day12 <- igg_raw_Day12 %>%
  filter(
    Vessel_centroid_inside_tumor_ROI == 1,
    PV_QC_status == "Analyzable"
  ) %>%
  select(
    File,
    Sample_ID,
    Treatment,
    Tumor,
    ROI,
    Vessel_ID,
    Vessel_centroid_inside_tumor_ROI,
    PV_QC_status,
    IgG_PV_Mean_bg_corrected,
    Leakage_status
  )

cat("\nNumber of vessels after IgG QC:", nrow(igg_data_Day12), "\n")

cat(
  "Vessels with missing IgG values:",
  sum(is.na(igg_data_Day12$IgG_PV_Mean_bg_corrected)),
  "\n"
)

# ---- 5. Check duplicate vessel identifiers ----

duplicate_igg_ids_Day12 <- igg_data_Day12 %>%
  count(Sample_ID, Vessel_ID) %>%
  filter(n > 1)

if (nrow(duplicate_igg_ids_Day12) > 0) {
  print(duplicate_igg_ids_Day12)
  stop("Duplicate Sample_ID / Vessel_ID combinations in Day 12 IgG data.")
}

# ---- 6. Calculate sample-specific IgG quartiles ----

igg_data_Day12 <- igg_data_Day12 %>%
  group_by(Sample_ID) %>%
  mutate(
    IgG_Q25 = if (all(is.na(IgG_PV_Mean_bg_corrected))) {
      NA_real_
    } else {
      quantile(
        IgG_PV_Mean_bg_corrected,
        probs = 0.25,
        na.rm = TRUE,
        names = FALSE
      )
    },
    IgG_Q75 = if (all(is.na(IgG_PV_Mean_bg_corrected))) {
      NA_real_
    } else {
      quantile(
        IgG_PV_Mean_bg_corrected,
        probs = 0.75,
        na.rm = TRUE,
        names = FALSE
      )
    },
    IgG_qualifier = case_when(
      is.na(IgG_PV_Mean_bg_corrected) ~ NA_character_,
      IgG_PV_Mean_bg_corrected <= IgG_Q25 ~ "Low",
      IgG_PV_Mean_bg_corrected >= IgG_Q75 ~ "High",
      TRUE ~ "Mid"
    )
  ) %>%
  ungroup()

# ---- 7. Assign BBB Stability and IgG sign ----

igg_data_Day12 <- igg_data_Day12 %>%
  mutate(
    `BBB Stability` = case_when(
      IgG_qualifier == "Low" ~ "High BBB",
      IgG_qualifier == "Mid" ~ "Mid BBB",
      IgG_qualifier == "High" ~ "Low BBB",
      TRUE ~ NA_character_
    ),
    IgG_sign = case_when(
      IgG_PV_Mean_bg_corrected < 0 ~ "Negative",
      IgG_PV_Mean_bg_corrected > 0 ~ "Positive",
      IgG_PV_Mean_bg_corrected == 0 ~ "Zero",
      TRUE ~ NA_character_
    ),
    `BBB Stability` = factor(
      `BBB Stability`,
      levels = c("High BBB", "Mid BBB", "Low BBB")
    ),
    IgG_qualifier = factor(
      IgG_qualifier,
      levels = c("Low", "Mid", "High")
    )
  )

cat("\nBBB Stability classification:\n")
print(table(igg_data_Day12$`BBB Stability`, useNA = "ifany"))

cat("\nIgG qualifier classification:\n")
print(table(igg_data_Day12$IgG_qualifier, useNA = "ifany"))

# ---- 8. Summarize IgG classifications ----

cat("\nIgG classification by sample:\n")

print(
  igg_data_Day12 %>%
    count(Sample_ID, Treatment, `BBB Stability`)
)

cat("\nIgG sign by qualifier:\n")

print(
  table(
    igg_data_Day12$IgG_qualifier,
    igg_data_Day12$IgG_sign,
    useNA = "ifany"
  )
)

# ---- 9. Read first-tab vessel measurements ----

vessel_gene_data_Day12 <- map_dfr(day12_excel_files, function(file) {
  
  filename <- basename(file)
  
  treatment <- str_extract(filename, "^(Avastin|Saline)")
  tumor_number <- str_extract(filename, "(?<=Tumor_)\\d+")
  roi_number <- str_extract(filename, "(?<=ROI_)\\d+")
  
  read_excel(file, sheet = 1) %>%
    mutate(
      Treatment = treatment,
      Tumor = paste0("Tumor-", tumor_number),
      ROI = paste0("ROI-", roi_number),
      Sample_ID = paste(Treatment, Tumor, ROI, sep = "-")
    )
})

# ---- 10. Check vessel IDs ----

duplicate_gene_ids_Day12 <- vessel_gene_data_Day12 %>%
  count(Sample_ID, Vessel_ID) %>%
  filter(n > 1)

if (nrow(duplicate_gene_ids_Day12) > 0) {
  print(duplicate_gene_ids_Day12)
  stop("Duplicate Sample_ID / Vessel_ID combinations in Day 12 vessel data.")
}

# ---- 11. Check overlap before joining ----

gene_ids_Day12 <- vessel_gene_data_Day12 %>%
  distinct(Sample_ID, Vessel_ID)

igg_ids_Day12 <- igg_data_Day12 %>%
  distinct(Sample_ID, Vessel_ID)

cat("\nFirst-tab vessels without a QC-passing IgG match:\n")

print(
  anti_join(
    gene_ids_Day12,
    igg_ids_Day12,
    by = c("Sample_ID", "Vessel_ID")
  )
)

cat("\nQC-passing IgG vessels without a first-tab match:\n")

print(
  anti_join(
    igg_ids_Day12,
    gene_ids_Day12,
    by = c("Sample_ID", "Vessel_ID")
  )
)

# ---- 12. Join IgG classifications to vessel measurements ----
# Retain vessels that pass IgG QC and match the first tab.

vessel_analysis_full_data_Day12 <- vessel_gene_data_Day12 %>%
  inner_join(
    igg_data_Day12 %>%
      select(
        Sample_ID,
        Vessel_ID,
        IgG_PV_Mean_bg_corrected,
        IgG_sign,
        IgG_qualifier,
        IgG_Q25,
        IgG_Q75,
        `BBB Stability`
      ),
    by = c("Sample_ID", "Vessel_ID")
  )

# ---- 13. Define marker columns and labels ----

day12_marker_cols <- c(
  "CD31_Mean_bg_corrected",
  "GLUT1_Mean_bg_corrected",
  "PLVAP_Mean_bg_corrected",
  "MFSD2A_Mean_bg_corrected",
  "ABCB1_Mean_bg_corrected"
)

day12_marker_labels <- c(
  CD31_Mean_bg_corrected = "CD31",
  GLUT1_Mean_bg_corrected = "GLUT1",
  PLVAP_Mean_bg_corrected = "PLVAP",
  MFSD2A_Mean_bg_corrected = "MFSD2A",
  ABCB1_Mean_bg_corrected = "ABCB1"
)

day12_marker_order <- c(
  "CD31", "GLUT1", "PLVAP", "MFSD2A", "ABCB1"
)

day12_treatment_levels <- c("Saline", "Avastin")
day12_bbb_levels <- c("High", "Mid", "Low")

# ---- 14. Select analysis variables ----

vessel_analysis_Day12 <- vessel_analysis_full_data_Day12 %>%
  select(
    Sample_ID,
    Treatment,
    Tumor,
    ROI,
    Vessel_ID,
    IgG_PV_Mean_bg_corrected,
    IgG_sign,
    IgG_qualifier,
    IgG_Q25,
    IgG_Q75,
    `BBB Stability`,
    all_of(day12_marker_cols),
    Area_calibrated
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = day12_treatment_levels
    ),
    Area_calibrated = suppressWarnings(
      as.numeric(Area_calibrated)
    ),
    across(
      all_of(day12_marker_cols),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

# ---- 15. Calculate area-normalized scores in wide format ----

day12_area_norm_map <- c(
  CD31_Mean_bg_corrected = "CD31_per_area",
  GLUT1_Mean_bg_corrected = "GLUT1_per_area",
  PLVAP_Mean_bg_corrected = "PLVAP_per_area",
  MFSD2A_Mean_bg_corrected = "MFSD2A_per_area",
  ABCB1_Mean_bg_corrected = "ABCB1_per_area"
)

for (marker_col in names(day12_area_norm_map)) {
  
  normalized_col <- unname(day12_area_norm_map[marker_col])
  
  vessel_analysis_Day12[[normalized_col]] <- ifelse(
    is.finite(vessel_analysis_Day12[[marker_col]]) &
      is.finite(vessel_analysis_Day12$Area_calibrated) &
      vessel_analysis_Day12$Area_calibrated > 0,
    vessel_analysis_Day12[[marker_col]] /
      vessel_analysis_Day12$Area_calibrated,
    NA_real_
  )
}

# Log2(x + 1), matching the Day 5 long-format workflow
day12_area_norm_cols <- unname(day12_area_norm_map)

for (normalized_col in day12_area_norm_cols) {
  
  log_col <- paste0(normalized_col, "_log2")
  
  vessel_analysis_Day12[[log_col]] <- ifelse(
    is.finite(vessel_analysis_Day12[[normalized_col]]) &
      vessel_analysis_Day12[[normalized_col]] >= 0,
    log2(vessel_analysis_Day12[[normalized_col]] + 1),
    NA_real_
  )
}

# ---- 16. Check calibrated area and normalized values ----

cat("\nArea QC summary:\n")

print(
  vessel_analysis_Day12 %>%
    summarise(
      n_vessels = n(),
      n_zero_area = sum(Area_calibrated == 0, na.rm = TRUE),
      n_negative_area = sum(Area_calibrated < 0, na.rm = TRUE),
      n_missing_area = sum(is.na(Area_calibrated))
    )
)

cat("\nArea-normalized marker QC:\n")

print(
  vessel_analysis_Day12 %>%
    summarise(
      across(
        all_of(day12_area_norm_cols),
        list(
          n_missing = ~sum(is.na(.x)),
          n_nonpositive = ~sum(.x <= 0, na.rm = TRUE)
        )
      )
    ),
  width = Inf
)

# ---- 17. Create long-format marker data ----

day12_all_vessels_long <- vessel_analysis_Day12 %>%
  pivot_longer(
    cols = all_of(day12_marker_cols),
    names_to = "Marker_column",
    values_to = "Expression"
  ) %>%
  mutate(
    Marker = unname(day12_marker_labels[Marker_column]),
    
    Expression = suppressWarnings(as.numeric(Expression)),
    
    Expression_per_area = if_else(
      is.finite(Expression) &
        is.finite(Area_calibrated) &
        Area_calibrated > 0,
      Expression / Area_calibrated,
      NA_real_
    ),
    
    # Separate derived measure; not the area-normalized score
    Estimated_integrated_signal = if_else(
      is.finite(Expression) &
        is.finite(Area_calibrated) &
        Area_calibrated > 0,
      Expression * Area_calibrated,
      NA_real_
    ),
    
    Expression_log2 = if_else(
      is.finite(Expression) & Expression >= 0,
      log2(Expression + 1),
      NA_real_
    ),
    
    Expression_per_area_log2 = if_else(
      is.finite(Expression_per_area) &
        Expression_per_area >= 0,
      log2(Expression_per_area + 1),
      NA_real_
    ),
    
    Estimated_integrated_signal_log2 = if_else(
      is.finite(Estimated_integrated_signal) &
        Estimated_integrated_signal >= 0,
      log2(Estimated_integrated_signal + 1),
      NA_real_
    ),
    
    Marker = factor(
      Marker,
      levels = day12_marker_order
    ),
    
    Treatment = factor(
      Treatment,
      levels = day12_treatment_levels
    ),
    
    BBB_State = factor(
      recode(
        as.character(`BBB Stability`),
        "High BBB" = "High",
        "Mid BBB" = "Mid",
        "Low BBB" = "Low"
      ),
      levels = day12_bbb_levels
    ),
    
    BBB_Stability_original = `BBB Stability`,
    IgG_quartile = IgG_qualifier
  ) %>%
  select(
    Sample_ID,
    Treatment,
    Tumor,
    ROI,
    Vessel_ID,
    Area_calibrated,
    IgG_PV_Mean_bg_corrected,
    IgG_sign,
    IgG_qualifier,
    IgG_quartile,
    IgG_Q25,
    IgG_Q75,
    BBB_Stability_original,
    BBB_State,
    Marker,
    Marker_column,
    Expression,
    Expression_per_area,
    Estimated_integrated_signal,
    Expression_log2,
    Expression_per_area_log2,
    Estimated_integrated_signal_log2
  )

# ---- 18. Validate long-format data ----

cat("\nLong-format dimensions:\n")
print(dim(day12_all_vessels_long))

cat("\nUnique vessels in wide-format dataset:\n")
print(
  vessel_analysis_Day12 %>%
    summarise(
      n_rows = n(),
      n_unique_vessels = n_distinct(
        paste(Sample_ID, Vessel_ID, sep = "__")
      ),
      n_samples = n_distinct(Sample_ID)
    )
)

cat("\nVessel counts by sample and BBB state:\n")

print(
  vessel_analysis_Day12 %>%
    count(Sample_ID, Treatment, `BBB Stability`)
)

cat("\nUsable area-normalized observations by marker and treatment:\n")

print(
  day12_all_vessels_long %>%
    group_by(Marker, Treatment) %>%
    summarise(
      n_vessels = sum(!is.na(Expression_per_area_log2)),
      n_samples = n_distinct(
        Sample_ID[!is.na(Expression_per_area_log2)]
      ),
      .groups = "drop"
    ),
  n = Inf
)

cat("\nArea-normalized values by BBB state:\n")

print(
  day12_all_vessels_long %>%
    group_by(Marker, BBB_State) %>%
    summarise(
      n_vessels = sum(!is.na(Expression_per_area_log2)),
      n_missing = sum(is.na(Expression_per_area_log2)),
      .groups = "drop"
    ),
  n = Inf
)

# ---- 19. Export Day 12 datasets ----

write_csv(
  igg_data_Day12,
  file.path(day12_area_norm_dir, "Day12_IgG_QC_and_BBB_Classification.csv")
)

write_csv(
  vessel_analysis_full_data_Day12,
  file.path(day12_area_norm_dir, "Day12_Full_QC_Matched_Vessel_Data.csv")
)

write_csv(
  vessel_analysis_Day12,
  file.path(day12_area_norm_dir, "Day12_Vessel_Analysis_Area_Normalized_Wide.csv")
)

write_csv(
  day12_all_vessels_long,
  file.path(day12_area_norm_dir, "Day12_All_QC_Passed_Vessels_Long_Format.csv")
)

write_xlsx(
  list(
    IgG_Classification = igg_data_Day12,
    Full_Matched_Data = vessel_analysis_full_data_Day12,
    Area_Normalized_Wide = vessel_analysis_Day12
  ),
  path = file.path(day12_area_norm_dir, "Day12_Area_Normalized_Data.xlsx")
)

cat("\nDay 12 datasets exported to:\n")
cat(normalizePath(day12_area_norm_dir), "\n")



#----

# ---- Day 12: Complete area-normalized Avastin analysis ----
# Uses day12_all_vessels_long.
# No print() commands.
# Biological replicate = mouse/sample, not individual vessel.

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(lme4)
library(lmerTest)
library(emmeans)
library(purrr)
library(tibble)

# ---- 1. Output folders and settings ----

day12_root_dir <- file.path("Arpan Multiplex", "Day12_Area_Normalized")
day12_results_dir <- file.path(day12_root_dir, "Treatment_Models")
day12_sensitivity_dir <- file.path(day12_root_dir, "Sample_Level_Sensitivity")
day12_area_dir <- file.path(day12_root_dir, "Vessel_Area")
day12_igg_dir <- file.path(day12_root_dir, "IgG_Quantile_Distributions")
day12_igg_spread_dir <- file.path(day12_root_dir, "IgG_Variability")
day12_plvap_dir <- file.path(day12_root_dir, "PLVAP_Spread")
day12_mfsd_plvap_dir <- file.path(day12_root_dir, "Correlations_MFSD2A_PLVAP")
day12_glut_mfsd_dir <- file.path(day12_root_dir, "Correlations_GLUT1_MFSD2A")

day12_dirs <- c(
  day12_results_dir, day12_sensitivity_dir, day12_area_dir,
  day12_igg_dir, day12_igg_spread_dir, day12_plvap_dir,
  day12_mfsd_plvap_dir, day12_glut_mfsd_dir
)

invisible(lapply(
  day12_dirs, dir.create,
  recursive = TRUE, showWarnings = FALSE
))

stopifnot(exists("day12_all_vessels_long"))

day12_treatment_colors <- c(
  "Saline" = "#0072B2",
  "Avastin" = "#D55E00"
)

day12_marker_order <- c("GLUT1", "PLVAP", "MFSD2A", "ABCB1", "CD31")
day12_main_markers <- c("GLUT1", "PLVAP", "MFSD2A", "ABCB1")
day12_bbb_order <- c("High BBB", "Mid BBB", "Low BBB")

day12_bbb_colors <- c(
  "High BBB" = "#009E73",
  "Mid BBB" = "#E6AB02",
  "Low BBB" = "#D55E00"
)

day12_quantile_colors <- c(
  Q25 = "#0072B2",
  Median = "black",
  Q75 = "#D55E00"
)

# All plots use Helvetica; font sizes are approximately 50% larger.
day12_base_theme <- theme_classic(base_size = 22.5, base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica"),
    plot.title = element_text(face = "bold", size = 25.5),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(face = "bold", size = 21),
    axis.text = element_text(size = 18, color = "black"),
    strip.text = element_text(face = "bold", size = 19.5),
    strip.background = element_rect(fill = "grey92", color = "grey75"),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 16.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )

# Consistent PNG/SVG export.
day12_save_plot <- function(plot, folder, filename, width, height) {
  ggsave(
    file.path(folder, paste0(filename, ".png")),
    plot, width = width, height = height, dpi = 400,
    bg = "white"
  )
  ggsave(
    file.path(folder, paste0(filename, ".svg")),
    plot, width = width, height = height,
    bg = "white"
  )
}

day12_df_method <- if (requireNamespace("pbkrtest", quietly = TRUE)) {
  "kenward-roger"
} else {
  "satterthwaite"
}

message("Day 12 mixed-model degrees-of-freedom method: ", day12_df_method)

# ---- 2. Standardize and validate the Day 12 data ----

day12_raw <- day12_all_vessels_long

day12_required_columns <- c(
  "Sample_ID", "Treatment", "Vessel_ID", "Marker",
  "Expression_per_area_log2", "Expression_log2",
  "Area_calibrated", "IgG_PV_Mean_bg_corrected"
)

day12_missing_columns <- setdiff(day12_required_columns, names(day12_raw))

if (length(day12_missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(day12_missing_columns, collapse = ", ")
  )
}

if (!"BBB_Stability_original" %in% names(day12_raw) &&
    !"BBB_State" %in% names(day12_raw)) {
  stop("Neither BBB_Stability_original nor BBB_State exists in day12_raw.")
}

day12_data <- day12_raw %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    Vessel_ID = as.character(Vessel_ID),
    Marker = toupper(trimws(as.character(Marker))),
    Treatment_original = as.character(Treatment),
    
    Treatment = case_when(
      str_detect(str_to_lower(Treatment_original), "saline|pbs|vehicle") ~ "Saline",
      str_detect(str_to_lower(Treatment_original), "avastin") ~ "Avastin",
      TRUE ~ NA_character_
    ),
    
    BBB_raw = if ("BBB_Stability_original" %in% names(day12_raw)) {
      as.character(BBB_Stability_original)
    } else {
      as.character(BBB_State)
    },
    
    BBB_State = case_when(
      str_detect(str_to_lower(BBB_raw), "high") ~ "High BBB",
      str_detect(str_to_lower(BBB_raw), "mid") ~ "Mid BBB",
      str_detect(str_to_lower(BBB_raw), "low") ~ "Low BBB",
      TRUE ~ NA_character_
    ),
    
    Expression_per_area_log2 = as.numeric(Expression_per_area_log2),
    Expression_log2 = as.numeric(Expression_log2),
    Area_calibrated = as.numeric(Area_calibrated),
    IgG_value = as.numeric(IgG_PV_Mean_bg_corrected)
  ) %>%
  filter(
    !is.na(Sample_ID),
    !is.na(Vessel_ID),
    Treatment %in% c("Saline", "Avastin")
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Saline", "Avastin")),
    BBB_State = factor(BBB_State, levels = day12_bbb_order),
    Marker = factor(Marker, levels = day12_marker_order)
  )

day12_sample_counts <- day12_data %>%
  distinct(Sample_ID, Treatment) %>%
  count(Treatment, name = "n_samples")

write_csv(
  day12_sample_counts,
  file.path(day12_root_dir, "Day12_Sample_Counts.csv")
)

day12_vessel_counts <- day12_data %>%
  distinct(Sample_ID, Treatment, Vessel_ID, BBB_State) %>%
  count(Treatment, BBB_State, name = "n_vessels")

write_csv(
  day12_vessel_counts,
  file.path(day12_root_dir, "Day12_Vessel_Counts_By_Treatment_BBB.csv")
)

if (n_distinct(day12_data$Treatment) < 2) {
  stop("Both Saline and Avastin samples must be present.")
}

# ---- 3. Main treatment mixed models by BBB state ----

day12_model_data <- day12_data %>%
  filter(
    as.character(Marker) %in% day12_main_markers,
    !is.na(BBB_State),
    is.finite(Expression_per_area_log2)
  ) %>%
  mutate(Marker = factor(as.character(Marker), levels = day12_main_markers))

day12_fit_one_model <- function(marker_name, bbb_name) {
  
  dat <- day12_model_data %>%
    filter(
      as.character(Marker) == marker_name,
      as.character(BBB_State) == bbb_name
    ) %>%
    droplevels()
  
  base_stats <- tibble(
    Marker = marker_name,
    BBB_State = bbb_name,
    n_vessels = nrow(dat),
    n_samples = n_distinct(dat$Sample_ID),
    log2_difference = NA_real_,
    SE = NA_real_,
    df = NA_real_,
    p_value = NA_real_,
    fold_change = NA_real_,
    FC_CI_lower = NA_real_,
    FC_CI_upper = NA_real_,
    model_status = "Insufficient data"
  )
  
  if (nrow(dat) < 4 || n_distinct(dat$Treatment) < 2 ||
      n_distinct(dat$Sample_ID) < 3) {
    return(list(model = NULL, stats = base_stats))
  }
  
  model <- tryCatch(
    lmerTest::lmer(
      Expression_per_area_log2 ~ Treatment + (1 | Sample_ID),
      data = dat, REML = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(model)) {
    base_stats$model_status <- "Model failed"
    return(list(model = NULL, stats = base_stats))
  }
  
  cr <- tryCatch({
    em <- emmeans::emmeans(model, ~ Treatment, lmer.df = day12_df_method)
    ct <- emmeans::contrast(
      em,
      method = list("Avastin - Saline" = c(-1, 1)),
      adjust = "none"
    )
    as.data.frame(summary(ct, infer = c(TRUE, TRUE), level = 0.95))
  }, error = function(e) NULL)
  
  if (is.null(cr) || nrow(cr) == 0) {
    base_stats$model_status <- "Contrast failed"
    return(list(model = model, stats = base_stats))
  }
  
  lo <- intersect(c("lower.CL", "asymp.LCL"), names(cr))[1]
  hi <- intersect(c("upper.CL", "asymp.UCL"), names(cr))[1]
  
  stats <- tibble(
    Marker = marker_name,
    BBB_State = bbb_name,
    n_vessels = nrow(dat),
    n_samples = n_distinct(dat$Sample_ID),
    log2_difference = cr$estimate[1],
    SE = cr$SE[1],
    df = if ("df" %in% names(cr)) cr$df[1] else NA_real_,
    p_value = cr$p.value[1],
    fold_change = 2^cr$estimate[1],
    FC_CI_lower = if (!is.na(lo)) 2^cr[[lo]][1] else NA_real_,
    FC_CI_upper = if (!is.na(hi)) 2^cr[[hi]][1] else NA_real_,
    model_status = "OK"
  )
  
  list(model = model, stats = stats)
}

day12_model_grid <- expand_grid(
  Marker = day12_main_markers,
  BBB_State = day12_bbb_order
)

day12_model_results <- purrr::pmap(
  day12_model_grid,
  function(Marker, BBB_State) day12_fit_one_model(Marker, BBB_State)
)

day12_model_stats <- bind_rows(lapply(day12_model_results, `[[`, "stats")) %>%
  group_by(BBB_State) %>%
  mutate(
    p_value_Holm = p.adjust(p_value, method = "holm"),
    significant_Holm = !is.na(p_value_Holm) & p_value_Holm < 0.05,
    significance = case_when(
      is.na(p_value_Holm) ~ NA_character_,
      p_value_Holm < 0.001 ~ "***",
      p_value_Holm < 0.01 ~ "**",
      p_value_Holm < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup() %>%
  mutate(
    Marker = factor(Marker, levels = day12_main_markers),
    BBB_State = factor(BBB_State, levels = day12_bbb_order)
  ) %>%
  arrange(BBB_State, Marker)

write_csv(
  day12_model_stats,
  file.path(day12_results_dir, "Day12_Area_Normalized_Treatment_Models_By_BBB_State.csv")
)

# Combined treatment plot: all markers and BBB states.

p_day12_treatment <- ggplot(
  day12_model_data,
  aes(Treatment, Expression_per_area_log2, color = Treatment, fill = Treatment)
) +
  geom_violin(alpha = 0.14, trim = TRUE, linewidth = 0.6) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.25, linewidth = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.28, size = 1.2) +
  facet_grid(Marker ~ BBB_State, scales = "free_y") +
  scale_color_manual(values = day12_treatment_colors) +
  scale_fill_manual(values = day12_treatment_colors) +
  labs(
    title = "Day 12 area-normalized vascular marker expression",
    subtitle = "Vessel-level distributions; mixed models account for sample clustering",
    x = NULL, y = "Area-normalized expression [log2(x + 1)]"
  ) +
  day12_base_theme +
  theme(legend.position = "top")

day12_save_plot(
  p_day12_treatment, day12_results_dir,
  "Day12_Area_Normalized_Treatment_By_BBB", 18, 14
)

# Separate 2 x 2 gene figure for each BBB state.

for (bbb in day12_bbb_order) {
  
  plot_dat <- day12_model_data %>%
    filter(as.character(BBB_State) == bbb)
  
  p <- ggplot(
    plot_dat,
    aes(Treatment, Expression_per_area_log2, color = Treatment, fill = Treatment)
  ) +
    geom_violin(alpha = 0.14, trim = TRUE, linewidth = 0.6) +
    geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.25, linewidth = 0.7) +
    geom_jitter(width = 0.12, alpha = 0.28, size = 1.2) +
    facet_wrap(~ Marker, scales = "free_y", ncol = 2) +
    scale_color_manual(values = day12_treatment_colors) +
    scale_fill_manual(values = day12_treatment_colors) +
    labs(
      title = paste("Day 12 treatment comparison:", bbb),
      subtitle = "Four area-normalized markers; vessels are clustered within samples",
      x = NULL, y = "Area-normalized expression [log2(x + 1)]"
    ) +
    day12_base_theme +
    theme(legend.position = "top")
  
  day12_save_plot(
    p, day12_results_dir,
    paste0("Day12_Treatment_2x2_", gsub(" ", "_", bbb)),
    13, 12
  )
}

# ---- 4. Sample-level sensitivity analysis ----
# One mean per mouse, marker, and BBB state; inference is at sample level.

day12_sensitivity_input <- day12_data %>%
  filter(
    as.character(Marker) %in% day12_main_markers,
    is.finite(Expression_per_area_log2)
  ) %>%
  mutate(Marker = as.character(Marker))

day12_summarize_samples <- function(dat, analysis_label) {
  dat %>%
    group_by(Marker, Treatment, Sample_ID) %>%
    summarise(
      Sample_mean_log2 = mean(Expression_per_area_log2, na.rm = TRUE),
      Sample_median_log2 = median(Expression_per_area_log2, na.rm = TRUE),
      n_vessels = sum(is.finite(Expression_per_area_log2)),
      .groups = "drop"
    ) %>%
    mutate(Analysis = analysis_label)
}

day12_sample_all <- day12_summarize_samples(
  day12_sensitivity_input, "All vessels"
)

day12_sample_by_bbb <- day12_sensitivity_input %>%
  filter(!is.na(BBB_State)) %>%
  group_by(BBB_State) %>%
  group_split() %>%
  lapply(function(dat) {
    day12_summarize_samples(dat, as.character(unique(dat$BBB_State))[1])
  }) %>%
  bind_rows()

day12_sample_summaries <- bind_rows(day12_sample_all, day12_sample_by_bbb)

write_csv(
  day12_sample_summaries,
  file.path(day12_sensitivity_dir, "Day12_Sample_Level_Means_All_and_By_BBB.csv")
)

day12_welch_test <- function(sample_data, analysis_label) {
  
  bind_rows(lapply(day12_main_markers, function(marker_name) {
    
    dat <- sample_data %>% filter(Marker == marker_name)
    saline <- dat %>% filter(Treatment == "Saline") %>% pull(Sample_mean_log2)
    avastin <- dat %>% filter(Treatment == "Avastin") %>% pull(Sample_mean_log2)
    
    saline <- saline[is.finite(saline)]
    avastin <- avastin[is.finite(avastin)]
    
    if (length(saline) < 2 || length(avastin) < 2) {
      return(tibble(
        Analysis = analysis_label, Marker = marker_name,
        n_Saline = length(saline), n_Avastin = length(avastin),
        Mean_Saline = if (length(saline)) mean(saline) else NA_real_,
        SD_Saline = if (length(saline) > 1) sd(saline) else NA_real_,
        Mean_Avastin = if (length(avastin)) mean(avastin) else NA_real_,
        SD_Avastin = if (length(avastin) > 1) sd(avastin) else NA_real_,
        Difference_Avastin_minus_Saline = NA_real_,
        CI_lower = NA_real_, CI_upper = NA_real_, p_value = NA_real_
      ))
    }
    
    test <- tryCatch(t.test(avastin, saline, var.equal = FALSE),
                     error = function(e) NULL)
    
    if (is.null(test)) {
      return(tibble(
        Analysis = analysis_label, Marker = marker_name,
        n_Saline = length(saline), n_Avastin = length(avastin),
        Mean_Saline = mean(saline), SD_Saline = sd(saline),
        Mean_Avastin = mean(avastin), SD_Avastin = sd(avastin),
        Difference_Avastin_minus_Saline = mean(avastin) - mean(saline),
        CI_lower = NA_real_, CI_upper = NA_real_, p_value = NA_real_
      ))
    }
    
    tibble(
      Analysis = analysis_label, Marker = marker_name,
      n_Saline = length(saline), n_Avastin = length(avastin),
      Mean_Saline = mean(saline), SD_Saline = sd(saline),
      Mean_Avastin = mean(avastin), SD_Avastin = sd(avastin),
      Difference_Avastin_minus_Saline = mean(avastin) - mean(saline),
      CI_lower = unname(test$conf.int[1]),
      CI_upper = unname(test$conf.int[2]),
      p_value = test$p.value
    )
  })) %>%
    mutate(
      p_value_Holm = p.adjust(p_value, method = "holm"),
      significance_Holm = case_when(
        is.na(p_value_Holm) ~ NA_character_,
        p_value_Holm < 0.001 ~ "***",
        p_value_Holm < 0.01 ~ "**",
        p_value_Holm < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      Analysis = analysis_label
    )
}

day12_welch_all <- day12_welch_test(day12_sample_all, "All vessels")

day12_welch_bbb <- day12_sample_by_bbb %>%
  group_by(Analysis) %>%
  group_split() %>%
  lapply(function(dat) {
    day12_welch_test(dat, as.character(unique(dat$Analysis))[1])
  }) %>%
  bind_rows()

day12_welch_results <- bind_rows(day12_welch_all, day12_welch_bbb)

write_csv(
  day12_welch_results,
  file.path(day12_sensitivity_dir, "Day12_Sample_Level_Welch_Statistics.csv")
)

# Combined sample-level plot.

p_day12_sample_means <- ggplot(
  day12_sample_summaries,
  aes(Treatment, Sample_mean_log2, color = Treatment)
) +
  geom_line(aes(group = Sample_ID), color = "grey65", linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 4) +
  facet_grid(Marker ~ Analysis, scales = "free_y") +
  scale_color_manual(values = day12_treatment_colors) +
  labs(
    title = "Day 12 sample-level sensitivity analysis",
    subtitle = "Each point represents one mouse/sample",
    x = NULL, y = "Mean area-normalized expression [log2(x + 1)]"
  ) +
  day12_base_theme +
  theme(legend.position = "top")

day12_save_plot(
  p_day12_sample_means, day12_sensitivity_dir,
  "Day12_Sample_Level_Means", 20, 14
)

# Boxplots using sample-level means, one 2 x 2 figure per BBB state.
# Points are the biological replicates; the boxplot is descriptive with n
# determined by the number of available mouse-level summaries.

for (bbb in day12_bbb_order) {
  
  plot_dat <- day12_sample_by_bbb %>%
    filter(Analysis == bbb) %>%
    mutate(Marker = factor(Marker, levels = day12_main_markers))
  
  p <- ggplot(
    plot_dat,
    aes(Treatment, Sample_mean_log2, color = Treatment, fill = Treatment)
  ) +
    geom_boxplot(
      width = 0.55, alpha = 0.20, outlier.shape = NA,
      linewidth = 0.8
    ) +
    geom_jitter(
      width = 0.08, size = 4, alpha = 0.9
    ) +
    facet_wrap(~ Marker, scales = "free_y", ncol = 2) +
    scale_color_manual(values = day12_treatment_colors) +
    scale_fill_manual(values = day12_treatment_colors) +
    labs(
      title = paste("Day 12 sample-level sensitivity:", bbb),
      subtitle = "Boxplots summarize sample means; each point is one mouse",
      x = NULL, y = "Sample mean expression [log2(x + 1)]"
    ) +
    day12_base_theme +
    theme(legend.position = "top")
  
  day12_save_plot(
    p, day12_sensitivity_dir,
    paste0("Day12_Sensitivity_Boxplots_2x2_", gsub(" ", "_", bbb)),
    13, 12
  )
}

# ---- 5. Calibrated vessel-area analysis ----

day12_area_data <- day12_data %>%
  filter(is.finite(Area_calibrated), Area_calibrated >= 0) %>%
  distinct(Sample_ID, Treatment, Vessel_ID, BBB_State, Area_calibrated) %>%
  mutate(Area_log2 = log2(Area_calibrated + 1))

write_csv(
  day12_area_data,
  file.path(day12_area_dir, "Day12_Calibrated_Vessel_Area_Data.csv")
)

day12_area_sample_summary <- day12_area_data %>%
  group_by(Sample_ID, Treatment, BBB_State) %>%
  summarise(
    n_vessels = n(),
    mean_area = mean(Area_calibrated, na.rm = TRUE),
    median_area = median(Area_calibrated, na.rm = TRUE),
    mean_log2_area = mean(Area_log2, na.rm = TRUE),
    median_log2_area = median(Area_log2, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  day12_area_sample_summary,
  file.path(day12_area_dir, "Day12_Area_Sample_Level_Summary.csv")
)

day12_area_model <- tryCatch(
  lmerTest::lmer(
    Area_log2 ~ Treatment * BBB_State + (1 | Sample_ID),
    data = filter(day12_area_data, !is.na(BBB_State)),
    REML = FALSE
  ),
  error = function(e) NULL
)

if (!is.null(day12_area_model)) {
  
  day12_area_anova <- as.data.frame(anova(day12_area_model)) %>%
    tibble::rownames_to_column("Effect")
  
  write_csv(
    day12_area_anova,
    file.path(day12_area_dir, "Day12_Area_Treatment_BBB_Interaction_ANOVA.csv")
  )
  
  day12_area_emm <- emmeans::emmeans(
    day12_area_model, ~ Treatment | BBB_State,
    lmer.df = day12_df_method
  )
  
  day12_area_contrasts <- as.data.frame(
    summary(
      emmeans::contrast(
        day12_area_emm,
        method = list("Avastin - Saline" = c(-1, 1)),
        adjust = "none"
      ),
      infer = c(TRUE, TRUE)
    )
  ) %>%
    group_by(BBB_State) %>%
    mutate(
      p_value_Holm = p.adjust(p.value, method = "holm"),
      fold_change = 2^estimate,
      percent_change = (fold_change - 1) * 100
    ) %>%
    ungroup()
  
  write_csv(
    day12_area_contrasts,
    file.path(day12_area_dir, "Day12_Area_Treatment_Contrasts_By_BBB.csv")
  )
}

p_day12_area <- ggplot(
  day12_area_data,
  aes(Treatment, Area_log2, color = Treatment, fill = Treatment)
) +
  geom_violin(alpha = 0.14, trim = TRUE, linewidth = 0.6) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.25, linewidth = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.3, size = 1.2) +
  facet_wrap(~ BBB_State, scales = "free_y", nrow = 1) +
  scale_color_manual(values = day12_treatment_colors) +
  scale_fill_manual(values = day12_treatment_colors) +
  labs(
    title = "Day 12 calibrated vessel area",
    subtitle = "Vessel-level distributions; mouse-level summaries exported separately",
    x = NULL, y = "Calibrated vessel area [log2(x + 1)]"
  ) +
  day12_base_theme +
  theme(legend.position = "none")

day12_save_plot(
  p_day12_area, day12_area_dir,
  "Day12_Calibrated_Vessel_Area", 16, 8
)

# Additional separate vessel-area plot for each BBB state.

for (bbb in day12_bbb_order) {
  
  plot_dat <- day12_area_data %>%
    filter(as.character(BBB_State) == bbb)
  
  if (nrow(plot_dat) == 0) next
  
  p <- ggplot(
    plot_dat,
    aes(Treatment, Area_log2, color = Treatment, fill = Treatment)
  ) +
    geom_violin(alpha = 0.14, trim = TRUE, linewidth = 0.6) +
    geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.25, linewidth = 0.7) +
    geom_jitter(width = 0.12, alpha = 0.3, size = 1.2) +
    scale_color_manual(values = day12_treatment_colors) +
    scale_fill_manual(values = day12_treatment_colors) +
    labs(
      title = paste("Day 12 calibrated vessel area:", bbb),
      x = NULL, y = "Calibrated vessel area [log2(x + 1)]"
    ) +
    day12_base_theme +
    theme(legend.position = "none")
  
  day12_save_plot(
    p, day12_area_dir,
    paste0("Day12_Calibrated_Vessel_Area_", gsub(" ", "_", bbb)),
    8, 7
  )
}

# ---- 6. Perivascular IgG distributions and quantile cutoffs ----
# One record per sample/vessel prevents duplication across markers.
# BBB state is shown using the available vessel classification.
# Quantiles are calculated separately within each mouse/sample.

day12_igg_data <- day12_data %>%
  filter(is.finite(IgG_value)) %>%
  distinct(Sample_ID, Treatment, Vessel_ID, IgG_value, BBB_State)

day12_igg_quantiles <- day12_igg_data %>%
  group_by(Sample_ID, Treatment) %>%
  summarise(
    n_vessels = n(),
    IgG_Q25 = quantile(IgG_value, 0.25, na.rm = TRUE, names = FALSE),
    IgG_median = median(IgG_value, na.rm = TRUE),
    IgG_Q75 = quantile(IgG_value, 0.75, na.rm = TRUE, names = FALSE),
    IgG_IQR = IQR(IgG_value, na.rm = TRUE),
    IgG_mean = mean(IgG_value, na.rm = TRUE),
    IgG_SD = sd(IgG_value, na.rm = TRUE),
    IgG_min = min(IgG_value, na.rm = TRUE),
    IgG_max = max(IgG_value, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  day12_igg_quantiles,
  file.path(day12_igg_dir, "Day12_IgG_Sample_Quantiles.csv")
)

# Assign sample slots for a consistent 2-treatment x 3-sample layout.

day12_sample_map <- day12_igg_quantiles %>%
  arrange(Treatment, Sample_ID) %>%
  group_by(Treatment) %>%
  mutate(Sample_slot = paste0("Sample ", row_number())) %>%
  ungroup()

if (any(table(day12_sample_map$Treatment) != 3)) {
  warning("At least one treatment does not have exactly three samples.")
}

day12_igg_data <- day12_igg_data %>%
  left_join(
    day12_sample_map %>% select(Treatment, Sample_ID, Sample_slot),
    by = c("Treatment", "Sample_ID")
  ) %>%
  mutate(
    Sample_slot = factor(
      Sample_slot, levels = c("Sample 1", "Sample 2", "Sample 3")
    )
  )

day12_igg_quantiles <- day12_igg_quantiles %>%
  left_join(
    day12_sample_map %>% select(Treatment, Sample_ID, Sample_slot),
    by = c("Treatment", "Sample_ID")
  ) %>%
  mutate(
    Sample_slot = factor(
      Sample_slot, levels = c("Sample 1", "Sample 2", "Sample 3")
    )
  )

day12_igg_data <- day12_igg_data %>%
  left_join(
    day12_igg_quantiles %>%
      select(Sample_ID, Treatment, IgG_Q25, IgG_median, IgG_Q75),
    by = c("Sample_ID", "Treatment")
  ) %>%
  mutate(
    BBB_group_quantile = case_when(
      IgG_value <= IgG_Q25 ~ "High BBB",
      IgG_value >= IgG_Q75 ~ "Low BBB",
      TRUE ~ "Mid BBB"
    ),
    BBB_group_quantile = factor(
      BBB_group_quantile, levels = day12_bbb_order
    )
  )

day12_igg_group_counts <- day12_igg_data %>%
  count(Treatment, Sample_ID, BBB_group_quantile, name = "n_vessels")

write_csv(
  day12_igg_group_counts,
  file.path(day12_igg_dir, "Day12_IgG_Quantile_Group_Counts.csv")
)

# Quantile line data: one set of cutoffs per sample.

day12_igg_cutoff_lines <- day12_igg_quantiles %>%
  select(Treatment, Sample_ID, Sample_slot, IgG_Q25, IgG_median, IgG_Q75) %>%
  pivot_longer(
    cols = c(IgG_Q25, IgG_median, IgG_Q75),
    names_to = "Quantile",
    values_to = "Cutoff"
  ) %>%
  mutate(
    Quantile = recode(
      Quantile,
      IgG_Q25 = "Q25",
      IgG_median = "Median",
      IgG_Q75 = "Q75"
    ),
    Quantile = factor(Quantile, levels = c("Q25", "Median", "Q75"))
  )

# Plot 1: density distributions, BBB-state colors, and dotted quantile lines.
# Facets use Treatment x Sample_slot, arranged as 2 rows x 3 columns.

p_day12_igg_density <- ggplot(
  day12_igg_data,
  aes(x = IgG_value)
) +
  geom_density(
    aes(color = BBB_group_quantile, fill = BBB_group_quantile),
    alpha = 0.28,
    linewidth = 1.1,
    adjust = 1,
    na.rm = TRUE
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "black",
    linewidth = 0.9
  ) +
  geom_vline(
    data = day12_igg_quantiles,
    aes(xintercept = IgG_Q25),
    inherit.aes = FALSE,
    color = "#0072B2",
    linetype = "dotted",
    linewidth = 1.1
  ) +
  geom_vline(
    data = day12_igg_quantiles,
    aes(xintercept = IgG_median),
    inherit.aes = FALSE,
    color = "black",
    linetype = "dotted",
    linewidth = 1.1
  ) +
  geom_vline(
    data = day12_igg_quantiles,
    aes(xintercept = IgG_Q75),
    inherit.aes = FALSE,
    color = "#D55E00",
    linetype = "dotted",
    linewidth = 1.1
  ) +
  facet_grid(Treatment ~ Sample_slot, scales = "free") +
  scale_color_manual(values = day12_bbb_colors, drop = FALSE) +
  scale_fill_manual(
    values = day12_bbb_colors,
    breaks = day12_bbb_order,
    labels = c("High BBB stability", "Mid BBB stability", "Low BBB stability"),
    drop = FALSE
  ) +
  labs(
    title = "Day 12 perivascular IgG density distributions",
    subtitle = paste(
      "Each panel represents one mouse.",
      "Blue, black, and red-orange dotted lines: Q25, median, and Q75."
    ),
    x = "Perivascular IgG mean (background-corrected)",
    y = "Density",
    color = "BBB state",
    fill = "BBB state"
  ) +
  day12_base_theme +
  theme(legend.position = "bottom")

day12_save_plot(
  p_day12_igg_density, day12_igg_dir,
  "Day12_IgG_Density_BBB_Quantiles_2x3", 18, 12
)

# Plot 2: histogram distributions with sample-specific Q25 and Q75 lines.

p_day12_igg_hist <- ggplot(day12_igg_data, aes(x = IgG_value)) +
  geom_histogram(
    bins = 35, fill = "grey75", color = "white", linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.9
  ) +
  geom_vline(
    data = day12_igg_quantiles,
    aes(xintercept = IgG_Q25),
    inherit.aes = FALSE, color = "#0072B2",
    linetype = "dotted", linewidth = 1.1
  ) +
  geom_vline(
    data = day12_igg_quantiles,
    aes(xintercept = IgG_median),
    inherit.aes = FALSE, color = "black",
    linetype = "dotted", linewidth = 1.1
  ) +
  geom_vline(
    data = day12_igg_quantiles,
    aes(xintercept = IgG_Q75),
    inherit.aes = FALSE, color = "#D55E00",
    linetype = "dotted", linewidth = 1.1
  ) +
  facet_grid(Treatment ~ Sample_slot, scales = "free_y") +
  labs(
    title = "Day 12 perivascular IgG histograms",
    subtitle = "Dotted lines show sample-specific Q25, median, and Q75",
    x = "Perivascular IgG mean (background-corrected)",
    y = "Number of vessels"
  ) +
  day12_base_theme

day12_save_plot(
  p_day12_igg_hist, day12_igg_dir,
  "Day12_IgG_Histograms_2x3", 18, 12
)

# Plot 3: sample-specific quartile cutoff dot plot.

p_day12_igg_cutoffs <- ggplot(
  day12_igg_cutoff_lines,
  aes(x = Cutoff, y = Sample_ID, color = Quantile, shape = Quantile)
) +
  geom_vline(
    xintercept = 0, linetype = "dotted", color = "grey35", linewidth = 0.9
  ) +
  geom_point(size = 5, position = position_dodge(width = 0.45)) +
  facet_grid(Treatment ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = day12_quantile_colors) +
  scale_shape_manual(values = c(Q25 = 16, Median = 18, Q75 = 17)) +
  labs(
    title = "Day 12 sample-specific IgG quantile cutoffs",
    subtitle = "Cutoffs calculated independently within each mouse",
    x = "Perivascular IgG mean (background-corrected)",
    y = NULL, color = "Quantile", shape = "Quantile"
  ) +
  day12_base_theme +
  theme(legend.position = "bottom")

day12_save_plot(
  p_day12_igg_cutoffs, day12_igg_dir,
  "Day12_IgG_Quantile_Cutoffs", 13, 9
)

# ---- 7. Mouse-level IgG variability: IQR and exact permutation ----

day12_igg_iqr <- day12_igg_quantiles %>%
  select(Sample_ID, Treatment, n_vessels, IgG_mean, IgG_median, IgG_SD, IgG_IQR)

write_csv(
  day12_igg_iqr,
  file.path(day12_igg_spread_dir, "Day12_IgG_Variability_By_Mouse.csv")
)

day12_igg_saline <- day12_igg_iqr %>%
  filter(Treatment == "Saline") %>% pull(IgG_IQR)

day12_igg_avastin <- day12_igg_iqr %>%
  filter(Treatment == "Avastin") %>% pull(IgG_IQR)

day12_igg_perm <- tibble()

if (length(day12_igg_saline) == 3 && length(day12_igg_avastin) == 3) {
  
  day12_igg_all_iqr <- c(day12_igg_saline, day12_igg_avastin)
  day12_igg_observed_diff <- mean(day12_igg_avastin) - mean(day12_igg_saline)
  day12_igg_assignments <- combn(seq_along(day12_igg_all_iqr), 3)
  
  day12_igg_perm_diffs <- apply(day12_igg_assignments, 2, function(idx) {
    mean(day12_igg_all_iqr[idx]) - mean(day12_igg_all_iqr[-idx])
  })
  
  day12_igg_perm <- tibble(
    n_Saline = length(day12_igg_saline),
    n_Avastin = length(day12_igg_avastin),
    Mean_IQR_Saline = mean(day12_igg_saline),
    Mean_IQR_Avastin = mean(day12_igg_avastin),
    Median_IQR_Saline = median(day12_igg_saline),
    Median_IQR_Avastin = median(day12_igg_avastin),
    Difference_Avastin_minus_Saline = day12_igg_observed_diff,
    Ratio_Avastin_to_Saline = mean(day12_igg_avastin) / mean(day12_igg_saline),
    Percent_change = (
      mean(day12_igg_avastin) / mean(day12_igg_saline) - 1
    ) * 100,
    Exact_permutation_p_two_sided = mean(
      abs(day12_igg_perm_diffs) >= abs(day12_igg_observed_diff) - 1e-12
    ),
    Exact_permutation_p_one_sided_lower = mean(
      day12_igg_perm_diffs <= day12_igg_observed_diff + 1e-12
    )
  )
}

write_csv(
  day12_igg_perm,
  file.path(day12_igg_spread_dir, "Day12_IgG_IQR_Exact_Permutation.csv")
)

p_day12_igg_iqr <- ggplot(
  day12_igg_iqr, aes(Treatment, IgG_IQR, color = Treatment)
) +
  geom_point(size = 5) +
  scale_color_manual(values = day12_treatment_colors) +
  labs(
    title = "Day 12 perivascular IgG variability",
    subtitle = "Each point represents one mouse; IQR measures within-mouse dispersion",
    x = NULL, y = "Within-mouse IgG interquartile range"
  ) +
  day12_base_theme +
  theme(legend.position = "none")

day12_save_plot(
  p_day12_igg_iqr, day12_igg_spread_dir,
  "Day12_IgG_IQR_By_Mouse", 9, 7
)

# ---- 8. PLVAP within-mouse variability ----

day12_plvap_data <- day12_data %>%
  filter(
    as.character(Marker) == "PLVAP",
    is.finite(Expression_per_area_log2),
    !is.na(BBB_State)
  ) %>%
  mutate(BBB_State = factor(as.character(BBB_State), levels = day12_bbb_order))

day12_plvap_spread <- day12_plvap_data %>%
  group_by(Sample_ID, Treatment, BBB_State) %>%
  summarise(
    n_vessels = n(),
    PLVAP_mean = mean(Expression_per_area_log2, na.rm = TRUE),
    PLVAP_median = median(Expression_per_area_log2, na.rm = TRUE),
    PLVAP_SD = sd(Expression_per_area_log2, na.rm = TRUE),
    PLVAP_IQR = IQR(Expression_per_area_log2, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  day12_plvap_spread,
  file.path(day12_plvap_dir, "Day12_PLVAP_Within_Mouse_Spread.csv")
)

p_day12_plvap_sd <- ggplot(
  day12_plvap_spread, aes(Treatment, PLVAP_SD, color = Treatment)
) +
  geom_point(size = 5) +
  facet_wrap(~ BBB_State, scales = "free_y", nrow = 1) +
  scale_color_manual(values = day12_treatment_colors) +
  labs(
    title = "PLVAP variability within each mouse",
    x = NULL, y = "Within-mouse SD of PLVAP expression [log2(x + 1)]"
  ) +
  day12_base_theme +
  theme(legend.position = "none")

p_day12_plvap_iqr <- ggplot(
  day12_plvap_spread, aes(Treatment, PLVAP_IQR, color = Treatment)
) +
  geom_point(size = 5) +
  facet_wrap(~ BBB_State, scales = "free_y", nrow = 1) +
  scale_color_manual(values = day12_treatment_colors) +
  labs(
    title = "PLVAP interquartile range within each mouse",
    x = NULL, y = "Within-mouse PLVAP IQR [log2(x + 1)]"
  ) +
  day12_base_theme +
  theme(legend.position = "none")

day12_save_plot(
  p_day12_plvap_sd, day12_plvap_dir,
  "Day12_PLVAP_Within_Mouse_SD", 16, 8
)

day12_save_plot(
  p_day12_plvap_iqr, day12_plvap_dir,
  "Day12_PLVAP_Within_Mouse_IQR", 16, 8
)

# ---- 9. Helper functions for matched-vessel correlations ----

day12_make_vessel_wide <- function(data, marker_x, marker_y, low_bbb_only = FALSE) {
  
  dat <- data %>%
    filter(
      as.character(Marker) %in% c(marker_x, marker_y),
      is.finite(Expression_per_area_log2)
    )
  
  if (low_bbb_only) {
    dat <- dat %>% filter(as.character(BBB_State) == "Low BBB")
  }
  
  dat %>%
    group_by(Sample_ID, Treatment, Vessel_ID, Marker) %>%
    summarise(Expression = mean(Expression_per_area_log2, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(Marker = as.character(Marker)) %>%
    pivot_wider(names_from = Marker, values_from = Expression) %>%
    filter(
      is.finite(.data[[marker_x]]),
      is.finite(.data[[marker_y]])
    )
}

day12_cor_test <- function(x, y, method = "spearman") {
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3 || dplyr::n_distinct(x) < 2 ||
      dplyr::n_distinct(y) < 2) {
    return(tibble(n = length(x), correlation = NA_real_, p_value = NA_real_))
  }
  
  test <- tryCatch(
    cor.test(x, y, method = method, exact = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(test)) {
    return(tibble(n = length(x), correlation = NA_real_, p_value = NA_real_))
  }
  
  tibble(
    n = length(x),
    correlation = unname(test$estimate),
    p_value = test$p.value
  )
}

day12_run_correlation <- function(
    marker_x, marker_y, output_dir, analysis_label, low_bbb_only = FALSE) {
  
  wide <- day12_make_vessel_wide(
    day12_data, marker_x, marker_y, low_bbb_only
  )
  
  if (nrow(wide) == 0) {
    message("No matched vessels for correlation: ", analysis_label)
    return(invisible(NULL))
  }
  
  pooled_stats <- wide %>%
    group_by(Treatment) %>%
    group_modify(~ {
      day12_cor_test(.x[[marker_x]], .x[[marker_y]]) %>%
        rename(n_vessels = n, spearman_rho = correlation, spearman_p = p_value)
    }) %>%
    ungroup()
  
  write_csv(
    pooled_stats,
    file.path(output_dir, paste0("Day12_", analysis_label, "_Pooled_Correlations.csv"))
  )
  
  mouse_stats <- wide %>%
    group_by(Sample_ID, Treatment) %>%
    group_modify(~ {
      day12_cor_test(.x[[marker_x]], .x[[marker_y]]) %>%
        rename(n_vessels = n, spearman_rho = correlation, spearman_p = p_value)
    }) %>%
    ungroup()
  
  write_csv(
    mouse_stats,
    file.path(output_dir, paste0("Day12_", analysis_label, "_Mouse_Level_Correlations.csv"))
  )
  
  label_data <- pooled_stats %>%
    mutate(
      label = paste0(
        "n = ", n_vessels,
        "\nSpearman \u03c1 = ",
        ifelse(is.na(spearman_rho), "NA", sprintf("%.2f", spearman_rho)),
        "\np = ",
        ifelse(is.na(spearman_p), "NA",
               ifelse(spearman_p < 0.001, "<0.001", sprintf("%.3f", spearman_p)))
      )
    )
  
  p <- ggplot(
    wide,
    aes(x = .data[[marker_x]], y = .data[[marker_y]], color = Sample_ID)
  ) +
    geom_point(alpha = 0.45, size = 2.5) +
    geom_smooth(
      aes(group = Treatment), method = "lm", se = TRUE,
      color = "black", linewidth = 1
    ) +
    geom_label(
      data = label_data,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = -0.05, vjust = 1.1,
      size = 4.8, family = "Helvetica", color = "black", fill = "white"
    ) +
    facet_wrap(~ Treatment, nrow = 1) +
    scale_color_discrete() +
    labs(
      title = paste0("Day 12 ", marker_x, "–", marker_y, " relationship"),
      subtitle = if (low_bbb_only) {
        "Low-BBB vessels only; pooled-vessel statistics are descriptive"
      } else {
        "All BBB states; pooled-vessel statistics are descriptive"
      },
      x = paste0(marker_x, " expression per area [log2(x + 1)]"),
      y = paste0(marker_y, " expression per area [log2(x + 1)]"),
      color = "Mouse/sample"
    ) +
    day12_base_theme +
    theme(legend.position = "bottom")
  
  day12_save_plot(
    p, output_dir, paste0("Day12_", analysis_label, "_Scatterplot"), 16, 10
  )
  
  # Additional sample-faceted scatterplot: 2 columns x 3 rows.
  # Each panel is one mouse; treatment is retained in the strip label.
  
  wide <- wide %>%
    left_join(
      day12_sample_map %>% select(Sample_ID, Treatment, Sample_slot),
      by = c("Sample_ID", "Treatment")
    ) %>%
    mutate(
      Sample_label = paste0(as.character(Treatment), "\n", Sample_ID)
    )
  
  p_sample <- ggplot(
    wide,
    aes(x = .data[[marker_x]], y = .data[[marker_y]], color = Treatment)
  ) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_smooth(
      method = "lm", se = FALSE, color = "black", linewidth = 1
    ) +
    facet_wrap(~ Sample_label, ncol = 2, scales = "free") +
    scale_color_manual(values = day12_treatment_colors) +
    labs(
      title = paste0("Day 12 ", marker_x, "–", marker_y, " by mouse"),
      subtitle = if (low_bbb_only) "Low-BBB vessels only" else "All BBB states",
      x = paste0(marker_x, " expression per area [log2(x + 1)]"),
      y = paste0(marker_y, " expression per area [log2(x + 1)]"),
      color = "Treatment"
    ) +
    day12_base_theme +
    theme(legend.position = "top")
  
  day12_save_plot(
    p_sample, output_dir,
    paste0("Day12_", analysis_label, "_Scatterplot_By_Sample_2col"),
    14, 18
  )
  
  invisible(list(
    vessel_data = wide,
    pooled_stats = pooled_stats,
    mouse_stats = mouse_stats,
    plot = p,
    sample_plot = p_sample
  ))
}

# ---- 10. MFSD2A–PLVAP correlation in Low-BBB vessels ----

day12_mfsd_plvap_results <- day12_run_correlation(
  marker_x = "MFSD2A",
  marker_y = "PLVAP",
  output_dir = day12_mfsd_plvap_dir,
  analysis_label = "MFSD2A_PLVAP_LowBBB",
  low_bbb_only = TRUE
)

# ---- 11. GLUT1–MFSD2A correlation across all BBB states ----

day12_glut_mfsd_results <- day12_run_correlation(
  marker_x = "GLUT1",
  marker_y = "MFSD2A",
  output_dir = day12_glut_mfsd_dir,
  analysis_label = "GLUT1_MFSD2A_All_Vessels",
  low_bbb_only = FALSE
)

# ---- 12. Export final analysis index ----

day12_analysis_summary <- tibble(
  Analysis = c(
    "Treatment mixed models by BBB state",
    "Sample-level Welch sensitivity analysis",
    "Sample-level sensitivity boxplots by BBB state",
    "Calibrated vessel area and separate BBB-state plots",
    "IgG distributions, density groups, and sample-specific quantiles",
    "Mouse-level IgG IQR and exact permutation",
    "PLVAP within-mouse SD and IQR",
    "MFSD2A-PLVAP correlation in Low-BBB vessels",
    "GLUT1-MFSD2A correlation across all BBB states",
    "Sample-faceted correlation scatterplots"
  ),
  Output_folder = c(
    day12_results_dir,
    day12_sensitivity_dir,
    day12_sensitivity_dir,
    day12_area_dir,
    day12_igg_dir, 
    day12_igg_spread_dir,
    day12_plvap_dir,
    day12_mfsd_plvap_dir,
    day12_glut_mfsd_dir,
    day12_root_dir
  )
)

write_csv(
  day12_analysis_summary,
  file.path(day12_root_dir, "Day12_Analysis_Output_Index.csv")
)

message("Day 12 analyses completed.")
message("Results saved under: ", normalizePath(day12_root_dir))
message("Mixed-model inference used: ", day12_df_method)



#IgG Density plot clean----

# ---- Day 12: IgG density by treatment, tumor, and BBB state ----

library(tidyverse)
library(ggplot2)

# ---- Output folder ----

igg_density_dir <- file.path(
  "Arpan Multiplex",
  "Day12_Area_Normalized",
  "IgG_Quantile_Distributions"
)

dir.create(
  igg_density_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ---- Prepare IgG data ----

igg_data <- day12_all_vessels_long %>%
  mutate(
    Treatment = case_when(
      str_detect(str_to_lower(as.character(Treatment)), "avastin") ~ "Avastin",
      str_detect(str_to_lower(as.character(Treatment)), "saline|pbs|vehicle") ~ "Saline",
      TRUE ~ as.character(Treatment)
    ),
    Sample_ID = as.character(Sample_ID),
    IgG_value = as.numeric(IgG_PV_Mean_bg_corrected)
  ) %>%
  filter(
    Treatment %in% c("Saline", "Avastin"),
    !is.na(Sample_ID),
    is.finite(IgG_value)
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("Saline", "Avastin")
    )
  )

# ---- Assign tumor labels ----

sample_lookup <- igg_data %>%
  distinct(Treatment, Sample_ID) %>%
  arrange(Treatment, Sample_ID) %>%
  group_by(Treatment) %>%
  mutate(
    Sample_slot = paste("Tumor", row_number())
  ) %>%
  ungroup()

igg_data <- igg_data %>%
  left_join(
    sample_lookup,
    by = c("Treatment", "Sample_ID")
  ) %>%
  mutate(
    Sample_slot = factor(
      Sample_slot,
      levels = c("Tumor 1", "Tumor 2", "Tumor 3")
    )
  )

# ---- Calculate sample-specific BBB cutoffs ----

sample_cutoffs <- igg_data %>%
  group_by(Treatment, Sample_ID, Sample_slot) %>%
  summarise(
    n_vessels = n(),
    Q25 = quantile(IgG_value, 0.25, na.rm = TRUE),
    Median = median(IgG_value, na.rm = TRUE),
    Q75 = quantile(IgG_value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# ---- Assign BBB state within each tumor ----

igg_data <- igg_data %>%
  left_join(
    sample_cutoffs %>%
      select(Treatment, Sample_ID, Q25, Q75),
    by = c("Treatment", "Sample_ID")
  ) %>%
  mutate(
    BBB_group = case_when(
      IgG_value <= Q25 ~ "High",
      IgG_value >= Q75 ~ "Low",
      TRUE ~ "Mid"
    ),
    BBB_group = factor(
      BBB_group,
      levels = c("High", "Mid", "Low")
    )
  )

# ---- Save cutoff and count tables ----

write_csv(
  sample_cutoffs,
  file.path(igg_density_dir, "IgG_sample_specific_quantiles.csv")
)

bbb_counts <- igg_data %>%
  count(Treatment, Sample_ID, Sample_slot, BBB_group, name = "n_vessels") %>%
  arrange(Treatment, Sample_slot, BBB_group)

write_csv(
  bbb_counts,
  file.path(igg_density_dir, "IgG_BBB_state_counts.csv")
)

# ---- BBB state colors ----

bbb_colors <- c(
  "High" = "#559985",
  "Mid" = "#B89A3E",
  "Low" = "#B8794F"
)

# ---- Density plot: six panels, three BBB states per panel ----

igg_density_by_state <- ggplot(
  igg_data,
  aes(
    x = IgG_value,
    color = BBB_group,
    fill = BBB_group
  )
) +
  geom_density(
    aes(group = BBB_group),
    alpha = 0.25,
    linewidth = 1.5,
    adjust = 1,
    na.rm = TRUE
  ) +
  geom_vline(
    data = sample_cutoffs,
    aes(xintercept = Q25),
    inherit.aes = FALSE,
    color = "#0072B2",
    linetype = "dashed",
    linewidth = 1.1
  ) +
  geom_vline(
    data = sample_cutoffs,
    aes(xintercept = Q75),
    inherit.aes = FALSE,
    color = "#D55E00",
    linetype = "dashed",
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dotted",
    linewidth = 1.05
  ) +
  facet_grid(
    rows = vars(Treatment),
    cols = vars(Sample_slot),
    labeller = labeller(
      Treatment = c(
        "Saline" = "Vehicle",
        "Avastin" = "Avastin"
      )
    )
  ) +
  scale_color_manual(
    values = bbb_colors,
    name = "BBB state"
  ) +
  scale_fill_manual(
    values = bbb_colors,
    name = "BBB state"
  ) +
  labs(
    title = "Day 12: Perivascular IgG Density by BBB State",
    subtitle = "BBB states defined using tumor-specific Q25 and Q75 cutoffs",
    x = "Perivascular IgG (background-corrected mean)",
    y = "Density"
  ) +
  theme_bw(
    base_size = 28,
    base_family = "Helvetica"
  ) +
  theme(
    plot.title = element_text(
      size = 36,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    plot.subtitle = element_text(
      size = 24,
      hjust = 0.5,
      color = "black"
    ),
    axis.title = element_text(
      size = 26,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(
      size = 20,
      color = "black"
    ),
    strip.text = element_text(
      size = 24,
      face = "bold",
      color = "black"
    ),
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 1.1
    ),
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 1.1
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.4
    ),
    panel.spacing = grid::unit(1, "lines"),
    legend.title = element_text(
      size = 22,
      face = "bold"
    ),
    legend.text = element_text(size = 20),
    legend.position = "bottom",
    plot.background = element_rect(
      fill = "white",
      color = NA
    )
  )

# ---- Export PNG ----

ggsave(
  filename = file.path(
    igg_density_dir,
    "Day12_IgG_Density_By_BBB_State.png"
  ),
  plot = igg_density_by_state,
  width = 18,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ---- Export SVG ----

ggsave(
  filename = file.path(
    igg_density_dir,
    "Day12_IgG_Density_By_BBB_State.svg"
  ),
  plot = igg_density_by_state,
  width = 18,
  height = 14,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# ---- Display plot and summaries ----

igg_density_by_state

sample_cutoffs %>%
  arrange(Treatment, Sample_slot)

bbb_counts
