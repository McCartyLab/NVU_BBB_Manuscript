
#Day 5 Area calibrated 

# ---- Setup: Build QC-passed Day 5 analysis data ----

library(tidyverse)

set.seed(123)

# ---- Output folders and analysis settings 

day5_new_analysis_dir <- file.path(
  "Arpan Multiplex",
  "Day5_Area_Normalized"
)

day5_setup_results_dir <- file.path(
  day5_new_analysis_dir,
  "Setup_and_QC"
)

dir.create(
  day5_setup_results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

day5_igg_file <- file.path(
  "Arpan Multiplex", "IgG Data", "IgG_data_processed.csv"
)

day5_marker_cols <- c(
  "MFSD2A_Mean_bg_corrected",
  "ABCB1_Mean_bg_corrected",
  "GLUT1_Mean_bg_corrected",
  "PLVAP_Mean_bg_corrected"
)

day5_marker_labels <- c(
  "MFSD2A_Mean_bg_corrected" = "MFSD2A",
  "ABCB1_Mean_bg_corrected" = "ABCB1",
  "GLUT1_Mean_bg_corrected" = "GLUT1",
  "PLVAP_Mean_bg_corrected" = "PLVAP"
)

day5_marker_order <- c("MFSD2A", "ABCB1", "GLUT1", "PLVAP")

day5_treatment_levels <- c("Saline", "Avastin")

day5_treatment_colors <- c(
  "Saline" = "#0072B2",
  "Avastin" = "#D55E00"
)

day5_bbb_levels <- c("High BBB", "Mid BBB", "Low BBB")

day5_excluded_vessel_ids <- "V1052"

day5_igg_required_cols <- c(
  "Sample_ID",
  "Vessel_ID",
  "IgG_PV_Mean_bg_corrected",
  "IgG_sign",
  "IgG_qualifier",
  "BBB Stability",
  "Vessel_centroid_inside_tumor_ROI",
  "PV_QC_status"
)

# Each sample corresponds to one mouse/tumor ROI.

day5_gene_files <- tibble(
  Sample_ID = c(
    "Saline-Tumor-1-ROI-1",
    "Saline-Tumor-2-ROI-2",
    "Saline-Tumor-3-ROI-3",
    "Avastin-Tumor-1-ROI-1",
    "Avastin-Tumor-2-ROI-2",
    "Avastin-Tumor-3-ROI-3"
  ),
  Treatment = c(
    "Saline", "Saline", "Saline",
    "Avastin", "Avastin", "Avastin"
  ),
  Tumor = c("T1", "T2", "T3", "T1", "T2", "T3"),
  ROI = c("ROI-1", "ROI-2", "ROI-3", "ROI-1", "ROI-2", "ROI-3"),
  File = c(
    file.path(
      "Arpan Multiplex",
      "Saline_Day5_Tumor_1_ROI_1_Vessel_measurements.csv"
    ),
    file.path(
      "Arpan Multiplex",
      "Saline_Day5_Tumor_2_ROI_2_Vessel_measurements.csv"
    ),
    file.path(
      "Arpan Multiplex",
      "Saline_Day5_Tumor_3_ROI_3_Vessel_measurements.csv"
    ),
    file.path(
      "Arpan Multiplex",
      "Avastin_Day5_Tumor_1_ROI_1_Vessel_measurements.csv"
    ),
    file.path(
      "Arpan Multiplex",
      "Avastin_Day5_Tumor_2_ROI_2_Vessel_measurements.csv"
    ),
    file.path(
      "Arpan Multiplex",
      "Avastin_Day5_Tumor_3_ROI_3_Analysis_Vessel_measurements.csv"
    )
  )
)

# ---- Validate input files and columns 

if (!file.exists(day5_igg_file)) {
  stop("IgG input file not found: ", day5_igg_file)
}

igg_raw <- readr::read_csv(
  day5_igg_file,
  show_col_types = FALSE
)

missing_igg_cols <- setdiff(day5_igg_required_cols, names(igg_raw))

if (length(missing_igg_cols) > 0) {
  stop(
    "Missing required columns in IgG file: ",
    paste(missing_igg_cols, collapse = ", ")
  )
}

# Area_calibrated is required for the new area-related calculations.

required_gene_cols <- c(
  "Vessel_ID",
  "Area_calibrated",
  day5_marker_cols
)

gene_file_checks <- day5_gene_files %>%
  mutate(
    file_exists = file.exists(File),
    missing_columns = NA_character_,
    file_check_pass = FALSE
  )

for (i in seq_len(nrow(gene_file_checks))) {
  
  if (!gene_file_checks$file_exists[i]) {
    gene_file_checks$missing_columns[i] <- "FILE_NOT_FOUND"
    next
  }
  
  header_data <- tryCatch(
    readr::read_csv(
      gene_file_checks$File[i],
      n_max = 0,
      show_col_types = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(header_data)) {
    gene_file_checks$missing_columns[i] <- "FILE_READ_ERROR"
    next
  }
  
  missing_cols <- setdiff(required_gene_cols, names(header_data))
  
  gene_file_checks$missing_columns[i] <- if (length(missing_cols) == 0) {
    ""
  } else {
    paste(missing_cols, collapse = "; ")
  }
  
  gene_file_checks$file_check_pass[i] <- length(missing_cols) == 0
}

write_csv(
  gene_file_checks,
  file.path(
    day5_setup_results_dir,
    "gene_expression_input_file_checks.csv"
  )
)

if (any(!gene_file_checks$file_check_pass)) {
  stop(
    "One or more gene-expression files are missing or lack required columns. ",
    "See Setup_and_QC/gene_expression_input_file_checks.csv."
  )
}

# ---- Standardize identifiers and check duplicate IgG records ----

igg_raw <- igg_raw %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    Vessel_ID = as.character(Vessel_ID),
    IgG_PV_Mean_bg_corrected = suppressWarnings(
      as.numeric(IgG_PV_Mean_bg_corrected)
    )
  )

igg_duplicate_records <- igg_raw %>%
  count(Sample_ID, Vessel_ID, name = "n_records") %>%
  filter(n_records > 1)

if (nrow(igg_duplicate_records) > 0) {
  
  write_csv(
    igg_duplicate_records,
    file.path(
      day5_setup_results_dir,
      "IgG_duplicate_vessel_records.csv"
    )
  )
  
  stop(
    "Duplicate Sample_ID/Vessel_ID records were found in the IgG file. ",
    "Review Setup_and_QC/IgG_duplicate_vessel_records.csv."
  )
}

# ---- Import marker measurements and calibrated vessel area 

gene_raw <- purrr::map_dfr(
  seq_len(nrow(day5_gene_files)),
  function(i) {
    
    metadata <- day5_gene_files[i, ]
    
    readr::read_csv(
      metadata$File,
      show_col_types = FALSE
    ) %>%
      mutate(
        Sample_ID = metadata$Sample_ID,
        Treatment = metadata$Treatment,
        Tumor = metadata$Tumor,
        ROI = metadata$ROI,
        Vessel_ID = as.character(Vessel_ID),
        Area_calibrated = suppressWarnings(
          as.numeric(Area_calibrated)
        ),
        across(
          all_of(day5_marker_cols),
          ~ suppressWarnings(as.numeric(.x))
        )
      ) %>%
      select(
        Sample_ID,
        Treatment,
        Tumor,
        ROI,
        Vessel_ID,
        Area_calibrated,
        all_of(day5_marker_cols)
      )
  }
) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = day5_treatment_levels
    )
  )

# ---- Exclude duplicate vessel V1052

gene_excluded_vessels <- gene_raw %>%
  filter(Vessel_ID %in% day5_excluded_vessel_ids)

write_csv(
  gene_excluded_vessels,
  file.path(
    day5_setup_results_dir,
    "gene_expression_excluded_duplicate_vessels.csv"
  )
)

gene_data <- gene_raw %>%
  filter(!Vessel_ID %in% day5_excluded_vessel_ids)

gene_duplicate_records <- gene_data %>%
  count(Sample_ID, Vessel_ID, name = "n_records") %>%
  filter(n_records > 1)

if (nrow(gene_duplicate_records) > 0) {
  
  write_csv(
    gene_duplicate_records,
    file.path(
      day5_setup_results_dir,
      "gene_expression_duplicate_vessel_records.csv"
    )
  )
  
  stop(
    "Duplicate Sample_ID/Vessel_ID records remain in gene-expression data."
  )
}

# ---- Calculate area-related marker measurements----
#
# Expression remains the original recorded mean intensity.
#
# Expression_per_area:
#   Mean_bg_corrected / Area_calibrated
#
# Estimated_integrated_signal:
#   Mean_bg_corrected * Area_calibrated
#
# These are distinct metrics. The per-area ratio is retained for your
# requested normalization, while the estimated integrated signal is
# retained for comparison. Neither calculation supplies missing
# background subtraction.

invalid_area_records <- gene_data %>%
  filter(
    is.na(Area_calibrated) |
      !is.finite(Area_calibrated) |
      Area_calibrated <= 0
  )

write_csv(
  invalid_area_records,
  file.path(
    day5_setup_results_dir,
    "vessels_with_invalid_Area_calibrated.csv"
  )
)

gene_data <- gene_data %>%
  mutate(
    across(
      all_of(day5_marker_cols),
      ~ if_else(is.finite(.x), .x, NA_real_)
    )
  )

for (marker_col in day5_marker_cols) {
  
  marker_name <- unname(day5_marker_labels[marker_col])
  
  area_normalized_col <- paste0(marker_name, "_per_area")
  integrated_signal_col <- paste0(marker_name, "_integrated_signal")
  
  gene_data[[area_normalized_col]] <- ifelse(
    is.finite(gene_data[[marker_col]]) &
      is.finite(gene_data$Area_calibrated) &
      gene_data$Area_calibrated > 0,
    gene_data[[marker_col]] / gene_data$Area_calibrated,
    NA_real_
  )
  
  gene_data[[integrated_signal_col]] <- ifelse(
    is.finite(gene_data[[marker_col]]) &
      is.finite(gene_data$Area_calibrated) &
      gene_data$Area_calibrated > 0,
    gene_data[[marker_col]] * gene_data$Area_calibrated,
    NA_real_
  )
}

write_csv(
  gene_data,
  file.path(
    day5_setup_results_dir,
    "Day5_Gene_Measurements_With_Area_Metrics.csv"
  )
)

# ---- Apply existing IgG vessel QC

igg_excluded_vessels <- igg_raw %>%
  filter(Vessel_ID %in% day5_excluded_vessel_ids)

write_csv(
  igg_excluded_vessels,
  file.path(
    day5_setup_results_dir,
    "IgG_excluded_duplicate_vessels.csv"
  )
)

igg_qc <- igg_raw %>%
  filter(!Vessel_ID %in% day5_excluded_vessel_ids) %>%
  mutate(
    centroid_inside = tolower(
      stringr::str_trim(
        as.character(Vessel_centroid_inside_tumor_ROI)
      )
    ) %in% c("1", "1.0", "true"),
    
    pv_analyzable = tolower(
      stringr::str_trim(as.character(PV_QC_status))
    ) == "analyzable",
    
    QC_pass = centroid_inside & pv_analyzable
  )

igg_qc_excluded <- igg_qc %>%
  filter(!QC_pass)

write_csv(
  igg_qc_excluded,
  file.path(
    day5_setup_results_dir,
    "IgG_vessels_excluded_by_QC.csv"
  )
)

igg_qc_summary <- igg_qc %>%
  group_by(Sample_ID) %>%
  summarise(
    n_before_QC = n(),
    n_centroid_inside = sum(centroid_inside, na.rm = TRUE),
    n_PV_analyzable = sum(pv_analyzable, na.rm = TRUE),
    n_QC_pass = sum(QC_pass, na.rm = TRUE),
    n_QC_excluded = sum(!QC_pass, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  igg_qc_summary,
  file.path(
    day5_setup_results_dir,
    "IgG_vessel_QC_summary.csv"
  )
)

# Keep the original BBB Stability assignment for comparison.
# The new quantile-derived BBB state will be assigned separately.

igg_qc_passed <- igg_qc %>%
  filter(QC_pass) %>%
  select(
    Sample_ID,
    Vessel_ID,
    IgG_PV_Mean_bg_corrected,
    IgG_sign,
    IgG_qualifier,
    `BBB Stability`
  ) %>%
  rename(
    BBB_Stability_original = `BBB Stability`
  )

# ---- Check sample matching 

gene_samples <- sort(unique(gene_data$Sample_ID))
igg_samples <- sort(unique(igg_qc_passed$Sample_ID))

gene_samples_missing_igg <- tibble(
  Sample_ID = setdiff(gene_samples, igg_samples)
)

igg_samples_without_gene <- tibble(
  Sample_ID = setdiff(igg_samples, gene_samples)
)

write_csv(
  gene_samples_missing_igg,
  file.path(
    day5_setup_results_dir,
    "gene_samples_missing_IgG.csv"
  )
)

write_csv(
  igg_samples_without_gene,
  file.path(
    day5_setup_results_dir,
    "IgG_samples_without_gene_expression.csv"
  )
)

if (nrow(gene_samples_missing_igg) > 0) {
  stop(
    "One or more gene-expression samples have no QC-passed IgG match. ",
    "See Setup_and_QC/gene_samples_missing_IgG.csv."
  )
}

if (nrow(igg_samples_without_gene) > 0) {
  message(
    "Some QC-passed IgG samples have no gene-expression counterpart; ",
    "they will be omitted from the matched analysis."
  )
}

gene_vessels_without_igg <- gene_data %>%
  anti_join(
    igg_qc_passed,
    by = c("Sample_ID", "Vessel_ID")
  )

igg_vessels_without_gene <- igg_qc_passed %>%
  anti_join(
    gene_data,
    by = c("Sample_ID", "Vessel_ID")
  )

write_csv(
  gene_vessels_without_igg,
  file.path(
    day5_setup_results_dir,
    "gene_vessels_without_IgG_match.csv"
  )
)

write_csv(
  igg_vessels_without_gene,
  file.path(
    day5_setup_results_dir,
    "IgG_vessels_without_gene_match.csv"
  )
)

# ---- Match gene measurements to QC-passed IgG records 

day5_all_vessels <- gene_data %>%
  inner_join(
    igg_qc_passed,
    by = c("Sample_ID", "Vessel_ID")
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = day5_treatment_levels
    ),
    Sample_ID = as.character(Sample_ID),
    IgG_PV_Mean_bg_corrected = suppressWarnings(
      as.numeric(IgG_PV_Mean_bg_corrected)
    )
  )

if (nrow(day5_all_vessels) == 0) {
  stop("No gene-expression vessels matched QC-passed IgG records.")
}

if (anyDuplicated(day5_all_vessels[c("Sample_ID", "Vessel_ID")])) {
  stop(
    "Duplicate vessel records remain after joining gene-expression and IgG data."
  )
}

# ---- Assign BBB states using IgG quartiles within each sample
#
# Q1: lowest 25% of IgG signal -> High BBB
# Q2 and Q3: middle 50%         -> Mid BBB
# Q4: highest 25% of IgG signal -> Low BBB
#
# ntile() creates approximately equal-sized quartiles within each sample.
# If IgG values are tied at a boundary, tied values may be split across
# quartiles. This is documented in the exported data and should be kept
# in mind when interpreting borderline vessels.

day5_all_vessels <- day5_all_vessels %>%
  group_by(Sample_ID) %>%
  mutate(
    IgG_quartile = if_else(
      is.finite(IgG_PV_Mean_bg_corrected),
      as.integer(
        ntile(
          rank(
            IgG_PV_Mean_bg_corrected,
            ties.method = "first",
            na.last = "keep"
          ),
          4
        )
      ),
      NA_integer_
    ),
    BBB_State = case_when(
      IgG_quartile == 1L ~ "High BBB",
      IgG_quartile %in% c(2L, 3L) ~ "Mid BBB",
      IgG_quartile == 4L ~ "Low BBB",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(
    BBB_State = factor(
      BBB_State,
      levels = day5_bbb_levels
    )
  )

# ---- Export vessels that could not receive a quantile BBB state

vessels_missing_igg_signal <- day5_all_vessels %>%
  filter(is.na(BBB_State))

write_csv(
  vessels_missing_igg_signal,
  file.path(
    day5_setup_results_dir,
    "Vessels_missing_IgG_signal_for_BBB_quantiles.csv"
  )
)

# ---- Summarize sample-specific BBB-state assignments

day5_bbb_quantile_summary <- day5_all_vessels %>%
  filter(!is.na(BBB_State)) %>%
  count(
    Sample_ID,
    Treatment,
    BBB_State,
    name = "n_vessels"
  ) %>%
  group_by(Sample_ID) %>%
  mutate(
    sample_total = sum(n_vessels),
    percent_of_sample = 100 * n_vessels / sample_total
  ) %>%
  ungroup()

write_csv(
  day5_bbb_quantile_summary,
  file.path(
    day5_setup_results_dir,
    "BBB_quantile_state_counts_by_sample.csv"
  )
)

# ---- Export the complete matched dataset
#
# Includes all QC-passed IgG-matched vessels, even if an IgG signal is
# missing and a quantile-derived BBB state could not be assigned.

write_csv(
  day5_all_vessels,
  file.path(
    day5_setup_results_dir,
    "Day5_All_QC_Passed_IgG_Matched_Vessels.csv"
  )
)

# ---- Prepare long-format marker data ----
#
# Expression:
#   Original recorded Mean_bg_corrected value
#
# Expression_per_area:
#   Mean_bg_corrected / Area_calibrated
#
# Estimated_integrated_signal:
#   Mean_bg_corrected * Area_calibrated
#
# The per-area value is the requested ratio. It is not a background-
# corrected intensity, because a measured background mean was unavailable.

day5_all_vessels_long <- day5_all_vessels %>%
  pivot_longer(
    cols = all_of(day5_marker_cols),
    names_to = "Marker_column",
    values_to = "Expression"
  ) %>%
  mutate(
    Marker = unname(day5_marker_labels[Marker_column]),
    
    Expression = suppressWarnings(as.numeric(Expression)),
    
    Expression_per_area = if_else(
      is.finite(Expression) &
        is.finite(Area_calibrated) &
        Area_calibrated > 0,
      Expression / Area_calibrated,
      NA_real_
    ),
    
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
      is.finite(Expression_per_area) & Expression_per_area >= 0,
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
      levels = day5_marker_order
    ),
    
    Treatment = factor(
      Treatment,
      levels = day5_treatment_levels
    ),
    
    BBB_State = factor(
      BBB_State,
      levels = day5_bbb_levels
    )
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
    BBB_Stability_original,
    IgG_quartile,
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

write_csv(
  day5_all_vessels_long,
  file.path(
    day5_setup_results_dir,
    "Day5_All_QC_Passed_Vessels_Long_Format.csv"
  )
)

write.csv(day5_all_vessels_long, "Arpan Multiplex/Day5_Area_Normalized/day5_all_vessels_long.csv")

# ---- Summarize area and marker-value availability

day5_measurement_summary <- day5_all_vessels_long %>%
  group_by(Treatment, Sample_ID, Marker) %>%
  summarise(
    n_matched_vessels = n(),
    n_valid_area = sum(
      is.finite(Area_calibrated) & Area_calibrated > 0,
      na.rm = TRUE
    ),
    n_valid_original_mean = sum(
      is.finite(Expression),
      na.rm = TRUE
    ),
    n_valid_per_area = sum(
      is.finite(Expression_per_area),
      na.rm = TRUE
    ),
    n_valid_integrated_signal = sum(
      is.finite(Estimated_integrated_signal),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_csv(
  day5_measurement_summary,
  file.path(
    day5_setup_results_dir,
    "Marker_and_Area_Data_Summary.csv"
  )
)

# ---- Report output locations and basic QC totals

message("Day 5 setup and QC completed.")
message("Results directory: ", day5_setup_results_dir)
message("Matched vessels: ", nrow(day5_all_vessels))
message(
  "Vessels with a quantile-derived BBB state: ",
  sum(!is.na(day5_all_vessels$BBB_State))
)
message(
  "Vessels without a quantile-derived BBB state: ",
  sum(is.na(day5_all_vessels$BBB_State))
)

print(day5_bbb_quantile_summary)
print(day5_measurement_summary)


#check for normality----

# ---- Check normality of area-normalized gene expression

library(tidyverse)
library(patchwork)

# ---- Output folders

day5_area_dir <- file.path(
  "Arpan Multiplex",
  "Day5_Area_Normalized"
)

normality_results_dir <- file.path(
  day5_area_dir, "Expression_Normality", "Results"
)

normality_plots_dir <- file.path(
  day5_area_dir, "Expression_Normality", "Plots"
)

dir.create(normality_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(normality_plots_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Prepare area-normalized expression data 
# Use the original Mean_bg_corrected / Area_calibrated ratio.
# This checks distributions; it does not change the underlying data.

marker_names <- c("MFSD2A", "ABCB1", "GLUT1", "PLVAP")

area_expression_cols <- paste0(marker_names, "_per_area")

missing_cols <- setdiff(area_expression_cols, names(gene_data))

if (length(missing_cols) > 0) {
  stop(
    "These area-normalized columns are missing from gene_data: ",
    paste(missing_cols, collapse = ", "),
    ". Run the revised setup script first."
  )
}

expression_area_long <- gene_data %>%
  select(
    Sample_ID,
    Treatment,
    all_of(area_expression_cols)
  ) %>%
  pivot_longer(
    cols = all_of(area_expression_cols),
    names_to = "Marker_column",
    values_to = "Expression_per_area"
  ) %>%
  mutate(
    Marker = stringr::str_remove(Marker_column, "_per_area$"),
    Expression_per_area = suppressWarnings(
      as.numeric(Expression_per_area)
    ),
    Expression_log2 = if_else(
      is.finite(Expression_per_area) & Expression_per_area >= 0,
      log2(Expression_per_area + 1),
      NA_real_
    ),
    Treatment = factor(Treatment, levels = c("Saline", "Avastin")),
    Marker = factor(Marker, levels = marker_names)
  ) %>%
  filter(
    is.finite(Expression_per_area),
    Expression_per_area >= 0
  )

if (nrow(expression_area_long) == 0) {
  stop("No finite, nonnegative area-normalized expression values found.")
}

write_csv(
  expression_area_long,
  file.path(normality_results_dir, "Area_normalized_expression_values.csv")
)

# ---- Calculate skewness 

calculate_skewness <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) < 3 || sd(x) == 0) {
    return(NA_real_)
  }
  
  mean((x - mean(x))^3) / sd(x)^3
}

# ---- Summary statistics by sample and marker 

expression_summary <- expression_area_long %>%
  group_by(Treatment, Sample_ID, Marker) %>%
  summarise(
    n = n(),
    mean_raw = mean(Expression_per_area),
    sd_raw = sd(Expression_per_area),
    median_raw = median(Expression_per_area),
    IQR_raw = IQR(Expression_per_area),
    skewness_raw = calculate_skewness(Expression_per_area),
    mean_log2 = mean(Expression_log2, na.rm = TRUE),
    sd_log2 = sd(Expression_log2, na.rm = TRUE),
    median_log2 = median(Expression_log2, na.rm = TRUE),
    skewness_log2 = calculate_skewness(Expression_log2),
    .groups = "drop"
  )

write_csv(
  expression_summary,
  file.path(normality_results_dir, "Expression_summary_by_sample.csv")
)

# ---- Shapiro-Wilk tests 
# Shapiro-Wilk accepts 3 to 5000 observations.
# For groups above 5000, report that the test was not run.
# Use plots and skewness as well as p-values to assess distributions.

run_shapiro <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) < 3) {
    return(tibble(
      n = length(x),
      W = NA_real_,
      p_value = NA_real_,
      note = "Fewer than 3 observations"
    ))
  }
  
  if (length(unique(x)) < 2) {
    return(tibble(
      n = length(x),
      W = NA_real_,
      p_value = NA_real_,
      note = "No variation"
    ))
  }
  
  if (length(x) > 5000) {
    return(tibble(
      n = length(x),
      W = NA_real_,
      p_value = NA_real_,
      note = "More than 5000 observations; test not run"
    ))
  }
  
  result <- shapiro.test(x)
  
  tibble(
    n = length(x),
    W = unname(result$statistic),
    p_value = result$p.value,
    note = ""
  )
}

normality_by_sample <- expression_area_long %>%
  group_by(Treatment, Sample_ID, Marker) %>%
  group_modify(~ {
    raw_test <- run_shapiro(.x$Expression_per_area) %>%
      rename(
        n_raw = n,
        W_raw = W,
        p_raw = p_value,
        note_raw = note
      )
    
    log_test <- run_shapiro(.x$Expression_log2) %>%
      rename(
        n_log2 = n,
        W_log2 = W,
        p_log2 = p_value,
        note_log2 = note
      )
    
    bind_cols(raw_test, log_test)
  }) %>%
  ungroup()

write_csv(
  normality_by_sample,
  file.path(normality_results_dir, "Normality_tests_by_sample.csv")
)

# ---- Normality tests pooled by treatment and marker

normality_by_treatment <- expression_area_long %>%
  group_by(Treatment, Marker) %>%
  group_modify(~ {
    raw_test <- run_shapiro(.x$Expression_per_area) %>%
      rename(
        n_raw = n,
        W_raw = W,
        p_raw = p_value,
        note_raw = note
      )
    
    log_test <- run_shapiro(.x$Expression_log2) %>%
      rename(
        n_log2 = n,
        W_log2 = W,
        p_log2 = p_value,
        note_log2 = note
      )
    
    bind_cols(raw_test, log_test)
  }) %>%
  ungroup()

write_csv(
  normality_by_treatment,
  file.path(normality_results_dir, "Normality_tests_by_treatment.csv")
)

# ---- Histograms: original area-normalized expression

p_hist_raw <- ggplot(
  expression_area_long,
  aes(x = Expression_per_area)
) +
  geom_histogram(bins = 40, fill = "#3B75AF", color = "white") +
  facet_grid(Marker ~ Sample_ID, scales = "free") +
  labs(
    title = "Area-normalized Expression: Original Scale",
    x = "Mean_bg_corrected / Area_calibrated",
    y = "Number of vessels"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    strip.text = element_text(face = "bold", size = 8)
  )

# ---- Histograms: log2-transformed expression

p_hist_log2 <- ggplot(
  expression_area_long,
  aes(x = Expression_log2)
) +
  geom_histogram(bins = 40, fill = "#D55E00", color = "white") +
  facet_grid(Marker ~ Sample_ID, scales = "free") +
  labs(
    title = "Area-normalized Expression: log2(x + 1)",
    x = "log2(Area-normalized expression + 1)",
    y = "Number of vessels"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    strip.text = element_text(face = "bold", size = 8)
  )

ggsave(
  file.path(normality_plots_dir, "Area_expression_histograms_raw_vs_log2.png"),
  p_hist_raw / p_hist_log2,
  width = 19,
  height = 15,
  dpi = 300
)

# ---- Q-Q plots: original expression 
p_qq_raw <- ggplot(
  expression_area_long,
  aes(sample = Expression_per_area)
) +
  stat_qq(alpha = 0.35, size = 0.7) +
  stat_qq_line() +
  facet_grid(Marker ~ Sample_ID, scales = "free") +
  labs(
    title = "Q-Q Plots: Original Area-normalized Expression",
    x = "Theoretical quantiles",
    y = "Observed quantiles"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 8)
  )

# ---- Q-Q plots: log2-transformed expression

p_qq_log2 <- ggplot(
  expression_area_long,
  aes(sample = Expression_log2)
) +
  stat_qq(alpha = 0.35, size = 0.7) +
  stat_qq_line() +
  facet_grid(Marker ~ Sample_ID, scales = "free") +
  labs(
    title = "Q-Q Plots: log2-transformed Area-normalized Expression",
    x = "Theoretical quantiles",
    y = "Observed quantiles"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 8)
  )

ggsave(
  file.path(normality_plots_dir, "Area_expression_QQ_raw_vs_log2.png"),
  p_qq_raw / p_qq_log2,
  width = 19,
  height = 15,
  dpi = 300
)

# ---- Print results 
cat("\nSummary statistics by sample and marker:\n")
print(expression_summary, n = Inf)

cat("\nNormality tests by sample:\n")
print(normality_by_sample, n = Inf)

cat("\nNormality tests by treatment:\n")
print(normality_by_treatment, n = Inf)

cat("\nResults saved under:\n", normality_results_dir, "\n")
cat("Plots saved under:\n", normality_plots_dir, "\n")




##use Log normalized values##




#Linear model fit for each BBB state----


#--- Day 5: High, Mid, and Low BBB-stability vessels
#--- Treatment effects: Kenward-Roger mixed models + Holm correction
#--- Holm correction across four markers within each BBB state
#--- Plot annotations display Holm-adjusted p-values and significance

library(dplyr)
library(tidyr)
library(purrr)
library(lme4)
library(lmerTest)
library(emmeans)
library(pbkrtest)
library(readr)
library(ggplot2)

#--- Output directories
output_dir <- file.path(
  "Arpan Multiplex",
  "Day5_Area_Normalized",
  "Linear Model"
)

results_dir <- file.path(output_dir, "Results")
plots_dir <- file.path(output_dir, "Plots")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

#--- Analysis settings
markers_to_test <- c("MFSD2A", "ABCB1", "GLUT1", "PLVAP")
bbb_states <- c("High", "Mid", "Low")

treatment_colors <- c(
  "Saline" = "#3C78A8",
  "Avastin" = "#D95F59"
)

#--- Check required columns
required_cols <- c(
  "Sample_ID",
  "Treatment",
  "BBB_State",
  "Marker",
  "Expression_per_area"
)

missing_cols <- setdiff(
  required_cols,
  names(day5_all_vessels_long)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in day5_all_vessels_long: ",
    paste(missing_cols, collapse = ", ")
  )
}

#--- Prepare expression values
day5_input <- day5_all_vessels_long %>%
  mutate(
    Expression_per_area = suppressWarnings(
      as.numeric(as.character(Expression_per_area))
    )
  )

#--- Apply QC filters when the corresponding columns are present
if ("Vessel_ID" %in% names(day5_input)) {
  day5_input <- day5_input %>%
    filter(as.character(Vessel_ID) != "V1052")
}

roi_col <- intersect(
  c(
    "vessel_centroid_inside_tumor_ROI",
    "Vessel_centroid_inside_tumor_ROI"
  ),
  names(day5_input)
)

if (length(roi_col) > 0) {
  day5_input <- day5_input %>%
    filter(
      as.character(.data[[roi_col[1]]]) == "1"
    )
}

pv_col <- intersect(
  c("PV_QC_status", "PV_QC_Status"),
  names(day5_input)
)

if (length(pv_col) > 0) {
  day5_input <- day5_input %>%
    filter(
      as.character(.data[[pv_col[1]]]) == "Analyzable"
    )
}

#--- Standardize treatment and BBB-state labels
day5_input <- day5_input %>%
  mutate(
    BBB_State_clean = tolower(
      trimws(as.character(BBB_State))
    ),
    Treatment_clean = trimws(as.character(Treatment)),
    Marker = trimws(as.character(Marker)),
    Sample_ID = trimws(as.character(Sample_ID))
  ) %>%
  mutate(
    BBB_State_standard = case_when(
      BBB_State_clean %in% c(
        "high", "high bbb", "high bbb stability",
        "high-bbb", "high-bbb stability"
      ) ~ "High",
      
      BBB_State_clean %in% c(
        "mid", "middle", "mid bbb", "mid bbb stability",
        "mid-bbb", "mid-bbb stability",
        "intermediate", "intermediate bbb"
      ) ~ "Mid",
      
      BBB_State_clean %in% c(
        "low", "low bbb", "low bbb stability",
        "low-bbb", "low-bbb stability"
      ) ~ "Low",
      
      TRUE ~ NA_character_
    ),
    Treatment = case_when(
      Treatment_clean %in% c(
        "Saline", "PBS", "-Avastin"
      ) ~ "Saline",
      
      Treatment_clean %in% c(
        "Avastin", "+Avastin"
      ) ~ "Avastin",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    BBB_State_standard %in% bbb_states,
    Marker %in% markers_to_test,
    Treatment %in% c("Saline", "Avastin"),
    !is.na(Sample_ID),
    Sample_ID != "",
    !is.na(Expression_per_area),
    is.finite(Expression_per_area),
    Expression_per_area >= 0
  ) %>%
  mutate(
    BBB_State_standard = factor(
      BBB_State_standard,
      levels = bbb_states
    ),
    Treatment = factor(
      Treatment,
      levels = c("Saline", "Avastin")
    ),
    Marker = factor(
      Marker,
      levels = markers_to_test
    ),
    Expression_log2 = log2(Expression_per_area + 1)
  ) %>%
  droplevels()

if (nrow(day5_input) == 0) {
  stop(
    "No data remain after filtering. Check BBB-state labels, ",
    "treatment labels, and QC columns."
  )
}

#--- Verify sample-to-treatment assignment
sample_treatment_check <- day5_input %>%
  distinct(Sample_ID, Treatment) %>%
  count(Sample_ID, name = "n_treatments")

if (any(sample_treatment_check$n_treatments > 1)) {
  stop(
    "At least one Sample_ID is assigned to multiple treatments. ",
    "Check Sample_ID uniqueness and treatment labels."
  )
}

#--- Export vessel counts and independent sample counts
sample_counts <- day5_input %>%
  distinct(
    BBB_State_standard, Marker, Treatment, Sample_ID
  ) %>%
  count(
    BBB_State_standard, Marker, Treatment,
    name = "n_independent_samples"
  )

vessel_counts <- day5_input %>%
  count(
    BBB_State_standard, Marker, Treatment,
    name = "n_vessels"
  )

sample_summary <- full_join(
  sample_counts,
  vessel_counts,
  by = c("BBB_State_standard", "Marker", "Treatment")
) %>%
  arrange(BBB_State_standard, Marker, Treatment)

write_csv(
  sample_summary,
  file.path(
    results_dir,
    "All_BBB_States_Sample_and_Vessel_Counts.csv"
  )
)

print(sample_summary, n = Inf)

#--- Fit one treatment-only mixed model per marker and BBB state
#--- Treatment contrast: Avastin minus Saline
#--- Kenward-Roger inference

fit_bbb_model <- function(state_name, marker_name) {
  
  dat <- day5_input %>%
    filter(
      as.character(BBB_State_standard) == state_name,
      as.character(Marker) == marker_name
    ) %>%
    droplevels()
  
  n_by_treatment <- dat %>%
    distinct(Sample_ID, Treatment) %>%
    count(Treatment)
  
  insufficient_samples <-
    n_distinct(dat$Treatment) != 2 ||
    any(n_by_treatment$n < 2) ||
    n_distinct(dat$Sample_ID) < 3
  
  if (insufficient_samples) {
    warning(
      "Insufficient independent samples for state ",
      state_name, ", marker ", marker_name
    )
    
    return(list(
      state = state_name,
      marker = marker_name,
      model = NULL,
      contrast = NULL,
      means = NULL,
      diagnostic = data.frame(
        BBB_State = state_name,
        Marker = marker_name,
        n_vessels = nrow(dat),
        n_samples = n_distinct(dat$Sample_ID),
        singular_fit = NA,
        model_status = "Not fitted: insufficient samples"
      )
    ))
  }
  
  model <- tryCatch(
    lmerTest::lmer(
      Expression_log2 ~ Treatment + (1 | Sample_ID),
      data = dat,
      REML = TRUE
    ),
    error = function(e) {
      warning(
        "Model failed for state ", state_name,
        ", marker ", marker_name, ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (is.null(model)) {
    return(list(
      state = state_name,
      marker = marker_name,
      model = NULL,
      contrast = NULL,
      means = NULL,
      diagnostic = data.frame(
        BBB_State = state_name,
        Marker = marker_name,
        n_vessels = nrow(dat),
        n_samples = n_distinct(dat$Sample_ID),
        singular_fit = NA,
        model_status = "Model fitting failed"
      )
    ))
  }
  
  #--- Estimated marginal means using Kenward-Roger inference
  emm <- tryCatch(
    emmeans::emmeans(
      model,
      ~ Treatment,
      lmer.df = "kenward-roger"
    ),
    error = function(e) {
      warning(
        "Kenward-Roger inference failed for state ",
        state_name, ", marker ", marker_name, ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (is.null(emm)) {
    return(list(
      state = state_name,
      marker = marker_name,
      model = model,
      contrast = NULL,
      means = NULL,
      diagnostic = data.frame(
        BBB_State = state_name,
        Marker = marker_name,
        n_vessels = nobs(model),
        n_samples = n_distinct(model.frame(model)$Sample_ID),
        singular_fit = lme4::isSingular(model, tol = 1e-4),
        model_status = "Kenward-Roger inference failed"
      )
    ))
  }
  
  #--- Contrast: Avastin minus Saline
  contrast_object <- emmeans::contrast(
    emm,
    method = list("Avastin - Saline" = c(-1, 1)),
    adjust = "none"
  )
  
  contrast_table <- as.data.frame(
    summary(
      contrast_object,
      infer = c(TRUE, TRUE),
      level = 0.95,
      adjust = "none"
    )
  ) %>%
    mutate(
      BBB_State = state_name,
      Marker = marker_name
    )
  
  #--- Standardize confidence-limit column names
  if (
    !"lower.CL" %in% names(contrast_table) &&
    "asymp.LCL" %in% names(contrast_table)
  ) {
    contrast_table$lower.CL <- contrast_table$asymp.LCL
  }
  
  if (
    !"upper.CL" %in% names(contrast_table) &&
    "asymp.UCL" %in% names(contrast_table)
  ) {
    contrast_table$upper.CL <- contrast_table$asymp.UCL
  }
  
  #--- Estimated marginal means
  means_table <- as.data.frame(
    summary(
      emm,
      infer = c(TRUE, TRUE)
    )
  ) %>%
    mutate(
      BBB_State = state_name,
      Marker = marker_name
    )
  
  #--- Model diagnostics
  diagnostic <- data.frame(
    BBB_State = state_name,
    Marker = marker_name,
    n_vessels = nobs(model),
    n_samples = n_distinct(model.frame(model)$Sample_ID),
    random_intercept_variance =
      as.numeric(VarCorr(model)$Sample_ID[1]),
    residual_variance = sigma(model)^2,
    singular_fit = lme4::isSingular(model, tol = 1e-4),
    model_status = "Fitted"
  )
  
  list(
    state = state_name,
    marker = marker_name,
    model = model,
    contrast = contrast_table,
    means = means_table,
    diagnostic = diagnostic
  )
}

#--- Run all 12 state-by-marker models
model_grid <- expand_grid(
  BBB_State = bbb_states,
  Marker = markers_to_test
)

all_model_results <- purrr::map2(
  model_grid$BBB_State,
  model_grid$Marker,
  fit_bbb_model
)

#--- Export model diagnostics
all_diagnostics <- bind_rows(
  lapply(all_model_results, function(x) x$diagnostic)
) %>%
  arrange(
    factor(BBB_State, levels = bbb_states),
    factor(Marker, levels = markers_to_test)
  )

write_csv(
  all_diagnostics,
  file.path(
    results_dir,
    "All_BBB_States_Model_Diagnostics.csv"
  )
)

print(all_diagnostics, n = Inf)

#--- Combine successful treatment contrasts
all_contrasts_raw <- bind_rows(
  lapply(all_model_results, function(x) x$contrast)
)

if (nrow(all_contrasts_raw) == 0) {
  stop("No treatment contrasts were successfully estimated.")
}

#--- Holm correction across four markers within each BBB state
all_contrasts <- all_contrasts_raw %>%
  group_by(BBB_State) %>%
  mutate(
    p_value_Holm = p.adjust(p.value, method = "holm")
  ) %>%
  ungroup() %>%
  mutate(
    Fold_change_transformed = 2^estimate,
    FC_CI_lower_transformed = 2^lower.CL,
    FC_CI_upper_transformed = 2^upper.CL,
    Significance = case_when(
      p_value_Holm < 0.001 ~ "***",
      p_value_Holm < 0.01 ~ "**",
      p_value_Holm < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    p_label = paste0(
      "Holm-adjusted p = ",
      format.pval(
        p_value_Holm,
        digits = 3,
        eps = 0.001
      )
    )
  ) %>%
  select(
    BBB_State, Marker, contrast, estimate, SE, df,
    lower.CL, upper.CL, p.value, p_value_Holm,
    Fold_change_transformed,
    FC_CI_lower_transformed,
    FC_CI_upper_transformed,
    Significance, p_label
  ) %>%
  arrange(
    factor(BBB_State, levels = bbb_states),
    factor(Marker, levels = markers_to_test)
  )

write_csv(
  all_contrasts,
  file.path(
    results_dir,
    "All_BBB_States_Treatment_Contrasts_KenwardRoger_Holm.csv"
  )
)

print(all_contrasts, n = Inf)

#--- Export estimated marginal means
all_means <- bind_rows(
  lapply(all_model_results, function(x) x$means)
) %>%
  arrange(BBB_State, Marker, Treatment)

write_csv(
  all_means,
  file.path(
    results_dir,
    "All_BBB_States_Estimated_Marginal_Means.csv"
  )
)

#--- Export state-specific result tables
for (state_name in bbb_states) {
  
  state_results <- all_contrasts %>%
    filter(BBB_State == state_name)
  
  write_csv(
    state_results,
    file.path(
      results_dir,
      paste0(
        "Day5_", state_name,
        "_BBB_Treatment_Contrasts_KenwardRoger_Holm.csv"
      )
    )
  )
}

#--- Prepare vessel-level plot data
#--- Individual vessels are descriptive; inference uses the mixed model
plot_data <- day5_input %>%
  mutate(
    BBB_State_standard = factor(
      as.character(BBB_State_standard),
      levels = bbb_states
    ),
    Treatment = factor(
      as.character(Treatment),
      levels = c("Saline", "Avastin")
    ),
    Marker = factor(
      as.character(Marker),
      levels = markers_to_test
    )
  )

#--- Calculate annotation positions for every state-marker panel
#--- Use the panel's observed range so labels adapt to different scales
panel_positions <- plot_data %>%
  group_by(BBB_State_standard, Marker) %>%
  summarise(
    y_min = min(Expression_log2, na.rm = TRUE),
    y_max = max(Expression_log2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_range = y_max - y_min,
    y_range = ifelse(
      !is.finite(y_range) | y_range == 0,
      pmax(abs(y_max), 1),
      y_range
    ),
    y_label = y_max + 0.12 * y_range,
    y_top = y_max + 0.40 * y_range
  )

#--- Prepare plot annotations from the model results
#--- Display Holm-adjusted p-values, not raw p-values
plot_stats <- all_contrasts %>%
  mutate(
    BBB_State_standard = factor(
      BBB_State,
      levels = bbb_states
    ),
    Marker = factor(
      Marker,
      levels = markers_to_test
    ),
    stat_label = paste0(
      "Holm-adjusted p = ",
      format.pval(
        p_value_Holm,
        digits = 3,
        eps = 0.001
      ),
      "\n",
      Significance
    )
  ) %>%
  select(
    BBB_State_standard,
    Marker,
    p_value_Holm,
    Significance,
    stat_label
  ) %>%
  left_join(
    panel_positions,
    by = c("BBB_State_standard", "Marker")
  ) %>%
  filter(
    !is.na(p_value_Holm),
    is.finite(y_label)
  )

#--- Shared plotting function
make_bbb_plot <- function(
    data,
    stats,
    combined = FALSE,
    state_name = NULL
) {
  
  p <- ggplot(
    data,
    aes(
      x = Treatment,
      y = Expression_log2,
      fill = Treatment
    )
  ) +
    geom_violin(
      trim = TRUE,
      alpha = 0.30,
      color = NA
    ) +
    geom_boxplot(
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.75
    ) +
    geom_point(
      aes(color = Treatment),
      position = position_jitter(
        width = 0.10,
        height = 0,
        seed = 123
      ),
      size = if (combined) 1.0 else 1.2,
      alpha = 0.35
    ) +
    geom_text(
      data = stats,
      aes(
        x = 1.5,
        y = y_label,
        label = stat_label
      ),
      inherit.aes = FALSE,
      size = 3.3,
      lineheight = 1.0,
      color = "black"
    ) +
    scale_fill_manual(
      values = treatment_colors,
      drop = FALSE
    ) +
    scale_color_manual(
      values = treatment_colors,
      drop = FALSE
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.35))
    ) +
    labs(
      title = if (combined) {
        "Day 5: Avastin Treatment Across BBB-Stability States"
      } else {
        paste0(
          "Day 5: ", state_name,
          " BBB-Stability Vessels"
        )
      },
      subtitle = if (combined) {
        paste(
          "Kenward-Roger mixed models;",
          "Holm correction across four markers within each BBB state"
        )
      } else {
        paste(
          "Avastin versus Saline;",
          "Kenward-Roger mixed models with Holm correction"
        )
      },
      x = NULL,
      y = "log2(Expression per area + 1)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      strip.text = element_text(
        face = "bold",
        size = if (combined) 12 else 13
      ),
      plot.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.margin = margin(10, 12, 10, 10)
    )
  
  if (combined) {
    p <- p +
      facet_grid(
        BBB_State_standard ~ Marker,
        scales = "free_y"
      )
  } else {
    p <- p +
      facet_wrap(
        ~ Marker,
        scales = "free_y",
        ncol = 2
      )
  }
  
  p
}

#--- Combined plot: High, Mid, and Low BBB states
p_all_states <- make_bbb_plot(
  data = plot_data,
  stats = plot_stats,
  combined = TRUE
)

ggsave(
  file.path(
    plots_dir,
    "Day5_All_BBB_States_Treatment_Expression.svg"
  ),
  p_all_states,
  width = 13,
  height = 10
)

ggsave(
  file.path(
    plots_dir,
    "Day5_All_BBB_States_Treatment_Expression.png"
  ),
  p_all_states,
  width = 13,
  height = 10,
  dpi = 300
)

print(p_all_states)

#--- Separate plots for each BBB state
for (state_name in bbb_states) {
  
  state_plot_data <- plot_data %>%
    filter(
      as.character(BBB_State_standard) == state_name
    )
  
  state_plot_stats <- plot_stats %>%
    filter(
      as.character(BBB_State_standard) == state_name
    )
  
  if (nrow(state_plot_data) == 0) next
  
  p_state <- make_bbb_plot(
    data = state_plot_data,
    stats = state_plot_stats,
    combined = FALSE,
    state_name = state_name
  )
  
  ggsave(
    file.path(
      plots_dir,
      paste0(
        "Day5_", state_name,
        "_BBB_Treatment_Expression.svg"
      )
    ),
    p_state,
    width = 9,
    height = 7
  )
  
  ggsave(
    file.path(
      plots_dir,
      paste0(
        "Day5_", state_name,
        "_BBB_Treatment_Expression.png"
      )
    ),
    p_state,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  print(p_state)
}

#--- Completion messages
message("All BBB-state Day 5 analyses completed.")
message("Models: Expression_log2 ~ Treatment + (1 | Sample_ID).")
message("Inference: Kenward-Roger treatment contrasts.")
message("Multiple testing: Holm correction across four markers within each BBB state.")
message("Plots display Holm-adjusted p-values and significance labels.")
message("Review sample counts, convergence, and singular-fit diagnostics.")





#--- Sample-level sensitivity analyses using Welch's t-test----


#--- Load packages

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr",
  "purrr",
  "tibble",
  "svglite"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(purrr)
library(tibble)

#--- Check input

if (!exists("analysis_data")) {
  stop(
    "Object 'analysis_data' was not found. ",
    "Run the Day 5 data preparation section first."
  )
}

required_columns <- c(
  "Sample_ID",
  "Treatment",
  "Tumor",
  "BBB_State",
  "Marker",
  "Expression_log2"
)

missing_columns <- setdiff(required_columns, names(analysis_data))

if (length(missing_columns) > 0) {
  stop(
    "Missing columns in analysis_data: ",
    paste(missing_columns, collapse = ", ")
  )
}

#--- Output directories

output_root <- file.path(
  "Arpan Multiplex",
  "Day5_Area_Normalized/Sensitivity"
)

results_dir <- file.path(output_root, "Results")
plots_dir <- file.path(output_root, "Plots")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

#--- Publication plot theme and colors

treatment_colors <- c(
  "Saline" = "#3C78A8",
  "Avastin" = "#D95F59"
)

publication_theme <- theme_classic(base_size = 17) +
  theme(
    plot.title = element_text(
      size = 21,
      face = "bold",
      hjust = 0.5,
      color = "#20252B"
    ),
    plot.subtitle = element_text(
      size = 14,
      hjust = 0.5,
      color = "#4D555E"
    ),
    axis.title = element_text(
      size = 17,
      face = "bold",
      color = "#20252B"
    ),
    axis.text = element_text(
      size = 15,
      color = "#20252B"
    ),
    strip.background = element_rect(
      fill = "#EEF1F4",
      color = NA
    ),
    strip.text = element_text(
      size = 17,
      face = "bold",
      color = "#20252B"
    ),
    legend.position = "top",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    panel.spacing = grid::unit(1.2, "lines"),
    plot.margin = margin(12, 18, 12, 12)
  )

marker_levels <- c("MFSD2A", "ABCB1", "GLUT1", "PLVAP")

bbb_states <- c("High BBB", "Mid BBB", "Low BBB")

#--- Prepare vessel-level input

sensitivity_input <- analysis_data %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    Treatment = as.character(Treatment),
    Tumor = as.character(Tumor),
    BBB_State = as.character(BBB_State),
    Marker = as.character(Marker),
    Expression_log2 = as.numeric(Expression_log2)
  ) %>%
  filter(
    Treatment %in% c("Saline", "Avastin"),
    Marker %in% marker_levels,
    is.finite(Expression_log2),
    !is.na(Sample_ID),
    !is.na(BBB_State)
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Saline", "Avastin")),
    Marker = factor(Marker, levels = marker_levels)
  )

if (nrow(sensitivity_input) == 0) {
  stop("No usable observations remain in sensitivity_input.")
}

#--- Verify one treatment and tumor label per biological sample

sample_treatment_check <- sensitivity_input %>%
  distinct(Sample_ID, Treatment) %>%
  count(Sample_ID, name = "n_treatments") %>%
  filter(n_treatments != 1)

if (nrow(sample_treatment_check) > 0) {
  stop(
    "At least one Sample_ID is assigned to multiple treatments. ",
    "Resolve sample labels before analysis."
  )
}

sample_tumor_check <- sensitivity_input %>%
  distinct(Sample_ID, Tumor) %>%
  count(Sample_ID, name = "n_tumor_labels") %>%
  filter(n_tumor_labels != 1)

if (nrow(sample_tumor_check) > 0) {
  stop(
    "At least one Sample_ID maps to multiple Tumor labels. ",
    "Check the sample identifiers before proceeding."
  )
}

#--- Define sample-level summaries

# Each biological sample contributes one mean per marker.
# The mean is calculated on the log2-transformed expression scale.
# The biological sample, not the individual vessel, is the unit of inference.

summarize_samples <- function(data, analysis_label) {
  data %>%
    group_by(Marker, Treatment, Sample_ID, Tumor) %>%
    summarise(
      Sample_mean_log2 = mean(Expression_log2, na.rm = TRUE),
      Sample_median_log2 = median(Expression_log2, na.rm = TRUE),
      n_vessels = sum(is.finite(Expression_log2)),
      .groups = "drop"
    ) %>%
    mutate(Analysis = analysis_label)
}

#--- Analysis 1: all vessels

sample_all_vessels <- summarize_samples(
  sensitivity_input,
  "All vessels"
)

#--- Analysis 2: High BBB vessels

sample_high_bbb <- sensitivity_input %>%
  filter(BBB_State == "High BBB") %>%
  summarize_samples("High BBB")

#--- Analysis 3: Mid BBB vessels

sample_mid_bbb <- sensitivity_input %>%
  filter(BBB_State == "Mid BBB") %>%
  summarize_samples("Mid BBB")

#--- Analysis 4: Low BBB vessels

sample_low_bbb <- sensitivity_input %>%
  filter(BBB_State == "Low BBB") %>%
  summarize_samples("Low BBB")

if (nrow(sample_high_bbb) == 0) {
  warning("No observations found for BBB_State == 'High BBB'.")
}

if (nrow(sample_mid_bbb) == 0) {
  warning("No observations found for BBB_State == 'Mid BBB'.")
}

if (nrow(sample_low_bbb) == 0) {
  warning("No observations found for BBB_State == 'Low BBB'.")
}

#--- Export sample-level data

write_csv(
  sample_all_vessels,
  file.path(results_dir, "Day5_all_vessels_sample_level_means.csv")
)

write_csv(
  sample_high_bbb,
  file.path(results_dir, "Day5_HighBBB_sample_level_means.csv")
)

write_csv(
  sample_mid_bbb,
  file.path(results_dir, "Day5_MidBBB_sample_level_means.csv")
)

write_csv(
  sample_low_bbb,
  file.path(results_dir, "Day5_LowBBB_sample_level_means.csv")
)

#--- Welch statistical testing function

run_welch_tests <- function(sample_data, analysis_label) {
  
  output <- lapply(marker_levels, function(marker_name) {
    
    dat <- sample_data %>%
      filter(as.character(Marker) == marker_name)
    
    saline <- dat %>%
      filter(Treatment == "Saline") %>%
      pull(Sample_mean_log2)
    
    avastin <- dat %>%
      filter(Treatment == "Avastin") %>%
      pull(Sample_mean_log2)
    
    saline <- saline[is.finite(saline)]
    avastin <- avastin[is.finite(avastin)]
    
    n_saline <- length(saline)
    n_avastin <- length(avastin)
    
    if (n_saline < 2 || n_avastin < 2) {
      return(tibble(
        Analysis = analysis_label,
        Marker = marker_name,
        n_Saline = n_saline,
        n_Avastin = n_avastin,
        Mean_Saline = if (n_saline > 0) mean(saline) else NA_real_,
        SD_Saline = if (n_saline > 1) sd(saline) else NA_real_,
        Mean_Avastin = if (n_avastin > 0) mean(avastin) else NA_real_,
        SD_Avastin = if (n_avastin > 1) sd(avastin) else NA_real_,
        Difference_Avastin_minus_Saline = NA_real_,
        Difference_CI_lower = NA_real_,
        Difference_CI_upper = NA_real_,
        Welch_p = NA_real_
      ))
    }
    
    welch <- t.test(
      avastin,
      saline,
      var.equal = FALSE,
      conf.level = 0.95
    )
    
    tibble(
      Analysis = analysis_label,
      Marker = marker_name,
      n_Saline = n_saline,
      n_Avastin = n_avastin,
      Mean_Saline = mean(saline),
      SD_Saline = sd(saline),
      Mean_Avastin = mean(avastin),
      SD_Avastin = sd(avastin),
      Difference_Avastin_minus_Saline = mean(avastin) - mean(saline),
      Difference_CI_lower = unname(welch$conf.int[1]),
      Difference_CI_upper = unname(welch$conf.int[2]),
      Welch_p = welch$p.value
    )
  })
  
  bind_rows(output) %>%
    mutate(
      # Holm correction across four markers within this analysis
      Welch_p_Holm = p.adjust(Welch_p, method = "holm"),
      Welch_significance = case_when(
        is.na(Welch_p_Holm) ~ NA_character_,
        Welch_p_Holm < 0.001 ~ "***",
        Welch_p_Holm < 0.01 ~ "**",
        Welch_p_Holm < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

#--- Run all four sensitivity analyses

stats_all_vessels <- run_welch_tests(
  sample_all_vessels,
  "All vessels"
)

stats_high_bbb <- run_welch_tests(
  sample_high_bbb,
  "High BBB"
)

stats_mid_bbb <- run_welch_tests(
  sample_mid_bbb,
  "Mid BBB"
)

stats_low_bbb <- run_welch_tests(
  sample_low_bbb,
  "Low BBB"
)

sensitivity_stats <- bind_rows(
  stats_all_vessels,
  stats_high_bbb,
  stats_mid_bbb,
  stats_low_bbb
) %>%
  mutate(
    Marker = factor(Marker, levels = marker_levels),
    Analysis = factor(
      Analysis,
      levels = c("All vessels", "High BBB", "Mid BBB", "Low BBB")
    )
  ) %>%
  arrange(Analysis, Marker)

#--- Export statistical results

write_csv(
  sensitivity_stats,
  file.path(results_dir, "Day5_sample_level_Welch_statistics_all_analyses.csv")
)

write_csv(
  stats_all_vessels,
  file.path(results_dir, "Day5_all_vessels_Welch_tests.csv")
)

write_csv(
  stats_high_bbb,
  file.path(results_dir, "Day5_HighBBB_Welch_tests.csv")
)

write_csv(
  stats_mid_bbb,
  file.path(results_dir, "Day5_MidBBB_Welch_tests.csv")
)

write_csv(
  stats_low_bbb,
  file.path(results_dir, "Day5_LowBBB_Welch_tests.csv")
)

#--- Prepare plot annotations

make_annotations <- function(sample_data, stats_data) {
  
  positions <- sample_data %>%
    group_by(Marker) %>%
    summarise(
      y_min = min(Sample_mean_log2, na.rm = TRUE),
      y_max = max(Sample_mean_log2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = y_max - y_min,
      y_range = if_else(
        is.finite(y_range) & y_range > 0,
        y_range,
        1
      ),
      y_position = y_max + 0.25 * y_range
    )
  
  stats_data %>%
    left_join(positions, by = "Marker") %>%
    mutate(
      label = paste0(
        "Welch Holm-adjusted p = ",
        ifelse(
          is.na(Welch_p_Holm),
          "NA",
          formatC(Welch_p_Holm, format = "f", digits = 3)
        )
      )
    )
}

#--- Plotting function

make_sample_boxplot <- function(
    sample_data,
    stats_data,
    analysis_label,
    file_prefix
) {
  
  plot_annotations <- make_annotations(sample_data, stats_data)
  
  # Individual marker plots
  for (marker_name in marker_levels) {
    
    dat <- sample_data %>%
      filter(as.character(Marker) == marker_name)
    
    ann <- plot_annotations %>%
      filter(as.character(Marker) == marker_name)
    
    if (nrow(dat) == 0) next
    
    p <- ggplot(
      dat,
      aes(
        x = Treatment,
        y = Sample_mean_log2,
        fill = Treatment
      )
    ) +
      geom_boxplot(
        width = 0.52,
        outlier.shape = NA,
        alpha = 0.50,
        color = "#292D32",
        linewidth = 0.8
      ) +
      geom_point(
        aes(color = Treatment),
        position = position_jitter(width = 0.055, height = 0),
        size = 4.2,
        alpha = 0.95
      ) +
      geom_text(
        data = ann,
        aes(
          x = 1.5,
          y = y_position,
          label = label
        ),
        inherit.aes = FALSE,
        size = 4.0,
        fontface = "bold",
        color = "#292D32"
      ) +
      scale_fill_manual(values = treatment_colors, drop = FALSE) +
      scale_color_manual(values = treatment_colors, drop = FALSE) +
      scale_x_discrete(drop = FALSE) +
      scale_y_continuous(
        expand = expansion(mult = c(0.08, 0.30))
      ) +
      labs(
        title = paste0(marker_name, " expression"),
        subtitle = paste0(
          "Day 5 | ",
          analysis_label,
          " | one point per biological sample"
        ),
        x = "Treatment",
        y = expression("Sample mean " * log[2]("Expression per area" + 1))
      ) +
      publication_theme +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 17, face = "bold")
      )
    
    # Export individual marker plot as SVG only
    ggsave(
      filename = file.path(
        plots_dir,
        paste0(file_prefix, "_", marker_name, ".svg")
      ),
      plot = p,
      width = 8.5,
      height = 7.5,
      device = svglite::svglite,
      bg = "white"
    )
  }
  
  # Combined four-marker figure
  combined_plot <- ggplot(
    sample_data,
    aes(
      x = Treatment,
      y = Sample_mean_log2,
      fill = Treatment
    )
  ) +
    geom_boxplot(
      width = 0.52,
      outlier.shape = NA,
      alpha = 0.50,
      color = "#292D32",
      linewidth = 0.65
    ) +
    geom_point(
      aes(color = Treatment),
      position = position_jitter(width = 0.055, height = 0),
      size = 3.5,
      alpha = 0.95
    ) +
    geom_text(
      data = plot_annotations,
      aes(
        x = 1.5,
        y = y_position,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3.5,
      fontface = "bold",
      color = "#292D32"
    ) +
    facet_wrap(~ Marker, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = treatment_colors, drop = FALSE) +
    scale_color_manual(values = treatment_colors, drop = FALSE) +
    scale_y_continuous(
      expand = expansion(mult = c(0.08, 0.30))
    ) +
    labs(
      title = "Day 5 Avastin treatment: sample-level sensitivity analysis",
      subtitle = paste0(
        analysis_label,
        " | each point represents one independent tumor sample"
      ),
      x = "Treatment",
      y = expression("Sample mean " * log[2]("Expression per area" + 1))
    ) +
    publication_theme +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 17, face = "bold")
    )
  
  # Export combined four-marker plot as SVG only
  ggsave(
    filename = file.path(plots_dir, paste0(file_prefix, "_combined.svg")),
    plot = combined_plot,
    width = 13,
    height = 10,
    device = svglite::svglite,
    bg = "white"
  )
  
  invisible(combined_plot)
}

#--- Generate plots: all vessels

plot_all_vessels <- make_sample_boxplot(
  sample_data = sample_all_vessels,
  stats_data = stats_all_vessels,
  analysis_label = "All vessels",
  file_prefix = "Day5_AllVessels_SampleLevel"
)

#--- Generate plots: High BBB

plot_high_bbb <- make_sample_boxplot(
  sample_data = sample_high_bbb,
  stats_data = stats_high_bbb,
  analysis_label = "High BBB vessels",
  file_prefix = "Day5_HighBBB_SampleLevel"
)

#--- Generate plots: Mid BBB

plot_mid_bbb <- make_sample_boxplot(
  sample_data = sample_mid_bbb,
  stats_data = stats_mid_bbb,
  analysis_label = "Mid BBB vessels",
  file_prefix = "Day5_MidBBB_SampleLevel"
)

#--- Generate plots: Low BBB

plot_low_bbb <- make_sample_boxplot(
  sample_data = sample_low_bbb,
  stats_data = stats_low_bbb,
  analysis_label = "Low BBB vessels",
  file_prefix = "Day5_LowBBB_SampleLevel"
)

#--- Export treatment-direction summary

direction_summary <- sensitivity_stats %>%
  mutate(
    Direction = case_when(
      is.na(Difference_Avastin_minus_Saline) ~ "Not estimable",
      Difference_Avastin_minus_Saline > 0 ~ "Higher with Avastin",
      Difference_Avastin_minus_Saline < 0 ~ "Lower with Avastin",
      TRUE ~ "No difference"
    )
  ) %>%
  select(
    Analysis,
    Marker,
    n_Saline,
    n_Avastin,
    Mean_Saline,
    SD_Saline,
    Mean_Avastin,
    SD_Avastin,
    Difference_Avastin_minus_Saline,
    Difference_CI_lower,
    Difference_CI_upper,
    Direction,
    Welch_p,
    Welch_p_Holm,
    Welch_significance
  )

write_csv(
  direction_summary,
  file.path(results_dir, "Day5_Welch_effect_direction_summary.csv")
)

#--- Print results

cat("\n--- All-vessel sample-level Welch statistics ---\n")
print(stats_all_vessels, n = Inf)

cat("\n--- High-BBB sample-level Welch statistics ---\n")
print(stats_high_bbb, n = Inf)

cat("\n--- Mid-BBB sample-level Welch statistics ---\n")
print(stats_mid_bbb, n = Inf)

cat("\n--- Low-BBB sample-level Welch statistics ---\n")
print(stats_low_bbb, n = Inf)

cat("\n--- Output locations ---\n")
cat("Results: ", normalizePath(results_dir), "\n")
cat("Plots: ", normalizePath(plots_dir), "\n")

cat("\nAll four Welch sensitivity analyses completed.\n")




#Area calibrated----
# ---- Day 5 Avastin vs Saline: Log2 Area Analysis
# Primary: all vessels pooled
# Secondary: treatment comparisons within BBB states
# Interaction: Treatment x BBB state
# Main inference: mixed models with sample-level clustering
# Sensitivity: mouse-level Welch tests
# Plots:
#   1. All-vessel violin plot with statistics
#   2. Three-panel BBB-state violin plot with statistics
#   3. All-vessel mouse-level boxplot with statistics
#   4. Three-panel BBB-state mouse-level boxplot with statistics
# Output: CSV tables and PNG/SVG plots
#
# Note: no print() calls are used; results are saved to CSV files.

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(readr)
library(tibble)

# ---- Settings

area_data <- day5_all_vessels_long

treatment_colors <- c(
  "Saline" = "#0072B2",
  "Avastin" = "#D55E00"
)

output_dir <- file.path(
  "Arpan Multiplex",
  "Day5_Area_Log2_Comprehensive_Analysis"
)

results_dir <- file.path(output_dir, "Results")
plots_dir <- file.path(output_dir, "Plots")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Configure emmeans degrees-of-freedom limits

emmeans_options <- emmeans::emm_options()

emmeans::emm_options(
  pbkrtest.limit = max(
    nrow(area_data),
    emmeans_options$pbkrtest.limit
  ),
  lmerTest.limit = max(
    nrow(area_data),
    emmeans_options$lmerTest.limit
  )
)

# ---- Check required columns

required_cols <- c(
  "Area_calibrated",
  "Treatment",
  "Sample_ID"
)

missing_cols <- setdiff(required_cols, names(area_data))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", "),
    "\nAvailable columns:\n",
    paste(names(area_data), collapse = ", ")
  )
}

bbb_candidates <- c(
  "BBB_group",
  "BBB_State",
  "BBB_state",
  "BBB_stability",
  "BBB_Stability",
  "BBB_category",
  "BBB_Category"
)

bbb_matches <- bbb_candidates[bbb_candidates %in% names(area_data)]

if (length(bbb_matches) == 0) {
  stop(
    "Could not identify the BBB-state column. Available columns:\n",
    paste(names(area_data), collapse = ", ")
  )
}

bbb_col <- bbb_matches[1]

message("Using BBB-state column: ", bbb_col)

# ---- Standardize labels and filter eligible vessels

area_data <- area_data %>%
  mutate(
    BBB_state = trimws(as.character(.data[[bbb_col]])),
    Treatment = case_when(
      grepl("saline|pbs|vehicle", Treatment, ignore.case = TRUE) ~ "Saline",
      grepl("avastin", Treatment, ignore.case = TRUE) ~ "Avastin",
      TRUE ~ as.character(Treatment)
    ),
    Area_calibrated = as.numeric(Area_calibrated),
    Sample_ID = as.character(Sample_ID)
  ) %>%
  filter(
    Treatment %in% c("Saline", "Avastin"),
    !is.na(Area_calibrated),
    is.finite(Area_calibrated),
    !is.na(Sample_ID),
    Sample_ID != ""
  ) %>%
  mutate(
    BBB_state = case_when(
      grepl("^High", BBB_state, ignore.case = TRUE) ~ "High",
      grepl("^Mid", BBB_state, ignore.case = TRUE) ~ "Mid",
      grepl("^Low", BBB_state, ignore.case = TRUE) ~ "Low",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(BBB_state)) %>%
  mutate(
    BBB_state = factor(
      BBB_state,
      levels = c("High", "Mid", "Low")
    ),
    Treatment = factor(
      Treatment,
      levels = c("Saline", "Avastin")
    )
  ) %>%
  droplevels()

if (nrow(area_data) == 0) {
  stop("No eligible vessels remain after filtering.")
}

# ---- Verify sample-level experimental design

sample_treatment_check <- area_data %>%
  distinct(Sample_ID, Treatment) %>%
  count(Sample_ID, name = "n_treatments") %>%
  filter(n_treatments != 1)

if (nrow(sample_treatment_check) > 0) {
  stop("A Sample_ID is assigned to multiple treatments.")
}

sample_summary <- area_data %>%
  distinct(Sample_ID, Treatment) %>%
  count(Treatment, name = "n_samples")

sample_state_summary <- area_data %>%
  distinct(Sample_ID, BBB_state, Treatment) %>%
  count(BBB_state, Treatment, name = "n_samples")

write_csv(
  sample_summary,
  file.path(results_dir, "Day5_Area_Sample_Summary.csv")
)

write_csv(
  sample_state_summary,
  file.path(results_dir, "Day5_Area_Sample_Summary_By_BBB.csv")
)

# ---- Inspect original area values

area_QC <- area_data %>%
  summarise(
    n_vessels = n(),
    n_samples = n_distinct(Sample_ID),
    n_zero = sum(Area_calibrated == 0),
    n_negative = sum(Area_calibrated < 0),
    minimum = min(Area_calibrated),
    median = median(Area_calibrated),
    maximum = max(Area_calibrated)
  )

write_csv(
  area_QC,
  file.path(results_dir, "Day5_Area_Log2_QC.csv")
)

if (any(area_data$Area_calibrated <= 0)) {
  stop(
    "Area_calibrated contains zero or negative values. ",
    "Review Day5_Area_Log2_QC.csv before choosing a transformation. ",
    "No arbitrary offset has been applied."
  )
}

area_data <- area_data %>%
  mutate(
    Area_log2 = log2(Area_calibrated)
  )

# ---- Descriptive vessel-level summaries

vessel_summary <- area_data %>%
  group_by(BBB_state, Treatment) %>%
  summarise(
    n_vessels = n(),
    n_samples = n_distinct(Sample_ID),
    mean_area = mean(Area_calibrated),
    median_area = median(Area_calibrated),
    sd_area = sd(Area_calibrated),
    mean_log2_area = mean(Area_log2),
    median_log2_area = median(Area_log2),
    .groups = "drop"
  )

write_csv(
  vessel_summary,
  file.path(results_dir, "Day5_Area_Vessel_Summary.csv")
)

# ---- Helpers

get_ci_columns <- function(result_df) {
  
  lower_candidates <- c("lower.CL", "asymp.LCL", "lower.HPD")
  upper_candidates <- c("upper.CL", "asymp.UCL", "upper.HPD")
  
  lower_match <- intersect(lower_candidates, names(result_df))
  upper_match <- intersect(upper_candidates, names(result_df))
  
  if (length(lower_match) == 0 || length(upper_match) == 0) {
    stop(
      "Could not identify confidence interval columns. Returned columns: ",
      paste(names(result_df), collapse = ", ")
    )
  }
  
  list(
    lower = lower_match[1],
    upper = upper_match[1]
  )
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, "< 0.001", paste0("= ", sprintf("%.3f", p)))
  )
}

format_stats_label <- function(p, p_holm, fold_change = NULL) {
  
  labels <- paste0(
    "p ", format_p(p),
    "\nHolm-adjusted p ", format_p(p_holm)
  )
  
  if (!is.null(fold_change)) {
    labels <- paste0(
      labels,
      "\nFold change = ",
      ifelse(
        is.na(fold_change),
        "NA",
        sprintf("%.2f", fold_change)
      )
    )
  }
  
  labels
}

extract_treatment_contrast <- function(model, dat, analysis_name) {
  
  emm <- emmeans::emmeans(
    model,
    ~ Treatment,
    lmer.df = "kenward-roger",
    pbkrtest.limit = nrow(dat),
    lmerTest.limit = nrow(dat)
  )
  
  contrast_df <- as.data.frame(
    summary(
      emmeans::contrast(
        emm,
        method = list("Avastin - Saline" = c(-1, 1)),
        adjust = "none"
      ),
      infer = c(TRUE, TRUE),
      level = 0.95
    )
  )
  
  ci_cols <- get_ci_columns(contrast_df)
  
  statistic_candidates <- c("t.ratio", "z.ratio")
  statistic_match <- intersect(
    statistic_candidates,
    names(contrast_df)
  )
  
  statistic_col <- if (length(statistic_match) > 0) {
    statistic_match[1]
  } else {
    NA_character_
  }
  
  tibble(
    Analysis = analysis_name,
    log2_difference = contrast_df$estimate[1],
    SE = contrast_df$SE[1],
    df = contrast_df$df[1],
    test_statistic = if (!is.na(statistic_col)) {
      contrast_df[[statistic_col]][1]
    } else {
      NA_real_
    },
    statistic_type = statistic_col,
    fold_change = 2^contrast_df$estimate[1],
    FC_CI_lower = 2^contrast_df[[ci_cols$lower]][1],
    FC_CI_upper = 2^contrast_df[[ci_cols$upper]][1],
    p_value = contrast_df$p.value[1],
    inference_method = if ("z.ratio" %in% names(contrast_df)) {
      "Asymptotic"
    } else {
      "t-based"
    }
  )
}

# ---- 1. PRIMARY ANALYSIS: ALL VESSELS POOLED

n_samples_by_treatment <- area_data %>%
  distinct(Sample_ID, Treatment) %>%
  count(Treatment, name = "n_samples")

if (
  !all(c("Saline", "Avastin") %in% n_samples_by_treatment$Treatment) ||
  any(n_samples_by_treatment$n_samples < 2)
) {
  stop("At least two independent samples per treatment are required.")
}

model_pooled <- lmerTest::lmer(
  Area_log2 ~ Treatment + (1 | Sample_ID),
  data = area_data,
  REML = TRUE
)

pooled_result <- extract_treatment_contrast(
  model_pooled,
  area_data,
  "All_vessels"
) %>%
  mutate(
    n_Saline_samples = n_distinct(
      area_data$Sample_ID[area_data$Treatment == "Saline"]
    ),
    n_Avastin_samples = n_distinct(
      area_data$Sample_ID[area_data$Treatment == "Avastin"]
    ),
    n_Saline_vessels = sum(area_data$Treatment == "Saline"),
    n_Avastin_vessels = sum(area_data$Treatment == "Avastin")
  ) %>%
  select(
    Analysis,
    n_Saline_samples,
    n_Avastin_samples,
    n_Saline_vessels,
    n_Avastin_vessels,
    everything()
  )

# ---- 2. SECONDARY ANALYSIS: EACH BBB STATE SEPARATELY

fit_bbb_model <- function(state_name) {
  
  dat <- area_data %>%
    filter(BBB_state == state_name) %>%
    droplevels()
  
  n_saline <- n_distinct(dat$Sample_ID[dat$Treatment == "Saline"])
  n_avastin <- n_distinct(dat$Sample_ID[dat$Treatment == "Avastin"])
  
  n_saline_vessels <- sum(dat$Treatment == "Saline")
  n_avastin_vessels <- sum(dat$Treatment == "Avastin")
  
  if (
    n_saline < 2 ||
    n_avastin < 2 ||
    n_distinct(dat$Treatment) < 2
  ) {
    return(
      tibble(
        Analysis = state_name,
        n_Saline_samples = n_saline,
        n_Avastin_samples = n_avastin,
        n_Saline_vessels = n_saline_vessels,
        n_Avastin_vessels = n_avastin_vessels,
        log2_difference = NA_real_,
        SE = NA_real_,
        df = NA_real_,
        test_statistic = NA_real_,
        statistic_type = NA_character_,
        fold_change = NA_real_,
        FC_CI_lower = NA_real_,
        FC_CI_upper = NA_real_,
        p_value = NA_real_,
        inference_method = NA_character_,
        model_status = "Insufficient sample-level replication"
      )
    )
  }
  
  model <- tryCatch(
    lmerTest::lmer(
      Area_log2 ~ Treatment + (1 | Sample_ID),
      data = dat,
      REML = TRUE
    ),
    error = function(e) e
  )
  
  if (inherits(model, "error")) {
    return(
      tibble(
        Analysis = state_name,
        n_Saline_samples = n_saline,
        n_Avastin_samples = n_avastin,
        n_Saline_vessels = n_saline_vessels,
        n_Avastin_vessels = n_avastin_vessels,
        log2_difference = NA_real_,
        SE = NA_real_,
        df = NA_real_,
        test_statistic = NA_real_,
        statistic_type = NA_character_,
        fold_change = NA_real_,
        FC_CI_lower = NA_real_,
        FC_CI_upper = NA_real_,
        p_value = NA_real_,
        inference_method = NA_character_,
        model_status = conditionMessage(model)
      )
    )
  }
  
  extract_treatment_contrast(
    model,
    dat,
    state_name
  ) %>%
    mutate(
      n_Saline_samples = n_saline,
      n_Avastin_samples = n_avastin,
      n_Saline_vessels = n_saline_vessels,
      n_Avastin_vessels = n_avastin_vessels,
      model_status = "OK"
    ) %>%
    select(
      Analysis,
      n_Saline_samples,
      n_Avastin_samples,
      n_Saline_vessels,
      n_Avastin_vessels,
      everything()
    )
}

bbb_results <- bind_rows(
  lapply(c("High", "Mid", "Low"), fit_bbb_model)
) %>%
  mutate(
    p_value_Holm_BBB = p.adjust(p_value, method = "holm")
  )

# ---- Combine pooled and BBB-state treatment comparisons

treatment_results <- bind_rows(
  pooled_result %>% mutate(model_status = "OK"),
  bbb_results
) %>%
  mutate(
    p_value_Holm_four_tests = p.adjust(p_value, method = "holm")
  )

write_csv(
  treatment_results,
  file.path(results_dir, "Day5_Area_Treatment_Comparisons.csv")
)

# ---- 3. BBB-ADJUSTED MODEL AND TREATMENT x BBB INTERACTION

model_additive <- lmerTest::lmer(
  Area_log2 ~ Treatment + BBB_state + (1 | Sample_ID),
  data = area_data,
  REML = FALSE
)

model_interaction <- lmerTest::lmer(
  Area_log2 ~ Treatment * BBB_state + (1 | Sample_ID),
  data = area_data,
  REML = FALSE
)

# ---- Global interaction likelihood-ratio test

interaction_test <- anova(
  model_additive,
  model_interaction,
  refit = FALSE
)

interaction_test_df <- as.data.frame(interaction_test) %>%
  rownames_to_column("Model")

write_csv(
  interaction_test_df,
  file.path(results_dir, "Day5_Area_BBB_Interaction_LRT.csv")
)

# ---- Type III fixed-effect tests

contrasts(area_data$Treatment) <- contr.sum(2)
contrasts(area_data$BBB_state) <- contr.sum(3)

model_interaction_KR <- lmerTest::lmer(
  Area_log2 ~ Treatment * BBB_state + (1 | Sample_ID),
  data = area_data,
  REML = TRUE
)

interaction_type3 <- as.data.frame(
  anova(
    model_interaction_KR,
    type = 3,
    ddf = "Kenward-Roger"
  )
) %>%
  rownames_to_column("Effect")

write_csv(
  interaction_type3,
  file.path(
    results_dir,
    "Day5_Area_BBB_Interaction_Kenward_Roger.csv"
  )
)

# ---- Treatment contrasts within each BBB state from interaction model

interaction_emm <- emmeans::emmeans(
  model_interaction_KR,
  ~ Treatment | BBB_state,
  lmer.df = "kenward-roger",
  pbkrtest.limit = nrow(area_data),
  lmerTest.limit = nrow(area_data)
)

interaction_contrasts_raw <- as.data.frame(
  summary(
    emmeans::contrast(
      interaction_emm,
      method = list("Avastin - Saline" = c(-1, 1)),
      by = "BBB_state",
      adjust = "none"
    ),
    infer = c(TRUE, TRUE),
    level = 0.95
  )
)

interaction_ci_cols <- get_ci_columns(interaction_contrasts_raw)

statistic_candidates <- c("t.ratio", "z.ratio")
statistic_match <- intersect(
  statistic_candidates,
  names(interaction_contrasts_raw)
)

if (length(statistic_match) == 0) {
  stop(
    "Could not identify t.ratio or z.ratio. Returned columns: ",
    paste(names(interaction_contrasts_raw), collapse = ", ")
  )
}

statistic_col <- statistic_match[1]

interaction_contrasts <- interaction_contrasts_raw %>%
  mutate(
    fold_change = 2^estimate,
    FC_CI_lower = 2^.data[[interaction_ci_cols$lower]],
    FC_CI_upper = 2^.data[[interaction_ci_cols$upper]],
    p_value_Holm = p.adjust(p.value, method = "holm"),
    significant_Holm = !is.na(p_value_Holm) & p_value_Holm < 0.05,
    test_statistic = .data[[statistic_col]],
    statistic_type = statistic_col,
    inference_method = if ("z.ratio" %in% names(interaction_contrasts_raw)) {
      "Asymptotic"
    } else {
      "t-based"
    }
  ) %>%
  rename(
    log2_difference = estimate,
    p_value = p.value
  ) %>%
  select(
    BBB_state,
    contrast,
    log2_difference,
    SE,
    df,
    statistic_type,
    test_statistic,
    p_value,
    p_value_Holm,
    significant_Holm,
    fold_change,
    FC_CI_lower,
    FC_CI_upper,
    inference_method,
    everything()
  )

write_csv(
  interaction_contrasts,
  file.path(
    results_dir,
    "Day5_Area_Interaction_Model_Contrasts.csv"
  )
)

# ---- 4. MOUSE-LEVEL SENSITIVITY ANALYSIS

# One mean per mouse/sample across all vessels

mouse_level_pooled <- area_data %>%
  group_by(Sample_ID, Treatment) %>%
  summarise(
    n_vessels = n(),
    mean_log2_area = mean(Area_log2),
    median_log2_area = median(Area_log2),
    .groups = "drop"
  ) %>%
  mutate(BBB_state = "All_vessels")

# One mean per mouse/sample within each BBB state

mouse_level_bbb <- area_data %>%
  group_by(Sample_ID, Treatment, BBB_state) %>%
  summarise(
    n_vessels = n(),
    mean_log2_area = mean(Area_log2),
    median_log2_area = median(Area_log2),
    .groups = "drop"
  ) %>%
  mutate(BBB_state = as.character(BBB_state))

mouse_level_all <- bind_rows(
  mouse_level_pooled,
  mouse_level_bbb
) %>%
  mutate(
    BBB_state = factor(
      BBB_state,
      levels = c("All_vessels", "High", "Mid", "Low")
    )
  )

welch_mouse_test <- function(dat, analysis_name) {
  
  saline <- dat %>%
    filter(Treatment == "Saline") %>%
    pull(mean_log2_area)
  
  avastin <- dat %>%
    filter(Treatment == "Avastin") %>%
    pull(mean_log2_area)
  
  if (length(saline) < 2 || length(avastin) < 2) {
    return(
      tibble(
        Analysis = analysis_name,
        n_Saline_mice = length(saline),
        n_Avastin_mice = length(avastin),
        log2_difference = NA_real_,
        fold_change = NA_real_,
        FC_CI_lower = NA_real_,
        FC_CI_upper = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  test <- t.test(avastin, saline, var.equal = FALSE)
  difference <- mean(avastin) - mean(saline)
  
  tibble(
    Analysis = analysis_name,
    n_Saline_mice = length(saline),
    n_Avastin_mice = length(avastin),
    log2_difference = difference,
    fold_change = 2^difference,
    FC_CI_lower = 2^unname(test$conf.int[1]),
    FC_CI_upper = 2^unname(test$conf.int[2]),
    p_value = test$p.value
  )
}

mouse_sensitivity_results <- bind_rows(
  lapply(
    c("All_vessels", "High", "Mid", "Low"),
    function(nm) {
      dat <- mouse_level_all %>%
        filter(BBB_state == nm)
      
      welch_mouse_test(dat, nm)
    }
  )
) %>%
  mutate(
    p_value_Holm = p.adjust(p_value, method = "holm")
  )

write_csv(
  mouse_level_all,
  file.path(results_dir, "Day5_Area_Mouse_Level_Data.csv")
)

write_csv(
  mouse_sensitivity_results,
  file.path(
    results_dir,
    "Day5_Area_Mouse_Level_Welch_Sensitivity.csv"
  )
)

# ---- 5. Prepare plot annotation tables

# Main violin plots use mixed-model treatment comparisons.
# Holm correction is across the pooled comparison and three BBB-state tests.

pooled_plot_stats <- treatment_results %>%
  filter(Analysis == "All_vessels") %>%
  transmute(
    label = format_stats_label(
      p_value,
      p_value_Holm_four_tests,
      fold_change
    )
  )

bbb_plot_stats <- treatment_results %>%
  filter(Analysis %in% c("High", "Mid", "Low")) %>%
  mutate(
    BBB_state = factor(
      Analysis,
      levels = c("High", "Mid", "Low")
    ),
    label = format_stats_label(
      p_value,
      p_value_Holm_four_tests,
      fold_change
    )
  ) %>%
  select(BBB_state, label)

# Sensitivity plots use Welch tests on mouse-level means.

mouse_pooled_plot_stats <- mouse_sensitivity_results %>%
  filter(Analysis == "All_vessels") %>%
  transmute(
    label = format_stats_label(
      p_value,
      p_value_Holm
    )
  )

mouse_bbb_plot_stats <- mouse_sensitivity_results %>%
  filter(Analysis %in% c("High", "Mid", "Low")) %>%
  mutate(
    BBB_state = factor(
      Analysis,
      levels = c("High", "Mid", "Low")
    ),
    label = format_stats_label(
      p_value,
      p_value_Holm
    )
  ) %>%
  select(BBB_state, label)

# ---- 6. Plot 1: standalone all-vessel violin plot

all_y_range <- range(area_data$Area_log2, na.rm = TRUE)
all_y_padding <- max(diff(all_y_range) * 0.16, 0.15)

pooled_plot_stats <- pooled_plot_stats %>%
  mutate(
    Treatment = factor("Saline", levels = c("Saline", "Avastin")),
    y_position = all_y_range[2] + all_y_padding
  )

p_all_vessels <- ggplot(
  area_data,
  aes(x = Treatment, y = Area_log2, fill = Treatment)
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.35,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  geom_jitter(
    aes(color = Treatment),
    width = 0.12,
    alpha = 0.22,
    size = 1.0
  ) +
  geom_text(
    data = pooled_plot_stats,
    aes(
      x = Treatment,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1,
    size = 4.0,
    lineheight = 1.05
  ) +
  scale_fill_manual(values = treatment_colors) +
  scale_color_manual(values = treatment_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0.04, 0.23))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Day 5 vessel area: all vessels",
    subtitle = paste(
      "Mixed-effects model; sample included as a random intercept",
      "| Log2-transformed area"
    ),
    x = NULL,
    y = "log2(Area_calibrated)",
    fill = "Treatment",
    color = "Treatment"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    plot.margin = margin(12, 55, 12, 12)
  )

# ---- 7. Plot 2: three-panel BBB-state violin plot

bbb_y_positions <- area_data %>%
  group_by(BBB_state) %>%
  summarise(
    y_max = max(Area_log2, na.rm = TRUE),
    y_range = diff(range(Area_log2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    y_position = y_max + pmax(y_range * 0.18, 0.15)
  )

bbb_plot_stats <- bbb_plot_stats %>%
  left_join(
    bbb_y_positions %>% select(BBB_state, y_position),
    by = "BBB_state"
  ) %>%
  mutate(
    Treatment = factor("Saline", levels = c("Saline", "Avastin"))
  )

p_bbb_vessels <- ggplot(
  area_data,
  aes(x = Treatment, y = Area_log2, fill = Treatment)
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.35,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  geom_jitter(
    aes(color = Treatment),
    width = 0.12,
    alpha = 0.22,
    size = 0.9
  ) +
  geom_text(
    data = bbb_plot_stats,
    aes(
      x = Treatment,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1,
    size = 3.5,
    lineheight = 1.0
  ) +
  facet_wrap(
    ~ BBB_state,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(values = treatment_colors) +
  scale_color_manual(values = treatment_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0.04, 0.25))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Day 5 vessel area by BBB stability",
    subtitle = paste(
      "Mixed-effects treatment comparisons;",
      "Holm correction across four comparisons"
    ),
    x = NULL,
    y = "log2(Area_calibrated)",
    fill = "Treatment",
    color = "Treatment"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    plot.margin = margin(12, 35, 12, 12)
  )

# ---- 8. Plot 3: standalone all-vessel mouse-level boxplot

mouse_pooled_data <- mouse_level_all %>%
  filter(BBB_state == "All_vessels") %>%
  droplevels()

mouse_pooled_range <- range(
  mouse_pooled_data$mean_log2_area,
  na.rm = TRUE
)

mouse_pooled_padding <- max(
  diff(mouse_pooled_range) * 0.22,
  0.15
)

mouse_pooled_plot_stats <- mouse_pooled_plot_stats %>%
  mutate(
    Treatment = factor("Saline", levels = c("Saline", "Avastin")),
    y_position = mouse_pooled_range[2] + mouse_pooled_padding
  )

p_mouse_all <- ggplot(
  mouse_pooled_data,
  aes(x = Treatment, y = mean_log2_area, fill = Treatment)
) +
  geom_boxplot(
    width = 0.48,
    alpha = 0.35,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = Treatment),
    width = 0.08,
    height = 0,
    size = 3.5,
    alpha = 0.9
  ) +
  geom_text(
    data = mouse_pooled_plot_stats,
    aes(
      x = Treatment,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1,
    size = 4.0,
    lineheight = 1.05
  ) +
  scale_fill_manual(values = treatment_colors) +
  scale_color_manual(values = treatment_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0.06, 0.25))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Day 5 vessel area: mouse-level sensitivity analysis",
    subtitle = paste(
      "Each point is one mouse's mean log2 vessel area;",
      "Welch's two-sample t-test"
    ),
    x = NULL,
    y = "Mean log2(Area_calibrated) per mouse",
    fill = "Treatment",
    color = "Treatment"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    plot.margin = margin(12, 55, 12, 12)
  )

# ---- 9. Plot 4: three-panel BBB-state mouse-level boxplot

mouse_bbb_data <- mouse_level_all %>%
  filter(BBB_state %in% c("High", "Mid", "Low")) %>%
  droplevels()

mouse_bbb_y_positions <- mouse_bbb_data %>%
  group_by(BBB_state) %>%
  summarise(
    y_max = max(mean_log2_area, na.rm = TRUE),
    y_range = diff(range(mean_log2_area, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    y_position = y_max + pmax(y_range * 0.25, 0.15)
  )

mouse_bbb_plot_stats <- mouse_bbb_plot_stats %>%
  left_join(
    mouse_bbb_y_positions %>% select(BBB_state, y_position),
    by = "BBB_state"
  ) %>%
  mutate(
    Treatment = factor("Saline", levels = c("Saline", "Avastin"))
  )

p_mouse_bbb <- ggplot(
  mouse_bbb_data,
  aes(x = Treatment, y = mean_log2_area, fill = Treatment)
) +
  geom_boxplot(
    width = 0.48,
    alpha = 0.35,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = Treatment),
    width = 0.08,
    height = 0,
    size = 2.8,
    alpha = 0.9
  ) +
  geom_text(
    data = mouse_bbb_plot_stats,
    aes(
      x = Treatment,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1,
    size = 3.4,
    lineheight = 1.0
  ) +
  facet_wrap(
    ~ BBB_state,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(values = treatment_colors) +
  scale_color_manual(values = treatment_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0.06, 0.27))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Day 5 vessel area: mouse-level sensitivity by BBB state",
    subtitle = paste(
      "Each point is one mouse's mean log2 vessel area;",
      "Welch tests with Holm correction across four comparisons"
    ),
    x = NULL,
    y = "Mean log2(Area_calibrated) per mouse",
    fill = "Treatment",
    color = "Treatment"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom",
    plot.margin = margin(12, 35, 12, 12)
  )

# ---- 10. Save plots as PNG and SVG

ggsave(
  file.path(plots_dir, "Day5_Area_All_Vessels_Violin.png"),
  p_all_vessels,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(plots_dir, "Day5_Area_All_Vessels_Violin.svg"),
  p_all_vessels,
  width = 8.5,
  height = 6.5
)

ggsave(
  file.path(plots_dir, "Day5_Area_By_BBB_State_Violin_3Panel.png"),
  p_bbb_vessels,
  width = 15,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(plots_dir, "Day5_Area_By_BBB_State_Violin_3Panel.svg"),
  p_bbb_vessels,
  width = 15,
  height = 6.5
)

ggsave(
  file.path(plots_dir, "Day5_Area_All_Vessels_Mouse_Level_Boxplot.png"),
  p_mouse_all,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(plots_dir, "Day5_Area_All_Vessels_Mouse_Level_Boxplot.svg"),
  p_mouse_all,
  width = 8.5,
  height = 6.5
)

ggsave(
  file.path(plots_dir, "Day5_Area_Mouse_Level_By_BBB_State_3Panel.png"),
  p_mouse_bbb,
  width = 15,
  height = 6.5,
  dpi = 300
)

ggsave(
  file.path(plots_dir, "Day5_Area_Mouse_Level_By_BBB_State_3Panel.svg"),
  p_mouse_bbb,
  width = 15,
  height = 6.5
)

message(
  "Analysis complete. Results and plots saved to: ",
  output_dir
)


#Avastin treatment was associated with reduced calibrated vessel area across BBB-stability states, with the strongest statistical evidence observed in the Low-BBB group.










#IgG----

# ---- Day 5 IgG distributions, quantile cutoffs, and BBB-state density plots

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

# ---- Output folders

output_dir <- file.path(
  "Arpan Multiplex",
  "Day5_Area_Normalized",
  "IgG_Quantile_Distributions"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Prepare IgG data

igg_data <- day5_all_vessels_long %>%
  mutate(
    Treatment = case_when(
      grepl("saline|pbs|vehicle", Treatment, ignore.case = TRUE) ~ "Saline",
      grepl("avastin", Treatment, ignore.case = TRUE) ~ "Avastin",
      TRUE ~ as.character(Treatment)
    ),
    IgG_value = as.numeric(IgG_PV_Mean_bg_corrected),
    Sample_ID = as.character(Sample_ID)
  ) %>%
  filter(
    !is.na(IgG_value),
    !is.na(Sample_ID),
    !is.na(Treatment),
    Treatment %in% c("Saline", "Avastin")
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Saline", "Avastin"))
  )

# ---- Confirm sample counts

sample_counts <- igg_data %>%
  group_by(Treatment, Sample_ID) %>%
  summarise(
    n_vessels = n(),
    .groups = "drop"
  ) %>%
  arrange(Treatment, Sample_ID)

print(sample_counts)

# ---- Establish 3-column layout within each treatment
# Each treatment should have three samples.
# Sample_slot determines the column; Sample_ID remains the true sample identifier.

sample_map <- sample_counts %>%
  arrange(Treatment, Sample_ID) %>%
  group_by(Treatment) %>%
  mutate(
    Sample_slot = paste0("Sample ", row_number())
  ) %>%
  ungroup()

if (any(table(sample_map$Treatment) != 3)) {
  warning("At least one treatment does not have exactly three samples. Check sample_counts.")
}

igg_data <- igg_data %>%
  left_join(
    sample_map %>% select(Treatment, Sample_ID, Sample_slot),
    by = c("Treatment", "Sample_ID")
  ) %>%
  mutate(
    Sample_slot = factor(
      Sample_slot,
      levels = c("Sample 1", "Sample 2", "Sample 3")
    )
  )

# ---- Calculate Q25, median, and Q75 independently per sample

igg_cutoffs <- igg_data %>%
  group_by(Treatment, Sample_ID, Sample_slot) %>%
  summarise(
    Q25 = quantile(IgG_value, 0.25, na.rm = TRUE, names = FALSE),
    Median = median(IgG_value, na.rm = TRUE),
    Q75 = quantile(IgG_value, 0.75, na.rm = TRUE, names = FALSE),
    Min_IgG = min(IgG_value, na.rm = TRUE),
    Max_IgG = max(IgG_value, na.rm = TRUE),
    n_vessels = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Saline", "Avastin")),
    Sample_label = paste(Treatment, Sample_ID, sep = " | ")
  )

write_csv(
  igg_cutoffs,
  file.path(output_dir, "Day5_IgG_Q25_Q50_Q75_by_sample.csv")
)

# ---- Assign BBB stability groups using sample-specific cutoffs
# High: IgG <= Q25
# Mid:  Q25 < IgG < Q75
# Low:  IgG >= Q75

igg_data <- igg_data %>%
  left_join(
    igg_cutoffs %>%
      select(Treatment, Sample_ID, Q25, Median, Q75),
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

# ---- Save counts by treatment, sample, and BBB group

bbb_counts <- igg_data %>%
  count(Treatment, Sample_ID, BBB_group, name = "n_vessels") %>%
  complete(
    Treatment,
    Sample_ID,
    BBB_group,
    fill = list(n_vessels = 0)
  ) %>%
  arrange(Treatment, Sample_ID, BBB_group)

write_csv(
  bbb_counts,
  file.path(output_dir, "Day5_IgG_BBB_group_counts_by_sample.csv")
)

print(bbb_counts)

# ---- Prepare quantile dots for histograms

# Calculate histogram heights to place cutoff dots near the baseline.
histogram_heights <- igg_data %>%
  group_by(Treatment, Sample_ID, Sample_slot) %>%
  summarise(
    max_bin_count = max(
      hist(IgG_value, breaks = 35, plot = FALSE)$counts
    ),
    .groups = "drop"
  ) %>%
  mutate(
    dot_y = pmax(max_bin_count * 0.06, 0.5)
  )

cutoff_dots <- igg_cutoffs %>%
  pivot_longer(
    cols = c(Q25, Median, Q75),
    names_to = "Quantile",
    values_to = "Cutoff"
  ) %>%
  left_join(
    histogram_heights,
    by = c("Treatment", "Sample_ID", "Sample_slot")
  ) %>%
  mutate(
    Quantile = factor(
      Quantile,
      levels = c("Q25", "Median", "Q75")
    )
  )

quantile_colors <- c(
  Q25 = "#0072B2",
  Median = "black",
  Q75 = "#D55E00"
)

quantile_shapes <- c(
  Q25 = 16,
  Median = 18,
  Q75 = 17
)

# ---- Plot 1: Histograms in a 3-column x 2-row grid 
# Top row: Saline; bottom row: Avastin.
# Each column represents the first, second, or third sample within treatment.

igg_histogram <- ggplot(
  igg_data,
  aes(x = IgG_value)
) +
  geom_histogram(
    bins = 35,
    fill = "grey75",
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "black",
    linewidth = 0.7
  ) +
  geom_vline(
    data = igg_cutoffs,
    aes(xintercept = Q25),
    inherit.aes = FALSE,
    color = "#0072B2",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    data = igg_cutoffs,
    aes(xintercept = Q75),
    inherit.aes = FALSE,
    color = "#D55E00",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_point(
    data = cutoff_dots,
    aes(
      x = Cutoff,
      y = dot_y,
      color = Quantile,
      shape = Quantile
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  facet_grid(
    Treatment ~ Sample_slot,
    scales = "free_y",
    labeller = labeller(
      Sample_slot = function(x) {
        lookup <- sample_map %>%
          distinct(Sample_slot, Sample_ID) %>%
          mutate(Sample_slot = as.character(Sample_slot))
        ids <- lookup$Sample_ID[match(x, lookup$Sample_slot)]
        ifelse(is.na(ids), x, ids)
      }
    )
  ) +
  scale_color_manual(
    values = quantile_colors,
    labels = c(
      Q25 = "Q25",
      Median = "Median",
      Q75 = "Q75"
    )
  ) +
  scale_shape_manual(
    values = quantile_shapes,
    labels = c(
      Q25 = "Q25",
      Median = "Median",
      Q75 = "Q75"
    )
  ) +
  labs(
    title = "Day 5 Perivascular IgG Distributions",
    subtitle = paste(
      "Saline top row; Avastin bottom row.",
      "Blue/red dashed lines: Q25/Q75; dots: Q25, median, Q75."
    ),
    x = "Perivascular IgG mean (background-corrected)",
    y = "Number of vessels",
    color = "Quantile",
    shape = "Quantile"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(igg_histogram)

ggsave(
  file.path(output_dir, "Day5_IgG_Histograms_3x2.png"),
  igg_histogram,
  width = 15,
  height = 9,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Day5_IgG_Histograms_3x2.svg"),
  igg_histogram,
  width = 15,
  height = 9,
  bg = "white"
)

# ---- Prepare background zones for density plots
# Zones show the sample-specific IgG ranges assigned to High, Mid, and Low.
# They are visual guides for the same quantile-based classification above.

density_zones <- bind_rows(
  igg_cutoffs %>%
    transmute(
      Treatment, Sample_ID, Sample_slot,
      BBB_group = "High",
      xmin = Min_IgG,
      xmax = Q25
    ),
  igg_cutoffs %>%
    transmute(
      Treatment, Sample_ID, Sample_slot,
      BBB_group = "Mid",
      xmin = Q25,
      xmax = Q75
    ),
  igg_cutoffs %>%
    transmute(
      Treatment, Sample_ID, Sample_slot,
      BBB_group = "Low",
      xmin = Q75,
      xmax = Max_IgG
    )
) %>%
  mutate(
    BBB_group = factor(
      BBB_group,
      levels = c("High", "Mid", "Low")
    )
  ) %>%
  filter(xmax > xmin)

bbb_colors <- c(
  High = "#009E73",
  Mid = "#E6AB02",
  Low = "#D55E00"
)

# ---- Plot 2: Density distributions with High/Mid/Low BBB zones 
# Each sample has its own density curve.
# The background bands indicate the sample-specific quantile-defined zones.

igg_density <- ggplot() +
  geom_rect(
    data = density_zones,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = 0,
      ymax = Inf,
      fill = BBB_group
    ),
    inherit.aes = FALSE,
    alpha = 0.16
  ) +
  geom_density(
    data = igg_data,
    aes(
      x = IgG_value,
      group = interaction(Treatment, Sample_ID)
    ),
    color = "black",
    fill = "grey65",
    alpha = 0.55,
    linewidth = 0.8,
    adjust = 1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "black",
    linewidth = 0.7
  ) +
  geom_vline(
    data = igg_cutoffs,
    aes(xintercept = Q25),
    inherit.aes = FALSE,
    color = "#0072B2",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    data = igg_cutoffs,
    aes(xintercept = Q75),
    inherit.aes = FALSE,
    color = "#D55E00",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_grid(
    Treatment ~ Sample_slot,
    labeller = labeller(
      Sample_slot = function(x) {
        lookup <- sample_map %>%
          distinct(Sample_slot, Sample_ID) %>%
          mutate(Sample_slot = as.character(Sample_slot))
        ids <- lookup$Sample_ID[match(x, lookup$Sample_slot)]
        ifelse(is.na(ids), x, ids)
      }
    )
  ) +
  scale_fill_manual(
    values = bbb_colors,
    breaks = c("High", "Mid", "Low"),
    labels = c(
      High = "High BBB stability",
      Mid = "Mid BBB stability",
      Low = "Low BBB stability"
    ),
    name = "Quantile-defined BBB group"
  ) +
  labs(
    title = "Day 5 Perivascular IgG Density Distributions",
    subtitle = paste(
      "Shaded regions show sample-specific quantile-defined BBB groups.",
      "Blue/red dashed lines: Q25/Q75; dotted line: zero."
    ),
    x = "Perivascular IgG mean (background-corrected)",
    y = "Density"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(igg_density)

ggsave(
  file.path(output_dir, "Day5_IgG_Density_BBB_Groups_3x2.png"),
  igg_density,
  width = 15,
  height = 9,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Day5_IgG_Density_BBB_Groups_3x2.svg"),
  igg_density,
  width = 15,
  height = 9,
  bg = "white"
)

# ---- Plot 3: Quantile cutoff dot plot 

cutoff_dotplot_data <- igg_cutoffs %>%
  select(
    Treatment, Sample_ID, Sample_slot,
    Q25, Median, Q75
  ) %>%
  pivot_longer(
    cols = c(Q25, Median, Q75),
    names_to = "Quantile",
    values_to = "Cutoff"
  ) %>%
  mutate(
    Quantile = factor(
      Quantile,
      levels = c("Q25", "Median", "Q75")
    )
  )

igg_quantile_dotplot <- ggplot(
  cutoff_dotplot_data,
  aes(
    x = Cutoff,
    y = Sample_ID,
    color = Quantile,
    shape = Quantile
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "grey35",
    linewidth = 0.7
  ) +
  geom_point(
    size = 4,
    position = position_dodge(width = 0.45)
  ) +
  facet_grid(
    Treatment ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_color_manual(values = quantile_colors) +
  scale_shape_manual(values = quantile_shapes) +
  labs(
    title = "Day 5 IgG Quantile Cutoffs by Sample",
    subtitle = "Cutoffs calculated independently within each sample",
    x = "Perivascular IgG mean (background-corrected)",
    y = NULL,
    color = "Quantile",
    shape = "Quantile"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.text = element_text(size = 11),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(igg_quantile_dotplot)

ggsave(
  file.path(output_dir, "Day5_IgG_Quantile_Cutoffs_Dotplot.png"),
  igg_quantile_dotplot,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Day5_IgG_Quantile_Cutoffs_Dotplot.svg"),
  igg_quantile_dotplot,
  width = 12,
  height = 8,
  bg = "white"
)

# ---- Print exact cutoffs for inspection 

igg_cutoffs %>%
  select(
    Treatment, Sample_ID, n_vessels,
    Q25, Median, Q75
  ) %>%
  arrange(Treatment, Sample_ID) %>%
  print()


#clean density plot----

# ---- Day 5: IgG density by treatment, tumor, and BBB state ----

library(tidyverse)
library(ggplot2)

# ---- Output folder ----

igg_density_dir <- file.path(
  "Arpan Multiplex",
  "Day5_Area_Normalized",
  "IgG_Quantile_Distributions"
)

dir.create(
  igg_density_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ---- Prepare IgG data ----

igg_data <- day5_all_vessels_long %>%
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
    title = "Day 5: Perivascular IgG Density by BBB State",
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
    "Day5_IgG_Density_By_BBB_State.png"
  ),
  plot = igg_density_by_state,
  width = 18,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ---- Export SVG ----

ggsave(
  filename = file.path(
    igg_density_dir,
    "Day5_IgG_Density_By_BBB_State.svg"
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



#IQR plot----

# ---- Day 5: Mouse-level IgG IQR plot and statistical analysis

library(dplyr)
library(ggplot2)
library(readr)

# ---- Output directory

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Calculate IgG spread separately for each mouse

igg_mouse_iqr <- igg_data %>%
  group_by(Treatment, Sample_ID, Sample_slot) %>%
  summarise(
    n_vessels = sum(is.finite(IgG_value)),
    Q25 = quantile(
      IgG_value[is.finite(IgG_value)],
      probs = 0.25,
      names = FALSE
    ),
    Median_IgG = median(
      IgG_value[is.finite(IgG_value)]
    ),
    Q75 = quantile(
      IgG_value[is.finite(IgG_value)],
      probs = 0.75,
      names = FALSE
    ),
    IQR_IgG = IQR(
      IgG_value[is.finite(IgG_value)]
    ),
    SD_IgG = sd(
      IgG_value[is.finite(IgG_value)]
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("Saline", "Avastin")
    ),
    Sample_slot = factor(
      Sample_slot,
      levels = c("Sample 1", "Sample 2", "Sample 3")
    )
  ) %>%
  arrange(Treatment, Sample_slot)

# ---- Inspect mouse-level results

print(igg_mouse_iqr, n = Inf)

# ---- Save mouse-level results

write_csv(
  igg_mouse_iqr,
  file.path(output_dir, "Day5_IgG_Mouse_Level_IQR.csv")
)

# ---- Descriptive statistics by treatment

igg_iqr_treatment_summary <- igg_mouse_iqr %>%
  group_by(Treatment) %>%
  summarise(
    n_mice = n(),
    Mean_IQR = mean(IQR_IgG, na.rm = TRUE),
    SD_IQR = sd(IQR_IgG, na.rm = TRUE),
    Median_IQR = median(IQR_IgG, na.rm = TRUE),
    Min_IQR = min(IQR_IgG, na.rm = TRUE),
    Max_IQR = max(IQR_IgG, na.rm = TRUE),
    Mean_SD_within_mouse = mean(SD_IgG, na.rm = TRUE),
    .groups = "drop"
  )

print(igg_iqr_treatment_summary)

write_csv(
  igg_iqr_treatment_summary,
  file.path(output_dir, "Day5_IgG_IQR_Treatment_Summary.csv")
)

# ---- Calculate treatment effect sizes

saline_iqr <- igg_mouse_iqr %>%
  filter(Treatment == "Saline") %>%
  pull(IQR_IgG)

avastin_iqr <- igg_mouse_iqr %>%
  filter(Treatment == "Avastin") %>%
  pull(IQR_IgG)

stopifnot(
  length(saline_iqr) == 3,
  length(avastin_iqr) == 3,
  all(is.finite(c(saline_iqr, avastin_iqr)))
)

mean_saline_iqr <- mean(saline_iqr)
mean_avastin_iqr <- mean(avastin_iqr)

observed_difference <- mean_avastin_iqr - mean_saline_iqr
observed_ratio <- mean_avastin_iqr / mean_saline_iqr
observed_percent_reduction <- 100 * (1 - observed_ratio)

# ---- Exact permutation test on mouse-level IQRs
# Test statistic: mean Avastin IQR - mean Saline IQR.
# Enumerate all 20 allocations of 3 of the 6 mice to Saline.
# This assumes mouse treatment labels are exchangeable under the null.

all_iqr <- igg_mouse_iqr$IQR_IgG
group_allocations <- combn(seq_along(all_iqr), 3)

permutation_differences <- apply(
  group_allocations,
  2,
  function(saline_indices) {
    avastin_indices <- setdiff(
      seq_along(all_iqr),
      saline_indices
    )
    
    mean(all_iqr[avastin_indices]) -
      mean(all_iqr[saline_indices])
  }
)

# Two-sided test
p_exact_two_sided <- mean(
  abs(permutation_differences) >=
    abs(observed_difference) - 1e-10
)

# One-sided alternative: Avastin has lower IQR than Saline.
# Use this as the primary directional test only if prespecified.
p_exact_one_sided <- mean(
  permutation_differences <= observed_difference + 1e-10
)

# ---- Combine effect sizes and permutation results

igg_iqr_test_results <- tibble(
  Mean_IQR_Saline = mean_saline_iqr,
  Mean_IQR_Avastin = mean_avastin_iqr,
  Median_IQR_Saline = median(saline_iqr),
  Median_IQR_Avastin = median(avastin_iqr),
  Difference_Avastin_minus_Saline = observed_difference,
  Ratio_Avastin_over_Saline = observed_ratio,
  Percent_Reduction_in_Mean_IQR = observed_percent_reduction,
  Exact_Permutation_P_Two_Sided = p_exact_two_sided,
  Exact_Permutation_P_One_Sided = p_exact_one_sided,
  N_Saline_Mice = length(saline_iqr),
  N_Avastin_Mice = length(avastin_iqr)
)

print(igg_iqr_test_results, width = Inf)

write_csv(
  igg_iqr_test_results,
  file.path(output_dir, "Day5_IgG_Mouse_IQR_Exact_Permutation_Test.csv")
)

# ---- Prepare treatment-level mean and SD for plotting

igg_iqr_plot_summary <- igg_mouse_iqr %>%
  group_by(Treatment) %>%
  summarise(
    Mean_IQR = mean(IQR_IgG, na.rm = TRUE),
    SD_IQR = sd(IQR_IgG, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("Saline", "Avastin")
    ),
    SD_lower = pmax(0, Mean_IQR - SD_IQR),
    SD_upper = Mean_IQR + SD_IQR
  )

# ---- Assign fixed x positions for consistent point and label alignment

igg_mouse_iqr_plot_data <- igg_mouse_iqr %>%
  mutate(
    x_position = ifelse(Treatment == "Saline", 1, 2),
    point_offset = case_when(
      Sample_slot == "Sample 1" ~ -0.09,
      Sample_slot == "Sample 2" ~ 0,
      Sample_slot == "Sample 3" ~ 0.09,
      TRUE ~ 0
    ),
    x_plot = x_position + point_offset
  )

igg_iqr_plot_summary <- igg_iqr_plot_summary %>%
  mutate(
    x_position = ifelse(Treatment == "Saline", 1, 2)
  )

# ---- Plot: mouse-level IQR

treatment_colors <- c(
  Saline = "#0072B2",
  Avastin = "#D55E00"
)

igg_mouse_iqr_plot <- ggplot() +
  # Each point is one independent mouse
  geom_point(
    data = igg_mouse_iqr_plot_data,
    aes(
      x = x_plot,
      y = IQR_IgG,
      color = Treatment
    ),
    size = 4
  ) +
  # Mouse labels
  geom_text(
    data = igg_mouse_iqr_plot_data,
    aes(
      x = x_plot,
      y = IQR_IgG,
      label = Sample_slot
    ),
    vjust = -1.0,
    size = 3.5,
    color = "black",
    show.legend = FALSE
  ) +
  # Mean +/- 1 SD across mice
  geom_errorbar(
    data = igg_iqr_plot_summary,
    aes(
      x = x_position,
      ymin = SD_lower,
      ymax = SD_upper
    ),
    inherit.aes = FALSE,
    width = 0.12,
    color = "black",
    linewidth = 0.7
  ) +
  # Horizontal dash marks treatment mean
  geom_segment(
    data = igg_iqr_plot_summary,
    aes(
      x = x_position - 0.18,
      xend = x_position + 0.18,
      y = Mean_IQR,
      yend = Mean_IQR
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.2
  ) +
  scale_color_manual(values = treatment_colors) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Saline", "Avastin"),
    limits = c(0.65, 2.35)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.04, 0.20))
  ) +
  labs(
    title = "Day 5: Perivascular IgG Variability Across Mice",
    subtitle = sprintf(
      "Mean IQR reduction with Avastin: %.1f%% | Exact permutation p = %.3f (two-sided)",
      observed_percent_reduction,
      p_exact_two_sided
    ),
    x = NULL,
    y = "Within-mouse IgG interquartile range (IQR)"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 13),
    legend.position = "none"
  )

print(igg_mouse_iqr_plot)

# ---- Save PNG

ggsave(
  filename = file.path(output_dir, "Day5_IgG_Mouse_Level_IQR.png"),
  plot = igg_mouse_iqr_plot,
  device = "png",
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)

# ---- Save SVG without requiring Hmisc

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = file.path(output_dir, "Day5_IgG_Mouse_Level_IQR.svg"),
    plot = igg_mouse_iqr_plot,
    device = svglite::svglite,
    width = 8,
    height = 6,
    bg = "white"
  )
} else {
  message(
    "Package 'svglite' is not installed. ",
    "The PNG was saved; install svglite separately if SVG output is needed." 
  )
}





#MFSD2A-PLVAP Index----

# ---- Day 5: MFSD2A-PLVAP state index and co-expression analysis ----

library(dplyr)
library(tidyr)
library(ggplot2)
library(lmerTest)
library(emmeans)

# ---- Output directory
day5_state_dir <- "Arpan Multiplex/Day5_Area_Normalized/Cell_State_Plasticity"
dir.create(day5_state_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helper functions

format_p <- function(p) {
  if (length(p) == 0 || is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  paste0("= ", sprintf("%.3f", p))
}

format_ci <- function(estimate, lower, upper, digits = 2) {
  paste0(
    sprintf(paste0("%.", digits, "f"), estimate),
    " (95% CI ",
    sprintf(paste0("%.", digits, "f"), lower),
    " to ",
    sprintf(paste0("%.", digits, "f"), upper),
    ")"
  )
}

# ---- 1. Prepare the Low-BBB vessel data

stopifnot(exists("day5_all_vessels_long"))

day5_state_input <- day5_all_vessels_long

required_cols <- c("Sample_ID", "Treatment", "Marker", "Expression")
missing_cols <- setdiff(required_cols, names(day5_state_input))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in day5_all_vessels_long: ",
    paste(missing_cols, collapse = ", ")
  )
}

bbb_col <- if ("BBB Stability" %in% names(day5_state_input)) {
  "BBB Stability"
} else if ("BBB_State" %in% names(day5_state_input)) {
  "BBB_State"
} else {
  stop("Could not find 'BBB Stability' or 'BBB_State' in day5_all_vessels_long.")
}

# Use existing area-normalized expression or calculate it
if (!"Expression_per_area" %in% names(day5_state_input)) {
  if (!"Area_calibrated" %in% names(day5_state_input)) {
    stop(
      "Neither Expression_per_area nor Area_calibrated is available."
    )
  }
  
  day5_state_input <- day5_state_input %>%
    mutate(
      Expression_per_area = if_else(
        is.finite(Area_calibrated) & Area_calibrated > 0,
        Expression / Area_calibrated,
        NA_real_
      )
    )
}

day5_state_input <- day5_state_input %>%
  mutate(
    Expression_per_area_log2 = if_else(
      is.finite(Expression_per_area) & Expression_per_area >= 0,
      log2(Expression_per_area + 1),
      NA_real_
    ),
    Treatment = factor(Treatment, levels = c("Saline", "Avastin"))
  )

# Select optional metadata safely
metadata_cols <- intersect(
  c("Tumor", "ROI", "Vessel_ID"),
  names(day5_state_input)
)

day5_lowbbb_state <- day5_state_input %>%
  filter(
    .data[[bbb_col]] %in% c("Low BBB", "Low"),
    Marker %in% c("MFSD2A", "PLVAP"),
    !is.na(Sample_ID),
    !is.na(Treatment),
    is.finite(Expression_per_area_log2)
  ) %>%
  select(
    Sample_ID, Treatment,
    all_of(metadata_cols),
    all_of(bbb_col), Marker,
    Expression, Expression_per_area,
    Expression_per_area_log2
  ) %>%
  distinct()

if (nrow(day5_lowbbb_state) == 0) {
  stop(
    "No Low-BBB vessels found. Check the values in column '",
    bbb_col, "'."
  )
}

cat("\nIndependent samples by treatment:\n")
print(
  day5_lowbbb_state %>%
    distinct(Sample_ID, Treatment) %>%
    count(Treatment, name = "n_samples")
)

cat("\nVessel counts by treatment and marker:\n")
print(
  day5_lowbbb_state %>%
    distinct(Sample_ID, Treatment, Marker, Vessel_ID) %>%
    count(Treatment, Marker, name = "n_vessels")
)

# ---- 2. Build paired MFSD2A and PLVAP vessel data

id_cols <- intersect(
  c("Sample_ID", "Treatment", "Tumor", "ROI", "Vessel_ID"),
  names(day5_lowbbb_state)
)

day5_state_wide <- day5_lowbbb_state %>%
  select(
    all_of(id_cols),
    Marker, Expression_per_area, Expression_per_area_log2
  ) %>%
  pivot_wider(
    names_from = Marker,
    values_from = c(Expression_per_area, Expression_per_area_log2),
    names_sep = "__"
  ) %>%
  filter(
    is.finite(Expression_per_area_log2__MFSD2A),
    is.finite(Expression_per_area_log2__PLVAP)
  ) %>%
  mutate(
    MFSD2A_PLVAP_index =
      Expression_per_area_log2__MFSD2A -
      Expression_per_area_log2__PLVAP,
    Treatment = factor(Treatment, levels = c("Saline", "Avastin"))
  )

if (nrow(day5_state_wide) == 0) {
  stop(
    "No vessels had both MFSD2A and PLVAP measurements after filtering."
  )
}

write.csv(
  day5_state_wide,
  file.path(day5_state_dir, "Day5_LowBBB_MFSD2A_PLVAP_StateIndex_ByVessel.csv"),
  row.names = FALSE
)

# ---- 3. Summarize the index and markers by independent sample

day5_state_mouse <- day5_state_wide %>%
  group_by(Sample_ID, Treatment) %>%
  summarise(
    n_vessels = n(),
    median_index = median(MFSD2A_PLVAP_index, na.rm = TRUE),
    mean_index = mean(MFSD2A_PLVAP_index, na.rm = TRUE),
    median_MFSD2A = median(Expression_per_area_log2__MFSD2A, na.rm = TRUE),
    median_PLVAP = median(Expression_per_area_log2__PLVAP, na.rm = TRUE),
    .groups = "drop"
  )

print(day5_state_mouse)

write.csv(
  day5_state_mouse,
  file.path(day5_state_dir, "Day5_LowBBB_StateIndex_MouseLevel.csv"),
  row.names = FALSE
)

# ---- 4. Exact two-sided permutation test using sample-level medians

exact_permutation_test <- function(data, outcome) {
  
  dat <- data %>%
    select(Sample_ID, Treatment, value = all_of(outcome)) %>%
    filter(is.finite(value)) %>%
    distinct()
  
  if (anyDuplicated(dat$Sample_ID)) {
    stop("Each Sample_ID must have exactly one value for the permutation test.")
  }
  
  dat$Treatment <- as.character(dat$Treatment)
  
  saline <- dat$value[dat$Treatment == "Saline"]
  avastin <- dat$value[dat$Treatment == "Avastin"]
  
  if (length(saline) != 3 || length(avastin) != 3) {
    warning(
      "Expected 3 samples per group; found Saline n = ",
      length(saline), ", Avastin n = ", length(avastin), "."
    )
  }
  
  if (length(saline) == 0 || length(avastin) == 0) {
    stop("Both treatment groups must have at least one sample.")
  }
  
  observed_difference <- mean(avastin) - mean(saline)
  all_values <- dat$value
  n_saline <- length(saline)
  
  allocations <- combn(seq_along(all_values), n_saline)
  
  perm_differences <- apply(allocations, 2, function(idx_saline) {
    idx_avastin <- setdiff(seq_along(all_values), idx_saline)
    mean(all_values[idx_avastin]) - mean(all_values[idx_saline])
  })
  
  p_exact <- mean(
    abs(perm_differences) >= abs(observed_difference) - 1e-12
  )
  
  tibble(
    outcome = outcome,
    n_Saline = length(saline),
    n_Avastin = length(avastin),
    mean_Saline = mean(saline),
    mean_Avastin = mean(avastin),
    difference_Avastin_minus_Saline = observed_difference,
    p_exact_two_sided = p_exact,
    n_permutations = length(perm_differences)
  )
}

day5_state_index_perm <- exact_permutation_test(
  day5_state_mouse, "median_index"
)

print(day5_state_index_perm)

# ---- 5. Mixed model for the composite state index
# Same basic model structure as for individual markers:
# index ~ Treatment + (1 | Sample_ID)

day5_state_lmm <- lmerTest::lmer(
  MFSD2A_PLVAP_index ~ Treatment + (1 | Sample_ID),
  data = day5_state_wide,
  REML = TRUE
)

print(summary(day5_state_lmm))
print(anova(day5_state_lmm))

day5_state_lmm_emm <- emmeans(
  day5_state_lmm,
  ~ Treatment,
  lmer.df = "kenward-roger"
)

day5_state_lmm_contrast <- as.data.frame(
  contrast(
    day5_state_lmm_emm,
    method = list("Avastin - Saline" = c(-1, 1))
  )
)

print(day5_state_lmm_contrast)

write.csv(
  day5_state_index_perm,
  file.path(day5_state_dir, "Day5_StateIndex_ExactPermutation_Statistics.csv"),
  row.names = FALSE
)

write.csv(
  day5_state_lmm_contrast,
  file.path(day5_state_dir, "Day5_StateIndex_MixedModel_Contrast.csv"),
  row.names = FALSE
)

# ---- 6. Prepare statistics label for the index plot

index_lmm_p <- day5_state_lmm_contrast$p.value[1]
index_lmm_est <- day5_state_lmm_contrast$estimate[1]
index_lmm_lower <- day5_state_lmm_contrast$lower.CL[1]
index_lmm_upper <- day5_state_lmm_contrast$upper.CL[1]
index_perm_p <- day5_state_index_perm$p_exact_two_sided[1]

index_stats_label <- paste0(
  "Mixed model, Avastin - Saline: ",
  format_ci(index_lmm_est, index_lmm_lower, index_lmm_upper),
  "; p ", format_p(index_lmm_p),
  "\nExact two-sided permutation on sample medians: p ",
  format_p(index_perm_p),
  "\nIndependent samples: Saline n = ",
  day5_state_index_perm$n_Saline[1],
  ", Avastin n = ", day5_state_index_perm$n_Avastin[1]
)

# ---- 7. Plot: state index with statistics

day5_index_plot <- ggplot() +
  geom_violin(
    data = day5_state_wide,
    aes(x = Treatment, y = MFSD2A_PLVAP_index, fill = Treatment),
    alpha = 0.13,
    trim = TRUE,
    color = NA
  ) +
  geom_point(
    data = day5_state_wide,
    aes(x = Treatment, y = MFSD2A_PLVAP_index, color = Treatment),
    position = position_jitter(width = 0.10, height = 0),
    alpha = 0.18,
    size = 1.5
  ) +
  geom_point(
    data = day5_state_mouse,
    aes(x = Treatment, y = median_index, shape = Sample_ID, color = Treatment),
    size = 4.2,
    stroke = 1.1
  ) +
  scale_color_manual(
    values = c("Saline" = "#0072B2", "Avastin" = "#D55E00")
  ) +
  scale_fill_manual(
    values = c("Saline" = "#0072B2", "Avastin" = "#D55E00")
  ) +
  labs(
    title = "MFSD2A-PLVAP state index in Low-BBB vessels",
    subtitle = "Faint points = vessels; large symbols = independent sample medians",
    x = NULL,
    y = "MFSD2A-PLVAP index\n(log2 area-normalized expression + 1)",
    color = "Treatment",
    fill = "Treatment",
    shape = "Mouse/sample",
    caption = index_stats_label
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(
      hjust = 0,
      size = 10.5,
      lineheight = 1.2
    ),
    plot.caption.position = "plot",
    legend.position = "right",
    plot.margin = margin(10, 12, 10, 12)
  )

print(day5_index_plot)

ggsave(
  file.path(day5_state_dir, "Day5_LowBBB_MFSD2A_PLVAP_StateIndex_WITH_STATS.png"),
  day5_index_plot, width = 9, height = 7, dpi = 400
)

ggsave(
  file.path(day5_state_dir, "Day5_LowBBB_MFSD2A_PLVAP_StateIndex_WITH_STATS.svg"),
  day5_index_plot, width = 9, height = 7
)

# ---- 8. Marker-specific mixed models

fit_day5_marker_lmm <- function(marker_name) {
  
  outcome <- switch(
    marker_name,
    "MFSD2A" = "Expression_per_area_log2__MFSD2A",
    "PLVAP" = "Expression_per_area_log2__PLVAP",
    stop("Unsupported marker: ", marker_name)
  )
  
  model <- lmerTest::lmer(
    reformulate(
      c("Treatment", "(1 | Sample_ID)"),
      response = outcome
    ),
    data = day5_state_wide,
    REML = TRUE
  )
  
  emm <- emmeans(
    model,
    ~ Treatment,
    lmer.df = "kenward-roger"
  )
  
  contrast_df <- as.data.frame(
    contrast(
      emm,
      method = list("Avastin - Saline" = c(-1, 1))
    )
  )
  
  contrast_df$Marker <- marker_name
  contrast_df
}

day5_marker_lmm_stats <- bind_rows(
  fit_day5_marker_lmm("MFSD2A"),
  fit_day5_marker_lmm("PLVAP")
) %>%
  mutate(
    p_Holm = p.adjust(p.value, method = "holm")
  )

print(day5_marker_lmm_stats)

write.csv(
  day5_marker_lmm_stats,
  file.path(day5_state_dir, "Day5_LowBBB_Marker_MixedModel_Contrasts.csv"),
  row.names = FALSE
)

# ---- 9. Sample-level marker medians and exact permutation tests

day5_marker_mouse <- day5_state_wide %>%
  group_by(Sample_ID, Treatment) %>%
  summarise(
    median_MFSD2A = median(Expression_per_area_log2__MFSD2A, na.rm = TRUE),
    median_PLVAP = median(Expression_per_area_log2__PLVAP, na.rm = TRUE),
    .groups = "drop"
  )

day5_marker_perm <- bind_rows(
  exact_permutation_test(day5_marker_mouse, "median_MFSD2A"),
  exact_permutation_test(day5_marker_mouse, "median_PLVAP")
) %>%
  mutate(p_exact_Holm = p.adjust(p_exact_two_sided, method = "holm"))

print(day5_marker_perm)

write.csv(
  day5_marker_perm,
  file.path(day5_state_dir, "Day5_LowBBB_Marker_MouseMedian_Permutation.csv"),
  row.names = FALSE
)

# ---- 10. Co-expression scatterplot sample summaries

day5_scatter_mouse <- day5_state_wide %>%
  group_by(Sample_ID, Treatment) %>%
  summarise(
    n_vessels = n(),
    MFSD2A_median = median(Expression_per_area_log2__MFSD2A, na.rm = TRUE),
    PLVAP_median = median(Expression_per_area_log2__PLVAP, na.rm = TRUE),
    .groups = "drop"
  )

# Descriptive sample-level Spearman correlation (six samples)
day5_mouse_spearman <- cor.test(
  day5_scatter_mouse$MFSD2A_median,
  day5_scatter_mouse$PLVAP_median,
  method = "spearman",
  exact = FALSE
)

print(day5_mouse_spearman)

# ---- 11. Build scatterplot statistics label

mfsd2a_stats <- day5_marker_lmm_stats %>%
  filter(Marker == "MFSD2A")

plvap_stats <- day5_marker_lmm_stats %>%
  filter(Marker == "PLVAP")

scatter_stats_label <- paste0(
  "Treatment effect from marker-specific mixed models (Avastin - Saline):",
  "\nMFSD2A estimate = ",
  format_ci(
    mfsd2a_stats$estimate[1],
    mfsd2a_stats$lower.CL[1],
    mfsd2a_stats$upper.CL[1]
  ),
  "; Holm-adjusted p ", format_p(mfsd2a_stats$p_Holm[1]),
  "\nPLVAP estimate = ",
  format_ci(
    plvap_stats$estimate[1],
    plvap_stats$lower.CL[1],
    plvap_stats$upper.CL[1]
  ),
  "; Holm-adjusted p ", format_p(plvap_stats$p_Holm[1]),
  "\nSample-level Spearman correlation: rho = ",
  sprintf("%.2f", unname(day5_mouse_spearman$estimate)),
  "; p ", format_p(day5_mouse_spearman$p.value),
  " (six sample-level medians)"
)

# ---- 12. Plot: MFSD2A versus PLVAP co-expression with statistics

day5_coplot <- ggplot() +
  geom_point(
    data = day5_state_wide,
    aes(
      x = Expression_per_area_log2__MFSD2A,
      y = Expression_per_area_log2__PLVAP,
      color = Treatment
    ),
    alpha = 0.20,
    size = 1.7
  ) +
  geom_point(
    data = day5_scatter_mouse,
    aes(
      x = MFSD2A_median,
      y = PLVAP_median,
      color = Treatment,
      shape = Sample_ID
    ),
    size = 4.5,
    stroke = 1.2
  ) +
  scale_color_manual(
    values = c("Saline" = "#0072B2", "Avastin" = "#D55E00")
  ) +
  labs(
    title = "MFSD2A-PLVAP co-expression in Low-BBB vessels",
    subtitle = "Faint points = vessels; large symbols = independent sample medians",
    x = "MFSD2A (log2 area-normalized expression + 1)",
    y = "PLVAP (log2 area-normalized expression + 1)",
    color = "Treatment",
    shape = "Mouse/sample",
    caption = scatter_stats_label
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(
      hjust = 0,
      size = 10,
      lineheight = 1.2
    ),
    plot.caption.position = "plot",
    legend.position = "right",
    plot.margin = margin(10, 12, 10, 12)
  )

print(day5_coplot)

ggsave(
  file.path(day5_state_dir, "Day5_LowBBB_MFSD2A_vs_PLVAP_Scatter_WITH_STATS.png"),
  day5_coplot, width = 10, height = 8, dpi = 400
)

ggsave(
  file.path(day5_state_dir, "Day5_LowBBB_MFSD2A_vs_PLVAP_Scatter_WITH_STATS.svg"),
  day5_coplot, width = 10, height = 8
)

# ---- 13. Descriptive vessel-level Spearman correlations by treatment
# These are descriptive only; vessels are nested within samples.

day5_spearman_by_treatment <- day5_state_wide %>%
  group_by(Treatment) %>%
  summarise(
    n_vessels = n(),
    spearman_rho = cor(
      Expression_per_area_log2__MFSD2A,
      Expression_per_area_log2__PLVAP,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )

print(day5_spearman_by_treatment)

write.csv(
  day5_spearman_by_treatment,
  file.path(day5_state_dir, "Day5_LowBBB_Coexpression_DescriptiveCorrelations.csv"),
  row.names = FALSE
)

# ---- 14. Final summary

cat("\nMouse/sample-level state-index medians:\n")
print(day5_state_mouse)

cat("\nExact two-sided permutation test for state index:\n")
print(day5_state_index_perm)

cat("\nMixed-model treatment contrast for state index:\n")
print(day5_state_lmm_contrast)

cat("\nMarker-specific mixed-model contrasts:\n")
print(day5_marker_lmm_stats)

cat("\nSample-level Spearman correlation:\n")
print(day5_mouse_spearman)

cat("\nAll output files saved to:\n")
cat(normalizePath(day5_state_dir), "\n")



#Correlations----


# ---- MFSD2A vs PLVAP scatterplots and treatment-specific correlations ----

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)
library(stringr)

# ---- Settings ----

scatter_dir <- "Arpan Multiplex/Day5_Area_Normalized/scatterplots"

dir.create(
  scatter_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(exists("day5_all_vessels_long"))

# ---- Identify required columns ----

dat <- day5_all_vessels_long

bbb_candidates <- c("BBB Stability", "BBB_State", "BBB_group", "BBB_Group")
bbb_col <- bbb_candidates[bbb_candidates %in% names(dat)][1]

if (length(bbb_col) == 0 || is.na(bbb_col)) {
  stop("Could not identify the BBB-state column.")
}

if (!"Sample_ID" %in% names(dat)) {
  stop("Sample_ID is missing.")
}

if (!"Treatment" %in% names(dat)) {
  stop("Treatment is missing.")
}

if (!"Marker" %in% names(dat)) {
  stop("Marker is missing.")
}

if (!"Vessel_ID" %in% names(dat)) {
  stop("Vessel_ID is required to pair MFSD2A and PLVAP within vessels.")
}

if (!"Expression_per_area" %in% names(dat)) {
  stop(
    "Expression_per_area is missing. ",
    "Use the area-normalized dataset for this analysis."
  )
}

# ---- Prepare and standardize data ----

dat <- dat %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    Vessel_ID = as.character(Vessel_ID),
    Treatment = case_when(
      str_to_lower(as.character(Treatment)) %in%
        c("saline", "pbs", "vehicle") ~ "Saline",
      str_to_lower(as.character(Treatment)) == "avastin" ~ "Avastin",
      TRUE ~ as.character(Treatment)
    ),
    BBB_raw = as.character(.data[[bbb_col]]),
    BBB_State = case_when(
      str_detect(str_to_lower(BBB_raw), "low") ~ "Low BBB",
      str_detect(str_to_lower(BBB_raw), "mid|intermediate") ~ "Mid BBB",
      str_detect(str_to_lower(BBB_raw), "high") ~ "High BBB",
      TRUE ~ NA_character_
    ),
    Marker = as.character(Marker),
    Expression_per_area = as.numeric(Expression_per_area)
  ) %>%
  filter(
    Treatment %in% c("Saline", "Avastin"),
    Marker %in% c("MFSD2A", "PLVAP"),
    !is.na(Sample_ID),
    !is.na(Vessel_ID),
    !is.na(Expression_per_area),
    is.finite(Expression_per_area),
    Expression_per_area >= 0
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Saline", "Avastin")),
    BBB_State = factor(
      BBB_State,
      levels = c("Low BBB", "Mid BBB", "High BBB")
    ),
    Expression_log2 = log2(Expression_per_area + 1)
  )

if (nrow(dat) == 0) {
  stop("No eligible MFSD2A or PLVAP observations were found.")
}

# ---- Pair MFSD2A and PLVAP within the same vessel ----
# Median is used if duplicate rows exist for a vessel-marker combination.

vessel_wide <- dat %>%
  group_by(Sample_ID, Treatment, BBB_State, Vessel_ID, Marker) %>%
  summarise(
    Expression_log2 = median(Expression_log2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Marker,
    values_from = Expression_log2
  ) %>%
  filter(
    !is.na(MFSD2A),
    !is.na(PLVAP),
    is.finite(MFSD2A),
    is.finite(PLVAP)
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Saline", "Avastin"))
  )

if (nrow(vessel_wide) == 0) {
  stop("No vessels have paired MFSD2A and PLVAP measurements.")
}

write_csv(
  vessel_wide,
  file.path(scatter_dir, "paired_vessel_MFSD2A_PLVAP.csv")
)

# ---- Correlation helper ----

get_correlation <- function(df, method = "spearman") {
  
  df <- df %>%
    filter(
      is.finite(MFSD2A),
      is.finite(PLVAP)
    )
  
  n <- nrow(df)
  
  if (n < 3 ||
      dplyr::n_distinct(df$MFSD2A) < 2 ||
      dplyr::n_distinct(df$PLVAP) < 2) {
    return(tibble(
      n = n,
      correlation = NA_real_,
      p_value = NA_real_
    ))
  }
  
  test <- suppressWarnings(
    cor.test(
      df$MFSD2A,
      df$PLVAP,
      method = method,
      exact = FALSE
    )
  )
  
  tibble(
    n = n,
    correlation = unname(test$estimate),
    p_value = test$p.value
  )
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# ---- 1. Pooled-vessel correlations by treatment ----
# These are descriptive: vessels from one mouse are not independent
# biological replicates.

pooled_overall <- vessel_wide %>%
  group_by(Treatment) %>%
  group_modify(~ {
    spearman <- get_correlation(.x, "spearman")
    pearson <- get_correlation(.x, "pearson")
    
    tibble(
      n_vessels = nrow(.x),
      spearman_rho = spearman$correlation,
      spearman_p = spearman$p_value,
      pearson_r = pearson$correlation,
      pearson_p = pearson$p_value
    )
  }) %>%
  ungroup()

write_csv(
  pooled_overall,
  file.path(scatter_dir, "pooled_correlations_by_treatment.csv")
)

print(pooled_overall)

# ---- 2. Pooled-vessel correlations by treatment and BBB state ----

pooled_by_bbb <- vessel_wide %>%
  filter(!is.na(BBB_State)) %>%
  group_by(Treatment, BBB_State) %>%
  group_modify(~ {
    spearman <- get_correlation(.x, "spearman")
    pearson <- get_correlation(.x, "pearson")
    
    tibble(
      n_vessels = nrow(.x),
      spearman_rho = spearman$correlation,
      spearman_p = spearman$p_value,
      pearson_r = pearson$correlation,
      pearson_p = pearson$p_value
    )
  }) %>%
  ungroup()

write_csv(
  pooled_by_bbb,
  file.path(scatter_dir, "pooled_correlations_by_treatment_BBB_state.csv")
)

print(pooled_by_bbb, n = Inf)

# ---- 3. Plot pooled correlations: all BBB states combined ----

overall_labels <- pooled_overall %>%
  mutate(
    label = paste0(
      "n = ", n_vessels,
      "\nSpearman \u03c1 = ",
      ifelse(is.na(spearman_rho), "NA", sprintf("%.2f", spearman_rho)),
      ", p = ", vapply(spearman_p, format_p, character(1)),
      "\nPearson r = ",
      ifelse(is.na(pearson_r), "NA", sprintf("%.2f", pearson_r)),
      ", p = ", vapply(pearson_p, format_p, character(1))
    )
  )

p_overall <- ggplot(
  vessel_wide,
  aes(x = MFSD2A, y = PLVAP, color = Sample_ID)
) +
  geom_point(alpha = 0.45, size = 1.8) +
  geom_smooth(
    aes(group = Treatment),
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 0.9,
    inherit.aes = TRUE
  ) +
  geom_label(
    data = overall_labels,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.04,
    vjust = 1.08,
    size = 3.5,
    color = "black",
    fill = "white"
  ) +
  facet_wrap(~ Treatment, nrow = 1) +
  labs(
    title = "MFSD2A–PLVAP relationship in Low-, Mid-, and High-BBB vessels",
    subtitle = "Vessels pooled within treatment; points are colored by mouse",
    x = "MFSD2A expression per area [log2(x + 1)]",
    y = "PLVAP expression per area [log2(x + 1)]",
    color = "Mouse"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(scatter_dir, "Scatter_MFSD2A_PLVAP_Pooled_Treatment.png"),
  p_overall, width = 11, height = 6.5, dpi = 400
)

ggsave(
  file.path(scatter_dir, "Scatter_MFSD2A_PLVAP_Pooled_Treatment.svg"),
  p_overall, width = 11, height = 6.5
)

print(p_overall)

# ---- 4. Plot correlations separately by BBB state and treatment ----

bbb_labels <- pooled_by_bbb %>%
  mutate(
    label = paste0(
      "n = ", n_vessels,
      "\nSpearman \u03c1 = ",
      ifelse(is.na(spearman_rho), "NA", sprintf("%.2f", spearman_rho)),
      ", p = ", vapply(spearman_p, format_p, character(1)),
      "\nPearson r = ",
      ifelse(is.na(pearson_r), "NA", sprintf("%.2f", pearson_r)),
      ", p = ", vapply(pearson_p, format_p, character(1))
    )
  )

p_bbb <- ggplot(
  vessel_wide %>% filter(!is.na(BBB_State)),
  aes(x = MFSD2A, y = PLVAP, color = Sample_ID)
) +
  geom_point(alpha = 0.45, size = 1.5) +
  geom_smooth(
    aes(group = Treatment),
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.8
  ) +
  geom_label(
    data = bbb_labels,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.04,
    vjust = 1.08,
    size = 3.0,
    color = "black",
    fill = "white"
  ) +
  facet_grid(BBB_State ~ Treatment, drop = TRUE) +
  labs(
    title = "MFSD2A–PLVAP correlations by BBB state",
    subtitle = "Pooled vessel-level correlations are descriptive",
    x = "MFSD2A expression per area [log2(x + 1)]",
    y = "PLVAP expression per area [log2(x + 1)]",
    color = "Mouse"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(scatter_dir, "Scatter_MFSD2A_PLVAP_By_BBB_State.png"),
  p_bbb, width = 12, height = 10, dpi = 400
)

ggsave(
  file.path(scatter_dir, "Scatter_MFSD2A_PLVAP_By_BBB_State.svg"),
  p_bbb, width = 12, height = 10
)

print(p_bbb)

# ---- 5. Correlation within each individual mouse ----
# Each mouse must have enough paired vessels for a meaningful estimate.

mouse_correlations <- vessel_wide %>%
  group_by(Treatment, Sample_ID) %>%
  group_modify(~ {
    spearman <- get_correlation(.x, "spearman")
    pearson <- get_correlation(.x, "pearson")
    
    tibble(
      n_vessels = nrow(.x),
      spearman_rho = spearman$correlation,
      spearman_p = spearman$p_value,
      pearson_r = pearson$correlation,
      pearson_p = pearson$p_value
    )
  }) %>%
  ungroup()

write_csv(
  mouse_correlations,
  file.path(scatter_dir, "correlations_by_individual_mouse.csv")
)

print(mouse_correlations, n = Inf)

# ---- 6. Individual mouse correlations within each BBB state ----

mouse_bbb_correlations <- vessel_wide %>%
  filter(!is.na(BBB_State)) %>%
  group_by(Treatment, Sample_ID, BBB_State) %>%
  group_modify(~ {
    spearman <- get_correlation(.x, "spearman")
    pearson <- get_correlation(.x, "pearson")
    
    tibble(
      n_vessels = nrow(.x),
      spearman_rho = spearman$correlation,
      spearman_p = spearman$p_value,
      pearson_r = pearson$correlation,
      pearson_p = pearson$p_value
    )
  }) %>%
  ungroup()

write_csv(
  mouse_bbb_correlations,
  file.path(scatter_dir, "correlations_by_mouse_and_BBB_state.csv")
)

print(mouse_bbb_correlations, n = Inf)

# ---- 7. Compare mouse-level correlations between treatments ----
# Use Fisher z-transformed Spearman rho for group summaries.
# Permute treatment labels across mice, not across vessels.
# With 3 mice per group, there are only 20 possible assignments.
#
# This is exploratory and assumes independent mice and exchangeable
# treatment labels under the null hypothesis.

exact_mouse_permutation <- function(df) {
  
  df <- df %>%
    filter(
      !is.na(Treatment),
      !is.na(Sample_ID),
      !is.na(spearman_rho)
    ) %>%
    distinct(Sample_ID, Treatment, spearman_rho)
  
  if (nrow(df) != 6 ||
      sum(df$Treatment == "Saline") != 3 ||
      sum(df$Treatment == "Avastin") != 3) {
    return(tibble(
      n_Saline = sum(df$Treatment == "Saline"),
      n_Avastin = sum(df$Treatment == "Avastin"),
      mean_rho_Saline = NA_real_,
      mean_rho_Avastin = NA_real_,
      difference_mean_rho = NA_real_,
      p_exact_two_sided = NA_real_,
      n_permutations = NA_integer_
    ))
  }
  
  # Bound correlations to avoid infinite Fisher z values.
  rho <- pmax(pmin(df$spearman_rho, 0.999999), -0.999999)
  z <- atanh(rho)
  
  saline_idx <- which(df$Treatment == "Saline")
  observed <- mean(z[df$Treatment == "Avastin"]) -
    mean(z[saline_idx])
  
  allocations <- combn(seq_len(nrow(df)), 3, simplify = FALSE)
  
  perm_stats <- vapply(allocations, function(saline_set) {
    avastin_set <- setdiff(seq_len(nrow(df)), saline_set)
    mean(z[avastin_set]) - mean(z[saline_set])
  }, numeric(1))
  
  p_exact <- mean(abs(perm_stats) >= abs(observed) - 1e-12)
  
  tibble(
    n_Saline = 3L,
    n_Avastin = 3L,
    mean_rho_Saline = mean(df$spearman_rho[df$Treatment == "Saline"]),
    mean_rho_Avastin = mean(df$spearman_rho[df$Treatment == "Avastin"]),
    difference_mean_rho =
      mean(df$spearman_rho[df$Treatment == "Avastin"]) -
      mean(df$spearman_rho[df$Treatment == "Saline"]),
    difference_mean_fisher_z = observed,
    p_exact_two_sided = p_exact,
    n_permutations = length(perm_stats)
  )
}

# Overall, each mouse's correlation uses all BBB states.

mouse_comparison_overall <- exact_mouse_permutation(mouse_correlations) %>%
  mutate(analysis = "All BBB states")

# Separate mouse-level correlation comparisons for each BBB state.

mouse_comparison_bbb <- mouse_bbb_correlations %>%
  group_by(BBB_State) %>%
  group_modify(~ exact_mouse_permutation(.x)) %>%
  ungroup() %>%
  mutate(analysis = as.character(BBB_State))

mouse_comparison <- bind_rows(
  mouse_comparison_overall,
  mouse_comparison_bbb
) %>%
  mutate(
    p_exact_Holm = p.adjust(p_exact_two_sided, method = "holm")
  )

write_csv(
  mouse_comparison,
  file.path(scatter_dir, "mouse_level_correlation_treatment_comparisons.csv")
)

print(mouse_comparison, n = Inf)

# ---- 8. Plot each mouse's Spearman correlation ----

p_mouse_rho <- ggplot(
  mouse_correlations,
  aes(x = Treatment, y = spearman_rho, color = Treatment)
) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 4
  ) +
  geom_line(
    aes(group = Sample_ID),
    color = "grey65",
    linewidth = 0.6,
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c("Saline" = "#0072B2", "Avastin" = "#D55E00")
  ) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  labs(
    title = "Mouse-level MFSD2A–PLVAP Spearman correlations",
    subtitle = "Each point represents one mouse; n = 3 mice per treatment",
    x = NULL,
    y = "Spearman correlation (\u03c1)",
    color = "Treatment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  file.path(scatter_dir, "Mouse_Level_Spearman_Correlations.png"),
  p_mouse_rho, width = 7, height = 6, dpi = 400
)

ggsave(
  file.path(scatter_dir, "Mouse_Level_Spearman_Correlations.svg"),
  p_mouse_rho, width = 7, height = 6
)

print(p_mouse_rho)

# ---- 9. Summary of generated outputs ----

message("Analysis complete. Outputs saved to: ", scatter_dir)



#PLVAP variance----

# ---- PLVAP spread and variability analysis: Avastin vs Saline ----

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(purrr)

# ---- Output folder ----

spread_dir <- "Arpan Multiplex/Day5_Area_Normalized/PLVAP spread"

dir.create(
  spread_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(exists("day5_all_vessels_long"))

# ---- Identify columns ----

dat <- day5_all_vessels_long

bbb_candidates <- c(
  "BBB Stability", "BBB_State", "BBB_group", "BBB_Group"
)

bbb_col <- bbb_candidates[bbb_candidates %in% names(dat)][1]

if (length(bbb_col) == 0 || is.na(bbb_col)) {
  stop("Could not identify the BBB-state column.")
}

required_cols <- c(
  "Sample_ID", "Treatment", "Marker",
  "Vessel_ID", "Expression_per_area"
)

missing_cols <- setdiff(required_cols, names(dat))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ---- Prepare PLVAP data ----
# The analysis uses log2(area-normalized expression + 1).

plvap <- dat %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    Vessel_ID = as.character(Vessel_ID),
    Treatment = case_when(
      str_to_lower(as.character(Treatment)) %in%
        c("saline", "pbs", "vehicle") ~ "Saline",
      str_to_lower(as.character(Treatment)) == "avastin" ~ "Avastin",
      TRUE ~ as.character(Treatment)
    ),
    BBB_raw = as.character(.data[[bbb_col]]),
    BBB_State = case_when(
      str_detect(str_to_lower(BBB_raw), "low") ~ "Low BBB",
      str_detect(str_to_lower(BBB_raw), "mid|intermediate") ~ "Mid BBB",
      str_detect(str_to_lower(BBB_raw), "high") ~ "High BBB",
      TRUE ~ NA_character_
    ),
    Expression_per_area = as.numeric(Expression_per_area)
  ) %>%
  filter(
    Marker == "PLVAP",
    Treatment %in% c("Saline", "Avastin"),
    !is.na(Sample_ID),
    !is.na(Vessel_ID),
    is.finite(Expression_per_area),
    Expression_per_area >= 0
  ) %>%
  group_by(Sample_ID, Treatment, BBB_State, Vessel_ID) %>%
  summarise(
    Expression_per_area = median(
      Expression_per_area, na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("Saline", "Avastin")
    ),
    BBB_State = factor(
      BBB_State,
      levels = c("Low BBB", "Mid BBB", "High BBB")
    ),
    PLVAP_log2 = log2(Expression_per_area + 1)
  )

if (nrow(plvap) == 0) {
  stop("No eligible PLVAP measurements were found.")
}

write_csv(
  plvap,
  file.path(spread_dir, "PLVAP_area_normalized_data.csv")
)

# ---- Helper: descriptive variability statistics ----

spread_stats <- function(x) {
  
  x <- x[is.finite(x)]
  n <- length(x)
  
  tibble(
    n_vessels = n,
    mean = if (n > 0) mean(x) else NA_real_,
    median = if (n > 0) median(x) else NA_real_,
    SD = if (n > 1) sd(x) else NA_real_,
    variance = if (n > 1) var(x) else NA_real_,
    IQR = if (n > 0) IQR(x) else NA_real_,
    MAD = if (n > 0) mad(x, constant = 1) else NA_real_,
    Q25 = if (n > 0) unname(quantile(x, 0.25)) else NA_real_,
    Q75 = if (n > 0) unname(quantile(x, 0.75)) else NA_real_
  )
}

# ---- Analysis 1: pooled descriptive spread ----
# Summaries are calculated across vessels pooled within each treatment.
# These are descriptive, not mouse-level inferential tests.

pooled_summary <- plvap %>%
  group_by(Treatment) %>%
  group_modify(~ spread_stats(.x$PLVAP_log2)) %>%
  ungroup() %>%
  mutate(BBB_State = "All BBB states")

bbb_pooled_summary <- plvap %>%
  filter(!is.na(BBB_State)) %>%
  group_by(BBB_State, Treatment) %>%
  group_modify(~ spread_stats(.x$PLVAP_log2)) %>%
  ungroup()

pooled_summary_all <- bind_rows(
  pooled_summary,
  bbb_pooled_summary
) %>%
  select(BBB_State, Treatment, everything())

write_csv(
  pooled_summary_all,
  file.path(spread_dir, "PLVAP_pooled_spread_summary.csv")
)

print(pooled_summary_all, n = Inf)

# ---- Calculate variance ratios and SD ratios ----
# Ratio < 1 means lower observed variability with Avastin.

variance_ratios <- pooled_summary_all %>%
  select(
    BBB_State, Treatment,
    n_vessels, SD, variance, IQR, MAD
  ) %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(n_vessels, SD, variance, IQR, MAD),
    names_sep = "_"
  ) %>%
  mutate(
    variance_ratio_Avastin_over_Saline =
      variance_Avastin / variance_Saline,
    SD_ratio_Avastin_over_Saline =
      SD_Avastin / SD_Saline,
    IQR_ratio_Avastin_over_Saline =
      IQR_Avastin / IQR_Saline,
    MAD_ratio_Avastin_over_Saline =
      MAD_Avastin / MAD_Saline
  )

write_csv(
  variance_ratios,
  file.path(spread_dir, "PLVAP_variance_ratios.csv")
)

print(variance_ratios, n = Inf)

# ---- Analysis 2: Brown-Forsythe test ----
# This is a median-centered Levene test.
# Pooled-vessel p-values are exploratory because vessels are nested
# within mice and are not independent biological replicates.

brown_forsythe <- function(df) {
  
  df <- df %>%
    filter(
      !is.na(Treatment),
      is.finite(PLVAP_log2)
    ) %>%
    mutate(Treatment = droplevels(Treatment))
  
  n_groups <- n_distinct(df$Treatment)
  
  if (n_groups != 2 ||
      any(table(df$Treatment) < 2)) {
    return(tibble(
      n_vessels = nrow(df),
      BF_F = NA_real_,
      df1 = NA_real_,
      df2 = NA_real_,
      p_value = NA_real_
    ))
  }
  
  df <- df %>%
    group_by(Treatment) %>%
    mutate(
      absolute_deviation = abs(
        PLVAP_log2 - median(PLVAP_log2, na.rm = TRUE)
      )
    ) %>%
    ungroup()
  
  fit <- try(
    aov(absolute_deviation ~ Treatment, data = df),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(tibble(
      n_vessels = nrow(df),
      BF_F = NA_real_,
      df1 = NA_real_,
      df2 = NA_real_,
      p_value = NA_real_
    ))
  }
  
  a <- summary(fit)[[1]]
  
  tibble(
    n_vessels = nrow(df),
    BF_F = a[1, "F value"],
    df1 = a[1, "Df"],
    df2 = a[2, "Df"],
    p_value = a[1, "Pr(>F)"]
  )
}

bf_overall <- brown_forsythe(plvap) %>%
  mutate(BBB_State = "All BBB states")

bf_by_bbb <- plvap %>%
  filter(!is.na(BBB_State)) %>%
  group_by(BBB_State) %>%
  group_modify(~ brown_forsythe(.x)) %>%
  ungroup()

bf_results <- bind_rows(
  bf_overall,
  bf_by_bbb
) %>%
  mutate(
    p_Holm = p.adjust(p_value, method = "holm")
  )

write_csv(
  bf_results,
  file.path(spread_dir, "PLVAP_Brown_Forsythe_results.csv")
)

print(bf_results, n = Inf)

# ---- Analysis 3: variability within each mouse ----
# Each row is one independent mouse/sample within a BBB state.

mouse_spread <- plvap %>%
  group_by(Sample_ID, Treatment, BBB_State) %>%
  group_modify(~ spread_stats(.x$PLVAP_log2)) %>%
  ungroup()

write_csv(
  mouse_spread,
  file.path(spread_dir, "PLVAP_spread_by_mouse.csv")
)

print(mouse_spread, n = Inf)

# ---- Exact permutation test on mouse-level spread ----
# Tests the difference in mean mouse-level variability between treatments.
# Permutes treatment labels across mice, not across vessels.
# With 3 mice per treatment there are only 20 allocations.
# Exact p-values therefore have limited resolution.

exact_spread_test <- function(df, metric) {
  
  df <- df %>%
    filter(
      !is.na(Treatment),
      !is.na(Sample_ID),
      is.finite(.data[[metric]])
    ) %>%
    distinct(Sample_ID, Treatment, .data[[metric]])
  
  # Verify that each mouse belongs to only one treatment.
  mouse_check <- df %>%
    distinct(Sample_ID, Treatment) %>%
    count(Sample_ID)
  
  if (any(mouse_check$n != 1)) {
    stop("At least one Sample_ID is associated with multiple treatments.")
  }
  
  n_sal <- sum(df$Treatment == "Saline")
  n_ava <- sum(df$Treatment == "Avastin")
  
  if (nrow(df) != 6 || n_sal != 3 || n_ava != 3) {
    return(tibble(
      metric = metric,
      n_Saline = n_sal,
      n_Avastin = n_ava,
      mean_Saline = NA_real_,
      mean_Avastin = NA_real_,
      difference_Avastin_minus_Saline = NA_real_,
      ratio_Avastin_over_Saline = NA_real_,
      p_exact_two_sided = NA_real_,
      n_permutations = NA_integer_
    ))
  }
  
  x <- df[[metric]]
  treatment <- as.character(df$Treatment)
  
  observed_sal <- mean(x[treatment == "Saline"])
  observed_ava <- mean(x[treatment == "Avastin"])
  observed_diff <- observed_ava - observed_sal
  
  allocations <- combn(seq_len(nrow(df)), 3, simplify = FALSE)
  
  perm_diffs <- vapply(allocations, function(sal_idx) {
    ava_idx <- setdiff(seq_len(nrow(df)), sal_idx)
    mean(x[ava_idx]) - mean(x[sal_idx])
  }, numeric(1))
  
  p_exact <- mean(
    abs(perm_diffs) >= abs(observed_diff) - 1e-12
  )
  
  tibble(
    metric = metric,
    n_Saline = n_sal,
    n_Avastin = n_ava,
    mean_Saline = observed_sal,
    mean_Avastin = observed_ava,
    difference_Avastin_minus_Saline = observed_diff,
    ratio_Avastin_over_Saline =
      if (observed_sal != 0) observed_ava / observed_sal else NA_real_,
    p_exact_two_sided = p_exact,
    n_permutations = length(perm_diffs)
  )
}

# Test the SD, IQR, and MAD across mice.
# The test is performed for all BBB states combined and separately
# within Low, Mid, and High BBB states.


# ---- Exact mouse-level spread comparisons ----

mouse_spread_comparisons <- mouse_spread_long %>%
  group_by(BBB_State, metric) %>%
  group_modify(~ {
    result <- exact_spread_test(.x, "spread_value")
    
    # Remove metric because group_modify() adds it automatically
    result %>% select(-any_of("metric"))
  }) %>%
  ungroup() %>%
  mutate(
    p_Holm_within_metric = ave(
      p_exact_two_sided,
      metric,
      FUN = function(x) p.adjust(x, method = "holm")
    )
  )

write_csv(
  mouse_spread_comparisons,
  file.path(spread_dir, "PLVAP_mouse_level_spread_tests.csv")
)

print(mouse_spread_comparisons, n = Inf)


# ---- Plot 1: pooled vessel-level PLVAP distributions ----

p_distribution <- ggplot(
  plvap,
  aes(x = Treatment, y = PLVAP_log2, fill = Treatment)
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.25,
    color = NA
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.6
  ) +
  geom_point(
    aes(color = Sample_ID),
    position = position_jitter(width = 0.12, height = 0),
    alpha = 0.28,
    size = 1.3
  ) +
  facet_wrap(~ BBB_State, nrow = 1, drop = TRUE) +
  scale_fill_manual(
    values = c("Saline" = "#0072B2", "Avastin" = "#D55E00")
  ) +
  labs(
    title = "PLVAP distribution by treatment",
    subtitle = "Area-normalized expression; log2(x + 1)",
    x = NULL,
    y = "PLVAP expression per area [log2(x + 1)]",
    fill = "Treatment",
    color = "Mouse"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(spread_dir, "PLVAP_distribution_by_treatment.png"),
  p_distribution, width = 12, height = 6, dpi = 400
)

ggsave(
  file.path(spread_dir, "PLVAP_distribution_by_treatment.svg"),
  p_distribution, width = 12, height = 6
)

print(p_distribution)

# ---- Plot 2: mouse-level standard deviations ----

p_mouse_sd <- ggplot(
  mouse_spread %>% filter(is.finite(SD)),
  aes(x = Treatment, y = SD, color = Treatment)
) +
  geom_point(
    position = position_jitter(width = 0.07, height = 0),
    size = 4
  ) +
  geom_line(
    aes(group = Sample_ID),
    color = "grey65",
    linewidth = 0.7
  ) +
  facet_wrap(~ BBB_State, nrow = 1, drop = TRUE) +
  scale_color_manual(
    values = c("Saline" = "#0072B2", "Avastin" = "#D55E00")
  ) +
  labs(
    title = "PLVAP variability within individual mice",
    subtitle = "Each point represents one mouse; lower SD indicates a narrower distribution",
    x = NULL,
    y = "Within-mouse SD of PLVAP expression",
    color = "Treatment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  file.path(spread_dir, "PLVAP_within_mouse_SD.png"),
  p_mouse_sd, width = 11, height = 6, dpi = 400
)

ggsave(
  file.path(spread_dir, "PLVAP_within_mouse_SD.svg"),
  p_mouse_sd, width = 11, height = 6
)

print(p_mouse_sd)

# ---- Plot 3: mouse-level IQR ----

p_mouse_iqr <- ggplot(
  mouse_spread %>% filter(is.finite(IQR)),
  aes(x = Treatment, y = IQR, color = Treatment)
) +
  geom_point(
    position = position_jitter(width = 0.07, height = 0),
    size = 4
  ) +
  geom_line(
    aes(group = Sample_ID),
    color = "grey65",
    linewidth = 0.7
  ) +
  facet_wrap(~ BBB_State, nrow = 1, drop = TRUE) +
  scale_color_manual(
    values = c("Saline" = "#0072B2", "Avastin" = "#D55E00")
  ) +
  labs(
    title = "PLVAP interquartile range within individual mice",
    subtitle = "Each point represents one mouse; lower IQR indicates a narrower middle 50% of values",
    x = NULL,
    y = "Within-mouse IQR of PLVAP expression",
    color = "Treatment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  file.path(spread_dir, "PLVAP_within_mouse_IQR.png"),
  p_mouse_iqr, width = 11, height = 6, dpi = 400
)

ggsave(
  file.path(spread_dir, "PLVAP_within_mouse_IQR.svg"),
  p_mouse_iqr, width = 11, height = 6
)

print(p_mouse_iqr)

# ---- Finished ----

message("PLVAP spread analysis complete.")
message("Results saved to: ", spread_dir)


#GLUT1 MFSD2A corr ----
# ---- Day 5: GLUT1-MFSD2A correlation analysis ----

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# ---- Output directory ----
correlation_dir <- file.path(
  "Arpan Multiplex",
  "Day5_Area_Normalized",
  "Correlations_GLUT1_MFSD2A"
)

dir.create(
  correlation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ---- Confirm required columns ----
required_columns <- c(
  "Sample_ID", "Treatment", "Vessel_ID",
  "Marker", "Expression_per_area_log2"
)

missing_columns <- setdiff(
  required_columns,
  colnames(day5_all_vessels_long)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ---- Extract GLUT1 and MFSD2A measurements ----
correlation_long <- day5_all_vessels_long %>%
  mutate(
    Marker_clean = toupper(trimws(as.character(Marker))),
    Treatment = factor(
      as.character(Treatment),
      levels = c("Saline", "Avastin")
    ),
    Expression_per_area_log2 = as.numeric(
      Expression_per_area_log2
    )
  ) %>%
  filter(
    Marker_clean %in% c("GLUT1", "MFSD2A"),
    !is.na(Sample_ID),
    !is.na(Vessel_ID),
    !is.na(Treatment),
    is.finite(Expression_per_area_log2)
  ) %>%
  select(
    Sample_ID, Treatment, Vessel_ID,
    Marker_clean, Expression_per_area_log2
  )

# ---- Match both markers within each vessel ----
correlation_wide <- correlation_long %>%
  group_by(Sample_ID, Treatment, Vessel_ID, Marker_clean) %>%
  summarise(
    Expression = mean(Expression_per_area_log2),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Marker_clean,
    values_from = Expression
  ) %>%
  filter(
    is.finite(GLUT1),
    is.finite(MFSD2A)
  ) %>%
  mutate(
    Treatment = factor(
      as.character(Treatment),
      levels = c("Saline", "Avastin")
    )
  )

if (nrow(correlation_wide) < 3) {
  stop("Fewer than three matched vessels have both marker measurements.")
}

# ---- Spearman correlation function ----
calculate_spearman <- function(data) {
  
  data <- data %>%
    filter(
      is.finite(GLUT1),
      is.finite(MFSD2A)
    )
  
  n <- nrow(data)
  
  if (
    n < 3 ||
    dplyr::n_distinct(data$GLUT1) < 2 ||
    dplyr::n_distinct(data$MFSD2A) < 2
  ) {
    return(tibble(
      n_vessels = n,
      Spearman_rho = NA_real_,
      p_value = NA_real_
    ))
  }
  
  test <- suppressWarnings(
    cor.test(
      data$GLUT1,
      data$MFSD2A,
      method = "spearman",
      exact = FALSE
    )
  )
  
  tibble(
    n_vessels = n,
    Spearman_rho = unname(test$estimate),
    p_value = test$p.value
  )
}

# ---- 1. Correlation within each mouse ----
# These are the mouse-specific correlations.
mouse_correlations <- correlation_wide %>%
  group_by(Treatment, Sample_ID) %>%
  group_modify(~ calculate_spearman(.x)) %>%
  ungroup() %>%
  arrange(Treatment, Sample_ID)

print(mouse_correlations, n = Inf)

write_csv(
  mouse_correlations,
  file.path(correlation_dir, "Day5_GLUT1_MFSD2A_Spearman_By_Mouse.csv")
)

# ---- 2. Pooled vessel-level correlations by treatment ----
# Descriptive analysis: vessels from the same mouse are not independent.
treatment_correlations <- correlation_wide %>%
  group_by(Treatment) %>%
  group_modify(~ calculate_spearman(.x)) %>%
  ungroup()

print(treatment_correlations, n = Inf)

write_csv(
  treatment_correlations,
  file.path(
    correlation_dir,
    "Day5_GLUT1_MFSD2A_Spearman_By_Treatment_Descriptive.csv"
  )
)

# ---- 3. Save matched vessel-level data ----
write_csv(
  correlation_wide,
  file.path(
    correlation_dir,
    "Day5_GLUT1_MFSD2A_Matched_Vessels.csv"
  )
)

# ---- 4. Scatter plot, faceted by treatment and mouse ----
treatment_colors <- c(
  "Saline" = "#0072B2",
  "Avastin" = "#D55E00"
)

p_mouse <- ggplot(
  correlation_wide,
  aes(x = GLUT1, y = MFSD2A, color = Treatment)
) +
  geom_point(alpha = 0.65, size = 1.8) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.8
  ) +
  facet_grid(Treatment ~ Sample_ID, scales = "free") +
  scale_color_manual(values = treatment_colors) +
  labs(
    title = "Day 5 GLUT1-MFSD2A Correlation",
    subtitle = "Area-normalized, log2-transformed expression in matched vessels",
    x = "GLUT1 expression per area (log2 scale)",
    y = "MFSD2A expression per area (log2 scale)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "none",
    panel.spacing = grid::unit(1, "lines")
  )

print(p_mouse)

ggsave(
  file.path(
    correlation_dir,
    "Day5_GLUT1_MFSD2A_Correlation_By_Mouse.png"
  ),
  p_mouse,
  width = 15,
  height = 8,
  dpi = 600
)

# ---- 5. Treatment-level descriptive scatter plot ----
p_treatment <- ggplot(
  correlation_wide,
  aes(x = GLUT1, y = MFSD2A, color = Treatment)
) +
  geom_point(alpha = 0.55, size = 1.8) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1
  ) +
  facet_wrap(~ Treatment, scales = "free") +
  scale_color_manual(values = treatment_colors) +
  labs(
    title = "GLUT1-MFSD2A Relationship by Treatment",
    subtitle = "Pooled vessel-level visualization; mice are the biological replicates",
    x = "GLUT1 expression per area (log2 scale)",
    y = "MFSD2A expression per area (log2 scale)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "none"
  )

print(p_treatment)

ggsave(
  file.path(
    correlation_dir,
    "Day5_GLUT1_MFSD2A_Correlation_By_Treatment.png"
  ),
  p_treatment,
  width = 11,
  height = 6,
  dpi = 600
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    file.path(
      correlation_dir,
      "Day5_GLUT1_MFSD2A_Correlation_By_Mouse.svg"
    ),
    p_mouse,
    width = 15,
    height = 8,
    device = svglite::svglite
  )
  
  ggsave(
    file.path(
      correlation_dir,
      "Day5_GLUT1_MFSD2A_Correlation_By_Treatment.svg"
    ),
    p_treatment,
    width = 11,
    height = 6,
    device = svglite::svglite
  )
}

message("Correlation analysis complete. Results saved to: ", correlation_dir)

# ---- Add correlation labels to plots ----

library(dplyr)
library(ggplot2)

# ---- Format correlation labels ----
format_cor_label <- function(rho, p, n) {
  if (is.na(rho) || is.na(p)) {
    return(paste0("n = ", n, "\nCorrelation unavailable"))
  }
  
  p_label <- if (p < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", formatC(p, format = "f", digits = 3))
  }
  
  paste0(
    "Spearman \u03c1 = ", sprintf("%.2f", rho),
    "\n", p_label,
    "\nn = ", n, " vessels"
  )
}

# ---- Mouse-level correlation labels ----
mouse_labels <- mouse_correlations %>%
  mutate(
    Sample_ID = as.character(Sample_ID),
    label = mapply(
      format_cor_label,
      Spearman_rho,
      p_value,
      n_vessels
    )
  )

# ---- Plot: individual mouse panels ----
p_mouse <- ggplot(
  correlation_wide,
  aes(x = GLUT1, y = MFSD2A, color = Treatment)
) +
  geom_point(alpha = 0.60, size = 1.6) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.8
  ) +
  geom_label(
    data = mouse_labels,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 3.2,
    color = "black",
    fill = "white",
    label.size = 0.2
  ) +
  facet_grid(Treatment ~ Sample_ID, scales = "free") +
  scale_color_manual(
    values = c(
      "Saline" = "#0072B2",
      "Avastin" = "#D55E00"
    )
  ) +
  labs(
    title = "Day 5 GLUT1-MFSD2A Correlation",
    subtitle = "Spearman correlations calculated separately for each mouse",
    x = "GLUT1 expression per area (log2 scale)",
    y = "MFSD2A expression per area (log2 scale)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "none",
    panel.spacing = grid::unit(1, "lines")
  )

print(p_mouse)

ggsave(
  file.path(
    correlation_dir,
    "Day5_GLUT1_MFSD2A_Correlation_By_Mouse.png"
  ),
  p_mouse,
  width = 16,
  height = 8,
  dpi = 600
)

# ---- Treatment-level correlation labels ----
treatment_labels <- treatment_correlations %>%
  mutate(
    Treatment = factor(
      as.character(Treatment),
      levels = c("Saline", "Avastin")
    ),
    label = mapply(
      format_cor_label,
      Spearman_rho,
      p_value,
      n_vessels
    )
  )

# ---- Plot: pooled vessel-level correlations by treatment ----
p_treatment <- ggplot(
  correlation_wide,
  aes(x = GLUT1, y = MFSD2A, color = Treatment)
) +
  geom_point(alpha = 0.45, size = 1.5) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1
  ) +
  geom_label(
    data = treatment_labels,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.1,
    size = 3.6,
    color = "black",
    fill = "white",
    label.size = 0.2
  ) +
  facet_wrap(~ Treatment, scales = "free") +
  scale_color_manual(
    values = c(
      "Saline" = "#0072B2",
      "Avastin" = "#D55E00"
    )
  ) +
  labs(
    title = "GLUT1-MFSD2A Relationship by Treatment",
    subtitle = "Pooled vessel-level Spearman correlations",
    x = "GLUT1 expression per area (log2 scale)",
    y = "MFSD2A expression per area (log2 scale)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "none"
  )

print(p_treatment)

ggsave(
  file.path(
    correlation_dir,
    "Day5_GLUT1_MFSD2A_Correlation_By_Treatment.png"
  ),
  p_treatment,
  width = 11,
  height = 6,
  dpi = 600
)

# ---- SVG exports ----
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    file.path(
      correlation_dir,
      "Day5_GLUT1_MFSD2A_Correlation_By_Mouse.svg"
    ),
    p_mouse,
    width = 16,
    height = 8,
    device = svglite::svglite
  )
  
  ggsave(
    file.path(
      correlation_dir,
      "Day5_GLUT1_MFSD2A_Correlation_By_Treatment.svg"
    ),
    p_treatment,
    width = 11,
    height = 6,
    device = svglite::svglite
  )
}

