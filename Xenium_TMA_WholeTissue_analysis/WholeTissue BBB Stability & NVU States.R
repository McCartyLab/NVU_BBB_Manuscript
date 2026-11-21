#BBB Stability Score across datasets
#NVU States



#NVU States----
#Read Tumor Chunk Obj
vascular_subtypes_whole_tissue <- readRDS("vascular_subtypes_large_tissue.rds")

# NVU State composition per sample ----

library(dplyr)
library(ggplot2)

# Define final NVU state nomenclature
NVU_state_labels <- c(
  "Canonical NVU (no tumor)"        = "N.T.: Intact",
  "Intact NVU (+ tumor)"            = "T.I.: Intact",
  "No astrocyte/end-feet (+ tumor)" = "T.I.: Astrocyte displaced",
  "Pericyte dropout (+ tumor)"      = "T.I.: Pericyte displaced",
  "Bare vessel (+ tumor)"           = "T.I.: Endothelial only",
  "Astro, no pericyte (no tumor)"   = "N.T.: Pericyte displaced",
  "Bare EC vessel (no tumor)"       = "N.T.: Endothelial only",
  "EC/PC only (no tumor)"           = "N.T.: Astrocyte displaced"
)

vascular_subtypes_whole_tissue$NVU_States <- recode(
  vascular_subtypes_whole_tissue$long_label,
  !!!NVU_state_labels
)

# Extract classified NVU states and rename
nvu_plot_df <- vascular_subtypes_whole_tissue[[]] %>%
  dplyr::select(sample, long_label) %>%
  filter(!is.na(long_label)) %>%
  mutate(
    sample = recode(sample, "Gr_3" = "Grade 3"),
    NVU_States = recode(long_label, !!!NVU_state_labels)
  ) %>%
  dplyr::count(sample, NVU_States) %>%
  group_by(sample) %>%
  mutate(
    percentage = 100 * n / sum(n)
  ) %>%
  ungroup()

# Set NVU state order
nvu_plot_df$NVU_States <- factor(
  nvu_plot_df$NVU_States,
  levels = c(
    "N.T.: Intact",
    "N.T.: Astrocyte displaced",
    "N.T.: Pericyte displaced",
    "N.T.: Endothelial only",
    "T.I.: Intact",
    "T.I.: Astrocyte displaced",
    "T.I.: Pericyte displaced",
    "T.I.: Endothelial only"
  )
)

# Set sample order
nvu_plot_df$sample <- factor(
  nvu_plot_df$sample,
  levels = c(
    "Normal",
    "Grade 3",
    "GBM",
    "Breast",
    "Kidney",
    "Lung",
    "Melanoma",
    "Thyroid"
  )
)


# 100% stacked barplot
chunks_nvu_states_barplot <- ggplot(
  nvu_plot_df,
  aes(
    x = sample,
    y = percentage,
    fill = NVU_States
  )
) +
  geom_col(
    width = 0.8,
    color = "black",
    linewidth = 0.25
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "NVU state (%)",
    fill = "NVU state"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 13,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 13
    ),
    axis.title.y = element_text(
      size = 15
    ),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(
      size = 12
    )
  )

# Set directory where you want to save
output_dir <- "Vascular subtype analysis/Figures/"
# Make sure the directory exists
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# Export the plot as SVG
ggsave(
  filename = file.path(output_dir, "Chunks NVU States Stacked Barplot.svg"),
  plot = chunks_nvu_states_barplot,
  width = 12,      # width in inches
  height = 4,     # height in inches
  dpi = 300
)




#EC BBB Stability 

#scRNA
scRNA_ec <- readRDS("Bejarano et al. 2025 Immunity/Bejarano_seurat_modules.rds")
scRNA_ec$disease_category <- dplyr::case_when(
  scRNA_ec$Disease.subtype == "Non.Tumor" ~ "Non-tumor",
  scRNA_ec$Disease.subtype == "IDH mut glioma" ~ "LGG",
  scRNA_ec$Disease.subtype == "IDH WT glioblastoma" ~ "GBM",
  scRNA_ec$Disease.subtype == "BrM" ~ "BrM",
  TRUE ~ NA_character_
)

#Normalize Data
scRNA_ec <- NormalizeData(
  scRNA_ec,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


#TMA core
TMA_full_EC.PC_subcl.EC <- readRDS("TMA_full_EC.PC_subcl.EC.rds")
TMA_EC_sub <- subset(TMA_full_EC.PC_subcl.EC,
                     vasc_sub == "EC")



TMA_EC_sub$disease_category <- dplyr::case_when(
  TMA_EC_sub$sample == "Normal" ~ "Non-tumor",
  TMA_EC_sub$sample == "G2" ~ "LGG",
  TMA_EC_sub$sample == "G3" ~ "LGG",
  TMA_EC_sub$sample == "GBM" ~ "GBM",
  TMA_EC_sub$sample %in% c(
    "BR",
    "COL",
    "END",
    "LU_NSC",
    "MEL",
    "RCC",
    "THY"
  ) ~ "BrM",
  TRUE ~ NA_character_
)


#Xenium chunk
chunk_EC_xenium <- subset(vascular_subtypes_whole_tissue, subset = new_clusters == "EC")



chunk_EC_xenium$disease_category <- dplyr::case_when(
  chunk_EC_xenium$sample == "Normal" ~ "Non-tumor",
  chunk_EC_xenium$sample == "Gr_3" ~ "LGG",
  chunk_EC_xenium$sample == "GBM" ~ "GBM",
  chunk_EC_xenium$sample %in% c(
    "Breast",
    "Kidney",
    "Lung",
    "Melanoma",
    "Thyroid"
  ) ~ "BrM",
  TRUE ~ NA_character_
)


# BBB Stability Score by disease and dataset----
# Violin = EC-level score distribution
# Points = mean score per biological sample/core
# Datasets shown as separate facets

library(dplyr)
library(ggplot2)

# 1. Prepare EC-level data
scRNA_plot <- scRNA_ec@meta.data %>%
  filter(
    !is.na(BBB_Stability_Score),
    !is.na(disease_category),
    !is.na(Patient)
  ) %>%
  mutate(
    Dataset = "scRNA",
    Plot_ID = Patient
  ) %>%
  select(
    Dataset,
    disease_category,
    Plot_ID,
    BBB_Stability_Score
  )

TMA_plot <- TMA_EC_sub@meta.data %>%
  filter(
    !is.na(BBB_Stability_Score),
    !is.na(disease_category),
    !is.na(sample_id)
  ) %>%
  mutate(
    Dataset = "Xenium Cores",
    Plot_ID = sample_id
  ) %>%
  select(
    Dataset,
    disease_category,
    Plot_ID,
    BBB_Stability_Score
  )

chunk_plot <- chunk_EC_xenium@meta.data %>%
  filter(
    !is.na(BBB_Stability_full),
    !is.na(disease_category),
    !is.na(sample)
  ) %>%
  mutate(
    Dataset = "Xenium Chunks",
    Plot_ID = sample,
    BBB_Stability_Score = BBB_Stability_full
  ) %>%
  select(
    Dataset,
    disease_category,
    Plot_ID,
    BBB_Stability_Score
  )

# 2. Combine datasets
BBB_plot_df <- bind_rows(
  scRNA_plot,
  TMA_plot,
  chunk_plot
) %>%
  mutate(
    disease_category = factor(
      disease_category,
      levels = c("Non-tumor", "LGG", "GBM", "BrM")
    ),
    Dataset = factor(
      Dataset,
      levels = c("scRNA", "Xenium Cores", "Xenium Chunks")
    )
  )

# 3. Calculate mean BBB Stability Score for each sample/core
BBB_sample_means <- BBB_plot_df %>%
  group_by(
    Dataset,
    disease_category,
    Plot_ID
  ) %>%
  summarise(
    mean_BBB_Stability = mean(
      BBB_Stability_Score,
      na.rm = TRUE
    ),
    n_EC = n(),
    .groups = "drop"
  )

# 4. Plot
all_BBB_stability_vlnplots <- ggplot(
  BBB_plot_df,
  aes(
    x = disease_category,
    y = BBB_Stability_Score
  )
) +
  geom_violin(
    fill = "grey90",
    color = "black",
    linewidth = 0.4,
    trim = FALSE,
    scale = "width"
  ) +
  geom_point(
    data = BBB_sample_means,
    aes(
      x = disease_category,
      y = mean_BBB_Stability
    ),
    size = 3,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.8,
    position = position_jitter(
      width = 0.08,
      height = 0
    )
  ) +
  facet_wrap(
    ~ Dataset,
    nrow = 1
  ) +
  labs(
    x = NULL,
    y = "BBB Stability Score"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 13,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(
      size = 15
    ),
    strip.text = element_text(
      size = 14,
      face = "bold"
    ),
    legend.position = "none"
  )




# Set directory where you want to save
output_dir <- "Vascular subtype analysis/Figures/"
# Make sure the directory exists
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# Export the plot as SVG
ggsave(
  filename = file.path(output_dir, "CrossDataset BBB Stability Scores.svg"),
  plot = all_BBB_stability_vlnplots,
  width = 8,      # width in inches
  height = 4,     # height in inches
  dpi = 300
)



# BBB Stability Score by Xenium Chunk sample----
# High = top 25% within sample
# Low  = bottom 25% within sample
# Dotted lines = sample-specific Q25 and Q75 cutoffs

library(dplyr)
library(ggplot2)

# Prepare plotting data
BBB_chunk_plot <- chunk_EC_xenium@meta.data %>%
  filter(
    !is.na(BBB_Stability_full),
    !is.na(BBB_EC_status_within_sample),
    !is.na(sample)
  ) %>%
  mutate(
    # Rename grade samples
    sample = recode(
      sample,
      "Gr_2" = "Grade 2",
      "Gr_3" = "Grade 3"
    ),
    
    BBB_EC_status_within_sample = factor(
      BBB_EC_status_within_sample,
      levels = c("High", "Mid", "Low")
    )
  )

# Calculate sample-specific 25th and 75th percentile cutoffs
BBB_chunk_cutoffs <- BBB_chunk_plot %>%
  group_by(sample) %>%
  summarise(
    Q25 = quantile(BBB_Stability_full, 0.25, na.rm = TRUE),
    Q75 = quantile(BBB_Stability_full, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Make sample order
sample_order <- c(
  "Normal",
  "Grade 3",
  "GBM",
  "Breast",
  "Kidney",
  "Lung",
  "Melanoma",
  "Thyroid"
)

BBB_chunk_plot <- BBB_chunk_plot %>%
  mutate(
    sample = factor(sample, levels = sample_order)
  )

BBB_chunk_cutoffs <- BBB_chunk_cutoffs %>%
  mutate(
    sample = factor(sample, levels = sample_order)
  )

# Plot
chunk_bbb_stability_sample <- ggplot(
  BBB_chunk_plot,
  aes(
    x = sample,
    y = BBB_Stability_full,
    fill = BBB_EC_status_within_sample
  )
) +
  geom_violin(
    scale = "width",
    trim = FALSE,
    color = "black",
    linewidth = 0.3
  ) +
  
  # Q25 cutoff
  geom_segment(
    data = BBB_chunk_cutoffs,
    aes(
      x = as.numeric(sample) - 0.38,
      xend = as.numeric(sample) + 0.38,
      y = Q25,
      yend = Q25
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  
  # Q75 cutoff
  geom_segment(
    data = BBB_chunk_cutoffs,
    aes(
      x = as.numeric(sample) - 0.38,
      xend = as.numeric(sample) + 0.38,
      y = Q75,
      yend = Q75
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  
  scale_fill_manual(
    values = c(
      "Low"  = "#D73027",   # red
      "Mid"  = "grey80",
      "High" = "#4575B4"    # blue
    ),
    labels = c(
      "High" = "High (top 25%)",
      "Mid"  = "Intermediate (25–75%)",
      "Low"  = "Low (bottom 25%)"
    ),
    name = "BBB Stability"
  ) +
  
  labs(
    x = NULL,
    y = "BBB Stability Score"
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 13,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(
      size = 13,
      face = "bold"
    ),
    legend.text = element_text(size = 11)
  )

# Set directory where you want to save
output_dir <- "Vascular subtype analysis/Figures/"
# Make sure the directory exists
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# Export the plot as SVG
ggsave(
  filename = file.path(output_dir, "Chunk BBB Stability_sample.svg"),
  plot = chunk_bbb_stability_sample,
  width = 8,      # width in inches
  height = 4,     # height in inches
  dpi = 300
)


# BBB Stability Score by scRNA-seq patient ----
# High = top 25% within sample
# Low  = bottom 25% within sample
# Dotted lines = sample-specific Q25 and Q75 cutoffs

library(dplyr)
library(ggplot2)

# Prepare plotting data
BBB_scRNA_plot <- scRNA_ec@meta.data %>%
  filter(
    !is.na(BBB_Stability_Score),
    !is.na(BBB_group),
    !is.na(Patient)
  ) %>%
  mutate(
    # Display labels only
    Patient_label = case_when(
      Patient == "NT-38" ~ "Non-Tumor",
      grepl("IDH_WT", Patient) ~ "GBM",
      grepl("IDH_mut", Patient) ~ "Grade 3",
      grepl("BrM", Patient) ~ "BrM",
      TRUE ~ Patient
    ),
    
    BBB_group = as.character(BBB_group),
    BBB_group = recode(
      BBB_group,
      "High_BBB_Stability" = "High",
      "Intermediate" = "Mid",
      "Low_BBB_Stability" = "Low"
    ),
    BBB_group = factor(
      BBB_group,
      levels = c("High", "Mid", "Low")
    )
  )

# Calculate patient-specific 25th and 75th percentile cutoffs
# using the ORIGINAL Patient IDs
BBB_scRNA_cutoffs <- BBB_scRNA_plot %>%
  group_by(Patient, Patient_label) %>%
  summarise(
    Q25 = quantile(BBB_Stability_Score, 0.25, na.rm = TRUE),
    Q75 = quantile(BBB_Stability_Score, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Factor display labels in desired order
scRNA_order <- c(
  "Non-Tumor",
  "Grade 3",
  "GBM",
  "BrM"
)

BBB_scRNA_plot <- BBB_scRNA_plot %>%
  mutate(
    Patient_label = factor(
      Patient_label,
      levels = scRNA_order
    )
  )

BBB_scRNA_cutoffs <- BBB_scRNA_cutoffs %>%
  mutate(
    Patient_label = factor(
      Patient_label,
      levels = scRNA_order
    )
  )

# Plot
scRNA_bbb_stability_patient <- ggplot(
  BBB_scRNA_plot,
  aes(
    x = Patient_label,
    y = BBB_Stability_Score,
    fill = BBB_group
  )
) +
  geom_violin(
    scale = "width",
    trim = FALSE,
    color = "black",
    linewidth = 0.3
  ) +
  
  # Q25 cutoff
  geom_segment(
    data = BBB_scRNA_cutoffs,
    aes(
      x = as.numeric(Patient_label) - 0.38,
      xend = as.numeric(Patient_label) + 0.38,
      y = Q25,
      yend = Q25
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  
  # Q75 cutoff
  geom_segment(
    data = BBB_scRNA_cutoffs,
    aes(
      x = as.numeric(Patient_label) - 0.38,
      xend = as.numeric(Patient_label) + 0.38,
      y = Q75,
      yend = Q75
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  
  scale_fill_manual(
    values = c(
      "Low"  = "#D73027",
      "Mid"  = "grey80",
      "High" = "#4575B4"
    ),
    labels = c(
      "High" = "High (top 25%)",
      "Mid"  = "Intermediate (25–75%)",
      "Low"  = "Low (bottom 25%)"
    ),
    breaks = c("High", "Mid", "Low"),
    name = "BBB Stability"
  ) +
  
  labs(
    x = NULL,
    y = "BBB Stability Score"
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 13,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(
      size = 13,
      face = "bold"
    ),
    legend.text = element_text(size = 11)
  )

# Save
output_dir <- "Vascular subtype analysis/Figures/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  filename = file.path(output_dir, "scRNA BBB Stability_patient.svg"),
  plot = scRNA_bbb_stability_patient,
  width = 8,
  height = 4,
  dpi = 300
)


# BBB Stability Score by TMA Tumor Core ----
# High = top 25% within sample
# Low  = bottom 25% within sample
# Dotted lines = sample-specific Q25 and Q75 cutoffs

# Prepare plotting data
BBB_TMA_plot <- TMA_EC_sub@meta.data %>%
  filter(
    !is.na(BBB_Stability_Score),
    !is.na(BBB_group),
    !is.na(sample)
  ) %>%
  mutate(
    # Display labels only
    sample_label = recode(
      sample,
      "Normal" = "Non-Tumor",
      "G2" = "Grade 2",
      "G3" = "Grade 3",
      "LU_NSC" = "NSCLC"
    ),
    
    BBB_group = as.character(BBB_group),
    BBB_group = recode(
      BBB_group,
      "High_BBB_Stability" = "High",
      "Intermediate" = "Mid",
      "Low_BBB_Stability" = "Low"
    ),
    BBB_group = factor(
      BBB_group,
      levels = c("High", "Mid", "Low")
    )
  )

# Calculate core-specific 25th and 75th percentile cutoffs
# using the ORIGINAL sample IDs
BBB_TMA_cutoffs <- BBB_TMA_plot %>%
  group_by(sample, sample_label) %>%
  summarise(
    Q25 = quantile(BBB_Stability_Score, 0.25, na.rm = TRUE),
    Q75 = quantile(BBB_Stability_Score, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Factor display labels in desired order
TMA_order <- c(
  "Non-Tumor",
  "Grade 2",
  "Grade 3",
  "GBM",
  "RCC",
  "THY",
  "BR",
  "NSCLC",
  "END",
  "MEL",
  "COL"
)

BBB_TMA_plot <- BBB_TMA_plot %>%
  mutate(
    sample_label = factor(
      sample_label,
      levels = TMA_order
    )
  )

BBB_TMA_cutoffs <- BBB_TMA_cutoffs %>%
  mutate(
    sample_label = factor(
      sample_label,
      levels = TMA_order
    )
  )

# Plot
TMA_bbb_stability_core <- ggplot(
  BBB_TMA_plot,
  aes(
    x = sample_label,
    y = BBB_Stability_Score,
    fill = BBB_group
  )
) +
  geom_violin(
    scale = "width",
    trim = FALSE,
    color = "black",
    linewidth = 0.3
  ) +
  
  # Q25 cutoff
  geom_segment(
    data = BBB_TMA_cutoffs,
    aes(
      x = as.numeric(sample_label) - 0.38,
      xend = as.numeric(sample_label) + 0.38,
      y = Q25,
      yend = Q25
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  
  # Q75 cutoff
  geom_segment(
    data = BBB_TMA_cutoffs,
    aes(
      x = as.numeric(sample_label) - 0.38,
      xend = as.numeric(sample_label) + 0.38,
      y = Q75,
      yend = Q75
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  
  scale_fill_manual(
    values = c(
      "Low"  = "#D73027",
      "Mid"  = "grey80",
      "High" = "#4575B4"
    ),
    labels = c(
      "High" = "High (top 25%)",
      "Mid"  = "Intermediate (25–75%)",
      "Low"  = "Low (bottom 25%)"
    ),
    breaks = c("High", "Mid", "Low"),
    name = "BBB Stability"
  ) +
  
  labs(
    x = NULL,
    y = "BBB Stability Score"
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 13,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(
      size = 13,
      face = "bold"
    ),
    legend.text = element_text(size = 11)
  )

# Save
ggsave(
  filename = file.path(output_dir, "TMA BBB Stability_core.svg"),
  plot = TMA_bbb_stability_core,
  width = 8,
  height = 4,
  dpi = 300
)







# DotPlot of Leading Edge Genes from scRNA in our Xenium Chunk----

# Genes to plot
BBB_genes <- c(
  "SLC2A1",
  "MFSD2A",
  "ABCG2",
  "ABCB1",
  "SLC7A5",
  "SLC16A1",
  "TFRC",
  "SLCO1A2",
  "ATP10A",
  "CLDN5",
  "TJP1",
  "LSR",
  "F11R",
  "OCLN",
  "TJP2",
  "ESAM",
  "CTNNB1",
  "PECAM1",
  "SLC1A1",
  "SLC2A3",
  "SLC38A3",
  "ANXA1",
  "CAV1",
  "SLC38A2",
  "TGFB2",
  "CTSS",
  "CTSH",
  "IFITM3"
)

# Check BBB status distribution
table(
  vascular_subtypes_whole_tissue$BBB_EC_status_within_sample,
  useNA = "ifany"
)

# Keep only High and Low BBB-status ECs
vascular_BBB <- subset(
  vascular_subtypes_whole_tissue,
  subset = BBB_EC_status_within_sample %in% c("High", "Low")
)

# Set desired order: Low → High on x-axis
vascular_BBB$BBB_EC_status_within_sample <- factor(
  vascular_BBB$BBB_EC_status_within_sample,
  levels = c("Low", "High")
)

# Check which genes are present
BBB_genes_present <- BBB_genes[
  BBB_genes %in% rownames(vascular_BBB)
]

BBB_genes_missing <- BBB_genes[
  !BBB_genes %in% rownames(vascular_BBB)
]

if (length(BBB_genes_missing) > 0) {
  cat("Missing genes:\n")
  print(BBB_genes_missing)
}

# Order genes by overall average expression across High + Low BBB-status ECs

# Calculate mean expression for each gene
expr_data <- FetchData(
  vascular_BBB,
  vars = BBB_genes_present
)

gene_mean_expression <- colMeans(
  expr_data,
  na.rm = TRUE
)

# Order genes from highest to lowest expression
BBB_genes_ordered <- names(
  sort(gene_mean_expression, decreasing = TRUE)
)

# Make the gene order explicit
BBB_status_dotplot <- DotPlot(
  vascular_BBB,
  group.by = "BBB_EC_status_within_sample",
  features = BBB_genes_ordered
) +
  coord_flip() +
  labs(
    title = "High BBB Stability Associated \n Genes: Xenium Chunk",
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

BBB_status_dotplot

ggsave(
  filename = file.path(
    output_dir,
    "Leading Edge Dotplot_Chunk.svg"
  ),
  plot = BBB_status_dotplot,
  width = 6,
  height = 8,
  dpi = 300
)



#EC BBB Stability per NVU State----
# EC BBB Stability Status Across NVU States
# High vs Low BBB stability

library(dplyr)
library(ggplot2)

# Define NVU state order
NVU_order <- c(
  "N.T.: Intact",
  "N.T.: Astrocyte displaced",
  "N.T.: Pericyte displaced",
  "N.T.: Endothelial only",
  "T.I.: Intact",
  "T.I.: Astrocyte displaced",
  "T.I.: Pericyte displaced",
  "T.I.: Endothelial only"
)

# Extract metadata
NVU_BBB_df <- vascular_subtypes_whole_tissue[[]] %>%
  filter(
    !is.na(NVU_States),
    BBB_EC_status_within_sample %in% c("High", "Low")
  ) %>%
  mutate(
    NVU_States = factor(NVU_States, levels = NVU_order),
    BBB_EC_status_within_sample = factor(
      BBB_EC_status_within_sample,
      levels = c("High", "Low")
    ),
    NVU_group = case_when(
      grepl("^N\\.T\\.:", NVU_States) ~ "N.T.",
      grepl("^T\\.I\\.:", NVU_States) ~ "T.I."
    ),
    NVU_state_clean = sub(
      "^(N\\.T\\.:|T\\.I\\.:) ",
      "",
      as.character(NVU_States)
    )
  )

# Summarize High/Low proportions within each NVU state
NVU_BBB_summary <- NVU_BBB_df %>%
  count(
    NVU_group,
    NVU_States,
    NVU_state_clean,
    BBB_EC_status_within_sample,
    name = "n"
  ) %>%
  group_by(NVU_States) %>%
  mutate(
    total_high_low = sum(n),
    percent = 100 * n / total_high_low
  ) %>%
  ungroup()

# Set state order within each facet
NVU_BBB_summary$NVU_state_clean <- factor(
  NVU_BBB_summary$NVU_state_clean,
  levels = c(
    "Intact",
    "Astrocyte displaced",
    "Pericyte displaced",
    "Endothelial only"
  )
)

# Plot
NVU_BBB_dotplot <- ggplot(
  NVU_BBB_summary,
  aes(
    x = percent,
    y = NVU_state_clean,
    size = n,
    group = BBB_EC_status_within_sample
  )
) +
  geom_line(
    aes(group = NVU_state_clean),
    linewidth = 0.6,
    color = "grey75"
  ) +
  geom_point(
    aes(
      shape = BBB_EC_status_within_sample,
      color = BBB_EC_status_within_sample
    ),
    stroke = 1
  ) +
  facet_grid(
    NVU_group ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_color_manual(
    values = c(
      "High" = "#4575B4",
      "Low" = "#D73027"
    )
  ) +
  scale_shape_manual(
    values = c(
      "High" = 16,
      "Low" = 17
    )
  ) +
  scale_size_continuous(
    range = c(3, 9),
    name = "ECs"
  ) +
  labs(
    title = "EC BBB Stability Across NVU States",
    x = "ECs with High or Low BBB Stability (%)",
    y = NULL,
    color = "BBB stability",
    shape = "BBB stability"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    legend.text = element_text(size = 10),
    panel.spacing = unit(1, "lines")
  )

NVU_BBB_dotplot

# Save
ggsave(
  "Vascular subtype analysis/Figures/EC_BBB_Stability_Across_NVU_States.svg",
  NVU_BBB_dotplot,
  width = 8,
  height = 6,
  device = "svg"
)

xf
# EC BBB Stability Across NVU States----

#BY Condition
# Create simplified condition labels
vascular_subtypes_whole_tissue$new_conditons <- dplyr::case_when(
  vascular_subtypes_whole_tissue$condition == "Normal" ~ "Non-Tumor",
  vascular_subtypes_whole_tissue$condition == "Primary Tumor" ~ "Gliomas",
  vascular_subtypes_whole_tissue$condition == "Secondary Tumor" ~ "BrM",
  TRUE ~ NA_character_
)

# Set factor order
vascular_subtypes_whole_tissue$new_conditons <- factor(
  vascular_subtypes_whole_tissue$new_conditons,
  levels = c("Non-Tumor", "Gliomas", "BrM")
)

# Check
table(
  vascular_subtypes_whole_tissue$new_conditons,
  useNA = "ifany"
)

# High vs Low BBB stability
# Faceted by condition AND tumor-infiltrated status

library(dplyr)
library(ggplot2)

# Extract metadata
NVU_BBB_df <- vascular_subtypes_whole_tissue[[]] %>%
  filter(
    !is.na(NVU_States),
    !is.na(BBB_EC_status_within_sample),
    new_conditons %in% c("Non-Tumor", "Gliomas", "BrM")
  ) %>%
  mutate(
    # Separate Non-Tumor vs Tumor-Infiltrated
    NVU_context = case_when(
      grepl("^N\\.T\\.:", NVU_States) ~ "Non-Tumor",
      grepl("^T\\.I\\.:", NVU_States) ~ "Tumor-Infiltrated",
      TRUE ~ NA_character_
    ),
    
    # Collapse N.T. and T.I. into the same four biological states
    NVU_state_clean = case_when(
      NVU_States %in% c(
        "N.T.: Intact",
        "T.I.: Intact"
      ) ~ "Intact",
      
      NVU_States %in% c(
        "N.T.: Astrocyte displaced",
        "T.I.: Astrocyte displaced"
      ) ~ "Astrocyte displaced",
      
      NVU_States %in% c(
        "N.T.: Pericyte displaced",
        "T.I.: Pericyte displaced"
      ) ~ "Pericyte displaced",
      
      NVU_States %in% c(
        "N.T.: Endothelial only",
        "T.I.: Endothelial only"
      ) ~ "Endothelial only",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(NVU_context),
    !is.na(NVU_state_clean)
  ) %>%
  mutate(
    NVU_context = factor(
      NVU_context,
      levels = c(
        "Non-Tumor",
        "Tumor-Infiltrated"
      )
    ),
    
    NVU_state_clean = factor(
      NVU_state_clean,
      levels = rev(c(
        "Intact",
        "Astrocyte displaced",
        "Pericyte displaced",
        "Endothelial only"
      ))
    ),
    
    BBB_EC_status_within_sample = factor(
      BBB_EC_status_within_sample,
      levels = c("High", "Mid", "Low")
    ),
    
    new_conditons = factor(
      new_conditons,
      levels = c(
        "Non-Tumor",
        "Gliomas",
        "BrM"
      )
    )
  )

# Summarize High / Mid / Low
# Mid remains in the denominator but is not plotted
NVU_BBB_summary <- NVU_BBB_df %>%
  count(
    new_conditons,
    NVU_context,
    NVU_state_clean,
    BBB_EC_status_within_sample,
    name = "n"
  ) %>%
  group_by(
    new_conditons,
    NVU_context,
    NVU_state_clean
  ) %>%
  mutate(
    total_ECs = sum(n),
    percent = 100 * n / total_ECs
  ) %>%
  ungroup() %>%
  filter(
    BBB_EC_status_within_sample %in% c("High", "Low")
  )

# Plot
NVU_BBB_dotplot <- ggplot(
  NVU_BBB_summary,
  aes(
    x = percent,
    y = NVU_state_clean,
    color = BBB_EC_status_within_sample,
    shape = BBB_EC_status_within_sample
  )
) +
  geom_line(
    aes(
      group = interaction(
        new_conditons,
        NVU_context,
        NVU_state_clean
      )
    ),
    color = "grey50",
    linewidth = 0.8
  ) +
  geom_point(
    size = 4.5,
    stroke = 1
  ) +
  facet_grid(
    new_conditons ~ NVU_context
  ) +
  scale_color_manual(
    values = c(
      "High" = "#4575B4",
      "Low" = "#D73027"
    ),
    labels = c(
      "High" = "High BBB stability",
      "Low" = "Low BBB stability"
    )
  ) +
  scale_shape_manual(
    values = c(
      "High" = 16,
      "Low" = 17
    ),
    labels = c(
      "High" = "High BBB stability",
      "Low" = "Low BBB stability"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "EC BBB Stability Across NVU States",
    x = "ECs with High or Low BBB Stability (%)",
    y = NULL,
    color = "BBB stability",
    shape = "BBB stability"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    legend.text = element_text(size = 10),
    panel.spacing = unit(1, "lines")
  )

NVU_BBB_dotplot

ggsave(
  "Vascular subtype analysis/Figures/EC_BBB_Stability_Across_NVU_States.svg",
  NVU_BBB_dotplot,
  width = 10,
  height = 7,
  device = "svg"
)

#Astrocyte displacement is associated with a marked shift toward low endothelial BBB stability in gliomas, whereas other forms of NVU disruption show more context-dependent relationships with endothelial barrier state.




#Density of Cell types per NVU State
whole_tissue <- readRDS("NVU_BBB_annotated_Seurat_scripts01to07.rds")


# Create new_clusters from the existing clusters column
whole_tissue$new_clusters <- as.character(whole_tissue$clusters)

# Split EC/PC cells according to vascular_subtype
whole_tissue$new_clusters[
  whole_tissue$clusters == "EC/PC" &
    !is.na(whole_tissue$vascular_subtype)
] <- whole_tissue$vascular_subtype[
  whole_tissue$clusters == "EC/PC" &
    !is.na(whole_tissue$vascular_subtype)
]

# Convert back to factor if desired
whole_tissue$new_clusters <- factor(whole_tissue$new_clusters)

# Check the result
table(whole_tissue$new_clusters, useNA = "ifany")

#Subset
# # Subset to EC, Pericyte, Tumor, and Ast-like cells
# whole_tissue_NVU <- subset(
#   whole_tissue,
#   subset = new_clusters %in% c("EC", "Pericyte", "Tumor", "Ast-like")
# )
# 
# # Check the result
# table(whole_tissue_NVU$new_clusters, useNA = "ifany")
# 
# # Check object size
#object.size(whole_tissue_NVU)

# Define final NVU state nomenclature
NVU_state_labels <- c(
  "Canonical NVU (no tumor)"        = "N.T.: Intact",
  "Intact NVU (+ tumor)"            = "T.I.: Intact",
  "No astrocyte/end-feet (+ tumor)" = "T.I.: Astrocyte displaced",
  "Pericyte dropout (+ tumor)"      = "T.I.: Pericyte displaced",
  "Bare vessel (+ tumor)"           = "T.I.: Endothelial only",
  "Astro, no pericyte (no tumor)"   = "N.T.: Pericyte displaced",
  "Bare EC vessel (no tumor)"       = "N.T.: Endothelial only",
  "EC/PC only (no tumor)"           = "N.T.: Astrocyte displaced"
)

# Create new NVU state column
whole_tissue$NVU_States <- dplyr::recode(
  whole_tissue$long_label,
  !!!NVU_state_labels
)

# Set final NVU state order
whole_tissue$NVU_States <- factor(
  whole_tissue$NVU_States,
  levels = c(
    "N.T.: Intact",
    "N.T.: Astrocyte displaced",
    "N.T.: Pericyte displaced",
    "N.T.: Endothelial only",
    "T.I.: Intact",
    "T.I.: Astrocyte displaced",
    "T.I.: Pericyte displaced",
    "T.I.: Endothelial only"
  )
)

# Check
table(whole_tissue$NVU_States, useNA = "ifany")


#metadata of non NAs
# Keep only cells with a defined NVU state
whole_tissue_NVU_meta <- whole_tissue@meta.data %>%
  dplyr::filter(!is.na(NVU_States))

# Check remaining cells
table(whole_tissue_NVU_meta$NVU_States, useNA = "ifany")

# Check cell types within each NVU state
table(
  whole_tissue_NVU_meta$NVU_States,
  whole_tissue_NVU_meta$new_clusters,
  useNA = "ifany"
)


#Subset metadata to just tumor infiltrated and nvu cell types + tumor

# Define the NVU states and cell types for the analysis
TI_states <- c(
  "T.I.: Intact",
  "T.I.: Astrocyte displaced",
  "T.I.: Pericyte displaced",
  "T.I.: Endothelial only"
)

NVU_celltypes <- c(
  "EC",
  "Pericyte",
  "Ast-like",
  "Tumor"
)

# Subset metadata to Tumor-Infiltrated NVU states and relevant cell types
whole_tissue_TI_meta <- whole_tissue_NVU_meta %>%
  dplyr::filter(
    NVU_States %in% TI_states,
    new_clusters %in% NVU_celltypes
  )

# Set ordering
whole_tissue_TI_meta$NVU_States <- factor(
  whole_tissue_TI_meta$NVU_States,
  levels = TI_states
)

whole_tissue_TI_meta$new_clusters <- factor(
  whole_tissue_TI_meta$new_clusters,
  levels = c("EC", "Pericyte", "Ast-like", "Tumor")
)

# Check
table(
  whole_tissue_TI_meta$NVU_States,
  whole_tissue_TI_meta$new_clusters
)



#Cell type per NVU State----
library(dplyr)
library(ggplot2)

# Define the NVU states to show
NVU_states_plot <- c(
  "N.T.: Intact",
  "T.I.: Intact",
  "T.I.: Astrocyte displaced",
  "T.I.: Pericyte displaced",
  "T.I.: Endothelial only"
)

# Define cell types
NVU_celltypes <- c(
  "EC",
  "Pericyte",
  "Ast-like",
  "Tumor"
)

# Create metadata subset
NVU_composition <- whole_tissue_NVU_meta %>%
  filter(
    NVU_States %in% NVU_states_plot,
    new_clusters %in% NVU_celltypes
  ) %>%
  count(NVU_States, new_clusters, name = "n") %>%
  group_by(NVU_States) %>%
  mutate(
    percentage = 100 * n / sum(n)
  ) %>%
  ungroup()

# Add zero-count combinations so missing cell types are explicitly shown
NVU_composition <- NVU_composition %>%
  tidyr::complete(
    NVU_States = NVU_states_plot,
    new_clusters = NVU_celltypes,
    fill = list(n = 0, percentage = 0)
  )

# Set ordering
NVU_composition$NVU_States <- factor(
  NVU_composition$NVU_States,
  levels = NVU_states_plot
)

NVU_composition$new_clusters <- factor(
  NVU_composition$new_clusters,
  levels = c("EC", "Pericyte", "Ast-like", "Tumor"),
  labels = c("Endothelial", "Pericyte", "Astrocyte", "Tumor")
)

# Dot plot
NVU_cell_composition_dotplot <- ggplot(
  NVU_composition,
  aes(
    x = NVU_States,
    y = new_clusters,
    size = n,
    color = percentage
  )
) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(
    name = "Cell count",
    range = c(2, 12)
  ) +
  scale_color_gradient(
    name = "Composition (%)",
    low = "grey90",
    high = "#2166AC"
  ) +
  labs(
    title = "Cellular Composition Across NVU States",
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      size = 16,
      angle = 30,
      hjust = 1
    ),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

NVU_cell_composition_dotplot

# Save
ggsave(
  "Vascular subtype analysis/Figures/NVU_Cell_Composition_DotPlot.svg",
  NVU_cell_composition_dotplot,
  width = 8,
  height = 4
)






