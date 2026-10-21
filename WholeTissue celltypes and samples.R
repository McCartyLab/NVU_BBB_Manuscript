#Xenium Whole-tissue 

#Sample composition



#Whole-Tissue----

whole_tissue <- readRDS("NVU_BBB_annotated_Seurat_scripts01to07.rds")

library(dplyr)

whole_tissue@meta.data <- whole_tissue@meta.data %>%
  mutate(
    new_clusters = if_else(
      clusters == "EC/PC" & !is.na(vascular_subtype),
      vascular_subtype,
      clusters
    )
  )

whole_tissue$new_clusters_condition <- paste0(whole_tissue$new_clusters,"_",whole_tissue$condition)


# UMAP - Large Tissue----

# Define cluster order
desired_order <- c(
  "EC",
  "Pericyte",
  "Ast-like",
  "Fibroblasts",
  "OPC/OL-like",
  "Immune",
  "VLMC",
  "Neu-like",
  "Tumor"
)

# Display labels for the figure
display_labels <- c(
  "EC"          = "Endothelial",
  "Pericyte"    = "Pericyte",
  "Ast-like"    = "Ast-like",
  "Fibroblasts" = "Fibroblasts",
  "OPC/OL-like" = "OPC/OL-like",
  "Immune"      = "Immune",
  "VLMC"        = "VLMC",
  "Neu-like"    = "Neu-like",
  "Tumor"       = "Tumor"
)

# Custom colors
cluster_colors <- c(
  "EC"          = "#e31a1c",  # red
  "Pericyte"    = "#fb6a4a",  # coral
  "Ast-like"    = "#228b22",  # dark green
  "Fibroblasts" = "#a6761d",  # brown
  "OPC/OL-like" = "#e69f00",  # orange
  "Immune"      = "#6a3d9a",  # purple
  "VLMC"        = "#00a6a6",  # teal
  "Neu-like"    = "#1f78b4",  # blue
  "Tumor"       = "gray50"    # gray
)

# Set identities
Idents(whole_tissue) <- factor(
  whole_tissue$new_clusters,
  levels = c(desired_order, "LowQual_drop")
)

# Get cells to plot without creating a new Seurat object
umap_cells <- rownames(whole_tissue[[]])[
  whole_tissue$new_clusters %in% desired_order
]

# Plot UMAP
ChunkUMAP <- DimPlot(
  whole_tissue,
  reduction = "umap",
  cells = umap_cells,
  cols = cluster_colors,
  label = FALSE
) +
  scale_color_manual(
    values = cluster_colors,
    labels = display_labels,
    breaks = desired_order
  ) +
  theme_santi.umap() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 16)
    )
  )

# Set directory where you want to save
output_dir <- "Vascular subtype analysis/Figures/"

# Make sure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Export UMAP as SVG
ggsave(
  filename = file.path(output_dir, "UMAP: Helv_Tumor Chunk.tiff"),
  plot = ChunkUMAP,
  width = 10,
  height = 9,
  dpi = 300
)


#Thyroid sample UMAP
# UMAP - Thyroid Sample Only
# Large Tissue

# Define cluster order
desired_order <- c(
  "EC",
  "Pericyte",
  "Ast-like",
  "Fibroblasts",
  "OPC/OL-like",
  "Immune",
  "VLMC",
  "Neu-like",
  "Tumor"
)

# Display labels
display_labels <- c(
  "EC"          = "Endothelial",
  "Pericyte"    = "Pericyte",
  "Ast-like"    = "Ast-like",
  "Fibroblasts" = "Fibroblasts",
  "OPC/OL-like" = "OPC/OL-like",
  "Immune"      = "Immune",
  "VLMC"        = "VLMC",
  "Neu-like"    = "Neu-like",
  "Tumor"       = "Tumor"
)

# Custom colors
cluster_colors <- c(
  "EC"          = "#e31a1c",
  "Pericyte"    = "#fb6a4a",
  "Ast-like"    = "#228b22",
  "Fibroblasts" = "#a6761d",
  "OPC/OL-like" = "#e69f00",
  "Immune"      = "#6a3d9a",
  "VLMC"        = "#00a6a6",
  "Neu-like"    = "#1f78b4",
  "Tumor"       = "gray50"
)

# Set identities
Idents(whole_tissue) <- factor(
  whole_tissue$new_clusters,
  levels = c(desired_order, "LowQual_drop")
)

# Get metadata vectors
sample_vec <- as.character(whole_tissue$sample)
cluster_vec <- as.character(whole_tissue$new_clusters)

# Get ONLY Thyroid cells with valid cell-type annotations
umap_cells <- rownames(whole_tissue[[]])[
  sample_vec == "Thyroid" &
    cluster_vec %in% desired_order
]

# Plot UMAP
ThyroidUMAP <- DimPlot(
  whole_tissue,
  reduction = "umap",
  cells = umap_cells,
  cols = cluster_colors,
  label = FALSE
) +
  scale_color_manual(
    values = cluster_colors,
    labels = display_labels,
    breaks = desired_order,
    drop = FALSE
  ) +
  theme_santi.umap() +
  guides(
    color = guide_legend(
      override.aes = list(size = 16)
    )
  )

# Output directory
output_dir <- "Vascular subtype analysis/Figures/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Export
ggsave(
  filename = file.path(
    output_dir,
    "UMAP: Thyroid.svg"
  ),
  plot = ThyroidUMAP,
  width = 9,
  height = 7,
  dpi = 300
)

# Display
ThyroidUMAP


# Whole Tissue: Cell Type Proportions by Sample
# Each sample sums to 100%

library(ggplot2)
library(scales)

# Extract only the metadata vectors we need
sample_vec <- as.character(whole_tissue$sample)
cluster_vec <- as.character(whole_tissue$new_clusters)

# Explicitly remove unwanted cells
keep <- (
  !is.na(sample_vec) &
    sample_vec != "" &
    !is.na(cluster_vec) &
    cluster_vec != "" &
    cluster_vec != "LowQual_drop"
)

# Build count table
celltype_counts <- as.data.frame(
  table(
    sample = sample_vec[keep],
    new_clusters = cluster_vec[keep],
    useNA = "no"
  ),
  stringsAsFactors = FALSE
)

# Remove zero-count combinations
celltype_counts <- celltype_counts[
  celltype_counts$Freq > 0,
  ,
  drop = FALSE
]

# Rename Grade 3
celltype_counts$sample[
  celltype_counts$sample == "Gr_3"
] <- "Grade 3"

# Sample order
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

# Whole Tissue cell-type order
celltype_order <- c(
  "EC",
  "Pericyte",
  "Ast-like",
  "Fibroblasts",
  "OPC/OL-like",
  "Immune",
  "VLMC",
  "Neu-like",
  "Tumor"
)

# Keep only expected samples and cell types
celltype_counts <- celltype_counts[
  celltype_counts$sample %in% sample_order &
    celltype_counts$new_clusters %in% celltype_order,
  ,
  drop = FALSE
]

# Set factors on the small plotting table
celltype_counts$sample <- factor(
  celltype_counts$sample,
  levels = sample_order
)

celltype_counts$new_clusters <- factor(
  celltype_counts$new_clusters,
  levels = celltype_order
)

# Whole Tissue UMAP colors
cluster_colors <- c(
  "EC"          = "#e31a1c",
  "Pericyte"    = "#fb6a4a",
  "Ast-like"    = "#228b22",
  "Fibroblasts" = "#a6761d",
  "OPC/OL-like" = "#e69f00",
  "Immune"      = "#6a3d9a",
  "VLMC"        = "#00a6a6",
  "Neu-like"    = "#1f78b4",
  "Tumor"       = "gray50"
)

# Display labels
display_labels <- c(
  "EC"          = "Endothelial",
  "Pericyte"    = "Pericyte",
  "Ast-like"    = "Ast-like",
  "Fibroblasts" = "Fibroblasts",
  "OPC/OL-like" = "OPC/OL-like",
  "Immune"      = "Immune",
  "VLMC"        = "VLMC",
  "Neu-like"    = "Neu-like",
  "Tumor"       = "Tumor"
)

# Plot
LargeTissue_sample_prop_plot <- ggplot(
  celltype_counts,
  aes(
    x = sample,
    fill = new_clusters,
    weight = Freq
  )
) +
  geom_bar(
    position = "fill",
    width = 0.8
  ) +
  scale_fill_manual(
    values = cluster_colors,
    labels = display_labels,
    breaks = celltype_order,
    drop = FALSE
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(scale = 100),
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Proportion of cells (%)",
    fill = "Cell type"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  theme(
    axis.text.x = element_text(
      size = 18,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 18,
      colour = "black"
    ),
    axis.title = element_text(
      size = 20,
      colour = "black",
      face = "bold"
    ),
    legend.position = "right",
    legend.text = element_text(
      size = 14
    ),
    legend.title = element_text(
      size = 16,
      face = "bold"
    )
  )

# Display
LargeTissue_sample_prop_plot

# Save
output_dir <- "Vascular subtype analysis/Figures/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(
  filename = file.path(
    output_dir,
    "Whole Tissue Sample Cell Type Proportions.svg"
  ),
  plot = LargeTissue_sample_prop_plot,
  width = 8,
  height = 4,
  dpi = 600
)


#Spatial plot of ECs in whole tissue samples----
# EC Spatial Plots + EC Percentage per Sample
# Excludes LowQual_drop from EC percentage denominator

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Colors
normal_color <- "#228B22"
glioma_color <- "#377EB8"
brm_color <- "#D65F00"   # brighter burnt orange

# Output directory
output_dir <- "Vascular subtype analysis/Figures/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get EC cells
ec_cells <- WhichCells(
  whole_tissue,
  expression = new_clusters == "EC"
)

# Spatial plots
p_normal <- ImageDimPlot(
  whole_tissue,
  fov = "Normal_fov",
  cells = ec_cells,
  cols = normal_color,
  size = 1.2,
  border.color = "black",
  border.size = 0.15,
  dark.background = FALSE
) +
  ggtitle("Normal") +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )

p_gr3 <- ImageDimPlot(
  whole_tissue,
  fov = "fov",
  cells = ec_cells,
  cols = glioma_color,
  size = 1.2,
  border.color = "black",
  border.size = 0.15,
  dark.background = FALSE
) +
  ggtitle("Grade 3") +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )

p_gbm <- ImageDimPlot(
  whole_tissue,
  fov = "fov.3",
  cells = ec_cells,
  cols = glioma_color,
  size = 1.2,
  border.color = "black",
  border.size = 0.15,
  dark.background = FALSE
) +
  ggtitle("GBM") +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )

p_breast <- ImageDimPlot(
  whole_tissue,
  fov = "Breast_fov",
  cells = ec_cells,
  cols = brm_color,
  size = 1.2,
  border.color = "black",
  border.size = 0.15,
  dark.background = FALSE
) +
  ggtitle("Breast") +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )

p_kidney <- ImageDimPlot(
  whole_tissue,
  fov = "Kidney_fov",
  cells = ec_cells,
  cols = brm_color,
  size = 1.2,
  border.color = "black",
  border.size = 0.15,
  dark.background = FALSE
) +
  ggtitle("Kidney") +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )

p_lung <- ImageDimPlot(
  whole_tissue,
  fov = "Lung_fov",
  cells = ec_cells,
  cols = brm_color,
  size = 1.2,
  border.color = "black",
  border.size = 0.15,
  dark.background = FALSE
) +
  ggtitle("Lung") +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )

p_melanoma <- ImageDimPlot(
  whole_tissue,
  fov = "Melanoma_fov",
  cells = ec_cells,
  cols = brm_color,
  size = 1.2,
  border.color = "black",
  border.size = 0.15,
  dark.background = FALSE
) +
  ggtitle("Melanoma") +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )

p_thyroid <- ImageDimPlot(
  whole_tissue,
  fov = "Thyroid_fov",
  cells = ec_cells,
  cols = brm_color,
  size = 1.2,
  border.color = "black",
  border.size = 0.15,
  dark.background = FALSE
) +
  ggtitle("Thyroid") +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )


# Calculate EC percentage per sample
sample_vec <- as.character(whole_tissue$sample)
cluster_vec <- as.character(whole_tissue$new_clusters)

# Exclude NA, empty labels, and LowQual_drop
keep <- (
  !is.na(sample_vec) &
    sample_vec != "" &
    !is.na(cluster_vec) &
    cluster_vec != "" &
    cluster_vec != "LowQual_drop"
)

sample_vec <- sample_vec[keep]
cluster_vec <- cluster_vec[keep]

# Sample order
sample_order <- c(
  "Normal",
  "Gr_3",
  "GBM",
  "Breast",
  "Kidney",
  "Lung",
  "Melanoma",
  "Thyroid"
)

# Total valid cells per sample
sample_counts <- table(sample_vec)

# EC cells per sample
ec_counts <- table(
  sample_vec[cluster_vec == "EC"]
)

# Calculate EC percentage
ec_percent <- 100 *
  as.numeric(ec_counts[sample_order]) /
  as.numeric(sample_counts[sample_order])

# Create plotting data frame
ec_plot_df <- data.frame(
  sample = sample_order,
  EC_cells = as.numeric(ec_counts[sample_order]),
  Total_cells = as.numeric(sample_counts[sample_order]),
  EC_percent = ec_percent,
  stringsAsFactors = FALSE
)

# Display names
ec_plot_df$sample <- c(
  "Normal",
  "Grade 3",
  "GBM",
  "Breast",
  "Kidney",
  "Lung",
  "Melanoma",
  "Thyroid"
)

# Condition
ec_plot_df$condition <- c(
  "Normal",
  "Glioma",
  "Glioma",
  "BrM",
  "BrM",
  "BrM",
  "BrM",
  "BrM"
)

# Order samples in plot
ec_plot_df$sample <- factor(
  ec_plot_df$sample,
  levels = rev(c(
    "Normal",
    "Grade 3",
    "GBM",
    "Breast",
    "Kidney",
    "Lung",
    "Melanoma",
    "Thyroid"
  ))
)

p_percent <- ggplot(
  ec_plot_df,
  aes(
    x = EC_percent,
    y = sample,
    fill = condition
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = normal_color,
      "Glioma" = glioma_color,
      "BrM" = brm_color
    )
  ) +
  scale_x_continuous(
    name = "Endothelial cells (%) per Sample",
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_discrete(
    name = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 12,
      face = "bold"
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  ) +
  ggtitle(" ")

# Combine into 3 x 3 layout
EC_spatial_3x3 <- (
  p_normal + p_gr3 + p_gbm +
    p_breast + p_kidney + p_lung +
    p_melanoma + p_thyroid + p_percent
) +
  plot_layout(
    ncol = 3,
    nrow = 3,
    widths = rep(1, 3),
    heights = rep(1, 3)
  )

# Display
EC_spatial_3x3


# Export percentage table
write.csv(
  ec_plot_df,
  file = file.path(
    output_dir,
    "EC_percentage_by_sample.csv"
  ),
  row.names = FALSE
)

# Export complete 3 x 3 figure
ggsave(
  filename = file.path(
    output_dir,
    "tiff_EC_spatial_3x3_with_percentage.tiff"
  ),
  plot = EC_spatial_3x3,
  width = 12,
  height = 12
)

# Print the percentage table for checking
ec_plot_df


# Marker gene DotPlot - Large Tissue----

# Keep LowQual_drop as a valid identity so it does not become NA
all_identity_levels <- c(
  desired_order,
  "LowQual_drop"
)

Idents(whole_tissue) <- factor(
  whole_tissue$new_clusters,
  levels = all_identity_levels
)

# Marker genes
marker_genes <- c(
  # Endothelial
  "VWF", "PECAM1", "CDH5",
  
  # Pericyte
  "PDGFRB", "RGS5", "NOTCH3",
  
  # Astrocyte
  "AQP4", "GJA1", "SLC1A3", "SLC1A2", "ALDH1L1",
  
  # Fibroblasts
  "DCN", "FBLN1", "LAMA2", "COL12A1",
  
  # Oligodendrocyte / OPC
  "MOG", "MYRF", "MAG", "CLDN11", "CNDP1",
  
  # Immune
  "PTPRC", "ITGB2", "ITGAX", "ITGAM",
  
  # VLMC
  "SOX10", "CAPN3", "ERBB3", "GPNMB",
  
  # Neuronal
  "SLC17A7", "SLC17A6", "C1QL3", "RORB", "FSTL4",
  
  # Tumor
  "VEGFA", "CAV1", "TOP2A", "MKI67"
)

# Publication-quality DotPlot
chunk_marker_dotplot <- DotPlot(
  whole_tissue,
  features = marker_genes,
  idents = desired_order,
  dot.scale = 6
) +
  RotatedAxis() +
  scale_y_discrete(
    labels = display_labels
  ) +
  theme_santi.umap() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(family = "Helvetica"),
    axis.title = element_text(family = "Helvetica"),
    legend.text = element_text(family = "Helvetica"),
    legend.title = element_text(family = "Helvetica")
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 24
    ),
    axis.text.y = element_text(
      size = 26
    ),
    legend.text = element_text(
      size = 24
    ),
    legend.title = element_text(
      size = 26
    )
  )

# Save DotPlot
ggsave(
  filename = file.path(output_dir, "Helv_Marker Dotplot: Tumor Chunk.svg"),
  plot = chunk_marker_dotplot,
  width = 24,
  height = 10,
  dpi = 300
)



# Large Tissue: Cell Type Proportions by Sample----
# Each cell type sums to 100% across samples

library(ggplot2)
library(scales)

# Extract ONLY the two metadata columns needed
sample_vec <- as.character(whole_tissue$sample)
cluster_vec <- as.character(whole_tissue$new_clusters)

# Explicitly remove unwanted cells
keep <- (
  !is.na(sample_vec) &
    sample_vec != "" &
    !is.na(cluster_vec) &
    cluster_vec != "" &
    cluster_vec != "LowQual_drop"
)

# Build a small count table
celltype_counts <- as.data.frame(
  table(
    sample = sample_vec[keep],
    new_clusters = cluster_vec[keep],
    useNA = "no"
  ),
  stringsAsFactors = FALSE
)

# Remove zero-count combinations
celltype_counts <- celltype_counts[
  celltype_counts$Freq > 0,
  ,
  drop = FALSE
]

# Explicitly rename the actual metadata label
celltype_counts$sample[
  celltype_counts$sample == "Gr_3"
] <- "Grade 3"

# Explicit sample order
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

# Explicit cell-type order
celltype_order <- c(
  "EC",
  "Pericyte",
  "Ast-like",
  "Fibroblasts",
  "OPC/OL-like",
  "Immune",
  "VLMC",
  "Neu-like",
  "Tumor"
)

# Keep ONLY the expected cell types and samples
celltype_counts <- celltype_counts[
  celltype_counts$sample %in% sample_order &
    celltype_counts$new_clusters %in% celltype_order,
  ,
  drop = FALSE
]

# Set factors ONLY on the small plotting data frame
celltype_counts$sample <- factor(
  celltype_counts$sample,
  levels = sample_order
)

celltype_counts$new_clusters <- factor(
  celltype_counts$new_clusters,
  levels = celltype_order
)

# Display labels
display_labels <- c(
  "EC"          = "Endothelial",
  "Pericyte"    = "Pericyte",
  "Ast-like"    = "Ast-like",
  "Fibroblasts" = "Fibroblasts",
  "OPC/OL-like" = "OPC/OL-like",
  "Immune"      = "Immune",
  "VLMC"        = "VLMC",
  "Neu-like"    = "Neu-like",
  "Tumor"       = "Tumor"
)

# Plot
LargeTissue_celltype_prop_plot <- ggplot(
  celltype_counts,
  aes(
    x = new_clusters,
    fill = sample,
    weight = Freq
  )
) +
  geom_bar(
    position = "fill",
    width = 0.8
  ) +
  scale_fill_brewer(
    palette = "Set3",
    drop = FALSE
  ) +
  scale_x_discrete(
    labels = display_labels,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(scale = 100),
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Proportion of cell type (%)",
    fill = "Sample"
  ) +
  theme_classic() +
  theme( text = element_text(family = "Helvetica"), axis.text = element_text(family = "Helvetica"), axis.title = element_text(family = "Helvetica"), legend.text = element_text(family = "Helvetica"), legend.title = element_text(family = "Helvetica") ) +
  theme(
    axis.text.x = element_text(
      size = 18,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 18,
      colour = "black"
    ),
    axis.title = element_text(
      size = 20,
      colour = "black",
      face = "bold"
    ),
    legend.position = "right",
    legend.text = element_text(
      size = 12
    ),
    legend.title = element_text(
      size = 14
    )
  )

# Display
LargeTissue_celltype_prop_plot

# Save
output_dir <- "Vascular subtype analysis/Figures/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(
  filename = file.path(
    output_dir,
    "Helv_Large Tissue Cell Type Proportions by Sample.svg"
  ),
  plot = LargeTissue_celltype_prop_plot,
  width = 9,
  height = 6,
  dpi = 600
)


#Stacked barplots of proportions
#sample by condition

cond_prop <- whole_tissue@meta.data %>%
  count(sample, condition) %>%   # count cells per group
  group_by(condition) %>%
  mutate(prop = n / sum(n) * 100) %>%       # compute % per sample
  ungroup()


library(RColorBrewer)
#export svg as TMA cond prop plot, 600 x 600

cond_prop <- cond_prop %>%
  mutate(sample = recode(sample, "Gr_3" = "Grade 3"))

whole_tissue_prop <-ggplot(cond_prop, aes(x = condition, y = prop, fill = sample)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_brewer(palette = "Set3") +
  scale_x_discrete(
    labels = c(
      "Normal" = "Normal",
      "Primary Tumor" = "Primary",
      "Secondary Tumor" = "Metastasis"
    )
  ) +
  labs(
    x = NULL,
    y = "Cell proportion (%)",
    fill = "Sample"
  ) +
  theme_classic() +
  theme( text = element_text(family = "Helvetica"), axis.text = element_text(family = "Helvetica"), axis.title = element_text(family = "Helvetica"), legend.text = element_text(family = "Helvetica"), legend.title = element_text(family = "Helvetica") ) +
  theme(
    axis.text.x = element_text(size = 18, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title = element_text(size = 20, colour = "black", face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# Save
output_dir <- "Vascular subtype analysis/Figures/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(
  filename = file.path(
    output_dir,
    "Helv_Large Tissue Cond Proportion.svg"
  ),
  plot = whole_tissue_prop,
  width = 6,
  height = 6,
  dpi = 600
)

