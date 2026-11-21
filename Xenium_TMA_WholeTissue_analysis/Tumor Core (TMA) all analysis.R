#Xenium Tumor Core (TMA)

#Cell type annotation, cell composition
#DESeq2 Pseudobulk DEG analysis
#MFSD2A analysis
#Trajectory analysis

library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyr)
library(purrr)
library(csv)
library(RColorBrewer)
library(ggthemes)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(arrow)

theme_santi.umap<-function(base_size = 80, face ="bold", base_family= "Helvetica") {
  theme_few() %+replace%
    theme(axis.ticks = element_line(),
          strip.placement = "outside",
          panel.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(color = "black", size = 0.5),
          plot.title = element_blank(),
          axis.text = element_text(size = 28),
          axis.title = element_text(size = 34),
          legend.text = element_text(size = 36),
          legend.key.size = unit(4, "line")
    )
}


# Load TMA vascular subtype object----
TMA_full_vascular_subtypes <- readRDS(
  "TMA_full_vascular_subtypes.rds"
)





# UMAP - TMA----

# Define cluster order
desired_order <- c(
  "EC",
  "Pericyte",
  "Ast-like",
  "Mixed",
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
  "Mixed"       = "Mixed",
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
  "Mixed"       = "#969696",  # neutral gray
  "OPC/OL-like" = "#e69f00",  # orange
  "Immune"      = "#6a3d9a",  # purple
  "VLMC"        = "#00a6a6",  # teal
  "Neu-like"    = "#1f78b4",  # blue
  "Tumor"       = "gray50"    # gray
)

# Set identities
Idents(TMA_full_vascular_subtypes) <- factor(
  TMA_full_vascular_subtypes$new_clusters,
  levels = desired_order
)

# Get cells to plot
umap_cells <- rownames(TMA_full_vascular_subtypes[[]])[
  TMA_full_vascular_subtypes$new_clusters %in% desired_order
]

# Plot UMAP
TMA_UMAP <- DimPlot(
  TMA_full_vascular_subtypes,
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

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Export UMAP as SVG
ggsave(
  filename = file.path(
    output_dir,
    "UMAP: TMA Vascular Subtypes.tiff"
  ),
  plot = TMA_UMAP,
  width = 9,
  height = 8,
  dpi = 300
)


# Marker gene DotPlot - TMA----

# Set identities
Idents(TMA_full_vascular_subtypes) <- factor(
  TMA_full_vascular_subtypes$new_clusters,
  levels = desired_order
)

# Marker genes
marker.list_TMA_full <- list(
  
  # Endothelial / Pericyte 
  "EC/PC" = c(
    # Endothelial
    "VWF", "PECAM1", "CDH5",
    
    # Pericyte
    "PDGFRB", "RGS5", "NOTCH3"
  ),
  
  # Astrocyte
  "Ast-like" = c(
    "AQP4", "GJA1", "SLC1A3", "SLC1A2", "ALDH1L1"
  ),
  
  # Oligodendrocyte / OPC
  "OPC/OL-like" = c(
    "MOG", "MAG", "CLDN11", "CNDP1"
  ),
  
  # Immune
  "Immune" = c(
    "PTPRC", "ITGB2", "ITGAX", "ITGAM"
  ),
  
  #VLMC
  "VLMC" = c( "DCN",
              "COL12A1",
              "LAMA2",
              "FBLN1"
  ),
  
  # Neuronal
  "Neu-like" = c(
    "NDST4", "TAC1", "FSTL4"
  ),
  
  # Tumor
  "Tumor" = c(
    "VEGFA", "CAV1", "TOP2A", "MKI67"
  )
)




# Flatten marker list
marker_genes_TMA <- unique(
  unlist(marker.list_TMA_full)
)


# Publication-quality DotPlot
TMA_marker_dotplot <- DotPlot(
  TMA_full_vascular_subtypes,
  features = marker_genes_TMA,
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
  filename = file.path(
    output_dir,
    "Helv_Marker Dotplot: TMA Vascular Subtypes.svg"
  ),
  plot = TMA_marker_dotplot,
  width = 24,
  height = 10,
  dpi = 300
)






#Prop by cell type plot----
# Cell type proportions by sample
# Each cell type sums to 100% across samples

celltype_prop <- TMA_full_vascular_subtypes@meta.data %>%
  count(sample, new_clusters) %>%
  group_by(new_clusters) %>%
  mutate(prop = n / sum(n) * 100) %>%
  ungroup()


# Sample name cleanup
celltype_prop <- celltype_prop %>%
  mutate(
    sample = recode(
      sample,
      "G2" = "Grade 2",
      "G3" = "Grade 3"
    )
  )


# Cell type order
celltype_order <- c(
  "EC",
  "Pericyte",
  "Ast-like",
  "Mixed",
  "OPC/OL-like",
  "Immune",
  "VLMC",
  "Neu-like",
  "Tumor"
)


# Set cell type order
celltype_prop <- celltype_prop %>%
  mutate(
    new_clusters = factor(
      new_clusters,
      levels = celltype_order
    )
  )


# Plot
TMA_celltype_prop_plot <- ggplot(
  celltype_prop,
  aes(
    x = new_clusters,
    y = prop,
    fill = sample
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  scale_fill_brewer(
    palette = "Set3"
  ) +
  scale_x_discrete(
    labels = c(
      "EC"          = "Endothelial",
      "Pericyte"    = "Pericyte",
      "Ast-like"    = "Ast-like",
      "Mixed"       = "Mixed",
      "OPC/OL-like" = "OPC/OL-like",
      "Immune"      = "Immune",
      "VLMC"        = "VLMC",
      "Neu-like"    = "Neu-like",
      "Tumor"       = "Tumor"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    x = NULL,
    y = "Proportion of cell type (%)",
    fill = "Sample"
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
      size = 12
    ),
    legend.title = element_text(
      size = 14
    )
  )


# Save as SVG
ggsave(
  filename = file.path(
    output_dir,
    "Helv_TMA Cell Type Proportions by Sample.svg"
  ),
  plot = TMA_celltype_prop_plot,
  width = 9,
  height = 6,
  dpi = 600
)



#TMA prop plot by condition----
library(dplyr)
library(ggplot2)

# Sample order
sample_order <- c(
  "Normal",
  "G2",
  "G3",
  "GBM",
  "BR",
  "COL",
  "END",
  "RCC",
  "NSCLC",
  "MEL",
  "THY"
)

# Prepare data and calculate exact stack positions
plot_df <- cond_prop_TMA %>%
  mutate(
    sample = as.character(sample),
    condition = factor(
      condition,
      levels = c("Normal", "Gliomas", "Brain Met")
    ),
    sample = factor(
      sample,
      levels = sample_order
    )
  ) %>%
  group_by(condition) %>%
  arrange(sample, .by_group = TRUE) %>%
  mutate(
    ymin = c(0, head(cumsum(prop), -1)),
    ymax = cumsum(prop)
  ) %>%
  ungroup()

# Plot
TMA_cond_prop_plot <- ggplot(
  plot_df
) +
  geom_rect(
    aes(
      xmin = as.numeric(condition) - 0.4,
      xmax = as.numeric(condition) + 0.4,
      ymin = ymin,
      ymax = ymax,
      fill = sample
    )
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#8DD3C7",
      "G2" = "#FFFFB3",
      "G3" = "#BEBADA",
      "GBM" = "#FB8072",
      "BR" = "#80B1D3",
      "COL" = "#FDB462",
      "END" = "#B3DE69",
      "RCC" = "#FCCDE5",
      "NSCLC" = "#D9D9D9",
      "MEL" = "#BC80BD",
      "THY" = "#CCEBC5"
    ),
    breaks = sample_order,
    labels = c(
      "Normal",
      "Grade 2",
      "Grade 3",
      "GBM",
      "BR",
      "COL",
      "END",
      "RCC",
      "NSCLC",
      "MEL",
      "THY"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = 1:3,
    labels = c(
      "Normal",
      "Gliomas",
      "Brain Met"
    ),
    limits = c(0.5, 3.5)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Cell proportion (%)",
    fill = "Sample"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(
      family = "Helvetica",
      size = 18,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      family = "Helvetica",
      size = 18,
      colour = "black"
    ),
    axis.title = element_text(
      family = "Helvetica",
      size = 20,
      colour = "black",
      face = "bold"
    ),
    legend.text = element_text(
      family = "Helvetica",
      size = 12
    ),
    legend.title = element_text(
      family = "Helvetica",
      size = 14
    ),
    legend.position = "right"
  )

TMA_cond_prop_plot

ggsave(
  filename = file.path(
    output_dir,
    "Helv_TMA Condition Proportions by Sample.svg"
  ),
  plot = TMA_cond_prop_plot,
  width = 6,
  height = 6,
  dpi = 600
)
# TMA EC BBB gene expression boxplots----

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)



# Load TMA EC/PC subcluster object
TMA_full_EC.PC_subcl.EC <- readRDS(
  "TMA_full_EC.PC_subcl.EC.rds"
)

# Create updated cluster annotation
TMA_full_EC.PC_subcl.EC@meta.data <- 
  TMA_full_EC.PC_subcl.EC@meta.data %>%
  mutate(
    new_clusters = if_else(
      clusters == "EC/PC" & !is.na(vasc_sub),
      vasc_sub,
      clusters
    )
  )

# Subset endothelial cells
TMA_EC_sub <- subset(
  TMA_full_EC.PC_subcl.EC,
  subset = new_clusters == "EC"
)

# Sample order
sample_order <- c(
  "Normal",
  "Grade 2",
  "Grade 3",
  "GBM",
  "BR",
  "COL",
  "END",
  "RCC",
  "NSCLC",
  "MEL",
  "THY"
)

# Genes to plot
genes_to_plot <- c(
  "PLVAP",
  "COL4A1",
  "MFSD2A",
  "CLDN5",
  "ABCB1",
  "ABCG2"
)

# Fetch expression data
plot_df <- FetchData(
  TMA_EC_sub,
  vars = c(genes_to_plot, "sample"),
  assay = "Xenium"
)

# Convert sample to character before modifying
plot_df$sample <- as.character(plot_df$sample)

# Rename G2/G3
plot_df$sample[plot_df$sample == "G2"] <- "Grade 2"
plot_df$sample[plot_df$sample == "G3"] <- "Grade 3"

# Remove cells with missing or empty sample labels
plot_df <- plot_df[
  !is.na(plot_df$sample) &
    plot_df$sample != "",
  ,
  drop = FALSE
]

# Keep only samples in the desired sample order
plot_df <- plot_df[
  plot_df$sample %in% sample_order,
  ,
  drop = FALSE
]

# Convert to ordered factor
plot_df$sample <- factor(
  plot_df$sample,
  levels = sample_order
)

# Check sample distribution
print(table(plot_df$sample, useNA = "ifany"))

# Function to generate boxplots
make_boxplot <- function(gene) {
  
  ggplot(
    plot_df,
    aes(
      x = sample,
      y = .data[[gene]]
    )
  ) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.45
    ) +
    geom_jitter(
      width = 0.18,
      height = 0,
      size = 0.6,
      alpha = 0.20,
      color = "black"
    ) +
    theme_classic() +
    labs(
      title = gene,
      x = NULL,
      y = "Expression"
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 20
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        color = "black",
        size = 12
      ),
      axis.text.y = element_text(
        color = "black",
        size = 12
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 14
      ),
      plot.margin = margin(
        t = 8,
        r = 8,
        b = 8,
        l = 8
      )
    )
}

# Generate individual plots
p_PLVAP  <- make_boxplot("PLVAP")
p_COL4A1 <- make_boxplot("COL4A1")
p_MFSD2A <- make_boxplot("MFSD2A")
p_CLDN5  <- make_boxplot("CLDN5")
p_ABCB1  <- make_boxplot("ABCB1")
p_ABCG2  <- make_boxplot("ABCG2")

# Combine into 2 × 3 figure
TMA_EC_BBB_genes <- (
  p_PLVAP + p_COL4A1 + p_MFSD2A +
    p_CLDN5 + p_ABCB1 + p_ABCG2
) +
  plot_layout(
    ncol = 3,
    nrow = 2
  )

# Display
TMA_EC_BBB_genes

# Set directory where you want to save
output_dir <- "Vascular subtype analysis/Figures/"

# Make sure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Export UMAP as SVG
ggsave(
  filename = file.path(output_dir, "TMA genes boxplots.svg"),
  plot = TMA_EC_BBB_genes,
  width = 18,
  height = 10,
  dpi = 300
)









# DESeq2 DEG analysis -----

# Install/load DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

library(Seurat)
library(DESeq2)
library(Matrix)
library(dplyr)
library(tibble)


# Update Seurat object
TMA_full <- UpdateSeuratObject(TMA_full)


# Rename Neurons -> Neu-like BEFORE any analysis
TMA_full$new_clusters <- as.character(TMA_full$new_clusters)

TMA_full$new_clusters[
  TMA_full$new_clusters == "Neurons"
] <- "Neu-like"


# Check cell/sample composition
table(
  TMA_full$sample,
  TMA_full$new_clusters
)

table(
  TMA_full$sample_id,
  TMA_full$new_clusters
)


# Define condition
# Normal = normal brain tissue
# Primary Tumor = primary gliomas (G2, G3, GBM)
# Secondary Tumor = brain metastases

TMA_full$condition <- dplyr::case_when(
  TMA_full$sample == "Normal" ~ "Normal",
  TMA_full$sample %in% c("G2", "G3", "GBM") ~ "Primary Tumor",
  TMA_full$sample %in% c(
    "BR", "COL", "END", "LU_NSC", "MEL", "RCC", "THY"
  ) ~ "Secondary Tumor",
  TRUE ~ NA_character_
)

TMA_full$condition <- factor(
  TMA_full$condition,
  levels = c(
    "Normal",
    "Primary Tumor",
    "Secondary Tumor"
  )
)


# Check condition assignment
table(
  TMA_full$sample,
  TMA_full$condition,
  useNA = "ifany"
)


# Loading gene panel
gene.list <- read.csv(
  "Xenium Human Brain Panel Gene List-Santi.csv"
)


# Filter out NVU custom genes
nvu_genes <- gene.list %>%
  filter(Santi.Annotations == "NVU") %>%
  pull(gene)

nvu_genes2 <- nvu_genes[
  !(nvu_genes %in% c(
    "ALDH1L1",
    "AQP4",
    "GJA1",
    "GFAP",
    "MLC1",
    "BDNF",
    "GABBR1",
    "NDP",
    "ATP1B2",
    "WNT7A",
    "ALB",
    "SLC1A2"
  ))
]

nvu_genes <- nvu_genes2


# Pseudobulk function
pseudobulk_counts <- function(obj, celltype,
                              celltype_col = "new_clusters",
                              sample_col = "sample_id",
                              cond_col = "condition") {
  
  # Select cells matching the cell type
  cells_use <- rownames(obj@meta.data)[
    obj@meta.data[[celltype_col]] == celltype
  ]
  
  if (length(cells_use) == 0) {
    stop(
      paste0(
        "No cells found for celltype: ",
        celltype
      )
    )
  }
  
  # Subset Seurat object
  sub <- subset(
    obj,
    cells = cells_use
  )
  
  # Seurat v5: use layer instead of slot
  counts <- GetAssayData(
    sub,
    layer = "counts"
  )
  
  # Align metadata to count matrix columns
  meta <- sub@meta.data[
    colnames(counts),
    ,
    drop = FALSE
  ] %>%
    dplyr::select(
      all_of(c(sample_col, cond_col))
    ) %>%
    dplyr::mutate(
      sample_id = .data[[sample_col]],
      condition = .data[[cond_col]]
    )
  
  # Remove cells with missing sample/condition information
  keep_cells <- !is.na(meta$sample_id) &
    !is.na(meta$condition)
  
  meta <- meta[
    keep_cells,
    ,
    drop = FALSE
  ]
  
  counts <- counts[
    ,
    keep_cells,
    drop = FALSE
  ]
  
  # Pseudobulk:
  # sum counts across cells within each biological sample
  pb <- t(
    rowsum(
      t(as.matrix(counts)),
      group = meta$sample_id
    )
  )
  
  # Sample-level metadata
  coldata <- meta %>%
    dplyr::distinct(
      sample_id,
      condition
    ) %>%
    dplyr::arrange(sample_id)
  
  # Make sure count matrix and metadata have identical sample order
  pb <- pb[
    ,
    coldata$sample_id,
    drop = FALSE
  ]
  
  # Row names must match count matrix column names
  rownames(coldata) <- coldata$sample_id
  
  return(
    list(
      counts = pb,
      coldata = coldata
    )
  )
}


# DESeq2 function
run_DESeq2 <- function(pb_counts, coldata, cond1, cond2) {
  
  # Set condition factor
  coldata$condition <- factor(
    coldata$condition,
    levels = c(
      "Normal",
      "Primary Tumor",
      "Secondary Tumor"
    )
  )
  
  # Check that both conditions are present
  tab <- table(coldata$condition)
  
  n1 <- ifelse(
    cond1 %in% names(tab),
    tab[cond1],
    0
  )
  
  n2 <- ifelse(
    cond2 %in% names(tab),
    tab[cond2],
    0
  )
  
  # Skip comparison if insufficient biological replicates
  if (n1 < 2 | n2 < 2) {
    message(
      "Skipping ",
      cond1,
      " vs ",
      cond2,
      " (need >=2 biological replicates each). Found: ",
      n1,
      " and ",
      n2
    )
    return(NULL)
  }
  
  # Keep only the two conditions being compared
  keep <- coldata$condition %in% c(
    cond1,
    cond2
  )
  
  coldata <- droplevels(
    coldata[
      keep,
      ,
      drop = FALSE
    ]
  )
  
  pb_counts <- pb_counts[
    ,
    rownames(coldata),
    drop = FALSE
  ]
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_counts),
    colData = coldata,
    design = ~ condition
  )
  
  # Filter low-count genes
  keep_genes <- rowSums(
    counts(dds) >= 10
  ) >= 3
  
  if (sum(keep_genes) == 0) {
    message(
      "Skipping ",
      cond1,
      " vs ",
      cond2,
      ": no genes passed the count filter."
    )
    return(NULL)
  }
  
  dds <- dds[
    keep_genes,
  ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(
    dds,
    contrast = c(
      "condition",
      cond1,
      cond2
    )
  )
  
  return(
    as.data.frame(res) %>%
      tibble::rownames_to_column("gene")
  )
}


# Output directories
outdir <- "Vascular subtype analysis/TMAs Condition Pairwise Comparisons_noLuReg/"

raw_dir <- file.path(
  outdir,
  "Raw Outputs"
)

filt_dir <- file.path(
  outdir,
  "Filtered Outputs"
)

nvu_dir <- file.path(
  outdir,
  "Custom NVU Filtered DEGs"
)

dir.create(
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  filt_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  nvu_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Function for exporting DEG results
write_deg_csv <- function(
    res_df,
    celltype,
    comparison,
    outdir
) {
  
  # Don't try to write NULL results
  if (is.null(res_df)) {
    return(invisible(NULL))
  }
  
  fname <- paste0(
    gsub("/", "_", celltype),
    "_",
    comparison,
    ".csv"
  )
  
  write.csv(
    res_df,
    file = file.path(
      outdir,
      fname
    ),
    row.names = FALSE
  )
}


# Function for filtering significant DEGs
filter_degs <- function(
    res_df,
    lfc_cutoff = 1,
    padj_cutoff = 0.05
) {
  
  # Return NULL if DESeq2 comparison was skipped
  if (is.null(res_df)) {
    return(NULL)
  }
  
  res_df %>%
    dplyr::filter(
      !is.na(padj),
      padj < padj_cutoff,
      abs(log2FoldChange) > lfc_cutoff
    )
}


# Function for filtering significant NVU genes
filter_sig_nvu <- function(
    res_df,
    nvu_genes,
    lfc_cutoff = 1,
    padj_cutoff = 0.05
) {
  
  # Return NULL if DESeq2 comparison was skipped
  if (is.null(res_df)) {
    return(NULL)
  }
  
  res_df %>%
    dplyr::filter(
      gene %in% nvu_genes,
      !is.na(padj),
      padj < padj_cutoff,
      abs(log2FoldChange) > lfc_cutoff
    )
}


# Cell types
celltypes <- c(
  "EC",
  "Pericyte",
  "VLMC",
  "Mixed",
  "Ast-like",
  "Immune",
  "OPC/OL-like",
  "Neu-like",
  "Tumor"
)


# Run DESeq2
results_all <- list()

for (ct in celltypes) {
  
  cat("\n========================================\n")
  cat("Processing:", ct, "\n")
  cat("========================================\n")
  
  # Generate pseudobulk counts
  pb <- pseudobulk_counts(
    TMA_full,
    ct
  )
  
  # Show number of biological samples by condition
  print(
    table(
      pb$coldata$condition
    )
  )
  
  # Show sample IDs
  print(
    table(
      pb$coldata$condition,
      pb$coldata$sample_id
    )
  )
  
  
  # Primary vs Normal
  cat("\nPrimary Tumor vs Normal\n")
  
  res_np <- run_DESeq2(
    pb$counts,
    pb$coldata,
    "Primary Tumor",
    "Normal"
  )
  
  
  # Secondary vs Normal
  cat("\nSecondary Tumor vs Normal\n")
  
  res_ns <- run_DESeq2(
    pb$counts,
    pb$coldata,
    "Secondary Tumor",
    "Normal"
  )
  
  
  # Secondary vs Primary
  cat("\nSecondary Tumor vs Primary Tumor\n")
  
  res_sp <- run_DESeq2(
    pb$counts,
    pb$coldata,
    "Secondary Tumor",
    "Primary Tumor"
  )
  
  
  # Store raw results in R
  results_all[[ct]] <- list(
    Primary_vs_Normal = res_np,
    Secondary_vs_Normal = res_ns,
    Secondary_vs_Primary = res_sp
  )
  
  
  # RAW results
  write_deg_csv(
    res_np,
    ct,
    "Primary_vs_Normal",
    raw_dir
  )
  
  write_deg_csv(
    res_ns,
    ct,
    "Secondary_vs_Normal",
    raw_dir
  )
  
  write_deg_csv(
    res_sp,
    ct,
    "Secondary_vs_Primary",
    raw_dir
  )
  
  
  # FILTERED results
  sig_np <- filter_degs(
    res_np
  )
  
  sig_ns <- filter_degs(
    res_ns
  )
  
  sig_sp <- filter_degs(
    res_sp
  )
  
  write_deg_csv(
    sig_np,
    ct,
    "Primary_vs_Normal_sig",
    filt_dir
  )
  
  write_deg_csv(
    sig_ns,
    ct,
    "Secondary_vs_Normal_sig",
    filt_dir
  )
  
  write_deg_csv(
    sig_sp,
    ct,
    "Secondary_vs_Primary_sig",
    filt_dir
  )
  
  
  # NVU SIGNIFICANT ONLY
  sig_nvu_np <- filter_sig_nvu(
    res_np,
    nvu_genes
  )
  
  sig_nvu_ns <- filter_sig_nvu(
    res_ns,
    nvu_genes
  )
  
  sig_nvu_sp <- filter_sig_nvu(
    res_sp,
    nvu_genes
  )
  
  write_deg_csv(
    sig_nvu_np,
    ct,
    "Primary_vs_Normal_NVU_sig",
    nvu_dir
  )
  
  write_deg_csv(
    sig_nvu_ns,
    ct,
    "Secondary_vs_Normal_NVU_sig",
    nvu_dir
  )
  
  write_deg_csv(
    sig_nvu_sp,
    ct,
    "Secondary_vs_Primary_NVU_sig",
    nvu_dir
  )
  
  
  # Report number of significant genes
  cat("\nSignificant DEGs:\n")
  
  cat(
    "  Primary vs Normal:",
    ifelse(
      is.null(sig_np),
      0,
      nrow(sig_np)
    ),
    "\n"
  )
  
  cat(
    "  Secondary vs Normal:",
    ifelse(
      is.null(sig_ns),
      0,
      nrow(sig_ns)
    ),
    "\n"
  )
  
  cat(
    "  Secondary vs Primary:",
    ifelse(
      is.null(sig_sp),
      0,
      nrow(sig_sp)
    ),
    "\n"
  )
}






#Comparing Samples
#Differential expression analysis was performed using DESeq2 on pseudobulked counts aggregated per sample. Gene expression differences between sample groups were assessed using a negative binomial generalized linear model with sample group as the design factor. P-values were adjusted for multiple testing using the Benjamini–Hochberg procedure. Genes with |log2 fold change| > 1 and FDR < 0.05 were considered significant.
# DESeq2: Sample Pairwise Comparisons -----

# Using:
# - Negative binomial distribution
# - Sample-level pseudobulk counts
# - Biological replicates = sample_id
# - FDR-adjusted p-values using Benjamini-Hochberg
# - Significant DEGs: |log2FC| > 1 and padj < 0.05


# Load packages
library(Seurat)
library(DESeq2)
library(dplyr)
library(tibble)
library(Matrix)


# Update Seurat object
TMA_full <- UpdateSeuratObject(TMA_full)


# Rename Neurons to Neu-like before any analysis
TMA_full$new_clusters <- as.character(TMA_full$new_clusters)

TMA_full$new_clusters[
  TMA_full$new_clusters == "Neurons"
] <- "Neu-like"


# Check updated cell-type annotations
table(
  TMA_full$new_clusters
)


# Add sample_group
# Each sample is treated as its own group for sample-to-sample comparisons
TMA_full$sample_group <- TMA_full$sample


# Check sample/sample_id composition
table(
  TMA_full$sample,
  TMA_full$new_clusters
)

table(
  TMA_full$sample_id,
  TMA_full$new_clusters
)


# Check mapping between sample and sample_id
table(
  TMA_full$sample,
  TMA_full$sample_id
)


# Pseudobulk function
pseudobulk_counts <- function(obj, celltype,
                              celltype_col = "new_clusters",
                              sample_col = "sample_id",
                              group_col = "sample_group") {
  
  # Select cells matching the cell type
  cells_use <- rownames(obj@meta.data)[
    obj@meta.data[[celltype_col]] == celltype
  ]
  
  if (length(cells_use) == 0) {
    stop(
      paste0(
        "No cells found for celltype: ",
        celltype
      )
    )
  }
  
  # Subset Seurat object
  sub <- subset(
    obj,
    cells = cells_use
  )
  
  # Seurat v5: use layer instead of slot
  counts <- GetAssayData(
    sub,
    layer = "counts"
  )
  
  # Align metadata to count matrix columns
  meta <- sub@meta.data[
    colnames(counts),
    ,
    drop = FALSE
  ] %>%
    dplyr::select(
      all_of(c(sample_col, group_col))
    ) %>%
    dplyr::mutate(
      sample_id = .data[[sample_col]],
      sample_group = .data[[group_col]]
    )
  
  # Remove cells with missing sample information
  keep_cells <- !is.na(meta$sample_id) &
    !is.na(meta$sample_group)
  
  meta <- meta[
    keep_cells,
    ,
    drop = FALSE
  ]
  
  counts <- counts[
    ,
    keep_cells,
    drop = FALSE
  ]
  
  # Pseudobulk:
  # sum counts across cells belonging to the same biological sample
  pb <- t(
    rowsum(
      t(as.matrix(counts)),
      group = meta$sample_id
    )
  )
  
  # Sample-level metadata
  coldata <- meta %>%
    dplyr::distinct(
      sample_id,
      sample_group
    ) %>%
    dplyr::arrange(sample_id)
  
  # Make sure pseudobulk columns match coldata
  pb <- pb[
    ,
    coldata$sample_id,
    drop = FALSE
  ]
  
  # Row names must match count matrix column names
  rownames(coldata) <- coldata$sample_id
  
  return(
    list(
      counts = pb,
      coldata = coldata
    )
  )
}


# DESeq2 function
run_DESeq2_group <- function(pb_counts,
                             coldata,
                             group_col = "sample_group",
                             group1,
                             group2,
                             min_reps = 2) {
  
  # Convert group to factor
  coldata[[group_col]] <- factor(
    coldata[[group_col]]
  )
  
  # Count biological replicates
  tab <- table(
    coldata[[group_col]]
  )
  
  n1 <- ifelse(
    group1 %in% names(tab),
    tab[group1],
    0
  )
  
  n2 <- ifelse(
    group2 %in% names(tab),
    tab[group2],
    0
  )
  
  # Skip comparisons with insufficient replicates
  if (n1 < min_reps | n2 < min_reps) {
    message(
      "Skipping ",
      group1,
      " vs ",
      group2,
      " (need >= ",
      min_reps,
      " biological replicates each). Found: ",
      n1,
      " and ",
      n2
    )
    
    return(NULL)
  }
  
  # Keep only the two groups being compared
  keep <- coldata[[group_col]] %in% c(
    group1,
    group2
  )
  
  coldata <- coldata[
    keep,
    ,
    drop = FALSE
  ]
  
  coldata[[group_col]] <- droplevels(
    coldata[[group_col]]
  )
  
  # Match pseudobulk samples to coldata
  pb_counts <- pb_counts[
    ,
    rownames(coldata),
    drop = FALSE
  ]
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_counts),
    colData = coldata,
    design = as.formula(
      paste("~", group_col)
    )
  )
  
  # Filter low-count genes
  keep_genes <- rowSums(
    counts(dds) >= 10
  ) >= 3
  
  if (sum(keep_genes) == 0) {
    message(
      "Skipping ",
      group1,
      " vs ",
      group2,
      ": no genes passed the count filter."
    )
    
    return(NULL)
  }
  
  dds <- dds[
    keep_genes,
  ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(
    dds,
    contrast = c(
      group_col,
      group1,
      group2
    )
  )
  
  return(
    as.data.frame(res) %>%
      tibble::rownames_to_column("gene")
  )
}


# Function for exporting DEG results
write_deg_csv <- function(res_df,
                          celltype,
                          comparison,
                          outdir) {
  
  # Do nothing for skipped comparisons
  if (is.null(res_df)) {
    return(invisible(NULL))
  }
  
  fname <- paste0(
    gsub("/", "_", celltype),
    "_",
    comparison,
    ".csv"
  )
  
  write.csv(
    res_df,
    file = file.path(
      outdir,
      fname
    ),
    row.names = FALSE
  )
}


# Function for filtering significant DEGs
filter_degs <- function(res_df,
                        lfc_cutoff = 1,
                        padj_cutoff = 0.05) {
  
  # Return NULL if comparison was skipped
  if (is.null(res_df)) {
    return(NULL)
  }
  
  res_df %>%
    dplyr::filter(
      !is.na(padj),
      padj < padj_cutoff,
      abs(log2FoldChange) > lfc_cutoff
    )
}


# Function for filtering significant NVU genes
filter_sig_nvu <- function(res_df,
                           nvu_genes,
                           lfc_cutoff = 1,
                           padj_cutoff = 0.05) {
  
  # Return NULL if comparison was skipped
  if (is.null(res_df)) {
    return(NULL)
  }
  
  res_df %>%
    dplyr::filter(
      gene %in% nvu_genes,
      !is.na(padj),
      padj < padj_cutoff,
      abs(log2FoldChange) > lfc_cutoff
    )
}


# Output directories
outdir <- paste0(
  "Vascular subtype analysis/",
  "TMAs Sample Pairwise Comparisons(noLuReg)/"
)

raw_dir <- file.path(
  outdir,
  "Raw Outputs"
)

filt_dir <- file.path(
  outdir,
  "Filtered Outputs"
)

nvu_dir <- file.path(
  outdir,
  "Custom NVU Filtered DEGs"
)


dir.create(
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  filt_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  nvu_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Pairwise sample comparisons
sample_groups <- sort(
  unique(
    na.omit(TMA_full$sample_group)
  )
)

comparisons <- t(
  combn(
    sample_groups,
    2
  )
)

comparisons <- as.data.frame(
  comparisons,
  stringsAsFactors = FALSE
)

colnames(comparisons) <- c(
  "group1",
  "group2"
)


# Force Normal to always be the second group
# This makes the interpretation:
#   group1 vs Normal
# rather than:
#   Normal vs group1

comparisons <- comparisons %>%
  mutate(
    g1 = ifelse(
      group1 == "Normal",
      group2,
      group1
    ),
    g2 = ifelse(
      group1 == "Normal",
      group1,
      group2
    )
  ) %>%
  select(
    g1,
    g2
  )


# Cell types
celltypes <- c(
  "EC",
  "Pericyte",
  "VLMC",
  "Mixed",
  "Ast-like",
  "Immune",
  "OPC/OL-like",
  "Neu-like",
  "Tumor"
)


# Store results in R
results_sample_pairwise <- list()


# Run DESeq2
for (ct in celltypes) {
  
  cat("\n========================================\n")
  cat("Processing cell type:", ct, "\n")
  cat("========================================\n")
  
  # Generate pseudobulk counts
  pb <- pseudobulk_counts(
    TMA_full,
    ct
  )
  
  # Show samples available for this cell type
  cat("\nSamples available:\n")
  print(
    pb$coldata
  )
  
  cat("\nSample groups:\n")
  print(
    table(
      pb$coldata$sample_group
    )
  )
  
  
  # Initialize results for this cell type
  results_sample_pairwise[[ct]] <- list()
  
  
  # Run every pairwise comparison
  for (i in seq_len(nrow(comparisons))) {
    
    g1 <- comparisons$g1[i]
    g2 <- comparisons$g2[i]
    
    comparison_name <- paste0(
      g1,
      "_vs_",
      g2
    )
    
    cat(
      "\nRunning:",
      g1,
      "vs",
      g2,
      "\n"
    )
    
    
    # Run DESeq2
    res <- run_DESeq2_group(
      pb_counts = pb$counts,
      coldata = pb$coldata,
      group_col = "sample_group",
      group1 = g1,
      group2 = g2
    )
    
    
    # Store result
    results_sample_pairwise[[ct]][[
      comparison_name
    ]] <- res
    
    
    # Skip everything else if comparison failed
    if (is.null(res)) {
      cat(
        "  -> Skipped\n"
      )
      next
    }
    
    
    # RAW results
    write_deg_csv(
      res,
      ct,
      comparison_name,
      raw_dir
    )
    
    
    # FILTERED results
    sig_res <- filter_degs(
      res,
      lfc_cutoff = 1,
      padj_cutoff = 0.05
    )
    
    write_deg_csv(
      sig_res,
      ct,
      paste0(
        comparison_name,
        "_sig"
      ),
      filt_dir
    )
    
    
    # NVU significant genes
    sig_nvu <- filter_sig_nvu(
      res,
      nvu_genes,
      lfc_cutoff = 1,
      padj_cutoff = 0.05
    )
    
    write_deg_csv(
      sig_nvu,
      ct,
      paste0(
        comparison_name,
        "_NVU_sig"
      ),
      nvu_dir
    )
    
    
    # Report results
    cat(
      "  -> Total genes:",
      nrow(res),
      "\n"
    )
    
    cat(
      "  -> Significant DEGs:",
      nrow(sig_res),
      "\n"
    )
    
    cat(
      "  -> Significant NVU genes:",
      nrow(sig_nvu),
      "\n"
    )
  }
}


cat("\n========================================\n")
cat("Sample pairwise DESeq2 analysis complete.\n")
cat("========================================\n")






#Creating ShinyApp to make heatmaps
#ShinyApp----
#export shinyapp
reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
           "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}


devtools::install_github("SGDDNB/ShinyCell")
library(ShinyCell)


levels(x = TMA_full) <- c('EC/PC', 'Ast-like', 'Immune', "OPC/OL-like", "Neu-like", "Tumor")


scConf = createConfig(TMA_full)
makeShinyApp(TMA_full, gex.assay = "Xenium", scConf, gene.mapping = TRUE,
             shiny.title = "TMA_full") 



#Export ShinyApp per cluster

seurat_obj <- TMA_full
new_clusters_to_split <- levels(seurat_obj$new_clusters)

base_outdir <- "TMA_ShinyCell_per_cluster_apps"
dir.create(base_outdir, recursive = TRUE, showWarnings = FALSE)

for (cl in new_clusters_to_split) {
  
  cat("Creating app for cluster:", cl, "\n")
  
  # Subset by cluster
  obj_sub <- subset(seurat_obj, subset = new_clusters == cl)
  
  # Create config
  scConf <- createConfig(obj_sub)
  
  # Build Shiny app (creates folder named "shinyApp")
  makeShinyApp(
    obj_sub,
    scConf = scConf,
    gex.assay = "Xenium",   # change if different
    gex.slot  = "counts",   # change if different
    shiny.title = paste0("TMA ShinyCell: ", cl)
  )
  
  # Destination folder
  new_folder <- file.path(base_outdir, paste0("ShinyApp_", cl))
  
  # Remove destination if it exists
  if (dir.exists(new_folder)) {
    unlink(new_folder, recursive = TRUE)
  }
  
  # Copy generated shinyApp
  dir.create(new_folder, recursive = TRUE)
  file.copy("shinyApp", new_folder, recursive = TRUE)
  
  # Remove original shinyApp so next iteration is clean
  unlink("shinyApp", recursive = TRUE)
  
  cat("Finished:", cl, "\n\n")
}






#Volcano Plots of Condition Comparisons EC/PC
library(ggplot2)
library(dplyr)
library(ggrepel)

plot_volcano <- function(df, title,
                         lfc_cutoff = 1,
                         padj_cutoff = 0.05,
                         top_n = 10) {
  
  # Remove NA padj
  df <- df %>%
    filter(!is.na(padj))
  
  # Add columns
  df <- df %>%
    mutate(
      negLog10Padj = -log10(padj),
      Significance = case_when(
        padj < padj_cutoff & log2FoldChange > lfc_cutoff  ~ "Up",
        padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
        TRUE ~ "Not Sig"
      )
    )
  
  # Select top genes for labeling
  top_genes <- df %>%
    arrange(padj)  %>%
    head(50)
  
  # Symmetric x-axis
  max_lfc <- max(abs(df$log2FoldChange), na.rm = TRUE)
  
  ggplot(df, aes(x = log2FoldChange, y = negLog10Padj, color = Significance)) +
    geom_point(alpha = 0.7, size = 2) +
    
    scale_color_manual(values = c(
      "Up" = "#D62728",
      "Down" = "#1F77B4",
      "Not Sig" = "grey70"
    )) +
    
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf
    ) +
    
    coord_cartesian(xlim = c(-max_lfc, max_lfc)) +
    
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    )
}


input_dir <- "TMA6 & TMA8 Met/More Filtered analysis/Condition Pairwise Comparisons/Custom NVU Filtered DEGs/"

all_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)




#Plotting volcano plots
volcano_plots <- lapply(names(EC_PC_list), function(name) {
  plot_volcano(EC_PC_list[[name]], title = name)
})

names(volcano_plots) <- names(EC_PC_list)

volcano_plots$Primary_vs_Normal_NVU_sig
volcano_plots$Secondary_vs_Normal_NVU_sig
volcano_plots$Secondary_vs_Primary_NVU_sig

# Define output directory
output_dir <- "TMA6 & TMA8 Met/More Filtered analysis/Condition Pairwise Comparisons/Custom NVU Filtered DEGs/Volcano Plots/"

# Create folder if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (i in seq_along(volcano_plots)) {
  ggsave(
    filename = file.path(output_dir,
                         paste0("Volcano_", names(volcano_plots)[i], ".svg")),
    plot = volcano_plots[[i]],
    width = 6,
    height = 5,
    dpi = 300
  )
}





#edited 4.7.26


#make volcano plot of Tumor Secondary vs Primary

DEG_list_met_vs_pri_tumor <- read_csv("TMA6 & TMA8 Met/More Filtered analysis/Condition Pairwise Comparisons_noLuReg/Custom NVU Filtered DEGs/Tumor_Secondary_vs_Primary_NVU_sig.csv")


# Create required columns
DEG_list_met_vs_pri_tumor_filtered <- DEG_list_met_vs_pri_tumor %>%
  filter(!gene %in% c("GFAP", "MLC1", "SLC1A2"))


df <- DEG_list_met_vs_pri_tumor_filtered %>%
  mutate(
    negLog10Padj = -log10(padj),
    Significance = case_when(
      log2FoldChange > 1 & padj < 0.05 ~ "Up",
      log2FoldChange < -1 & padj < 0.05 ~ "Down",
      TRUE ~ "Not Sig"
    )
  )



# Select top genes to label (top 10 by smallest padj)
top_genes <- df %>%
  arrange(padj) %>%
  slice_head(n = 50)

# Plot
ggplot(df, aes(x = log2FoldChange, y = negLog10Padj, color = Significance)) +
  geom_point(alpha = 0.7, size = 2) +
  
  scale_color_manual(values = c(
    "Up" = "#D62728",
    "Down" = "#1F77B4",
    "Not Sig" = "grey70"
  )) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 5,
    max.overlaps = Inf
  ) +
  
  coord_cartesian(
    xlim = c(
      -max(abs(df$log2FoldChange)),
      max(abs(df$log2FoldChange))
    )
  ) +
  
  labs(
    title = " ",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  )










library(ggplot2)
library(dplyr)
library(ggrepel)

input_dir <- "TMA6 & TMA8 Met/More Filtered analysis/Condition Pairwise Comparisons/Custom NVU Filtered DEGs/"
all_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

main_output_dir <- "TMA6 & TMA8 Met/More Filtered analysis/Condition Pairwise Comparisons/Custom NVU Filtered DEGs/Volcano_Plots/"

if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir, recursive = TRUE)
}

plot_volcano <- function(df, cell_type, comparison,
                         lfc_cutoff = 1,
                         padj_cutoff = 0.05,
                         top_n = 10) {
  
  df <- df %>%
    filter(!is.na(padj)) %>%
    mutate(
      negLog10Padj = -log10(padj),
      Significance = case_when(
        padj < padj_cutoff & log2FoldChange > lfc_cutoff  ~ "Up",
        padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
        TRUE ~ "Not Sig"
      )
    )
  
  top_genes <- df %>%
    arrange(padj) %>%
    head(50)
  
  max_lfc <- max(abs(df$log2FoldChange), na.rm = TRUE)
  
  ggplot(df, aes(log2FoldChange, negLog10Padj, color = Significance)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c(
      "Up" = "#D62728",
      "Down" = "#1F77B4",
      "Not Sig" = "grey70"
    )) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    geom_text_repel(data = top_genes, aes(label = gene), size = 5) +
    coord_cartesian(xlim = c(-max_lfc, max_lfc)) +
    labs(
      title = paste(cell_type, "–", comparison),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    )
}

for (file in all_files) {
  
  df <- read.csv(file)
  
  filename <- basename(file)
  filename_no_ext <- sub("\\.csv$", "", filename)
  
  # Split only at FIRST underscore
  split_name <- strsplit(filename_no_ext, "_", fixed = TRUE)[[1]]
  
  cell_type <- split_name[1]
  
  # Reconstruct comparison name (everything after first underscore)
  comparison <- paste(split_name[-1], collapse = "_")
  
  # Clean comparison label
  comparison <- gsub("_NVU_sig", "", comparison)
  comparison <- gsub("_", " ", comparison)
  
  # Create cell-type folder
  cell_output_dir <- file.path(main_output_dir, cell_type)
  
  if (!dir.exists(cell_output_dir)) {
    dir.create(cell_output_dir)
  }
  
  # Generate plot
  p <- plot_volcano(df, cell_type, comparison)
  
  # Save plot
  ggsave(
    filename = file.path(cell_output_dir,
                         paste0("Volcano_", filename_no_ext, ".svg")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}











# Correlation matrices of NVU genes across metastases per cell type

# NVU Gene Correlation Heatmaps per Cell Type (Upper Triangle)

library(Seurat)
library(dplyr)
library(pheatmap)

# 1️⃣ NVU gene list
gene.list <- read.csv("Xenium Human Brain Panel Gene List-Santi.csv")
nvu_genes <- gene.list %>%
  filter(Santi.Annotations == "NVU") %>%
  pull(gene) %>%
  toupper()  # ensure consistency with pseudobulk gene names

# 2️⃣ Function: pseudobulk per cell type

obj <- TMA_full

Idents(obj) <- "sample"

obj <- RenameIdents(
  obj,
  "G2" = "Grade 2",
  "G3" = "Grade 3",
  "LU_NSC" = "NSCLC"
)

obj$sample <- Idents(obj)

pseudobulk_counts <- function(obj, celltype,
                              celltype_col="new_clusters",
                              sample_col="sample") {
  
  # Subset cells of this cell type
  cells_use <- rownames(obj@meta.data)[obj@meta.data[[celltype_col]] == celltype]
  sub <- subset(obj, cells = cells_use)
  
  counts <- GetAssayData(sub, layer = "counts")
  
  meta <- sub@meta.data[colnames(counts), , drop = FALSE]
  
  # Aggregate counts per sample
  pb <- t(rowsum(t(as.matrix(counts)), group = meta[[sample_col]]))
  
  return(pb)
}

# 3️⃣ Iterate over each cell type
cell_types <- unique(obj$new_clusters)  # adjust cluster column name if needed

for (ct in cell_types) {
  
  cat("Processing cell type:", ct, "\n")
  
  # Pseudobulk counts per sample
  pb <- pseudobulk_counts(obj, celltype = ct)
  
  # Select only NVU genes
  genes_use <- intersect(rownames(pb), nvu_genes)
  if (length(genes_use) == 0) {
    cat("No NVU genes found for", ct, "- skipping\n")
    next
  }
  
  pb_nvu <- pb[genes_use, ]
  
  # Optional: log2 transform + scaling across genes
  pb_norm <- log2(pb_nvu + 1)
  pb_scaled <- t(scale(t(pb_norm)))  # z-score per gene across samples
  
  
  # 4️⃣ Correlation matrix (samples/metastases)
  cor_matrix <- cor(pb_scaled, method = "pearson")
  
  # Define your output directory
  outdir <- "Manuscript Figures/Manuscript Figs (no LU_REG)/Correlation Heatmaps"  # e.g., "results/heatmaps"
  
  # Make sure the directory exists
  if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # 5️⃣ Heatmap
  # Construct the full file path for this cluster
  out_file <- file.path(outdir, paste0("NVU_Correlation_Heatmap_", gsub("/", "_", ct), ".png"))
  
  # Save the heatmap
  png(filename = out_file, width = 6, height = 6, units = "in", res = 300)
  p <- pheatmap(cor_matrix,
                color = colorRampPalette(c("blue", "white", "red"))(100),
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                main = paste0("NVU Gene Correlation: ", ct),
                fontsize_row = 16,
                fontsize_col = 16)
  dev.off()
  
}


#Pseudobulk gene expression matrices were generated for each cell type by summing counts per sample. Analysis was restricted to genes in the custom NVU panel. Counts were log2-transformed and scaled (z-score per gene across samples) to normalize expression patterns. Pearson correlation matrices were computed between all metastatic samples within each cell type, and visualized as heatmaps to compare relative NVU gene expression patterns across metastases.

conflicts_prefer(base::intersect)








#MFSD2A plots----


#VlnPlots of Select genes

TMA_EC.PC<- subset(TMA_full, subset = new_clusters == "EC/PC")

TMA_EC.PC$sample <- factor(TMA_EC.PC$sample, levels = c("Normal", "G2", "G3", "GBM", "BR", "COL", "END", "RCC", "LU_NSC", "MEL", "THY")) #"LU_REG"

DotPlot(TMA_EC.PC, features = c("MFSD2A"), group.by = "sample")

#Pseudo-bulk and wilcox
df <- FetchData(
  TMA_EC.PC,
  vars = c("MFSD2A", "sample", "sample_id", "condition")
)

agg <- df %>%
  group_by(sample, sample_id, condition) %>%
  summarise(
    mean_expr = mean(MFSD2A),
    pct_expr = mean(MFSD2A > 0),
    .groups = "drop"
  )

library(lme4)

fit <- lmer(mean_expr ~ condition + (1 | sample), data = agg)

library(lmerTest)
anova(fit)

library(emmeans)
emmeans(fit, pairwise ~ condition)

condition_colors <- c(
  "Normal" = "#2ca02c",
  "Primary Tumor" = "#1f77b4",
  "Secondary Tumor" = "#ff7f0e"
)

agg$Tumor <- ifelse(
  agg$condition == "Normal", "Normal",
  ifelse(agg$condition == "Primary Tumor", "Primary", "Metastatic")
)

agg$Tumor <- factor(agg$Tumor, levels = c("Normal", "Primary", "Metastatic"))

condition_colors <- c(
  "Normal" = "#2ca02c",
  "Primary" = "#1f77b4",
  "Metastatic" = "#ff7f0e"
)


ggplot(agg, aes(x = Tumor, y = mean_expr, color = Tumor)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  scale_color_manual(values = condition_colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),     # axis labels
    axis.text = element_text(size = 16),      # tick labels
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(
    x = "Condition",
    y = "MFSD2A Mean Expression",
    color = "Condition"
  )

ggplot(agg, aes(condition, pct_expr, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  theme_classic()


ggplot(agg, aes(sample, mean_expr, color = Tumor)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  scale_color_manual(values = condition_colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),     # axis labels
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 16), 
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(
    x = " ",
    y = "MFSD2A Mean Expression",
    color = "Condition"
  )




#seeing what cell types are epxressing MFSD2A

TMA_full <- readRDS("TMA_full_vascular_subtypes.rds")

TMA_full$sample <- factor(TMA_full$sample, levels = c("Normal", "G2", "G3", "GBM", "BR", "COL", "END", "RCC", "LU_NSC", "MEL", "THY")) # "LU_REG"
TMA_full$new_clusters <- factor(TMA_full$new_clusters, levels = rev(c("Neu-like", "Ast-like", "Immune", "EC", "Pericyte", "VLMC", "Mixed", "OPC/OL-like", "Tumor")))

# Extract normalized expression (slot="data") or counts (slot="counts")
expr_mat <- GetAssayData(TMA_full, assay = "Xenium", layer = "data")  # or "counts" for raw

MFSD2A_expr <- expr_mat["MFSD2A", ]  # vector, names = cell barcodes

# Combine with metadata
df <- TMA_full@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  left_join(
    data.frame(cell = names(MFSD2A_expr), MFSD2A = as.numeric(MFSD2A_expr)),
    by = "cell"
  )

agg <- df %>%
  group_by(sample, new_clusters) %>%
  summarise(
    sum_MFSD2A = sum(MFSD2A),
    .groups = "drop"
  )

agg <- agg %>%
  group_by(sample) %>%
  mutate(pct_MFSD2A = sum_MFSD2A / sum(sum_MFSD2A) * 100) %>%
  ungroup()


cell_colors <- c(
  "Neu-like" = "#1f78b4",    # blue
  "Tumor" = "gray",       # red
  "OPC/OL-like" = "orange", # green
  "EC/PC" = "red",       # orange
  "Immune" = "#6a3d9a",      # purple
  "Ast-like" = "darkgreen"
)

theme_santi.umap2<-function(base_size = 80, face ="bold", base_family= "Arial") {
  theme_few() %+replace%
    theme(axis.ticks = element_line(),
          strip.placement = "outside",
          panel.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(color = "black", size = 0.5),
          plot.title = element_blank(),
          axis.text = element_text(size = 28),
          axis.title = element_text(size = 34),
          legend.text = element_text(size = 36),
          legend.key.size = unit(2, "line")
    )
}

#1200 x 400
ggplot(agg, aes(x = sample, y = pct_MFSD2A, fill = new_clusters)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_colors) +
  theme_santi.umap2() +
  labs(
    x = " ",
    y = "MFSD2A expression percentage",
    fill = "Cell Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 18),     # axis labels
        axis.text = element_text(size = 16),      # tick labels
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) 




#dotplot of MFSD2A expression in celltypes in tma
library(Seurat)
library(dplyr)
library(ggplot2)

# Extract MFSD2A expression and metadata
mf_data <- FetchData(
  TMA_full,
  vars = c("MFSD2A", "new_clusters", "sample")
)

# Summarize within each cell type × sample
mf_summary <- mf_data %>%
  group_by(new_clusters, sample) %>%
  summarise(
    pct_MFSD2A_pos = mean(MFSD2A > 0) * 100,
    avg_MFSD2A = mean(MFSD2A),
    .groups = "drop"
  )

# Recode sample names
mf_summary$sample <- recode(
  mf_summary$sample,
  "G2" = "Grade 2",
  "G3" = "Grade 3"
)

# Order samples
mf_summary$sample <- factor(
  mf_summary$sample,
  levels = c(
    "Normal",
    "Grade 2",
    "Grade 3",
    "GBM",
    "BR",
    "COL",
    "END",
    "RCC",
    "NSCLC",
    "MEL",
    "THY"
  )
)

# Order cell types
mf_summary$new_clusters <- factor(
  mf_summary$new_clusters,
  levels = rev(c(
    "Neu-like",
    "Ast-like",
    "Immune",
    "EC",
    "Pericyte",
    "VLMC",
    "Mixed",
    "OPC/OL-like",
    "Tumor"
  ))
)

# Dot heatmap
mfsd2a_dot<- ggplot(
  mf_summary,
  aes(
    x = sample,
    y = new_clusters
  )
) +
  geom_point(
    aes(
      size = pct_MFSD2A_pos,
      color = avg_MFSD2A
    )
  ) +
  scale_size(
    range = c(0, 10),
    name = "MFSD2A+ cells (%)"
  ) +
  scale_color_gradient(
    low = "white",
    high = "darkblue",
    name = "Mean MFSD2A\ntranscripts/cell"
  ) +
  labs(
    x = " ",
    y = " "
  ) +
  ggtitle("MFSD2A Expression in Tumor Core") +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


# Set directory where you want to save
output_dir <- "Vascular subtype analysis/Figures/"

# Make sure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Export UMAP as SVG
ggsave(
  filename = "Vascular subtype analysis/Figures/tma_mfsd2a_dot.svg",
  plot = mfsd2a_dot,
  width = 9.5,
  height = 5
)



#Volcanoes of genes of interest

TMA_EC.PC

TMA_EC.PC$sample <- factor(TMA_EC.PC$sample, levels = c("Normal", "G2", "G3", "GBM", "BR", "COL", "END", "RCC", "LU_NSC", "MEL", "THY")) # "LU_REG"

Idents(TMA_EC.PC) <- "sample"

TMA_EC.PC <- RenameIdents(
  TMA_EC.PC,
  "G2" = "Grade 2",
  "G3" = "Grade 3",
  "LU_NSC" = "NSCLC"
)

TMA_EC.PC$sample <- Idents(TMA_EC.PC)

TMA_EC.PC$sample <- factor(TMA_EC.PC$sample, levels = c("Normal", "Grade 2", "Grade 3", "GBM", "BR", "COL", "END", "RCC", "NSCLC", "MEL", "THY")) # "LU_REG"



VlnPlot(TMA_EC.PC, group.by = "sample", features = c("PLVAP", "COL4A1", "MFSD2A", "CLDN5", "ABCB1", "ABCG2"), pt.size = 0.1)

#export 
VlnPlot(
  TMA_EC.PC,
  group.by = "sample",
  features = c("PLVAP", "COL4A1", "MFSD2A", "CLDN5", "ABCB1", "ABCG2"),
  pt.size = 0.1,
  ncol = 3   # 👈 number of columns
)







#Ast-like volcano plots


input_dir <- "TMA6 & TMA8 Met/More Filtered analysis/Condition Pairwise Comparisons_noLuReg/Filtered Outputs/"
all_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
#just keep ast-like
all_files <- list.files(
  input_dir,
  pattern = "^Ast-like_.*\\.csv$",
  full.names = TRUE
)

#Load data from Human Protein Atlas
hpa_astro <- read_tsv("Manuscript Figures/Manuscript Figs (no LU_REG)/Ast_Volcano_Plots/rna_single_cell_type.tsv") %>% 
  filter(`Cell type` == "astrocytes") %>% 
  filter(nCPM > 5)

#now filter by genes in the NVU custom + ast markers
nvu_ast_genes <- gene.list %>% 
  filter(Santi.Annotations %in% c("NVU", "Astrocyte"))

xenium_genes <- nvu_ast_genes$gene

astro_expressed <- xenium_genes[xenium_genes %in% hpa_astro$`Gene name`]




main_output_dir <- "Manuscript Figures/Manuscript Figs (no LU_REG)/Ast_Volcano_Plots/"

if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir, recursive = TRUE)
}

plot_volcano <- function(df, cell_type, comparison,
                         lfc_cutoff = 1,
                         padj_cutoff = 0.05,
                         top_n = 10) {
  
  df <- df %>%
    filter(!is.na(padj)) %>%
    mutate(
      negLog10Padj = -log10(padj),
      Significance = case_when(
        padj < padj_cutoff & log2FoldChange > lfc_cutoff  ~ "Up",
        padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
        TRUE ~ "Not Sig"
      )
    )
  
  top_genes <- df %>%
    arrange(padj) %>%
    head(50)
  
  max_lfc <- max(abs(df$log2FoldChange), na.rm = TRUE)
  
  ggplot(df, aes(log2FoldChange, negLog10Padj, color = Significance)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c(
      "Up" = "#D62728",
      "Down" = "#1F77B4",
      "Not Sig" = "grey70"
    )) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    geom_text_repel(data = top_genes, aes(label = gene), size = 5) +
    coord_cartesian(xlim = c(-max_lfc, max_lfc)) +
    labs(
      title = paste(cell_type, "–", comparison),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14)
    )
}

for (file in all_files) {
  
  df <- read.csv(file)
  
  # Filter DEGs to only astrocyte-expressed genes
  df <- df %>% filter(gene %in% astro_expressed)
  
  filename <- basename(file)
  filename_no_ext <- sub("\\.csv$", "", filename)
  
  # Split only at FIRST underscore
  split_name <- strsplit(filename_no_ext, "_", fixed = TRUE)[[1]]
  
  cell_type <- split_name[1]
  
  # Reconstruct comparison name (everything after first underscore)
  comparison <- paste(split_name[-1], collapse = "_")
  
  # Clean comparison label
  comparison <- gsub("_NVU_sig", "", comparison)
  comparison <- gsub("_", " ", comparison)
  
  # Create cell-type folder
  cell_output_dir <- file.path(main_output_dir, cell_type)
  if (!dir.exists(cell_output_dir)) {
    dir.create(cell_output_dir)
  }
  
  # Generate plot with filtered genes
  p <- plot_volcano(df, cell_type, comparison)
  
  # Save plot
  ggsave(
    filename = file.path(cell_output_dir,
                         paste0("Volcano_", filename_no_ext, ".svg")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}






#Tae Hoon Genes of Interest Genes

#load genes
install.packages("readxl")
library(readxl)
T_H_genes <- read_excel("Tae Hoon Data/TMA8 EC Custom Candidate Gene List 2025 12.xlsx")

genes <- T_H_genes$`Gene List`

#Pipeline
run_pipeline <- function(gene, seurat_obj) {
  
  message("Processing: ", gene)
  
  # Subset
  sub_obj <- subset(seurat_obj, subset = new_clusters == "EC/PC")
  
  # Fetch data
  df <- FetchData(
    sub_obj,
    vars = c(gene, "sample", "sample_id", "condition")
  )
  
  colnames(df)[1] <- "expr"
  
  # Aggregate
  agg <- df %>%
    group_by(sample, sample_id, condition) %>%
    summarise(
      mean_expr = mean(expr, na.rm = TRUE),
      pct_expr = mean(expr > 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Tumor labeling
  agg$Tumor <- ifelse(
    agg$condition == "Normal", "Normal",
    ifelse(agg$condition == "Primary Tumor", "Primary", "Metastatic")
  )
  
  agg$Tumor <- factor(agg$Tumor, levels = c("Normal", "Primary", "Metastatic"))
  
  # Mixed model
  fit <- lmer(mean_expr ~ condition + (1 | sample), data = agg)
  
  stats <- list(
    anova = anova(fit),
    emmeans = emmeans(fit, pairwise ~ condition)
  )
  
  # Plot 1: Tumor groups
  p1 <- ggplot(agg, aes(x = Tumor, y = mean_expr, color = Tumor)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2) +
    theme_classic() +
    ggtitle(paste0(gene, " Mean Expression"))
  
  # Plot 2: mean expression per sample
  p2 <-  ggplot(agg, aes(sample, mean_expr, color = Tumor)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 2) + 
    scale_color_manual(values = condition_colors) + theme_classic() + 
    theme( 
      legend.position = "none", 
      axis.title = element_text(size = 18),
      axis.text.x = element_text(size = 16, angle = 45, hjust = 1), 
      axis.text.y = element_text(size = 16), plot.title = element_text(size = 16, face = "bold") 
    ) + 
    labs( 
      x = " ",
      y = " Mean Expression",
      color = "Condition"
    ) +
    ggtitle(paste0(gene, " Mean Expression"))
  
  return(list(
    gene = gene,
    stats = stats,
    agg = agg,
    plots = list(p1 = p1, p2 = p2)
  ))
}



#run all genes
results <- lapply(genes, run_pipeline, seurat_obj = TMA_EC.PC)

names(results) <- genes



#save plots
out_dir <- "T_H_Genes"

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

for (res in results) {
  if (is.null(res)) next  # skip failed genes
  
  g <- res$gene
  
  ggsave(
    filename = file.path(out_dir, paste0(g, "_mean_expr_condition.png")),
    plot = res$plots$p1,
    width = 5,
    height = 4
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(g, "_mean_expr_sample.png")),
    plot = res$plots$p2,
    width = 5,
    height = 4
  )
}







#celltype composition

# Get expression matrix ONCE (important for speed)
expr_mat <- GetAssayData(TMA_full, assay = "Xenium", layer = "data")

# Clean genes
genes <- trimws(T_H_genes$`Gene List`)

for (gene in genes) {
  
  message("Processing: ", gene)
  
  # Skip if gene not found
  if (!gene %in% rownames(expr_mat)) {
    message("Skipping ", gene, " (not found)")
    next
  }
  
  # Extract expression
  gene_expr <- expr_mat[gene, ]
  gene_expr <- as.numeric(gene_expr)
  names(gene_expr) <- colnames(expr_mat)
  
  # Attach to metadata (robust, no join)
  meta <- TMA_full@meta.data
  meta$expr <- gene_expr[rownames(meta)]
  
  # Skip if no expression
  if (all(is.na(meta$expr)) || sum(meta$expr, na.rm = TRUE) == 0) {
    message("Skipping ", gene, " (no expression)")
    next
  }
  
  # Aggregate
  agg_ct <- meta %>%
    tibble::rownames_to_column("cell") %>%
    group_by(sample, new_clusters) %>%
    summarise(
      sum_expr = sum(expr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(sample) %>%
    mutate(
      total = sum(sum_expr, na.rm = TRUE),
      pct_expr = ifelse(total > 0, sum_expr / total * 100, 0)
    ) %>%
    ungroup()
  
  # Plot
  p <- ggplot(agg_ct, aes(x = sample, y = pct_expr, fill = new_clusters)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cell_colors) +
    theme_santi.umap2() +
    labs(
      x = " ",
      y = paste0(gene, " expression (%)"),
      fill = "Cell Type"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16)
    )
  
  # Save
  
  out_dir <- "T_H_Genes_2"
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  
  ggsave(
    filename = file.path(out_dir, paste0(gene, "_celltype_composition.png")),
    plot = p,
    width = 8,
    height = 4
  )
}



# TMA EC NVU Remodeling Continuum----

# TMA EC NVU Remodeling Continuum
# Complete analysis and supplementary figure

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(UCell)
library(grid)
library(viridis)
library(patchwork)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("slingshot")

library(slingshot)


# Subset endothelial cells from Normal, Grade 2, Grade 3, and GBM

TMA_full_vascular_subtypes <- readRDS(
  "TMA_full_vascular_subtypes.rds"
)

TMA_EC_primary <- subset(
  TMA_full_vascular_subtypes,
  subset =
    new_clusters == "EC" &
    sample %in% c(
      "Normal",
      "G2",
      "G3",
      "GBM"
    )
)

# Define disease-state labels

TMA_EC_primary$sample_plot <- factor(
  dplyr::recode(
    as.character(TMA_EC_primary$sample),
    "Normal" = "Normal",
    "G2" = "Grade 2",
    "G3" = "Grade 3",
    "GBM" = "GBM"
  ),
  levels = c(
    "Normal",
    "Grade 2",
    "Grade 3",
    "GBM"
  )
)

TMA_EC_primary$disease_stage <- dplyr::recode(
  as.character(TMA_EC_primary$sample),
  "Normal" = 0,
  "G2" = 1,
  "G3" = 2,
  "GBM" = 3
)

# Define NVU and BBB modules

nvu_modules <- list(
  Tight_Junctions = c(
    "CLDN5",
    "OCLN",
    "TJP1",
    "TJP2",
    "TJP3",
    "F11R",
    "JAM2",
    "JAM3",
    "LSR"
  ),
  
  BBB_Transport = c(
    "ABCB1",
    "ABCG2",
    "SLC2A1",
    "SLCO1A2",
    "SLCO1C1",
    "LRP1",
    "MFSD2A"
  ),
  
  Angiogenesis = c(
    "ANGPT2",
    "VEGFA",
    "DLL4",
    "NDP",
    "KDR",
    "TEK"
  ),
  
  Vascular_Inflam. = c(
    "ICAM1",
    "VCAM1",
    "IL1B",
    "IL6",
    "CCL2",
    "S1PR3"
  ),
  
  ECM_Remodeling = c(
    "COL4A1",
    "LAMB1",
    "FN1",
    "VCL",
    "ITGB1",
    "MMP2",
    "MMP9",
    "VIM"
  )
)

# Calculate UCell module scores

TMA_EC_primary <- AddModuleScore_UCell(
  TMA_EC_primary,
  features = nvu_modules
)

# Calculate BBB Stability Score

TMA_EC_primary$BBB_Stability <-
  (
    TMA_EC_primary$Tight_Junctions_UCell +
      TMA_EC_primary$BBB_Transport_UCell
  ) -
  (
    TMA_EC_primary$Angiogenesis_UCell +
      TMA_EC_primary$Vascular_Inflam._UCell +
      TMA_EC_primary$ECM_Remodeling_UCell
  )

# Define genes used for trajectory construction

genes_use <- unique(
  unlist(nvu_modules)
)

genes_use <- genes_use[
  genes_use %in%
    rownames(TMA_EC_primary)
]

# Scale NVU genes

TMA_EC_primary <- ScaleData(
  TMA_EC_primary,
  features = genes_use,
  verbose = FALSE
)

# PCA using NVU/BBB genes

TMA_EC_primary <- RunPCA(
  TMA_EC_primary,
  features = genes_use,
  verbose = FALSE
)

# UMAP

TMA_EC_primary <- RunUMAP(
  TMA_EC_primary,
  dims = 1:10,
  reduction = "pca",
  verbose = FALSE
)

# Cluster endothelial cells for Slingshot

TMA_EC_primary <- FindNeighbors(
  TMA_EC_primary,
  dims = 1:10,
  verbose = FALSE
)

TMA_EC_primary <- FindClusters(
  TMA_EC_primary,
  resolution = 0.3,
  verbose = FALSE
)

cluster_assignments <- Idents(
  TMA_EC_primary
)

# Identify the cluster containing the largest number of Normal ECs

cluster_sample_table <- table(
  cluster_assignments,
  TMA_EC_primary$sample_plot
)

normal_cluster <- rownames(
  cluster_sample_table
)[
  which.max(
    cluster_sample_table[, "Normal"]
  )
]

print(normal_cluster)

# Run Slingshot

umap <- Embeddings(
  TMA_EC_primary,
  "umap"
)

sce <- slingshot(
  umap,
  clusterLabels = cluster_assignments,
  start.clus = normal_cluster
)

# Extract NVU Remodeling Continuum

TMA_EC_primary$NVU_Remodeling_Continuum <-
  slingPseudotime(sce)[, 1]

# Extract Slingshot trajectory curve

curve <- slingCurves(sce)[[1]]

curve_df <- data.frame(
  UMAP_1 = curve$s[
    curve$ord,
    1
  ],
  UMAP_2 = curve$s[
    curve$ord,
    2
  ]
)

# Create directional arrow using final portion of trajectory

n_curve <- nrow(
  curve_df
)

arrow_end_index <- n_curve

arrow_start_index <- max(
  1,
  n_curve - max(
    5,
    round(
      n_curve * 0.08
    )
  )
)

sling_arrow_df <- data.frame(
  x_start = curve_df$UMAP_1[
    arrow_start_index
  ],
  y_start = curve_df$UMAP_2[
    arrow_start_index
  ],
  x_end = curve_df$UMAP_1[
    arrow_end_index
  ],
  y_end = curve_df$UMAP_2[
    arrow_end_index
  ]
)

# Prepare sample-level dataset
# Each sample_id is one independent biological replicate

continuum_sample_df <- TMA_EC_primary@meta.data %>%
  dplyr::filter(
    !is.na(
      NVU_Remodeling_Continuum
    ),
    !is.na(sample_id),
    !is.na(sample_plot)
  ) %>%
  dplyr::group_by(
    sample_id,
    sample_plot,
    disease_stage
  ) %>%
  dplyr::summarise(
    NVU_Remodeling_Continuum = median(
      NVU_Remodeling_Continuum,
      na.rm = TRUE
    ),
    Tight_Junctions = median(
      Tight_Junctions_UCell,
      na.rm = TRUE
    ),
    BBB_Transport = median(
      BBB_Transport_UCell,
      na.rm = TRUE
    ),
    Angiogenesis = median(
      Angiogenesis_UCell,
      na.rm = TRUE
    ),
    Vascular_Inflam. = median(
      Vascular_Inflam._UCell,
      na.rm = TRUE
    ),
    ECM_Remodeling = median(
      ECM_Remodeling_UCell,
      na.rm = TRUE
    ),
    BBB_Stability = median(
      BBB_Stability,
      na.rm = TRUE
    ),
    n_EC = dplyr::n(),
    .groups = "drop"
  )

# Confirm independent sample numbers

print(
  table(
    continuum_sample_df$sample_plot
  )
)

# Disease-stage association with continuum

continuum_kw <- kruskal.test(
  NVU_Remodeling_Continuum ~ sample_plot,
  data = continuum_sample_df
)

continuum_stage_cor <- cor.test(
  continuum_sample_df$disease_stage,
  continuum_sample_df$NVU_Remodeling_Continuum,
  method = "spearman",
  exact = FALSE
)

# Pairwise disease-state comparisons

pairwise_continuum <- pairwise.wilcox.test(
  continuum_sample_df$NVU_Remodeling_Continuum,
  continuum_sample_df$sample_plot,
  p.adjust.method = "BH"
)

# Correlations between continuum and NVU modules

module_names <- c(
  "Tight_Junctions",
  "BBB_Transport",
  "Angiogenesis",
  "Vascular_Inflam.",
  "ECM_Remodeling",
  "BBB_Stability"
)

module_labels <- c(
  "Tight Junctions",
  "BBB Transport",
  "Angiogenesis",
  "Vascular Inflammation",
  "ECM Remodeling",
  "BBB Stability"
)

module_results <- lapply(
  module_names,
  function(module) {
    
    test <- cor.test(
      continuum_sample_df$NVU_Remodeling_Continuum,
      continuum_sample_df[[module]],
      method = "spearman",
      exact = FALSE
    )
    
    data.frame(
      Feature = module,
      Spearman_rho = unname(
        test$estimate
      ),
      P_value = test$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(
      P_value,
      method = "BH"
    )
  ) %>%
  arrange(P_value)

# Disease-stage correlations with each module

disease_order_results <- lapply(
  module_names,
  function(module) {
    
    test <- cor.test(
      continuum_sample_df$disease_stage,
      continuum_sample_df[[module]],
      method = "spearman",
      exact = FALSE
    )
    
    data.frame(
      Feature = module,
      Spearman_rho = unname(
        test$estimate
      ),
      P_value = test$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(
      P_value,
      method = "BH"
    )
  ) %>%
  arrange(P_value)

# Define disease-state colors

condition_colors <- c(
  "Normal" = "grey40",
  "Grade 2" = "#1B9E77",
  "Grade 3" = "#D95F02",
  "GBM" = "#7570B3"
)

# Extract UMAP coordinates

umap_df <- as.data.frame(
  Embeddings(
    TMA_EC_primary,
    "umap"
  )
)

colnames(umap_df)[1:2] <- c(
  "UMAP_1",
  "UMAP_2"
)

umap_df$sample_plot <-
  TMA_EC_primary$sample_plot

umap_df$NVU_Remodeling_Continuum <-
  TMA_EC_primary$NVU_Remodeling_Continuum

# Panel A
# UMAP colored by disease state

pA <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = sample_plot
  )
) +
  geom_point(
    size = 1,
    alpha = 0.7
  ) +
  scale_color_manual(
    values = condition_colors,
    name = "Disease state"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  ) +
  labs(
    title = "A",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(
    base_size = 16
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      face = "bold"
    ),
    axis.title = element_text(
      size = 16
    ),
    axis.text = element_text(
      size = 13
    ),
    legend.title = element_text(
      size = 14
    ),
    legend.text = element_text(
      size = 13
    )
  )

# Panel B
# UMAP colored by NVU Remodeling Continuum
# with actual Slingshot trajectory

pB <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = NVU_Remodeling_Continuum
  )
) +
  geom_point(
    size = 0.1,
    alpha = 0.7
  ) +
  scale_color_viridis_c(
    option = "plasma",
    name = "NVU Remodeling\nContinuum"
  ) +
  geom_path(
    data = curve_df,
    aes(
      x = UMAP_1,
      y = UMAP_2
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.2
  ) +
  geom_segment(
    data = sling_arrow_df,
    aes(
      x = x_start,
      y = y_start,
      xend = x_end,
      yend = y_end
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.2,
    arrow = arrow(
      length = unit(
        0.25,
        "cm"
      ),
      type = "closed"
    )
  ) +
  labs(
    title = "B",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(
    base_size = 16
  ) +
  theme(
    plot.title = element_text(
      size = 20,
      face = "bold"
    ),
    axis.title = element_text(
      size = 16
    ),
    axis.text = element_text(
      size = 13
    ),
    legend.title = element_text(
      size = 14
    ),
    legend.text = element_text(
      size = 13
    )
  )


# Panel C
# Sample-level NVU Remodeling Continuum across disease states

continuum_sample_df$sample_plot <- factor(
  continuum_sample_df$sample_plot,
  levels = c(
    "GBM",
    "Grade 3",
    "Grade 2",
    "Normal"
  )
)

pC <- ggplot(
  continuum_sample_df,
  aes(
    x = NVU_Remodeling_Continuum,
    y = sample_plot
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.6
  ) +
  geom_quasirandom(
    aes(
      group = sample_plot
    ),
    width = 0.18,
    size = 2.8,
    alpha = 0.85
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  labs(
    x = "NVU Remodeling Continuum",
    y = NULL
  ) +
  annotate(
    "text",
    x = Inf,
    y = 4,
    label = paste0(
      "Spearman \u03c1 = ",
      round(
        unname(
          continuum_stage_cor$estimate
        ),
        2
      ),
      "\nP = ",
      format.pval(
        continuum_stage_cor$p.value,
        digits = 2,
        eps = 0.001
      )
    ),
    hjust = 1.05,
    vjust = 1.2,
    size = 4.5
  ) +
  theme_classic(
    base_size = 16
  ) +
  theme(
    axis.title.x = element_text(
      size = 17
    ),
    axis.text.x = element_text(
      size = 14
    ),
    axis.text.y = element_text(
      size = 14
    )
  )

# Prepare module plotting dataframe

module_plot_df <- continuum_sample_df %>%
  dplyr::select(
    sample_id,
    sample_plot,
    NVU_Remodeling_Continuum,
    all_of(module_names)
  ) %>%
  pivot_longer(
    cols = all_of(module_names),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  mutate(
    Module = factor(
      Module,
      levels = module_names,
      labels = module_labels
    )
  )

# Calculate statistics for Panel D

module_stats <- lapply(
  module_names,
  function(module) {
    
    test <- cor.test(
      continuum_sample_df$NVU_Remodeling_Continuum,
      continuum_sample_df[[module]],
      method = "spearman",
      exact = FALSE
    )
    
    data.frame(
      Module = module,
      rho = unname(
        test$estimate
      ),
      P = test$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(
      P,
      method = "BH"
    ),
    Module = factor(
      Module,
      levels = module_names,
      labels = module_labels
    )
  )

module_plot_df <- module_plot_df %>%
  left_join(
    module_stats,
    by = "Module"
  )

module_annotations <- module_plot_df %>%
  group_by(Module) %>%
  summarise(
    rho = dplyr::first(rho),
    FDR = dplyr::first(FDR),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "\u03c1 = ",
      sprintf(
        "%.2f",
        rho
      ),
      "\nFDR = ",
      format.pval(
        FDR,
        digits = 2,
        eps = 0.001
      )
    )
  )

# Panel D
# Individual NVU modules versus continuum

pD <- ggplot(
  module_plot_df,
  aes(
    x = NVU_Remodeling_Continuum,
    y = Score
  )
) +
  geom_point(
    aes(
      color = sample_plot
    ),
    size = 2.3,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 0.6
  ) +
  facet_wrap(
    ~Module,
    scales = "free_y",
    ncol = 3
  ) +
  geom_text(
    data = module_annotations,
    aes(
      x = Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.2,
    size = 3.8
  ) +
  scale_color_manual(
    values = condition_colors
  ) +
  labs(
    x = "NVU Remodeling Continuum",
    y = "Module Score",
    color = NULL
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    axis.title = element_text(
      size = 15
    ),
    axis.text = element_text(
      size = 11
    ),
    legend.text = element_text(
      size = 12
    ),
    legend.position = "bottom",
    panel.spacing = unit(
      1.1,
      "lines"
    ),
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )

# Combine panels

top_row <- pA +
  pB +
  pC +
  patchwork::plot_layout(
    widths = c(
      1,
      1,
      1
    )
  )

p_NVU_continuum <- (
  top_row /
    pD
) +
  patchwork::plot_layout(
    heights = c(
      1,
      1.15
    )
  )

print(p_NVU_continuum)

# Save supplementary figure

# Define output directory

# Set directory where you want to save
output_dir <- "Vascular subtype analysis/Figures/"

# Make sure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Export UMAP as SVG
ggsave(
  filename = "Vascular subtype analysis/Figures/Trajectory_TMAs.svg",
  plot = p_NVU_continuum,
  width = 24,
  height = 24
)


# Save statistical results

write.csv(
  continuum_sample_df,
  "NVU_Remodeling_Continuum_sample_level_data.csv",
  row.names = FALSE
)

write.csv(
  module_results,
  "NVU_Remodeling_Continuum_module_correlations.csv",
  row.names = FALSE
)

write.csv(
  disease_order_results,
  "NVU_Remodeling_Continuum_disease_ordering.csv",
  row.names = FALSE
)

# Print key statistics

cat("\n")
cat("NVU Remodeling Continuum statistics\n")
cat("\n")

cat(
  "Kruskal-Wallis P = ",
  format.pval(
    continuum_kw$p.value,
    digits = 3,
    eps = 0.001
  ),
  "\n"
)

cat(
  "Disease-stage Spearman rho = ",
  round(
    unname(
      continuum_stage_cor$estimate
    ),
    3
  ),
  "\n"
)

cat(
  "Disease-stage Spearman P = ",
  format.pval(
    continuum_stage_cor$p.value,
    digits = 3,
    eps = 0.001
  ),
  "\n"
)

print(module_results)


#Exporting each plot separately
# Export individual supplementary panels as SVG

output_dir <- "Vascular subtype analysis/Figures"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Panel A
ggsave(
  filename = file.path(
    output_dir,
    "Trajectory_TMAs_Panel_A.svg"
  ),
  plot = pA,
  width = 5,
  height = 4,
  units = "in"
)

# Panel B
ggsave(
  filename = file.path(
    output_dir,
    "Trajectory_TMAs_Panel_B.svg"
  ),
  plot = pB,
  width = 5,
  height = 4,
  units = "in"
)

# Panel C
ggsave(
  filename = file.path(
    output_dir,
    "Trajectory_TMAs_Panel_C.svg"
  ),
  plot = pC,
  width = 10,
  height = 4,
  units = "in"
)

# Panel D
ggsave(
  filename = file.path(
    output_dir,
    "Trajectory_TMAs_Panel_D.svg"
  ),
  plot = pD,
  width = 12,
  height = 7,
  units = "in"
)

