#Unified framework of BBB Disruption

#Use TMA & Big Tissue Data to show module scores
#Use scRNA seq data to show MAST analysis and GSEA pathways involved in BBB Disruption (High BBB Stability --> Low BBB Stability) for LGG, GBM, Mets

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



#TMA Data
TMA_full_EC.PC <- read_rds("TMA_full_EC.PC_modules.rds")

#New Categories
TMA_full_EC.PC$disease_category <- dplyr::case_when(
  TMA_full_EC.PC$sample == "Normal" ~ "Non-tumor",
  TMA_full_EC.PC$sample %in% c("G2", "G3") ~ "LGG",
  TMA_full_EC.PC$sample == "GBM" ~ "GBM",
  TMA_full_EC.PC$sample %in% c("BR", "COL", "END", "LU_NSC", "MEL", "RCC", "THY") ~ "BrM",
  TRUE ~ NA_character_
)

table(
  TMA_full_EC.PC$disease_category,
  TMA_full_EC.PC$sample
)

#Subcluster EC/PC cluster
DimPlot(TMA_full_EC.PC)

## Re-run dimensional reduction and clustering
## on EC/PC compartment

DefaultAssay(TMA_full_EC.PC) <- "Xenium"

# Normalize if necessary
TMA_full_EC.PC <- NormalizeData(
  TMA_full_EC.PC,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# Variable features
TMA_full_EC.PC <- FindVariableFeatures(
  TMA_full_EC.PC,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

# Scale
TMA_full_EC.PC <- ScaleData(
  TMA_full_EC.PC,
  features = VariableFeatures(TMA_full_EC.PC),
  verbose = FALSE
)

# PCA
TMA_full_EC.PC <- RunPCA(
  TMA_full_EC.PC,
  features = VariableFeatures(TMA_full_EC.PC),
  npcs = 30,
  verbose = FALSE
)

# Inspect PCA
ElbowPlot(TMA_full_EC.PC, ndims = 30)

# Neighbors
TMA_full_EC.PC <- FindNeighbors(
  TMA_full_EC.PC,
  dims = 1:15
)

# Clusters
TMA_full_EC.PC <- FindClusters(
  TMA_full_EC.PC,
  resolution = 0.3
)

# UMAP
TMA_full_EC.PC <- RunUMAP(
  TMA_full_EC.PC,
  dims = 1:15
)

# Plot
DimPlot(
  TMA_full_EC.PC,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
)

#Markers
## Find marker genes for EC/PC subclusters

TMA_ECPC_markers <- FindAllMarkers(
  TMA_full_EC.PC,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25
)

# Top markers per cluster
TMA_ECPC_top_markers <- TMA_ECPC_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20
  )

TMA_ECPC_top_markers

## Vascular identity marker panel
## EC vs Pericyte vs VSMC vs Mixed

vascular_markers <- list(
  
  # Endothelial cells
  EC = c(
    "PECAM1",
    "CDH5",
    "CLDN5",
    "VWF",
    "FLT1",
    "KDR",
    "TEK",
    "ESAM",
    "PLVAP",
    "ENG",
    "NOS3",
    "GJA4",
    "NRP1",
    "NOTCH1",
    "ADGRL4",
    "PODXL"
  ),
  
  # Pericytes
  Pericyte = c(
    "PDGFRB",
    "RGS5",
    "NOTCH3",
    "KCNJ8",
    "ABCC9",
    "CSPG4"
  ),
  
  # VLMC / ECM-associated vascular cells
  VLMC = c(
    "DCN",
    "COL12A1",
    "LAMA2",
    "FBLN1"
  )
)

# Keep only genes present in the object
vascular_markers <- lapply(
  vascular_markers,
  function(x) intersect(x, rownames(TMA_full_EC.PC))
)

vascular_markers

## Vascular identity DotPlot

DotPlot(
  TMA_full_EC.PC,
  group.by = "Xenium_snn_res.0.3",
  features = vascular_markers,
  dot.scale = 6
) +
  RotatedAxis()


#reannonate subclusters
## Rename vascular clusters

vascular_annotation <- c(
  "0" = "Pericyte",
  "1" = "EC",
  "2" = "VLMC",
  "3" = "Mixed",
  "4" = "Mixed",
  "5" = "EC",
  "6" = "VLMC",
  "7" = "EC",
  "8" = "Mixed"
)

# Set identities to the clustering at resolution 0.3
Idents(TMA_full_EC.PC) <- "Xenium_snn_res.0.3"

# Rename clusters
TMA_full_EC.PC <- RenameIdents(
  TMA_full_EC.PC,
  `0` = "Pericyte",
  `1` = "EC",
  `2` = "VLMC",
  `3` = "Mixed",
  `4` = "Mixed",
  `5` = "EC",
  `6` = "VLMC",
  `7` = "EC",
  `8` = "Mixed"
)

# Check
table(Idents(TMA_full_EC.PC))

TMA_full_EC.PC$vasc_sub <- Idents(TMA_full_EC.PC)


# Visualize
DimPlot(
  TMA_full_EC.PC,
  label = TRUE,
  repel = TRUE
)

write_rds(TMA_full_EC.PC, "TMA_full_EC.PC_subcl.EC.rds")






TMA_full_EC.PC_subcl.EC <- readRDS("TMA_full_EC.PC_subcl.EC.rds")

#Subset new EC subcluster

TMA_EC_sub <- subset(TMA_full_EC.PC_subcl.EC,
                     vasc_sub == "EC")


# Extract metadata
TMA_meta <- TMA_EC_sub@meta.data

# Define modules
nvu_heatmap_modules <- c(
  "Tight_Junctions",
  "BBB_Transport",
  "Angiogenesis",
  "Vascular_Inflam.",
  "ECM_Remodeling"
)


# Step 1: Average EC cells within each patient/sample
TMA_patient_summary <- TMA_meta %>%
  group_by(
    sample_id,
    sample,
    disease_category
  ) %>%
  summarise(
    across(
      all_of(nvu_heatmap_modules),
      ~ mean(.x, na.rm = TRUE)
    ),
    BBB_Stability_Score = mean(
      BBB_Stability_Score,
      na.rm = TRUE
    ),
    n_cells = n(),
    .groups = "drop"
  )


# Check patient-level summary
head(TMA_patient_summary)


# Step 2: Average patients within disease category
TMA_disease_summary <- TMA_patient_summary %>%
  group_by(
    disease_category
  ) %>%
  summarise(
    across(
      all_of(nvu_heatmap_modules),
      ~ mean(.x, na.rm = TRUE)
    ),
    BBB_Stability_Score = mean(
      BBB_Stability_Score,
      na.rm = TRUE
    ),
    n_patients = n(),
    .groups = "drop"
  )


# View final heatmap values
TMA_disease_summary

#Prepare for heatmap
TMA_disease_summary <- TMA_disease_summary %>%
  mutate(
    Dataset = "TMA",
    Disease = disease_category
  )





#Large Tissue
whole_tissue <- readRDS("NVU_BBB_annotated_Seurat_scripts01to07.rds")

ec_whole_tissue <- subset(whole_tissue,
                          subset = vascular_subtype == "EC")



ec_whole_tissue$disease_category <- dplyr::case_when(
  ec_whole_tissue$sample == "Normal" ~ "Non-tumor",
  ec_whole_tissue$sample == "Gr_3" ~ "LGG",
  ec_whole_tissue$sample == "GBM" ~ "GBM",
  ec_whole_tissue$sample %in% c(
    "Breast",
    "Kidney",
    "Lung",
    "Melanoma",
    "Thyroid"
  ) ~ "BrM",
  TRUE ~ NA_character_
)


table(
  ec_whole_tissue$disease_category,
  ec_whole_tissue$sample
)


# Extract metadata
large_meta <- ec_whole_tissue@meta.data


# Define modules
nvu_heatmap_modules <- c(
  "Tight_Junctions_UCell",
  "BBB_Transport_UCell",
  "Angiogenesis_UCell",
  "Vascular_Inflam_UCell",
  "ECM_Remodeling_UCell"
)


# Step 1:
# Average EC/PC cells within each tissue specimen
whole_tissue_patient_summary <- large_meta %>%
  dplyr::group_by(
    sample,
    disease_category
  ) %>%
  dplyr::summarise(
    across(
      all_of(nvu_heatmap_modules),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_cells = n(),
    .groups = "drop"
  )

# Check number of biological replicates

whole_tissue_patient_summary %>%
  dplyr::count(disease_category)


whole_tissue_disease_summary <- whole_tissue_patient_summary %>%
  dplyr::group_by(
    disease_category
  ) %>%
  dplyr::summarise(
    across(
      all_of(nvu_heatmap_modules),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_samples = n(),
    .groups = "drop"
  )

#rename columns to take out UCell
whole_tissue_disease_summary <- whole_tissue_disease_summary %>%
  dplyr::rename(
    Tight_Junctions = Tight_Junctions_UCell,
    BBB_Transport = BBB_Transport_UCell,
    Angiogenesis = Angiogenesis_UCell,
    Vascular_Inflam. = Vascular_Inflam_UCell,
    ECM_Remodeling = ECM_Remodeling_UCell
  )

#Add Dataset label
whole_tissue_disease_summary <- whole_tissue_disease_summary %>%
  mutate(
    Dataset = "Xenium_whole_tissue",
    Disease = disease_category
  )









#Bejarano scRNA seq dataset

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



table(
  scRNA_ec$disease_category,
  scRNA_ec$Disease.subtype
)


nvu_heatmap_modules <- c(
  "Tight_Junctions",
  "BBB_Transport",
  "Angiogenesis",
  "Vascular_Inflam.",
  "ECM_Remodeling"
)


scRNA_patient_summary <- scRNA_ec@meta.data %>%
  dplyr::group_by(
    Patient,
    disease_category
  ) %>%
  dplyr::summarise(
    across(
      all_of(nvu_heatmap_modules),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_cells = n(),
    .groups = "drop"
  )


scRNA_patient_summary

#average patients within disease categories
scRNA_disease_summary <- scRNA_patient_summary %>%
  dplyr::group_by(
    disease_category
  ) %>%
  dplyr::summarise(
    across(
      all_of(nvu_heatmap_modules),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_patients = n(),
    .groups = "drop"
  )


scRNA_disease_summary

#add datatset column
scRNA_disease_summary <- scRNA_disease_summary %>%
  mutate(
    Dataset = "scRNA",
    Disease = disease_category
  )


#Combine datasets
NVU_cross_platform <- bind_rows(
  TMA_disease_summary,
  whole_tissue_disease_summary,
  scRNA_disease_summary
)


NVU_cross_platform
write.csv(NVU_cross_platform, "NVU_cross_platform.csv")


#Heatmap----

nvu_heatmap_modules <- c(
  "Tight_Junctions",
  "BBB_Transport",
  "Angiogenesis",
  "Vascular_Inflam.",
  "ECM_Remodeling"
)


NVU_heatmap_df <- NVU_cross_platform %>%
  dplyr::select(
    Dataset,
    Disease,
    all_of(nvu_heatmap_modules)
  ) %>%
  pivot_longer(
    cols = all_of(nvu_heatmap_modules),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  group_by(
    Dataset,
    Module
  ) %>%
  mutate(
    Zscore = as.numeric(scale(Score))
  ) %>%
  ungroup()


# enforce module order
NVU_heatmap_df$Module <- factor(
  NVU_heatmap_df$Module,
  levels = nvu_heatmap_modules
)



NVU_heatmap_df$Disease <- factor(
  NVU_heatmap_df$Disease,
  levels = c(
    "Non-tumor",
    "LGG",
    "GBM",
    "BrM"
  )
)


NVU_heatmap_df$Dataset <- factor(
  NVU_heatmap_df$Dataset,
  levels = c(
    "TMA",
    "Xenium_whole_tissue",
    "scRNA"
  )
)


id="r1xq2t"
NVU_heatmap_df$Row <- factor(
  paste(
    NVU_heatmap_df$Disease,
    NVU_heatmap_df$Dataset,
    sep = " | "
  ),
  levels = c(
    paste("Non-tumor",
          c("TMA","Xenium_whole_tissue","scRNA"),
          sep=" | "),
    paste("LGG",
          c("TMA","Xenium_whole_tissue","scRNA"),
          sep=" | "),
    paste("GBM",
          c("TMA","Xenium_whole_tissue","scRNA"),
          sep=" | "),
    paste("BrM",
          c("TMA","Xenium_whole_tissue","scRNA"),
          sep=" | ")
  )
)

## Rename and order datasets

NVU_heatmap_df$Dataset <- recode(
  NVU_heatmap_df$Dataset,
  "scRNA" = "scRNA-seq",
  "TMA" = "Xenium TMA",
  "Xenium_whole_tissue" = "Xenium whole tissue"
)

NVU_heatmap_df$Dataset <- factor(
  NVU_heatmap_df$Dataset,
  levels = c(
    "scRNA-seq",
    "Xenium TMA",
    "Xenium whole tissue"
  )
)


## Order diseases

NVU_heatmap_df$Disease <- factor(
  NVU_heatmap_df$Disease,
  levels = rev(c(
    "Non-tumor",
    "LGG",
    "GBM",
    "BrM"
  )
  ))

## Heatmap

ggplot(
  NVU_heatmap_df,
  aes(
    x = Module,
    y = Dataset,
    fill = Zscore
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  facet_grid(
    Disease ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Relative\nmodule score"
  ) +
  theme_classic(
    base_size = 14
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      #face = "bold",
      size = 18
    ),
    axis.text.y = element_text(
      size = 13
    ),
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 90,
      face = "bold",
      size = 16
    ),
    strip.placement = "outside",
    panel.spacing = unit(0.6, "lines"),
    legend.position = "right"
  )

library(svglite)

ggsave(
  filename = "NVU dereg. framework/Module_summary_htmap.svg",
  plot = last_plot(),
  device = "svg",
  width = 8,
  height = 8,
  units = "in"
)



#PCA plot----

pca_modules <- c(
  "Tight_Junctions",
  "BBB_Transport",
  "Angiogenesis",
  "Vascular_Inflam.",
  "ECM_Remodeling"
)

pca_df_scaled <- NVU_cross_platform %>%
  dplyr::select(
    Dataset,
    Disease,
    all_of(pca_modules)
  ) %>%
  group_by(Dataset) %>%
  mutate(
    across(
      all_of(pca_modules),
      ~ as.numeric(scale(.x))
    )
  ) %>%
  ungroup()


pca_res <- prcomp(
  pca_df_scaled[,pca_modules],
  scale. = FALSE
)


pca_scores <- as.data.frame(
  pca_res$x
) %>%
  bind_cols(
    pca_df_scaled %>% 
      dplyr::select(Dataset,Disease)
  )


#variance explained
summary(pca_res) 



#check PC1
pca_res$rotation
#Tight Junctions and BBB Transport drive PC1

#Rename labels
pca_scores$Dataset <- dplyr::recode(
  pca_scores$Dataset,
  "scRNA" = "scRNA-seq",
  "TMA" = "Xenium TMA",
  "Xenium_whole_tissue" = "Xenium whole tissue"
)

pca_scores$Dataset <- factor(
  pca_scores$Dataset,
  levels = c(
    "scRNA-seq",
    "Xenium TMA",
    "Xenium whole tissue"
  )
)

## Disease consensus (average across all datasets)

disease_centroids <- pca_scores %>%
  group_by(Disease) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  ) %>%
  mutate(
    Disease = factor(
      Disease,
      levels = c(
        "Non-tumor",
        "LGG",
        "GBM",
        "BrM"
      )
    )
  )

## PCA plot
ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = Disease,
    shape = Dataset
  )
) +
  
  # Reference lines
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey60"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey60"
  ) +
  
  # Individual datasets
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  
  # Disease consensus across datasets
  geom_point(
    data = disease_centroids,
    aes(
      PC1,
      PC2,
      fill = Disease
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 6,
    color = "black",
    stroke = 1.2
  ) +
  
  # Disease labels
  geom_text_repel(
    data = disease_centroids,
    aes(
      PC1,
      PC2,
      label = Disease
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4,
    box.padding = 0.35,
    point.padding = 0.4
  ) +
  
  scale_color_manual(
    values = c(
      "Non-tumor" = "#4DAF4A",
      "LGG" = "#377EB8",
      "GBM" = "#E69F00",
      "BrM" = "#D73027"
    ),
    guide = "none"
  ) +
  
  scale_fill_manual(
    values = c(
      "Non-tumor" = "#4DAF4A",
      "LGG" = "#377EB8",
      "GBM" = "#E69F00",
      "BrM" = "#D73027"
    )
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "PC1 (84% variance)",
    y = "PC2 (10.7% variance)",
    fill = "Samples",
    shape = "Dataset"
  ) +
  
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical"
    
  )


ggsave(
  filename = "New NVU dereg. framework/PCA_summary_plot.svg",
  plot = last_plot(),
  device = "svg",
  width = 6,
  height = 8,
  units = "in"
)


#Loadings plot----


# Extract PC1 and PC2 loadings
loading_df <- as.data.frame(
  pca_res$rotation[, c("PC1", "PC2")]
)

loading_df$Module <- rownames(loading_df)

# Calculate scaling factor so arrows fit with PCA scores
arrow_scale <- max(abs(pca_scores$PC1), abs(pca_scores$PC2)) /
  max(abs(loading_df$PC1), abs(loading_df$PC2)) * 0.6

loading_df <- loading_df %>%
  mutate(
    PC1_arrow = PC1 * arrow_scale,
    PC2_arrow = PC2 * arrow_scale
  )

loading_df <- loading_df %>%
  mutate(
    Module = recode(
      Module,
      "Tight_Junctions" = "Tight Junctions",
      "BBB_Transport" = "BBB Transport",
      "Angiogenesis" = "Angiogenesis",
      "Vascular_Inflam." = "Vascular Inflammation",
      "ECM_Remodeling" = "ECM-Remodeling"
    )
  )

ggplot() +
  
  # PCA module loading arrows
  geom_segment(
    data = loading_df,
    aes(
      x = 0,
      y = 0,
      xend = PC1_arrow,
      yend = PC2_arrow
    ),
    arrow = arrow(
      length = unit(0.25, "cm")
    ),
    linewidth = 1
  ) +
  
  # Module labels
  geom_text_repel(
    data = loading_df,
    aes(
      x = PC1_arrow,
      y = PC2_arrow,
      label = Module
    ),
    size = 5,
    box.padding = 0.8,
    point.padding = 0.5,
    force = 2,
    min.segment.length = 0
  ) +
  
  # Origin line
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  coord_equal() +
  
  theme_classic(
    base_size = 14
  ) +
  
  labs(
    x = "PC1:\
    BBB stability ↔ vascular remodeling",
    y = "PC2:\ 
    Vascular inflammatory remodeling"
  )



ggsave(
  filename = "NVU dereg. framework/PCA_loadings_summary_plot.svg",
  plot = last_plot(),
  device = "svg",
  width = 5,
  height = 8,
  units = "in"
)


#BBB Stability Score Boxplots----
# BBB Stability Score Distributions ----
# BBB STABILITY SCORE BY DISEASE AND DATASET

library(dplyr)
library(ggplot2)

# Prepare cell-level data

# scRNA-seq
scRNA_cell_df <- scRNA_ec[[]] %>%
  filter(!is.na(BBB_Stability_Score)) %>%
  mutate(
    Dataset = "scRNA-seq",
    Sample = Patient
  ) %>%
  select(
    Dataset,
    Sample,
    disease_category,
    BBB_Stability_Score
  )

# Xenium TMA
TMA_cell_df <- TMA_EC_sub[[]] %>%
  filter(!is.na(BBB_Stability_Score)) %>%
  mutate(
    Dataset = "Xenium TMA",
    Sample = sample_id
  ) %>%
  select(
    Dataset,
    Sample,
    disease_category,
    BBB_Stability_Score
  )

# Xenium whole tissue
whole_tissue_cell_df <- ec_whole_tissue[[]] %>%
  filter(!is.na(BBB_Stability)) %>%
  mutate(
    Dataset = "Xenium whole tissue",
    Sample = sample,
    BBB_Stability_Score = BBB_Stability
  ) %>%
  select(
    Dataset,
    Sample,
    disease_category,
    BBB_Stability_Score
  )

# Combine all cell-level data
BBB_cell_df <- bind_rows(
  scRNA_cell_df,
  TMA_cell_df,
  whole_tissue_cell_df
)

# Calculate mean BBB Stability Score for each biological sample
# scRNA-seq = Patient
# TMA = sample_id
# Whole tissue = sample

BBB_sample_means <- BBB_cell_df %>%
  group_by(Dataset, Sample, disease_category) %>%
  summarise(
    Sample_mean = mean(BBB_Stability_Score),
    .groups = "drop"
  )

# Set disease order
BBB_cell_df$disease_category <- factor(
  BBB_cell_df$disease_category,
  levels = c("Non-tumor", "LGG", "GBM", "BrM")
)

BBB_sample_means$disease_category <- factor(
  BBB_sample_means$disease_category,
  levels = c("Non-tumor", "LGG", "GBM", "BrM")
)

# Set dataset order
BBB_cell_df$Dataset <- factor(
  BBB_cell_df$Dataset,
  levels = c(
    "scRNA-seq",
    "Xenium TMA",
    "Xenium whole tissue"
  )
)

BBB_sample_means$Dataset <- factor(
  BBB_sample_means$Dataset,
  levels = c(
    "scRNA-seq",
    "Xenium TMA",
    "Xenium whole tissue"
  )
)

# Plot distributions with sample-level means
ggplot(
  BBB_cell_df,
  aes(
    x = disease_category,
    y = BBB_Stability_Score
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.5
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_point(
    data = BBB_sample_means,
    aes(
      x = disease_category,
      y = Sample_mean
    ),
    size = 2.8,
    position = position_jitter(width = 0.08),
    inherit.aes = FALSE
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
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )



#MAST/GSEA results----
#From Bejarano et al. scRNA seq

#LGG
LGG_GSEA_df <- read.csv("Bejarano et al. 2025 Immunity/Analysis/Bej_MAST_gsea_BBB_LGG_df_normalized.csv")
#Positive
LGG_GSEA_df %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  dplyr::select(Description, NES, p.adjust)


LGG_high_selected <- c(
  "transport across blood-brain barrier",
  "organic anion transport",
  "xenobiotic transport",
  "lipid transport",
  "cell-cell junction organization",
  "transepithelial transport",
  "protein localization to cell junction"
)


#negative NES
LGG_GSEA_df %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  dplyr::select(Description, NES, p.adjust)


LGG_low_selected <- c(
  "cytoplasmic translation",
  "ribosome biogenesis",
  "regulation of endothelial cell chemotaxis"
)




#GBM
GBM_GSEA_df <- read.csv("Bejarano et al. 2025 Immunity/Analysis/Bej_MAST_gsea_BBB_GBM_df_normalized.csv")

#Positive
GBM_GSEA_df %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  dplyr::select(Description, NES, p.adjust)



GBM_high_selected <- c(
  "transport across blood-brain barrier",
  "organic anion transport",
  "monocarboxylic acid transport",
  "lipid transport",
  "detoxification",
  "response to toxic substance"
)


#negative NES
GBM_GSEA_df %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  dplyr::select(Description, NES, p.adjust)


GBM_low_selected <- character(0)


#Brain Mets
BrM_GSEA_df <- read.csv("Bejarano et al. 2025 Immunity/Analysis/Bej_MAST_gsea_BBB_BrM_df_normalized.csv")
#Positive
BrM_GSEA_df %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  dplyr::select(Description, NES, p.adjust)


BrM_high_selected <- c(
  "transport across blood-brain barrier",
  "organic anion transport",
  "xenobiotic transport",
  "import across plasma membrane",
  "response to toxic substance",
  "regulation of cytokine production involved in immune response",
  "antigen processing and presentation",
  "interferon-mediated signaling pathway",
  "regulation of adaptive immune response"
)


#negative NES
BrM_GSEA_df %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  dplyr::select(Description, NES, p.adjust)




BrM_low_selected <- character(0)


# MAST / GSEA SUMMARY HEATMAP----
# Individual GO terms + biological program

library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(stringr)
library(svglite)


## 1. CURATE INDIVIDUAL GO TERMS

# LGG

LGG_high_selected <- c(
  "transport across blood-brain barrier",
  "organic anion transport",
  "xenobiotic transport",
  "lipid transport",
  "cell-cell junction organization",
  "transepithelial transport",
  "protein localization to cell junction"
)

LGG_low_selected <- c(
  "cytoplasmic translation",
  "ribosome biogenesis",
  "regulation of endothelial cell chemotaxis"
)


# GBM

GBM_high_selected <- c(
  "transport across blood-brain barrier",
  "organic anion transport",
  "monocarboxylic acid transport",
  "lipid transport",
  "detoxification",
  "response to toxic substance"
)

GBM_low_selected <- character(0)


# Brain metastasis

BrM_high_selected <- c(
  "transport across blood-brain barrier",
  "organic anion transport",
  "xenobiotic transport",
  "import across plasma membrane",
  "response to toxic substance",
  "regulation of cytokine production involved in immune response",
  "antigen processing and presentation",
  "interferon-mediated signaling pathway",
  "regulation of adaptive immune response"
)

BrM_low_selected <- character(0)


## 2. COMBINE ALL SELECTED TERMS

all_selected <- unique(c(
  LGG_high_selected,
  LGG_low_selected,
  GBM_high_selected,
  GBM_low_selected,
  BrM_high_selected,
  BrM_low_selected
))


## 3. EXTRACT NES VALUES

LGG_heat <- LGG_GSEA_df %>%
  filter(
    Description %in% all_selected
  ) %>%
  select(
    Description,
    NES
  ) %>%
  mutate(
    Tumor = "LGG"
  )


GBM_heat <- GBM_GSEA_df %>%
  filter(
    Description %in% all_selected
  ) %>%
  select(
    Description,
    NES
  ) %>%
  mutate(
    Tumor = "GBM"
  )


BrM_heat <- BrM_GSEA_df %>%
  filter(
    Description %in% all_selected
  ) %>%
  select(
    Description,
    NES
  ) %>%
  mutate(
    Tumor = "BrM"
  )


heat_df <- bind_rows(
  LGG_heat,
  GBM_heat,
  BrM_heat
)


## 4. DEFINE BBB DIRECTION

high_terms <- unique(c(
  LGG_high_selected,
  GBM_high_selected,
  BrM_high_selected
))

low_terms <- unique(c(
  LGG_low_selected,
  GBM_low_selected,
  BrM_low_selected
))


heat_df <- heat_df %>%
  mutate(
    BBB_State = case_when(
      Description %in% high_terms ~ "BBB preserved",
      Description %in% low_terms ~ "BBB disrupted",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(BBB_State)
  )


## 5. DEFINE BIOLOGICAL PROGRAMS

# Only programs represented by the CURRENT selected pathways

program_order <- c(
  "BBB / tight junction",
  "BBB transport",
  "Migration / chemotaxis",
  "Immune surveillance / interferon",
  "Chemokine / inflammatory activation",
  "Translation / ribosome",
  "Detoxification / stress"
)


heat_df <- heat_df %>%
  mutate(
    Program = case_when(
      
      # BBB transport
      Description %in% c(
        "transport across blood-brain barrier",
        "organic anion transport",
        "xenobiotic transport",
        "lipid transport",
        "monocarboxylic acid transport",
        "transepithelial transport",
        "import across plasma membrane"
      ) ~ "BBB transport",
      
      
      # BBB / tight junction
      Description %in% c(
        "cell-cell junction organization",
        "protein localization to cell junction"
      ) ~ "BBB / tight junction",
      
      
      # Migration / chemotaxis
      Description %in% c(
        "regulation of endothelial cell chemotaxis"
      ) ~ "Migration / chemotaxis",
      
      
      # Immune surveillance / interferon
      Description %in% c(
        "antigen processing and presentation",
        "interferon-mediated signaling pathway",
        "regulation of adaptive immune response"
      ) ~ "Immune surveillance / interferon",
      
      
      # Chemokine / inflammatory activation
      Description %in% c(
        "regulation of cytokine production involved in immune response"
      ) ~ "Chemokine / inflammatory activation",
      
      
      # Translation / ribosome
      Description %in% c(
        "cytoplasmic translation",
        "ribosome biogenesis"
      ) ~ "Translation / ribosome",
      
      
      # Detoxification / stress
      Description %in% c(
        "detoxification",
        "response to toxic substance"
      ) ~ "Detoxification / stress",
      
      
      TRUE ~ NA_character_
    ),
    
    Program = factor(
      Program,
      levels = program_order
    )
  ) %>%
  filter(
    !is.na(Program)
  )


## 6. CREATE HEATMAP MATRIX

heat_matrix <- heat_df %>%
  select(
    Description,
    Tumor,
    BBB_State,
    NES
  ) %>%
  pivot_wider(
    names_from = c(
      Tumor,
      BBB_State
    ),
    values_from = NES
  ) %>%
  column_to_rownames(
    "Description"
  ) %>%
  as.matrix()


## 7. ORDER COLUMNS

desired_columns <- c(
  "LGG_BBB preserved",
  "GBM_BBB preserved",
  "BrM_BBB preserved",
  "LGG_BBB disrupted",
  "GBM_BBB disrupted",
  "BrM_BBB disrupted"
)


# Add missing columns as NA
# GBM and BrM currently have no selected BBB-disrupted pathways

missing_columns <- setdiff(
  desired_columns,
  colnames(heat_matrix)
)

if (length(missing_columns) > 0) {
  
  heat_matrix <- cbind(
    heat_matrix,
    matrix(
      NA,
      nrow = nrow(heat_matrix),
      ncol = length(missing_columns),
      dimnames = list(
        rownames(heat_matrix),
        missing_columns
      )
    )
  )
}


heat_matrix <- heat_matrix[
  ,
  desired_columns,
  drop = FALSE
]


## 8. CREATE ROW INFORMATION

row_info <- heat_df %>%
  select(
    Description,
    BBB_State,
    Program
  ) %>%
  distinct()


## 9. ORDER ROWS

row_order_df <- row_info %>%
  mutate(
    BBB_order = case_when(
      BBB_State == "BBB preserved" ~ 1,
      BBB_State == "BBB disrupted" ~ 2
    )
  ) %>%
  arrange(
    BBB_order,
    Program,
    Description
  )


ordered_rows <- row_order_df$Description[
  row_order_df$Description %in%
    rownames(heat_matrix)
]


heat_matrix_plot <- heat_matrix[
  ordered_rows,
  ,
  drop = FALSE
]


## 10. MATCH ROW ANNOTATION

row_annotation <- row_order_df %>%
  filter(
    Description %in% rownames(heat_matrix_plot)
  ) %>%
  select(
    Description,
    Program
  )


rownames(row_annotation) <- row_annotation$Description


row_annotation <- row_annotation[
  rownames(heat_matrix_plot),
  ,
  drop = FALSE
]


## 11. WRAP PATHWAY NAMES

wrapped_names <- stringr::str_wrap(
  rownames(heat_matrix_plot),
  width = 42
)


rownames(heat_matrix_plot) <- wrapped_names

rownames(row_annotation) <- wrapped_names


## 12. NES COLOR SCALE

max_abs <- max(
  abs(heat_matrix_plot),
  na.rm = TRUE
)


col_fun <- circlize::colorRamp2(
  c(
    -max_abs,
    0,
    max_abs
  ),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)


## 13. BIOLOGICAL PROGRAM COLORS

program_colors <- c(
  "BBB / tight junction" = "#7B3294",
  "BBB transport" = "#008837",
  "Migration / chemotaxis" = "#E08214",
  "Immune surveillance / interferon" = "#C51B7D",
  "Chemokine / inflammatory activation" = "#DE77AE",
  "Translation / ribosome" = "#762A83",
  "Detoxification / stress" = "#8C510A"
)


program_colors <- program_colors[
  program_order
]


## 14. LEFT ROW ANNOTATION

row_ha <- rowAnnotation(
  
  Program = row_annotation$Program,
  
  col = list(
    Program = program_colors
  ),
  
  annotation_name_side = "top",
  
  annotation_legend_param = list(
    Program = list(
      
      title = "Biological programs",
      
      title_gp = gpar(
        fontsize = 13,
        fontface = "bold"
      ),
      
      labels_gp = gpar(
        fontsize = 11
      ),
      
      nrow = 3,
      ncol = 3,
      
      grid_height = unit(
        5,
        "mm"
      ),
      
      grid_width = unit(
        5,
        "mm"
      )
    )
  )
)


## 15. DETERMINE BBB ROW SPLIT

ordered_direction <- row_order_df$BBB_State[
  match(
    gsub(
      "\n",
      " ",
      rownames(heat_matrix_plot)
    ),
    row_order_df$Description
  )
]


high_n <- sum(
  ordered_direction == "BBB preserved",
  na.rm = TRUE
)


low_n <- sum(
  ordered_direction == "BBB disrupted",
  na.rm = TRUE
)


## 16. COLUMN SPLIT

column_split <- factor(
  c(
    rep(
      "BBB preserved",
      3
    ),
    rep(
      "BBB disrupted",
      3
    )
  ),
  levels = c(
    "BBB preserved",
    "BBB disrupted"
  )
)


## 17. CLEAN COLUMN NAMES

colnames(heat_matrix_plot) <- c(
  "LGG",
  "GBM",
  "BrM",
  "LGG",
  "GBM",
  "BrM"
)


## 18. HEATMAP

ht <- ComplexHeatmap::Heatmap(
  
  heat_matrix_plot,
  
  name = "NES",
  
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  left_annotation = row_ha,
  
  row_split = factor(
    c(
      rep(
        "BBB preserved",
        high_n
      ),
      rep(
        "BBB disrupted",
        low_n
      )
    ),
    levels = c(
      "BBB preserved",
      "BBB disrupted"
    )
  ),
  
  row_title = NULL,
  
  column_split = column_split,
  
  column_title = c(
    "BBB preserved",
    "BBB disrupted"
  ),
  
  column_title_side = "top",
  
  column_title_gp = gpar(
    fontsize = 12,
    fontface = "bold"
  ),
  
  row_gap = unit(
    2.5,
    "mm"
  ),
  
  column_gap = unit(
    3,
    "mm"
  ),
  
  border = TRUE,
  
  rect_gp = gpar(
    col = "grey100",
    lwd = 0.5
  ),
  
  row_names_gp = gpar(
    fontsize = 12
  ),
  
  row_names_max_width = unit(
    9,
    "cm"
  ),
  
  column_names_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  column_names_rot = 0,
  
  heatmap_legend_param = list(
    title = "Normalized enrichment score",
    direction = "horizontal"
  )
)


## 19. DRAW AND SAVE AS SVG

MAST_GSEA_svg_file <- "New NVU dereg. framework/MAST_GSEA_summary_plot.svg"


svglite(
  file = MAST_GSEA_svg_file,
  width = 9,
  height = 16
)


draw(
  ht,
  
  heatmap_legend_side = "bottom",
  
  annotation_legend_side = "bottom",
  
  merge_legends = TRUE,
  
  padding = unit(
    c(
      5,
      5,
      5,
      5
    ),
    "mm"
  )
)


dev.off()



#Subset heatmap Pos NES----
# Individual GO terms + biological program
# Main figure: pathways positively associated with higher BBB Stability Score

library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(stringr)
library(svglite)


## 1. CURATE INDIVIDUAL GO TERMS

# LGG

LGG_high_selected <- c(
  "transport across blood-brain barrier",
  "organic anion transport",
  "xenobiotic transport",
  "lipid transport",
  "cell-cell junction organization",
  "transepithelial transport",
  "protein localization to cell junction"
)


# GBM

GBM_high_selected <- c(
  "transport across blood-brain barrier",
  "organic anion transport",
  "monocarboxylic acid transport",
  "lipid transport",
  "detoxification",
  "response to toxic substance"
)


# Brain metastasis

BrM_high_selected <- c(
  "transport across blood-brain barrier",
  "organic anion transport",
  "xenobiotic transport",
  "import across plasma membrane",
  "response to toxic substance",
  "regulation of cytokine production involved in immune response",
  "antigen processing and presentation",
  "interferon-mediated signaling pathway",
  "regulation of adaptive immune response"
)


## 2. COMBINE ALL SELECTED TERMS

all_selected <- unique(c(
  LGG_high_selected,
  GBM_high_selected,
  BrM_high_selected
))


## 3. EXTRACT NES VALUES

LGG_heat <- LGG_GSEA_df %>%
  filter(
    Description %in% all_selected
  ) %>%
  select(
    Description,
    NES
  ) %>%
  mutate(
    Tumor = "LGG"
  )


GBM_heat <- GBM_GSEA_df %>%
  filter(
    Description %in% all_selected
  ) %>%
  select(
    Description,
    NES
  ) %>%
  mutate(
    Tumor = "GBM"
  )


BrM_heat <- BrM_GSEA_df %>%
  filter(
    Description %in% all_selected
  ) %>%
  select(
    Description,
    NES
  ) %>%
  mutate(
    Tumor = "BrM"
  )


heat_df <- bind_rows(
  LGG_heat,
  GBM_heat,
  BrM_heat
)


## 4. DEFINE BBB DIRECTION

# All pathways shown in the main figure are positively
# associated with higher BBB Stability Score

heat_df <- heat_df %>%
  mutate(
    BBB_State = "High BBB Stability"
  )


## 5. DEFINE BIOLOGICAL PROGRAMS

program_order <- c(
  "BBB / tight junction",
  "BBB transport",
  "Immune surveillance / interferon",
  "Inflammatory signaling",
  "Detoxification / stress"
)


heat_df <- heat_df %>%
  mutate(
    Program = case_when(
      
      # BBB transport
      Description %in% c(
        "transport across blood-brain barrier",
        "organic anion transport",
        "xenobiotic transport",
        "lipid transport",
        "monocarboxylic acid transport",
        "transepithelial transport",
        "import across plasma membrane"
      ) ~ "BBB transport",
      
      
      # BBB / tight junction
      Description %in% c(
        "cell-cell junction organization",
        "protein localization to cell junction"
      ) ~ "BBB / tight junction",
      
      
      # Immune surveillance / interferon
      Description %in% c(
        "antigen processing and presentation",
        "interferon-mediated signaling pathway",
        "regulation of adaptive immune response"
      ) ~ "Immune surveillance / interferon",
      
      
      # Inflammatory signaling
      Description %in% c(
        "regulation of cytokine production involved in immune response"
      ) ~ "Inflammatory signaling",
      
      
      # Detoxification / stress
      Description %in% c(
        "detoxification",
        "response to toxic substance"
      ) ~ "Detoxification / stress",
      
      
      TRUE ~ NA_character_
    ),
    
    Program = factor(
      Program,
      levels = program_order
    )
  ) %>%
  filter(
    !is.na(Program)
  )


## 6. CREATE HEATMAP MATRIX

heat_matrix <- heat_df %>%
  select(
    Description,
    Tumor,
    NES
  ) %>%
  pivot_wider(
    names_from = Tumor,
    values_from = NES
  ) %>%
  column_to_rownames(
    "Description"
  ) %>%
  as.matrix()


## 7. ORDER COLUMNS

desired_columns <- c(
  "LGG",
  "GBM",
  "BrM"
)


# Add missing columns as NA if necessary

missing_columns <- setdiff(
  desired_columns,
  colnames(heat_matrix)
)

if (length(missing_columns) > 0) {
  
  heat_matrix <- cbind(
    heat_matrix,
    matrix(
      NA,
      nrow = nrow(heat_matrix),
      ncol = length(missing_columns),
      dimnames = list(
        rownames(heat_matrix),
        missing_columns
      )
    )
  )
}


heat_matrix <- heat_matrix[
  ,
  desired_columns,
  drop = FALSE
]


## 8. CREATE ROW INFORMATION

row_info <- heat_df %>%
  select(
    Description,
    BBB_State,
    Program
  ) %>%
  distinct()


## 9. ORDER ROWS

row_order_df <- row_info %>%
  mutate(
    Program_order = match(
      Program,
      program_order
    )
  ) %>%
  arrange(
    Program_order,
    Description
  )


ordered_rows <- row_order_df$Description[
  row_order_df$Description %in%
    rownames(heat_matrix)
]


heat_matrix_plot <- heat_matrix[
  ordered_rows,
  ,
  drop = FALSE
]


## 10. MATCH ROW ANNOTATION

row_annotation <- row_order_df %>%
  filter(
    Description %in% rownames(heat_matrix_plot)
  ) %>%
  select(
    Description,
    Program
  )


rownames(row_annotation) <- row_annotation$Description


row_annotation <- row_annotation[
  rownames(heat_matrix_plot),
  ,
  drop = FALSE
]


## 11. WRAP PATHWAY NAMES

wrapped_names <- stringr::str_wrap(
  rownames(heat_matrix_plot),
  width = 42
)


rownames(heat_matrix_plot) <- wrapped_names

rownames(row_annotation) <- wrapped_names


## 12. NES COLOR SCALE

# All displayed NES values are positive.
# Use a sequential scale rather than a diverging scale.

max_nes <- max(
  heat_matrix_plot,
  na.rm = TRUE
)


col_fun <- circlize::colorRamp2(
  c(
    0,
    max_nes
  ),
  c(
    "#F7FBFF",
    "#2166AC"
  )
)


## 13. BIOLOGICAL PROGRAM COLORS

program_colors <- c(
  "BBB / tight junction" = "#7B3294",
  "BBB transport" = "#008837",
  "Immune surveillance / interferon" = "#C51B7D",
  "Inflammatory signaling" = "#DE77AE",
  "Detoxification / stress" = "#8C510A"
)


program_colors <- program_colors[
  program_order
]


## 14. LEFT ROW ANNOTATION

row_ha <- rowAnnotation(
  
  Program = row_annotation$Program,
  
  col = list(
    Program = program_colors
  ),
  
  annotation_name_side = "bottom",
  
  annotation_legend_param = list(
    Program = list(
      
      title = "Biological programs",
      
      title_gp = gpar(
        fontsize = 13,
        fontface = "bold"
      ),
      
      labels_gp = gpar(
        fontsize = 11
      ),
      
      nrow = 2,
      ncol = 3,
      
      grid_height = unit(
        5,
        "mm"
      ),
      
      grid_width = unit(
        5,
        "mm"
      )
    )
  )
)


## 15. COLUMN SPLIT

# All columns represent pathways positively associated
# with higher BBB Stability Score

column_split <- factor(
  rep(
    "High BBB Stability",
    3
  ),
  levels = "High BBB Stability"
)


## 16. CLEAN COLUMN NAMES

colnames(heat_matrix_plot) <- c(
  "LGG",
  "GBM",
  "BrM"
)


## 17. HEATMAP

ht <- ComplexHeatmap::Heatmap(
  
  heat_matrix_plot,
  
  name = "NES",
  
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  left_annotation = row_ha,
  
  row_split = row_annotation$Program,
  
  row_title = NULL,
  
  column_split = column_split,
  
  column_title = "High BBB Stability Associated
  Pathways: scRNA seq",
  
  column_title_side = "top",
  
  # Title of the split column group
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),
  
  row_gap = unit(
    2.5,
    "mm"
  ),
  
  column_gap = unit(
    3,
    "mm"
  ),
  
  border = TRUE,
  
  rect_gp = gpar(
    col = "grey100",
    lwd = 0.5
  ),
  
  # Larger pathway names
  row_names_gp = gpar(
    fontsize = 15
  ),
  row_names_max_width = unit(
    9,
    "cm"
  ),
  
  column_names_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  column_names_rot = 0,
  
  heatmap_legend_param = list(
    title = "Normalized enrichment score",
    direction = "horizontal"
  )
)


## 18. DRAW AND SAVE AS SVG

MAST_GSEA_svg_file <- "Vascular subtype analysis/Figures/Pos_MAST_GSEA_summary_plot.svg"


svglite(
  file = MAST_GSEA_svg_file,
  width = 7,
  height = 11
)


draw(
  ht,
  
  heatmap_legend_side = "bottom",
  
  annotation_legend_side = "bottom",
  
  merge_legends = TRUE,
  
  padding = unit(
    c(
      5,
      5,
      5,
      5
    ),
    "mm"
  )
)


dev.off()





# Leading Edge----
# Compare leading-edge genes across LGG, GBM, and BrM
# Separate High BBB and Low BBB analyses
# Based ONLY on the FINAL SELECTED GSEA PATHWAYS


## 20. DEFINE FINAL PATHWAY -> PROGRAM MAPPING

# This mapping contains ONLY pathways included in the
# finalized GSEA pathway selections used in the heatmap.

program_map <- tribble(
  
  ~Description, ~Program,
  
  # BBB transport
  "transport across blood-brain barrier",
  "BBB transport",
  
  "organic anion transport",
  "BBB transport",
  
  "xenobiotic transport",
  "BBB transport",
  
  "lipid transport",
  "BBB transport",
  
  "monocarboxylic acid transport",
  "BBB transport",
  
  "transepithelial transport",
  "BBB transport",
  
  "import across plasma membrane",
  "BBB transport",
  
  
  # BBB / tight junction
  "cell-cell junction organization",
  "BBB / tight junction",
  
  "protein localization to cell junction",
  "BBB / tight junction",
  
  
  # Migration / chemotaxis
  "regulation of endothelial cell chemotaxis",
  "Migration / chemotaxis",
  
  
  # Immune surveillance / interferon
  "antigen processing and presentation",
  "Immune surveillance / interferon",
  
  "interferon-mediated signaling pathway",
  "Immune surveillance / interferon",
  
  "regulation of adaptive immune response",
  "Immune surveillance / interferon",
  
  
  # Chemokine / inflammatory activation
  "regulation of cytokine production involved in immune response",
  "Chemokine / inflammatory activation",
  
  
  # Translation / ribosome
  "cytoplasmic translation",
  "Translation / ribosome",
  
  "ribosome biogenesis",
  "Translation / ribosome",
  
  
  # Detoxification / stress
  "detoxification",
  "Detoxification / stress",
  
  "response to toxic substance",
  "Detoxification / stress"
)


## 21. DEFINE FINAL SELECTED PATHWAYS

# Keep this synchronized with the MAST/GSEA heatmap.

final_high_pathways <- list(
  
  LGG = LGG_high_selected,
  
  GBM = GBM_high_selected,
  
  BrM = BrM_high_selected
)


final_low_pathways <- list(
  
  LGG = LGG_low_selected,
  
  GBM = GBM_low_selected,
  
  BrM = BrM_low_selected
)


## 22. CHECK PATHWAY MAPPING

all_final_pathways <- unique(
  c(
    unlist(final_high_pathways),
    unlist(final_low_pathways)
  )
)


unmapped_pathways <- setdiff(
  all_final_pathways,
  program_map$Description
)


if (
  length(unmapped_pathways) > 0
) {
  
  warning(
    "The following final GSEA pathways are not assigned to a biological program: ",
    paste(
      unmapped_pathways,
      collapse = ", "
    )
  )
  
}


extra_program_paths <- setdiff(
  program_map$Description,
  all_final_pathways
)


if (
  length(extra_program_paths) > 0
) {
  
  message(
    "Program map contains pathways not currently selected: ",
    paste(
      extra_program_paths,
      collapse = ", "
    )
  )
  
}


## 23. EXTRACT LEADING-EDGE GENES

extract_LE_program <- function(
    gsea_df,
    selected_pathways,
    direction = c(
      "positive",
      "negative"
    )
){
  
  direction <- match.arg(
    direction
  )
  
  
  if (
    length(selected_pathways) == 0 ||
    nrow(gsea_df) == 0
  ) {
    
    return(
      tibble()
    )
  }
  
  
  gsea_df %>%
    
    # ONLY finalized selected pathways
    filter(
      Description %in% selected_pathways
    ) %>%
    
    # Positive = higher BBB Stability
    # Negative = lower BBB Stability
    filter(
      if (
        direction == "positive"
      ) {
        NES > 0
      } else {
        NES < 0
      }
    ) %>%
    
    select(
      Description,
      NES,
      core_enrichment
    ) %>%
    
    filter(
      !is.na(core_enrichment),
      core_enrichment != ""
    ) %>%
    
    separate_rows(
      core_enrichment,
      sep = "/"
    ) %>%
    
    rename(
      Gene = core_enrichment
    ) %>%
    
    left_join(
      program_map,
      by = "Description"
    ) %>%
    
    filter(
      !is.na(Program),
      !grepl(
        "^MT[-.]",
        Gene
      )
    ) %>%
    
    distinct()
}


## 24. HIGH BBB / BBB-PRESERVED LEADING EDGE

high_LE_program <- bind_rows(
  
  extract_LE_program(
    LGG_GSEA_df,
    final_high_pathways$LGG,
    direction = "positive"
  ) %>%
    mutate(
      Tumor = "LGG"
    ),
  
  
  extract_LE_program(
    GBM_GSEA_df,
    final_high_pathways$GBM,
    direction = "positive"
  ) %>%
    mutate(
      Tumor = "GBM"
    ),
  
  
  extract_LE_program(
    BrM_GSEA_df,
    final_high_pathways$BrM,
    direction = "positive"
  ) %>%
    mutate(
      Tumor = "BrM"
    )
) %>%
  mutate(
    BBB_state = "BBB preserved"
  )


## 25. LOW BBB / BBB-DISRUPTED LEADING EDGE

low_LE_program <- bind_rows(
  
  extract_LE_program(
    LGG_GSEA_df,
    final_low_pathways$LGG,
    direction = "negative"
  ) %>%
    mutate(
      Tumor = "LGG"
    ),
  
  
  extract_LE_program(
    GBM_GSEA_df,
    final_low_pathways$GBM,
    direction = "negative"
  ) %>%
    mutate(
      Tumor = "GBM"
    ),
  
  
  extract_LE_program(
    BrM_GSEA_df,
    final_low_pathways$BrM,
    direction = "negative"
  ) %>%
    mutate(
      Tumor = "BrM"
    )
) %>%
  mutate(
    BBB_state = "BBB disrupted"
  )


## 26. CHECK FINAL LE PATHWAYS

cat(
  "\nHIGH BBB leading-edge pathways:\n"
)

print(
  high_LE_program %>%
    distinct(
      Tumor,
      Program,
      Description
    ) %>%
    arrange(
      Tumor,
      Program,
      Description
    )
)


cat(
  "\nLOW BBB leading-edge pathways:\n"
)

print(
  low_LE_program %>%
    distinct(
      Tumor,
      Program,
      Description
    ) %>%
    arrange(
      Tumor,
      Program,
      Description
    )
)


## 27. RANK LEADING-EDGE GENES

high_program_ranked <- high_LE_program %>%
  
  group_by(
    Program,
    Gene
  ) %>%
  
  summarise(
    
    Tumor_Count = n_distinct(
      Tumor
    ),
    
    Pathway_Count = n_distinct(
      Description
    ),
    
    Tumors = paste(
      sort(
        unique(
          Tumor
        )
      ),
      collapse = " + "
    ),
    
    .groups = "drop"
  )


low_program_ranked <- low_LE_program %>%
  
  group_by(
    Program,
    Gene
  ) %>%
  
  summarise(
    
    Tumor_Count = n_distinct(
      Tumor
    ),
    
    Pathway_Count = n_distinct(
      Description
    ),
    
    Tumors = paste(
      sort(
        unique(
          Tumor
        )
      ),
      collapse = " + "
    ),
    
    .groups = "drop"
  )


## 28. CHECKS

cat(
  "\nHIGH BBB LE genes by program and tumor:\n"
)

print(
  table(
    high_LE_program$Program,
    high_LE_program$Tumor
  )
)


cat(
  "\nLOW BBB LE genes by program and tumor:\n"
)

print(
  table(
    low_LE_program$Program,
    low_LE_program$Tumor
  )
)


cat(
  "\nHIGH BBB pathway contributions:\n"
)

print(
  high_LE_program %>%
    dplyr::count(
      Program,
      Description,
      Tumor
    )
)


cat(
  "\nLOW BBB pathway contributions:\n"
)

print(
  low_LE_program %>%
    count(
      Program,
      Description,
      Tumor
    )
)


## 29. PREPARE RANKED TABLES FOR PLOTTING

high_program_ranked_plot <- high_LE_program %>%
  
  filter(
    !grepl(
      "^MT[-.]",
      Gene
    )
  ) %>%
  
  group_by(
    Program,
    Gene
  ) %>%
  
  summarise(
    
    Tumor_Count = n_distinct(
      Tumor
    ),
    
    Pathway_Count = n_distinct(
      Description
    ),
    
    Tumors = paste(
      sort(
        unique(
          Tumor
        )
      ),
      collapse = " + "
    ),
    
    .groups = "drop"
  )


low_program_ranked_plot <- low_LE_program %>%
  
  filter(
    !grepl(
      "^MT[-.]",
      Gene
    )
  ) %>%
  
  group_by(
    Program,
    Gene
  ) %>%
  
  summarise(
    
    Tumor_Count = n_distinct(
      Tumor
    ),
    
    Pathway_Count = n_distinct(
      Description
    ),
    
    Tumors = paste(
      sort(
        unique(
          Tumor
        )
      ),
      collapse = " + "
    ),
    
    .groups = "drop"
  )


## 30. ORDER PROGRAMS

program_order <- c(
  "BBB / tight junction",
  "BBB transport",
  "Migration / chemotaxis",
  "Immune surveillance / interferon",
  "Chemokine / inflammatory activation",
  "Translation / ribosome",
  "Detoxification / stress"
)


high_program_ranked_plot$Program <- factor(
  high_program_ranked_plot$Program,
  levels = program_order
)


low_program_ranked_plot$Program <- factor(
  low_program_ranked_plot$Program,
  levels = program_order
)


## 31. SELECT TOP 15 GENES PER PROGRAM

high_program_ranked_plot <- high_program_ranked_plot %>%
  
  group_by(
    Program
  ) %>%
  
  arrange(
    desc(Tumor_Count),
    desc(Pathway_Count),
    Gene
  ) %>%
  
  slice_head(
    n = 15
  ) %>%
  
  ungroup()


low_program_ranked_plot <- low_program_ranked_plot %>%
  
  group_by(
    Program
  ) %>%
  
  arrange(
    desc(Tumor_Count),
    desc(Pathway_Count),
    Gene
  ) %>%
  
  slice_head(
    n = 15
  ) %>%
  
  ungroup()


## 32. PLOTTING FUNCTION

plot_LE_program <- function(
    rank_df,
    title
){
  
  if (
    nrow(rank_df) == 0
  ) {
    
    return(
      ggplot() +
        theme_void() +
        labs(
          title = title
        )
    )
  }
  
  
  plot_df <- rank_df %>%
    
    separate_rows(
      Tumors,
      sep = " \\+ "
    ) %>%
    
    mutate(
      
      Tumors = factor(
        Tumors,
        levels = c(
          "LGG",
          "GBM",
          "BrM"
        )
      ),
      
      Gene = str_remove(
        Gene,
        "^MT[-.]"
      )
    ) %>%
    
    group_by(
      Program
    ) %>%
    
    mutate(
      
      Gene = factor(
        Gene,
        levels = unique(
          Gene[
            order(
              Tumor_Count,
              Pathway_Count
            )
          ]
        )
      )
    ) %>%
    
    ungroup()
  
  
  ggplot(
    plot_df,
    aes(
      x = Tumors,
      y = Gene,
      size = Pathway_Count
    )
  ) +
    
    geom_point(
      alpha = 0.8
    ) +
    
    facet_wrap(
      ~Program,
      scales = "free_y",
      ncol = 1
    ) +
    
    scale_size_continuous(
      name = "Number of\ncontributing pathways",
      range = c(
        2,
        8
      )
    ) +
    
    labs(
      x = NULL,
      y = NULL,
      title = title
    ) +
    
    theme_bw(
      base_size = 12
    ) +
    
    theme(
      
      strip.text.x = element_text(
        face = "bold",
        size = 14
      ),
      
      strip.background = element_rect(
        linewidth = 0.5
      ),
      
      axis.text.y = element_text(
        size = 16
      ),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 20
      ),
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 16
      ),
      
      legend.position = "bottom"
    )
}


## 33. BBB PRESERVED LEADING EDGE

high_LE_dotplot <- plot_LE_program(
  high_program_ranked_plot,
  "Leading-edge genes\n defining BBB preserved programs"
)


## 34. BBB DISRUPTED LEADING EDGE

low_LE_dotplot <- plot_LE_program(
  low_program_ranked_plot,
  "Leading-edge genes\n defining BBB disrupted programs"
)


## 35. DISPLAY

high_LE_dotplot

low_LE_dotplot


## 36. EXPORT BBB PRESERVED LEADING-EDGE DOTPLOT

ggsave(
  filename = "New NVU dereg. framework/LE_BBB_preserved_programs.svg",
  plot = high_LE_dotplot,
  device = svglite,
  width = 5,
  height = 22,
  units = "in"
)


## 37. EXPORT BBB DISRUPTED LEADING-EDGE DOTPLOT

ggsave(
  filename = "New NVU dereg. framework/LE_BBB_disrupted_programs.svg",
  plot = low_LE_dotplot,
  device = svglite,
  width = 5,
  height = 15,
  units = "in"
)





# Xenium validation of FINAL GSEA leading-edge genes ----

library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)


# 1. DEFINE FINAL LEADING-EDGE GENE SETS

# These genes come directly from the finalized GSEA pathway analysis.

high_LE_genes <- high_LE_program %>%
  
  distinct(
    Gene
  ) %>%
  
  pull(
    Gene
  ) %>%
  
  unique()


low_LE_genes <- low_LE_program %>%
  
  distinct(
    Gene
  ) %>%
  
  pull(
    Gene
  ) %>%
  
  unique()


# Remove mitochondrial genes

high_LE_genes <- high_LE_genes[
  !grepl(
    "^MT[-.]",
    high_LE_genes
  )
]


low_LE_genes <- low_LE_genes[
  !grepl(
    "^MT[-.]",
    low_LE_genes
  )
]


# 2. PREPARE XENIUM EC OBJECT

Idents(
  ec_whole_tissue
) <- ec_whole_tissue$long_label


# Remove cells without an NVU annotation

ec_whole_tissue_noNA <- subset(
  
  ec_whole_tissue,
  
  cells = colnames(
    ec_whole_tissue
  )[
    
    !is.na(
      ec_whole_tissue$long_label
    )
    
  ]
)


# Keep High and Low BBB cells only

ec_whole_tissue_noNA <- subset(
  
  ec_whole_tissue_noNA,
  
  BBB_EC_status_within_sample %in% c(
    "High",
    "Low"
  )
)


# 3. KEEP ONLY GENES PRESENT ON THE XENIUM PANEL

xenium_genes <- rownames(
  ec_whole_tissue_noNA[["Xenium"]]
)


high_LE_genes_xenium <- high_LE_genes[
  high_LE_genes %in%
    xenium_genes
]


low_LE_genes_xenium <- low_LE_genes[
  low_LE_genes %in%
    xenium_genes
]


cat(
  "\nBBB-preserved LE genes:",
  length(high_LE_genes),
  "\n"
)


cat(
  "BBB-preserved LE genes on Xenium panel:",
  length(high_LE_genes_xenium),
  "\n"
)


cat(
  "\nBBB-disrupted LE genes:",
  length(low_LE_genes),
  "\n"
)


cat(
  "BBB-disrupted LE genes on Xenium panel:",
  length(low_LE_genes_xenium),
  "\n"
)


cat(
  "\nBBB-preserved Xenium genes:\n"
)

print(
  high_LE_genes_xenium
)


cat(
  "\nBBB-disrupted Xenium genes:\n"
)

print(
  low_LE_genes_xenium
)

# PART 1 — XENIUM DOTPLOT BY BBB STABILITY GROUP

# 4. COMBINE FINAL LE GENES

shared_LE_genes <- intersect(
  high_LE_genes_xenium,
  low_LE_genes_xenium
)


high_only_LE_genes <- setdiff(
  high_LE_genes_xenium,
  shared_LE_genes
)


low_only_LE_genes <- setdiff(
  low_LE_genes_xenium,
  shared_LE_genes
)


genes_to_plot <- c(
  high_only_LE_genes,
  low_only_LE_genes
)

# 5. EXTRACT DOTPLOT DATA

if (
  length(genes_to_plot) > 0
) {
  
  LE_BBB_dot_data <- DotPlot(
    
    ec_whole_tissue_noNA,
    
    features = genes_to_plot,
    
    group.by = "BBB_EC_status_within_sample",
    
    assay = "Xenium"
    
  )$data
  
  
  # Add biological LE category
  
  LE_BBB_dot_data <- LE_BBB_dot_data %>%
    
    mutate(
      
      Program = case_when(
        
        features.plot %in% high_only_LE_genes ~
          "BBB preserved LE",
        
        features.plot %in% low_only_LE_genes ~
          "BBB disrupted LE",
        
        TRUE ~
          NA_character_
        
      )
      
    )
  
  
  # Order BBB groups
  
  LE_BBB_dot_data$id <- factor(
    
    LE_BBB_dot_data$id,
    
    levels = c(
      "Low",
      "High"
    )
    
  )
  
  
  # Order programs
  
  LE_BBB_dot_data$Program <- factor(
    
    LE_BBB_dot_data$Program,
    
    levels = c(
      "BBB preserved LE",
      "BBB disrupted LE"
    )
    
  )
  
  
  # Gene ordering
  
  LE_BBB_dot_data$Gene <- factor(
    
    LE_BBB_dot_data$features.plot,
    
    levels = genes_to_plot
    
  )
  
  
  # 6. PLOT — FACETED BY LE STATE
  
  Xenium_LE_BBB_dotplot <- ggplot(
    
    LE_BBB_dot_data,
    
    aes(
      x = Gene,
      y = id,
      size = pct.exp,
      color = avg.exp.scaled
    )
    
  ) +
    
    geom_point(
      alpha = 0.85
    ) +
    
    facet_wrap(
      ~Program,
      scales = "free_x",
      nrow = 1
    ) +
    
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = "Scaled\nexpression"
    ) +
    
    scale_size(
      range = c(
        1,
        8
      ),
      name = "% cells"
    ) +
    
    labs(
      x = NULL,
      y = "BBB stability group",
      title = "Xenium validation of GSEA leading-edge genes"
    ) +
    
    theme_bw(
      base_size = 12
    ) +
    
    theme(
      
      strip.text = element_text(
        face = "bold",
        size = 14
      ),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 11
      ),
      
      axis.text.y = element_text(
        size = 12
      ),
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      
      legend.position = "right"
      
    )
  
  
  Xenium_LE_BBB_dotplot
  
  
  # Save
  
  ggsave(
    
    filename =
      "New NVU dereg. framework/Xenium_FINAL_LE_BBB_stability_dotplot.svg",
    
    plot = Xenium_LE_BBB_dotplot,
    
    device = svglite,
    
    width = 12,
    
    height = 3.5,
    
    units = "in"
    
  )
  
}


# PART 2 — LE GENES ACROSS NVU STATES
# FACETED BY LGG / GBM / BrM

# 7. DEFINE NVU STATE ORDER

nvu_order <- c(
  
  "Canonical NVU (no tumor)",
  
  "EC/PC only (no tumor)",
  
  "Astro, no pericyte (no tumor)",
  
  "Bare EC vessel (no tumor)",
  
  "Intact NVU (+ tumor)",
  
  "Pericyte dropout (+ tumor)",
  
  "No astrocyte/end-feet (+ tumor)",
  
  "Bare vessel (+ tumor)"
  
)


ec_whole_tissue_noNA$long_label <- factor(
  
  ec_whole_tissue_noNA$long_label,
  
  levels = nvu_order
  
)


# 8. DEFINE TUMOR TYPE

ec_whole_tissue_tumor <- subset(
  
  ec_whole_tissue_noNA,
  
  subset = sample != "Normal"
)


ec_whole_tissue_tumor$Tumor <- case_when(
  
  ec_whole_tissue_tumor$sample == "Gr_3" ~
    "LGG",
  
  ec_whole_tissue_tumor$sample == "GBM" ~
    "GBM",
  
  ec_whole_tissue_tumor$sample %in% c(
    "Breast",
    "Kidney",
    "Lung",
    "Melanoma",
    "Thyroid"
  ) ~
    "BrM",
  
  TRUE ~
    NA_character_
  
)


ec_whole_tissue_tumor$Tumor <- factor(
  
  ec_whole_tissue_tumor$Tumor,
  
  levels = c(
    "LGG",
    "GBM",
    "BrM"
  )
  
)

# 9. FUNCTION TO EXTRACT DOTPLOT DATA

get_NVU_LE_dotplot_data <- function(
    object,
    genes,
    tumor_name
){
  
  if (
    length(genes) == 0 ||
    ncol(object) == 0
  ) {
    
    return(
      tibble()
    )
    
  }
  
  
  dp <- DotPlot(
    
    object,
    
    features = genes,
    
    group.by = "long_label",
    
    assay = "Xenium"
    
  )
  
  
  dp$data %>%
    
    mutate(
      
      Tumor = tumor_name,
      
      Gene = factor(
        features.plot,
        levels = genes
      ),
      
      id = factor(
        id,
        levels = rev(nvu_order)
      )
      
    )
  
}

# 10. SPLIT BY TUMOR TYPE

ec_LGG <- subset(
  
  ec_whole_tissue_tumor,
  
  subset = Tumor == "LGG"
)


ec_GBM <- subset(
  
  ec_whole_tissue_tumor,
  
  subset = Tumor == "GBM"
)


ec_BrM <- subset(
  
  ec_whole_tissue_tumor,
  
  subset = Tumor == "BrM"
)


# BBB-PRESERVED LE GENES ACROSS NVU STATES


high_NVU_dot_data <- bind_rows(
  
  get_NVU_LE_dotplot_data(
    ec_LGG,
    high_LE_genes_xenium,
    "LGG"
  ),
  
  get_NVU_LE_dotplot_data(
    ec_GBM,
    high_LE_genes_xenium,
    "GBM"
  ),
  
  get_NVU_LE_dotplot_data(
    ec_BrM,
    high_LE_genes_xenium,
    "BrM"
  )
  
)


high_NVU_dot_data$Tumor <- factor(
  
  high_NVU_dot_data$Tumor,
  
  levels = c(
    "LGG",
    "GBM",
    "BrM"
  )
  
)


if (
  nrow(high_NVU_dot_data) > 0
) {
  
  Xenium_BBB_preserved_NVU_plot <- ggplot(
    
    high_NVU_dot_data,
    
    aes(
      x = Gene,
      y = id,
      size = pct.exp,
      color = avg.exp.scaled
    )
    
  ) +
    
    geom_point() +
    
    facet_wrap(
      ~Tumor,
      nrow = 1
    ) +
    
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = "Scaled\nexpression"
    ) +
    
    scale_size(
      range = c(
        1,
        8
      ),
      name = "% cells"
    ) +
    
    labs(
      x = NULL,
      y = NULL,
      title = "BBB-preserved leading-edge genes across NVU states"
    ) +
    
    theme_bw(
      base_size = 12
    ) +
    
    theme(
      
      strip.text = element_text(
        face = "bold",
        size = 14
      ),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 9
      ),
      
      axis.text.y = element_text(
        size = 10
      ),
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      
      legend.position = "right"
      
    )
  
  
  Xenium_BBB_preserved_NVU_plot
  
  
  ggsave(
    
    filename =
      "New NVU dereg. framework/Xenium_FINAL_LE_BBB_preserved_NVU_states.svg",
    
    plot = Xenium_BBB_preserved_NVU_plot,
    
    device = svglite,
    
    width = 20,
    
    height = 4,
    
    units = "in"
    
  )
  
}


# BBB-DISRUPTED LE GENES ACROSS NVU STATES


low_NVU_dot_data <- bind_rows(
  
  get_NVU_LE_dotplot_data(
    ec_LGG,
    low_LE_genes_xenium,
    "LGG"
  ),
  
  get_NVU_LE_dotplot_data(
    ec_GBM,
    low_LE_genes_xenium,
    "GBM"
  ),
  
  get_NVU_LE_dotplot_data(
    ec_BrM,
    low_LE_genes_xenium,
    "BrM"
  )
  
)


low_NVU_dot_data$Tumor <- factor(
  
  low_NVU_dot_data$Tumor,
  
  levels = c(
    "LGG",
    "GBM",
    "BrM"
  )
  
)


if (
  nrow(low_NVU_dot_data) > 0
) {
  
  Xenium_BBB_disrupted_NVU_plot <- ggplot(
    
    low_NVU_dot_data,
    
    aes(
      x = Gene,
      y = id,
      size = pct.exp,
      color = avg.exp.scaled
    )
    
  ) +
    
    geom_point() +
    
    facet_wrap(
      ~Tumor,
      nrow = 1
    ) +
    
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = "Scaled\nexpression"
    ) +
    
    scale_size(
      range = c(
        1,
        8
      ),
      name = "% cells"
    ) +
    
    labs(
      x = NULL,
      y = NULL,
      title = "BBB-disrupted leading-edge genes across NVU states"
    ) +
    
    theme_bw(
      base_size = 12
    ) +
    
    theme(
      
      strip.text = element_text(
        face = "bold",
        size = 14
      ),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 9
      ),
      
      axis.text.y = element_text(
        size = 10
      ),
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      
      legend.position = "right"
      
    )
  
  
  Xenium_BBB_disrupted_NVU_plot
  
  
  ggsave(
    
    filename =
      "New NVU dereg. framework/Xenium_FINAL_LE_BBB_disrupted_NVU_states.svg",
    
    plot = Xenium_BBB_disrupted_NVU_plot,
    
    device = svglite,
    
    width = 8,
    
    height = 4,
    
    units = "in"
    
  )
  
}



