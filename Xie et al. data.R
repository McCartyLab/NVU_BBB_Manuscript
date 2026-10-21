#Xie et al. data


library(dplyr)
library(tidyr)
library(purrr)
library(csv)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggrepel)
library(tidyverse)
library(stringr)
library(Seurat)
library(SeuratObject)
library(hdf5r)
library(scCustomize)
library(arrow)



#load data
Xie_vascular_seurat <- readRDS("Xie et al. 2024 Neuron/GBM_vascular_cells.Rds")



#Plotting data
DimPlot(Xie_vascular_seurat, group.by = "tumour.type")
DimPlot(Xie_vascular_seurat, group.by = "cell.type")

DimPlot(Xie_vascular_immune_seurat, group.by = "tumour.type")



#subset endothelial cells
Xie_ec <- subset(Xie_vascular_seurat,
                 subset = cell.type == "Endothelial")


Idents(Xie_ec) <- "tumour.type"

DimPlot(Xie_ec, group.by = "tumour.type")

Xie_high_vs_control <- FindMarkers(
  Xie_ec,
  ident.1 = "Higher-grade glioma",
  ident.2 = "Control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25
)


DotPlot(Xie_ec,
        idents = c("Higher-grade glioma", "Control"),
        features = c("PLVAP", "LAMB1", "COL4A1",
                     "LDLR", "SLC16A1")) +
  coord_flip() +
  ggtitle("EC DEGs: GBM vs Normal (Xie et al. 2024)")





# Xie et al. EC scRNA-seq
# NVU Module Scores and BBB Stability Score
# Publication-ready plots

library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

# 1. Define NVU modules

nvu_modules <- list(
  Tight_Junctions = c(
    "CLDN5", "OCLN", "TJP1", "TJP2", "TJP3",
    "F11R", "JAM2", "JAM3", "LSR"
  ),
  
  BBB_Transport = c(
    "ABCB1", "ABCG2", "SLC2A1", "SLCO1A2",
    "SLCO1C1", "LRP1", "MFSD2A"
  ),
  
  Angiogenesis = c(
    "ANGPT2", "VEGFA", "DLL4", "NDP", "KDR", "TEK"
  ),
  
  Vascular_Inflam = c(
    "ICAM1", "VCAM1", "IL1B", "IL6", "CCL2", "S1PR3"
  ),
  
  ECM_Remodeling = c(
    "COL4A1", "LAMB1", "FN1", "VCL",
    "ITGB1", "MMP2", "MMP9", "VIM"
  )
)


# 2. Calculate UCell module scores

Xie_ec <- AddModuleScore_UCell(
  Xie_ec,
  features = nvu_modules
)


# Rename UCell metadata columns
colnames(Xie_ec@meta.data)[
  colnames(Xie_ec@meta.data) %in%
    paste0(names(nvu_modules), "_UCell")
] <- names(nvu_modules)


# 3. Calculate BBB Stability Score

Xie_ec$BBB_Stability_Score <-
  (
    Xie_ec$Tight_Junctions +
      Xie_ec$BBB_Transport
  ) -
  (
    Xie_ec$Angiogenesis +
      Xie_ec$Vascular_Inflam +
      Xie_ec$ECM_Remodeling
  )


# 4. Calculate sample-level scores

Xie_meta <- Xie_ec@meta.data

Xie_ec_summary <- Xie_meta %>%
  group_by(sampleID, tumour.type) %>%
  summarise(
    Tight_Junctions = mean(
      Tight_Junctions,
      na.rm = TRUE
    ),
    
    BBB_Transport = mean(
      BBB_Transport,
      na.rm = TRUE
    ),
    
    Angiogenesis = mean(
      Angiogenesis,
      na.rm = TRUE
    ),
    
    Vascular_Inflam = mean(
      Vascular_Inflam,
      na.rm = TRUE
    ),
    
    ECM_Remodeling = mean(
      ECM_Remodeling,
      na.rm = TRUE
    ),
    
    BBB_Stability = mean(
      BBB_Stability_Score,
      na.rm = TRUE
    ),
    
    n_cells = n(),
    
    .groups = "drop"
  )

# 5. Set disease group order

Xie_ec_summary <- Xie_ec_summary %>%
  mutate(
    tumour.type = factor(
      tumour.type,
      levels = c(
        "Control",
        "Lower-grade glioma",
        "Higher-grade glioma"
      ),
      labels = c(
        "Normal",
        "LGG",
        "GBM"
      )
    )
  )

# 6. Convert module scores to long format

Xie_ec_long <- Xie_ec_summary %>%
  pivot_longer(
    cols = c(
      Tight_Junctions,
      BBB_Transport,
      Angiogenesis,
      Vascular_Inflam,
      ECM_Remodeling
    ),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  mutate(
    Module = factor(
      Module,
      levels = c(
        "Tight_Junctions",
        "BBB_Transport",
        "Angiogenesis",
        "Vascular_Inflam",
        "ECM_Remodeling"
      ),
      labels = c(
        "Tight Junctions",
        "BBB Transport",
        "Angiogenesis",
        "Vascular Inflammation",
        "ECM Remodeling"
      )
    )
  )

# 7. Publication-ready NVU module plot

p_modules <- ggplot(
  Xie_ec_long,
  aes(
    x = tumour.type,
    y = Score
  )
) +
  
  # Individual samples
  geom_quasirandom(
    aes(group = tumour.type),
    width = 0.15,
    size = 1.8,
    alpha = 0.45,
    color = "black"
  ) +
  
  # Boxplot on top
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ) +
  # One panel per module
  facet_wrap(
    ~Module,
    scales = "free_y",
    nrow = 1
  ) +
  
  labs(
    title = "Xiet et al. EC scRNA-seq",
    x = NULL,
    y = "Mean UCell score"
  ) +
  
  theme_classic(base_size = 13) +
  theme( text = element_text(family = "Helvetica"), axis.text = element_text(family = "Helvetica"), axis.title = element_text(family = "Helvetica"), legend.text = element_text(family = "Helvetica"), legend.title = element_text(family = "Helvetica") ) +
  theme(
    plot.title = element_text(
      size = 28,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      size = 14,
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 14,
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 14,
      face = "bold"
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    
    panel.spacing = unit(
      1.2,
      "lines"
    ),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )


p_modules

# Save
output_dir <- "Vascular subtype analysis/Figures/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(
  filename = file.path(
    output_dir,
    "Xie_p_modules.jpeg"
  ),
  plot = p_modules,
  width = 12,
  height = 6,
  dpi = 600
)




# 8. Publication-ready BBB Stability Score plot

p_BBB <- ggplot(
  Xie_ec_summary,
  aes(
    x = tumour.type,
    y = BBB_Stability
  )
) +
  
  # Individual samples
  geom_quasirandom(
    aes(group = tumour.type),
    width = 0.15,
    size = 1.8,
    alpha = 0.45,
    color = "black"
  ) +
  
  # Boxplot on top
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ) +
  
  labs(
    title = "Xiet et al. EC scRNA-seq",
    x = NULL,
    y = "BBB Stability Score"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      size = 28,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      size = 14,
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 14,
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 14,
      face = "bold"
    ),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )


p_BBB


# Save
output_dir <- "Vascular subtype analysis/Figures/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(
  filename = file.path(
    output_dir,
    "Xie_p_BBB_Stability.jpeg"
  ),
  plot = p_BBB,
  width = 6,
  height = 5.5,
  dpi = 600
)






#Final fig----
# 7. Publication-ready NVU module + BBB Stability plot

library(ggplot2)
library(ggbeeswarm)
library(dplyr)

# ------------------------------------------------------------
# Prepare module data
# ------------------------------------------------------------

Xie_ec_long <- Xie_ec_long %>%
  mutate(
    Module = factor(
      Module,
      levels = c(
        "Tight Junctions",
        "BBB Transport",
        "Angiogenesis",
        "Vascular Inflammation",
        "ECM Remodeling"
      ),
      labels = c(
        "Tight Junctions",
        "BBB Transport",
        "Angiogenesis",
        "Vascular Inflammation",
        "ECM Remodeling"
      )
    )
  )


# ------------------------------------------------------------
# Prepare BBB Stability data
# ------------------------------------------------------------

Xie_BBB_long <- Xie_ec_summary %>%
  transmute(
    tumour.type = tumour.type,
    Score = BBB_Stability,
    Module = "BBB Stability"
  ) %>%
  mutate(
    Module = factor(
      Module,
      levels = c(
        "Tight Junctions",
        "BBB Transport",
        "Angiogenesis",
        "Vascular Inflammation",
        "ECM Remodeling",
        "BBB Stability"
      )
    )
  )


# ------------------------------------------------------------
# Combine module + BBB Stability data
# ------------------------------------------------------------

Xie_combined <- bind_rows(
  Xie_ec_long %>%
    select(
      tumour.type,
      Score,
      Module
    ),
  Xie_BBB_long
) %>%
  mutate(
    Module = factor(
      Module,
      levels = c(
        "Tight Junctions",
        "BBB Transport",
        "Angiogenesis",
        "Vascular Inflammation",
        "ECM Remodeling",
        "BBB Stability"
      )
    )
  )


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p_Xie_NVU <- ggplot(
  Xie_combined,
  aes(
    x = tumour.type,
    y = Score
  )
) +
  
  # Individual samples
  geom_quasirandom(
    aes(group = tumour.type),
    width = 0.15,
    size = 1.5,
    alpha = 0.45,
    color = "black"
  ) +
  
  # Boxplot
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ) +
  
  # 2 columns × 3 rows
  facet_wrap(
    ~Module,
    ncol = 2,
    nrow = 3,
    scales = "free"
  ) +
  
  labs(
    title = "Xie et al. EC scRNA-seq",
    x = NULL,
    y = "Mean UCell score"
  ) +
  
  theme_classic(base_size = 13) +
  theme( text = element_text(family = "Helvetica"), axis.text = element_text(family = "Helvetica"), axis.title = element_text(family = "Helvetica"), legend.text = element_text(family = "Helvetica"), legend.title = element_text(family = "Helvetica") ) +
  theme(
    plot.title = element_text(
      size = 32,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      size = 20,
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 16,
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 20,
      face = "bold"
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 22,
      face = "bold",
      color = "black"
    ),
    
    panel.spacing = unit(
      1.2,
      "lines"
    ),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )


p_Xie_NVU


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

output_dir <- "Vascular subtype analysis/Figures/"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

ggsave(
  filename = file.path(
    output_dir,
    "Helv_Xie_NVU_modules_BBB_Stability.jpeg"
  ),
  plot = p_Xie_NVU,
  width = 14,
  height = 13,
  dpi = 600
)






#quantiles

#Identify top and bottom 25%
library(dplyr)


# Calculate subtype-specific quartiles and assign groups
Xie_meta <- Xie_meta %>%
  group_by(tumour.type) %>%
  mutate(
    q25 = quantile(BBB_Stability_Score, 0.25, na.rm = TRUE),
    q75 = quantile(BBB_Stability_Score, 0.75, na.rm = TRUE),
    BBB_group = case_when(
      BBB_Stability_Score >= q75 ~ "High_BBB_Stability",
      BBB_Stability_Score <= q25 ~ "Low_BBB_Stability",
      TRUE ~ "Intermediate"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(-q25, -q75)

# Add the grouping back to the Seurat object
Xie_ec$BBB_group <- Xie_meta$BBB_group

# Check numbers
table(Xie_ec$tumour.type,
      Xie_ec$BBB_group)


DimPlot(Xie_ec, group.by = "BBB_group")

VlnPlot(Xie_ec, features = c("BBB_Stability_Score"))





#Check MAST results


LE_dotplot <- read.csv("LE_dotplot_GBM.csv")

#GBM
Idents(Xie_ec) <- "tumour.type"
Xie_ec_GBM <- subset(Xie_ec, ident = "Higher-grade glioma")


#Genes up in high bbb stability cells
high_genes_to_plot <- LE_dotplot %>%
  filter(NES > 0) %>%
  pull(Gene) %>%
  unique() %>%
  na.omit()

DotPlot(
  Xie_ec_GBM,
  group.by = "BBB_group",
  features = high_genes_to_plot
) +
  coord_flip()


#Genes up in low bbb stability cells
low_genes_to_plot <- LE_dotplot %>%
  filter(NES < 0) %>%
  pull(Gene) %>%
  unique() %>%
  na.omit()

DotPlot(
  Xie_ec_GBM,
  group.by = "BBB_group",
  features = low_genes_to_plot
) +
  coord_flip()



#Genes high and low in Xenium ECs

# Xiet et al. EC scRNA-seq: BBB-associated genes
# Tumor-associated vs. Normal-associated genes

library(Seurat)
library(ggplot2)
library(dplyr)

# Tumor-associated genes
tumor_associated_genes <- unique(c(
  "PLVAP",
  "LAMB1",
  "KDR",
  "MMP2",
  "MMP9",
  "IL1B",
  "SLC4A4",
  "MCAM",
  "FZD4",
  "COL4A1"
))

# Normal-associated genes
normal_associated_genes <- unique(c(
  "SLC2A3",
  "SLC1A1",
  "LDLR",
  "TFRC",
  "CLDN5",
  "OCLN",
  "IL6",
  "VCAM1",
  "F11R",
  "DLL4",
  "MFSD2A",
  "ABCG2",
  "ABCB1",
  "SLCO1A2",
  "SLC16A1"
))

# Remove any genes appearing in both groups
normal_associated_genes <- setdiff(
  normal_associated_genes,
  tumor_associated_genes
)

# Combined gene list
genes_to_plot <- c(
  tumor_associated_genes,
  normal_associated_genes
)

# Generate DotPlot
p <- DotPlot(
  Xie_ec_GBM,
  group.by = "BBB_group",
  features = genes_to_plot,
  dot.scale = 7
)

# Extract DotPlot data
dot_data <- p$data

# Add gene group
dot_data$GeneGroup <- ifelse(
  as.character(dot_data$features) %in% tumor_associated_genes,
  "Tumor-associated",
  "Normal-associated"
)

# Set facet order
dot_data$GeneGroup <- factor(
  dot_data$GeneGroup,
  levels = c(
    "Tumor-associated",
    "Normal-associated"
  )
)

# Preserve gene order
dot_data$features <- factor(
  as.character(dot_data$features),
  levels = rev(genes_to_plot)
)

# Replace underscores with spaces in BBB_group labels 
dot_data$id <- str_replace_all( as.character(dot_data$id), "_", " " )
# Publication-ready DotPlot
Xie_dotplot<- ggplot(
  dot_data,
  aes(
    x = id,
    y = features,
    size = pct.exp,
    color = avg.exp.scaled
  )
) +
  geom_point() +
  
  # Separate tumor-associated and normal-associated genes
  facet_grid(
    rows = vars(GeneGroup),
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_size(
    range = c(1.5, 7),
    name = "Percent\nexpressed"
  ) +
  
  scale_color_gradient(
    low = "grey90",
    high = "#B2182B",
    name = "Average\nexpression"
  ) +
  
  labs(
    title = "Xie et al. EC scRNA-seq",
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      size = 14,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 12,
      color = "black"
    ),
    
    strip.text.y = element_text(
      size = 15,
      face = "bold"
    ),
    
    strip.background = element_blank(),
    
    panel.spacing = unit(
      1,
      "lines"
    ),
    
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )

# Save
output_dir <- "Vascular subtype analysis/Figures/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(
  filename = file.path(
    output_dir,
    "Xie_dotplot.svg"
  ),
  plot = Xie_dotplot,
  width = 8,
  height = 8,
  dpi = 600
)






#Check VEGF and Wnt pathway signaling in high and low BBB stability ECs----

install.packages("remotes")

remotes::install_github(
  "jinworks/CellChat",
  dependencies = TRUE
)

library(CellChat)

CellChatDB <- CellChatDB.human

# See available pathway names
sort(unique(CellChatDB$interaction$pathway_name))

#pull pathways of interest
vegf_lr <- CellChatDB$interaction %>%
  filter(pathway_name == "VEGF")

wnt_lr <- CellChatDB$interaction %>%
  filter(pathway_name == "WNT")

vegf_lr
wnt_lr

#pull gene names
vegf_genes <- extractGeneSubsetFromPair(
  pairLR = vegf_lr,
  complex_input = CellChatDB$complex,
  geneInfo = CellChatDB$geneInfo,
  combined = TRUE
)

wnt_genes <- extractGeneSubsetFromPair(
  pairLR = wnt_lr,
  complex_input = CellChatDB$complex,
  geneInfo = CellChatDB$geneInfo,
  combined = TRUE
)

vegf_genes <- unique(vegf_genes)
wnt_genes <- unique(wnt_genes)

vegf_genes
wnt_genes

#check if genes in Xie data, should be all present
vegf_genes_xie <- intersect(
  vegf_genes,
  rownames(Xie_ec)
)

wnt_genes_xie <- intersect(
  wnt_genes,
  rownames(Xie_ec)
)

cat("VEGF genes in Xie:", length(vegf_genes_xie), "\n")
cat("WNT genes in Xie:", length(wnt_genes_xie), "\n")

vegf_genes_xie
wnt_genes_xie


#Score pathways using Module
library(UCell)

Xie_ec <- AddModuleScore_UCell(
  Xie_ec,
  features = list(
    VEGF_CellChat = vegf_genes_xie,
    WNT_CellChat = wnt_genes_xie
  )
)


#plot
plot_df <- Xie_ec@meta.data %>%
  select(
    BBB_group,
    tumour.type,
    VEGF_CellChat_UCell,
    WNT_CellChat_UCell
  ) %>%
  mutate(
    Disease = recode(
      tumour.type,
      "Control" = "Normal",
      "Lower-grade glioma" = "LGG",
      "Higher-grade glioma" = "GBM"
    ),
    Disease = factor(
      Disease,
      levels = c("Normal", "LGG", "GBM")
    )
  ) %>%
  pivot_longer(
    cols = c(
      VEGF_CellChat_UCell,
      WNT_CellChat_UCell
    ),
    names_to = "Pathway",
    values_to = "Score"
  ) %>%
  mutate(
    Pathway = recode(
      Pathway,
      VEGF_CellChat_UCell = "VEGF",
      WNT_CellChat_UCell = "WNT"
    ),
    BBB_group = factor(
      BBB_group,
      levels = c(
        "High_BBB_Stability",
        "Intermediate",
        "Low_BBB_Stability"
      )
    )
  )

ggplot(
  plot_df,
  aes(
    x = BBB_group,
    y = Score
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA
  ) +
  facet_grid(
    Pathway ~ Disease,
    scales = "free_y"
  ) +
  ggtitle("ECs: Xie et al. scRNA") +
  theme_classic() +
  labs(
    x = NULL,
    y = "UCell score"
  )


