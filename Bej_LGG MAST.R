#Running MAST regression, GSEA, and leading edge analysis on Bejarano scRNA 
#LGG only


library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(tidyverse)
library(stringr)
library(Seurat)
library(SeuratObject)
install.packages("Matrix")
library(Matrix)



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(
  c(
    "MAST",
    "clusterProfiler",
    "org.Hs.eg.db",
    "enrichplot"
  ),
  ask = FALSE,
  update = FALSE
)

library(MAST)

library(RSQLite)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)




#Load RDS

Bej_seurat <- readRDS(
  "Bejarano et al. 2025 Immunity/Bejarano_seurat_modules.rds"
)


#Clean object

Bej_seurat_clean <- Bej_seurat

Bej_seurat_clean@reductions <- list()
Bej_seurat_clean@graphs <- list()


rownames(Bej_seurat_clean@meta.data) <- Cells(Bej_seurat_clean)



# Check LGG labels first

unique(Bej_seurat_clean$Disease.subtype)


# SUBSET LGG
# Replace with exact metadata labels from above

Bej_LGG <- subset(
  Bej_seurat_clean,
  subset = Disease.subtype %in% c(
    "IDH mut glioma"
  )
)



# Check patient distribution

table(Bej_LGG$Patient)


#***New***
#Normalize Data
Bej_LGG <- NormalizeData(
  Bej_LGG,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# BBB stability distribution

ggplot(
  Bej_LGG@meta.data,
  aes(
    x = Patient,
    y = BBB_Stability_Score
  )
) +
  geom_boxplot() +
  theme_classic()



#------------------------------------------------
# MAST
#------------------------------------------------


# Detection rate

detected_rate_LGG <- Matrix::colSums(
  GetAssayData(
    Bej_LGG,
    assay = "RNA",
    layer = "counts"
  ) > 0
) / nrow(Bej_LGG)


names(detected_rate_LGG) <- Cells(Bej_LGG)


Bej_LGG$detected_rate <- detected_rate_LGG



# Expression matrix

expr_LGG <- GetAssayData(
  Bej_LGG,
  assay = "RNA",
  layer = "data"
)



# Filter genes

keep_genes_LGG <- Matrix::rowSums(expr_LGG > 0) >
  (0.01 * ncol(expr_LGG))


expr_LGG <- expr_LGG[keep_genes_LGG,]



dim(expr_LGG)



# Dense matrix

expr_dense_LGG <- as.matrix(expr_LGG)



# Metadata

cdata_LGG <- Bej_LGG@meta.data


rownames(cdata_LGG) <- Cells(Bej_LGG)


cdata_LGG <- cdata_LGG[
  colnames(expr_dense_LGG),
  ,
  drop = FALSE
]


fdata_LGG <- data.frame(
  primerid = rownames(expr_dense_LGG),
  row.names = rownames(expr_dense_LGG)
)


cdata_LGG$wellKey <- rownames(cdata_LGG)



# Create MAST object
sca_LGG <- FromMatrix(
  exprsArray = expr_dense_LGG,
  cData = cdata_LGG,
  fData = fdata_LGG,
  check_sanity = FALSE
)



# factors

colData(sca_LGG)$Patient <- factor(
  colData(sca_LGG)$Patient
)



# Run MAST

zlm_BBB_LGG <- zlm(
  ~ BBB_Stability_Score + detected_rate,
  sca = sca_LGG,
  method = "bayesglm",
  ebayes = FALSE
)



# Extract results
extract_MAST <- function(model, coefficient){
  res <- summary(
    model,
    doLRT = coefficient
  )
  results <- merge(
    res$datatable[
      res$datatable$component == "H",
    ],
    res$datatable[
      res$datatable$component == "logFC",
    ],
    by = c("primerid","contrast")
  )
  results$FDR <- p.adjust(
    results$`Pr(>Chisq).x`,
    method = "BH"
  )
  results %>%
    arrange(FDR)
}


BBB_LGG_results <- extract_MAST(
  zlm_BBB_LGG,
  "BBB_Stability_Score"
)




# Volcano

ggplot(
  BBB_LGG_results,
  aes(
    x = coef.y,
    y = -log10(FDR)
  )
) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Association with BBB Stability",
    y = "-log10(FDR)"
  )



#------------------------------------------------
# GSEA
#------------------------------------------------


BBB_LGG_ranked <- BBB_LGG_results %>%
  arrange(desc(coef.y))


geneList_BBB_LGG <- BBB_LGG_ranked$coef.y


names(geneList_BBB_LGG) <- BBB_LGG_ranked$primerid


geneList_BBB_LGG <- sort(
  geneList_BBB_LGG,
  decreasing = TRUE
)



gsea_BBB_LGG <- gseGO(
  geneList = geneList_BBB_LGG,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH"
)



gsea_BBB_LGG_df <- as.data.frame(
  gsea_BBB_LGG
)

write.csv(gsea_BBB_LGG_df, "Bejarano et al. 2025 Immunity/Analysis/Bej_MAST_gsea_BBB_LGG_df_normalized.csv")


#------------------------------------------------
# Select LGG pathways AFTER looking at results
#------------------------------------------------


#Split High and Low Pathways
#high use pvalue instead of q value or pvaladj
LGG_highBBB_pathways <- gsea_BBB_LGG_df %>%
  filter(
    NES > 0,
    pvalue < 0.05
  ) %>%
  arrange(desc(NES))

#select key non-redundant terms

#selected key non-redundant terms

LGG_high_selected <- c(
  "transport across blood-brain barrier",
  "vascular transport",
  "positive regulation of bicellular tight junction assembly",
  "regulation of establishment of endothelial barrier",
  "xenobiotic transport",
  "neutral amino acid transport",
  "endothelial cell development",
  "blood vessel remodeling"
)



#low 
LGG_lowBBB_pathways <- gsea_BBB_LGG_df %>%
  filter(
    NES < 0,
    pvalue < 0.05
  ) %>%
  arrange(NES)


LGG_low_selected <- c(
  "regulation of endothelial cell migration",
  "positive regulation of endothelial cell chemotaxis",
  "positive regulation of vascular endothelial growth factor signaling pathway",
  "extracellular matrix constituent secretion",
  "regulation of extracellular matrix assembly",
  "negative regulation of cell-cell adhesion mediated by cadherin"
)



#------------------------------------------------
# Leading edge analysis
#------------------------------------------------


LGG_selected_pathways <- c(
  LGG_high_selected,
  LGG_low_selected
)



LGG_LE_combined <- gsea_BBB_LGG_df %>%
  filter(
    Description %in% LGG_selected_pathways
  ) %>%
  dplyr::select(
    Description,
    NES,
    core_enrichment
  ) %>%
  separate_rows(
    core_enrichment,
    sep="/"
  ) %>%
  rename(
    Gene = core_enrichment
  )



## 1. Select representative GSEA pathways

LGG_selected_pathways <- c(
  LGG_high_selected,
  LGG_low_selected
)

gsea_plot_df_LGG <- gsea_BBB_LGG_df %>%
  filter(
    Description %in% LGG_selected_pathways
  ) %>%
  mutate(
    Direction = ifelse(
      NES > 0,
      "High BBB stability",
      "Low BBB stability"
    )
  ) %>%
  arrange(desc(NES))


## 2. NES barplot

ggplot(
  gsea_plot_df_LGG,
  aes(
    x = reorder(Description, NES),
    y = NES,
    fill = Direction
  )
) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "High BBB stability" = "firebrick",
      "Low BBB stability" = "steelblue"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "Normalized Enrichment Score (NES)",
    fill = NULL
  )


## 3. Extract leading-edge genes from selected pathways

LGG_LE_combined <- gsea_BBB_LGG_df %>%
  filter(
    Description %in% LGG_selected_pathways
  ) %>%
  dplyr::select(
    Description,
    NES,
    core_enrichment
  ) %>%
  tidyr::separate_rows(
    core_enrichment,
    sep = "/"
  ) %>%
  dplyr::rename(
    Gene = core_enrichment
  )


## Remove mitochondrial genes

mt_genes <- c(
  "MT-CO1",
  "MT-CO2",
  "MT-CO3",
  "MT-ND1",
  "MT-ND2",
  "MT-ND4",
  "MT-ND5",
  "MT-ATP6",
  "MT-CYB"
)


LGG_LE_combined <- LGG_LE_combined %>%
  filter(
    !Gene %in% mt_genes
  )


## 4. Count gene recurrence across pathways

LE_gene_frequency_LGG <- LGG_LE_combined %>%
  dplyr::count(
    Gene,
    name = "Pathway_Count"
  )


## 5. Build plotting dataframe

LE_dotplot_LGG <- LGG_LE_combined %>%
  group_by(
    Description,
    Gene
  ) %>%
  summarise(
    NES = unique(NES)[1],
    .groups = "drop"
  ) %>%
  left_join(
    LE_gene_frequency_LGG,
    by = "Gene"
  )


## 6. Order pathways

high_pathways_LGG <- gsea_BBB_LGG_df %>%
  filter(
    Description %in% LGG_selected_pathways,
    NES > 0
  ) %>%
  arrange(
    desc(NES)
  ) %>%
  pull(Description)


low_pathways_LGG <- gsea_BBB_LGG_df %>%
  filter(
    Description %in% LGG_selected_pathways,
    NES < 0
  ) %>%
  arrange(
    NES
  ) %>%
  pull(Description)


LE_dotplot_LGG$Description <- factor(
  LE_dotplot_LGG$Description,
  levels = c(
    rev(low_pathways_LGG),
    rev(high_pathways_LGG)
  )
)


## 7. Order genes

gene_order_LGG <- LE_dotplot_LGG %>%
  group_by(
    Gene
  ) %>%
  summarise(
    meanNES = mean(NES),
    .groups = "drop"
  ) %>%
  arrange(
    desc(meanNES)
  ) %>%
  pull(Gene)


LE_dotplot_LGG$Gene <- factor(
  LE_dotplot_LGG$Gene,
  levels = gene_order_LGG
)


## 8. Automatically calculate divider position

n_high_genes_LGG <- LE_dotplot_LGG %>%
  filter(
    NES > 0
  ) %>%
  pull(Gene) %>%
  unique() %>%
  length()


## 9. Plot leading-edge dot plot

ggplot(
  LE_dotplot_LGG,
  aes(
    x = Gene,
    y = Description,
    size = Pathway_Count,
    color = NES
  )
) +
  geom_point(
    alpha = 0.9
  ) +
  geom_vline(
    xintercept = n_high_genes_LGG + 0.5,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_color_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    name = "NES"
  ) +
  scale_size_continuous(
    range = c(2,8),
    name = "Leading-edge\npathway overlap"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 10
    ),
    legend.position = "right"
  ) +
  labs(
    x = "Leading-edge genes",
    y = NULL
  ) +
  ggtitle(
    "LGG High vs Low BBB Stability Pathways & Genes"
  )
