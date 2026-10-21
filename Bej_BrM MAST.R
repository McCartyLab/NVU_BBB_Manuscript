############################################################
## Brain metastasis BBB stability GSEA + leading edge analysis
############################################################


library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


############################################################
## 1. Subset brain metastasis endothelial cells
############################################################


#Load RDS
Bej_seurat <- readRDS("Bejarano et al. 2025 Immunity/Bejarano_seurat_modules.rds")

Bej_seurat_clean <- Bej_seurat

Bej_seurat_clean@reductions <- list()
Bej_seurat_clean@graphs <- list()

#align metadata rownames to seurat cell names
rownames(Bej_seurat@meta.data) <- Cells(Bej_seurat)


Bej_BrM <- subset(
  Bej_seurat_clean,
  subset = Disease.subtype == "BrM"
)

#***New***
#Normalize Data
Bej_BrM <- NormalizeData(
  Bej_BrM,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


############################################################
## 2. Detection rate
############################################################

detected_rate_BrM <- Matrix::colSums(
  GetAssayData(
    Bej_BrM,
    assay = "RNA",
    layer = "counts"
  ) > 0
) / nrow(Bej_BrM)


names(detected_rate_BrM) <- Cells(Bej_BrM)

Bej_BrM$detected_rate <- detected_rate_BrM



############################################################
## 3. Extract expression matrix
############################################################

expr_BrM <- GetAssayData(
  Bej_BrM,
  assay = "RNA",
  layer = "data"
)


keep_genes_BrM <- Matrix::rowSums(expr_BrM > 0) >
  (0.01 * ncol(expr_BrM))


expr_BrM <- expr_BrM[keep_genes_BrM, ]



############################################################
## 4. Create MAST object
############################################################

expr_dense_BrM <- as.matrix(expr_BrM)


cdata_BrM <- Bej_BrM@meta.data

rownames(cdata_BrM) <- Cells(Bej_BrM)

cdata_BrM <- cdata_BrM[
  colnames(expr_dense_BrM),
  ,
  drop = FALSE
]


fdata_BrM <- data.frame(
  primerid = rownames(expr_dense_BrM),
  row.names = rownames(expr_dense_BrM)
)


cdata_BrM$wellKey <- rownames(cdata_BrM)


sca_BrM <- FromMatrix(
  exprsArray = expr_dense_BrM,
  cData = cdata_BrM,
  fData = fdata_BrM,
  check_sanity = FALSE
)



############################################################
## 5. MAST model
############################################################

zlm_BBB_BrM <- zlm(
  ~ BBB_Stability_Score + detected_rate,
  sca = sca_BrM,
  method = "bayesglm",
  ebayes = FALSE
)

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


BBB_BrM_results <- extract_MAST(
  zlm_BBB_BrM,
  "BBB_Stability_Score"
)



############################################################
## 6. Rank genes for GSEA
############################################################

BBB_BrM_ranked <- BBB_BrM_results %>%
  arrange(desc(coef.y))


geneList_BBB_BrM <- BBB_BrM_ranked$coef.y

names(geneList_BBB_BrM) <- BBB_BrM_ranked$primerid


geneList_BBB_BrM <- sort(
  geneList_BBB_BrM,
  decreasing = TRUE
)



############################################################
## 7. GO GSEA
############################################################

gsea_BBB_BrM <- gseGO(
  geneList = geneList_BBB_BrM,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH"
)


gsea_BBB_BrM_df <- as.data.frame(
  gsea_BBB_BrM
)

write.csv(gsea_BBB_BrM_df, "Bejarano et al. 2025 Immunity/Analysis/Bej_MAST_gsea_BBB_BrM_df_normalized.csv")


############################################################
## 8. Select BrM representative pathways
############################################################

BrM_highBBB_pathways <- gsea_BBB_BrM_df %>%
  filter(
    NES > 0,
    pvalue < 0.05
  ) %>%
  arrange(desc(NES))



BrM_lowBBB_pathways <- gsea_BBB_BrM_df %>%
  filter(
    NES < 0,
    pvalue < 0.05
  ) %>%
  arrange(NES)
############################################################
## Clean BrM representative pathways + leading-edge dot plot
############################################################

############################################################
## 1. Select representative BrM pathways
############################################################

BrM_high_selected <- c(
  "maintenance of blood-brain barrier",
  "establishment of blood-brain barrier",
  "regulation of establishment of endothelial barrier",
  "tight junction assembly",
  "positive regulation of bicellular tight junction assembly",
  "negative regulation of vascular permeability",
  "transport across blood-brain barrier",
  "vascular transport",
  "organic anion transport",
  "transferrin transport",
  "detoxification",
  "stress response to metal ion"
)


BrM_low_selected <- c(
  "positive regulation of endothelial cell chemotaxis",
  "regulation of chemokine-mediated signaling pathway",
  "oxidative phosphorylation",
  "mitochondrial electron transport, ubiquinol to cytochrome c",
  "response to hydroperoxide",
  "nucleotide biosynthetic process"
)


BrM_selected_pathways <- c(
  BrM_high_selected,
  BrM_low_selected
)



############################################################
## 2. NES barplot
############################################################

gsea_plot_BrM <- gsea_BBB_BrM_df %>%
  filter(
    Description %in% BrM_selected_pathways
  ) %>%
  mutate(
    Direction = ifelse(
      NES > 0,
      "High BBB stability",
      "Low BBB stability"
    )
  )


ggplot(
  gsea_plot_BrM,
  aes(
    x = reorder(Description, NES),
    y = NES,
    fill = Direction
  )
) +
  geom_col(
    width = 0.75
  ) +
  geom_hline(
    yintercept = 0
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "High BBB stability" = "firebrick",
      "Low BBB stability" = "steelblue"
    )
  ) +
  theme_classic(
    base_size = 14
  ) +
  labs(
    x = NULL,
    y = "Normalized Enrichment Score (NES)",
    fill = NULL
  )



############################################################
## 3. Extract leading-edge genes
############################################################

BrM_LE_combined <- gsea_BBB_BrM_df %>%
  filter(
    Description %in% BrM_selected_pathways
  ) %>%
  dplyr::select(
    Description,
    NES,
    core_enrichment
  ) %>%
  mutate(
    core_enrichment = as.character(core_enrichment)
  ) %>%
  separate_rows(
    core_enrichment,
    sep = "/"
  ) %>%
  rename(
    Gene = core_enrichment
  ) %>%
  mutate(
    Gene = as.character(Gene)
  ) %>%
  filter(
    Gene != "",
    !grepl("^MT-", Gene)
  )


# Check structure
str(BrM_LE_combined)



############################################################
## 4. Count gene recurrence across pathways
############################################################

LE_gene_frequency_BrM <- BrM_LE_combined %>%
  dplyr::group_by(
    Gene
  ) %>%
  dplyr::summarise(
    Pathway_Count = n(),
    .groups = "drop"
  )



############################################################
## 5. Build dotplot dataframe
############################################################

LE_dotplot_BrM <- BrM_LE_combined %>%
  group_by(
    Description,
    Gene
  ) %>%
  summarise(
    NES = dplyr::first(NES),
    .groups = "drop"
  ) %>%
  left_join(
    LE_gene_frequency_BrM,
    by = "Gene"
  )



############################################################
## 6. Order pathways
############################################################

high_pathways_BrM <- gsea_BBB_BrM_df %>%
  filter(
    Description %in% BrM_selected_pathways,
    NES > 0
  ) %>%
  arrange(
    desc(NES)
  ) %>%
  pull(Description)


low_pathways_BrM <- gsea_BBB_BrM_df %>%
  filter(
    Description %in% BrM_selected_pathways,
    NES < 0
  ) %>%
  arrange(
    NES
  ) %>%
  pull(Description)


LE_dotplot_BrM$Description <- factor(
  LE_dotplot_BrM$Description,
  levels = c(
    rev(low_pathways_BrM),
    rev(high_pathways_BrM)
  )
)



############################################################
## 7. Order genes
############################################################

gene_order_BrM <- LE_dotplot_BrM %>%
  group_by(
    Gene
  ) %>%
  summarise(
    meanNES = mean(NES),
    overlap = max(Pathway_Count),
    .groups = "drop"
  ) %>%
  arrange(
    desc(meanNES),
    desc(overlap)
  ) %>%
  pull(Gene)


LE_dotplot_BrM$Gene <- factor(
  LE_dotplot_BrM$Gene,
  levels = gene_order_BrM
)



############################################################
## 8. Leading-edge dot plot
############################################################

ggplot(
  LE_dotplot_BrM,
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
  ggtitle("Mets High vs Low BBB Stability Pathways & Genes")


