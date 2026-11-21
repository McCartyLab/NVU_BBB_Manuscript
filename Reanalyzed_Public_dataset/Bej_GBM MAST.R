#Running MAST regression, GSEA, and leading edge analysis on Bejarano scRNA 
#GBM only


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





#MAST Model---

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("MAST", force = T)
library(Matrix)
library(MAST)

remove.packages("arrow")
install.packages("arrow")

library(arrow)
library(MAST)


#Load RDS
Bej_seurat <- readRDS("Bejarano et al. 2025 Immunity/Bejarano_seurat_modules.rds")

Bej_seurat_clean <- Bej_seurat

Bej_seurat_clean@reductions <- list()
Bej_seurat_clean@graphs <- list()

#align metadata rownames to seurat cell names
rownames(Bej_seurat@meta.data) <- Cells(Bej_seurat)

Bej_GBM <- subset(Bej_seurat_clean, 
                  subset = Disease.subtype == "IDH WT glioblastoma")


ggplot(
  Bej_GBM@meta.data,
  aes(
    x = Patient,
    y = BBB_Stability_Score
  )
) +
  geom_boxplot() +
  theme_classic() 

#***New***
#Normalize Data
Bej_GBM <- NormalizeData(
  Bej_GBM,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)



#Calculate detection rate
detected_rate_GBM <- Matrix::colSums(
  GetAssayData(
    Bej_GBM,
    assay = "RNA",
    layer = "counts"
  ) > 0
) / nrow(Bej_GBM)


names(detected_rate_GBM) <- Cells(Bej_GBM)

Bej_GBM$detected_rate <- detected_rate_GBM


#extract expression Matrix
expr_GBM <- GetAssayData(
  Bej_GBM,
  assay = "RNA",
  layer = "data"
)


#filter genes
keep_genes_GBM <- Matrix::rowSums(expr_GBM > 0) >
  (0.01 * ncol(expr_GBM))


expr_GBM <- expr_GBM[keep_genes_GBM,]


dim(expr_GBM)


#convert to dense matrix
expr_dense_GBM <- as.matrix(expr_GBM)


#Create MAST sca object
cdata_GBM <- Bej_GBM@meta.data


rownames(cdata_GBM) <- Cells(Bej_GBM)


cdata_GBM <- cdata_GBM[
  colnames(expr_dense_GBM),
  ,
  drop = FALSE
]


fdata_GBM <- data.frame(
  primerid = rownames(expr_dense_GBM),
  row.names = rownames(expr_dense_GBM)
)


cdata_GBM$wellKey <- rownames(cdata_GBM)


sca_GBM <- FromMatrix(
  exprsArray = expr_dense_GBM,
  cData = cdata_GBM,
  fData = fdata_GBM,
  check_sanity = FALSE
)

#set factors
colData(sca_GBM)$Patient <- factor(
  colData(sca_GBM)$Patient
)

#Run MAST model, w/o patient as factor bc only 2
zlm_BBB_GBM <- zlm(
  ~ BBB_Stability_Score + detected_rate,
  sca = sca_GBM,
  method = "bayesglm",
  ebayes = FALSE
)


#extract MAST results
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


BBB_GBM_results <- extract_MAST(
  zlm_BBB_GBM,
  "BBB_Stability_Score"
)



#volcano plot
ggplot(
  BBB_GBM_results,
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


#rank genes for GSEA
BBB_GBM_ranked <- BBB_GBM_results %>%
  arrange(desc(coef.y))


geneList_BBB_GBM <- BBB_GBM_ranked$coef.y

names(geneList_BBB_GBM) <- BBB_GBM_ranked$primerid


geneList_BBB_GBM <- sort(
  geneList_BBB_GBM,
  decreasing = TRUE
)


#Run GO GSEA
library(clusterProfiler)
library(org.Hs.eg.db)


gsea_BBB_GBM <- gseGO(
  geneList = geneList_BBB_GBM,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH"
)


gsea_BBB_GBM_df <- as.data.frame(gsea_BBB_GBM)

write.csv(gsea_BBB_GBM_df, "Bejarano et al. 2025 Immunity/Analysis/Bej_MAST_gsea_BBB_GBM_df_normalized.csv")








#Split High and Low Pathways
  #high use pvalue instead of q value or pvaladj
GBM_highBBB_pathways <- gsea_BBB_GBM_df %>%
  filter(
    NES > 0,
    pvalue < 0.05
  ) %>%
  arrange(desc(NES))

  #select key non-redundant terms

GBM_high_selected <- c(
  "transport across blood-brain barrier",
  "vascular transport",
  "organic anion transport",
  "detoxification",
  "stress response to metal ion",
  "positive regulation of bicellular tight junction assembly",
  "interferon-mediated signaling pathway"
)

#"Leading-edge analysis identified candidate endothelial programs associated with increased BBB stability, including transporter expression, barrier maintenance, metal ion detoxification, and endothelial stress-response pathways."
  
#low 
GBM_lowBBB_pathways <- gsea_BBB_GBM_df %>%
  filter(
    NES < 0,
    pvalue < 0.05
  ) %>%
  arrange(NES)

#selected key non-redundant terms

GBM_low_selected <- c(
  "oxidative phosphorylation",
  "mitochondrial electron transport, cytochrome c to oxygen",
  "endothelial cell chemotaxis",
  "regulation of extracellular matrix assembly",
  "negative regulation of cell-cell adhesion mediated by cadherin",
  "positive regulation of neuroinflammatory response"
)





##Leading Edge Analysis----

# High BBB Stability

GBM_high_LE <- gsea_BBB_GBM_df %>%
  filter(
    Description %in% GBM_high_selected #Use only selected representative pathways
  ) %>%
  dplyr::select(
    Description,
    NES,
    core_enrichment
  ) %>%
  separate_rows(
    core_enrichment,
    sep = "/"
  ) %>%
  rename(
    Gene = core_enrichment
  )


# Summarize recurring leading edge genes

GBM_high_LE_summary <- GBM_high_LE %>%
  group_by(Gene) %>%
  summarise(
    n_pathways = n(),
    pathways = paste(unique(Description), collapse = "; "),
    mean_NES = mean(NES)
  ) %>%
  arrange(
    desc(n_pathways),
    desc(mean_NES)
  )


head(
  GBM_high_LE_summary,
  30
)

GBM_high_specific <- GBM_high_LE_summary %>%
  filter(n_pathways >= 2)


#cnet plot for high

# IDs of selected pathways
GBM_high_cnet_ids <- gsea_BBB_GBM_df %>%
  filter(
    Description %in% GBM_high_selected
  ) %>%
  pull(ID)


# manually subset the gseaResult object
GBM_high_gsea <- gsea_BBB_GBM

GBM_high_gsea@result <- GBM_high_gsea@result %>%
  filter(ID %in% GBM_high_cnet_ids)

cnetplot(
  GBM_high_gsea,
  foldChange = geneList_BBB,
  showCategory = nrow(GBM_high_gsea@result)
)



# Low BBB Stability

GBM_low_LE <- gsea_BBB_GBM_df %>%
  filter(
    Description %in% GBM_low_selected #Use only selected representative pathways
  ) %>%
  dplyr::select(
    Description,
    NES,
    core_enrichment
  ) %>%
  separate_rows(
    core_enrichment,
    sep = "/"
  ) %>%
  rename(
    Gene = core_enrichment
  )


GBM_low_LE_summary <- GBM_low_LE %>%
  group_by(Gene) %>%
  summarise(
    n_pathways = n(),
    pathways = paste(unique(Description), collapse = "; "),
    mean_NES = mean(NES)
  ) %>%
  arrange(
    desc(n_pathways),
    desc(mean_NES)
  )


head(
  GBM_low_LE_summary,
  30
)

GBM_low_specific <- GBM_low_LE_summary %>%
  filter(n_pathways >= 2)

library(clusterProfiler)
library(enrichplot)
library(dplyr)

GBM_low_cnet <- gsea_BBB_GBM_df %>%
  filter(
    Description %in% GBM_low_selected
  )

cnetplot(
  gsea_BBB_GBM,
  showCategory = GBM_low_selected,
  foldChange = geneList_BBB
)

#High BBB stability endothelial cells maintain barrier transport and detoxification programs, whereas low BBB stability cells exhibit metabolic remodeling, junctional disruption, and migratory/inflammatory activation


#



library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)

## 1. Select representative GSEA pathways

GBM_selected_pathways <- c(
  GBM_high_selected,
  GBM_low_selected
)

gsea_plot_df <- gsea_BBB_GBM_df %>%
  filter(
    Description %in% GBM_selected_pathways
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
  gsea_plot_df,
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


## Leading-edge pathway-gene dot plot


# Extract leading-edge genes from selected pathways

GBM_LE_combined <- gsea_BBB_GBM_df %>%
  filter(
    Description %in% GBM_selected_pathways
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

## Remove mitochondrial genes from leading-edge visualization

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


GBM_LE_combined <- GBM_LE_combined %>%
  filter(
    !Gene %in% mt_genes
  )

## Count gene recurrence across pathways

LE_gene_frequency <- GBM_LE_combined %>%
  dplyr::count(
    Gene,
    name = "Pathway_Count"
  )


## Build plotting dataframe

LE_dotplot <- GBM_LE_combined %>%
  group_by(
    Description,
    Gene
  ) %>%
  summarise(
    NES = unique(NES)[1],
    .groups = "drop"
  ) %>%
  left_join(
    LE_gene_frequency,
    by = "Gene"
  )

## Order pathways

# Separate high and low BBB pathways

high_pathways <- gsea_BBB_GBM_df %>%
  filter(
    Description %in% GBM_selected_pathways,
    NES > 0
  ) %>%
  arrange(
    desc(NES)
  ) %>%
  pull(Description)


low_pathways <- gsea_BBB_GBM_df %>%
  filter(
    Description %in% GBM_selected_pathways,
    NES < 0
  ) %>%
  arrange(
    NES
  ) %>%
  pull(Description)


# Put low first in factor levels because ggplot places first levels at bottom
LE_dotplot$Description <- factor(
  LE_dotplot$Description,
  levels = c(
    rev(low_pathways),
    rev(high_pathways)
  )
)


## Order genes left to right

# Order genes by NES contribution:
# High BBB genes first, then Low BBB genes
gene_order <- LE_dotplot %>%
  group_by(Gene) %>%
  summarise(
    meanNES = mean(NES),
    .groups = "drop"
  ) %>%
  arrange(desc(meanNES)) %>%
  pull(Gene)


# Reverse so high BBB genes appear on the left
# and low BBB genes appear on the right
LE_dotplot$Gene <- factor(
  LE_dotplot$Gene,
  levels = (gene_order)
)

## Plot

ggplot(
  LE_dotplot,
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
    xintercept = 47.5,
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
  ggtitle("GBM High vs Low BBB Stability Pathways & Genes")


