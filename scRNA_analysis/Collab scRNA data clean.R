#scRNA seq from collaborators

all_Harmony <-readRDS("Collab scRNA seq Object /all_Harmony.rds")


# Create a new column with the modified values
all_Harmony$sample_type <- ifelse(
  all_Harmony$sample == "Cutaneous_Melanoma", 
  "Melanoma", 
  all_Harmony$sample
)

#Final annnotions based of MarkerGenes
# Vector of new cluster labels (must be named with current cluster IDs as characters)
"0" = "Immune - CD8 T-Cell", #cytotoxic CD8 T-cell
"1" = "Tumor - Prolif.", #melanoma-like tumor cell (prolif sign.)
"2" = "Immune - TAMs", #Myeloid - Activated monocyte/macrophage (TAMs)
"3" = "Immune - TAMs", #TAMs
"4" = "Immune - TAMs" , #Differentiated tissue macrophages / TAM-like macrophages
"5" = "Pericyte", #Pericyte/Vascular Smooth Muscle Cells
"6" = "Tumor - Prolif.", #Prolif MKI67 Tumor cells
"7" = "Neurons", #not ast-like
"8" = "Tumor", #Tumor Melanoma cells?
"9" = "Immune - TRegs", #T cells → CD4+ Regulatory T cells (Tregs)
"10" = "Oligodendrocytes", 
"11" = "Tumor - Epithelial", #Tumor (epithelial carcinoma)
"12" = "Immune - B-lineage", #B-lineage immune cells
"13" = "Oligodendrocytes", 
"14" = "Endothelial", 
"15" = "Astrocytes",
"16" = "Immune - CD8 T-Cells", #Clycling CD8 T Cells
"17" = "Immune - Neutrophils", #Activated inflammatory myeloid cells / monocytes / neutrophils  
"18" = "Immune - Exhausted T-Cells" #Activated/exhausted T cell subset



broad.annotations <- c(
  "0" = "Immune",
  "1" = "Tumor",
  "2" = "Immune",
  "3" = "Immune",
  "4" = "Immune",
  "5" = "Pericyte",
  "6" = "Tumor",
  "7" = "Neural",
  "8" = "Tumor",
  "9" = "Immune",
  "10" = "Oligodendrocytes",
  "11" = "Tumor",
  "12" = "Immune",
  "13" = "Oligodendrocytes",
  "14" = "Endothelial",
  "15" = "Astrocytes",
  "16" = "Immune",
  "17" = "Immune",
  "18" = "Immune"
)

Idents(all_Harmony) <- "seurat_clusters"
all_Harmony <- RenameIdents(all_Harmony, broad.annotations)
all_Harmony$broad_categories <- Idents(all_Harmony)

DimPlot(all_Harmony, group.by = "broad_categories", label = T)


fine.annotations <- c(
  #IMMUNE CELLS
  # T cells
  "0" = "CD8+ Cytotoxic T cells", #KILL INFECTED/TUMOR CELLS
  "16" = "CD8+ Cycling T cells", #PROLIFERATING CD8+ T CELS
  "18" = "CD8+ Exhausted T cells", #DYSFUNCTIONAL CD8+ T CELLS
  "9" = "CD4+ Regulatory T cells", #IMMUNE SUPPRESSION, TOLERANCE
  
  # Myeloid
  "2" = "TAMs - Activated", #ACTIVATED MACROPHAGES IN TUMOR ENV.
  "3" = "TAMs - Inflammatory", #PRO-INFLAMM. MACROPHAGES
  "4" = "TAMs - M2-like", #IMMUNOSUPPRESSIVE, TISSUE REPAIR MACROPHAGES
  "17" = "Inflammatory Monocytes / Neutrophil-like", #INFLAMM. MONOCYTES AND NEUTROPHILE-LIKE CELLS
  
  # B cells
  "12" = "B cells", #ANTIBODY-PRODUCING ADAPTIVE IMMUNE CELLS
  
  # Tumor
  "1" = "Tumor - Proliferative",
  "6" = "Tumor - Highly Proliferative",
  "8" = "Tumor - Invasive / Dedifferentiated",
  "11" = "Tumor - Epithelial",
  
  # Neural
  "7" = "Neurons",
  "10" = "Oligodendrocytes - Mature",
  "13" = "Oligodendrocytes - Immature",
  "15" = "Astrocytes",
  
  # Stromal
  "5" = "Pericytes/vSMC",
  "14" = "Endothelial"
)


Idents(all_Harmony) <- "seurat_clusters"
all_Harmony <- RenameIdents(all_Harmony, fine.annotations)
all_Harmony$fine_annotations <- Idents(all_Harmony)

DimPlot(all_Harmony, group.by = "fine_annotations", label = T)



Idents(all_Harmony) <- "broad_categories"

# Define your clusters
clusters <- c("Neural", "Tumor", "Oligodendrocytes", "Endothelial", "Pericyte", "Immune", "Astrocytes")

# Create a custom color vector
cluster_colors <- c(
  "Neural" = "#1f78b4",    # blue
  "Tumor" = "gray",       # red
  "Oligodendrocytes" = "orange", # green
  "Endothelial" = "red",       
  "Pericyte" = "#B8860B",
  "Immune" = "#6a3d9a",      # purple
  "Astrocytes" = "darkgreen"
)

# Make sure your identities are set correctly
Idents(all_Harmony) <- factor(Idents(all_Harmony), levels = clusters)

theme_santi.umap<-function(base_size = 80, face ="bold", base_family= "Arial") {
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

# Specify the order you want in the legend
desired_order <- c("Tumor", "Neural", "Oligodendrocytes", "Immune", "Astrocytes", "Pericyte", "Endothelial")

# Reorder the factor levels
Idents(all_Harmony) <- factor(Idents(all_Harmony), levels = rev(desired_order))
# Plot UMAP with custom colors
DimPlot(all_Harmony, reduction = "umap", cols = cluster_colors, label = F) +
  ggtitle("Integrated Big Tissue Mets UMAP") +
  theme_santi.umap() +
  guides(color = guide_legend(override.aes = list(size = 16)))


DimPlot(all_Harmony, reduction = "umap", split.by = "sample_type", cols = cluster_colors, label = F) +
  ggtitle("Integrated Big Tissue Mets UMAP") +
  theme_santi.umap() +
  guides(color = guide_legend(override.aes = list(size = 16)))


#DotPlot Marker Genes
Idents(all_Harmony) <- "seurat_clusters"

cluster.markers <- list(
  
  # 0 — Cytotoxic CD8 T cells
  "0" = c("CD8A","CD8B","GZMB","PRF1","NKG7"),
  
  # 1 — Melanoma tumor (MITF-like)
  "1" = c("MITF","SOX10","PMEL","MLANA","TYRP1"),
  
  # 2 — Activated TAMs
  "2" = c("CD68","C1QA","C1QB","TYROBP","APOE"),
  
  # 3 — TAMs (inflammatory skew)
  "3" = c("LST1","FCGR3A","CTSS","ITGAM","CD14"),
  
  # 4 — Differentiated / M2-like macrophages
  "4" = c("CD163","MRC1","TREM2","MSR1","C1QC"),
  
  # 5 — Pericyte / VSMC
  "5" = c("PDGFRB","RGS5","CSPG4","ACTA2","TAGLN"),
  
  # 6 — Proliferating tumor
  "6" = c("MKI67","TOP2A","UBE2C","BIRC5","CDK1"),
  
  # 7 — Neuron-like
  "7" = c("RBFOX3","MAP2","SNAP25","SYT1","SLC17A7"),
  
  # 8 — Tumor (invasive / dedifferentiated melanoma)
  "8" = c("AXL","NGFR","ITGA5","ZEB1","FN1"),
  
  # 9 — Tregs
  "9" = c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2"),
  
  # 10 — Oligodendrocytes
  "10" = c("MBP","PLP1","MOG","MOBP","MAG"),
  
  # 11 — Epithelial carcinoma tumor
  "11" = c("EPCAM","KRT8","KRT18","KRT19","MUC1"),
  
  # 12 — B cells
  "12" = c("MS4A1","CD79A","CD79B","BANK1","CD74"),
  
  # 13 — Oligodendrocytes (OPC-like)
  "13" = c("PDGFRA","CSPG4","OLIG1","OLIG2","SOX10"),
  
  # 14 — Endothelial
  "14" = c("PECAM1","VWF","CLDN5","KDR","FLT1"),
  
  # 15 — Astrocytes
  "15" = c("GFAP","AQP4","ALDH1L1","SLC1A2","SLC1A3"),
  
  # 16 — Cycling CD8 T cells
  "16" = c("CD8A","MKI67","TOP2A","STMN1","PCNA"),
  
  # 17 — Inflammatory monocytes / neutrophil-like
  "17" = c("S100A8","S100A9","FCN1","LYZ","CXCL8"),
  
  # 18 — Exhausted / activated T cells
  "18" = c("PDCD1","LAG3","HAVCR2","TOX","CXCL13")
)


dotplot.genes <- unique(unlist(cluster.markers))


DotPlot(
  all_Harmony,
  features = dotplot.genes,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  scale_color_gradientn(colors = c("lightgrey", "#3b82f6", "#b91c1c")) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90")
  )



library(Seurat)
library(dplyr)

all_Harmony@meta.data <- all_Harmony@meta.data %>%
  select(
    orig.ident,
    nCount_RNA,
    nFeature_RNA,
    sample_id,
    sample,
    percent.mt,
    doublet_call,
    seurat_clusters,
    broad_categories,
    sample_type
  ) %>%
  rename(cell_type = broad_categories)

# Verify the remaining metadata columns
colnames(all_Harmony@meta.data)


write_rds(all_Harmony, "Collab scRNA seq Object /scRNAseq_annotated_obj.rds")


Idents(all_Harmony) <- "broad_categories"

# Define your clusters
clusters <- c("Neural", "Tumor", "Oligodendrocytes", "Endothelial", "Pericyte", "Immune", "Astrocytes")

# Create a custom color vector
cluster_colors <- c(
  "Neural" = "#1f78b4",    # blue
  "Tumor" = "gray",       # red
  "Oligodendrocytes" = "orange", # green
  "Endothelial" = "red",       
  "Pericyte" = "#B8860B",
  "Immune" = "#6a3d9a",      # purple
  "Astrocytes" = "darkgreen"
)

# Make sure your identities are set correctly
Idents(all_Harmony) <- factor(Idents(all_Harmony), levels = clusters)

theme_santi.umap<-function(base_size = 80, face ="bold", base_family= "Arial") {
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

# Specify the order you want in the legend
desired_order <- c("Tumor", "Neural", "Oligodendrocytes", "Immune", "Astrocytes", "Pericyte", "Endothelial")

# Reorder the factor levels
Idents(all_Harmony) <- factor(Idents(all_Harmony), levels = rev(desired_order))
# Plot UMAP with custom colors
DimPlot(all_Harmony, reduction = "umap", cols = cluster_colors, label = F) +
  ggtitle("Integrated Big Tissue Mets UMAP") +
  theme_santi.umap() +
  guides(color = guide_legend(override.aes = list(size = 16)))


DimPlot(all_Harmony, reduction = "umap", split.by = "sample_type", cols = cluster_colors, label = F) +
  ggtitle("Integrated Big Tissue Mets UMAP") +
  theme_santi.umap() +
  guides(color = guide_legend(override.aes = list(size = 16)))


#DotPlot Marker Genes
Idents(all_Harmony) <- "seurat_clusters"

cluster.markers <- list(
  
  # 0 — Cytotoxic CD8 T cells
  "0" = c("CD8A","CD8B","GZMB","PRF1","NKG7"),
  
  # 1 — Melanoma tumor (MITF-like)
  "1" = c("MITF","SOX10","PMEL","MLANA","TYRP1"),
  
  # 2 — Activated TAMs
  "2" = c("CD68","C1QA","C1QB","TYROBP","APOE"),
  
  # 3 — TAMs (inflammatory skew)
  "3" = c("LST1","FCGR3A","CTSS","ITGAM","CD14"),
  
  # 4 — Differentiated / M2-like macrophages
  "4" = c("CD163","MRC1","TREM2","MSR1","C1QC"),
  
  # 5 — Pericyte / VSMC
  "5" = c("PDGFRB","RGS5","CSPG4","ACTA2","TAGLN"),
  
  # 6 — Proliferating tumor
  "6" = c("MKI67","TOP2A","UBE2C","BIRC5","CDK1"),
  
  # 7 — Neuron-like
  "7" = c("RBFOX3","MAP2","SNAP25","SYT1","SLC17A7"),
  
  # 8 — Tumor (invasive / dedifferentiated melanoma)
  "8" = c("AXL","NGFR","ITGA5","ZEB1","FN1"),
  
  # 9 — Tregs
  "9" = c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2"),
  
  # 10 — Oligodendrocytes
  "10" = c("MBP","PLP1","MOG","MOBP","MAG"),
  
  # 11 — Epithelial carcinoma tumor
  "11" = c("EPCAM","KRT8","KRT18","KRT19","MUC1"),
  
  # 12 — B cells
  "12" = c("MS4A1","CD79A","CD79B","BANK1","CD74"),
  
  # 13 — Oligodendrocytes (OPC-like)
  "13" = c("PDGFRA","CSPG4","OLIG1","OLIG2","SOX10"),
  
  # 14 — Endothelial
  "14" = c("PECAM1","VWF","CLDN5","KDR","FLT1"),
  
  # 15 — Astrocytes
  "15" = c("GFAP","AQP4","ALDH1L1","SLC1A2","SLC1A3"),
  
  # 16 — Cycling CD8 T cells
  "16" = c("CD8A","MKI67","TOP2A","STMN1","PCNA"),
  
  # 17 — Inflammatory monocytes / neutrophil-like
  "17" = c("S100A8","S100A9","FCN1","LYZ","CXCL8"),
  
  # 18 — Exhausted / activated T cells
  "18" = c("PDCD1","LAG3","HAVCR2","TOX","CXCL13")
)


dotplot.genes <- unique(unlist(cluster.markers))


DotPlot(
  all_Harmony,
  features = dotplot.genes,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  scale_color_gradientn(colors = c("lightgrey", "#3b82f6", "#b91c1c")) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90")
  )



#DESeq2 for Tumor cells - good cell numbers----

#Tumor Cells
Tumor_obj <- subset(all_Harmony, idents = c("Tumor"))

DimPlot(Tumor_obj, group.by = "sample_type")

Tumor_obj_Br_Mel <- subset(Tumor_obj, subset = sample_type == c("Breast", "Melanoma"))


theme_santi.umap<-function(base_size = 80, face ="bold", base_family= "Arial") {
  theme_few() %+replace%
    theme(axis.ticks = element_line(),
          strip.placement = "outside",
          panel.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(color = "black", size = 0.5),
          plot.title = element_blank(),
          axis.text = element_text(size = 28),
          axis.title = element_text(size = 34),
          legend.text = element_text(size = 22),
          legend.key.size = unit(4, "line")
    )
}

sample_colors <- c(
  "Melanoma" = "#264653",  # deep teal
  "Breast"   = "#F4A261"   # soft orange
)


TumorUMAP <- DimPlot(
  Tumor_obj_Br_Mel,
  group.by = "sample_type",
  cols = sample_colors,
  label = FALSE
) +
  ggtitle(" ") + #scRNA Dataset Tumor Cells
  theme_santi.umap() +
  guides(color = guide_legend(override.aes = list(size = 12))) +  # 12 is huge
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.title = element_blank()
  )


#outdir <- "All Big Tissue /New Filtering/ UMAPs"
outdir <- "Manuscript Figures/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Save as SVG
ggsave(
  filename = file.path(outdir, "scRNA_Dataset_Tumor_Cells_UMAP.svg"),
  plot = TumorUMAP,
  width = 6,
  height = 3.5,
  units = "in"
)








#DESeq2 for tumor cells
# 1. Subset vascular cells
Tumor_obj

# 2. Get raw counts
counts <- GetAssayData(Tumor_obj, assay = "RNA", layer = "counts")

# 3. Prepare metadata
meta <- Tumor_obj@meta.data

# Create one row per sample
coldata <- meta %>%
  dplyr::select(sample_id, sample_type) %>%
  distinct() %>%
  as.data.frame(stringsAsFactors = FALSE)

rownames(coldata) <- coldata$sample_id
coldata$sample_id <- NULL

# 4. Aggregate counts per sample (pseudobulk)
samples <- rownames(coldata)

pseudobulk_counts <- matrix(
  0,
  nrow = nrow(counts),
  ncol = length(samples)
)
rownames(pseudobulk_counts) <- rownames(counts)
colnames(pseudobulk_counts) <- samples

for (s in samples) {
  cells <- colnames(Tumor_obj)[Tumor_obj$sample_id == s]
  pseudobulk_counts[, s] <- rowSums(counts[, cells, drop = FALSE])
}

# Check order matches
all(colnames(pseudobulk_counts) == rownames(coldata))

# 5. Create DESeq2 object
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts,
  colData = coldata,
  design = ~ sample_type
)

# 6. Filter lowly expressed genes (optional but recommended)
keep <- rowSums(counts(dds) >= 10) >= 2  # gene has ≥10 counts in ≥2 samples
dds <- dds[keep, ]

# 7. Run DESeq2
dds <- DESeq(dds)

# 8. Extract results
#Output folder
dir.create("DESeq2_Tumor_Pseudobulk", showWarnings = FALSE)

#Define comparisons
comparisons <- list(
  c("Breast", "Melanoma"),
  c("Breast", "Thyroid"),
  c("Thyroid", "Melanoma")
)


for (comp in comparisons) {
  
  group1 <- comp[1]
  group2 <- comp[2]
  
  # Run DESeq2 results
  res <- results(dds, contrast = c("sample_type", group1, group2))
  
  # Order by adjusted p-value
  res <- res[order(res$padj), ]
  
  # Convert to dataframe
  res_df <- as.data.frame(res)
  
  # Create filename
  filename <- paste0(
    "DESeq2_Tumor_Pseudobulk/",
    "Tumor_pseudobulk_",
    group1,
    "_vs_",
    group2,
    "_DEGs.csv"
  )
  
  # Write file
  write.csv(res_df, file = filename)
}



#Check DEG overlap between scRNA & Xenium

#Br vs. Mel
#Read files
scRNA_Br_Mel_deg <- read.csv("DESeq2_Tumor_Pseudobulk/Tumor_pseudobulk_Breast_vs_Melanoma_DEGs.csv")
scRNA_Br_Mel_deg$gene <- (scRNA_Br_Mel_deg$X)

xenium_Br_Mel_deg <- read.csv("Xenium Tumor Cells DEGs/Tumor_Breast_vs_Melanoma_NVU_sig.csv")
xenium_Br_Mel_deg$gene <- (xenium_Br_Mel_deg$gene)

#Merge and Plot
Br_Mel_merged_deg <- inner_join(scRNA_Br_Mel_deg, xenium_Br_Mel_deg, by = "gene")

colnames(Br_Mel_merged_deg)[colnames(Br_Mel_merged_deg) == "log2FoldChange"] <- "scRNA_logFC"
colnames(Br_Mel_merged_deg)[colnames(Br_Mel_merged_deg) == "avg_log2FC"] <- "xenium_logFC"

#highlight top DEGs
highlight_genes <- Br_Mel_merged_deg %>%
  filter(padj < 0.05 & p_val_adj < 0.05) %>%
  arrange(padj) %>%
  head(15)

#Pearson correlation score
# Calculate Pearson correlation if enough points
if(nrow(Br_Mel_merged_deg) >= 2) {
  
  cor_test <- cor.test(
    Br_Mel_merged_deg$xenium_logFC,
    Br_Mel_merged_deg$scRNA_logFC,
    method = "pearson"
  )
  
  r_val <- cor_test$estimate
  r_squared <- round(r_val^2, 2)  # R²
  
  r_label <- paste0("R² = ", r_squared)
  
} else {
  r_label <- "R² = NA"
}


#plot figure
p_fc <- ggplot(highlight_genes, aes(x = xenium_logFC, y = scRNA_logFC)) +
  geom_point(alpha = 0.8, color = "black", size =3) +
  # add gene labels
  geom_text_repel(
    data = highlight_genes,
    aes(label = gene),
    size = 8,
    max.overlaps = 80
  ) +
  geom_smooth(method = "lm", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = r_label,
    hjust = 1.1,
    vjust = 1.5,
    size = 10
  ) +
  labs(
    x = "Xenium log2 Fold Change",
    y = "scRNA log2 Fold Change",
    #title = "Cross-platform DEG Concordance"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 24),
    axis.title = element_text (size = 28)
  )

p_fc

library(VennDiagram)
sc_genes <- (subset(scRNA_Br_Mel_deg, padj < 0.05 & abs(log2FoldChange) > 1))$X  # padj < 0.05 & abs(log2FoldChange) > 0.5
xen_genes <- subset(xenium_Br_Mel_deg, p_val_adj < 0.05)$gene

venn.plot <- venn.diagram(
  x = list(
    scRNA = sc_genes,
    Xenium = xen_genes
  ),
  filename = NULL,
  fill = c("steelblue", "orange"),
  alpha = 0.5,
  cex = 3,
  cat.cex = 3
)

grid::grid.draw(venn.plot)
# Save as SVG 800 x 800








#GO Enrichment----
library(clusterProfiler)
library(org.Hs.eg.db)

#Split data
br_up <- subset(scRNA_Br_Mel_deg, padj < 0.05 & log2FoldChange > 1)
mel_up <- subset(scRNA_Br_Mel_deg, padj < 0.05 & log2FoldChange < -1)

br_gene_df <- bitr(br_up$gene,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db,
                   drop = TRUE)



mel_gene_df <- bitr(mel_up$gene,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db,
                    drop = TRUE)


ego_br <- enrichGO(
  gene = br_gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

ego_mel <- enrichGO(
  gene = mel_gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

p_br <- dotplot(ego_br, showCategory = 5) +
  ggtitle("Breast Met Enriched pathways") +
  theme_bw()

p_mel <- dotplot(ego_mel, showCategory = 5) +
  ggtitle("Melanoma Met Enriched pathways") +
  theme_bw()

library(patchwork)
p_br | p_mel



#CNET Plot

ego_br_sim <- simplify(
  ego_br,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)


selected_terms_br <- c(
  "tight junction assembly",
  "cell-cell junction organization",
  "positive regulation of vascular permeability",
  #"vascular endothelial growth factor receptor-1 signaling pathway",
  #"leukocyte migration involved in inflammatory response",
  #"response to tumor necrosis factor",
  "angiogenesis involved in wound healing",
  "maintenance of blood-brain barrier"
)

ego_br_subset <- ego_br_sim %>%
  filter(Description %in% selected_terms_br)

fc_br <- br_up$log2FoldChange
names(fc_br) <- br_up$gene

# remove duplicates (important)
fc_br <- fc_br[!duplicated(names(fc_br))]

cnet_Br <- cnetplot(
  ego_br_subset,
  foldChange = fc_br,
  node_label = "all",
  color_category = "darkgreen",
  color_edge = "darkgrey",
  size_category = 1.5,
  size_item = 1,
  size_edge = 0.8
) 

cnet_Br <- cnet_Br +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  ggtitle("Selected GO Terms in Breast Metastasis") +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"))

cnet_Br

#Melanoma
ego_mel_sim <- simplify(
  ego_mel,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

selected_terms_mel <- c(
  "regulation of cell junction assembly",
  "ameboidal-type cell migration",
  "actin filament organization",
  "cell-substrate adhesion",
  "neural crest cell migration"
  #"regulation of neuron projection development"
)

ego_mel_subset <- ego_mel_sim %>%
  filter(Description %in% selected_terms_mel)


fc_mel <- mel_up$log2FoldChange
names(fc_mel) <- mel_up$gene

# remove duplicates (important)
fc_mel <- fc_mel[!duplicated(names(fc_mel))]


library(ggraph)
library(igraph)

cnet_Mel <- cnetplot(
  ego_mel_subset,
  foldChange = fc_mel,   # fold change vector for your genes
  node_label = "all",
  color_category = "darkgreen",  # GO terms
  color_edge = "darkgrey",
  size_category = 1.5,            # GO node size
  size_item = 1,                # gene node size
  size_edge = 0.8
)

# Add color scale for fold changes (blue = negative, red = positive)
cnet_Mel <- cnet_Mel +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  ggtitle("Selected GO Terms in Melanoma Metastasis") +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(size = 18, face = "bold"))

cnet_Mel



cnet_Br / cnet_Mel
combined_cnet <- cnet_Br / cnet_Mel + plot_layout(heights = c(1, 1))

ggsave(
  filename = file.path(outdir, "combined_cnet.svg"),
  plot = combined_cnet,
  width = 10,
  height = 12
)


#Methods (summary)
# 
# Tumor cells were subset from integrated single-cell RNA-seq data processed using Seurat. To enable robust differential expression analysis, counts were aggregated at the sample level to generate pseudobulk profiles. Differential gene expression between tumor types (Breast vs. Melanoma) was performed using DESeq2 with a design formula of ~ sample_type, and genes with low counts were filtered prior to analysis.
# 
# To validate transcriptional signatures across platforms, differentially expressed genes (DEGs) were compared with spatial transcriptomics data (Xenium), and concordance was assessed using Pearson correlation and overlap analysis.
# 
# Functional enrichment analysis of upregulated genes in each tumor type was performed using clusterProfiler with Gene Ontology (biological process) annotations. Redundant GO terms were reduced using semantic similarity–based simplification, and selected pathways were visualized using dot plots and gene–concept network plots (cnetplots), incorporating log2 fold-change values.
# 
# Results (summary)
# 
# Pseudobulk differential expression analysis identified distinct transcriptional programs between breast and melanoma tumor cells, with strong cross-platform concordance observed between single-cell and spatial datasets.
# 
# Gene ontology analysis revealed divergent tumor-intrinsic strategies for interacting with the neurovascular unit (NVU) and blood–brain barrier (BBB):
#   
#   Breast tumor cells were enriched for pathways related to tight junction organization, vascular permeability, angiogenesis, and BBB maintenance, including tight junction assembly, positive regulation of vascular permeability, and angiogenesis involved in wound healing. These findings suggest a phenotype characterized by disruption and remodeling of the BBB, potentially facilitating metastatic entry through increased vascular leakiness and inflammatory signaling.
# Melanoma tumor cells, in contrast, showed enrichment for pathways associated with cytoskeletal remodeling, migration, and neural lineage programs, including ameboidal-type cell migration, actin filament organization, cell-substrate adhesion, and neural crest cell migration. Additionally, enrichment of regulation of cell junction assembly suggests active interaction with endothelial barriers. Together, these results indicate a strategy centered on active transmigration and adaptation to the neural microenvironment, rather than overt disruption of BBB integrity.
# 
# Comparative network visualization further highlighted this contrast, with breast tumor cells emphasizing vascular and barrier dysfunction, whereas melanoma tumor cells emphasized invasive and neural mimicry programs.
# 
# 🔥 (Optional one-line takeaway for figure caption)
# 
# Breast and melanoma brain metastases exhibit distinct tumor-intrinsic programs, with breast tumors promoting BBB disruption and vascular remodeling, while melanoma tumors preferentially engage cytoskeletal invasion and neural-like interactions with the NVU.
# 



#GO Enrichment of 3 way comparisons----
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(patchwork)

padj_cutoff <- 0.05
fc_cutoff <- 1


#reading DESeq2 results
scRNA_Br_Mel_deg <- read.csv("DESeq2_Tumor_Pseudobulk/Tumor_pseudobulk_Breast_vs_Melanoma_DEGs.csv")
scRNA_Br_Thy_deg <- read.csv("DESeq2_Tumor_Pseudobulk/Tumor_pseudobulk_Breast_vs_Thyroid_DEGs.csv")
scRNA_Thy_Mel_deg <- read.csv("DESeq2_Tumor_Pseudobulk/Tumor_pseudobulk_Thyroid_vs_Melanoma_DEGs.csv")

#Breast Enriched Genes
br_vs_mel_up <- scRNA_Br_Mel_deg %>%
  filter(padj < padj_cutoff & log2FoldChange > fc_cutoff) %>%
  pull(X)

br_vs_thy_up <- scRNA_Br_Thy_deg %>%
  filter(padj < padj_cutoff & log2FoldChange > fc_cutoff) %>%
  pull(X)

breast_genes <- union(br_vs_mel_up, br_vs_thy_up) #includes both lists, less strict than intersect


#Melanoma Enriched Genes
mel_vs_br_up <- scRNA_Br_Mel_deg %>%
  filter(padj < padj_cutoff & log2FoldChange < -fc_cutoff) %>%
  pull(X)

mel_vs_thy_up <- scRNA_Thy_Mel_deg %>%
  filter(padj < padj_cutoff & log2FoldChange < -fc_cutoff) %>%
  pull(X)

melanoma_genes <- union(mel_vs_br_up, mel_vs_thy_up)


#Thyroid Enriched Genes
thy_vs_br_up <- scRNA_Br_Thy_deg %>%
  filter(padj < padj_cutoff & log2FoldChange < -fc_cutoff) %>%
  pull(X)

thy_vs_mel_up <- scRNA_Thy_Mel_deg %>%
  filter(padj < padj_cutoff & log2FoldChange > fc_cutoff) %>%
  pull(X)

thyroid_genes <- union(thy_vs_br_up, thy_vs_mel_up)


#Convert Genes Symbols
convert_entrez <- function(genes) {
  bitr(genes,
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = org.Hs.eg.db,
       drop = TRUE)$ENTREZID
}

breast_entrez <- convert_entrez(breast_genes)
melanoma_entrez <- convert_entrez(melanoma_genes)
thyroid_entrez <- convert_entrez(thyroid_genes)


gene_clusters <- list(
  Breast = breast_entrez,
  Melanoma = melanoma_entrez,
  Thyroid = thyroid_entrez
)

cc_go <- compareCluster(
  geneCluster = gene_clusters,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

p_compare <- dotplot(cc_go, showCategory = 20) +
  ggtitle("Tumor-Type Specific GO Biological Processes") +
  theme_bw()

p_compare



cc_go_sim <- simplify(
  cc_go,
  cutoff = 0.7,       # similarity threshold
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

dotplot(cc_go_sim, showCategory = 20) +
  ggtitle("Simplified GO BP Enrichment") +
  theme_bw()




cnetplot(cc_go_sim,
         showCategory = 3,
         color_category = "#E5C494",  # GO node color
         color_item = "#B3B3B3",      # gene node color
         color_edge = "grey",
         size_category = 1.5,
         size_item = 1,
         size_edge = 0.5,
         node_label = "all",
         pie = "equal") +
  ggtitle("GO Term Network Across Tumor Types") +
  theme_bw()



emapplot(cc_go_sim, showCategory = 30) +
  ggtitle("GO Term Enrichment Map") +
  theme_bw()


heatplot(cc_go_sim, showCategory = 20) +
  ggtitle("GO Term Enrichment Heatmap") +
  theme_bw()

upsetplot(cc_go_sim)






#cnet per sample

library(clusterProfiler)
library(ggplot2)

# Extract enrichGO results per cluster
breast_genes  <- cc_go_sim@compareClusterResult %>% dplyr::filter(Cluster == "Breast") %>% pull(geneID)
melanoma_genes <- cc_go_sim@compareClusterResult %>% dplyr::filter(Cluster == "Melanoma") %>% pull(geneID)
thyroid_genes <- cc_go_sim@compareClusterResult %>% dplyr::filter(Cluster == "Thyroid") %>% pull(geneID)

# Convert geneIDs (strings separated by "/") into a vector for each cluster
breast_genes <- unique(unlist(strsplit(breast_genes, "/")))
melanoma_genes <- unique(unlist(strsplit(melanoma_genes, "/")))
thyroid_genes <- unique(unlist(strsplit(thyroid_genes, "/")))

# Map SYMBOL -> ENTREZID if needed
library(org.Hs.eg.db)
breast_entrez <- bitr(breast_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
melanoma_entrez <- bitr(melanoma_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
thyroid_entrez <- bitr(thyroid_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID

# Run enrichGO for each cluster individually
breast_go <- enrichGO(gene = breast_entrez,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      readable = TRUE)

melanoma_go <- enrichGO(gene = melanoma_entrez,
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        readable = TRUE)

thyroid_go <- enrichGO(gene = thyroid_entrez,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       readable = TRUE)

# Make cnetplots
# Basic cnetplot for enrichResult (Breast)
cnetplot(breast_go, showCategory = 10) +
  ggtitle("Breast Tumor GO Cnetplot")

# Melanoma
cnetplot(melanoma_go, showCategory = 10) +
  ggtitle("Melanoma Tumor GO Cnetplot")

# Thyroid
cnetplot(thyroid_go, showCategory = 3) +
  ggtitle("Thyroid Tumor GO Cnetplot")





