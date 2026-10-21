#Xenium WholeTissue DEG analysis



#Whole Tissue-----

#FindMarkers- Condition Pairwise Comparisons----
obj <- vascular_subtypes_whole_tissue

# Make sure identities are condition
Idents(obj) <- "condition"

# LOAD NVU PANEL

gene.list <- read.csv("Xenium Human Brain Panel Gene List-Santi.csv")

nvu_genes <- gene.list %>%
  filter(Santi.Annotations == "NVU") %>%
  pull(gene) %>%
  unique()

nvu_genes2 <- nvu_genes[!(nvu_genes %in% c("ALDH1L1", "AQP4" ,"GJA1","GFAP", "MLC1", "BDNF", "GABBR1", "NDP", "ATP1B2", "WNT7A", "ALB", "SLC1A2"))]

nvu_genes <- toupper(nvu_genes2)



# CELL TYPES

celltypes <- c("EC", "Pericyte", "VLMC", "LowQual_drop")

# DEFINE ONLY THE 3 COMPARISONS YOU WANT

comparisons <- data.frame(
  group1 = c("Primary Tumor",
             "Secondary Tumor",
             "Secondary Tumor"),
  
  group2 = c("Normal",
             "Normal",
             "Primary Tumor")
)

# OUTPUT DIRECTORIES
outdir   <- "Vascular subtype analysis/Big_Tissue_Celltype_FindMarkers_Condition/"
raw_dir  <- file.path(outdir, "Raw")
filt_dir <- file.path(outdir, "Filtered")
nvu_dir  <- file.path(outdir, "NVU_Filtered")

dir.create(raw_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nvu_dir,  recursive = TRUE, showWarnings = FALSE)

# FILTERING FUNCTIONS
filter_degs <- function(res_df,
                        lfc_cutoff = 1,
                        padj_cutoff = 0.05) {
  
  res_df %>%
    filter(!is.na(p_val_adj),
           p_val_adj < padj_cutoff,
           abs(avg_log2FC) > lfc_cutoff)
}

filter_sig_nvu <- function(res_df,
                           nvu_genes,
                           lfc_cutoff = 1,
                           padj_cutoff = 0.05) {
  
  res_df %>%
    mutate(gene = toupper(gene)) %>%
    filter(gene %in% nvu_genes,
           !is.na(p_val_adj),
           p_val_adj < padj_cutoff,
           abs(avg_log2FC) > lfc_cutoff)
}

# MAIN LOOP
for (ct in celltypes) {
  
  cat("Processing cell type:", ct, "\n")
  
  # Subset to cell type
  ct_obj <- subset(obj, subset = new_clusters == ct)
  
  if (ncol(ct_obj) < 50) {
    message("Skipping ", ct, " (too few cells)")
    next
  }
  
  for (i in 1:nrow(comparisons)) {
    
    g1 <- comparisons$group1[i]
    g2 <- comparisons$group2[i]
    
    cat("   ", g1, "vs", g2, "\n")
    
    deg <- tryCatch(
      FindMarkers(
        ct_obj,
        ident.1 = g1,
        ident.2 = g2,
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0
      ),
      error = function(e) {
        message("      SKIP: ", g1, " vs ", g2)
        NULL
      }
    )
    
    if (!is.null(deg) && nrow(deg) > 0) {
      
      deg <- deg %>%
        rownames_to_column("gene")
      
      # RAW
      write.csv(
        deg,
        file = file.path(
          raw_dir,
          paste0(gsub("/", "_", ct), "_",
                 g1, "_vs_", g2, ".csv")
        ),
        row.names = FALSE
      )
      
      # SIGNIFICANT
      sig <- filter_degs(deg)
      
      write.csv(
        sig,
        file = file.path(
          filt_dir,
          paste0(gsub("/", "_", ct), "_",
                 g1, "_vs_", g2, "_sig.csv")
        ),
        row.names = FALSE
      )
      
      # NVU SIGNIFICANT
      sig_nvu <- filter_sig_nvu(deg, nvu_genes)
      
      write.csv(
        sig_nvu,
        file = file.path(
          nvu_dir,
          paste0(gsub("/", "_", ct), "_",
                 g1, "_vs_", g2, "_NVU_sig.csv")
        ),
        row.names = FALSE
      )
    }
  }
}






#FindMarkers- Sample Pairwise Comparisons----
obj <- vascular_subtypes_whole_tissue

# Make sure identities are tissue groups
# (Change this if needed)

Idents(obj) <- obj$sample

# LOAD NVU PANEL
# gene.list <- read.csv("Xenium Human Brain Panel Gene List-Santi.csv")
# 
# nvu_genes <- gene.list %>%
#   filter(Santi.Annotations == "NVU") %>%
#   pull(gene) %>%      # change column name if needed
#   unique()
# 
# # Optional safety step (avoid case mismatch issues)
# nvu_genes <- toupper(nvu_genes2)

library(conflicted)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("rename", "dplyr")

conflict_prefer("intersect", "base")
conflict_prefer("setdiff", "base")


# CELL TYPES
celltypes <- c("EC", "Pericyte", "VLMC", "LowQual_drop")


# BUILD PAIRWISE COMPARISONS
# Normal always second

groups <- levels(Idents(obj))

# Separate Normal
non_normal <- setdiff(groups, "Normal")

# Tumor vs Tumor
tumor_comparisons <- t(combn(non_normal, 2))
tumor_comparisons <- as.data.frame(tumor_comparisons)
colnames(tumor_comparisons) <- c("group1", "group2")

# Tumor vs Normal (Normal always second)
normal_comparisons <- data.frame(
  group1 = non_normal,
  group2 = "Normal"
)

# Combine
comparisons <- bind_rows(tumor_comparisons, normal_comparisons)


# OUTPUT DIRECTORIES

outdir   <- "Vascular subtype analysis/Big_Tissue_Celltype_FindMarkers/"
raw_dir  <- file.path(outdir, "Raw")
filt_dir <- file.path(outdir, "Filtered")
nvu_dir  <- file.path(outdir, "NVU_Filtered")

dir.create(raw_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nvu_dir,  recursive = TRUE, showWarnings = FALSE)

# FILTERING FUNCTIONS

filter_degs <- function(res_df,
                        lfc_cutoff = 1,
                        padj_cutoff = 0.05) {
  
  res_df %>%
    filter(!is.na(p_val_adj),
           p_val_adj < padj_cutoff,
           abs(avg_log2FC) > lfc_cutoff)
}

filter_sig_nvu <- function(res_df,
                           nvu_genes,
                           lfc_cutoff = 1,
                           padj_cutoff = 0.05) {
  
  res_df %>%
    mutate(gene = toupper(gene)) %>%
    filter(gene %in% nvu_genes,
           !is.na(p_val_adj),
           p_val_adj < padj_cutoff,
           abs(avg_log2FC) > lfc_cutoff)
}

# MAIN LOOP
for (ct in celltypes) {
  
  cat("Processing cell type:", ct, "\n")
  
  # Subset to cell type
  ct_obj <- subset(obj, subset = new_clusters == ct)
  
  if (ncol(ct_obj) < 50) {
    message("Skipping ", ct, " (too few cells)")
    next
  }
  
  for (i in 1:nrow(comparisons)) {
    
    g1 <- comparisons$group1[i]
    g2 <- comparisons$group2[i]
    
    cat("   ", g1, "vs", g2, "\n")
    
    deg <- tryCatch(
      FindMarkers(
        ct_obj,
        ident.1 = g1,
        ident.2 = g2,
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0
      ),
      error = function(e) {
        message("      SKIP: ", g1, " vs ", g2)
        NULL
      }
    )
    
    if (!is.null(deg) && nrow(deg) > 0) {
      
      deg <- deg %>%
        rownames_to_column("gene")
      
      # RAW
      write.csv(
        deg,
        file = file.path(
          raw_dir,
          paste0(gsub("/", "_", ct), "_",
                 g1, "_vs_", g2, ".csv")
        ),
        row.names = FALSE
      )
      
      # SIGNIFICANT
      sig <- filter_degs(deg)
      
      write.csv(
        sig,
        file = file.path(
          filt_dir,
          paste0(gsub("/", "_", ct), "_",
                 g1, "_vs_", g2, "_sig.csv")
        ),
        row.names = FALSE
      )
      
      # NVU SIGNIFICANT
      sig_nvu <- filter_sig_nvu(deg, nvu_genes)
      
      write.csv(
        sig_nvu,
        file = file.path(
          nvu_dir,
          paste0(gsub("/", "_", ct), "_",
                 g1, "_vs_", g2, "_NVU_sig.csv")
        ),
        row.names = FALSE
      )
    }
  }
}






#Adding new column
whole_tissue$cluster_sample <- paste0(whole_tissue$clusters,"_",whole_tissue$sample)
Idents(whole_tissue) <- whole_tissue$cluster_sample

whole_tissue$cluster_sample <-
  as.factor(whole_tissue$cluster_sample)


DotPlot(whole_tissue, features = c("MMP9", "SLC1A2"),
        idents = c("EC/PC_Breast", "EC/PC_Normal"))



#Creating ShinyApp to make heatmaps
#ShinyApp----
#export shinyapp
reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
           "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}


devtools::install_github("SGDDNB/ShinyCell")
library(ShinyCell)

Idents(whole_tissue) <- whole_tissue$clusters
levels(x = whole_tissue) <- c('EC/PC', 'Ast-like', 'Immune', "Fibroblasts", "OPC/OL-like", "Neu-like", "Tumor")


scConf <- createConfig(whole_tissue)
makeShinyApp(whole_tissue, gex.assay = "Xenium", scConf, gene.mapping = TRUE,
             shiny.title = "whole_tissue") 



#Exporting a shinyapp for each cell type
seurat_obj <- whole_tissue
clusters_to_split <- levels(seurat_obj$clusters)

base_outdir <- "ShinyCell_per_cluster_apps"
dir.create(base_outdir, recursive = TRUE, showWarnings = FALSE)

for (cl in clusters_to_split) {
  
  cat("Creating app for cluster:", cl, "\n")
  
  obj_sub <- subset(seurat_obj, subset = clusters == cl)
  
  scConf <- createConfig(obj_sub)
  
  # Build app (creates folder called "shinyApp")
  makeShinyApp(
    obj_sub,
    scConf = scConf,
    gex.assay = "Xenium",
    gex.slot  = "counts",
    shiny.title = paste0("ShinyCell: ", cl)
  )
  
  # Destination folder
  new_folder <- file.path(base_outdir, paste0("ShinyApp_", cl))
  
  # Remove destination if exists
  if (dir.exists(new_folder)) {
    unlink(new_folder, recursive = TRUE)
  }
  
  # Copy shinyApp folder
  dir.create(new_folder, recursive = TRUE)
  file.copy("shinyApp", new_folder, recursive = TRUE)
  
  # Remove original shinyApp so next iteration starts clean
  unlink("shinyApp", recursive = TRUE)
  
  cat("Finished:", cl, "\n\n")
}








#Condition comparison volcano plots

# Input and output directories
input_dir <- "Vascular subtype analysis/Big_Tissue_Celltype_FindMarkers_Condition/NVU_Filtered/"  # folder with FindMarkers CSV files
output_dir <- "Vascular subtype analysis/Big_Tissue_Celltype_FindMarkers_Condition/NVU_Filtered/Volcano/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Volcano plot function
plot_volcano <- function(df, cell_type, comparison,
                         lfc_cutoff = 0.25,  # typical Seurat logFC threshold
                         padj_cutoff = 0.05,
                         top_n = 10) {
  
  df <- df %>%
    filter(!is.na(p_val_adj)) %>%
    mutate(
      negLog10Padj = -log10(p_val_adj),
      Significance = case_when(
        p_val_adj < padj_cutoff & avg_log2FC > lfc_cutoff  ~ "Up",
        p_val_adj < padj_cutoff & avg_log2FC < -lfc_cutoff ~ "Down",
        TRUE ~ "Not Sig"
      )
    )
  
  # Top genes by smallest adjusted p-value
  top_genes <- df %>%
    arrange(p_val_adj) %>%
    head(50)
  
  max_lfc <- max(abs(df$avg_log2FC), na.rm = TRUE)
  
  ggplot(df, aes(x = avg_log2FC, y = negLog10Padj, color = Significance)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c(
      "Up" = "#D62728",
      "Down" = "#1F77B4",
      "Not Sig" = "grey70"
    )) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    geom_text_repel(data = top_genes, aes(label = gene), size = 5, max.overlaps = 15) +
    coord_cartesian(xlim = c(-max_lfc, max_lfc)) +
    labs(
      title = paste(cell_type, "–", comparison),
      x = "Average Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    )
}

# Process all CSV files
all_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

for (file in all_files) {
  
  df <- read.csv(file)
  
  # File name without extension
  filename <- basename(file)
  filename_no_ext <- sub("\\.csv$", "", filename)
  
  # Split at first underscore for cell type
  split_name <- strsplit(filename_no_ext, "_", fixed = TRUE)[[1]]
  cell_type <- split_name[1]
  
  # Reconstruct comparison name (everything after first underscore)
  comparison <- paste(split_name[-1], collapse = "_")
  comparison <- gsub("_vs_", " vs ", comparison)
  
  # Create output folder per cell type
  cell_output_dir <- file.path(output_dir, cell_type)
  if (!dir.exists(cell_output_dir)) dir.create(cell_output_dir, recursive = TRUE)
  
  # Generate volcano plot
  p <- plot_volcano(df, cell_type, comparison)
  
  # Save plot as SVG
  ggsave(
    filename = file.path(cell_output_dir,
                         paste0("Volcano_", filename_no_ext, ".svg")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}







#All Primary and Secondary vs. Normal comparisons----
## Top 10 Up + Top 10 Down NVU genes per comparison

# 1️⃣  Folder containing all comparisons
deg_dir <- "Vascular subtype analysis/Big_Tissue_Celltype_FindMarkers/Raw/"
#gene_panel_file <- "Xenium Human Brain Panel Gene List-Santi.csv"

# 1️⃣ List all comparison files
deg_files <- list.files(deg_dir, pattern = "EC_.*\\.csv$", full.names = TRUE)

# Keep only comparisons vs Normal
deg_files_vs_normal <- deg_files[grepl("_vs_Normal\\.csv$", basename(deg_files))]

# 2️⃣ Function to read and tidy each file
read_deg <- function(file) {
  comp_name <- str_remove(basename(file), ".csv")
  
  read.csv(file) %>%
    select(gene, avg_log2FC, p_val_adj) %>%
    rename(log2FC = avg_log2FC, padj = p_val_adj) %>%
    mutate(comparison = comp_name)
}

# 3️⃣ Read and combine all comparisons
deg_all <- map_df(deg_files_vs_normal, read_deg)

# 4️⃣ Load NVU gene panel and keep only NVU genes
# gene.list <- read.csv(gene_panel_file)
# nvu_genes <- gene.list %>%
#   filter(Santi.Annotations == "NVU") %>%
#   pull(gene) %>%
#   unique() %>%
#   toupper()


deg_filtered <- deg_all %>%
  filter(gene %in% nvu_genes)

# 5️⃣ Select top 10 up and bottom 10 down per comparison
deg_top20 <- deg_filtered %>%
  filter(padj < 0.05) %>%
  group_by(comparison) %>%
  arrange(desc(log2FC)) %>%
  slice_head(n = 10) %>%  # top 10 up
  bind_rows(
    deg_filtered %>%
      filter(padj < 0.05) %>%
      group_by(comparison) %>%
      arrange(log2FC) %>%
      slice_head(n = 10)    # top 10 down
  ) %>%
  ungroup() %>%
  distinct(gene, comparison, .keep_all = TRUE)

genes_to_plot <- unique(deg_top20$gene)

# 6️⃣ Pivot for heatmap
heatmap_df <- deg_filtered %>%
  filter(gene %in% genes_to_plot) %>%
  select(gene, comparison, log2FC) %>%
  distinct() %>%
  pivot_wider(names_from = comparison, values_from = log2FC, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_df[,-1])
rownames(heatmap_matrix) <- heatmap_df$gene

# 7️⃣ Scale rows (genes)
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

# 8️⃣ Clean x-axis names
colnames(heatmap_matrix_scaled) <- str_remove(colnames(heatmap_matrix_scaled), "EC_")
colnames(heatmap_matrix_scaled) <- str_remove(colnames(heatmap_matrix_scaled), "_vs_Normal")

conflicts_prefer(ComplexHeatmap::pheatmap)
conflicts_prefer(pheatmap::pheatmap)



# 9️⃣ Flip heatmap orientation
# Genes = columns (x-axis)
# Comparisons = rows (y-axis)

colnames(heatmap_matrix_scaled) <- colnames(heatmap_matrix_scaled) %>%
  str_replace("^EC_Gr_3$", "Grade 3") %>%
  str_remove("^EC_")

heatmap_matrix_scaled_t <- t(heatmap_matrix_scaled)

# 10️⃣ Plot heatmap
EC_summary_heatmap <- pheatmap(
  heatmap_matrix_scaled_t,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 16,
  fontsize_col = 12,
  angle_col = 45,
  main = "NVU DEGs per Tumor vs Non-Tumor"
)

#Differential Expression Heatmap Analysis: Differentially expressed genes (DEGs) for endothelial/pericyte (EC/PC) cell types were obtained for each tumor subtype compared to normal tissue using the FindMarkers function. Genes present in a custom NVU (neurovascular unit) gene panel were excluded. For visualization, the top 5 upregulated and top 5 downregulated DEGs per comparison (adjusted p-value < 0.05) were selected. A log2 fold-change matrix was constructed with genes as rows and tumor comparisons as columns. Heatmaps were generated using the pheatmap R package, with row-scaling to highlight relative expression patterns across comparisons.

# Set directory where you want to save
output_dir <- "Vascular subtype analysis/Figures/"
# Make sure the directory exists
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# Export the plot as SVG
ggsave(
  filename = file.path(output_dir, "Chunks_Tumor_vs_NonTumor_DEGs_per_sample.svg"),
  plot = EC_summary_heatmap,
  width = 12,      # width in inches
  height = 6,     # height in inches
  dpi = 300
)

