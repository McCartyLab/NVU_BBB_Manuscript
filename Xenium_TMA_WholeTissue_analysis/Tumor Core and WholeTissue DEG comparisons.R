#Tumor Core and WholeTissue DEG comparisons


#Condition comparison volcano plots

# Input and output directories
input_dir <- "Vascular subtype analysis/TMA_Celltype_FindMarkers_Condition/NVU_Filtered/"  # folder with FindMarkers CSV files
output_dir <- "Vascular subtype analysis/TMA_Celltype_FindMarkers_Condition/NVU_Filtered/Volcano/"

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




#EC/PC cluster----

#Scatterplot of TMA and Big Tissue DEGs

#using volcano plot because TMA has no data

library(ggplot2)
library(dplyr)
library(ggrepel)

Big_Primary_vs_Normal <- read.csv("Vascular subtype analysis/Big_Tissue_Celltype_FindMarkers_Condition/NVU_Filtered/EC_Primary Tumor_vs_Normal_NVU_sig.csv")

# Prepare data
volcano_data <- Big_Primary_vs_Normal %>%
  mutate(
    gene = as.character(gene),
    neg_log10_padj = -log10(p_val_adj),
    sig = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down",
      TRUE ~ "Not Sig"
    )
  )

# Select top genes to label
labels <- volcano_data %>%
  filter(sig != "Not Sig") %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 41)

# Plot
EC_DEGs_Primary_vs_Normal_BigTissue <- ggplot(volcano_data, aes(x = avg_log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 2) +
  
  geom_text_repel(
    data = labels,
    aes(label = gene),
    size = 5,
    max.overlaps = 100
  ) +
  
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "Not Sig" = "gray"
  )) +
  
  theme_bw() +
  xlab("log2 Fold Change") +
  ylab("-log10(adj p-value)") +
  ggtitle("Gliomas vs Non-Tumor: Tissue Chunks") +
  
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  theme(
    legend.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 16, face = "bold")
  )

EC_DEGs_Primary_vs_Normal_BigTissue        


# Set directory where you want to save
output_dir <- "Vascular subtype analysis/Figures/"
# Make sure the directory exists
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# Export the plot as SVG
ggsave(
  filename = file.path(output_dir, "Glioma_vs_NonTumor_DEGs: Tumor Chunk.svg"),
  plot = EC_DEGs_Primary_vs_Normal_BigTissue,
  width = 8,      # width in inches
  height = 6,     # height in inches
  dpi = 300
)



##Secondary vs Normal----

# 1️⃣ Read data
Big_Secondary_vs_Normal <- read.csv("Vascular subtype analysis/Big_Tissue_Celltype_FindMarkers_Condition/NVU_Filtered/EC_Secondary Tumor_vs_Normal_NVU_sig.csv")
TMA_Secondary_vs_Normal <- read.csv("Vascular subtype analysis/TMAs Condition Pairwise Comparisons_noLuReg/Custom NVU Filtered DEGs/EC_Secondary_vs_Normal_NVU_sig.csv")

# 2️⃣ Prepare TMA DEGs
TMA_deg_Secondary <- TMA_Secondary_vs_Normal %>%
  select(gene, log2FC_TMA = log2FoldChange, padj_TMA = padj) %>%
  mutate(sig_TMA = padj_TMA < 0.05)

# 3️⃣ Prepare Big Tissue DEGs
Big_deg_Secondary <- Big_Secondary_vs_Normal %>%
  select(gene, log2FC_Big = avg_log2FC, padj_Big = p_val_adj) %>%
  mutate(sig_Big = padj_Big < 0.05)

# 4️⃣ Full join to keep all genes
deg_merged_Secondary <- full_join(TMA_deg_Secondary, Big_deg_Secondary, by = "gene")

# 5️⃣ Define coloring: overlapping genes in red, others black
deg_merged_Secondary <- deg_merged_Secondary %>%
  mutate(
    gene_color = ifelse(!is.na(log2FC_TMA) & !is.na(log2FC_Big), "Overlap", "Other"),
    gene_color = factor(gene_color, levels = c("Overlap","Other"))
  )

# 6️⃣ Replace NA with 0 for plotting only
deg_plot_Secondary <- deg_merged_Secondary %>%
  mutate(
    log2FC_TMA_plot = ifelse(is.na(log2FC_TMA), 0, log2FC_TMA),
    log2FC_Big_plot = ifelse(is.na(log2FC_Big), 0, log2FC_Big)
  )

# 7️⃣ Identify labels

# Overlapping genes
label_overlap_Secondary <- deg_plot_Secondary %>% filter(gene_color == "Overlap")

# Top and bottom 5 DEGs for TMA (non-overlapping)
label_TMA_Secondary <- deg_plot_Secondary %>%
  filter(gene_color == "Other" & !is.na(log2FC_TMA)) %>%
  arrange(desc(log2FC_TMA)) %>%
  slice_head(n = 5) %>%
  bind_rows(
    deg_plot_Secondary %>%
      filter(gene_color == "Other" & !is.na(log2FC_TMA)) %>%
      arrange(log2FC_TMA) %>%
      slice_head(n = 5)
  ) %>%
  distinct(gene, .keep_all = TRUE)

# Top and bottom 5 DEGs for Big Tissue (non-overlapping)
label_Big_Secondary <- deg_plot_Secondary %>%
  filter(gene_color == "Other" & !is.na(log2FC_Big)) %>%
  arrange(desc(log2FC_Big)) %>%
  slice_head(n = 5) %>%
  bind_rows(
    deg_plot_Secondary %>%
      filter(gene_color == "Other" & !is.na(log2FC_Big)) %>%
      arrange(log2FC_Big) %>%
      slice_head(n = 5)
  ) %>%
  distinct(gene, .keep_all = TRUE)

# Combine all labels and remove duplicates
labels_Secondary <- bind_rows(label_overlap_Secondary, label_TMA_Secondary, label_Big_Secondary) %>%
  distinct(gene, .keep_all = TRUE)

# 7️⃣ Calculate Pearson correlation for **overlapping gene**
overlap_deg <- deg_plot_Secondary %>% filter(gene_color == "Overlap")

cor_test <- cor.test(overlap_deg$log2FC_TMA, overlap_deg$log2FC_Big, method = "pearson")
r_val <- cor_test$estimate
r_squared <- round(r_val^2, 2)  # R²

# 8️⃣ Scatterplot with jitter, labels, and R²
EC_DEGs_Secondary_vs_Normal <- ggplot(deg_plot_Secondary, aes(x = log2FC_TMA_plot, y = log2FC_Big_plot)) +
  geom_jitter(aes(color = gene_color), width = 0.1, height = 0.1, alpha = 0.6, size = 4) +
  geom_smooth(data = overlap_deg, method = "lm", color = "black", se = FALSE) +
  geom_text_repel(data = labels_Secondary,
                  aes(label = gene),
                  size = 6,
                  max.overlaps = 100) +
  theme_bw() +
  xlab("TMA log2FC") +
  ylab("Big Tissue log2FC") +
  ggtitle("DEG Comparison: TMA vs Big Tissue (Secondary vs Normal)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c(
    "Overlap" = "red",  # genes present in both datasets
    "Other" = "black"   # genes present in only one dataset
  )) +
  annotate("text",
           x = min(deg_plot_Secondary$log2FC_TMA_plot, na.rm = TRUE),
           y = max(deg_plot_Secondary$log2FC_Big_plot, na.rm = TRUE),
           label = paste0("R² = ", r_squared),
           hjust = 0, vjust = 1.2, size = 5) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),     # axis labels
    axis.text = element_text(size = 16),      # tick labels
    plot.title = element_text(size = 16, face = "bold")
  )



EC_DEGs_Secondary_vs_Normal <- ggplot(deg_plot_Secondary, aes(x = log2FC_TMA_plot, y = log2FC_Big_plot)) +
  geom_jitter(
    data = subset(deg_plot_Secondary, gene_color == "Overlap"),
    aes(color = gene_color),
    width = 0.1, height = 0.1, alpha = 0.6, size = 4
  ) +
  geom_jitter(
    data = subset(deg_plot_Secondary, gene_color == "Other"),
    aes(color = gene_color),
    width = 0.1, height = 0.1, alpha = 0.6, size = 2
  ) +
  geom_smooth(data = overlap_deg, method = "lm", color = "black", se = FALSE) +
  geom_text_repel(data = labels_Secondary,
                  aes(label = gene),
                  size = 6,
                  max.overlaps = 100) +
  theme_bw() +
  xlab("TMA log2FC") +
  ylab("Big Tissue log2FC") +
  ggtitle("BrM vs Non-Tumor DEGs correlation across dataset") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c(
    "Overlap" = "red",  # genes present in both datasets
    "Other" = "black"   # genes present in only one dataset
  )) +
  annotate("text",
           x = min(deg_plot_Secondary$log2FC_TMA_plot, na.rm = TRUE),
           y = max(deg_plot_Secondary$log2FC_Big_plot, na.rm = TRUE),
           label = paste0("R² = ", r_squared),
           hjust = 0, vjust = 1.2, size = 5) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),     # axis labels
    axis.text = element_text(size = 16),      # tick labels
    plot.title = element_text(size = 16, face = "bold")
  )

# Set directory where you want to save
output_dir <- "Vascular subtype analysis/Figures/"
# Make sure the directory exists
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# Export the plot as SVG
ggsave(
  filename = file.path(output_dir, "BrM_vs_Non-tumor_DEGs_correlation.svg"),
  plot = EC_DEGs_Secondary_vs_Normal,
  width = 8,      # width in inches
  height = 6,     # height in inches
  dpi = 300
)



##Secondary vs Primary----

# 1️⃣ Read data
TMA_Secondary_vs_Primary <- read.csv("Vascular subtype analysis/TMAs Condition Pairwise Comparisons_noLuReg/Custom NVU Filtered DEGs/EC_Secondary_vs_Primary_NVU_sig.csv")
Big_Secondary_vs_Primary <- read.csv("Vascular subtype analysis/Big_Tissue_Celltype_FindMarkers_Condition/NVU_Filtered/EC_Secondary Tumor_vs_Primary Tumor_NVU_sig.csv")

# 2️⃣ Prepare TMA DEGs
TMA_deg_Secondary_Primary <- TMA_Secondary_vs_Primary %>%
  select(gene, log2FC_TMA = log2FoldChange, padj_TMA = padj) %>%
  mutate(sig_TMA = padj_TMA < 0.05)

# 3️⃣ Prepare Big Tissue DEGs
Big_deg_Secondary_Primary <- Big_Secondary_vs_Primary %>%
  select(gene, log2FC_Big = avg_log2FC, padj_Big = p_val_adj) %>%
  mutate(sig_Big = padj_Big < 0.05)

# 4️⃣ Full join to keep all genes
deg_merged_Secondary_Primary <- full_join(TMA_deg_Secondary_Primary, Big_deg_Secondary_Primary, by = "gene")

# 5️⃣ Define coloring: overlapping genes in red, others black
deg_merged_Secondary_Primary <- deg_merged_Secondary_Primary %>%
  mutate(
    gene_color = ifelse(!is.na(log2FC_TMA) & !is.na(log2FC_Big), "Overlap", "Other"),
    gene_color = factor(gene_color, levels = c("Overlap","Other"))
  )

# 6️⃣ Replace NA with 0 for plotting only
deg_plot_Secondary_Primary <- deg_merged_Secondary_Primary %>%
  mutate(
    log2FC_TMA_plot = ifelse(is.na(log2FC_TMA), 0, log2FC_TMA),
    log2FC_Big_plot = ifelse(is.na(log2FC_Big), 0, log2FC_Big)
  )

# 7️⃣ Identify labels

# Overlapping genes
label_overlap_Secondary_Primary <- deg_plot_Secondary_Primary %>% filter(gene_color == "Overlap")

# Top and bottom 5 DEGs for TMA (non-overlapping)
label_TMA_Secondary_Primary <- deg_plot_Secondary_Primary %>%
  filter(gene_color == "Other" & !is.na(log2FC_TMA)) %>%
  arrange(desc(log2FC_TMA)) %>%
  slice_head(n = 5) %>%
  bind_rows(
    deg_plot_Secondary_Primary %>%
      filter(gene_color == "Other" & !is.na(log2FC_TMA)) %>%
      arrange(log2FC_TMA) %>%
      slice_head(n = 5)
  ) %>%
  distinct(gene, .keep_all = TRUE)

# Top and bottom 5 DEGs for Big Tissue (non-overlapping)
label_Big_Secondary_Primary <- deg_plot_Secondary_Primary %>%
  filter(gene_color == "Other" & !is.na(log2FC_Big)) %>%
  arrange(desc(log2FC_Big)) %>%
  slice_head(n = 5) %>%
  bind_rows(
    deg_plot_Secondary_Primary %>%
      filter(gene_color == "Other" & !is.na(log2FC_Big)) %>%
      arrange(log2FC_Big) %>%
      slice_head(n = 5)
  ) %>%
  distinct(gene, .keep_all = TRUE)

# Combine all labels and remove duplicates
labels_Secondary_Primary <- bind_rows(label_overlap_Secondary_Primary, label_TMA_Secondary_Primary, label_Big_Secondary_Primary) %>%
  distinct(gene, .keep_all = TRUE)

# 7️⃣ Calculate Pearson correlation for overlapping genes
overlap_deg <- deg_plot_Secondary_Primary %>% filter(gene_color == "Overlap")

cor_test <- cor.test(overlap_deg$log2FC_TMA, overlap_deg$log2FC_Big, method = "pearson")
r_val <- cor_test$estimate
r_squared <- round(r_val^2, 2)  # R²

# 8️⃣ Scatterplot with jitter, labels, and R²
EC_DEGs_Secondary_vs_Primary <- ggplot(deg_plot_Secondary_Primary, aes(x = log2FC_TMA_plot, y = log2FC_Big_plot)) +
  geom_jitter(aes(color = gene_color), width = 0.1, height = 0.1, alpha = 0.6, size = 4) +
  geom_smooth(data = overlap_deg, method = "lm", color = "black", se = FALSE) +
  geom_text_repel(data = labels_Secondary_Primary,
                  aes(label = gene),
                  size = 6,
                  max.overlaps = 100) +
  theme_bw() +
  xlab("TMA log2FC") +
  ylab("Big Tissue log2FC") +
  ggtitle("DEG Comparison: TMA vs Big Tissue (Secondary vs Primary)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c(
    "Overlap" = "red",  # genes present in both datasets
    "Other" = "black"   # genes present in only one dataset
  )) +
  annotate("text",
           x = min(deg_plot_Secondary_Primary$log2FC_TMA_plot, na.rm = TRUE),
           y = max(deg_plot_Secondary_Primary$log2FC_Big_plot, na.rm = TRUE),
           label = paste0("R² = ", r_squared),
           hjust = 0, vjust = 1.2, size = 5) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),     # axis labels
    axis.text = element_text(size = 16),      # tick labels
    plot.title = element_text(size = 16, face = "bold")
  )


EC_DEGs_Secondary_vs_Primary <- ggplot(deg_plot_Secondary_Primary, aes(x = log2FC_TMA_plot, y = log2FC_Big_plot)) +
  geom_jitter(
    data = subset(deg_plot_Secondary_Primary, gene_color == "Overlap"),
    aes(color = gene_color),
    width = 0.1, height = 0.1, alpha = 0.6, size = 4
  ) +
  geom_jitter(
    data = subset(deg_plot_Secondary_Primary, gene_color == "Other"),
    aes(color = gene_color),
    width = 0.1, height = 0.1, alpha = 0.6, size = 2
  )+
  geom_smooth(data = overlap_deg, method = "lm", color = "black", se = FALSE) +
  geom_text_repel(data = labels_Secondary_Primary,
                  aes(label = gene),
                  size = 6,
                  max.overlaps = 100) +
  theme_bw() +
  xlab("TMA log2FC") +
  ylab("Big Tissue log2FC") +
  ggtitle("BrM vs Gliomas DEGs correlation across datasets") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c(
    "Overlap" = "red",  # genes present in both datasets
    "Other" = "black"   # genes present in only one dataset
  )) +
  annotate("text",
           x = min(deg_plot_Secondary_Primary$log2FC_TMA_plot, na.rm = TRUE),
           y = max(deg_plot_Secondary_Primary$log2FC_Big_plot, na.rm = TRUE),
           label = paste0("R² = ", r_squared),
           hjust = 0, vjust = 1.2, size = 5) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),     # axis labels
    axis.text = element_text(size = 16),      # tick labels
    plot.title = element_text(size = 16, face = "bold")
  )


# Set directory where you want to save
output_dir <- "Vascular subtype analysis/Figures/"
# Make sure the directory exists
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
# Export the plot as SVG
ggsave(
  filename = file.path(output_dir, "BrM_vs_Gliomas_DEG_Comparison.svg"),
  plot = EC_DEGs_Secondary_vs_Primary,
  width = 8,      # width in inches
  height = 6,     # height in inches
  dpi = 300
)



