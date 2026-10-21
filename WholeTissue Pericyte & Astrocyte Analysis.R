#Pericyte & Astrocyte DEG analysis

#Astrocyte BBB Support Analysis

# Pericyte DEG analysis across NVU states----

library(Seurat)
library(dplyr)

# Get pericyte cells
pericyte_cells <- WhichCells(
  whole_tissue,
  expression = new_clusters == "Pericyte"
)

pericyte_obj <- subset(
  whole_tissue,
  cells = pericyte_cells
)

# Set NVU state identity
Idents(pericyte_obj) <- "NVU_States"

# Check available pericytes
table(Idents(pericyte_obj), useNA = "ifany")

# N.T. Intact vs T.I. Intact
pericyte_NT_vs_TI <- FindMarkers(
  pericyte_obj,
  ident.1 = "T.I.: Intact",
  ident.2 = "N.T.: Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

# T.I. Intact vs T.I. Astrocyte displaced
pericyte_TI_vs_astro <- FindMarkers(
  pericyte_obj,
  ident.1 = "T.I.: Astrocyte displaced",
  ident.2 = "T.I.: Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

# Add gene names
pericyte_NT_vs_TI$gene <- rownames(pericyte_NT_vs_TI)
pericyte_TI_vs_astro$gene <- rownames(pericyte_TI_vs_astro)

# Order by adjusted p-value
pericyte_NT_vs_TI <- pericyte_NT_vs_TI %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))

pericyte_TI_vs_astro <- pericyte_TI_vs_astro %>%
  arrange(p_val_adj, desc(abs(avg_log2FC)))

# Inspect top genes
head(pericyte_NT_vs_TI, 20)
head(pericyte_TI_vs_astro, 20)


#NVU genes

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

#filter out genes non-NVU genes to avoid spillover
# Pericyte DEGs restricted to curated NVU genes

pericyte_NT_vs_TI_NVU <- pericyte_NT_vs_TI %>%
  filter(gene %in% nvu_genes) %>%
  arrange(desc(avg_log2FC))

pericyte_TI_vs_astro_NVU <- pericyte_TI_vs_astro %>%
  filter(gene %in% nvu_genes) %>%
  arrange(desc(avg_log2FC))

pericyte_NT_vs_TI_NVU
pericyte_TI_vs_astro_NVU

#volcano plots---

# Pericyte NVU Volcano Plots
# Curated NVU gene set
# N.T. Intact vs T.I. Intact
# T.I. Intact vs T.I. Astrocyte Displaced

library(dplyr)
library(ggplot2)
library(ggrepel)

# Prepare volcano data
prepare_volcano <- function(deg_df, gene_set) {
  
  df <- deg_df %>%
    filter(gene %in% gene_set) %>%
    mutate(
      neg_log10_padj = -log10(
        pmax(p_val_adj, .Machine$double.xmin)
      ),
      sig_up = p_val_adj < 0.05 & avg_log2FC >= 1,
      sig_down = p_val_adj < 0.05 & avg_log2FC <= -1,
      sig = sig_up | sig_down
    )
  
  # Top 10 upregulated genes
  top_up <- df %>%
    filter(avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 10)
  
  # Top 10 downregulated genes
  top_down <- df %>%
    filter(avg_log2FC < 0) %>%
    arrange(avg_log2FC) %>%
    slice_head(n = 10)
  
  # Add labels
  df <- df %>%
    mutate(
      label = ifelse(
        gene %in% c(top_up$gene, top_down$gene),
        gene,
        NA_character_
      )
    )
  
  return(df)
}


# Prepare N.T. Intact vs T.I. Intact
volcano_NT_TI <- prepare_volcano(
  pericyte_NT_vs_TI,
  nvu_genes
)


# Prepare T.I. Intact vs T.I. Astrocyte Displaced
volcano_TI_astro <- prepare_volcano(
  pericyte_TI_vs_astro,
  nvu_genes
)


# Volcano plot function
make_volcano <- function(df, title) {
  
  ggplot(df, aes(
    x = avg_log2FC,
    y = neg_log10_padj
  )) +
    
    # Non-significant genes
    geom_point(
      data = df %>% filter(!sig),
      size = 2.5,
      alpha = 0.6,
      color = "grey70"
    ) +
    
    # Significant upregulated genes
    geom_point(
      data = df %>% filter(sig_up),
      size = 3,
      alpha = 0.9,
      color = "red"
    ) +
    
    # Significant downregulated genes
    geom_point(
      data = df %>% filter(sig_down),
      size = 3,
      alpha = 0.9,
      color = "blue"
    ) +
    
    # Top 10 up/down gene labels
    geom_text_repel(
      data = df %>% filter(!is.na(label)),
      aes(label = label),
      size = 4,
      box.padding = 0.4,
      point.padding = 0.2,
      max.overlaps = Inf,
      min.segment.length = 0,
      seed = 123
    ) +
    
    # log2FC reference lines
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dotted",
      linewidth = 0.7
    ) +
    
    # FDR reference line
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      linewidth = 0.7
    ) +
    
    labs(
      title = title,
      x = "Average log2 fold change",
      y = expression(-log[10]("adjusted p-value"))
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
}


# Create volcano plots
volcano_NT_TI_plot <- make_volcano(
  volcano_NT_TI,
  "Pericytes: N.T. Intact vs T.I. Intact"
)

volcano_TI_astro_plot <- make_volcano(
  volcano_TI_astro,
  "Pericytes: T.I. Intact vs Astrocyte Displaced"
)


# Display plots
volcano_NT_TI_plot
volcano_TI_astro_plot


# Save plots
ggsave(
  "Vascular subtype analysis/Figures/Pericyte_NVU_NT_vs_TI_Volcano.svg",
  volcano_NT_TI_plot,
  width = 8,
  height = 7
)

ggsave(
  "Vascular subtype analysis/Figures/Pericyte_NVU_TI_vs_Astro_Volcano.svg",
  volcano_TI_astro_plot,
  width = 8,
  height = 7
)


#Per New_Condition----
# Cell type per NVU State----

library(dplyr)
library(tidyr)
library(ggplot2)

# Define condition order
condition_order <- c(
  "Non-Tumor",
  "Gliomas",
  "BrM"
)

# Define NVU states to show
NVU_states_plot <- c(
  "N.T.: Intact",
  "T.I.: Intact",
  "T.I.: Astrocyte displaced",
  "T.I.: Pericyte displaced",
  "T.I.: Endothelial only"
)

# Define cell types
NVU_celltypes <- c(
  "EC",
  "Pericyte",
  "Ast-like",
  "Tumor"
)


# Add simplified condition label
whole_tissue_NVU_meta <- whole_tissue_NVU_meta %>%
  mutate(
    new_conditons = case_when(
      condition == "Normal" ~ "Non-Tumor",
      condition == "Primary Tumor" ~ "Gliomas",
      condition == "Secondary Tumor" ~ "BrM",
      TRUE ~ NA_character_
    )
  )


# Check condition assignments
table(
  whole_tissue_NVU_meta$new_conditons,
  useNA = "ifany"
)


# Calculate cellular composition within each NVU state
NVU_composition <- whole_tissue_NVU_meta %>%
  filter(
    NVU_States %in% NVU_states_plot,
    new_clusters %in% NVU_celltypes,
    new_conditons %in% condition_order
  ) %>%
  count(
    new_conditons,
    NVU_States,
    new_clusters,
    name = "n"
  ) %>%
  group_by(
    new_conditons,
    NVU_States
  ) %>%
  mutate(
    percentage = 100 * n / sum(n)
  ) %>%
  ungroup()


# Add zero-count combinations so all cell types appear
NVU_composition <- NVU_composition %>%
  complete(
    new_conditons = condition_order,
    NVU_States = NVU_states_plot,
    new_clusters = NVU_celltypes,
    fill = list(
      n = 0,
      percentage = 0
    )
  )


# Set factor ordering
NVU_composition$new_conditons <- factor(
  NVU_composition$new_conditons,
  levels = condition_order
)

NVU_composition$NVU_States <- factor(
  NVU_composition$NVU_States,
  levels = NVU_states_plot
)

NVU_composition$new_clusters <- factor(
  NVU_composition$new_clusters,
  levels = NVU_celltypes,
  labels = c(
    "Endothelial",
    "Pericyte",
    "Astrocyte",
    "Tumor"
  )
)


# Plot cellular composition
NVU_cell_composition_dotplot <- ggplot(
  NVU_composition,
  aes(
    x = NVU_States,
    y = new_clusters,
    size = n,
    color = percentage
  )
) +
  geom_point(alpha = 0.9) +
  
  scale_size_continuous(
    name = "Cell count",
    range = c(2, 12)
  ) +
  
  scale_color_gradient(
    name = "Composition (%)",
    low = "grey90",
    high = "#2166AC"
  ) +
  
  facet_wrap(
    ~ new_conditons,
    ncol = 1
  ) +
  
  labs(
    title = "Cellular Composition Across NVU States",
    x = NULL,
    y = NULL
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    strip.text = element_text(
      size = 16,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 15,
      angle = 30,
      hjust = 1,
    ),
    axis.text.y = element_text(
      size = 15
    ),
    legend.title = element_text(
      size = 14
    ),
    legend.text = element_text(
      size = 12
    )
  )

NVU_cell_composition_dotplot


# Save composition plot
ggsave(
  "Vascular subtype analysis/Figures/NVU_Cell_Composition_DotPlot_by_Condition.svg",
  NVU_cell_composition_dotplot,
  width = 7.2,
  height = 8
)


# Pericyte DEG analysis across NVU states and conditions----

# IMPORTANT:
# Work from metadata/cell IDs rather than creating Seurat subsets.
# This avoids the memory exhaustion problem encountered with subset().


# Get metadata directly from the existing Seurat object
whole_tissue_meta <- whole_tissue[[]]


# Make condition labels
whole_tissue_meta$new_conditons <- dplyr::case_when(
  whole_tissue_meta$condition == "Normal" ~ "Non-Tumor",
  whole_tissue_meta$condition == "Primary Tumor" ~ "Gliomas",
  whole_tissue_meta$condition == "Secondary Tumor" ~ "BrM",
  TRUE ~ NA_character_
)


# Curated NVU gene set
nvu_genes <- unique(c(
  "CLDN5",
  "OCLN",
  "TJP1",
  "TJP2",
  "TJP3",
  "F11R",
  "JAM2",
  "JAM3",
  "LSR",
  "ABCB1",
  "ABCG2",
  "SLC2A1",
  "SLCO1A2",
  "SLCO1C1",
  "LRP1",
  "MFSD2A",
  "ANGPT2",
  "VEGFA",
  "DLL4",
  "NDP",
  "KDR",
  "TEK",
  "ICAM1",
  "VCAM1",
  "IL1B",
  "IL6",
  "CCL2",
  "S1PR3",
  "COL4A1",
  "LAMB1",
  "FN1",
  "VCL",
  "ITGB1",
  "MMP2",
  "MMP9",
  "VIM",
  "PLVAP",
  "ESAM",
  "ENG",
  "CTNNB1",
  "VWF",
  "TFRC",
  "ADGRL4",
  "KCNJ8",
  "SLC1A3",
  "SLC2A3",
  "SLC16A1",
  "SLC7A5",
  "ATP10A",
  "FZD4",
  "TSC2",
  "MCAM",
  "ADGRA2",
  "ABCC5",
  "FABP5"
))


# Check required metadata
stopifnot(
  all(c(
    "new_conditons",
    "NVU_States",
    "new_clusters"
  ) %in% colnames(whole_tissue_meta))
)


# Conditions
conditions_to_analyze <- c(
  "Non-Tumor",
  "Gliomas",
  "BrM"
)


# Initialize result lists
pericyte_NT_vs_TI_by_condition <- list(
  "Non-Tumor" = NULL,
  "Gliomas" = NULL,
  "BrM" = NULL
)

pericyte_TI_vs_astro_by_condition <- list(
  "Non-Tumor" = NULL,
  "Gliomas" = NULL,
  "BrM" = NULL
)


# 1. NON-TUMOR: T.I. Intact vs N.T. Intact

cells_TI_intact_NonTumor <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Non-Tumor" &
    whole_tissue_meta$new_clusters == "Pericyte" &
    whole_tissue_meta$NVU_States == "T.I.: Intact"
]

cells_NT_intact_NonTumor <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Non-Tumor" &
    whole_tissue_meta$new_clusters == "Pericyte" &
    whole_tissue_meta$NVU_States == "N.T.: Intact"
]

message(
  "Non-Tumor: N.T. Intact = ",
  length(cells_NT_intact_NonTumor),
  "; T.I. Intact = ",
  length(cells_TI_intact_NonTumor)
)

old_idents <- Idents(whole_tissue)

Idents(
  whole_tissue,
  cells = cells_TI_intact_NonTumor
) <- "T.I. Intact"

Idents(
  whole_tissue,
  cells = cells_NT_intact_NonTumor
) <- "N.T. Intact"

pericyte_NT_vs_TI_by_condition[["Non-Tumor"]] <- FindMarkers(
  whole_tissue,
  ident.1 = "T.I. Intact",
  ident.2 = "N.T. Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

Idents(whole_tissue) <- old_idents

# 2. GLIOMAS: T.I. Intact vs N.T. Intact

cells_TI_intact_Gliomas <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Gliomas" &
    whole_tissue_meta$new_clusters == "Pericyte" &
    whole_tissue_meta$NVU_States == "T.I.: Intact"
]

cells_NT_intact_Gliomas <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Gliomas" &
    whole_tissue_meta$new_clusters == "Pericyte" &
    whole_tissue_meta$NVU_States == "N.T.: Intact"
]

message(
  "Gliomas: N.T. Intact = ",
  length(cells_NT_intact_Gliomas),
  "; T.I. Intact = ",
  length(cells_TI_intact_Gliomas)
)

old_idents <- Idents(whole_tissue)

Idents(
  whole_tissue,
  cells = cells_TI_intact_Gliomas
) <- "T.I. Intact"

Idents(
  whole_tissue,
  cells = cells_NT_intact_Gliomas
) <- "N.T. Intact"

pericyte_NT_vs_TI_by_condition[["Gliomas"]] <- FindMarkers(
  whole_tissue,
  ident.1 = "T.I. Intact",
  ident.2 = "N.T. Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

Idents(whole_tissue) <- old_idents

# 3. BRM: T.I. Intact vs N.T. Intact

cells_TI_intact_BrM <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "BrM" &
    whole_tissue_meta$new_clusters == "Pericyte" &
    whole_tissue_meta$NVU_States == "T.I.: Intact"
]

cells_NT_intact_BrM <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "BrM" &
    whole_tissue_meta$new_clusters == "Pericyte" &
    whole_tissue_meta$NVU_States == "N.T.: Intact"
]

message(
  "BrM: N.T. Intact = ",
  length(cells_NT_intact_BrM),
  "; T.I. Intact = ",
  length(cells_TI_intact_BrM)
)

old_idents <- Idents(whole_tissue)

Idents(
  whole_tissue,
  cells = cells_TI_intact_BrM
) <- "T.I. Intact"

Idents(
  whole_tissue,
  cells = cells_NT_intact_BrM
) <- "N.T. Intact"

pericyte_NT_vs_TI_by_condition[["BrM"]] <- FindMarkers(
  whole_tissue,
  ident.1 = "T.I. Intact",
  ident.2 = "N.T. Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

Idents(whole_tissue) <- old_idents

# 4. GLIOMAS: T.I. Astrocyte displaced vs T.I. Intact

cells_TI_astro_Gliomas <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Gliomas" &
    whole_tissue_meta$new_clusters == "Pericyte" &
    whole_tissue_meta$NVU_States == "T.I.: Astrocyte displaced"
]

cells_TI_intact_astro_Gliomas <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Gliomas" &
    whole_tissue_meta$new_clusters == "Pericyte" &
    whole_tissue_meta$NVU_States == "T.I.: Intact"
]

message(
  "Gliomas: T.I. Intact = ",
  length(cells_TI_intact_astro_Gliomas),
  "; T.I. Astrocyte displaced = ",
  length(cells_TI_astro_Gliomas)
)

old_idents <- Idents(whole_tissue)

Idents(
  whole_tissue,
  cells = cells_TI_astro_Gliomas
) <- "T.I. Astrocyte displaced"

Idents(
  whole_tissue,
  cells = cells_TI_intact_astro_Gliomas
) <- "T.I. Intact"

pericyte_TI_vs_astro_by_condition[["Gliomas"]] <- FindMarkers(
  whole_tissue,
  ident.1 = "T.I. Astrocyte displaced",
  ident.2 = "T.I. Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

Idents(whole_tissue) <- old_idents

# 5. BRM: T.I. Astrocyte displaced vs T.I. Intact

cells_TI_astro_BrM <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "BrM" &
    whole_tissue_meta$new_clusters == "Pericyte" &
    whole_tissue_meta$NVU_States == "T.I.: Astrocyte displaced"
]

cells_TI_intact_astro_BrM <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "BrM" &
    whole_tissue_meta$new_clusters == "Pericyte" &
    whole_tissue_meta$NVU_States == "T.I.: Intact"
]

message(
  "BrM: T.I. Intact = ",
  length(cells_TI_intact_astro_BrM),
  "; T.I. Astrocyte displaced = ",
  length(cells_TI_astro_BrM)
)

old_idents <- Idents(whole_tissue)

Idents(
  whole_tissue,
  cells = cells_TI_astro_BrM
) <- "T.I. Astrocyte displaced"

Idents(
  whole_tissue,
  cells = cells_TI_intact_astro_BrM
) <- "T.I. Intact"

pericyte_TI_vs_astro_by_condition[["BrM"]] <- FindMarkers(
  whole_tissue,
  ident.1 = "T.I. Astrocyte displaced",
  ident.2 = "T.I. Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

Idents(whole_tissue) <- old_idents

# Filter DEG results to curated NVU genes

pericyte_NT_vs_TI_NVU_by_condition <- lapply(
  pericyte_NT_vs_TI_by_condition,
  function(x) {
    
    if (is.null(x)) {
      return(NULL)
    }
    
    x %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::filter(gene %in% nvu_genes)
  }
)


pericyte_TI_vs_astro_NVU_by_condition <- lapply(
  pericyte_TI_vs_astro_by_condition,
  function(x) {
    
    if (is.null(x)) {
      return(NULL)
    }
    
    x %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::filter(gene %in% nvu_genes)
  }
)


# Inspect results
pericyte_NT_vs_TI_NVU_by_condition

pericyte_TI_vs_astro_NVU_by_condition

# Volcano plot function

prepare_volcano <- function(deg_df) {
  
  if (is.null(deg_df) || nrow(deg_df) == 0) {
    return(NULL)
  }
  
  deg_df %>%
    dplyr::mutate(
      neg_log10_padj = -log10(
        pmax(p_val_adj, .Machine$double.xmin)
      ),
      significance = dplyr::case_when(
        p_val_adj < 0.05 &
          avg_log2FC >= 1 ~ "Up",
        
        p_val_adj < 0.05 &
          avg_log2FC <= -1 ~ "Down",
        
        TRUE ~ "Not significant"
      )
    )
}


make_volcano <- function(
    deg_df,
    plot_title
) {
  
  if (is.null(deg_df) || nrow(deg_df) == 0) {
    
    return(
      ggplot() +
        theme_void() +
        ggtitle(plot_title) +
        theme(
          plot.title = element_text(
            size = 16,
            face = "bold",
            hjust = 0.5
          )
        )
    )
  }
  
  
  plot_df <- prepare_volcano(deg_df)
  
  
  # Top 10 genes by positive log2FC
  top_up <- plot_df %>%
    dplyr::filter(avg_log2FC > 0) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::slice_head(n = 10)
  
  
  # Top 10 genes by negative log2FC
  top_down <- plot_df %>%
    dplyr::filter(avg_log2FC < 0) %>%
    dplyr::arrange(avg_log2FC) %>%
    dplyr::slice_head(n = 10)
  
  
  labels_df <- dplyr::bind_rows(
    top_up,
    top_down
  )
  
  
  ggplot(
    plot_df,
    aes(
      x = avg_log2FC,
      y = neg_log10_padj
    )
  ) +
    
    geom_point(
      aes(color = significance),
      size = 2.5,
      alpha = 0.8
    ) +
    
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dotted",
      linewidth = 0.6
    ) +
    
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    ggrepel::geom_text_repel(
      data = labels_df,
      aes(label = gene),
      size = 4,
      max.overlaps = 20,
      box.padding = 0.4,
      point.padding = 0.2,
      show.legend = FALSE
    ) +
    
    scale_color_manual(
      values = c(
        "Up" = "#D73027",
        "Down" = "#4575B4",
        "Not significant" = "grey70"
      ),
      labels = c(
        "Up" = "Upregulated",
        "Down" = "Downregulated",
        "Not significant" = "Not significant"
      )
    ) +
    
    labs(
      title = plot_title,
      x = "Average log2 fold change",
      y = expression(-log[10]("adjusted P")),
      color = "Significance"
    ) +
    
    theme_classic() +
    
    theme(
      # Plot title
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      
      # Axis titles
      axis.title.x = element_text(
        size = 14,
        face = "bold"
      ),
      
      axis.title.y = element_text(
        size = 14,
        face = "bold"
      ),
      
      # Axis tick labels
      axis.text.x = element_text(
        size = 12
      ),
      
      axis.text.y = element_text(
        size = 12
      ),
      
      # Legend
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      ),
      
      legend.key.size = grid::unit(
        0.7,
        "cm"
      ),
      
      # Plot margins
      plot.margin = margin(
        10,
        10,
        10,
        10
      )
    )
}

# Create volcano plots

# T.I. Intact vs N.T. Intact
volcano_NT_vs_TI <- lapply(
  conditions_to_analyze,
  function(cond) {
    
    make_volcano(
      pericyte_NT_vs_TI_NVU_by_condition[[cond]],
      paste0(
        cond,
        ": Pericytes\n",
        "T.I. Intact vs N.T. Intact"
      )
    )
  }
)

names(volcano_NT_vs_TI) <- conditions_to_analyze


# T.I. Astrocyte displaced vs T.I. Intact
volcano_TI_vs_astro <- lapply(
  conditions_to_analyze,
  function(cond) {
    
    make_volcano(
      pericyte_TI_vs_astro_NVU_by_condition[[cond]],
      paste0(
        cond,
        ": Pericytes\n",
        "T.I. Astrocyte displaced vs T.I. Intact"
      )
    )
  }
)

names(volcano_TI_vs_astro) <- conditions_to_analyze

# Display plots

volcano_NT_vs_TI$`Non-Tumor`
volcano_NT_vs_TI$Gliomas
volcano_NT_vs_TI$BrM

volcano_TI_vs_astro$`Non-Tumor`
volcano_TI_vs_astro$Gliomas
volcano_TI_vs_astro$BrM

# Save volcano plots

ggsave(
  "Vascular subtype analysis/Figures/Pericyte_NT_vs_TI_Non-Tumor_volcano.svg",
  volcano_NT_vs_TI$`Non-Tumor`,
  width = 7,
  height = 6
)

ggsave(
  "Vascular subtype analysis/Figures/Pericyte_NT_vs_TI_Gliomas_volcano.svg",
  volcano_NT_vs_TI$Gliomas,
  width = 7,
  height = 6
)

ggsave(
  "Vascular subtype analysis/Figures/Pericyte_NT_vs_TI_BrM_volcano.svg",
  volcano_NT_vs_TI$BrM,
  width = 7,
  height = 6
)

ggsave(
  "Vascular subtype analysis/Figures/Pericyte_TI_vs_Astro_Non-Tumor_volcano.svg",
  volcano_TI_vs_astro$`Non-Tumor`,
  width = 7,
  height = 6
)

ggsave(
  "Vascular subtype analysis/Figures/Pericyte_TI_vs_Astro_Gliomas_volcano.svg",
  volcano_TI_vs_astro$Gliomas,
  width = 7,
  height = 6
)

ggsave(
  "Vascular subtype analysis/Figures/Pericyte_TI_vs_Astro_BrM_volcano.svg",
  volcano_TI_vs_astro$BrM,
  width = 7,
  height = 6
)


# Ast-like DEG analysis across NVU states and conditions----

# IMPORTANT:
# Work from metadata/cell IDs rather than creating Seurat subsets.
# This avoids the memory exhaustion problem encountered with subset().


# Get metadata directly from the existing Seurat object
whole_tissue_meta <- whole_tissue[[]]


# Make condition labels
whole_tissue_meta$new_conditons <- dplyr::case_when(
  whole_tissue_meta$condition == "Normal" ~ "Non-Tumor",
  whole_tissue_meta$condition == "Primary Tumor" ~ "Gliomas",
  whole_tissue_meta$condition == "Secondary Tumor" ~ "BrM",
  TRUE ~ NA_character_
)


# Curated Ast-like NVU gene set
nvu_genes_astlike <- unique(c(
  "ALDH1L1",
  "AQP4",
  "APOE",
  "BCAN",
  "C1QL3",
  "GFAP",
  "GJA1",
  "MLC1",
  "PTPRZ1",
  "SLC1A1",
  "SLC1A2",
  "SLC1A3",
  "SLC4A4",
  "SOX9",
  "VIM",
  "ANGPT1",
  "ANGPT2",
  "CCL2",
  "COL4A1",
  "FZD4",
  "ICAM1",
  "IGFBP3",
  "IGFBP4",
  "IGFBP5",
  "IL1B",
  "IL6",
  "LRP1",
  "MFSD2A",
  "NPNT",
  "SLC2A1",
  "SLC2A3",
  "SLC16A1",
  "SLC7A5",
  "TGFB1",
  "TGFB2",
  "TGFBI",
  "TGFBR1",
  "TGFBR2",
  "THBS1",
  "VCAM1",
  "VEGFA",
  "WIF1",
  "WNT7A",
  "CEMIP",
  "CEMIP2",
  "CXCL14",
  "FABP5",
  "HILPDA",
  "SFRP2",
  "SERPINA3",
  "STAT3",
  "SULF1",
  "VCAN"
))

nvu_genes_astlike <- intersect(
  nvu_genes_astlike,
  rownames(whole_tissue)
)


# Check required metadata
stopifnot(
  all(c(
    "new_conditons",
    "NVU_States",
    "new_clusters"
  ) %in% colnames(whole_tissue_meta))
)


# Conditions
conditions_to_analyze <- c(
  "Non-Tumor",
  "Gliomas",
  "BrM"
)


# Initialize result lists
astlike_NT_vs_TI_by_condition <- list(
  "Non-Tumor" = NULL,
  "Gliomas" = NULL,
  "BrM" = NULL
)

astlike_TI_vs_pericyte_by_condition <- list(
  "Non-Tumor" = NULL,
  "Gliomas" = NULL,
  "BrM" = NULL
)


# 1. NON-TUMOR: T.I. Intact vs N.T. Intact

cells_TI_intact_NonTumor <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Non-Tumor" &
    whole_tissue_meta$new_clusters == "Ast-like" &
    whole_tissue_meta$NVU_States == "T.I.: Intact"
]

cells_NT_intact_NonTumor <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Non-Tumor" &
    whole_tissue_meta$new_clusters == "Ast-like" &
    whole_tissue_meta$NVU_States == "N.T.: Intact"
]

message(
  "Non-Tumor: N.T. Intact = ",
  length(cells_NT_intact_NonTumor),
  "; T.I. Intact = ",
  length(cells_TI_intact_NonTumor)
)

old_idents <- Idents(whole_tissue)

Idents(
  whole_tissue,
  cells = cells_TI_intact_NonTumor
) <- "T.I. Intact"

Idents(
  whole_tissue,
  cells = cells_NT_intact_NonTumor
) <- "N.T. Intact"

astlike_NT_vs_TI_by_condition[["Non-Tumor"]] <- FindMarkers(
  whole_tissue,
  ident.1 = "T.I. Intact",
  ident.2 = "N.T. Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

Idents(whole_tissue) <- old_idents


# 2. GLIOMAS: T.I. Intact vs N.T. Intact

cells_TI_intact_Gliomas <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Gliomas" &
    whole_tissue_meta$new_clusters == "Ast-like" &
    whole_tissue_meta$NVU_States == "T.I.: Intact"
]

cells_NT_intact_Gliomas <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Gliomas" &
    whole_tissue_meta$new_clusters == "Ast-like" &
    whole_tissue_meta$NVU_States == "N.T.: Intact"
]

message(
  "Gliomas: N.T. Intact = ",
  length(cells_NT_intact_Gliomas),
  "; T.I. Intact = ",
  length(cells_TI_intact_Gliomas)
)

old_idents <- Idents(whole_tissue)

Idents(
  whole_tissue,
  cells = cells_TI_intact_Gliomas
) <- "T.I. Intact"

Idents(
  whole_tissue,
  cells = cells_NT_intact_Gliomas
) <- "N.T. Intact"

astlike_NT_vs_TI_by_condition[["Gliomas"]] <- FindMarkers(
  whole_tissue,
  ident.1 = "T.I. Intact",
  ident.2 = "N.T. Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

Idents(whole_tissue) <- old_idents


# 3. BRM: T.I. Intact vs N.T. Intact

cells_TI_intact_BrM <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "BrM" &
    whole_tissue_meta$new_clusters == "Ast-like" &
    whole_tissue_meta$NVU_States == "T.I.: Intact"
]

cells_NT_intact_BrM <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "BrM" &
    whole_tissue_meta$new_clusters == "Ast-like" &
    whole_tissue_meta$NVU_States == "N.T.: Intact"
]

message(
  "BrM: N.T. Intact = ",
  length(cells_NT_intact_BrM),
  "; T.I. Intact = ",
  length(cells_TI_intact_BrM)
)

old_idents <- Idents(whole_tissue)

Idents(
  whole_tissue,
  cells = cells_TI_intact_BrM
) <- "T.I. Intact"

Idents(
  whole_tissue,
  cells = cells_NT_intact_BrM
) <- "N.T. Intact"

astlike_NT_vs_TI_by_condition[["BrM"]] <- FindMarkers(
  whole_tissue,
  ident.1 = "T.I. Intact",
  ident.2 = "N.T. Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

Idents(whole_tissue) <- old_idents


# 4. GLIOMAS: T.I. Pericyte displaced vs T.I. Intact

cells_TI_pericyte_Gliomas <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Gliomas" &
    whole_tissue_meta$new_clusters == "Ast-like" &
    whole_tissue_meta$NVU_States == "T.I.: Pericyte displaced"
]

cells_TI_intact_pericyte_Gliomas <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "Gliomas" &
    whole_tissue_meta$new_clusters == "Ast-like" &
    whole_tissue_meta$NVU_States == "T.I.: Intact"
]

message(
  "Gliomas: T.I. Intact = ",
  length(cells_TI_intact_pericyte_Gliomas),
  "; T.I. Pericyte displaced = ",
  length(cells_TI_pericyte_Gliomas)
)

old_idents <- Idents(whole_tissue)

Idents(
  whole_tissue,
  cells = cells_TI_pericyte_Gliomas
) <- "T.I. Pericyte displaced"

Idents(
  whole_tissue,
  cells = cells_TI_intact_pericyte_Gliomas
) <- "T.I. Intact"

astlike_TI_vs_pericyte_by_condition[["Gliomas"]] <- FindMarkers(
  whole_tissue,
  ident.1 = "T.I. Pericyte displaced",
  ident.2 = "T.I. Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

Idents(whole_tissue) <- old_idents


# 5. BRM: T.I. Pericyte displaced vs T.I. Intact

cells_TI_pericyte_BrM <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "BrM" &
    whole_tissue_meta$new_clusters == "Ast-like" &
    whole_tissue_meta$NVU_States == "T.I.: Pericyte displaced"
]

cells_TI_intact_pericyte_BrM <- rownames(whole_tissue_meta)[
  whole_tissue_meta$new_conditons == "BrM" &
    whole_tissue_meta$new_clusters == "Ast-like" &
    whole_tissue_meta$NVU_States == "T.I.: Intact"
]

message(
  "BrM: T.I. Intact = ",
  length(cells_TI_intact_pericyte_BrM),
  "; T.I. Pericyte displaced = ",
  length(cells_TI_pericyte_BrM)
)

old_idents <- Idents(whole_tissue)

Idents(
  whole_tissue,
  cells = cells_TI_pericyte_BrM
) <- "T.I. Pericyte displaced"

Idents(
  whole_tissue,
  cells = cells_TI_intact_pericyte_BrM
) <- "T.I. Intact"

astlike_TI_vs_pericyte_by_condition[["BrM"]] <- FindMarkers(
  whole_tissue,
  ident.1 = "T.I. Pericyte displaced",
  ident.2 = "T.I. Intact",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.10
)

Idents(whole_tissue) <- old_idents


# Filter DEG results to curated Ast-like NVU genes

astlike_NT_vs_TI_NVU_by_condition <- lapply(
  astlike_NT_vs_TI_by_condition,
  function(x) {
    
    if (is.null(x)) {
      return(NULL)
    }
    
    x %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::filter(gene %in% nvu_genes_astlike)
  }
)


astlike_TI_vs_pericyte_NVU_by_condition <- lapply(
  astlike_TI_vs_pericyte_by_condition,
  function(x) {
    
    if (is.null(x)) {
      return(NULL)
    }
    
    x %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::filter(gene %in% nvu_genes_astlike)
  }
)


# Inspect results
astlike_NT_vs_TI_NVU_by_condition

astlike_TI_vs_pericyte_NVU_by_condition


# Volcano plot function
# Same overall structure as the Pericyte function,
# with the plotting calculation made explicit.

prepare_volcano <- function(deg_df) {
  
  if (is.null(deg_df) || nrow(deg_df) == 0) {
    return(NULL)
  }
  
  deg_df %>%
    dplyr::mutate(
      neg_log10_padj = -log10(
        pmax(p_val_adj, .Machine$double.xmin)
      ),
      significance = dplyr::case_when(
        p_val_adj < 0.05 &
          avg_log2FC >= 1 ~ "Up",
        
        p_val_adj < 0.05 &
          avg_log2FC <= -1 ~ "Down",
        
        TRUE ~ "Not significant"
      )
    )
}


make_volcano <- function(
    deg_df,
    plot_title
) {
  
  if (is.null(deg_df) || nrow(deg_df) == 0) {
    
    return(
      ggplot() +
        theme_void() +
        ggtitle(plot_title) +
        theme(
          plot.title = element_text(
            size = 16,
            face = "bold",
            hjust = 0.5
          )
        )
    )
  }
  
  
  plot_df <- prepare_volcano(deg_df)
  
  
  # Top 10 genes by positive log2FC
  top_up <- plot_df %>%
    dplyr::filter(avg_log2FC > 0) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::slice_head(n = 10)
  
  
  # Top 10 genes by negative log2FC
  top_down <- plot_df %>%
    dplyr::filter(avg_log2FC < 0) %>%
    dplyr::arrange(avg_log2FC) %>%
    dplyr::slice_head(n = 10)
  
  
  labels_df <- dplyr::bind_rows(
    top_up,
    top_down
  )
  
  
  ggplot(
    plot_df,
    aes(
      x = avg_log2FC,
      y = neg_log10_padj
    )
  ) +
    
    geom_point(
      aes(color = significance),
      size = 2.5,
      alpha = 0.8
    ) +
    
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dotted",
      linewidth = 0.6
    ) +
    
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    ggrepel::geom_text_repel(
      data = labels_df,
      aes(label = gene),
      size = 4,
      max.overlaps = 20,
      box.padding = 0.4,
      point.padding = 0.2,
      show.legend = FALSE
    ) +
    
    scale_color_manual(
      values = c(
        "Up" = "#D73027",
        "Down" = "#4575B4",
        "Not significant" = "grey70"
      ),
      labels = c(
        "Up" = "Upregulated",
        "Down" = "Downregulated",
        "Not significant" = "Not significant"
      )
    ) +
    
    labs(
      title = plot_title,
      x = "Average log2 fold change",
      y = expression(-log[10]("adjusted P")),
      color = "Significance"
    ) +
    
    theme_classic() +
    
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      
      axis.title.x = element_text(
        size = 14,
        face = "bold"
      ),
      
      axis.title.y = element_text(
        size = 14,
        face = "bold"
      ),
      
      axis.text.x = element_text(
        size = 12
      ),
      
      axis.text.y = element_text(
        size = 12
      ),
      
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 12
      ),
      
      legend.key.size = grid::unit(
        0.7,
        "cm"
      ),
      
      plot.margin = margin(
        10,
        10,
        10,
        10
      )
    )
}


# Create volcano plots

# T.I. Intact vs N.T. Intact
volcano_astlike_NT_vs_TI <- lapply(
  conditions_to_analyze,
  function(cond) {
    
    make_volcano(
      astlike_NT_vs_TI_NVU_by_condition[[cond]],
      paste0(
        cond,
        ": Ast-like\n",
        "T.I. Intact vs N.T. Intact"
      )
    )
  }
)

names(volcano_astlike_NT_vs_TI) <- conditions_to_analyze


# T.I. Pericyte displaced vs T.I. Intact
volcano_astlike_TI_vs_pericyte <- lapply(
  conditions_to_analyze,
  function(cond) {
    
    make_volcano(
      astlike_TI_vs_pericyte_NVU_by_condition[[cond]],
      paste0(
        cond,
        ": Ast-like\n",
        "T.I. Pericyte displaced vs T.I. Intact"
      )
    )
  }
)

names(volcano_astlike_TI_vs_pericyte) <- conditions_to_analyze


# Display plots

volcano_astlike_NT_vs_TI$`Non-Tumor`
volcano_astlike_NT_vs_TI$Gliomas
volcano_astlike_NT_vs_TI$BrM

volcano_astlike_TI_vs_pericyte$`Non-Tumor`
volcano_astlike_TI_vs_pericyte$Gliomas
volcano_astlike_TI_vs_pericyte$BrM


# Save volcano plots

ggsave(
  "Vascular subtype analysis/Figures/Ast-like_NT_vs_TI_Non-Tumor_volcano.svg",
  volcano_astlike_NT_vs_TI$`Non-Tumor`,
  width = 7,
  height = 6
)

ggsave(
  "Vascular subtype analysis/Figures/Ast-like_NT_vs_TI_Gliomas_volcano.svg",
  volcano_astlike_NT_vs_TI$Gliomas,
  width = 7,
  height = 6
)

ggsave(
  "Vascular subtype analysis/Figures/Ast-like_NT_vs_TI_BrM_volcano.svg",
  volcano_astlike_NT_vs_TI$BrM,
  width = 7,
  height = 6
)

ggsave(
  "Vascular subtype analysis/Figures/Ast-like_TI_vs_Pericyte_Non-Tumor_volcano.svg",
  volcano_astlike_TI_vs_pericyte$`Non-Tumor`,
  width = 7,
  height = 6
)

ggsave(
  "Vascular subtype analysis/Figures/Ast-like_TI_vs_Pericyte_Gliomas_volcano.svg",
  volcano_astlike_TI_vs_pericyte$Gliomas,
  width = 7,
  height = 6
)

ggsave(
  "Vascular subtype analysis/Figures/Ast-like_TI_vs_Pericyte_BrM_volcano.svg",
  volcano_astlike_TI_vs_pericyte$BrM,
  width = 7,
  height = 6
)


#Ast-like BBB Support Module genes----
# Ast-like cells: BBB support and ECM remodeling gene expression
# Z-scored average expression across samples
# Samples ordered: Normal -> Grade 3 -> GBM -> metastatic samples

library(Seurat)
library(dplyr)
library(ggplot2)

# Define modules
ast_modules <- list(
  "BBB Support" = c(
    "AQP4",
    "GJA1",
    "APOE",
    "ANGPT1",
    "MLC1"
  ),
  "ECM Remodeling" = c(
    "MMP2",
    "MMP9",
    "VCAN"
  )
)

# All genes
ast_genes <- unlist(ast_modules)

# Get metadata
whole_tissue_meta <- whole_tissue[[]]

# Identify Ast-like cells
ast_meta <- whole_tissue_meta %>%
  tibble::rownames_to_column("cell") %>%
  filter(new_clusters == "Ast-like") %>%
  mutate(
    sample = if_else(
      sample == "Gr_3",
      "Grade 3",
      as.character(sample)
    )
  )

# Keep genes present in the object
ast_genes <- ast_genes[ast_genes %in% rownames(whole_tissue)]

# Extract only the required genes and Ast-like cells
# This avoids creating a large Seurat subset
expr_mat <- GetAssayData(
  whole_tissue,
  assay = "Xenium",
  layer = "data"
)[ast_genes, ast_meta$cell, drop = FALSE]

# Keep metadata aligned with expression matrix
ast_meta <- ast_meta %>%
  filter(cell %in% colnames(expr_mat))

expr_mat <- expr_mat[, ast_meta$cell, drop = FALSE]

# Calculate average expression and percent expressing per sample
dotplot_df <- lapply(
  seq_along(ast_modules),
  function(i) {
    
    module_name <- names(ast_modules)[i]
    module_genes <- ast_modules[[i]]
    module_genes <- module_genes[module_genes %in% rownames(expr_mat)]
    
    lapply(
      module_genes,
      function(gene) {
        
        gene_expr <- as.numeric(expr_mat[gene, ])
        
        tibble(
          sample = ast_meta$sample,
          expression = gene_expr
        ) %>%
          group_by(sample) %>%
          summarise(
            avg_expression = mean(expression, na.rm = TRUE),
            pct_expressed = 100 * mean(expression > 0, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          mutate(
            gene = gene,
            module = module_name
          )
      }
    ) %>%
      bind_rows()
  }
) %>%
  bind_rows()

# Z-score average expression separately for each gene across samples
dotplot_df <- dotplot_df %>%
  group_by(gene) %>%
  mutate(
    avg_expression_z = as.numeric(scale(avg_expression))
  ) %>%
  ungroup()

# Gene order
dotplot_df$gene <- factor(
  dotplot_df$gene,
  levels = ast_genes
)

# Module order
dotplot_df$module <- factor(
  dotplot_df$module,
  levels = c(
    "BBB Support",
    "ECM Remodeling"
  )
)

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

# Keep only samples present in the data
sample_order <- sample_order[
  sample_order %in% unique(as.character(dotplot_df$sample))
]

dotplot_df$sample <- factor(
  dotplot_df$sample,
  levels = rev(sample_order)
)

# Plot
p_Ast_modules_by_sample <- ggplot(
  dotplot_df,
  aes(
    x = gene,
    y = sample
  )
) +
  geom_point(
    aes(
      size = pct_expressed,
      color = avg_expression_z
    )
  ) +
  facet_grid(
    . ~ module,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_size(
    range = c(1, 7),
    name = "Percent expressed"
  ) +
  scale_color_gradient2(
    name = "Z-scored\naverage expression",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(
    title = "Ast-like BBB support and ECM remodeling programs",
    x = NULL,
    y = "Sample"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(
      size = 16,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 16),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_Ast_modules_by_sample)

# Save as SVG
ggsave(
  "Vascular subtype analysis/Figures/Ast_like_BBB_ECM_modules_by_sample_Zscore.svg",
  p_Ast_modules_by_sample,
  width = 12,
  height = 5,
  device = "svg"
)


# Ast-like UCell Module Scores----
# Calculate BBB Support and ECM Remodeling UCell scores
# directly on whole_tissue, then focus all downstream analysis
# on Ast-like cells.

library(Seurat)
library(UCell)
library(dplyr)
library(ggplot2)

# Output directory
dir.create(
  "Vascular subtype analysis/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# Define gene modules
BBB_Support_genes <- c(
  "AQP4",
  "GJA1",
  "APOE",
  "ANGPT1",
  "MLC1"
)

ECM_Remodeling_genes <- c(
  "MMP2",
  "MMP9",
  "VCAN"
)

# Check genes
cat("Checking genes...\n")

missing_BBB <- setdiff(
  BBB_Support_genes,
  rownames(whole_tissue)
)

missing_ECM <- setdiff(
  ECM_Remodeling_genes,
  rownames(whole_tissue)
)

if (length(missing_BBB) > 0) {
  stop(
    "Missing BBB Support genes: ",
    paste(missing_BBB, collapse = ", ")
  )
}

if (length(missing_ECM) > 0) {
  stop(
    "Missing ECM Remodeling genes: ",
    paste(missing_ECM, collapse = ", ")
  )
}

cat("All genes found.\n\n")

# Define UCell modules
ast_modules <- list(
  BBB_Support = BBB_Support_genes,
  ECM_Remodeling = ECM_Remodeling_genes
)

# Calculate UCell scores directly on whole_tissue
#
# This creates:
#   BBB_Support_UCell
#   ECM_Remodeling_UCell

cat("Calculating UCell module scores on whole_tissue...\n")

whole_tissue <- UCell::AddModuleScore_UCell(
  whole_tissue,
  features = ast_modules,
  name = NULL
)

cat("\nUCell scoring complete.\n")

# Check the resulting metadata columns
cat("\nUCell-related metadata columns:\n")

print(
  colnames(whole_tissue[[]])[
    grepl(
      "BBB|ECM|Support|Remodel",
      colnames(whole_tissue[[]]),
      ignore.case = TRUE
    )
  ]
)

# Rename UCell columns if necessary
#
# With name = NULL, these should normally be:
#   BBB_Support
#   ECM_Remodeling

if ("BBB_Support" %in% colnames(whole_tissue[[]])) {
  
  whole_tissue$BBB_Support_UCell <- whole_tissue$BBB_Support
  
}

if ("ECM_Remodeling" %in% colnames(whole_tissue[[]])) {
  
  whole_tissue$ECM_Remodeling_UCell <- whole_tissue$ECM_Remodeling
  
}

# Confirm the two desired score columns exist
if (!all(
  c(
    "BBB_Support_UCell",
    "ECM_Remodeling_UCell"
  ) %in% colnames(whole_tissue[[]])
)) {
  
  stop(
    "Could not find the expected UCell score columns. ",
    "Check the printed UCell-related metadata columns above."
  )
}

# Create composite Ast-like BBB Support Score
#
# Positive values:
#   relatively stronger BBB Support
#
# Negative values:
#   relatively stronger ECM Remodeling

whole_tissue$Ast_BBB_Support_Score <-
  whole_tissue$BBB_Support_UCell -
  whole_tissue$ECM_Remodeling_UCell

cat("\nFinal score columns:\n")

print(
  c(
    "BBB_Support_UCell",
    "ECM_Remodeling_UCell",
    "Ast_BBB_Support_Score"
  )
)

# Check overall score distributions
cat("\nOverall BBB Support score:\n")
print(
  summary(
    whole_tissue$BBB_Support_UCell
  )
)

cat("\nOverall ECM Remodeling score:\n")
print(
  summary(
    whole_tissue$ECM_Remodeling_UCell
  )
)

cat("\nOverall Ast BBB Support score:\n")
print(
  summary(
    whole_tissue$Ast_BBB_Support_Score
  )
)


# Get metadata only
#
# From this point onward, we work only with the metadata.
# We do NOT create another Seurat object.

whole_tissue_meta <- whole_tissue[[]]

# Focus on Ast-like cells
astlike_scores <- whole_tissue_meta %>%
  tibble::rownames_to_column("cell") %>%
  filter(
    new_clusters == "Ast-like"
  ) %>%
  select(
    cell,
    sample,
    condition,
    NVU_States,
    new_clusters,
    BBB_Support_UCell,
    ECM_Remodeling_UCell,
    Ast_BBB_Support_Score
  )

cat(
  "\nNumber of Ast-like cells:",
  nrow(astlike_scores),
  "\n"
)

# Check NVU states among Ast-like cells
cat("\nAst-like cells by NVU state:\n")

print(
  table(
    astlike_scores$NVU_States,
    useNA = "ifany"
  )
)

#NVU States for Ast-like cells
# Keep only the three NVU states of interest
#
# Actual labels present in whole_tissue:
#   N.T.: Intact
#   T.I.: Intact
#   T.I.: Pericyte displaced

nvu_states_ast <- c(
  "N.T.: Intact",
  "T.I.: Intact",
  "T.I.: Pericyte displaced"
)

astlike_scores <- astlike_scores %>%
  filter(
    NVU_States %in% nvu_states_ast
  )

cat(
  "\nAst-like cells in the three selected NVU states:",
  nrow(astlike_scores),
  "\n"
)

# Check cells by sample and NVU state
cat("\nAst-like cell counts by sample and NVU state:\n")

astlike_cell_counts <- astlike_scores %>%
  count(
    sample,
    condition,
    NVU_States,
    name = "n_cells"
  )

print(astlike_cell_counts)

# Calculate sample-level means
#
# Each sample is the biological replicate.
# Cells are first averaged within each sample/NVU state.

astlike_sample_summary <- astlike_scores %>%
  group_by(
    sample,
    condition,
    NVU_States
  ) %>%
  summarise(
    
    n_cells = n(),
    
    BBB_Support_mean =
      mean(
        BBB_Support_UCell,
        na.rm = TRUE
      ),
    
    ECM_Remodeling_mean =
      mean(
        ECM_Remodeling_UCell,
        na.rm = TRUE
      ),
    
    Ast_BBB_Support_Score_mean =
      mean(
        Ast_BBB_Support_Score,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

# Set NVU state order
astlike_sample_summary <- astlike_sample_summary %>%
  mutate(
    NVU_States = factor(
      NVU_States,
      levels = nvu_states_ast
    )
  )

# Print sample-level results
cat("\nSample-level Ast-like results:\n")

print(
  astlike_sample_summary
)

# Number of biological replicates per state
cat("\nNumber of samples per NVU state:\n")

print(
  astlike_sample_summary %>%
    count(NVU_States)
)

# Save sample-level summary
write.csv(
  astlike_sample_summary,
  "Vascular subtype analysis/Astlike_BBB_Support_sample_summary.csv",
  row.names = FALSE
)

# Save Ast-like cell-level scores
write.csv(
  astlike_scores,
  "Vascular subtype analysis/Astlike_BBB_Support_cell_scores.csv",
  row.names = FALSE
)

# Create plotting labels
#
# These are shortened labels for visualization only.
# The original NVU_States values are retained above.

astlike_sample_summary_plot <- astlike_sample_summary %>%
  mutate(
    
    NVU_State_Label = recode(
      as.character(NVU_States),
      
      "N.T.: Intact" =
        "N.T. Intact",
      
      "T.I.: Intact" =
        "T.I. Intact",
      
      "T.I.: Pericyte displaced" =
        "T.I. Pericyte displaced"
    ),
    
    NVU_State_Label = factor(
      NVU_State_Label,
      levels = c(
        "N.T. Intact",
        "T.I. Intact",
        "T.I. Pericyte displaced"
      )
    )
  )
# Ast-like NVU state plots
#
# Each point = one biological sample.
# Each line connects the same biological sample across NVU states.

library(ggplot2)

#----
# Ast-like BBB Support
# Sample-level MEAN scores across the three NVU states
#
# Each point = one biological sample/NVU state
# Each line = one biological sample
# Color = sample
# n_cells = number of Ast-like cells contributing to each point

library(dplyr)
library(ggplot2)

# Define NVU state order
nvu_state_levels <- c(
  "N.T. Intact",
  "T.I. Intact",
  "T.I. Pericyte displaced"
)

# Prepare cell-level data
astlike_BBB_mean_plot_data <- astlike_scores %>%
  mutate(
    NVU_State_Label = recode(
      as.character(NVU_States),
      "N.T.: Intact" = "N.T. Intact",
      "T.I.: Intact" = "T.I. Intact",
      "T.I.: Pericyte displaced" = "T.I. Pericyte displaced"
    ),
    NVU_State_Label = factor(
      NVU_State_Label,
      levels = nvu_state_levels
    )
  ) %>%
  filter(
    !is.na(BBB_Support_UCell),
    !is.na(NVU_State_Label)
  )


# Calculate sample-level MEAN BBB Support
astlike_BBB_sample_mean <- astlike_BBB_mean_plot_data %>%
  group_by(sample, condition, NVU_State_Label) %>%
  summarise(
    n_cells = n(),
    BBB_Support_mean = mean(
      BBB_Support_UCell,
      na.rm = TRUE
    ),
    BBB_Support_median = median(
      BBB_Support_UCell,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# Order samples
sample_order_mean <- astlike_BBB_sample_mean %>%
  distinct(sample, condition) %>%
  arrange(condition, sample) %>%
  pull(sample)

astlike_BBB_sample_mean <- astlike_BBB_sample_mean %>%
  mutate(
    sample = factor(
      sample,
      levels = sample_order_mean
    )
  )


# Main plot
#
# All samples shown together.
# Lines connect the same sample across NVU states.
# Color identifies the sample.

p_Ast_BBB_Support_sample_mean <- ggplot(
  astlike_BBB_sample_mean,
  aes(
    x = NVU_State_Label,
    y = BBB_Support_mean,
    group = sample,
    color = sample
  )
) +
  geom_line(
    linewidth = 0.9,
    alpha = 0.75
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    title = "Ast-like BBB Support across NVU states",
    subtitle = "Mean BBB Support UCell score within each biological sample",
    x = "NVU state",
    y = "Mean BBB Support UCell score",
    color = "Sample"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

print(p_Ast_BBB_Support_sample_mean)


# Save plot
ggsave(
  "Vascular subtype analysis/Figures/Astlike_BBB_Support_mean_by_sample_NVU.svg",
  p_Ast_BBB_Support_sample_mean,
  width = 8,
  height = 6
)


# Save sample-level summary
write.csv(
  astlike_BBB_sample_mean,
  "Vascular subtype analysis/Astlike_BBB_Support_mean_by_sample_NVU2.csv",
  row.names = FALSE
)


# Print summary
print(
  astlike_BBB_sample_mean %>%
    select(
      sample,
      condition,
      NVU_State_Label,
      n_cells,
      BBB_Support_mean,
      BBB_Support_median
    ),
  n = Inf
)

cat("\nAnalysis complete.\n")
cat("\nEach point = one biological sample/NVU state.\n")
cat("Each line connects the same sample across NVU states.\n")
cat("Color identifies the biological sample.\n")
cat("\nMean BBB Support scores saved.\n")
cat("SVG plot: Astlike_BBB_Support_mean_by_sample_NVU.svg\n")


#Summary plot-----
# Ast-like BBB Support by condition and NVU state
# Cell-level distributions pooled within condition
# NVU state is the primary biological stratification

library(dplyr)
library(ggplot2)

# Prepare Ast-like cell-level data
astlike_condition_scores <- whole_tissue[[]] %>%
  tibble::rownames_to_column("cell") %>%
  filter(
    new_clusters == "Ast-like",
    NVU_States %in% c(
      "N.T.: Intact",
      "T.I.: Intact",
      "T.I.: Pericyte displaced"
    )
  ) %>%
  select(
    cell,
    sample,
    condition,
    NVU_States,
    BBB_Support_UCell
  ) %>%
  mutate(
    NVU_State_Label = recode(
      NVU_States,
      "N.T.: Intact" = "N.T. Intact",
      "T.I.: Intact" = "T.I. Intact",
      "T.I.: Pericyte displaced" = "T.I. Pericyte displaced"
    ),
    condition = factor(
      recode(
        as.character(condition),
        "Normal" = "Normal",
        "Primary Tumor" = "Gliomas",
        "Secondary Tumor" = "BrM"
      ),
      levels = c(
        "Normal",
        "Gliomas",
        "BrM"
      )
    ),
    NVU_State_Label = factor(
      NVU_State_Label,
      levels = c(
        "N.T. Intact",
        "T.I. Intact",
        "T.I. Pericyte displaced"
      )
    )
  ) %>%
  filter(!is.na(condition))

# Count cells and samples contributing to each group
astlike_condition_counts <- astlike_condition_scores %>%
  group_by(condition, NVU_State_Label) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(sample),
    .groups = "drop"
  )

print(astlike_condition_counts)

# Condition-level summary
astlike_condition_summary <- astlike_condition_scores %>%
  group_by(condition, NVU_State_Label) %>%
  summarise(
    n_cells = n(),
    n_samples = n_distinct(sample),
    BBB_Support_mean = mean(BBB_Support_UCell, na.rm = TRUE),
    BBB_Support_median = median(BBB_Support_UCell, na.rm = TRUE),
    BBB_Support_sd = sd(BBB_Support_UCell, na.rm = TRUE),
    .groups = "drop"
  )

print(astlike_condition_summary)

# Plot
p_Ast_BBB_Support_condition <- ggplot(
  astlike_condition_scores,
  aes(
    x = NVU_State_Label,
    y = BBB_Support_UCell,
    fill = condition
  )
) +
  geom_violin(
    position = position_dodge(width = 0.85),
    scale = "width",
    trim = TRUE,
    alpha = 0.65,
    color = "black"
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.85),
    width = 0.15,
    outlier.shape = NA,
    color = "black"
  ) +
  labs(
    title = "Astrocytic BBB Support across NVU states",
    subtitle = "Astrocyte-level distributions by biological condition",
    x = "NVU state",
    y = "Astrocytic BBB Support Index",
    fill = "Condition"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.title.x = element_text(
      size = 17
    ),
    axis.title.y = element_text(
      size = 17
    ),
    axis.text.x = element_text(
      size = 15,
      angle = 30,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 15
    ),
    plot.title = element_text(
      face = "bold",
      size = 19
    ),
    plot.subtitle = element_text(
      size = 16
    ),
    legend.title = element_text(
      size = 16
    ),
    legend.text = element_text(
      size = 15
    ),
    legend.position = "right"
  )

print(p_Ast_BBB_Support_condition)

# Save figure
ggsave(
  "Vascular subtype analysis/Figures/Astlike_BBB_Support_condition_NVU.svg",
  p_Ast_BBB_Support_condition,
  width = 9,
  height = 5.5,
  device = "svg"
)
# Save data
write.csv(
  astlike_condition_scores,
  "Vascular subtype analysis/Astlike_BBB_Support_condition_NVU_cell_scores.csv",
  row.names = FALSE
)

write.csv(
  astlike_condition_summary,
  "Vascular subtype analysis/Astlike_BBB_Support_condition_NVU_summary.csv",
  row.names = FALSE
)





#EC and Ast Scatterplot----
# EC BBB Stability vs Ast-like BBB Support
# Sample-level analysis across NVU states

library(dplyr)
library(ggplot2)

# Metadata only -- avoids creating large Seurat subsets
whole_tissue_meta <- whole_tissue[[]] %>%
  tibble::rownames_to_column("cell")

# Define the NVU states used in the comparison
nvu_states_scatter <- c(
  "N.T.: Intact",
  "T.I.: Intact",
  "T.I.: Pericyte displaced"
)

# Calculate sample-level EC BBB Stability
ec_sample_scores <- whole_tissue_meta %>%
  filter(
    new_clusters == "EC",
    NVU_States %in% nvu_states_scatter,
    !is.na(BBB_Stability_full)
  ) %>%
  group_by(sample, condition, NVU_States) %>%
  summarise(
    EC_BBB_Stability = mean(BBB_Stability_full, na.rm = TRUE),
    n_EC = n(),
    .groups = "drop"
  )

# Calculate sample-level Ast-like BBB Support
ast_sample_scores <- whole_tissue_meta %>%
  filter(
    new_clusters == "Ast-like",
    NVU_States %in% nvu_states_scatter,
    !is.na(Ast_BBB_Support_Score)
  ) %>%
  group_by(sample, condition, NVU_States) %>%
  summarise(
    Ast_BBB_Support = mean(Ast_BBB_Support_Score, na.rm = TRUE),
    n_Ast = n(),
    .groups = "drop"
  )

# Join EC and Ast-like scores
ec_ast_sample_scores <- ec_sample_scores %>%
  inner_join(
    ast_sample_scores,
    by = c("sample", "condition", "NVU_States")
  ) %>%
  mutate(
    NVU_State_Label = recode(
      NVU_States,
      "N.T.: Intact" = "N.T. Intact",
      "T.I.: Intact" = "T.I. Intact",
      "T.I.: Pericyte displaced" = "T.I. Pericyte displaced"
    ),
    NVU_State_Label = factor(
      NVU_State_Label,
      levels = c(
        "N.T. Intact",
        "T.I. Intact",
        "T.I. Pericyte displaced"
      )
    )
  )

# Inspect the sample-level data
print(ec_ast_sample_scores)

# Count observations by NVU state
print(
  ec_ast_sample_scores %>%
    count(NVU_State_Label)
)

# Scatterplot
p_EC_Ast_BBB_relationship <- ggplot(
  ec_ast_sample_scores,
  aes(
    x = EC_BBB_Stability,
    y = Ast_BBB_Support,
    color = NVU_State_Label
  )
) +
  geom_point(
    size = 3.5,
    alpha = 0.85
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    linewidth = 0.7
  ) +
  labs(
    title = "Coordinated endothelial and astrocytic BBB support",
    x = "EC BBB Stability score",
    y = "Ast-like BBB Support score",
    color = "NVU state"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_EC_Ast_BBB_relationship)

# Save figure
ggsave(
  "Vascular subtype analysis/Figures/EC_vs_Ast_BBB_Support_sample_NVU.svg",
  p_EC_Ast_BBB_relationship,
  width = 6.5,
  height = 5.5,
  device = "svg"
)

# Save sample-level data
write.csv(
  ec_ast_sample_scores,
  "Vascular subtype analysis/EC_vs_Ast_BBB_Support_sample_NVU.csv",
  row.names = FALSE
)



#Final plot----
# EC BBB Stability vs Ast-like BBB Support
# Sample-level observations across NVU states

library(dplyr)
library(ggplot2)

# Rename condition labels for plotting
ec_ast_sample_scores <- ec_ast_sample_scores %>%
  mutate(
    condition_plot = as.character(condition),
    condition_plot = case_when(
      condition_plot == "Normal" ~ "Normal",
      condition_plot == "Primary Tumor" ~ "Gliomas",
      condition_plot == "Secondary Tumor" ~ "BrM",
      TRUE ~ condition_plot
    ),
    condition_plot = factor(
      condition_plot,
      levels = c(
        "Normal",
        "Gliomas",
        "BrM"
      )
    )
  )

# Define consistent NVU state colors
nvu_state_colors <- c(
  "N.T. Intact" = "#4D4D4D",
  "T.I. Intact" = "#E69F00",
  "T.I. Pericyte displaced" = "#56B4E9"
)

# Plot
p_EC_Ast_BBB_relationship <- ggplot(
  ec_ast_sample_scores,
  aes(
    x = EC_BBB_Stability,
    y = Ast_BBB_Support,
    color = NVU_State_Label,
    shape = condition_plot
  )
) +
  geom_point(
    size = 4,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = nvu_state_colors
  ) +
  scale_shape_manual(
    values = c(
      "Normal" = 16,
      "Gliomas" = 17,
      "BrM" = 15
    )
  ) +
  labs(
    title = "Coordinated endothelial and astrocytic BBB support",
    x = "EC BBB Stability score",
    y = "Ast-like BBB Support score",
    color = "NVU state",
    shape = "Condition"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      size = 16
    ),
    axis.text = element_text(
      size = 15
    ),
    legend.title = element_text(
      size = 15
    ),
    legend.text = element_text(
      size = 14
    ),
    legend.position = "right"
  )

print(p_EC_Ast_BBB_relationship)

# Save as SVG
ggsave(
  "Vascular subtype analysis/Figures/EC_vs_Ast_BBB_Support_sample_NVU.svg",
  p_EC_Ast_BBB_relationship,
  width = 8,
  height = 10,
  device = "svg"
)


