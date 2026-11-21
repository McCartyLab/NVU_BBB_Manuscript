

#EC Module & Tumor Module scatterplots----
# EC module vs Tumor module sample-level scatterplots

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)

# Define EC modules
ec_modules <- c(
  "Tight_Junctions",
  "BBB_Transport",
  "Angiogenesis",
  "Vascular_Inflam.",
  "ECM_Remodeling"
)

# Define Tumor modules
tumor_modules <- c(
  "Stem_like",
  "Hypoxia",
  "Mesenchymal",
  "Invasive",
  "Proliferative"
)

# Get EC cells
ec_cells <- whole_tissue@meta.data %>%
  filter(
    new_clusters == "EC",
    !is.na(sample),
    sample != ""
  )

# Get Tumor cells
tumor_cells <- whole_tissue@meta.data %>%
  filter(
    new_clusters == "Tumor",
    !is.na(sample),
    sample != ""
  )

# Calculate sample-level EC module means
ec_sample <- ec_cells %>%
  group_by(sample, condition) %>%
  summarise(
    across(
      all_of(ec_modules),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Calculate sample-level Tumor module means
tumor_sample <- tumor_cells %>%
  group_by(sample, condition) %>%
  summarise(
    across(
      all_of(tumor_modules),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Rename EC module columns
ec_sample <- ec_sample %>%
  rename_with(
    ~ paste0("EC_", .x),
    all_of(ec_modules)
  )

# Rename Tumor module columns
tumor_sample <- tumor_sample %>%
  rename_with(
    ~ paste0("Tumor_", .x),
    all_of(tumor_modules)
  )

# Keep only samples that have BOTH EC and Tumor cells
sample_module_df <- inner_join(
  ec_sample,
  tumor_sample,
  by = "sample",
  suffix = c("_EC", "_Tumor")
)

# Use condition from EC/tumor metadata
sample_module_df$condition <- dplyr::coalesce(
  sample_module_df$condition_EC,
  sample_module_df$condition_Tumor
)

# Remove duplicated condition columns
sample_module_df <- sample_module_df %>%
  select(
    sample,
    condition,
    everything(),
    -condition_EC,
    -condition_Tumor
  )

# Explicitly remove Normal
sample_module_df <- sample_module_df %>%
  filter(sample != "Normal")

# Check samples included
print(
  sample_module_df %>%
    select(sample, condition)
)

cat(
  "\nNumber of samples:",
  nrow(sample_module_df),
  "\n"
)

# Create scatterplots
plot_list <- list()

for (ec_module in ec_modules) {
  
  for (tumor_module in tumor_modules) {
    
    x_col <- paste0("EC_", ec_module)
    y_col <- paste0("Tumor_", tumor_module)
    
    p <- ggplot(
      sample_module_df,
      aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        color = condition
      )
    ) +
      geom_point(
        size = 7,
        alpha = 0.9
      ) +
      geom_text_repel(
        aes(label = sample),
        color = "black",
        size = 8,
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.2,
        segment.color = "grey50",
        segment.linewidth = 0.3
      ) +
      theme_classic() +
      labs(
        title = paste(
          ec_module,
          "vs",
          tumor_module
        ),
        x = ec_module,
        y = tumor_module,
        color = NULL
      ) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 20
        ),
        axis.title = element_text(
          face = "bold",
          size = 22
        ),
        axis.text = element_text(
          color = "black",
          size = 16
        ),
        legend.position = "none",
        legend.text = element_text(
          size = 9
        )
      )
    
    plot_list[[
      paste(
        ec_module,
        tumor_module,
        sep = "_vs_"
      )
    ]] <- p
  }
}

# Combine all 25 plots
EC_vs_Tumor_modules <- wrap_plots(
  plot_list,
  ncol = 5
)

# Display
EC_vs_Tumor_modules


# Save
output_dir <- "Vascular subtype analysis/Figures/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(
  filename = file.path(
    output_dir,
    "EC_vs_Tumor_modules.svg"
  ),
  plot = EC_vs_Tumor_modules,
  width = 24,
  height = 24,
  dpi = 600
)




#Gene Dimplot----
spatialgenesplot <- ImageFeaturePlot(whole_tissue, fov =  c("Lung_fov", "Breast_fov", "Thyroid_fov", "Melanoma_fov", "Kidney_fov"),
                                     features = c("SEMA3A", "F11R", "FN1", "TGFBR2", "VCAM1"),
                                     combine = T,
                                     coord.fixed = T)