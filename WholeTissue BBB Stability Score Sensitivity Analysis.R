#BBB Module Score Sensitivity Analysis


#Sensitivity analysis - BBB Stability Modules----
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)

#Excluded LSR, MFSD2A, SLCO1C1

# Original vs independent BBB Stability Score

df_bbb <- vascular_subtypes_whole_tissue@meta.data %>%
  dplyr::select(
    Original = BBB_Stability_full,
    Independent = BBB_Stability_indep
  ) %>%
  dplyr::filter(
    !is.na(Original),
    !is.na(Independent)
  )

# Number of cells
nrow(df_bbb)

# Pearson correlation
bbb_pearson <- cor.test(
  df_bbb$Original,
  df_bbb$Independent,
  method = "pearson"
)

bbb_pearson

# Spearman correlation
bbb_spearman <- cor.test(
  df_bbb$Original,
  df_bbb$Independent,
  method = "spearman"
)

bbb_spearman

# Correlation coefficient only
cor(
  df_bbb$Original,
  df_bbb$Independent,
  method = "pearson"
)

cor(
  df_bbb$Original,
  df_bbb$Independent,
  method = "spearman"
)

#plot sensitivity analysis

p_bbb_validation <- ggplot(
  df_bbb,
  aes(
    x = Original,
    y = Independent
  )
) +
  geom_point(
    alpha = 0.15,
    size = 0.5
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  labs(
    x = "Original BBB Stability Score",
    y = "Independent BBB Stability Score",
    title = "Sensitivity analysis of BBB Stability Score"
  ) +
  theme_classic()

p_bbb_validation


# Compare modules affected by gene exclusion

module_sensitivity <- list(
  
  "Tight Junctions" = c(
    "Tight_Junctions",
    "Tight_Junctions_indep_UCell"
  ),
  
  "BBB Transport" = c(
    "BBB_Transport",
    "BBB_Transport_indep_UCell"
  )
)

module_sensitivity_results <- purrr::map_dfr(
  names(module_sensitivity),
  
  function(module) {
    
    cols <- module_sensitivity[[module]]
    
    df <- vascular_subtypes_whole_tissue@meta.data %>%
      dplyr::select(
        Original = all_of(cols[1]),
        Independent = all_of(cols[2])
      ) %>%
      dplyr::filter(
        !is.na(Original),
        !is.na(Independent)
      )
    
    pearson <- cor.test(
      df$Original,
      df$Independent,
      method = "pearson"
    )
    
    spearman <- cor.test(
      df$Original,
      df$Independent,
      method = "spearman"
    )
    
    tibble(
      Module = module,
      N = nrow(df),
      Pearson_r = unname(pearson$estimate),
      Pearson_p = pearson$p.value,
      Spearman_rho = unname(spearman$estimate),
      Spearman_p = spearman$p.value
    )
  }
)

module_sensitivity_results


#plot
module_plot_data <- purrr::map_dfr(
  
  names(module_sensitivity),
  
  function(module) {
    
    cols <- module_sensitivity[[module]]
    
    vascular_subtypes_whole_tissue@meta.data %>%
      dplyr::select(
        Original = all_of(cols[1]),
        Independent = all_of(cols[2])
      ) %>%
      dplyr::filter(
        !is.na(Original),
        !is.na(Independent)
      ) %>%
      dplyr::mutate(
        Module = module
      )
  }
)

p_module_sensitivity <- ggplot(
  module_plot_data,
  aes(
    x = Original,
    y = Independent
  )
) +
  geom_point(
    alpha = 0.15,
    size = 0.5
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~Module,
    scales = "free"
  ) +
  labs(
    x = "Original module score",
    y = "Independent module score"
  ) +
  theme_classic()

p_module_sensitivity