# Statistical analysis - Negative ionization mode - Study II
Giovana Domeneghini Mercali

# Load libraries

``` r
library(factoextra) # visualizing PCA results
library(glue) # for easy pasting
library(plotly) # quick interactive plots
library(proxyC) # more efficient large matrices
library(data.table) # easy transposing
library(janitor) # for cleaning names and checking duplicates
library(notame) # for collapsing ions coming from the same metabolite
library(doParallel) # for parallelizing notame specifically
library(patchwork) # for making multi-panel plots
library(rstatix) # for additional univariate functionality
library(readxl)
library(svglite)
library(scales)  # para função `seq_gradient_pal()
library(openxlsx)
library(pheatmap) # for heínsatmaps

# this is at the end hoping that the default select will be that from dplyr
library(tidyverse) # for everything
```

# Read in data

Negative mode.

``` r
# read in metabolomics data
feature_table_long_log2 <- read_csv("Processed_data/feature_table_long_log2_neg.csv") 
```

    Rows: 52104 Columns: 7
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (4): sample_name, treatment, sample_or_qc, mz_rt
    dbl (3): run_order, 85.0285_0.3638, rel_abund

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
feature_table_wide_log2 <- read_csv("Processed_data/feature_table_wide_log2_neg.csv") 
```

    Rows: 24 Columns: 2176
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr    (3): sample_name, treatment, sample_or_qc
    dbl (2173): run_order, 85.0285_0.3638, 87.0077_0.3573, 96.9591_0.3664, 111.0...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# read in meta data - names of the experiments were changed to be more organized.
metadata <- read_csv("Data/metadata.csv",
                  trim_ws = TRUE)
```

    Rows: 24 Columns: 6
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (3): sample_name, short_sample_name, treatment
    dbl (3): run_order, soluble_solids, mass

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# read in meta data
annotation <- read_csv("Data/Annotation_neg.csv", 
                  locale=locale(encoding = "latin1"))
```

    Rows: 38 Columns: 2
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (2): Compound, mz_rt

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

# Boxplots

``` r
feature_table_long_log2_final <- feature_table_long_log2
```

``` r
boxplot_log2_transformed <- feature_table_long_log2_final %>%
  ggplot(aes(x = sample_name, y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_classic(base_size = 12) +
   theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)  # Adiciona borda completa
  ) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("UJ" = "#6A3D9A",
                                "P76" = "#2E86C1",
                                "P85" = "#EB984E",
                                "P95" = "#42846B",
                                "QC" = "#B2BABB")) +
  labs(y = "Relative abundance", x = "Sample", 
  title = "LC-MS (+) feature abundances by Sample",
       subtitle = "Data is log2 transformed")


boxplot_log2_transformed
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-4-1.png)

We can also look at this data by run order to see if we have any overall
run order effects visible.

``` r
boxplot_log2_transformed_run_order <-feature_table_long_log2_final %>%
  mutate(sample_name = fct_reorder(sample_name, run_order)) %>%
  ggplot(aes(x = sample_name , y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_classic(base_size = 12) +
   theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)  # Adiciona borda completa
  ) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("UJ" = "#6A3D9A",
                                "P76" = "#2E86C1",
                                "P85" = "#EB984E",
                                "P95" = "#42846B",
                                "QC" = "#B2BABB")) +
  labs(y = "Relative abundance", x = "Sample")

boxplot_log2_transformed_run_order
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-5-1.png)

``` r
boxplot <- boxplot_log2_transformed / boxplot_log2_transformed_run_order + plot_annotation(tag_levels = c("A"))

ggsave(plot = boxplot,
       filename = "Plots/boxplot_final_neg.jpeg", width = 19, height = 19, units = "cm")

ggsave(plot = boxplot,
       filename = "Plots/boxplot_final_neg.svg")
```

    Saving 7 x 5 in image

# PCAs

## With QCs

``` r
feature_table_wide_log2_final <- feature_table_long_log2_final %>%
  pivot_wider(names_from = mz_rt,
              values_from = rel_abund)
```

``` r
pca_qc <- prcomp(feature_table_wide_log2_final[,-c(1:4)], # remove metadata
                 scale = FALSE, # we did our own scaling
                 center = TRUE) # true is the default

summary(pca_qc)
```

    Importance of components:
                              PC1     PC2     PC3     PC4     PC5     PC6     PC7
    Standard deviation     3.0127 2.62045 2.40670 2.22571 2.14204 2.03224 1.97492
    Proportion of Variance 0.1107 0.08374 0.07063 0.06041 0.05595 0.05036 0.04756
    Cumulative Proportion  0.1107 0.19442 0.26505 0.32546 0.38141 0.43178 0.47934
                               PC8     PC9   PC10    PC11    PC12    PC13    PC14
    Standard deviation     1.91990 1.86510 1.8337 1.81722 1.78212 1.74800 1.73934
    Proportion of Variance 0.04495 0.04242 0.0410 0.04027 0.03873 0.03726 0.03689
    Cumulative Proportion  0.52429 0.56671 0.6077 0.64798 0.68671 0.72397 0.76086
                              PC15    PC16    PC17    PC18    PC19    PC20   PC21
    Standard deviation     1.71367 1.69093 1.67376 1.63816 1.63277 1.57942 1.4991
    Proportion of Variance 0.03581 0.03487 0.03416 0.03272 0.03251 0.03042 0.0274
    Cumulative Proportion  0.79667 0.83154 0.86570 0.89843 0.93094 0.96136 0.9888
                              PC22      PC23      PC24
    Standard deviation     0.96013 9.257e-05 2.087e-14
    Proportion of Variance 0.01124 0.000e+00 0.000e+00
    Cumulative Proportion  1.00000 1.000e+00 1.000e+00

Look at how much variance is explained by each PC.

``` r
importance_qc <- summary(pca_qc)$importance %>%
  as.data.frame()

head(importance_qc)
```

                                PC1      PC2     PC3      PC4      PC5      PC6
    Standard deviation     3.012716 2.620447 2.40670 2.225709 2.142045 2.032242
    Proportion of Variance 0.110680 0.083740 0.07063 0.060410 0.055950 0.050360
    Cumulative Proportion  0.110680 0.194420 0.26505 0.325460 0.381410 0.431780
                                PC7      PC8      PC9     PC10    PC11     PC12
    Standard deviation     1.974924 1.919902 1.865096 1.833717 1.81722 1.782118
    Proportion of Variance 0.047560 0.044950 0.042420 0.041000 0.04027 0.038730
    Cumulative Proportion  0.479340 0.524290 0.566710 0.607710 0.64798 0.686710
                               PC13     PC14     PC15     PC16     PC17     PC18
    Standard deviation     1.748001 1.739336 1.713666 1.690934 1.673758 1.638163
    Proportion of Variance 0.037260 0.036890 0.035810 0.034870 0.034160 0.032720
    Cumulative Proportion  0.723970 0.760860 0.796670 0.831540 0.865700 0.898430
                               PC19     PC20     PC21      PC22         PC23
    Standard deviation     1.632768 1.579415 1.499065 0.9601325 9.257306e-05
    Proportion of Variance 0.032510 0.030420 0.027400 0.0112400 0.000000e+00
    Cumulative Proportion  0.930940 0.961360 0.988760 1.0000000 1.000000e+00
                                   PC24
    Standard deviation     2.087285e-14
    Proportion of Variance 0.000000e+00
    Cumulative Proportion  1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_qc)
```

    Warning in geom_bar(stat = "identity", fill = barfill, color = barcolor, :
    Ignoring empty aesthetic: `width`.

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-10-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_qc)
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-11-1.png)

Make a scores plot but prettier.

``` r
# create a df of pca_qc$x
scores_UJ_qc <- as.data.frame(pca_qc$x)

# bind meta-data
scores_qc <- bind_cols(feature_table_wide_log2_final[,1:4], # first 4 columns
                       scores_UJ_qc)
```

Plot.

``` r
# create objects indicating percent variance explained by PC1 and PC2
PC1_percent_qc <- round((importance_qc[2,1])*100, # index 2nd row, 1st column, times 100
                         1) # round to 1 decimal
PC2_percent_qc <- round((importance_qc[2,2])*100, 1) 

# plot
# aes(text) is for setting tooltip with plotly later to indicate hover text
(scores_qc_plot <- scores_qc %>%
  ggplot(aes(x = PC1, y = PC2, fill = treatment, text = glue("Sample: {sample_name},
                                                           Treatment: {treatment}"))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(shape = 21, size = 3, color = "black") +
  theme_classic(base_size = 12) +
   theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)  # Adiciona borda completa
  ) +
    scale_fill_manual(values = c("UJ" = "#6A3D9A",
                                "P76" = "#2E86C1",
                                "P85" = "#EB984E",
                                "P95" = "#42846B",
                                "QC" = "#B2BABB")) +
  labs(x = glue("PC1: {PC1_percent_qc}%"), 
       y = glue("PC2: {PC2_percent_qc}%"), 
       title = "", tag = "A")) + 
  coord_cartesian(xlim = c(-20, 20), ylim = c(-10, 10))
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-13-1.png)

``` r
ggsave(plot = scores_qc_plot,
       filename = "Plots/pca_qc_final_neg.jpeg", width = 140, height = 120, units = "mm")
```

Make a loadings plot (points are features) even though it might not be
that useful.

``` r
fviz_pca_var(pca_qc)
```

    Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ℹ Please use `linewidth` instead.
    ℹ The deprecated feature was likely used in the ggpubr package.
      Please report the issue at <https://github.com/kassambara/ggpubr/issues>.

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-15-1.png)

See what I mean? Not that useful. There are some functions in PCAtools
that label only the points that most contribute to each PC. Could also
do this manually if its of interest.

## Without QCs

``` r
feature_table_wide_log2_final_noqc <- feature_table_wide_log2_final %>%
  filter(sample_or_qc == "Sample")


pca_noqc <- prcomp(feature_table_wide_log2_final_noqc[,-c(1:4)], # remove metadata
                 scale = FALSE, # we did our own scaling
                 center = TRUE) # true is the default

summary(pca_noqc)
```

    Importance of components:
                              PC1     PC2     PC3     PC4     PC5     PC6     PC7
    Standard deviation     3.7180 3.00273 2.76372 2.64469 2.53689 2.43080 2.39067
    Proportion of Variance 0.1485 0.09684 0.08204 0.07513 0.06913 0.06347 0.06139
    Cumulative Proportion  0.1485 0.24532 0.32736 0.40249 0.47162 0.53508 0.59647
                               PC8     PC9    PC10    PC11    PC12    PC13    PC14
    Standard deviation     2.34619 2.26675 2.20977 2.15906 2.14990 2.09665 2.07790
    Proportion of Variance 0.05912 0.05519 0.05245 0.05007 0.04965 0.04722 0.04638
    Cumulative Proportion  0.65559 0.71078 0.76323 0.81330 0.86294 0.91016 0.95653
                              PC15      PC16
    Standard deviation     2.01166 2.149e-14
    Proportion of Variance 0.04347 0.000e+00
    Cumulative Proportion  1.00000 1.000e+00

Look at how much variance is explained by each PC.

``` r
importance_noqc <- summary(pca_noqc)$importance %>%
  as.data.frame()

head(importance_noqc)
```

                                PC1      PC2      PC3      PC4      PC5      PC6
    Standard deviation     3.718038 3.002728 2.763719 2.644688 2.536892 2.430802
    Proportion of Variance 0.148480 0.096840 0.082040 0.075130 0.069130 0.063470
    Cumulative Proportion  0.148480 0.245320 0.327360 0.402490 0.471620 0.535080
                                PC7      PC8      PC9     PC10     PC11     PC12
    Standard deviation     2.390674 2.346185 2.266745 2.209767 2.159063 2.149905
    Proportion of Variance 0.061390 0.059120 0.055190 0.052450 0.050070 0.049650
    Cumulative Proportion  0.596470 0.655590 0.710780 0.763230 0.813300 0.862940
                              PC13    PC14    PC15         PC16
    Standard deviation     2.09665 2.07790 2.01166 2.148864e-14
    Proportion of Variance 0.04722 0.04638 0.04347 0.000000e+00
    Cumulative Proportion  0.91016 0.95653 1.00000 1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_noqc)
```

    Warning in geom_bar(stat = "identity", fill = barfill, color = barcolor, :
    Ignoring empty aesthetic: `width`.

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-18-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_noqc)
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-19-1.png)

Make a scores plot but prettier.

``` r
# create a df of pca_qc$x
scores_UJ_noqc <- as.data.frame(pca_noqc$x)

# bind meta-data
scores_noqc <- bind_cols(feature_table_wide_log2_final_noqc[,1:4], # metadata
                         scores_UJ_noqc)
```

Plot.

``` r
# create objects indicating percent variance explained by PC1 and PC2
PC1_percent_noqc <- round((importance_noqc[2,1])*100, # index 2nd row, 1st column, times 100
                         1) # round to 1 decimal
PC2_percent_noqc <- round((importance_noqc[2,2])*100, 1) 

# plot
# aes(text) is for setting tooltip with plotly later to indicate hover text
(scores_noqc_plot <- scores_noqc %>%
  ggplot(aes(x = PC1, y = PC2, fill = treatment, text = glue("Sample: {sample_name},
                                                             Treatment: {treatment}"))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(shape = 21, size = 3, color = "black") +
  theme_classic(base_size = 12) +
   theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)  # Adiciona borda completa
  ) +
    scale_fill_manual(values = c("UJ" = "#6A3D9A",
                                "P76" = "#2E86C1",
                                "P85" = "#EB984E",
                                "P95" = "#42846B")) +
  coord_cartesian(xlim = c(-20, 20), ylim = c(-10, 10)) +
  labs(x = glue("PC1: {PC1_percent_noqc}%"), 
       y = glue("PC2: {PC2_percent_noqc}%"), 
       title = "", tag = "A"))
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-21-1.png)

``` r
ggsave(plot = scores_noqc_plot,
       filename = "Plots/pca_noqc_final_neg.jpeg", width = 140, height = 120, units = "mm")
```

Make a loadings plot (points are features) even though it might not be
that useful.

``` r
fviz_pca_var(pca_noqc)
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-23-1.png)

See what I mean? Not that useful. There are some functions in PCAtools
that label only the points that most contribute to each PC. Could also
do this manually if its of interest.

# Add annotations

``` r
feature_table_long_log2_final <- left_join(feature_table_long_log2_final, annotation, by = "mz_rt")
  
unique(feature_table_long_log2_final$Compound)
```

    [1] NA

``` r
feature_table_long_log2_final <- feature_table_long_log2_final %>%
  mutate(id = if_else(is.na(Compound), as.character(mz_rt), as.character(Compound)))
```

``` r
write_csv(feature_table_long_log2_final,
          file = "Processed_data/feature_table_long_log2_final_annotated_neg.csv")
```

# Heatmap of all features

``` r
# find which features have no variance across samples
# having samples with no variance won't work for a heatmap
all_comparisons_novariance <- feature_table_long_log2_final %>%
  filter(sample_or_qc == "Sample") %>% # samples only
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)

# what features have no variance?
all_comparisons_novariance$id
```

    character(0)

``` r
# wrangle
feature_table_for_heatmap <- feature_table_wide_log2_final %>%
  filter(sample_or_qc == "Sample") %>% # samples only
  select(-run_order, -sample_or_qc) %>% # remove columns we don't need
  column_to_rownames(var = "sample_name") # move the sample names to rownames for plotting
```

``` r
heatmap <- pheatmap(feature_table_for_heatmap[,-1], # remove treatment column
         cluster_rows = TRUE, # HCA on samples
         cluster_cols = TRUE, # HCA on features 
         scale = "column",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         legend_breaks = c(-4, 0, 4),
         legend_labels = c('Low concentration', 'Medium concentration', 'High concentration'),
         border_color = F,
         show_colnames = FALSE, # scale for each feature
         filename = "Plots/heatmap_neg.png")
```

## ANOVA across groups

``` r
anova <- feature_table_long_log2_final %>%
  filter(!sample_or_qc == "QC") %>% # remove QCs
  dplyr::select(sample_name, treatment, mz_rt, rel_abund) %>%
  group_by(mz_rt) %>%
 # Remove grupos onde todos os valores de rel_abund são iguais
  filter(n_distinct(rel_abund) > 1) %>%
  anova_test(rel_abund ~ treatment)

# adjust pvalues for multiple testing
anova_padjusted <- p.adjust(anova$p, method = "BH") %>%
  as.data.frame() %>%
  rename(p_adj = 1)

anova_padj <- bind_cols(as.data.frame(anova), anova_padjusted) %>%
  rename(padj = 9)

# extract out only the significantly different features
anova_sig <- anova_padj %>%
  as.data.frame() %>%
  filter(padj <= 0.05)

# how many features are significantly different between the groups?
nrow(anova_sig)
```

    [1] 9

# Heatmap of features that are different between groups

``` r
feature_table_long_log2_final_anova <- inner_join(feature_table_long_log2_final, anova_sig, by = "mz_rt") %>%
   select(sample_name, treatment, sample_or_qc, mz_rt, rel_abund, id)
  
# unique(feature_table_long_log2_final$id)
```

``` r
# find which features have no variance across samples
# having samples with no variance won't work for a heatmap
all_comparisons_novariance <- feature_table_long_log2_final_anova %>%
  filter(sample_or_qc == "Sample") %>% # samples only
  dplyr::select(sample_name, treatment, mz_rt, rel_abund, id) %>%
  group_by(id) %>%
  summarize(stdev = sd(rel_abund),
            mean = mean(rel_abund)) %>%
  filter(stdev == 0)

# what features have no variance?
all_comparisons_novariance$id
```

    character(0)

``` r
# wrangle

feature_table_wide_log2_final_anova <- feature_table_long_log2_final_anova %>%
  select(-"id")  %>%
  pivot_wider(names_from = "mz_rt",
              values_from = "rel_abund")


feature_table_anova_for_heatmap <- feature_table_wide_log2_final_anova %>%
  filter(sample_or_qc == "Sample") %>% # samples only
  select(-sample_or_qc) %>% # remove columns we don't need
  column_to_rownames(var = "sample_name") # move the sample names to rownames for plotting


heatmap <- pheatmap(feature_table_anova_for_heatmap[,-1], # remove treatment column
         cluster_rows = TRUE, # HCA on samples
         cluster_cols = TRUE, # HCA on features 
         scale = "column",
         show_colnames = FALSE, # scale for each feature
         legend_breaks = c(-4, 0, 4),
         legend_labels = c('Low concentration', 'Medium concentration', 'High concentration'),
         border_color = F,
         filename = "Plots/heatmap_anova_neg.png")
```

``` r
ncol(feature_table_for_heatmap)
```

    [1] 2173

``` r
ncol(feature_table_anova_for_heatmap)
```

    [1] 10

# Volcano plots

``` r
volcano_function <- function(data, group1, group2, 
                             padj_threshold = 0.05, 
                             fc_threshold = 1.0) {
  
  # Identify features with zero variance
  no_variance <- data %>%
    filter(treatment %in% c(group1, group2)) %>%
    select(sample_name, treatment, mz_rt, rel_abund) %>%
    group_by(mz_rt) %>%
    summarize(
      stdev = sd(rel_abund, na.rm = TRUE),
      mean = mean(rel_abund, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(stdev == 0)
  
  features_to_remove <- no_variance$mz_rt
  
  
  # Perform t-test for each feature
  tests <- data %>%
    filter(treatment %in% c(group1, group2)) %>%
    filter(!mz_rt %in% features_to_remove) %>%
    select(sample_name, treatment, mz_rt, rel_abund) %>%
    group_by(mz_rt) %>%
    t_test(
      rel_abund ~ treatment,
      paired = FALSE,
      detailed = TRUE
    ) %>%
    ungroup() %>%
    
    # Benjamini-Hochberg correction across all features
    mutate(
      p.adj = p.adjust(p, method = "BH")
    )
  
  
  # Calculate log2 fold-change
  fc <- data %>%
    filter(treatment %in% c(group1, group2)) %>%
    group_by(treatment, mz_rt) %>%
    summarize(
      mean = mean(rel_abund, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = treatment,
      values_from = mean
    ) %>%
    mutate(
      log2_FC = .data[[group1]] - .data[[group2]]
    )
  
  
  # Merge fold-change and statistical results
  results <- right_join(fc, tests, by = "mz_rt") %>%
    mutate(
      neglog10padj = -log10(p.adj)
    )
  
  
  # Significant features higher in group 1
  higher_in_group1 <- results %>%
    filter(
      p.adj <= padj_threshold,
      log2_FC >= fc_threshold
    )
  
  
  # Significant features higher in group 2
  higher_in_group2 <- results %>%
    filter(
      p.adj <= padj_threshold,
      log2_FC <= -fc_threshold
    )
  
  
  # Messages
  message(
    "Total significant features (FDR <= ", 
    padj_threshold, "): ",
    nrow(results %>% filter(p.adj <= padj_threshold))
  )
  
  message(
    "Features with >=2-fold higher in ", group1,
    " (FDR <= ", padj_threshold, "): ",
    nrow(higher_in_group1)
  )
  
  message(
    "Features with >=2-fold higher in ", group2,
    " (FDR <= ", padj_threshold, "): ",
    nrow(higher_in_group2)
  )
  
  
  return(
    list(
      full_results = results,
      significant_group1 = higher_in_group1,
      significant_group2 = higher_in_group2
    )
  )
}
```

## UJ x P76

``` r
# Run function
volcano_UJ_P76 <- volcano_function(data = feature_table_long_log2_final, group1 = "UJ", group2 = "P76")
```

    Total significant features (FDR <= 0.05): 0

    Features with >=2-fold higher in UJ (FDR <= 0.05): 0

    Features with >=2-fold higher in P76 (FDR <= 0.05): 0

Se o menor p.adj da análise é 0.3748, isso significa que nenhuma feature
sobreviveu ao FDR. Isso normalmente acontece por um destes motivos: -
Poucas réplicas biológicas (baixa potência estatística); - Muitas
features (2171 testes aumentam a penalização); - Diferenças pequenas ou
variabilidade alta.

``` r
volcano_UJ_P76$full_results %>%
  summarise(
    min_p = min(p, na.rm = TRUE),
    min_padj = min(p.adj, na.rm = TRUE),
    n_sig = sum(p.adj < 0.05, na.rm = TRUE)
  )
```

    # A tibble: 1 × 3
         min_p min_padj n_sig
         <dbl>    <dbl> <int>
    1 0.000322    0.376     0

``` r
volcano_UJ_P76$full_results %>%
  select(mz_rt, p, p.adj, log2_FC) %>%
  arrange(p) %>%
  slice_head(n = 20)
```

    # A tibble: 20 × 4
       mz_rt                  p p.adj log2_FC
       <chr>              <dbl> <dbl>   <dbl>
     1 352.2618_4.8851 0.000322 0.376 -0.525 
     2 338.0956_2.7311 0.000395 0.376  0.149 
     3 615.2445_3.7984 0.000695 0.376 -0.585 
     4 695.1885_0.3253 0.000888 0.376 -0.362 
     5 453.2811_5.8416 0.000959 0.376 -0.287 
     6 281.1759_4.4884 0.00104  0.376 -0.311 
     7 375.148_7.058   0.00181  0.543 -0.453 
     8 311.0495_1.0681 0.00231  0.543  0.323 
     9 301.0355_3.4857 0.00248  0.543 -0.200 
    10 431.22_6.6004   0.00409  0.543 -0.288 
    11 478.293_5.9141  0.00431  0.543 -0.222 
    12 816.038_0.4991  0.00476  0.543 -0.239 
    13 530.3005_5.8866 0.0049   0.543 -0.381 
    14 607.1657_3.0708 0.00539  0.543 -0.269 
    15 383.047_0.3382  0.00562  0.543 -0.157 
    16 593.2723_5.8319 0.00593  0.543 -0.299 
    17 463.0878_2.8776 0.00604  0.543 -0.0630
    18 310.9692_0.3992 0.0063   0.543 -0.123 
    19 475.1966_4.4596 0.00663  0.543 -0.320 
    20 452.2775_5.8374 0.00668  0.543 -0.304 

``` r
volcano_UJ_P76$full_results %>%
  select(mz_rt, p, p.adj, log2_FC, neglog10padj) %>%
  arrange(p.adj)
```

    # A tibble: 2,171 × 5
       mz_rt                  p p.adj log2_FC neglog10padj
       <chr>              <dbl> <dbl>   <dbl>        <dbl>
     1 281.1759_4.4884 0.00104  0.376  -0.311        0.424
     2 338.0956_2.7311 0.000395 0.376   0.149        0.424
     3 352.2618_4.8851 0.000322 0.376  -0.525        0.424
     4 453.2811_5.8416 0.000959 0.376  -0.287        0.424
     5 615.2445_3.7984 0.000695 0.376  -0.585        0.424
     6 695.1885_0.3253 0.000888 0.376  -0.362        0.424
     7 133.0142_1.8948 0.0185   0.543   0.199        0.265
     8 167.0352_2.8411 0.0107   0.543   0.450        0.265
     9 213.9771_3.0912 0.00953  0.543  -0.465        0.265
    10 301.0355_3.4857 0.00248  0.543  -0.200        0.265
    # ℹ 2,161 more rows

``` r
(UJ_P76_volcano <- ggplot(
  volcano_UJ_P76$full_results,
  aes(x = log2_FC, y = neglog10padj, text = mz_rt)
) +
  
  # All features
  geom_point() +
  
  # Significant features higher in UJ
  geom_point(
    data = volcano_UJ_P76$significant_group1,
    aes(x = log2_FC, y = neglog10padj),
    color = "tomato3"
  ) +
  
  # Significant features higher in P76
  geom_point(
    data = volcano_UJ_P76$significant_group2,
    aes(x = log2_FC, y = neglog10padj),
    color = "royalblue4"
  ) +
  
  # Threshold lines
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = -1,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  
  # X-axis limits
  coord_cartesian(
    xlim = c(-2, 2)
  ) +
  
  # Same formatting as Code 2
  theme_classic(base_size = 8) +
  
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  ) +
  
  # Labels
  labs(
    x = expression(Log[2]~"(fold change)"),
    y = expression(-Log[10]~"(adjusted p-value)"),
    title = "UJ vs P76",
    tag = "A"
  )
)
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-41-1.png)

Saving volcano plot

``` r
ggsave(plot = UJ_P76_volcano,
       filename = "Plots/UJ_P76_volcano_neg.svg", width = 190, height = 110, units = "mm")
```

## UJ x P85

``` r
# Run function
volcano_UJ_P85 <- volcano_function(
  data = feature_table_long_log2_final,
  group1 = "UJ",
  group2 = "P85"
)
```

    Total significant features (FDR <= 0.05): 0

    Features with >=2-fold higher in UJ (FDR <= 0.05): 0

    Features with >=2-fold higher in P85 (FDR <= 0.05): 0

``` r
volcano_UJ_P85$full_results %>%
  select(mz_rt, p, p.adj, log2_FC, neglog10padj) %>%
  arrange(p.adj)
```

    # A tibble: 2,171 × 5
       mz_rt                   p p.adj log2_FC neglog10padj
       <chr>               <dbl> <dbl>   <dbl>        <dbl>
     1 329.1524_4.7039 0.0000586 0.127  -0.203        0.895
     2 510.0885_3.0601 0.000127  0.138  -0.250        0.861
     3 298.1227_1.8221 0.000582  0.316  -0.273        0.500
     4 375.148_7.058   0.000472  0.316  -0.557        0.500
     5 384.0503_0.3557 0.00107   0.332   0.248        0.479
     6 401.0884_2.6237 0.00103   0.332  -0.954        0.479
     7 556.2865_7.8875 0.000937  0.332  -0.593        0.479
     8 301.0355_3.4857 0.00162   0.352  -0.337        0.454
     9 304.0827_2.8256 0.00139   0.352   0.351        0.454
    10 370.07_2.9495   0.00154   0.352   0.634        0.454
    # ℹ 2,161 more rows

``` r
(UJ_P85_volcano <- ggplot(
  volcano_UJ_P85$full_results,
  aes(x = log2_FC, y = neglog10padj, text = mz_rt)
) +
  
  # All features
  geom_point() +
  
  # Significant features higher in UJ
  geom_point(
    data = volcano_UJ_P85$significant_group1,
    aes(x = log2_FC, y = neglog10padj),
    color = "tomato3"
  ) +
  
  # Significant features higher in P85
  geom_point(
    data = volcano_UJ_P85$significant_group2,
    aes(x = log2_FC, y = neglog10padj),
    color = "royalblue4"
  ) +
  
  # Threshold lines
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = -1,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  
  # X-axis limits
  coord_cartesian(
    xlim = c(-2, 2)
  ) +
  
  # Formatting
  theme_classic(base_size = 8) +
  
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  ) +
  
  # Labels
  labs(
    x = expression(Log[2]~"(fold change)"),
    y = expression(-Log[10]~"(adjusted p-value)"),
    title = "UJ vs P85",
    tag = "B"
  )
)
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-45-1.png)

``` r
ggsave(
  plot = UJ_P85_volcano,
  filename = "Plots/UJ_P85_volcano_neg.svg",
  width = 190,
  height = 110,
  units = "mm"
)
```

## UJ x P95

``` r
# Run function
volcano_UJ_P95 <- volcano_function(
  data = feature_table_long_log2_final,
  group1 = "UJ",
  group2 = "P95"
)
```

    Total significant features (FDR <= 0.05): 2

    Features with >=2-fold higher in UJ (FDR <= 0.05): 2

    Features with >=2-fold higher in P95 (FDR <= 0.05): 0

``` r
# View results
volcano_UJ_P95$full_results %>%
  select(mz_rt, p, p.adj, log2_FC, neglog10padj) %>%
  arrange(p.adj)
```

    # A tibble: 2,171 × 5
       mz_rt                   p  p.adj log2_FC neglog10padj
       <chr>               <dbl>  <dbl>   <dbl>        <dbl>
     1 370.07_2.9495   0.0000425 0.0461   1.33         1.34 
     2 380.0988_2.9488 0.0000386 0.0461   1.41         1.34 
     3 372.0679_2.9462 0.000113  0.0727   2.05         1.14 
     4 789.2977_3.6511 0.000134  0.0727  -0.548        1.14 
     5 304.0827_2.8256 0.000213  0.0925   0.588        1.03 
     6 340.0594_2.8314 0.000617  0.167    0.519        0.776
     7 352.2439_4.7097 0.000502  0.167   -0.409        0.776
     8 615.2445_3.7984 0.000568  0.167   -0.523        0.776
     9 353.1449_2.0239 0.00108   0.261    0.237        0.584
    10 593.1505_2.8414 0.00129   0.280   -0.221        0.553
    # ℹ 2,161 more rows

``` r
# Volcano plot

(UJ_P95_volcano <- ggplot(
  volcano_UJ_P95$full_results,
  aes(x = log2_FC, y = neglog10padj, text = mz_rt)
) +
  
  # All features
  geom_point() +
  
  # Significant features higher in UJ
  geom_point(
    data = volcano_UJ_P95$significant_group1,
    aes(x = log2_FC, y = neglog10padj),
    color = "tomato3"
  ) +
  
  # Significant features higher in P95
  geom_point(
    data = volcano_UJ_P95$significant_group2,
    aes(x = log2_FC, y = neglog10padj),
    color = "royalblue4"
  ) +
  
  # Threshold lines
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = -1,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  
  # X-axis limits
  coord_cartesian(
    xlim = c(-2, 2)
  ) +
  
  # Formatting
  theme_classic(base_size = 8) +
  
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  ) +
  
  # Labels
  labs(
    x = expression(Log[2]~"(fold change)"),
    y = expression(-Log[10]~"(adjusted p-value)"),
    title = "UJ vs P95",
    tag = "C"
  )
)
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-49-1.png)

``` r
ggsave(
  plot = UJ_P95_volcano,
  filename = "Plots/UJ_P95_volcano_neg.svg",
  width = 190,
  height = 110,
  units = "mm"
)
```

## P76 x P85

``` r
# Run function
volcano_P76_P85 <- volcano_function(
  data = feature_table_long_log2_final,
  group1 = "P76",
  group2 = "P85"
)
```

    Total significant features (FDR <= 0.05): 0

    Features with >=2-fold higher in P76 (FDR <= 0.05): 0

    Features with >=2-fold higher in P85 (FDR <= 0.05): 0

``` r
# View results
volcano_P76_P85$full_results %>%
  select(mz_rt, p, p.adj, log2_FC, neglog10padj) %>%
  arrange(p.adj)
```

    # A tibble: 2,171 × 5
       mz_rt                p p.adj log2_FC neglog10padj
       <chr>            <dbl> <dbl>   <dbl>        <dbl>
     1 1007.6805_4.2466 0.705 0.998  0.103      0.000891
     2 1013.2329_2.4503 0.694 0.998  0.0629     0.000891
     3 1015.2316_2.4787 0.521 0.998  0.175      0.000891
     4 1025.0717_0.6768 0.609 0.998  0.0784     0.000891
     5 1025.0724_0.5695 0.664 0.998  0.127      0.000891
     6 1025.5761_0.6758 0.474 0.998 -0.0665     0.000891
     7 1026.0754_0.6632 0.495 0.998  0.0443     0.000891
     8 1029.0441_0.4987 0.328 0.998  0.151      0.000891
     9 1033.0587_0.6602 0.651 0.998 -0.0488     0.000891
    10 1036.0645_0.7069 0.655 0.998  0.0365     0.000891
    # ℹ 2,161 more rows

``` r
# Volcano plot

(P76_P85_volcano <- ggplot(
  volcano_P76_P85$full_results,
  aes(x = log2_FC, y = neglog10padj, text = mz_rt)
) +
  
  # All features
  geom_point() +
  
  # Significant features higher in P76
  geom_point(
    data = volcano_P76_P85$significant_group1,
    aes(x = log2_FC, y = neglog10padj),
    color = "tomato3"
  ) +
  
  # Significant features higher in P85
  geom_point(
    data = volcano_P76_P85$significant_group2,
    aes(x = log2_FC, y = neglog10padj),
    color = "royalblue4"
  ) +
  
  # Threshold lines
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = -1,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  
  # X-axis limits
  coord_cartesian(
    xlim = c(-2, 2)
  ) +
  
  # Formatting
  theme_classic(base_size = 8) +
  
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  ) +
  
  # Labels
  labs(
    x = expression(Log[2]~"(fold change)"),
    y = expression(-Log[10]~"(adjusted p-value)"),
    title = "P76 vs P85"
  )
)
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-53-1.png)

``` r
ggsave(
  plot = P76_P85_volcano,
  filename = "Plots/P76_P85_volcano_neg.svg",
  width = 190,
  height = 110,
  units = "mm"
)
```

## P85 x P95

``` r
# Run function
volcano_P85_P95 <- volcano_function(
  data = feature_table_long_log2_final,
  group1 = "P85",
  group2 = "P95"
)
```

    Total significant features (FDR <= 0.05): 0

    Features with >=2-fold higher in P85 (FDR <= 0.05): 0

    Features with >=2-fold higher in P95 (FDR <= 0.05): 0

``` r
# View results
volcano_P85_P95$full_results %>%
  select(mz_rt, p, p.adj, log2_FC, neglog10padj) %>%
  arrange(p.adj)
```

    # A tibble: 2,171 × 5
       mz_rt                  p p.adj log2_FC neglog10padj
       <chr>              <dbl> <dbl>   <dbl>        <dbl>
     1 353.1449_2.0239 0.000132 0.287   0.269        0.543
     2 340.0594_2.8314 0.00118  0.497   0.359        0.304
     3 370.07_2.9495   0.00142  0.497   0.692        0.304
     4 372.0679_2.9462 0.00141  0.497   1.36         0.304
     5 380.0988_2.9488 0.00183  0.497   0.846        0.304
     6 517.077_2.583   0.0018   0.497  -0.385        0.304
     7 523.3187_5.1789 0.00096  0.497   0.476        0.304
     8 642.1671_2.9987 0.000526 0.497  -0.144        0.304
     9 659.1355_2.605  0.00226  0.545   0.684        0.263
    10 593.2723_5.8319 0.00275  0.597   0.467        0.224
    # ℹ 2,161 more rows

``` r
# Volcano plot

(P85_P95_volcano <- ggplot(
  volcano_P85_P95$full_results,
  aes(x = log2_FC, y = neglog10padj, text = mz_rt)
) +
  
  # All features
  geom_point() +
  
  # Significant features higher in P85
  geom_point(
    data = volcano_P85_P95$significant_group1,
    aes(x = log2_FC, y = neglog10padj),
    color = "tomato3"
  ) +
  
  # Significant features higher in P95
  geom_point(
    data = volcano_P85_P95$significant_group2,
    aes(x = log2_FC, y = neglog10padj),
    color = "royalblue4"
  ) +
  
  # Threshold lines
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = -1,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  
  # X-axis limits
  coord_cartesian(
    xlim = c(-2, 2)
  ) +
  
  # Formatting
  theme_classic(base_size = 8) +
  
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    )
  ) +
  
  # Labels
  labs(
    x = expression(Log[2]~"(fold change)"),
    y = expression(-Log[10]~"(adjusted p-value)"),
    title = "P85 vs P95",
    tag = "C"
  )
)
```

![](Neg-data---Part-V---Statistical-analysis_files/figure-commonmark/unnamed-chunk-57-1.png)

``` r
ggsave(
  plot = P85_P95_volcano,
  filename = "Plots/P85_P95_volcano_neg.svg",
  width = 190,
  height = 110,
  units = "mm"
)
```
