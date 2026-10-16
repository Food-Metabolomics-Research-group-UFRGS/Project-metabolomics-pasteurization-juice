# Metabolomic Data analysis after drift correction - Negative ionization mode - Study II
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
library(scales)  # para função `seq_gradient_pal()`
library(pheatmap) # for heínsatmaps

# this is at the end hoping that the default select will be that from dplyr
library(tidyverse) # for everything
```

# Read in data

Negative mode.

``` r
# read in metabolomics data
data_long <- read_csv("Notame/corrected_feature_log2_long.csv") 
```

    Rows: 52128 Columns: 10
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (4): sample_name, short_sample_name, treatment, mz_rt
    dbl (6): run_order, soluble_solids, mass, mz, rt, rel_abund

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
data_wide <- read_csv("Notame/corrected_feature_log2_wide.csv")
```

    Rows: 24 Columns: 2178
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr    (3): sample_name, short_sample_name, treatment
    dbl (2175): run_order, soluble_solids, mass, 85.0285_0.3638, 87.0077_0.3573,...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
data_wide <- data_wide  %>%
mutate(sample_or_qc = if_else(str_detect(sample_name, "QC"), true = "QC", false = "Sample")) %>%
  dplyr::select(sample_name, short_sample_name, treatment, run_order, sample_or_qc, everything()) %>%
  select(-mass, -soluble_solids, -sample_name) |>
  rename("sample_name" = short_sample_name)
```

# Assessing data quality

## Log2 transformed

``` r
boxplot_log2_transformed <- data_long %>%
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
  title = "LC-MS (-) feature abundances by Sample",
       subtitle = "Data is log2 transformed")


boxplot_log2_transformed 
```

![](Neg-data---Part-IV---Data-analysis---After-DC_files/figure-commonmark/unnamed-chunk-4-1.png)

We can also look at this data by run order to see if we have any overall
run order effects visible.

``` r
boxplot_log2_transformed_run_order <- data_long %>%
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
  labs(y = "Relative abundance", x = "Sample", 
  title = "LC-MS (-) feature abundances by Sample",
       subtitle = "Data is log2 transformed")

boxplot_log2_transformed_run_order
```

![](Neg-data---Part-IV---Data-analysis---After-DC_files/figure-commonmark/unnamed-chunk-5-1.png)

## PCAs all data

### With QCs

``` r
pca_qc <- prcomp(data_wide[,-c(1:7)], # remove metadata
                 scale = FALSE, # we did our own scaling
                 center = TRUE) # true is the default

summary(pca_qc)
```

    Importance of components:
                              PC1     PC2     PC3     PC4     PC5     PC6     PC7
    Standard deviation     3.0121 2.62007 2.40581 2.22517 2.14164 2.02977 1.97368
    Proportion of Variance 0.1107 0.08379 0.07064 0.06043 0.05598 0.05029 0.04755
    Cumulative Proportion  0.1107 0.19452 0.26517 0.32560 0.38158 0.43187 0.47942
                               PC8     PC9    PC10   PC11    PC12    PC13    PC14
    Standard deviation     1.91872 1.86416 1.83238 1.8148 1.78196 1.74750 1.73897
    Proportion of Variance 0.04493 0.04242 0.04098 0.0402 0.03876 0.03727 0.03691
    Cumulative Proportion  0.52435 0.56677 0.60775 0.6479 0.68670 0.72397 0.76088
                              PC15   PC16    PC17    PC18   PC19    PC20   PC21
    Standard deviation     1.71339 1.6886 1.67331 1.63784 1.6319 1.57890 1.4982
    Proportion of Variance 0.03583 0.0348 0.03418 0.03274 0.0325 0.03043 0.0274
    Cumulative Proportion  0.79671 0.8315 0.86569 0.89843 0.9309 0.96136 0.9888
                              PC22      PC23      PC24
    Standard deviation     0.95979 9.242e-05 2.098e-14
    Proportion of Variance 0.01124 0.000e+00 0.000e+00
    Cumulative Proportion  1.00000 1.000e+00 1.000e+00

Look at how much variance is explained by each PC.

``` r
importance_qc <- summary(pca_qc)$importance %>%
  as.data.frame()

head(importance_qc)
```

                                PC1      PC2     PC3      PC4      PC5      PC6
    Standard deviation     3.012098 2.620073 2.40581 2.225169 2.141636 2.029768
    Proportion of Variance 0.110740 0.083790 0.07064 0.060430 0.055980 0.050290
    Cumulative Proportion  0.110740 0.194520 0.26517 0.325600 0.381580 0.431870
                               PC7      PC8      PC9     PC10     PC11     PC12
    Standard deviation     1.97368 1.918716 1.864163 1.832381 1.814765 1.781961
    Proportion of Variance 0.04755 0.044930 0.042420 0.040980 0.040200 0.038760
    Cumulative Proportion  0.47942 0.524350 0.566770 0.607750 0.647940 0.686700
                               PC13     PC14     PC15     PC16     PC17     PC18
    Standard deviation     1.747501 1.738972 1.713388 1.688557 1.673313 1.637835
    Proportion of Variance 0.037270 0.036910 0.035830 0.034800 0.034180 0.032740
    Cumulative Proportion  0.723970 0.760880 0.796710 0.831520 0.865690 0.898430
                               PC19     PC20    PC21      PC22         PC23
    Standard deviation     1.631852 1.578898 1.49817 0.9597894 9.242203e-05
    Proportion of Variance 0.032500 0.030430 0.02740 0.0112400 0.000000e+00
    Cumulative Proportion  0.930930 0.961360 0.98876 1.0000000 1.000000e+00
                                   PC24
    Standard deviation     2.097934e-14
    Proportion of Variance 0.000000e+00
    Cumulative Proportion  1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_qc)
```

    Warning in geom_bar(stat = "identity", fill = barfill, color = barcolor, :
    Ignoring empty aesthetic: `width`.

![](Neg-data---Part-IV---Data-analysis---After-DC_files/figure-commonmark/unnamed-chunk-8-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_qc)
```

![](Neg-data---Part-IV---Data-analysis---After-DC_files/figure-commonmark/unnamed-chunk-9-1.png)

Make a scores plot but prettier.

``` r
# create a df of pca_qc$x
scores_raw_qc <- as.data.frame(pca_qc$x)

# bind meta-data
scores_qc <- bind_cols(data_wide[,1:4], # first 4 columns
                       scores_raw_qc)
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
       y = glue("PC2: {PC2_percent_qc}%")))
```

![](Neg-data---Part-IV---Data-analysis---After-DC_files/figure-commonmark/unnamed-chunk-11-1.png)

Make a loadings plot (points are features) even though it might not be
that useful.

``` r
fviz_pca_var(pca_qc)
```

    Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ℹ Please use `linewidth` instead.
    ℹ The deprecated feature was likely used in the ggpubr package.
      Please report the issue at <https://github.com/kassambara/ggpubr/issues>.

![](Neg-data---Part-IV---Data-analysis---After-DC_files/figure-commonmark/unnamed-chunk-12-1.png)

See what I mean? Not that useful. There are some functions in PCAtools
that label only the points that most contribute to each PC. Could also
do this manually if its of interest.

### Without QCs

``` r
data_wide_noqc <- data_wide %>%
  filter(sample_or_qc == "Sample")


pca_noqc <- prcomp(data_wide_noqc[,-c(1:7)], # remove metadata
                 scale = FALSE, # we did our own scaling
                 center = TRUE) # true is the default

summary(pca_noqc)
```

    Importance of components:
                              PC1     PC2     PC3     PC4     PC5     PC6     PC7
    Standard deviation     3.7171 3.00204 2.76185 2.64324 2.53534 2.42963 2.38904
    Proportion of Variance 0.1485 0.09689 0.08201 0.07512 0.06911 0.06347 0.06136
    Cumulative Proportion  0.1485 0.24544 0.32745 0.40256 0.47167 0.53514 0.59650
                               PC8     PC9    PC10   PC11    PC12    PC13    PC14
    Standard deviation     2.34397 2.26656 2.20926 2.1587 2.14936 2.09471 2.07485
    Proportion of Variance 0.05907 0.05523 0.05247 0.0501 0.04967 0.04717 0.04628
    Cumulative Proportion  0.65557 0.71080 0.76327 0.8134 0.86304 0.91022 0.95650
                             PC15      PC16
    Standard deviation     2.0115 2.138e-14
    Proportion of Variance 0.0435 0.000e+00
    Cumulative Proportion  1.0000 1.000e+00

Look at how much variance is explained by each PC.

``` r
importance_noqc <- summary(pca_noqc)$importance %>%
  as.data.frame()

head(importance_noqc)
```

                                PC1      PC2      PC3      PC4      PC5      PC6
    Standard deviation     3.717116 3.002041 2.761849 2.643243 2.535335 2.429633
    Proportion of Variance 0.148550 0.096890 0.082010 0.075120 0.069110 0.063470
    Cumulative Proportion  0.148550 0.245440 0.327450 0.402560 0.471670 0.535140
                                PC7      PC8      PC9     PC10     PC11     PC12
    Standard deviation     2.389037 2.343974 2.266562 2.209256 2.158668 2.149358
    Proportion of Variance 0.061360 0.059070 0.055230 0.052470 0.050100 0.049670
    Cumulative Proportion  0.596500 0.655570 0.710800 0.763270 0.813370 0.863040
                               PC13     PC14     PC15         PC16
    Standard deviation     2.094715 2.074848 2.011511 2.137668e-14
    Proportion of Variance 0.047170 0.046280 0.043500 0.000000e+00
    Cumulative Proportion  0.910220 0.956500 1.000000 1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_noqc)
```

    Warning in geom_bar(stat = "identity", fill = barfill, color = barcolor, :
    Ignoring empty aesthetic: `width`.

![](Neg-data---Part-IV---Data-analysis---After-DC_files/figure-commonmark/unnamed-chunk-15-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_noqc)
```

![](Neg-data---Part-IV---Data-analysis---After-DC_files/figure-commonmark/unnamed-chunk-16-1.png)

Make a scores plot but prettier.

``` r
# create a df of pca_qc$x
scores_raw_noqc <- as.data.frame(pca_noqc$x)

# bind meta-data
scores_noqc <- bind_cols(data_wide_noqc[,1:4], # metadata
                         scores_raw_noqc)
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
  labs(x = glue("PC1: {PC1_percent_noqc}%"), 
       y = glue("PC2: {PC2_percent_noqc}%")))
```

![](Neg-data---Part-IV---Data-analysis---After-DC_files/figure-commonmark/unnamed-chunk-18-1.png)

# Saving data for next steps

``` r
feature_table_wide_log2 <- data_wide 

write_csv(feature_table_wide_log2,
          file = "Processed_data/feature_table_wide_log2_neg.csv")
```

``` r
feature_table_long_log2 <- feature_table_wide_log2 %>%
  pivot_longer(cols = 6:ncol(.),
               names_to = "mz_rt",
               values_to = "rel_abund")

write_csv(feature_table_long_log2,
          file = "Processed_data/feature_table_long_log2_neg.csv")
```
