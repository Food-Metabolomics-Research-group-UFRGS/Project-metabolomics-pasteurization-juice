# Metabolomic first pass data analysis - Negative ionization mode - Study II
Giovana Domeneghini Mercali

### Load libraries

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

# this is at the end hoping that the default select will be that from dplyr
library(tidyverse) # for everything
```

## Read in data

Negative mode.

``` r
# read in metabolomics data
metab <- read_csv("Processed_data/neg_post_filtering.csv",
                  trim_ws = TRUE,
                  na = "0") # read in 0s to be NAs.
```

    Rows: 2174 Columns: 28
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr  (1): mz_rt
    dbl (27): row_ID, mz, rt, QC7, QC1, QC2, P_95_1, QC4, P_95_3, P_C_4, P_C_2, ...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# read in meta data
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

Take a quick look at our data.

``` r
# look at first 5 rows, first 5 columns 
metab[1:5,1:5]
```

    # A tibble: 5 × 5
      row_ID mz_rt              mz    rt   QC7
       <dbl> <chr>           <dbl> <dbl> <dbl>
    1    316 85.0285_0.3638   85.0 0.364   790
    2    215 87.0077_0.3573   87.0 0.357  1065
    3    221 96.9591_0.3664   97.0 0.366  4023
    4    193 111.0084_0.3444 111.  0.344 15875
    5   1046 111.0084_1.8865 111.  1.89  12735

``` r
# look at first 5 rows, all columns 
metadata[1:5,]
```

    # A tibble: 5 × 6
      sample_name short_sample_name treatment run_order soluble_solids  mass
      <chr>       <chr>             <chr>         <dbl>          <dbl> <dbl>
    1 P_76_1      P76-1             P76              21            9.6 0.496
    2 P_76_2      P76-2             P76               9            9.6 0.5  
    3 P_76_3      P76-3             P76               8            9.6 0.498
    4 P_76_4      P76-4             P76              14            9.6 0.506
    5 P_85_1      P85-1             P85               3            9.6 0.501

``` r
# check dimensions
dim(metab)
```

    [1] 2174   28

``` r
dim(metadata)
```

    [1] 24  6

``` r
# make metab sample columns all numeric
metab <- metab %>%
  mutate((across(.cols = 5:ncol(.), .fns = as.numeric)))

glimpse(metab)
```

    Rows: 2,174
    Columns: 28
    $ row_ID <dbl> 316, 215, 221, 193, 1046, 412, 502, 185, 74, 627, 432, 759, 206…
    $ mz_rt  <chr> "85.0285_0.3638", "87.0077_0.3573", "96.9591_0.3664", "111.0084…
    $ mz     <dbl> 85.02847, 87.00770, 96.95914, 111.00840, 111.00843, 112.01182, …
    $ rt     <dbl> 0.3637782, 0.3572833, 0.3663571, 0.3443712, 1.8864620, 0.543339…
    $ QC7    <dbl> 790, 1065, 4023, 15875, 12735, 1304, 1304, 1119, 15902, 1416, 9…
    $ QC1    <dbl> 809, 856, 3786, 19371, 12961, 1154, 1164, 1228, 17241, 1400, 11…
    $ QC2    <dbl> 774, 893, 4332, 18671, 16214, 1270, 1400, 1183, 17306, 1486, 11…
    $ P_95_1 <dbl> 841, 891, 4065, 17873, 12992, 1261, 1317, 1261, 16094, 1598, 10…
    $ QC4    <dbl> 851, 993, 3422, 18856, 12752, 1198, 1582, 1601, 19240, 1107, 94…
    $ P_95_3 <dbl> 806, 1104, 3708, 16800, 12781, 1035, 1099, 1035, 16632, 1436, 1…
    $ P_C_4  <dbl> 741, 1071, 3558, 17851, 14162, 1319, 1319, 1348, 17133, 1476, 1…
    $ P_C_2  <dbl> 645, 923, 3807, 17497, 12609, 1173, 1272, 1266, 17240, 1364, 11…
    $ P_85_3 <dbl> 778, 1313, 3522, 18747, 13912, 1109, 1290, 1301, 17056, 1468, 1…
    $ P_76_3 <dbl> 829, 870, 3384, 17912, 13545, 1100, 1391, 1271, 17596, 1556, 11…
    $ P_76_4 <dbl> 730, 1057, 3689, 17238, 12118, 1223, 1288, 1286, 16946, 1549, 1…
    $ P_C_3  <dbl> 685, 969, 3587, 17940, 13704, 916, 1350, 1289, 17228, 1640, 108…
    $ QC8    <dbl> 682, 1182, 3334, 17525, 12188, 1276, 1516, 1160, 18047, 1606, 9…
    $ QC3    <dbl> 666, 1084, 3486, 18477, 12961, 1090, 1291, 1267, 17826, 1500, 1…
    $ P_76_1 <dbl> 738, 905, 3599, 18832, 12958, 1222, 1222, 1206, 17125, 1393, 13…
    $ P_85_2 <dbl> 897, 919, 3300, 16722, 12208, 1179, 1263, 1377, 16835, 1821, 10…
    $ P_76_2 <dbl> 763, 847, 3575, 18562, 11861, 1038, 1234, 1181, 16700, 1475, 13…
    $ P_95_4 <dbl> 849, 1126, 3174, 18064, 12708, 1301, 1280, 1339, 17176, 1327, 1…
    $ QC5    <dbl> 762, 1066, 3807, 18425, 13532, 1269, 1320, 1269, 17290, 1524, 1…
    $ P_85_4 <dbl> 727, 1250, 3495, 16879, 12092, 1233, 1184, 1194, 16556, 1474, 9…
    $ P_85_1 <dbl> 635, 1273, 3306, 18149, 13450, 996, 1260, 1125, 16862, 1974, 13…
    $ P_C_1  <dbl> 916, 980, 3266, 18448, 13654, 1215, 1634, 1272, 17611, 1642, 99…
    $ QC6    <dbl> 815, 909, 3616, 17053, 12709, 1135, 1494, 1140, 17857, 1392, 11…
    $ P_95_2 <dbl> 654, 879, 3600, 17788, 14482, 1205, 1471, 1239, 16350, 1442, 10…

## Wrangle sample names

Here, the samples are in columns and the features are in rows. Samples
are coded so that the first number is the treatment code, and the last
code is the run order. We are going to transpose the data so that
samples are in rows and features are in columns, and we will also import
the metadata about the samples.

``` r
metab_t <- metab %>%
  select(-row_ID, -mz, -rt) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "mz_rt")

# make the first row the column names
colnames(metab_t) <- metab_t[1,]

# then remove the first row and rename the column that should be called sample_name
metab_t <- metab_t[-1,] %>%
  rename(sample_name = mz_rt) %>%
  mutate((across(.cols = 2:ncol(.), .fns = as.numeric)))

metab_t[1:10, 1:10]
```

       sample_name 85.0285_0.3638 87.0077_0.3573 96.9591_0.3664 111.0084_0.3444
    2          QC7            790           1065           4023           15875
    3          QC1            809            856           3786           19371
    4          QC2            774            893           4332           18671
    5       P_95_1            841            891           4065           17873
    6          QC4            851            993           3422           18856
    7       P_95_3            806           1104           3708           16800
    8        P_C_4            741           1071           3558           17851
    9        P_C_2            645            923           3807           17497
    10      P_85_3            778           1313           3522           18747
    11      P_76_3            829            870           3384           17912
       111.0084_1.8865 112.0118_0.5433 112.0119_0.6551 112.0119_0.38 115.0035_0.334
    2            12735            1304            1304          1119          15902
    3            12961            1154            1164          1228          17241
    4            16214            1270            1400          1183          17306
    5            12992            1261            1317          1261          16094
    6            12752            1198            1582          1601          19240
    7            12781            1035            1099          1035          16632
    8            14162            1319            1319          1348          17133
    9            12609            1173            1272          1266          17240
    10           13912            1109            1290          1301          17056
    11           13545            1100            1391          1271          17596

Add in the metadata and make a new column that will indicate whether a
sample is a “sample” or a “QC”. The metadata we have for the samples we
are no longer using (process blanks etc) are removed.

``` r
metab_plus <- left_join(metab_t, metadata, by = "sample_name") %>% # keeps only those samples in metab_t
  mutate(sample_or_qc = if_else(str_detect(sample_name, "QC"), true = "QC", false = "Sample")) %>%
  dplyr::select(sample_name, short_sample_name, treatment, run_order, mass, sample_or_qc, everything()) %>% # and move metadata to the front
mutate((across(.cols = 5:5, .fns = as.numeric)))
```

Go from wide to long data.

``` r
metab_plus_long <- metab_plus %>%
  pivot_longer(cols = -c(sample_name, short_sample_name, treatment, run_order, sample_or_qc, mass),  # remove metadata
               names_to = "mz_rt",
               values_to = "rel_abund")

glimpse(metab_plus_long)
```

    Rows: 52,200
    Columns: 8
    $ sample_name       <chr> "QC7", "QC7", "QC7", "QC7", "QC7", "QC7", "QC7", "QC…
    $ short_sample_name <chr> "QC-7", "QC-7", "QC-7", "QC-7", "QC-7", "QC-7", "QC-…
    $ treatment         <chr> "QC", "QC", "QC", "QC", "QC", "QC", "QC", "QC", "QC"…
    $ run_order         <dbl> 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, …
    $ mass              <dbl> 0.5041, 0.5041, 0.5041, 0.5041, 0.5041, 0.5041, 0.50…
    $ sample_or_qc      <chr> "QC", "QC", "QC", "QC", "QC", "QC", "QC", "QC", "QC"…
    $ mz_rt             <chr> "85.0285_0.3638", "87.0077_0.3573", "96.9591_0.3664"…
    $ rel_abund         <dbl> 790, 1065, 4023, 15875, 12735, 1304, 1304, 1119, 159…

Also add separate columns for mz and rt, and making both numeric.

``` r
metab_plus_long <- metab_plus_long %>%
  separate_wider_delim(cols = mz_rt,
                       delim = "_",
                       names = c("mz", "rt"),
                       cols_remove = FALSE) %>%
  mutate(across(.cols = 7:8, .fns = as.numeric))  %>% # convert mz and rt to numeric
  drop_na(mz, rt) # removing mz and rt values that are NA; maybe not actually necessary
```

    Warning: There were 2 warnings in `mutate()`.
    The first warning was:
    ℹ In argument: `across(.cols = 7:8, .fns = as.numeric)`.
    Caused by warning:
    ! NAs introduced by coercion
    ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.

``` r
# how did that go?
head(metab_plus_long)
```

    # A tibble: 6 × 10
      sample_name short_sample_name treatment run_order  mass sample_or_qc    mz
      <chr>       <chr>             <chr>         <dbl> <dbl> <chr>        <dbl>
    1 QC7         QC-7              QC               23 0.504 QC            85.0
    2 QC7         QC-7              QC               23 0.504 QC            87.0
    3 QC7         QC-7              QC               23 0.504 QC            97.0
    4 QC7         QC-7              QC               23 0.504 QC           111. 
    5 QC7         QC-7              QC               23 0.504 QC           111. 
    6 QC7         QC-7              QC               23 0.504 QC           112. 
    # ℹ 3 more variables: rt <dbl>, mz_rt <chr>, rel_abund <dbl>

## Correction - mass

``` r
metab_corrected_long <- metab_plus_long |>
 mutate(rel_abund = if_else(!is.na(rel_abund), true = rel_abund*0.1/(mass), false = rel_abund)) |>
  select(-mass, -sample_name) |>
  rename("sample_name" = short_sample_name)

head(metab_corrected_long)
```

    # A tibble: 6 × 8
      sample_name treatment run_order sample_or_qc    mz    rt mz_rt       rel_abund
      <chr>       <chr>         <dbl> <chr>        <dbl> <dbl> <chr>           <dbl>
    1 QC-7        QC               23 QC            85.0 0.364 85.0285_0.…      157.
    2 QC-7        QC               23 QC            87.0 0.357 87.0077_0.…      211.
    3 QC-7        QC               23 QC            97.0 0.366 96.9591_0.…      798.
    4 QC-7        QC               23 QC           111.  0.344 111.0084_0…     3149.
    5 QC-7        QC               23 QC           111.  1.89  111.0084_1…     2526.
    6 QC-7        QC               23 QC           112.  0.543 112.0118_0…      259.

## Data summaries

What mass range do I have?

``` r
range(metab_corrected_long$mz)
```

    [1]   85.0285 1287.3252

What retention time range do I have?

``` r
range(metab_corrected_long$rt)
```

    [1] 0.3003 7.9723

How many samples are in each of my meta-data groups?

``` r
# make wide data to make some calculations easier
metab_corrected_wide <- metab_corrected_long %>%
  dplyr::select(-mz, -rt) %>%
  pivot_wider(names_from = mz_rt,
              values_from = rel_abund)

# by sample vs QC
metab_corrected_wide %>%
  count(sample_or_qc)
```

    # A tibble: 2 × 2
      sample_or_qc     n
      <chr>        <int>
    1 QC               8
    2 Sample          16

What does my data coverage across mz and rt look like?

``` r
metab_corrected_long %>%
  group_by(mz_rt) %>% # so we only have one point per feature
  ggplot(aes(x = rt, y = mz)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Retention time (min)",
       y = "Mass to charge ratio (m/z)",
       title = "m/z by retention time plot (all features)",
       subtitle = "C18 reversed phase, negative ionization mode")
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-13-1.png)

Distribution of masses

``` r
metab_corrected_long %>%
  group_by(mz_rt) %>%
  ggplot(aes(x = mz)) +
  geom_histogram(bins = 100) +
  theme_minimal() +
  labs(x = "Mass to charge ratio (m/z)",
       y = "Number of features",
       title = "Distribution of features by mass")
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-14-1.png)

Distribution of retention times

``` r
metab_corrected_long %>%
  group_by(mz_rt) %>%
  ggplot(aes(x = rt)) +
  geom_density() +
  theme_minimal() +
  labs(x = "Retention time",
       y = "Number of features",
       title = "Distribution of features by retention time")
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-15-1.png)

## Missing data

### Surveying missingness

How many missing values are there for each feature? In this dataset,
missing values are coded as zero.

``` r
# all data including QCs
# how many missing values are there for each feature (row)
na_by_feature <- rowSums(is.na(metab)) %>%
  as.data.frame() %>%
  rename(missing_values = 1)


na_by_feature %>%
  ggplot(aes(x = missing_values)) +
  geom_histogram(bins = 67) + # since 67 samples
  theme_minimal() + 
  labs(title = "Number of missing values for each feature",
       x = "Number of missing values",
       y = "How many features have that many missing values")
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-16-1.png)

How many features have no missing values?

``` r
na_by_feature %>%
  count(missing_values == 0)
```

      missing_values == 0    n
    1                TRUE 2174

How many missing values are there for each sample?

``` r
# all data including QCs
# how many missing values are there for each feature (row)
na_by_sample <- colSums(is.na(metab)) %>%
  as.data.frame() %>%
  rename(missing_values = 1) %>%
  rownames_to_column(var = "feature") %>%
  filter(!feature %in% c("mz_rt", "rt", "mz", "row_ID"))

na_by_sample %>%
  ggplot(aes(x = missing_values)) +
  geom_histogram(bins = 100) + # there are 3473 features
  theme_minimal() + 
  labs(title = "Number of missing values for each sample",
       x = "Number of missing values",
       y = "How many samples have that many missing values")
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-18-1.png)

Which features have a lot of missing values?

``` r
contains_NAs_feature <- metab_corrected_long %>%
  group_by(mz_rt) %>%
  count(is.na(rel_abund)) %>%
  filter(`is.na(rel_abund)` == TRUE) %>% 
  arrange(desc(n))

head(contains_NAs_feature)
```

    # A tibble: 0 × 3
    # Groups:   mz_rt [0]
    # ℹ 3 variables: mz_rt <chr>, is.na(rel_abund) <lgl>, n <int>

Which samples have a lot of missing values?

``` r
contains_NAs_sample <- metab_corrected_long %>%
  group_by(sample_name) %>%
  count(is.na(rel_abund)) %>%
  filter(`is.na(rel_abund)` == TRUE) %>%
  arrange(desc(n))

head(contains_NAs_sample)
```

    # A tibble: 0 × 3
    # Groups:   sample_name [0]
    # ℹ 3 variables: sample_name <chr>, is.na(rel_abund) <lgl>, n <int>

Are there any missing values in the QCs? (There shouldn’t be.)

``` r
metab_QC <- metab %>%
  dplyr::select(contains("QC"))

na_by_sample <- colSums(is.na(metab_QC)) %>%
  as.data.frame() %>%
  rename(missing_values = 1) %>%
  rownames_to_column(var = "feature") %>%
  filter(!feature == "mz_rt")

sum(na_by_sample$missing_values) # nope
```

    [1] 0

### Imputing missing values

This is an optional step but some downstream analyses don’t handle
missingness well. Here we are imputing missing data with half the lowest
value observed for that feature.

``` r
# grab only the feature data and leave metadata
metab_corrected_wide_imputed <- metab_corrected_wide %>%
  dplyr::select(-c(1:4)) # the metadata columns

metab_corrected_wide_imputed[] <- lapply(metab_corrected_wide_imputed,
                                 function(x) ifelse(is.na(x), min(x, na.rm = TRUE)/2, x))

metab_corrected_wide_imputed_notame <- metab_corrected_wide_imputed

# bind back the metadata
metab_corrected_wide_imputed <- bind_cols(metab_corrected_wide[,1:4], metab_corrected_wide_imputed)
```

For notame later:

``` r
data_notame <- bind_cols(metab_corrected_wide[,1], metab_corrected_wide_imputed_notame) %>%
 as.data.frame()
```

Did imputing work?

``` r
# count missing values
metab_corrected_wide_imputed %>%
  dplyr::select(-c(1:4)) %>% # where the metadata is
  is.na() %>%
  sum()
```

    [1] 0

Create long imputed dataset.

``` r
metab_corrected_long_imputed <- metab_corrected_wide_imputed %>%
  pivot_longer(cols = 5:ncol(.),
               names_to = "mz_rt",
               values_to = "rel_abund")

head(metab_corrected_long_imputed)
```

    # A tibble: 6 × 6
      sample_name treatment run_order sample_or_qc mz_rt           rel_abund
      <chr>       <chr>         <dbl> <chr>        <chr>               <dbl>
    1 QC-7        QC               23 QC           85.0285_0.3638       157.
    2 QC-7        QC               23 QC           87.0077_0.3573       211.
    3 QC-7        QC               23 QC           96.9591_0.3664       798.
    4 QC-7        QC               23 QC           111.0084_0.3444     3149.
    5 QC-7        QC               23 QC           111.0084_1.8865     2526.
    6 QC-7        QC               23 QC           112.0118_0.5433      259.

Let’s also make separate mz and rt columns.

``` r
metab_corrected_long_imputed <- metab_corrected_long_imputed %>%
  separate_wider_delim(cols = mz_rt,
                       delim = "_",
                       names = c("mz", "rt"),
                       cols_remove = FALSE)

metab_corrected_long_imputed$mz <- as.numeric(metab_corrected_long_imputed$mz)
metab_corrected_long_imputed$rt <- as.numeric(metab_corrected_long_imputed$rt)
```

## Feature clustering with `notame`

We want to cluster features that likely come from the same metabolite
together, and we can do this using the package `notame`. You can learn
more
[here](http://127.0.0.1:24885/library/notame/doc/feature_clustering.html).

``` r
browseVignettes("notame")
```

Let’s make a m/z by retention time plot again before we start.

``` r
(before_notame <- metab_corrected_long_imputed %>%
  group_by(mz_rt) %>% # so we only have one point per feature
  ggplot(aes(x = rt, y = mz)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Retention time (min)",
       y = "Mass to charge ratio (m/z)",
       title = "m/z by retention time plot before notame, 2997 features",
       subtitle = "C18 reverse phase, negative ionization mode"))
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-28-1.png)

### Wrangling data

Transpose the wide data for notame and wrangle to the right format.
Below is info from the documentation:

- Data should be a data frame containing the abundances of features in
  each sample, one row per sample, each feature in a separate column
- Features should be a data frame containing information about the
  features, namely feature name (should be the same as the column name
  in data), mass and retention time

``` r
# change to results to numeric
data_notame <- data_notame %>%
  mutate(across(-sample_name, as.numeric))

tibble(data_notame)
```

    # A tibble: 24 × 2,175
       sample_name `85.0285_0.3638` `87.0077_0.3573` `96.9591_0.3664`
       <chr>                  <dbl>            <dbl>            <dbl>
     1 QC-7                    157.             211.             798.
     2 QC-1                    160.             170.             751.
     3 QC-2                    154.             177.             859.
     4 P95-1                   170.             180.             823.
     5 QC-4                    169.             197.             679.
     6 P95-3                   161.             220.             739.
     7 UJ-4                    143.             206.             686.
     8 UJ-2                    124.             177.             731.
     9 P85-3                   157.             265.             710.
    10 P76-3                   167.             175.             680.
    # ℹ 14 more rows
    # ℹ 2,171 more variables: `111.0084_0.3444` <dbl>, `111.0084_1.8865` <dbl>,
    #   `112.0118_0.5433` <dbl>, `112.0119_0.6551` <dbl>, `112.0119_0.38` <dbl>,
    #   `115.0035_0.334` <dbl>, `115.0036_0.7684` <dbl>, `115.0036_0.5518` <dbl>,
    #   `115.0038_0.9791` <dbl>, `116.0069_0.3527` <dbl>, `117.0193_0.8034` <dbl>,
    #   `128.0354_0.7648` <dbl>, `129.0192_1.0996` <dbl>, `129.0193_0.8865` <dbl>,
    #   `129.0193_0.6671` <dbl>, `129.0193_0.3243` <dbl>, …

Create df with features.

``` r
features <- metab_corrected_long_imputed %>%
  dplyr::select(mz_rt, mz, rt) %>%
  mutate(across(c(mz, rt), as.numeric)) %>%
  as.data.frame() %>%
  distinct()

glimpse(features)
```

    Rows: 2,174
    Columns: 3
    $ mz_rt <chr> "85.0285_0.3638", "87.0077_0.3573", "96.9591_0.3664", "111.0084_…
    $ mz    <dbl> 85.0285, 87.0077, 96.9591, 111.0084, 111.0084, 112.0118, 112.011…
    $ rt    <dbl> 0.3638, 0.3573, 0.3664, 0.3444, 1.8865, 0.5433, 0.6551, 0.3800, …

``` r
class(features)
```

    [1] "data.frame"

### Find connections

Set `cache = TRUE` for this chunk since its a bit slow especially if you
have a lot of features (this step took ~10 min).

``` r
connection <- find_connections(data = data_notame,
                               features = features,
                               corr_thresh = 0.9,
                               rt_window = 1/60,
                               name_col = "mz_rt",
                               mz_col = "mz",
                               rt_col = "rt")
```

    Warning: executing %dopar% sequentially: no parallel backend registered

    [1] 100
    [1] 200
    [1] 300
    [1] 400
    [1] 500
    [1] 600
    [1] 700
    [1] 800
    [1] 900
    [1] 1000
    [1] 1100
    [1] 1200
    [1] 1300
    [1] 1400
    [1] 1500
    [1] 1600
    [1] 1700
    [1] 1800
    [1] 1900
    [1] 2000
    [1] 2100

``` r
head(connection)
```

                    x               y       cor rt_diff mz_diff
    1 350.2464_4.7068  360.275_4.7078 0.9053289  0.0010 10.0286
    2  540.3297_5.884 541.3328_5.8869 0.9545825  0.0029  1.0031

### Clustering

Now that we have found all of the features that are connected based on
the parameters we have set, we now need to find clusters.

``` r
clusters <- find_clusters(connections = connection, 
                          d_thresh = 0.8)
```

    2 components found

Assign a cluster ID to each feature to keep, and the feature that is
picked is the one with the highest median peak intensity across the
samples.

``` r
# assign a cluster ID to all features
# clusters are named after feature with highest median peak height
features_clustered <- assign_cluster_id(data_notame, 
                                        clusters, 
                                        features, 
                                        name_col = "mz_rt")
```

Export out a list of your clusters this way you can use this later
during metabolite ID.

``` r
# export clustered feature list this way you have it
write_csv(features_clustered,
          "Processed_data/Clusters_neg.csv")
```

Pull data out from the clusters and see how many features we
removed/have now.

``` r
# lets see how many features are removed when we only keep one feature per cluster
pulled <- pull_clusters(data_notame, features_clustered, name_col = "mz_rt")

cluster_data <- pulled$cdata
cluster_features <- pulled$cfeatures

# how many features did we originally have after filtering?
nrow(metab)
```

    [1] 2174

``` r
# how many features got removed during clustering?
nrow(metab) - nrow(cluster_features)
```

    [1] 2

``` r
# what percentage of the original features were removed?
((nrow(metab) - nrow(cluster_features))/nrow(metab)) * 100
```

    [1] 0.09199632

Reduce our dataset to include only our new clusters. `cluster_data`
contains only the retained clusters, while `cluster_features` tells you
also which features are a part of each cluster.

``` r
# combined metadata_plus with cluster_features

metab_imputed_clustered_wide <- left_join(metab_corrected_wide_imputed[,1:4], cluster_data,
                                          by = "sample_name") 

dim(metab_imputed_clustered_wide) # we have 2474 features since 4 metadata columns
```

    [1]   24 2176

``` r
# make a long/tidy df
metab_imputed_clustered_long <- metab_imputed_clustered_wide %>%
  pivot_longer(cols = 5:ncol(.),
               names_to = "mz_rt",
               values_to = "rel_abund") %>%
  separate_wider_delim(cols = mz_rt, # make separate columns for mz and rt too
                       delim = "_",
                       names = c("mz", "rt"),
                       cols_remove = FALSE) %>%
  mutate(across(.cols = c("mz", "rt"), .fns = as.numeric)) # make mz and rt numeric
```

Write out the final clustered dataset.

``` r
write_csv(metab_imputed_clustered_long,
  "Processed_data/features_clustered_long_neg.csv")

write_csv(metab_imputed_clustered_wide,
  "Processed_data/features_clustered_wide_neg.csv")
```

``` r
metab_imputed_clustered_wide_t <- metab_imputed_clustered_wide %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "mz_rt")

write_csv(metab_imputed_clustered_wide_t,
  "Processed_data/features_manuscript_neg.csv")
```

Let’s look at a m/z by retention time plot again after clustering.

``` r
(after_notame <- metab_imputed_clustered_long %>%
  group_by(mz_rt) %>% # so we only have one point per feature
  ggplot(aes(x = rt, y = mz)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Retention time (min)",
       y = "Mass to charge ratio (m/z)",
       title = "m/z by retention time plot after notame, 2289 features",
       subtitle = "C18 reverse phase, negative ionization mode"))
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-40-1.png)

``` r
before_notame / after_notame
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-41-1.png)

## Assessing data quality

Let’s make sure that our data is of good quality.

### Untransformed data

First we are going to convert the type of some of the columns to match
what we want (e.g., run order converted to numeric, treatment to
factor).

``` r
tibble(metab_imputed_clustered_long)
```

    # A tibble: 52,128 × 8
       sample_name treatment run_order sample_or_qc    mz    rt mz_rt      rel_abund
       <chr>       <chr>         <dbl> <chr>        <dbl> <dbl> <chr>          <dbl>
     1 QC-7        QC               23 QC            85.0 0.364 85.0285_0…      157.
     2 QC-7        QC               23 QC            87.0 0.357 87.0077_0…      211.
     3 QC-7        QC               23 QC            97.0 0.366 96.9591_0…      798.
     4 QC-7        QC               23 QC           111.  0.344 111.0084_…     3149.
     5 QC-7        QC               23 QC           111.  1.89  111.0084_…     2526.
     6 QC-7        QC               23 QC           112.  0.543 112.0118_…      259.
     7 QC-7        QC               23 QC           112.  0.655 112.0119_…      259.
     8 QC-7        QC               23 QC           112.  0.38  112.0119_…      222.
     9 QC-7        QC               23 QC           115.  0.334 115.0035_…     3155.
    10 QC-7        QC               23 QC           115.  0.768 115.0036_…      281.
    # ℹ 52,118 more rows

``` r
# make run_order numeric
metab_imputed_clustered_long$run_order <- as.numeric(metab_imputed_clustered_long$run_order)

# make treatment and sample_or_qc a factor (i.e., categorical)
metab_imputed_clustered_long$treatment <- as.factor(metab_imputed_clustered_long$treatment)
metab_imputed_clustered_long$sample_or_qc <- as.factor(metab_imputed_clustered_long$sample_or_qc)

# did it work?
tibble(metab_imputed_clustered_long)
```

    # A tibble: 52,128 × 8
       sample_name treatment run_order sample_or_qc    mz    rt mz_rt      rel_abund
       <chr>       <fct>         <dbl> <fct>        <dbl> <dbl> <chr>          <dbl>
     1 QC-7        QC               23 QC            85.0 0.364 85.0285_0…      157.
     2 QC-7        QC               23 QC            87.0 0.357 87.0077_0…      211.
     3 QC-7        QC               23 QC            97.0 0.366 96.9591_0…      798.
     4 QC-7        QC               23 QC           111.  0.344 111.0084_…     3149.
     5 QC-7        QC               23 QC           111.  1.89  111.0084_…     2526.
     6 QC-7        QC               23 QC           112.  0.543 112.0118_…      259.
     7 QC-7        QC               23 QC           112.  0.655 112.0119_…      259.
     8 QC-7        QC               23 QC           112.  0.38  112.0119_…      222.
     9 QC-7        QC               23 QC           115.  0.334 115.0035_…     3155.
    10 QC-7        QC               23 QC           115.  0.768 115.0036_…      281.
    # ℹ 52,118 more rows

Let’s make a boxplot to see how the metabolite abundance looks across
different samples.

``` r
metab_imputed_clustered_long %>%
  ggplot(aes(x = sample_name, y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  labs(title = "LC-MS (-) feature abundances by sample",
       subtitle = "Data is unscaled",
       x = "Sample",
       y = "Relative abundance")
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-43-1.png)

Can’t really see anything because data is skewed.

### Transformed data

#### Log2 transformed

We will log2 transform our data.

``` r
metab_imputed_clustered_long_log2 <- metab_imputed_clustered_long %>%
  mutate(rel_abund = log2(rel_abund))
```

And then plot.

``` r
boxplot_log2_transformed <- metab_imputed_clustered_long_log2 %>%
  ggplot(aes(x = sample_name, y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Relative abundance", x = "Sample", 
  title = "LC-MS (-) feature abundances by Sample",
       subtitle = "Data is log2 transformed")

boxplot_log2_transformed 
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-45-1.png)

We can also look at this data by run order to see if we have any overall
run order effects visible.

``` r
boxplot_log2_transformed_run_order <- metab_imputed_clustered_long_log2 %>%
  mutate(sample_name = fct_reorder(sample_name, run_order)) %>%
  ggplot(aes(x = sample_name , y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "Relative abundance", x = "Sample")

boxplot_log2_transformed_run_order
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-46-1.png)

#### Log10 transformed

We will log10 transform our data.

``` r
metab_imputed_clustered_long_log10 <- metab_imputed_clustered_long %>%
  mutate(rel_abund = log10(rel_abund))
```

We can look at this data where we group by species.

``` r
metab_imputed_clustered_long_log10 %>%
  ggplot(aes(x = sample_name , y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "LC-MS (-) feature abundances by sample",
       subtitle = "Data is log10 transformed",
       y = "Relative abundance", x = "Sample")
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-48-1.png)

We can also look at this data by run order to see if we have any overall
run order effects visible.

``` r
metab_imputed_clustered_long_log10 %>%
  mutate(sample_name = fct_reorder(sample_name, run_order)) %>%
  ggplot(aes(x = sample_name , y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "LC-MS (-) feature abundances by sample",
       subtitle = "Data is log10 transformed",
       y = "Relative abundance", x = "Sample")
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-49-1.png)

#### Autoscaled

Scales to unit variance, where the mean for each feature is 0 with a
standard deviation of 1. This works well when all metabolites are
considered to be equivalently important though measurement errors can be
inflated. We’ve never actually scaled data this way for a project but
whatever here it is.

``` r
metab_imputed_clustered_wide <- metab_imputed_clustered_long %>%
  select(-mz, -rt) %>%
  pivot_wider(names_from = "mz_rt",
              values_from = "rel_abund")

# autoscaling on the now zero-centered matrix
autoscaled <- 
    bind_cols(metab_imputed_clustered_wide[1:4], # metadata
            lapply(metab_imputed_clustered_wide[5:ncol(metab_imputed_clustered_wide)], # metab data
                   scale)) # scale to mean 0 sd 1

autoscaled[1:10,1:10]
```

    # A tibble: 10 × 10
       sample_name treatment run_order sample_or_qc `85.0285_0.3638`[,1]
       <chr>       <fct>         <dbl> <fct>                       <dbl>
     1 QC-7        QC               23 QC                         0.292 
     2 QC-1        QC                1 QC                         0.527 
     3 QC-2        QC                2 QC                         0.0948
     4 P95-1       P95              13 Sample                     1.13  
     5 QC-4        QC               12 QC                         1.05  
     6 P95-3       P95               5 Sample                     0.540 
     7 UJ-4        UJ               16 Sample                    -0.576 
     8 UJ-2        UJ                6 Sample                    -1.76  
     9 P85-3       P85              10 Sample                     0.293 
    10 P76-3       P76               8 Sample                     0.908 
    # ℹ 5 more variables: `87.0077_0.3573` <dbl[,1]>, `96.9591_0.3664` <dbl[,1]>,
    #   `111.0084_0.3444` <dbl[,1]>, `111.0084_1.8865` <dbl[,1]>,
    #   `112.0118_0.5433` <dbl[,1]>

``` r
autoscaled_long <- autoscaled %>%
  pivot_longer(cols = 5:ncol(.),
               names_to = "mz_rt",
               values_to = "rel_abund") %>%
  mutate(treatment = as.factor(treatment))

autoscaled_long %>%
  mutate(sample_name = fct_reorder(sample_name, run_order)) %>%
  ggplot(aes(x = sample_name , y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "LC-MS (-) Feature Abundances by Sample",
       subtitle = "Data is autoscaled",
       y = "Relative abundance")
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-51-1.png)

This is weird we won’t use this.

#### Pareto scaled

Pareto scaling scales but keeps the fidelity of the original differences
in absolute measurement value more than autoscaling. Often data is
Pareto scaled after log transformation

``` r
metab_imputed_clustered_wide_log2 <- metab_imputed_clustered_long_log2 %>%
  select(-mz, -rt) %>%
  pivot_wider(names_from = "mz_rt",
              values_from = "rel_abund")

metab_imputed_clustered_wide_log2_metabs <- 
  metab_imputed_clustered_wide_log2[,5:ncol(metab_imputed_clustered_wide_log2)]

pareto_scaled <- 
  IMIFA::pareto_scale(metab_imputed_clustered_wide_log2_metabs, center = FALSE)

pareto_scaled <- bind_cols(metab_imputed_clustered_wide_log2[,1:4], pareto_scaled)
```

``` r
pareto_scaled_long <- pareto_scaled %>%
  pivot_longer(cols = 5:ncol(.),
               names_to = "mz_rt",
               values_to = "rel_abund")

pareto_scaled_long %>%
  # mutate(short_sample_name = fct_reorder(short_sample_name, treatment)) %>%
  ggplot(aes(x = sample_name, y = rel_abund, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "LC-MS (-) feature abundances by sample",
       subtitle = "Data is Pareto scaled",
       y = "Relative abundance", x = "Sample")
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-52-1.png)

I think pareto scaling is making everything look super the same. I am
going to use log2 transformed data for the rest of this analysis.

## PCAs all data

### With QCs

``` r
pca_qc <- prcomp(metab_imputed_clustered_wide_log2[,-c(1:4)], # remove metadata
                 scale = FALSE, # we did our own scaling
                 center = TRUE) # true is the default

summary(pca_qc)
```

    Importance of components:
                              PC1     PC2     PC3     PC4     PC5     PC6     PC7
    Standard deviation     2.9002 2.46970 2.25522 2.14817 2.03509 1.98713 1.91532
    Proportion of Variance 0.1007 0.07301 0.06088 0.05524 0.04957 0.04727 0.04391
    Cumulative Proportion  0.1007 0.17369 0.23457 0.28981 0.33938 0.38665 0.43056
                              PC8     PC9    PC10    PC11   PC12    PC13    PC14
    Standard deviation     1.9041 1.87010 1.82409 1.79977 1.7770 1.75552 1.74788
    Proportion of Variance 0.0434 0.04186 0.03983 0.03877 0.0378 0.03689 0.03657
    Cumulative Proportion  0.4740 0.51582 0.55565 0.59442 0.6322 0.66911 0.70567
                              PC15    PC16    PC17    PC18    PC19    PC20    PC21
    Standard deviation     1.72677 1.72335 1.69823 1.68870 1.63900 1.62667 1.61216
    Proportion of Variance 0.03569 0.03555 0.03452 0.03413 0.03216 0.03167 0.03111
    Cumulative Proportion  0.74137 0.77692 0.81144 0.84557 0.87773 0.90940 0.94051
                              PC22   PC23      PC24
    Standard deviation     1.59331 1.5593 2.117e-14
    Proportion of Variance 0.03039 0.0291 0.000e+00
    Cumulative Proportion  0.97090 1.0000 1.000e+00

Look at how much variance is explained by each PC.

``` r
importance_qc <- summary(pca_qc)$importance %>%
  as.data.frame()

head(importance_qc)
```

                                PC1      PC2      PC3      PC4     PC5     PC6
    Standard deviation     2.900187 2.469703 2.255218 2.148167 2.03509 1.98713
    Proportion of Variance 0.100680 0.073010 0.060880 0.055240 0.04957 0.04727
    Cumulative Proportion  0.100680 0.173690 0.234570 0.289810 0.33938 0.38665
                                PC7      PC8      PC9     PC10     PC11    PC12
    Standard deviation     1.915315 1.904127 1.870101 1.824092 1.799766 1.77697
    Proportion of Variance 0.043910 0.043400 0.041860 0.039830 0.038770 0.03780
    Cumulative Proportion  0.430560 0.473960 0.515820 0.555650 0.594420 0.63222
                               PC13    PC14     PC15     PC16     PC17     PC18
    Standard deviation     1.755522 1.74788 1.726767 1.723346 1.698233 1.688698
    Proportion of Variance 0.036890 0.03657 0.035690 0.035550 0.034520 0.034130
    Cumulative Proportion  0.669110 0.70567 0.741370 0.776920 0.811440 0.845570
                               PC19     PC20     PC21     PC22     PC23
    Standard deviation     1.638998 1.626668 1.612162 1.593305 1.559278
    Proportion of Variance 0.032160 0.031670 0.031110 0.030390 0.029100
    Cumulative Proportion  0.877730 0.909400 0.940510 0.970900 1.000000
                                   PC24
    Standard deviation     2.116953e-14
    Proportion of Variance 0.000000e+00
    Cumulative Proportion  1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_qc)
```

    Warning in geom_bar(stat = "identity", fill = barfill, color = barcolor, :
    Ignoring empty aesthetic: `width`.

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-55-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_qc)
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-56-1.png)

Make a scores plot but prettier.

``` r
# create a df of pca_qc$x
scores_raw_qc <- as.data.frame(pca_qc$x)

# bind meta-data
scores_qc <- bind_cols(metab_imputed_clustered_wide_log2[,1:4], # first 4 columns
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
  geom_point(shape = 21, color = "black") +
  theme_minimal() +
  labs(x = glue("PC1: {PC1_percent_qc}%"), 
       y = glue("PC2: {PC2_percent_qc}%"), 
       title = "PCA scores plot colored by treatment, negative mode"))
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-58-1.png)

Make a loadings plot (points are features) even though it might not be
that useful.

``` r
fviz_pca_var(pca_qc)
```

    Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ℹ Please use `linewidth` instead.
    ℹ The deprecated feature was likely used in the ggpubr package.
      Please report the issue at <https://github.com/kassambara/ggpubr/issues>.

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-59-1.png)

See what I mean? Not that useful. There are some functions in PCAtools
that label only the points that most contribute to each PC. Could also
do this manually if its of interest.

### Without QCs

``` r
metab_imputed_clustered_wide_log2_noqc <- metab_imputed_clustered_wide_log2 %>%
  filter(sample_or_qc == "Sample")


pca_noqc <- prcomp(metab_imputed_clustered_wide_log2_noqc[,-c(1:4)], # remove metadata
                 scale = FALSE, # we did our own scaling
                 center = TRUE) # true is the default

summary(pca_noqc)
```

    Importance of components:
                              PC1     PC2     PC3     PC4     PC5    PC6     PC7
    Standard deviation     3.5752 2.87362 2.63977 2.47725 2.46620 2.3541 2.33993
    Proportion of Variance 0.1456 0.09403 0.07935 0.06988 0.06926 0.0631 0.06235
    Cumulative Proportion  0.1456 0.23958 0.31893 0.38880 0.45806 0.5212 0.58351
                               PC8     PC9    PC10    PC11   PC12    PC13    PC14
    Standard deviation     2.27941 2.22610 2.17288 2.14696 2.1225 2.07624 2.07051
    Proportion of Variance 0.05916 0.05643 0.05376 0.05249 0.0513 0.04909 0.04882
    Cumulative Proportion  0.64267 0.69910 0.75286 0.80535 0.8566 0.90573 0.95455
                              PC15      PC16
    Standard deviation     1.99794 2.173e-14
    Proportion of Variance 0.04545 0.000e+00
    Cumulative Proportion  1.00000 1.000e+00

Look at how much variance is explained by each PC.

``` r
importance_noqc <- summary(pca_noqc)$importance %>%
  as.data.frame()

head(importance_noqc)
```

                                PC1      PC2      PC3      PC4      PC5      PC6
    Standard deviation     3.575216 2.873618 2.639773 2.477252 2.466199 2.354089
    Proportion of Variance 0.145550 0.094030 0.079350 0.069880 0.069260 0.063100
    Cumulative Proportion  0.145550 0.239580 0.318930 0.388800 0.458060 0.521160
                                PC7     PC8      PC9     PC10     PC11     PC12
    Standard deviation     2.339926 2.27941 2.226099 2.172877 2.146963 2.122455
    Proportion of Variance 0.062350 0.05916 0.056430 0.053760 0.052490 0.051300
    Cumulative Proportion  0.583510 0.64267 0.699100 0.752860 0.805350 0.856640
                               PC13     PC14     PC15         PC16
    Standard deviation     2.076241 2.070512 1.997938 2.173377e-14
    Proportion of Variance 0.049090 0.048820 0.045450 0.000000e+00
    Cumulative Proportion  0.905730 0.954550 1.000000 1.000000e+00

Generate a scree plot.

``` r
fviz_eig(pca_noqc)
```

    Warning in geom_bar(stat = "identity", fill = barfill, color = barcolor, :
    Ignoring empty aesthetic: `width`.

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-62-1.png)

Generate a scores plot (points are samples) quickly with `fviz_pca_ind`.

``` r
fviz_pca_ind(pca_noqc)
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-63-1.png)

Make a scores plot but prettier.

``` r
# create a df of pca_qc$x
scores_raw_noqc <- as.data.frame(pca_noqc$x)

# bind meta-data
scores_noqc <- bind_cols(metab_imputed_clustered_wide_log2_noqc[,1:4], # metadata
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
  geom_point(shape = 21, color = "black") +
  theme_minimal() +
  labs(x = glue("PC1: {PC1_percent_noqc}%"), 
       y = glue("PC2: {PC2_percent_noqc}%"), 
       title = "PCA scores plot colored by treatment, ositive mode"))
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-65-1.png)

Make a loadings plot (points are features) even though it might not be
that useful.

``` r
fviz_pca_var(pca_noqc)
```

![](Neg-Data---Part-II---Data-analysis_files/figure-commonmark/unnamed-chunk-66-1.png)

See what I mean? Not that useful. There are some functions in PCAtools
that label only the points that most contribute to each PC. Could also
do this manually if its of interest.
