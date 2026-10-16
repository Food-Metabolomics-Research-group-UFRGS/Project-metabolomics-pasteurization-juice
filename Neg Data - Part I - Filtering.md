# Filtering - Untargeted Metabolomics Data - Negative Mode - Study II
Giovana Domeneghini Mercali
2023-08-30

# Loading libraries

``` r
library(tidyverse)
library(janitor) # if you want to clean_names()
```

Once you get deconvoluted data from MZmine or similar programs, you need
to wrangle your data in such a way that you can conduct your analysis on
it.

# Read in data

First we want to read in our raw data. The code here is to read in data
directly from MZmine, so if you are using a different program for
deconvolution, you might need to make some light adjustments.

``` r
metadata_juice <- read_csv(file = "Data/Feature_list_all_data_neg.csv",
                      col_names = TRUE, # has headers
                      na = "0")
```

    Rows: 3887 Columns: 29
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    dbl (29): row ID, row m/z, row retention time, ProcessBlank1.mzML Peak heigh...

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# How many features do we have?
dim(metadata_juice)
```

    [1] 3887   29

``` r
# Look at beginning of the dataset

metadata_juice[1:8, 1:8]
```

    # A tibble: 8 × 8
      `row ID` `row m/z` `row retention time` `ProcessBlank1.mzML Peak height`
         <dbl>     <dbl>                <dbl>                            <dbl>
    1      316      85.0                0.364                               28
    2      215      87.0                0.357                               47
    3      142      89.0                0.286                              356
    4      685      89.0                0.841                              114
    5      221      97.0                0.366                              129
    6     4677     102.                 8.87                               907
    7      193     111.                 0.344                               71
    8     1046     111.                 1.89                                NA
    # ℹ 4 more variables: `ProcessBlank2.mzML Peak height` <dbl>,
    #   `QC7.mzML Peak height` <dbl>, `QC1.mzML Peak height` <dbl>,
    #   `QC2.mzML Peak height` <dbl>

``` r
# str(metadata_juice)
# names(metadata_juice)
#length(metadata_juice)
```

Note there is no metadata included in this file. Just m/z, retention
time, and a column for each sample, where values are peak heights. We
are using peak height instead of peak area because it is less dependent
on bad peak shape which you get sometimes with Metabolomics.

# RT filter

You might have deconvoluted more data than you plan to use in your
analysis. For example, you may want to exclude the first bit and last
bit of your run, since you do not expect to have good reproducibility in
those areas.

Here, we are filtering to only include features that elute between
0.3-8.0 min of this 10 min run.

``` r
metadata_juice_RTfilt <- metadata_juice %>%
  filter(between(`row retention time`, 0.3, 8.0))

# Did it work?
range(metadata_juice_RTfilt$`row retention time`)
```

    [1] 0.3003378 7.9830290

``` r
# How many features do we have now?
dim(metadata_juice_RTfilt)
```

    [1] 3648   29

# Cleaning up data

## Create mz_rt

This creates a unique identifier for each feature using its
mass-to-charge ratio (m/z) and retention time (RT).

``` r
MZ_RT <- metadata_juice_RTfilt %>%
  mutate(mz = round(metadata_juice_RTfilt$`row m/z`, digits = 4), # decrease number of decimals for m/z & rt
         rt = round(metadata_juice_RTfilt$`row retention time`, digits = 4),
         .before = 1,
         .keep = "unused") %>%
  unite(mz_rt, c(mz, rt), remove = TRUE) %>% # Combine m/z & rt with _ in between
  rename(mz = `row m/z`,
         rt = `row retention time`,
         row_ID = `row ID`) %>% # rename to be more R friendly
  select(row_ID, everything()) # reorder and move row_ID to front

MZ_RT[1:8, 1:8]
```

    # A tibble: 8 × 8
      row_ID mz_rt            mz    rt ProcessBlank1.mzML P…¹ ProcessBlank2.mzML P…²
       <dbl> <chr>         <dbl> <dbl>                  <dbl>                  <dbl>
    1    316 85.0285_0.36…  85.0 0.364                     28                     62
    2    215 87.0077_0.35…  87.0 0.357                     47                     69
    3    685 89.0234_0.84…  89.0 0.841                    114                    160
    4    221 96.9591_0.36…  97.0 0.366                    129                    335
    5    193 111.0084_0.3… 111.  0.344                     71                     89
    6   1046 111.0084_1.8… 111.  1.89                      NA                     74
    7    412 112.0118_0.5… 112.  0.543                      7                     33
    8    502 112.0119_0.6… 112.  0.655                      7                     65
    # ℹ abbreviated names: ¹​`ProcessBlank1.mzML Peak height`,
    #   ²​`ProcessBlank2.mzML Peak height`
    # ℹ 2 more variables: `QC7.mzML Peak height` <dbl>,
    #   `QC1.mzML Peak height` <dbl>

## Clean up file names

We are using `gsub()` to replace strings (i.e. characters) in our sample
names. Here is some useful info about
[gsub](http://www.endmemo.com/r/gsub.php), and two different tutorials
[here](https://www.youtube.com/watch?v=4ZokHoF99DY) and
[here](https://www.youtube.com/watch?v=r4Sh7H6wzPA). You will likely
need to change this code to suit your purposes.

``` r
# remove stuff from the end of file names, ".mzML Peak height"
newcolumnnames <- gsub(".mzML.*","", colnames(MZ_RT))
colnames(MZ_RT) <- newcolumnnames
```

What are our sample names?

``` r
colnames(MZ_RT)
```

     [1] "row_ID"        "mz_rt"         "mz"            "rt"           
     [5] "ProcessBlank1" "ProcessBlank2" "QC7"           "QC1"          
     [9] "QC2"           "P_95_1"        "QC4"           "P_95_3"       
    [13] "P_C_4"         "P_C_2"         "P_85_3"        "P_76_3"       
    [17] "P_76_4"        "P_C_3"         "QC8"           "QC3"          
    [21] "P_76_1"        "P_85_2"        "P_76_2"        "P_95_4"       
    [25] "QC5"           "P_85_4"        "P_85_1"        "P_C_1"        
    [29] "QC6"           "P_95_2"       

## Remove unwanted samples

Removing conditioning QCs and solvent blanks.

``` r
MZ_RT <- MZ_RT %>%
  select(-contains("Conditioning"), -contains("Solvent"), -contains("CQC"))

colnames(MZ_RT)
```

     [1] "row_ID"        "mz_rt"         "mz"            "rt"           
     [5] "ProcessBlank1" "ProcessBlank2" "QC7"           "QC1"          
     [9] "QC2"           "P_95_1"        "QC4"           "P_95_3"       
    [13] "P_C_4"         "P_C_2"         "P_85_3"        "P_76_3"       
    [17] "P_76_4"        "P_C_3"         "QC8"           "QC3"          
    [21] "P_76_1"        "P_85_2"        "P_76_2"        "P_95_4"       
    [25] "QC5"           "P_85_4"        "P_85_1"        "P_C_1"        
    [29] "QC6"           "P_95_2"       

# Start filtering

### CV function

Since base R does not have a function to calculate coefficient of
variance, let’s write one.

``` r
cv <- function(x){
        (sd(x)/mean(x))
}
```

## Counting QCs

Subset QCs and filter features to keep only those that are present in
100% of QCs. You could change this parameter based on your data.

``` r
# check dimensions of current df
dim(MZ_RT)
```

    [1] 3648   30

``` r
MZ_RT_QCs <- MZ_RT %>%
  select(mz_rt, contains("QC")) %>% # select QCs
  filter(rowSums(is.na(.)) <= 1) # remove rows that have 1 or more NAs
```

``` r
# check dimensions of QCs filtered df
dim(MZ_RT_QCs)
```

    [1] 3642    9

``` r
# how many features got removed with this filtering?
nrow(MZ_RT) - nrow(MZ_RT_QCs)
```

    [1] 6

## Filter on QC CV

Here we are removing features that have a CV of more than 30% in the
QCs. The rationale is that if a feature cannot be reproducibly measured
in samples that are all the same, it should not be included in our
analysis.

``` r
# calculate CV row-wise (1 means row-wise)
QC_CV <- apply(MZ_RT_QCs[, 2:ncol(MZ_RT_QCs)], 1, cv)


# bind the CV vector back to the QC df
MZ_RT_QCs_CV <- cbind(MZ_RT_QCs, QC_CV)

# filter for keeping features with QC_CV <= 0.30 (or 30%)
MZ_RT_QCs_CVfilt <- MZ_RT_QCs_CV %>%
  filter(QC_CV <= 0.30)
```

How many features did I remove with this CV filtering?

``` r
nrow(MZ_RT_QCs) - nrow(MZ_RT_QCs_CVfilt)
```

    [1] 276

## Merge back the rest of the data

MZ_RT_QCs_CVfilt only contains the QCs, We want to keep only the rows
that are present in this df, and then merge back all of the other
samples present in MZ_RT. We will do this by creating a vector that has
the mz_rt features we want to keep, and then using `filter()` and `%in%`
to keep only features that are a part of this list.

``` r
dim(MZ_RT_QCs_CVfilt)
```

    [1] 3366   10

``` r
dim(MZ_RT)
```

    [1] 3648   30

``` r
# make a character vector of the mz_rt features we want to keep
# i.e., the ones that passed our previous filtering steps

features_to_keep <- as.character(MZ_RT_QCs_CVfilt$mz_rt)

MZ_RT_filt <- MZ_RT %>%
  filter(mz_rt %in% features_to_keep)

dim(MZ_RT_filt)
```

    [1] 3366   30

You should have the same number of features in MZ_RT_QCs_CVfilt as you
do in your new filtered df MZ_RT_filt.

``` r
all.equal(nrow(MZ_RT_QCs_CVfilt), nrow(MZ_RT_filt))
```

    [1] TRUE

## Process blanks

We want to remove features that are present in our process blanks as
they are not coming from compounds present in our samples. In this
dataset, the sample (there is only two, typically you would have at
least 3 process blanks to average) representing this process blank (a
sample that includes all the extraction materials, minus the sample,
here the juice was replaced by mass with water) has “PB” in the sample
name.

``` r
# grab the name of the column/sample that is the process blank
grep("ProcessBlank", colnames(MZ_RT_filt), value = TRUE)
```

    [1] "ProcessBlank1" "ProcessBlank2"

Calculate the average value across the QCs, then remove features that
are not at least 10x higher in the QCs than in the process blank. To do
this we will use
[`apply()`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/apply).

`apply(X, MARGIN, FUN,...)` where X is your df, MARGIN is 1 for
row-wise, and 2 for col-wise, and FUN is your function

``` r
# pull avg peak height across QCs
avg_height_QC <- apply(MZ_RT_QCs_CVfilt[, 2:ncol(MZ_RT_QCs_CVfilt)], 1, mean)

# bind back to rest of data
MZ_RT_filt_QC_avg <- cbind(MZ_RT_filt, avg_height_QC)

# check dimensions
dim(MZ_RT_filt_QC_avg)
```

    [1] 3366   31

Pull the name of your process blank, and make a new column that
indicates how many fold higher your peak height is in your average QC vs
your process blank.

``` r
# pull name of process blank 
grep("ProcessBlank", colnames(MZ_RT_filt), value = TRUE)
```

    [1] "ProcessBlank1" "ProcessBlank2"

``` r
# make a new column that has a value of how many fold higher peak height is in QCs as compared to PB
# then you can avg your PBs together and do the same thing
MZ_RT_filt_QC_avg$ProcessBlank1[is.na(MZ_RT_filt_QC_avg$ProcessBlank1)] <- 0
MZ_RT_filt_QC_avg$ProcessBlank2[is.na(MZ_RT_filt_QC_avg$ProcessBlank2)] <- 0

MZ_RT_filt_PB <- MZ_RT_filt_QC_avg %>% 
  mutate(avg_height_PB = (ProcessBlank1 + ProcessBlank2)/2) %>%
  mutate(fold_higher_in_QC = avg_height_QC/avg_height_PB) %>%
  select(row_ID, mz_rt, mz, rt, avg_height_QC, avg_height_PB, fold_higher_in_QC)

head(MZ_RT_filt_PB)
```

      row_ID           mz_rt        mz        rt avg_height_QC avg_height_PB
    1    316  85.0285_0.3638  85.02847 0.3637782      683.2315            45
    2    215  87.0077_0.3573  87.00770 0.3572833      894.2347            58
    3    685  89.0234_0.8409  89.02345 0.8408943      761.3564           137
    4    221  96.9591_0.3664  96.95914 0.3663571     3311.7877           232
    5    193 111.0084_0.3444 111.00840 0.3443712    16028.1181            80
    6   1046 111.0084_1.8865 111.00843 1.8864620    11783.5660            37
      fold_higher_in_QC
    1         15.182923
    2         15.417839
    3          5.557346
    4         14.274947
    5        200.351477
    6        318.474758

``` r
dim(MZ_RT_filt_PB)
```

    [1] 3366    7

We want to keep features that are at least 10x higher in QCs than
process blanks, and we also want to keep Infs, because an Inf indicates
that a feature absent in the process blanks (i.e., you get an Inf
because you’re trying to divide by zero).

``` r
# keep features that are present at least 10x higher in QCs vs PB
# or, keep NAs because those are absent in blank
PB_features_to_keep <- MZ_RT_filt_PB %>%
  filter(fold_higher_in_QC > 10 | is.infinite(fold_higher_in_QC)) 

dim(PB_features_to_keep)
```

    [1] 2174    7

How many features did we remove?

``` r
nrow(MZ_RT_filt_QC_avg) - nrow(PB_features_to_keep)
```

    [1] 1192

Wow removed a lot of garbage! This is great.

Bind back metdata.

``` r
MZ_RT_filt_PBremoved <- MZ_RT_filt_QC_avg %>%
  filter(mz_rt %in% PB_features_to_keep$mz_rt)
```

## Duplicates

Do we have any duplicate features?

``` r
get_dupes(MZ_RT_filt_PBremoved, mz_rt)
```

    No duplicate combinations found of: mz_rt

     [1] mz_rt         dupe_count    row_ID        mz            rt           
     [6] ProcessBlank1 ProcessBlank2 QC7           QC1           QC2          
    [11] P_95_1        QC4           P_95_3        P_C_4         P_C_2        
    [16] P_85_3        P_76_3        P_76_4        P_C_3         QC8          
    [21] QC3           P_76_1        P_85_2        P_76_2        P_95_4       
    [26] QC5           P_85_4        P_85_1        P_C_1         QC6          
    [31] P_95_2        avg_height_QC
    <0 rows> (or 0-length row.names)

No duplicates found. In case there is duplicates, let’s remove them. The
code below is not currently running.

``` r
MZ_RT_filt_PBremoved_nodupes <- MZ_RT_filt_PBremoved %>%
  distinct(mz_rt, .keep_all = TRUE)

# how many duplicates did we remove?
nrow(MZ_RT_filt_PBremoved) - nrow(MZ_RT_filt_PBremoved_nodupes)

# check for dupes again, should be none
get_dupes(MZ_RT_filt_PBremoved_nodupes, mz_rt)

dim(MZ_RT_filt_PBremoved_nodupes)
```

Remove samples that we don’t need anymore.

``` r
MZ_RT_filt_PBremoved_extraremoved <- MZ_RT_filt_PBremoved %>%
  select(-ProcessBlank1, -ProcessBlank2, -avg_height_QC)

dim(MZ_RT_filt_PBremoved_extraremoved)
```

    [1] 2174   28

# Save your file

Now you have a list of features present in your samples after filtering
for CV in QCs, and removing all the extraneous columns we added to help
us do this, along with removing any process blanks.

``` r
write_csv(MZ_RT_filt_PBremoved_extraremoved,
          "Processed_data/neg_post_filtering.csv")
```
