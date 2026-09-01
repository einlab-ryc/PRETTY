# 2.2 Preprocess and QC Olink


## Purpose

This script runs the Olink QC pipeline on the raw NPX data. It removes
samples flagged with QC warnings, outliers, and unmatched pairs (using
`sid_pair` to drop a flagged sample’s partner too), then filters out
assays with high missingness. It outputs `olink_npx_qced_sid.tsv`, the
138-sample, 326-protein dataset used by every downstream proteomics
script, and appends the `olink_included`/`rnaseq_included` flags back
into `metadata.tsv`.

## Reproducibility

This script runs entirely from the public repository. Only
`olink_proteomics_npx.csv` and `metadata.tsv` are needed.

``` r
library(OlinkAnalyze)
library(tidyverse)
```

In the Olink report, we have the following QC relevant information, that
we will take as starting point:

- Olink Panel: Explore 384 Inflammation
- Normalization method: Intensity normalization (except PNLIPRP2)
- Number of samples: 187
- Samples passed QC: 171
- Number of proteins: 368
- Number of detected proteins: 338

## Read NPX data

This file is the same raw Olink export deposited in PRIDE. Sample
identifiers here are already the pseudonymous `sid` (`PTY_XXX`), except
for Olink’s own kit control samples.

``` r
pdata <- read_NPX("../data/olink_proteomics_npx.csv")
```

Remove Olink’s internal kit control samples (not participant data):

``` r
pdata <- pdata |> filter(!grepl("CONTROL", SampleID))
```

`SampleID` already contains `sid` for every remaining row. We keep the
`SampleID` column intact (OlinkAnalyze functions like `olink_qc_plot()`
and `olink_ttest()` require a column with this exact name internally)
and add `sid` as an alias for our own filtering/joining below.

``` r
pdata <- pdata |> mutate(sid = SampleID)
```

Add information about the groups:

``` r
mdata_groups <- read_tsv("../data/metadata.tsv") |>
    select(sid, group, nadm, mace, death_all)
pdata <- pdata |> left_join(mdata_groups, by = "sid")
```

## Number of measurements

First, we count the number of measurements in each group:

``` r
pdata |> count(group) |> mutate(pts = n / 368)
```

    # A tibble: 2 × 3
      group       n   pts
      <chr>   <int> <dbl>
    1 Control 34960    95
    2 Event   33856    92

We have 95 patients in the control group and 92 patients in the events
group. There are 3 patients that do not have a correct PS-matched pair
in the control group. We will remove these samples from further
analyses.

To identify a flagged sample’s matched partner during QC, we need the
case-control pairing structure. This is provided as the `sid_pair`
column in `metadata.tsv`.

``` r
id_matching_pairs <- read_tsv("../data/metadata.tsv") |>
    select(sid, sid_pair)
```

## Samples quality control

### Filter 1: QC Warning

Remove samples with Warning QC status:

``` r
sample_warn <- pdata |>
    filter(QC_Warning == "WARN") |>
    count(sid) |>
    mutate(percentage = round(n/n_distinct(pdata$Assay)*100)) |>
    arrange(desc(percentage))
sample_warn_ids <- sample_warn |> pull(sid)
sample_warn
```

    # A tibble: 16 × 3
       sid         n percentage
       <chr>   <int>      <dbl>
     1 PTY_081   288         78
     2 PTY_089   288         78
     3 PTY_133   288         78
     4 PTY_187   288         78
     5 PTY_156   253         69
     6 PTY_068   178         48
     7 PTY_131   178         48
     8 PTY_060   143         39
     9 PTY_059   115         31
    10 PTY_071   115         31
    11 PTY_182   115         31
    12 PTY_105   110         30
    13 PTY_169   110         30
    14 PTY_080    80         22
    15 PTY_166    80         22
    16 PTY_063    63         17

Check pairs of these samples:

``` r
sample_warn_ids_pairs <- id_matching_pairs |> filter(sid %in% sample_warn_ids) |> pull(sid_pair)
sample_warn_remove_all <- unique(c(sample_warn_ids, sample_warn_ids_pairs))
```

``` r
# Remove samples
pdata_filter1 <- pdata |> filter(!sid %in% sample_warn_remove_all)
```

This is the total number of remaining samples:

``` r
pdata_filter1 |> count(group) |> mutate(pts = n / 368)
```

    # A tibble: 2 × 3
      group       n   pts
      <chr>   <int> <dbl>
    1 Control 29440    80
    2 Event   28336    77

### Filter 2: Outliers

Identify outliers:

``` r
qc <- olink_qc_plot(pdata_filter1)
```

    `check_log` not provided. Running `check_npx()`.
    ℹ It is recommended that the user runs `check_npx()` to get a full picture of
      the results from the data validity check!

``` r
outliers <- qc$data |> filter(Outlier == 1) |> pull(SampleID)
outliers_sid <- pdata |> select(SampleID, sid) |> distinct() |> filter(SampleID %in% outliers) |> pull(sid)
```

5 outliers were identified. We are going to remove these 5 and their
paired participants.

Identify pairs of these samples:

``` r
outliers_ids_pairs <- id_matching_pairs |> filter(sid %in% outliers_sid) |> pull(sid_pair)
outliers_remove_all <- unique(c(outliers_sid, outliers_ids_pairs))
```

Remove:

``` r
pdata_filter2 <- pdata_filter1 |> filter(!sid %in% outliers_remove_all)
```

This is the total number of remaining samples:

``` r
pdata_filter2 |> count(group) |> mutate(pts = n / 368)
```

    # A tibble: 2 × 3
      group       n   pts
      <chr>   <int> <dbl>
    1 Control 27600    75
    2 Event   26496    72

### Filter 3: Unmatched samples

Identify unmatched samples (no valid PS-matched pair):

``` r
unmatched_sid <- id_matching_pairs |> filter(is.na(sid_pair)) |> pull(sid)
```

Remove:

``` r
pdata_filter3 <- pdata_filter2 |> filter(!sid %in% unmatched_sid)
```

This is the total number of remaining samples:

``` r
pdata_filter3 |> count(group) |> mutate(pts = n / 368)
```

    # A tibble: 2 × 3
      group       n   pts
      <chr>   <int> <dbl>
    1 Control 25392    69
    2 Event   25392    69

## Assays quality control

### Filter 4: Missingness

Identify assays with high missingness:

``` r
assays_todel <- pdata_filter3 |>
    mutate(Missing = ifelse(NPX > LOD, FALSE, TRUE)) |>
    group_by(Assay, group) |>
    summarise(Missing = round(mean(Missing), 4 ) > 0.25) |>
    group_by(Assay) |>
    summarise(Missing_status = all(Missing)) |>
    filter(Missing_status == TRUE) |>
    pull(Assay)
```

    `summarise()` has regrouped the output.
    ℹ Summaries were computed grouped by Assay and group.
    ℹ Output is grouped by Assay.
    ℹ Use `summarise(.groups = "drop_last")` to silence this message.
    ℹ Use `summarise(.by = c(Assay, group))` for per-operation grouping
      (`?dplyr::dplyr_by`) instead.

``` r
assays_todel
```

     [1] "ACTN4"   "AOC1"    "ARNT"    "ARTN"    "CASP2"   "CEP164"  "DGKZ"   
     [8] "EIF5A"   "FGF2"    "FKBP1B"  "FXYD5"   "GBP2"    "IL11"    "IL13"   
    [15] "IL15RA"  "IL17A"   "IL1A"    "IL2"     "IL20"    "IL22RA1" "IL24"   
    [22] "IL2RB"   "IL33"    "IL4"     "IL5"     "JCHAIN"  "JUN"     "LTO1"   
    [29] "NCLN"    "NFATC3"  "NRTN"    "PNPT1"   "PREB"    "PRKCQ"   "RAB37"  
    [36] "RGS8"    "SCRN1"   "TANK"    "TNFAIP8" "TPT1"    "VASH1"   "WAS"    

Remove:

``` r
pdata_filter4 <- pdata_filter3 |> filter(!Assay %in% assays_todel)
```

### Filter 5: Assay Warning

Count assay warnings:

``` r
pdata_filter4 |> count(Assay_Warning)
```

    # A tibble: 1 × 2
      Assay_Warning     n
      <chr>         <int>
    1 PASS          44988

There are no flagged assays.

## Final dataset

``` r
pdata_final <- pdata_filter4
```

The final dataset contains 138 samples and 326 proteins, distributed in
the following groups:

``` r
pdata_final |> select(sid, group) |> distinct() |> count(group)
```

    # A tibble: 2 × 2
      group       n
      <chr>   <int>
    1 Control    69
    2 Event      69

``` r
pdata_final |> select(sid, nadm) |> distinct() |> count(nadm)
```

    # A tibble: 2 × 2
       nadm     n
      <dbl> <int>
    1     0    91
    2     1    47

``` r
pdata_final |> select(sid, mace) |> distinct() |> count(mace)
```

    # A tibble: 2 × 2
       mace     n
      <dbl> <int>
    1     0   116
    2     1    22

``` r
pdata_final |> select(sid, death_all) |> distinct() |> count(death_all)
```

    # A tibble: 2 × 2
      death_all     n
          <dbl> <int>
    1         0   123
    2         1    15

Save to file:

``` r
# pdata_final |> write_tsv("../data/olink_npx_qced_sid.tsv")
```

Add qc columns to `metadata.tsv`:

``` r
olink_qced_sids <- pdata_final |> pull(sid) |> unique()
rnaseq_sids <- colnames(read_tsv("../data/rnaseq_transcriptomics_counts.tsv"))[-1]

metadata <- read_tsv("../data/metadata.tsv") |>
    mutate(
        olink_included = as.integer(sid %in% olink_qced_sids),
        rnaseq_included = as.integer(sid %in% rnaseq_sids)
    )

# metadata |> write_tsv("../data/metadata.tsv")
```

## Session info

``` r
sessioninfo::session_info()
```

    ─ Session info ───────────────────────────────────────────────────────────────
     setting  value
     version  R version 4.6.0 (2026-04-24)
     os       Ubuntu 26.04 LTS
     system   x86_64, linux-gnu
     ui       X11
     language (EN)
     collate  es_ES.UTF-8
     ctype    es_ES.UTF-8
     tz       Australia/Melbourne
     date     2026-08-31
     pandoc   3.8.3 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
     quarto   1.9.38 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/quarto

    ─ Packages ───────────────────────────────────────────────────────────────────
     package      * version date (UTC) lib source
     arrow          24.0.0  2026-04-29 [1] CRAN (R 4.6.0)
     assertthat     0.2.1   2019-03-21 [1] CRAN (R 4.6.0)
     bit            4.6.0   2025-03-06 [1] CRAN (R 4.6.0)
     bit64          4.8.2   2026-05-19 [1] CRAN (R 4.6.0)
     blob           1.3.0   2026-01-14 [1] CRAN (R 4.6.0)
     cli            3.6.6   2026-04-09 [1] CRAN (R 4.6.0)
     crayon         1.5.3   2024-06-20 [1] CRAN (R 4.6.0)
     data.table     1.18.4  2026-05-06 [1] CRAN (R 4.6.0)
     DBI            1.3.0   2026-02-25 [1] CRAN (R 4.6.0)
     dbplyr         2.5.2   2026-02-13 [1] CRAN (R 4.6.0)
     digest         0.6.39  2025-11-19 [1] CRAN (R 4.6.0)
     dplyr        * 1.2.1   2026-04-03 [1] CRAN (R 4.6.0)
     duckdb         1.5.2   2026-04-13 [1] CRAN (R 4.6.0)
     evaluate       1.0.5   2025-08-27 [1] CRAN (R 4.6.0)
     farver         2.1.2   2024-05-13 [1] CRAN (R 4.6.0)
     fastmap        1.2.0   2024-05-15 [1] CRAN (R 4.6.0)
     forcats      * 1.0.1   2025-09-25 [1] CRAN (R 4.6.0)
     generics       0.1.4   2025-05-09 [1] CRAN (R 4.6.0)
     ggplot2      * 4.0.3   2026-04-22 [1] CRAN (R 4.6.0)
     glue           1.8.1   2026-04-17 [1] CRAN (R 4.6.0)
     gtable         0.3.6   2024-10-25 [1] CRAN (R 4.6.0)
     hms            1.1.4   2025-10-17 [1] CRAN (R 4.6.0)
     htmltools      0.5.9   2025-12-04 [1] CRAN (R 4.6.0)
     jsonlite       2.0.0   2025-03-27 [1] CRAN (R 4.6.0)
     knitr          1.51    2025-12-20 [1] CRAN (R 4.6.0)
     lifecycle      1.0.5   2026-01-08 [1] CRAN (R 4.6.0)
     lubridate    * 1.9.5   2026-02-04 [1] CRAN (R 4.6.0)
     magrittr       2.0.5   2026-04-04 [1] CRAN (R 4.6.0)
     OlinkAnalyze * 5.0.0   2026-03-28 [1] CRAN (R 4.6.0)
     otel           0.2.0   2025-08-29 [1] CRAN (R 4.6.0)
     pillar         1.11.1  2025-09-17 [1] CRAN (R 4.6.0)
     pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.6.0)
     purrr        * 1.2.2   2026-04-10 [1] CRAN (R 4.6.0)
     R6             2.6.1   2025-02-15 [1] CRAN (R 4.6.0)
     RColorBrewer   1.1-3   2022-04-03 [1] CRAN (R 4.6.0)
     readr        * 2.2.0   2026-02-19 [1] CRAN (R 4.6.0)
     rlang          1.2.0   2026-04-06 [1] CRAN (R 4.6.0)
     rmarkdown      2.31    2026-03-26 [1] CRAN (R 4.6.0)
     rstudioapi     0.18.0  2026-01-16 [1] CRAN (R 4.6.0)
     S7             0.2.2   2026-04-22 [1] CRAN (R 4.6.0)
     scales         1.4.0   2025-04-24 [1] CRAN (R 4.6.0)
     sessioninfo    1.2.3   2025-02-05 [1] CRAN (R 4.6.0)
     stringi        1.8.7   2025-03-27 [1] CRAN (R 4.6.0)
     stringr      * 1.6.0   2025-11-04 [1] CRAN (R 4.6.0)
     tibble       * 3.3.1   2026-01-11 [1] CRAN (R 4.6.0)
     tidyr        * 1.3.2   2025-12-19 [1] CRAN (R 4.6.0)
     tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.6.0)
     tidyverse    * 2.0.0   2023-02-22 [1] CRAN (R 4.6.0)
     timechange     0.4.0   2026-01-29 [1] CRAN (R 4.6.0)
     tzdb           0.5.0   2025-03-15 [1] CRAN (R 4.6.0)
     utf8           1.2.6   2025-06-08 [1] CRAN (R 4.6.0)
     vctrs          0.7.3   2026-04-11 [1] CRAN (R 4.6.0)
     vroom          1.7.1   2026-03-31 [1] CRAN (R 4.6.0)
     withr          3.0.2   2024-10-28 [1] CRAN (R 4.6.0)
     xfun           0.60    2026-07-09 [1] CRAN (R 4.6.0)
     yaml           2.3.12  2025-12-10 [1] CRAN (R 4.6.0)

     [1] /opt/Rlibs
     [2] /usr/local/lib/R/site-library
     [3] /usr/lib/R/site-library
     [4] /usr/lib/R/library
     * ── Packages attached to the search path.

    ──────────────────────────────────────────────────────────────────────────────
