# 2.1 Preprocess metadata


## Purpose

This script builds `metadata.tsv`, the public metadata file every other
script reads. It takes the raw CoRIS extract, the Olink/RNA-seq ID
mappings, and the case-control pairing file, assigns pseudonymous IDs
(sid, sid_pair), cleans some variables, and computes the time-to-event
variables. Variables with re-identification risk (age, HIV transmission
route, country of origin, ART start period) are kept as empty
placeholder columns, so downstream scripts still run without those
values being exposed.

## Reproducibility

This script needs five private, controlled-access source files and can’t
be run from the public repository. It’s included so the derivation of
`metadata.tsv` is auditable, not so it can be re-executed.

``` r
library(haven)
library(lubridate)
library(readxl)
library(tidyverse)
```

## Source files

``` r
id_group_pbmc <- read_xlsx("../restricted_data/id_group_PBMC.xlsx")
id_matching <- read_xlsx("../restricted_data/id_matching_controls_cases.xlsx") # columns: ID_event, ID_control
PRETTY_metadata <- read_dta("../restricted_data/PRETTY_metadata.dta") # 187 samples, 68 variables
```

## Fix database issues

- `death_all`: `NA` means alive/no event (`death_all == 0`).
  Additionally correct a database inconsistency for participant
  `id == 9935`, who is recorded as not deceased (`death_all == 0`) being
  a control, when all the controls are recorded as `NA`.

``` r
# Use labels from Stata file in the group category
PRETTY_metadata$group <- haven::as_factor(PRETTY_metadata$group)
PRETTY_metadata |> count(group, death_all)
```

    # A tibble: 4 × 3
      group   death_all     n
      <fct>   <dbl+lbl> <int>
    1 Event    0 [No]      69
    2 Event    1 [Si]      23
    3 Control  0 [No]       1
    4 Control NA           94

``` r
PRETTY_metadata |> 
    filter(group == "Control", death_all == 0) |> 
    select(id, group, death_all)
```

    # A tibble: 1 × 3
         id group   death_all
      <dbl> <fct>   <dbl+lbl>
    1  9935 Control 0 [No]   

``` r
PRETTY_metadata <- PRETTY_metadata |>
    mutate(death_all = if_else(is.na(death_all) | id == 9935, 
                               0, death_all))
# Remove Stata labels
PRETTY_metadata$death_all <- as.numeric(PRETTY_metadata$death_all)
```

``` r
PRETTY_metadata |> count(group, death_all)
```

    # A tibble: 3 × 3
      group   death_all     n
      <fct>       <dbl> <int>
    1 Event           0    69
    2 Event           1    23
    3 Control         0    95

## Add Olink and RNA-seq sample IDs

``` r
PRETTY_metadata$olink_id <- PRETTY_metadata$id_muestra
```

``` r
id_group_pbmc <- id_group_pbmc |>
    mutate(rnaseq_id = gsub("-", "_", IDTUBO)) |>
    select(id_sample, rnaseq_id)

PRETTY_metadata <- PRETTY_metadata |>
    left_join(id_group_pbmc, by = c("id_muestra" = "id_sample"))
```

## Pseudonymize: assign `sid`

``` r
PRETTY_metadata <- PRETTY_metadata |>
    mutate(sid = paste0("PTY_", sprintf("%03d", row_number())))
```

## Matched-pair identifier (`sid_pair`)

``` r
id_matching_long <- bind_rows(
    id_matching |> mutate(id = ID_event, id_pair = ID_control,
                          .keep = "none"),
    id_matching |> mutate(id = ID_control, id_pair = ID_event,
                             .keep = "none"),
)

PRETTY_metadata <- PRETTY_metadata |>
    left_join(id_matching_long, by = "id") |>
    mutate(sid_pair = sid[match(id_pair, id)]) |>
    select(-id_pair)
```

## Time to event (`time`)

``` r
PRETTY_metadata <- PRETTY_metadata |>
    mutate(
        time = if_else(
            group == "Event",
            floor(interval(f_muestra24, fechaENOS) / months(1)),
            floor(interval(f_muestra24, art_endfu) / months(1))
        )
    )
```

## Months from ART initiation to sampling (`months_art_sampling`)

``` r
PRETTY_metadata <- PRETTY_metadata |>
    mutate(
        months_art_sampling = time_length(interval(art_ini, dmy(Fechadeobtención)), unit = "month")
    )
```

## De-identified placeholder columns

`age`, `mode` (HIV acquisition risk factor), `paisorigen` (country of
origin), and `artini_period` (year of ART initiation) carry elevated
re-identification risk. Therefore, they are included as `NA` columns so
other scripts can still run their full table structure publicly.

``` r
PRETTY_metadata <- PRETTY_metadata |>
    mutate(age = NA, mode = NA, paisorigen = NA, artini_period = NA)
```

## Recode `group`, `sex`, etc. to the labels used publicly

``` r
PRETTY_metadata <- PRETTY_metadata |>
    haven::as_factor() |>
    mutate(
        sex = fct_recode(sex, "Male" = "Hombre", "Female" = "Mujer"),
        edulvl = fct_recode(
            edulvl,
            "No studies" = "Sin estudios o primaria incompleta",
            "Primary school" = "Primaria actual completa (hasta 5º EGB)",
            "Secondary school" = "Secundaria obligatoria completa (EGB/ESO)",
            "High school" = "Bachillerato completo (BUP/COU/FP)",
            "University" = "Universitarios o superiores",
            "Other" = "5"
        ),
        art_reg = fct_recode(
            art_reg,
            "2 NRTI + 1 NNRTI" = "2NRTI+1NNRTI",
            "2 NRTI + 1 PI" = "2NRTI+1IP",
            "2 NRTI + 1 INSTI" = "2NRTI+1II",
            "Other" = "otro"
        )
    )
```

## Final public metadata file

``` r
metadata_public <- PRETTY_metadata |>
    select(
        sid, group, nadm, mace, death_all,
        age, sex, mode, paisorigen, edulvl,
        smoke, hta,
        art_reg,
        cvmax, nadircd4,
        CD4_0, CD8_0, CD4CD8_0,
        CD4_24m, CD8_24m, CD4CD8_24,
        months_art_sampling, artini_period,
        sid_pair, time
    )
```

``` r
# metadata_public |> write_tsv("../data/metadata.tsv")
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
     cellranger     1.1.0   2016-07-27 [1] CRAN (R 4.6.0)
     cli            3.6.6   2026-04-09 [1] CRAN (R 4.6.0)
     crayon         1.5.3   2024-06-20 [1] CRAN (R 4.6.0)
     digest         0.6.39  2025-11-19 [1] CRAN (R 4.6.0)
     dplyr        * 1.2.1   2026-04-03 [1] CRAN (R 4.6.0)
     evaluate       1.0.5   2025-08-27 [1] CRAN (R 4.6.0)
     farver         2.1.2   2024-05-13 [1] CRAN (R 4.6.0)
     fastmap        1.2.0   2024-05-15 [1] CRAN (R 4.6.0)
     forcats      * 1.0.1   2025-09-25 [1] CRAN (R 4.6.0)
     generics       0.1.4   2025-05-09 [1] CRAN (R 4.6.0)
     ggplot2      * 4.0.3   2026-04-22 [1] CRAN (R 4.6.0)
     glue           1.8.1   2026-04-17 [1] CRAN (R 4.6.0)
     gtable         0.3.6   2024-10-25 [1] CRAN (R 4.6.0)
     haven        * 2.5.5   2025-05-30 [1] CRAN (R 4.6.0)
     hms            1.1.4   2025-10-17 [1] CRAN (R 4.6.0)
     htmltools      0.5.9   2025-12-04 [1] CRAN (R 4.6.0)
     jsonlite       2.0.0   2025-03-27 [1] CRAN (R 4.6.0)
     knitr          1.51    2025-12-20 [1] CRAN (R 4.6.0)
     lifecycle      1.0.5   2026-01-08 [1] CRAN (R 4.6.0)
     lubridate    * 1.9.5   2026-02-04 [1] CRAN (R 4.6.0)
     magrittr       2.0.5   2026-04-04 [1] CRAN (R 4.6.0)
     otel           0.2.0   2025-08-29 [1] CRAN (R 4.6.0)
     pillar         1.11.1  2025-09-17 [1] CRAN (R 4.6.0)
     pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.6.0)
     purrr        * 1.2.2   2026-04-10 [1] CRAN (R 4.6.0)
     R6             2.6.1   2025-02-15 [1] CRAN (R 4.6.0)
     RColorBrewer   1.1-3   2022-04-03 [1] CRAN (R 4.6.0)
     readr        * 2.2.0   2026-02-19 [1] CRAN (R 4.6.0)
     readxl       * 1.5.0   2026-05-16 [1] CRAN (R 4.6.0)
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
     withr          3.0.2   2024-10-28 [1] CRAN (R 4.6.0)
     xfun           0.60    2026-07-09 [1] CRAN (R 4.6.0)
     yaml           2.3.12  2025-12-10 [1] CRAN (R 4.6.0)

     [1] /opt/Rlibs
     [2] /usr/local/lib/R/site-library
     [3] /usr/lib/R/site-library
     [4] /usr/lib/R/library
     * ── Packages attached to the search path.

    ──────────────────────────────────────────────────────────────────────────────
