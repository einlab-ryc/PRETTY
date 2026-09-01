# 10. External validation (CNICS cohort, Abelman et al. 2026)


## Purpose

This script reproduces the external validation in the CNICS cohort,
using LGALS9 NPX values from Abelman et al. It produces Fig. 6 and Table
3, the unadjusted and adjusted hazard ratios by quartile.

## Reproducibility

This script fetches public data from the BioStudies repository at render
time (<https://doi.org/10.6019/S-BSST2240>), so an internet connection
is required, but no data is stored in this repository. The sex-adjusted
sensitivity model uses a private, non-publishable CNICS extract and is
included for transparency on the additional analysis without
redistributing that data.

``` r
library(broom)
library(ggpubr)
library(ggsurvfit)
library(gtsummary)
library(patchwork)
library(survival)
library(tidyverse)
```

``` r
biostudies_base <- "https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/240/S-BSST2240/Files"

pdata <- read_csv(file.path(biostudies_base, "npxdat_2025_10_14.csv")) |>
    dplyr::select(ID, ends_with(c("LGALS9", "CCL13", "DAPP1", "SKAP2"))) |>
    dplyr::select(ID, starts_with("npx"))

mdata <- read_csv(file.path(biostudies_base, "other_2025_10_14.csv"))
mdata2 <- read_csv(file.path(biostudies_base, "table1_2025_10_14.csv"))
```

``` r
all_data <- pdata |>
    left_join(mdata, by = "ID") |>
    left_join(mdata2, by = "ID") |>
    mutate(CD4_Nadir_200 = ifelse(CD4_Nadir < 200, 1, 0)) |>
    pivot_longer(cols = starts_with("npx"), names_to = "protein", values_to = "npx",
                 names_prefix = "npx") |>
    group_by(protein) |>
    mutate(quantile = factor(ntile(npx, 4))) |>
    ungroup()
```

``` r
gal9_data <- all_data |>
    filter(protein == "LGALS9") |>
    mutate(quantile = relevel(factor(quantile), ref = "1"))
```

## Table 3: unadjusted hazard ratios (reproducible from public data)

``` r
# Unadjusted (categorical, per-quartile HR)
gal9_cox_crude <- coxph(Surv(followup, death) ~ quantile, data = gal9_data)

# Linear-trend (one HR per quartile increase) -- for the Global P value row
gal9_trend_crude <- coxph(Surv(followup, death) ~ as.numeric(quantile), data = gal9_data)
```

``` r
n_events <- gal9_data |> group_by(quantile) |> dplyr::count(death) |> filter(death == 1)
n_events
```

    # A tibble: 4 × 3
    # Groups:   quantile [4]
      quantile death     n
      <fct>    <dbl> <int>
    1 1            1     8
    2 2            1    13
    3 3            1    26
    4 4            1    46

``` r
table3_unadjusted <- bind_rows(
    tidy(gal9_cox_crude, exponentiate = TRUE, conf.int = TRUE) |>
        filter(str_detect(term, "^quantile")) |>
        mutate(quartile = str_remove(term, "quantile"), crude_estimate = estimate,
               crude_conf.low = conf.low, crude_conf.high = conf.high, .keep = "none"),
    tibble(quartile = "1", crude_estimate = 1, crude_conf.low = NA, crude_conf.high = NA)
) |>
    left_join(
        n_events |> mutate(quartile = as.character(quantile), n_events = n, .keep = "none"),
        by = "quartile"
    ) |>
    arrange(quartile)

table3_unadjusted <- table3_unadjusted |>
    bind_rows(tibble(
        quartile = "Global P value",
        crude_estimate = tidy(gal9_trend_crude)$p.value,
        crude_conf.low = NA, crude_conf.high = NA, n_events = NA
    ))

table3_unadjusted
```

    # A tibble: 5 × 6
      quartile       crude_estimate crude_conf.low crude_conf.high quantile n_events
      <chr>                   <dbl>          <dbl>           <dbl> <fct>       <int>
    1 1               1                     NA               NA    1               8
    2 2               1.63                   0.677            3.94 2              13
    3 3               3.12                   1.41             6.90 3              26
    4 4               6.33                   2.99            13.4  4              46
    5 Global P value  0.00000000355         NA               NA    <NA>           NA

## Table 3: adjusted hazard ratios (not reproducible from public data)

``` r
gal9_data_analytic <- read_csv("../restricted_data/CNICS_sex_extract.csv") |>
    left_join(gal9_data, by = "ID")
```

    Rows: 922 Columns: 3
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (1): sex
    dbl (2): ID, MSM

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# sex has four categories: Female, FTM, Male, MTF, but we want to analyse sex
# assigned at birth so we will reduce them to Female and Male
gal9_data_analytic$sex <- case_when(
    gal9_data_analytic$sex == "FTM" ~ "Female",
    gal9_data_analytic$sex == "MTF" ~ "Male",
    TRUE ~ gal9_data_analytic$sex
)
```

``` r
# Adjusted (categorical, per-quartile HR)
gal9_cox_adj <- coxph(Surv(followup, death) ~ quantile + factor(vacs3) +
                          `Race/Ethnicity` + CD4_Nadir_200 + sex, data = gal9_data_analytic)

# Linear-trend, adjusted -- for the Global P value row
gal9_trend_adj <- coxph(Surv(followup, death) ~ as.numeric(quantile) + factor(vacs3) +
                            `Race/Ethnicity` + CD4_Nadir_200 + sex, data = gal9_data_analytic)
```

``` r
table3 <- table3_unadjusted |>
    left_join(
        bind_rows(
            tidy(gal9_cox_adj, exponentiate = TRUE, conf.int = TRUE) |>
                filter(str_detect(term, "^quantile")) |>
                mutate(quartile = str_remove(term, "quantile"), adj_estimate = estimate,
                       adj_conf.low = conf.low, adj_conf.high = conf.high, .keep = "none"),
            tibble(quartile = "1", adj_estimate = 1, adj_conf.low = NA, adj_conf.high = NA)
        ),
        by = "quartile"
    ) |>
    arrange(quartile)

table3 <- table3 |>
    filter(quartile != "Global P value") |>
    bind_rows(tibble(
        quartile = "Global P value",
        crude_estimate = tidy(gal9_trend_crude)$p.value,
        crude_conf.low = NA, crude_conf.high = NA, n_events = NA,
        adj_estimate = tidy(gal9_trend_adj)$p.value[1],
        adj_conf.low = NA, adj_conf.high = NA
    ))

table3
```

    # A tibble: 5 × 9
      quartile       crude_estimate crude_conf.low crude_conf.high quantile n_events
      <chr>                   <dbl>          <dbl>           <dbl> <fct>       <int>
    1 1               1                     NA               NA    1               8
    2 2               1.63                   0.677            3.94 2              13
    3 3               3.12                   1.41             6.90 3              26
    4 4               6.33                   2.99            13.4  4              46
    5 Global P value  0.00000000355         NA               NA    <NA>           NA
    # ℹ 3 more variables: adj_estimate <dbl>, adj_conf.low <dbl>,
    #   adj_conf.high <dbl>

``` r
# table3 |> write_tsv("../results/table3.tsv")
```

## Fig. 6

``` r
spread_y <- function(y, min_sep = 0.03, y_min = 0, y_max = 1) {
    if (length(y) <= 1) return(y)
    ord <- order(y); ys <- y[ord]
    for (i in 2:length(ys)) if (ys[i] - ys[i - 1] < min_sep) ys[i] <- ys[i - 1] + min_sep
    overflow <- max(ys) - y_max
    if (overflow > 0) {
        ys <- ys - overflow
        for (i in (length(ys) - 1):1) if (ys[i + 1] - ys[i] < min_sep) ys[i] <- ys[i + 1] - min_sep
    }
    ys <- pmin(pmax(ys, y_min), y_max)
    out <- y; out[ord] <- ys; out
}

plot_cuminc_pub <- function(fit, title, xlab = "Years", ylab = "Time-to-event prob.",
                            risktable_stats = "cum.event", x_nudge = 0.6, min_sep = 0.03,
                            right_expand = 0.2, base_size = 10, p_text = NULL) {
    p <- ggcuminc(fit) +
        labs(x = xlab, y = ylab, color = "Quartile") +
        ggtitle(title) +
        theme_bw(base_size = base_size) +
        theme(plot.title = element_text(face = "bold", hjust = 0),
              axis.title = element_text(face = "bold"), legend.position = "none")
    
    lab_df <- p$data |>
        group_by(strata) |>
        slice_max(time, n = 1, with_ties = FALSE) |>
        ungroup() |>
        mutate(strata_chr = as.character(strata), qnum = str_extract(strata_chr, "\\d+"),
               label = paste0("Q", qnum), x_label = time + x_nudge)
    lab_df <- lab_df |> mutate(y_label = spread_y(estimate, min_sep = min_sep, y_min = 0, y_max = 1))
    
    p <- p +
        geom_text(data = lab_df, aes(x = x_label, y = y_label, label = label, color = strata),
                  inherit.aes = FALSE, fontface = "bold", hjust = 0, size = 3, show.legend = FALSE) +
        coord_cartesian(clip = "off") +
        scale_x_continuous(expand = expansion(mult = c(0.02, right_expand)))
    
    if (!is.null(p_text)) {
        p <- p + annotate("text", x = Inf, y = 0.02, label = p_text, hjust = 1.1, vjust = 0,
                          size = 3.2, fontface = "plain")
    }
    p <- p + add_risktable(risktable_stats = risktable_stats)
    ggsurvfit_build(p)
}
```

``` r
fit_gal9 <- survfit2(Surv(followup, as.factor(death)) ~ quantile, data = gal9_data)

p1 <- plot_cuminc_pub(fit_gal9, "Association between LGALS9 NPX concentrations and time to event",
                      p_text = "Global P value = <0.001")
```

    Plotting outcome "1".

``` r
p2 <- gal9_data |>
    mutate(death = factor(death, labels = c("Survived", "Deceased"))) |>
    ggplot(aes(x = death, y = npx, fill = death)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    labs(x = NULL, y = "LGALS9 NPX", title = "LGALS9 plasma NPX by mortality status") +
    theme_bw() +
    theme(legend.position = "none") +
    stat_compare_means(comparisons = list(c("Survived", "Deceased")), method = "t.test")

w1 <- wrap_elements(full = p1)
w2 <- wrap_elements(full = p2)

fig6 <- (w1 | w2) +
    plot_layout(widths = c(1.6, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 12), plot.tag.position = c(0.01, 0.99))

# svg("../results/fig6_CNICS_validation.svg", width = 10, height = 7)
fig6
```

![](10_validation_cnics_files/figure-commonmark/unnamed-chunk-12-1.png)

``` r
# dev.off()
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
     date     2026-09-01
     pandoc   3.8.3 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
     quarto   1.9.38 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/quarto

    ─ Packages ───────────────────────────────────────────────────────────────────
     package      * version date (UTC) lib source
     abind          1.4-8   2024-09-12 [1] CRAN (R 4.6.0)
     backports      1.5.1   2026-04-03 [1] CRAN (R 4.6.0)
     bit            4.6.0   2025-03-06 [1] CRAN (R 4.6.0)
     bit64          4.8.2   2026-05-19 [1] CRAN (R 4.6.0)
     broom        * 1.0.13  2026-05-14 [1] CRAN (R 4.6.0)
     car            3.1-5   2026-02-03 [1] CRAN (R 4.6.0)
     carData        3.0-6   2026-01-30 [1] CRAN (R 4.6.0)
     cli            3.6.6   2026-04-09 [1] CRAN (R 4.6.0)
     crayon         1.5.3   2024-06-20 [1] CRAN (R 4.6.0)
     curl           7.1.0   2026-04-22 [1] CRAN (R 4.6.0)
     digest         0.6.39  2025-11-19 [1] CRAN (R 4.6.0)
     dplyr        * 1.2.1   2026-04-03 [1] CRAN (R 4.6.0)
     evaluate       1.0.5   2025-08-27 [1] CRAN (R 4.6.0)
     farver         2.1.2   2024-05-13 [1] CRAN (R 4.6.0)
     fastmap        1.2.0   2024-05-15 [1] CRAN (R 4.6.0)
     forcats      * 1.0.1   2025-09-25 [1] CRAN (R 4.6.0)
     Formula        1.2-5   2023-02-24 [1] CRAN (R 4.6.0)
     generics       0.1.4   2025-05-09 [1] CRAN (R 4.6.0)
     ggplot2      * 4.0.3   2026-04-22 [1] CRAN (R 4.6.0)
     ggpubr       * 0.6.3   2026-02-24 [1] CRAN (R 4.6.0)
     ggsignif       0.6.4   2022-10-13 [1] CRAN (R 4.6.0)
     ggsurvfit    * 1.2.0   2025-09-13 [1] CRAN (R 4.6.0)
     glue           1.8.1   2026-04-17 [1] CRAN (R 4.6.0)
     gtable         0.3.6   2024-10-25 [1] CRAN (R 4.6.0)
     gtsummary    * 2.5.0   2025-12-05 [1] CRAN (R 4.6.0)
     hms            1.1.4   2025-10-17 [1] CRAN (R 4.6.0)
     htmltools      0.5.9   2025-12-04 [1] CRAN (R 4.6.0)
     jsonlite       2.0.0   2025-03-27 [1] CRAN (R 4.6.0)
     knitr          1.51    2025-12-20 [1] CRAN (R 4.6.0)
     labeling       0.4.3   2023-08-29 [1] CRAN (R 4.6.0)
     lattice        0.22-9  2026-02-09 [1] CRAN (R 4.6.0)
     lifecycle      1.0.5   2026-01-08 [1] CRAN (R 4.6.0)
     lubridate    * 1.9.5   2026-02-04 [1] CRAN (R 4.6.0)
     magrittr       2.0.5   2026-04-04 [1] CRAN (R 4.6.0)
     Matrix         1.7-5   2026-03-21 [4] CRAN (R 4.5.3)
     otel           0.2.0   2025-08-29 [1] CRAN (R 4.6.0)
     patchwork    * 1.3.2   2025-08-25 [1] CRAN (R 4.6.0)
     pillar         1.11.1  2025-09-17 [1] CRAN (R 4.6.0)
     pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.6.0)
     purrr        * 1.2.2   2026-04-10 [1] CRAN (R 4.6.0)
     R6             2.6.1   2025-02-15 [1] CRAN (R 4.6.0)
     RColorBrewer   1.1-3   2022-04-03 [1] CRAN (R 4.6.0)
     readr        * 2.2.0   2026-02-19 [1] CRAN (R 4.6.0)
     rlang          1.2.0   2026-04-06 [1] CRAN (R 4.6.0)
     rmarkdown      2.31    2026-03-26 [1] CRAN (R 4.6.0)
     rstatix        0.7.3   2025-10-18 [1] CRAN (R 4.6.0)
     rstudioapi     0.18.0  2026-01-16 [1] CRAN (R 4.6.0)
     S7             0.2.2   2026-04-22 [1] CRAN (R 4.6.0)
     scales         1.4.0   2025-04-24 [1] CRAN (R 4.6.0)
     sessioninfo    1.2.3   2025-02-05 [1] CRAN (R 4.6.0)
     stringi        1.8.7   2025-03-27 [1] CRAN (R 4.6.0)
     stringr      * 1.6.0   2025-11-04 [1] CRAN (R 4.6.0)
     survival     * 3.8-6   2026-01-16 [4] CRAN (R 4.5.2)
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
