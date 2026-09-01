# 4. Differential expression analysis of proteins


## Purpose

This script runs differential expression analysis on the Olink
proteomics data across the four outcomes (SNAE, malignancy,
cardiovascular event, mortality) using `olink_ttest`. It produces the
`dep_*.csv` files used by the mixomics and correlation scripts, and
generates the FDR-adjusted volcano plots (Fig. 1).

## Reproducibility

This script runs entirely from the public repository. Only
`olink_npx_qced_sid.tsv` is needed.

``` r
library(ggpubr)
library(ggrepel)
library(OlinkAnalyze)
library(tidyverse)
```

## Read data

``` r
pdata <- read_tsv("../data/olink_npx_qced_sid.tsv") |>
    mutate(
        group = factor(group, levels = c("Event", "Control")),
        nadm = factor(nadm, levels = c(1, 0)),
        mace = factor(mace, levels = c(1, 0)),
        death_all = factor(death_all, levels = c(1, 0))
    )
```

## Shared helpers

Colours/alphas and the differential-expression classification are
identical across every outcome. The FDR volcano plot is also shared, but
each outcome has small differences from the others (whether the FDR=0.05
reference line is drawn, y-axis expansion, label nudge/force, and
whether one assay’s label needs a separately-nudged position to avoid
overlap). Those are exposed as parameters rather than hardcoded, so each
call site states exactly how it differs.

``` r
cols   <- c("Up" = "#ffad73", "Down" = "#26b3ff", "n.s." = "grey")
alphas <- c("Up" = 1, "Down" = 1, "n.s." = 0.5)
highlight_assays <- c("LGALS9", "DAPP1", "SKAP2", "CCL13")

compute_dep <- function(pdata, outcome_var) {
    olink_ttest(pdata, outcome_var) |>
        mutate(prot_exp = case_when(
            estimate >= 0 & Adjusted_pval <= 0.05 ~ "Up",
            estimate <= 0 & Adjusted_pval <= 0.05 ~ "Down",
            TRUE ~ "n.s."
        ))
}

plot_volcano_fdr <- function(dep_df, title, xlim_range,
                             show_hline = TRUE,
                             y_expand = NULL,
                             label_nudge_y = 0.05,
                             label_force = 15,
                             split_label = NULL) {
    p <- dep_df |>
        ggplot(aes(x = estimate, y = -log10(Adjusted_pval), fill = prot_exp, alpha = prot_exp)) +
        geom_point(shape = 21, colour = "black")
    
    if (show_hline) {
        p <- p + geom_hline(yintercept = -log10(0.05), linetype = "dashed")
    }
    if (!is.null(y_expand)) {
        p <- p + scale_y_continuous(expand = expansion(mult = y_expand))
    }
    
    # Main label layer -- excludes split_label$assay when one assay needs a
    # separately-nudged position (e.g. SKAP2 in the NADM plot).
    label_data <- if (!is.null(split_label)) dep_df |> filter(Assay != split_label$assay) else dep_df
    p <- p +
        geom_label_repel(
            data = label_data,
            aes(label = if_else(Assay %in% highlight_assays, Assay, "")),
            size = 2, show.legend = FALSE, max.overlaps = Inf,
            nudge_y = label_nudge_y, min.segment.length = 0,
            force = label_force, label.padding = unit(0.1, "lines")
        )
    
    if (!is.null(split_label)) {
        p <- p +
            geom_label_repel(
                data = dep_df |> filter(Assay == split_label$assay),
                aes(label = if_else(Assay %in% highlight_assays, Assay, "")),
                size = 2, show.legend = FALSE, max.overlaps = Inf,
                nudge_y = split_label$nudge_y, nudge_x = split_label$nudge_x,
                min.segment.length = 0, force = 15, label.padding = unit(0.1, "lines")
            )
    }
    
    p +
        scale_fill_manual(values = cols, breaks = c("Up", "Down", "n.s."),
                          labels = c("Up", "Down", "n.s.")) +
        scale_alpha_manual(values = alphas) +
        xlim(xlim_range) +
        theme_bw() +
        theme(legend.position = "bottom", legend.background = element_rect(colour = "black"),
              axis.title = element_text(size = 9)) +
        labs(fill = "Differential abundance", title = title,
             x = expression(log[2](FC)), y = expression(-log[10]("FDR-adjusted P value"))) +
        guides(alpha = "none", fill = guide_legend(override.aes = list(alpha = 1)))
}
```

## Group

``` r
dep_group <- compute_dep(pdata, "group")
```

    `check_log` not provided. Running `check_npx()`.
    ℹ It is recommended that the user runs `check_npx()` to get a full picture of
      the results from the data validity check!
    T-test is performed on Event - Control.

``` r
dep_group
```

    # A tibble: 326 × 17
       Assay    OlinkID  UniProt Panel    estimate   Event Control statistic p.value
       <chr>    <chr>    <chr>   <chr>       <dbl>   <dbl>   <dbl>     <dbl>   <dbl>
     1 MAP2K6   OID20552 P52564  Inflamm…   -0.744 -0.327   0.417      -3.68 3.34e-4
     2 CLIP2    OID20559 Q9UDT6  Inflamm…   -1.04  -0.452   0.587      -3.63 4.08e-4
     3 MGLL     OID20711 Q99685  Inflamm…   -0.823 -0.368   0.455      -3.33 1.13e-3
     4 IKBKG    OID20544 Q9Y6K9  Inflamm…   -0.638 -0.292   0.346      -3.22 1.62e-3
     5 DBNL     OID20681 Q9UJU6  Inflamm…   -0.813 -0.314   0.498      -3.11 2.27e-3
     6 BANK1    OID20594 Q8NDB2  Inflamm…   -0.823 -0.374   0.450      -3.10 2.37e-3
     7 LGALS9   OID20781 O00182  Inflamm…    0.328  0.282  -0.0458      3.09 2.48e-3
     8 EDAR     OID20525 Q9UNE0  Inflamm…   -0.396 -0.0398  0.356      -3.06 2.68e-3
     9 ARHGEF12 OID20554 Q9NZN5  Inflamm…   -0.638 -0.0409  0.597      -3.05 2.79e-3
    10 TRIM5    OID20526 Q9C035  Inflamm…   -0.426 -0.0741  0.352      -3.02 3.06e-3
    # ℹ 316 more rows
    # ℹ 8 more variables: parameter <dbl>, conf.low <dbl>, conf.high <dbl>,
    #   method <chr>, alternative <chr>, Adjusted_pval <dbl>, Threshold <chr>,
    #   prot_exp <chr>

``` r
# dep_group |> write_tsv("../results/tableS2_dep_group.tsv")
```

``` r
p1 <- plot_volcano_fdr(dep_group, "PWH with vs without SNAEs", xlim_range = c(-1.1, 1.1),
                       y_expand = c(0.02, 0.15), label_nudge_y = 0.1, label_force = 1)
p1
```

![](04_prots_diffexp_files/figure-commonmark/unnamed-chunk-5-1.png)

## NADM

``` r
dep_nadm <- compute_dep(pdata, "nadm")
```

    `check_log` not provided. Running `check_npx()`.
    ℹ It is recommended that the user runs `check_npx()` to get a full picture of
      the results from the data validity check!
    T-test is performed on 1 - 0.

``` r
# dep_nadm |> write_tsv("../results/tableS3_dep_nadm.tsv")
```

``` r
p2 <- plot_volcano_fdr(dep_nadm, "PWH with vs without neoplasias", xlim_range = c(-1.7, 1.7),
                       split_label = list(assay = "SKAP2", nudge_y = 0.08, nudge_x = -0.3))
p2
```

![](04_prots_diffexp_files/figure-commonmark/unnamed-chunk-7-1.png)

## MACE

``` r
dep_mace <- compute_dep(pdata, "mace")
```

    `check_log` not provided. Running `check_npx()`.
    ℹ It is recommended that the user runs `check_npx()` to get a full picture of
      the results from the data validity check!
    T-test is performed on 1 - 0.

``` r
# dep_mace |> write_tsv("../results/tableS4_dep_mace.tsv")
```

``` r
p3 <- plot_volcano_fdr(dep_mace, "PWH with vs without CV events", xlim_range = c(-1.2, 1.2),
                       show_hline = FALSE)
p3
```

![](04_prots_diffexp_files/figure-commonmark/unnamed-chunk-9-1.png)

## Death_all

``` r
dep_death_all <- compute_dep(pdata, "death_all")
```

    `check_log` not provided. Running `check_npx()`.
    ℹ It is recommended that the user runs `check_npx()` to get a full picture of
      the results from the data validity check!
    T-test is performed on 1 - 0.

``` r
# dep_death_all |> write_tsv("../results/tableS5_dep_death_all.tsv")
```

``` r
p4 <- plot_volcano_fdr(dep_death_all, "PWH with vs without mortality", xlim_range = c(-2.8, 2.8))
p4
```

![](04_prots_diffexp_files/figure-commonmark/unnamed-chunk-11-1.png)

``` r
# svg("../results/fig1_AtoD.svg", width = 7, height = 5)
ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"),
          common.legend = TRUE, legend = "bottom")
```

![](04_prots_diffexp_files/figure-commonmark/unnamed-chunk-12-1.png)

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
     date     2026-08-31
     pandoc   3.8.3 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
     quarto   1.9.38 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/quarto

    ─ Packages ───────────────────────────────────────────────────────────────────
     package      * version   date (UTC) lib source
     abind          1.4-8     2024-09-12 [1] CRAN (R 4.6.0)
     arrow          24.0.0    2026-04-29 [1] CRAN (R 4.6.0)
     assertthat     0.2.1     2019-03-21 [1] CRAN (R 4.6.0)
     backports      1.5.1     2026-04-03 [1] CRAN (R 4.6.0)
     bit            4.6.0     2025-03-06 [1] CRAN (R 4.6.0)
     bit64          4.8.2     2026-05-19 [1] CRAN (R 4.6.0)
     blob           1.3.0     2026-01-14 [1] CRAN (R 4.6.0)
     broom          1.0.13    2026-05-14 [1] CRAN (R 4.6.0)
     car            3.1-5     2026-02-03 [1] CRAN (R 4.6.0)
     carData        3.0-6     2026-01-30 [1] CRAN (R 4.6.0)
     cli            3.6.6     2026-04-09 [1] CRAN (R 4.6.0)
     cowplot        1.2.0     2025-07-07 [1] CRAN (R 4.6.0)
     crayon         1.5.3     2024-06-20 [1] CRAN (R 4.6.0)
     DBI            1.3.0     2026-02-25 [1] CRAN (R 4.6.0)
     dbplyr         2.5.2     2026-02-13 [1] CRAN (R 4.6.0)
     digest         0.6.39    2025-11-19 [1] CRAN (R 4.6.0)
     dplyr        * 1.2.1     2026-04-03 [1] CRAN (R 4.6.0)
     duckdb         1.5.2     2026-04-13 [1] CRAN (R 4.6.0)
     evaluate       1.0.5     2025-08-27 [1] CRAN (R 4.6.0)
     farver         2.1.2     2024-05-13 [1] CRAN (R 4.6.0)
     fastmap        1.2.0     2024-05-15 [1] CRAN (R 4.6.0)
     forcats      * 1.0.1     2025-09-25 [1] CRAN (R 4.6.0)
     Formula        1.2-5     2023-02-24 [1] CRAN (R 4.6.0)
     generics       0.1.4     2025-05-09 [1] CRAN (R 4.6.0)
     ggplot2      * 4.0.3     2026-04-22 [1] CRAN (R 4.6.0)
     ggpubr       * 0.6.3     2026-02-24 [1] CRAN (R 4.6.0)
     ggrepel      * 0.9.8     2026-03-17 [1] CRAN (R 4.6.0)
     ggsignif       0.6.4     2022-10-13 [1] CRAN (R 4.6.0)
     glue           1.8.1     2026-04-17 [1] CRAN (R 4.6.0)
     gridExtra      2.3       2017-09-09 [1] CRAN (R 4.6.0)
     gtable         0.3.6     2024-10-25 [1] CRAN (R 4.6.0)
     hms            1.1.4     2025-10-17 [1] CRAN (R 4.6.0)
     htmltools      0.5.9     2025-12-04 [1] CRAN (R 4.6.0)
     jsonlite       2.0.0     2025-03-27 [1] CRAN (R 4.6.0)
     knitr          1.51      2025-12-20 [1] CRAN (R 4.6.0)
     labeling       0.4.3     2023-08-29 [1] CRAN (R 4.6.0)
     lifecycle      1.0.5     2026-01-08 [1] CRAN (R 4.6.0)
     lubridate    * 1.9.5     2026-02-04 [1] CRAN (R 4.6.0)
     magrittr       2.0.5     2026-04-04 [1] CRAN (R 4.6.0)
     OlinkAnalyze * 5.0.0     2026-03-28 [1] CRAN (R 4.6.0)
     otel           0.2.0     2025-08-29 [1] CRAN (R 4.6.0)
     pillar         1.11.1    2025-09-17 [1] CRAN (R 4.6.0)
     pkgconfig      2.0.3     2019-09-22 [1] CRAN (R 4.6.0)
     purrr        * 1.2.2     2026-04-10 [1] CRAN (R 4.6.0)
     R6             2.6.1     2025-02-15 [1] CRAN (R 4.6.0)
     RColorBrewer   1.1-3     2022-04-03 [1] CRAN (R 4.6.0)
     Rcpp           1.1.1-1.1 2026-04-24 [1] CRAN (R 4.6.0)
     readr        * 2.2.0     2026-02-19 [1] CRAN (R 4.6.0)
     rlang          1.2.0     2026-04-06 [1] CRAN (R 4.6.0)
     rmarkdown      2.31      2026-03-26 [1] CRAN (R 4.6.0)
     rstatix        0.7.3     2025-10-18 [1] CRAN (R 4.6.0)
     rstudioapi     0.18.0    2026-01-16 [1] CRAN (R 4.6.0)
     S7             0.2.2     2026-04-22 [1] CRAN (R 4.6.0)
     scales         1.4.0     2025-04-24 [1] CRAN (R 4.6.0)
     sessioninfo    1.2.3     2025-02-05 [1] CRAN (R 4.6.0)
     stringi        1.8.7     2025-03-27 [1] CRAN (R 4.6.0)
     stringr      * 1.6.0     2025-11-04 [1] CRAN (R 4.6.0)
     tibble       * 3.3.1     2026-01-11 [1] CRAN (R 4.6.0)
     tidyr        * 1.3.2     2025-12-19 [1] CRAN (R 4.6.0)
     tidyselect     1.2.1     2024-03-11 [1] CRAN (R 4.6.0)
     tidyverse    * 2.0.0     2023-02-22 [1] CRAN (R 4.6.0)
     timechange     0.4.0     2026-01-29 [1] CRAN (R 4.6.0)
     tzdb           0.5.0     2025-03-15 [1] CRAN (R 4.6.0)
     utf8           1.2.6     2025-06-08 [1] CRAN (R 4.6.0)
     vctrs          0.7.3     2026-04-11 [1] CRAN (R 4.6.0)
     vroom          1.7.1     2026-03-31 [1] CRAN (R 4.6.0)
     withr          3.0.2     2024-10-28 [1] CRAN (R 4.6.0)
     xfun           0.60      2026-07-09 [1] CRAN (R 4.6.0)
     yaml           2.3.12    2025-12-10 [1] CRAN (R 4.6.0)

     [1] /opt/Rlibs
     [2] /usr/local/lib/R/site-library
     [3] /usr/lib/R/site-library
     [4] /usr/lib/R/library
     * ── Packages attached to the search path.

    ──────────────────────────────────────────────────────────────────────────────
