# 6. ELISA validation panels (Figure 4)


## Purpose

This script builds the ELISA validation panels for the four biomarkers,
LGALS9, CCL13, DAPP1, and SKAP2. It compares cases against controls,
compares clinical event subtypes, compares against death, and correlates
each biomarker with its Olink NPX value. It produces Fig. 4 (all
measurements) and Fig. S4 (sensitivity analysis with Tukey outliers
excluded).

## Reproducibility

This script runs entirely from the public repository. Only
`elisa_validation.tsv`, `olink_proteomics_npx.csv`, and `metadata.tsv`
are needed.

``` r
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
```

``` r
BASE_FONT    <- 8    # legible at the 21x14.85cm final print size (iterated by rendering)
TITLE_SIZE   <- BASE_FONT * 0.85
BRACKET_SIZE <- BASE_FONT * 0.28
POINT_SIZE   <- 0.7

# Colour-blind-friendly palette (Okabe & Ito 2008)
CB_COLORS <- c(
    Control                = "#0072B2",  # blue
    Cases                  = "#D55E00",  # vermillion
    Cancer                 = "#E69F00",  # orange
    Cardiovascular         = "#009E73",  # bluish green
    `Death (all-cause)`    = "#CC79A7"   # reddish purple
)
```

``` r
format_p <- function(p, digits = 4) {
    threshold <- 10^(-digits)
    ifelse(p < threshold, sprintf("<%s", format(threshold, scientific = FALSE)),
           sprintf("%.4f", p))
}
```

## Load data

``` r
elisa <- read_tsv("../data/elisa_validation.tsv")
olink <- read_csv("../data/olink_proteomics_npx.csv")
metadata <- read_tsv("../data/metadata.tsv") |>
    mutate(group = if_else(group == "Event", 1, 0))
```

One row per protein x sample (unique per sample – no group expansion, so
no double-counting anywhere downstream).

``` r
protein_long <- bind_rows(lapply(c("DAPP1", "CCL13", "SKAP2", "LGALS9"), function(p) {
    conc <- elisa |> filter(assay == p) |> dplyr::select(sid, conc = mean_conc)
    npx  <- olink |> filter(Assay == p) |> dplyr::select(sid = SampleID, npx = NPX)
    
    conc |>
        inner_join(npx, by = "sid") |>
        inner_join(metadata |> dplyr::select(sid, group, nadm, mace, death_all), by = "sid") |>
        mutate(protein = p) |>
        filter(!is.na(conc))
}))
protein_dfs <- split(protein_long, protein_long$protein)
```

``` r
flag_iqr <- function(x) {
    q <- quantile(x, c(0.25, 0.75), na.rm = TRUE); iqr <- diff(q)
    x < (q[1] - 1.5 * iqr) | x > (q[2] + 1.5 * iqr)
}

# Outlier-excluded, deduplicated per-sample datasets -- one exclusion rule per
# protein applied once, not per comparison (the panel comparisons below share
# overlapping membership by design). Same Tukey 1.5xIQR rule for all 4
# proteins -- no protein-specific exception.
scatter_dfs <- protein_dfs
scatter_dfs_clean <- lapply(protein_dfs, function(d) d %>% filter(!flag_iqr(conc)))
```

## ELISA vs. NPX correlation, all 4 proteins

``` r
cor_tab <- bind_rows(lapply(names(protein_dfs), function(p) {
    d <- protein_dfs[[p]]
    ct <- cor.test(d$conc, d$npx, method = "spearman")
    tibble(protein = p, n = nrow(d), spearman_rho = unname(ct$estimate), p_value = ct$p.value)
}))
```

    Warning in cor.test.default(d$conc, d$npx, method = "spearman"): cannot compute
    exact p-value with ties
    Warning in cor.test.default(d$conc, d$npx, method = "spearman"): cannot compute
    exact p-value with ties

``` r
cor_tab
```

    # A tibble: 4 × 4
      protein     n spearman_rho  p_value
      <chr>   <int>        <dbl>    <dbl>
    1 CCL13     136       0.657  3.92e-18
    2 DAPP1     137       0.0706 4.12e- 1
    3 LGALS9    137       0.451  4.63e- 8
    4 SKAP2     136      -0.0544 5.29e- 1

## Column definitions

``` r
# Column 1: Cases (cancer, cardiovascular event, or non-accidental death) vs. Control.
make_cases_df <- function(df) {
    df %>% filter(group %in% c(0, 1)) %>%
        mutate(clinical_group = factor(ifelse(group == 0, "Control", "Cases"),
                                       levels = c("Control", "Cases")))
}

# Column 2: Cancer vs. Cardiovascular vs. Control -- a genuine 3-way partition.
make_component_df <- function(df) {
    df %>%
        mutate(clinical_group = case_when(
            group == 0             ~ "Control",
            nadm == 1 & mace == 0  ~ "Cancer",
            mace == 1 & nadm == 0  ~ "Cardiovascular",
            TRUE                   ~ NA_character_
        )) %>%
        filter(!is.na(clinical_group)) %>%
        mutate(clinical_group = factor(clinical_group,
                                       levels = c("Control", "Cancer", "Cardiovascular")))
}

# Column 3: Death (all-cause) vs. Control.
make_death_df <- function(df) {
    df %>% filter(group == 0 | death_all == 1) %>%
        mutate(clinical_group = factor(ifelse(group == 0, "Control", "Death (all-cause)"),
                                       levels = c("Control", "Death (all-cause)")))
}

comparison_builders <- list(
    Cases     = make_cases_df,
    Component = make_component_df,
    Death     = make_death_df
)
comparison_headers <- c(
    Cases     = "Cases vs. Control",
    Component = "Cancer vs. Cardiovascular vs. Control",
    Death     = "Death (all-cause) vs. Control",
    Scatter   = "ELISA vs. Olink NPX"
)
```

## Panel figures

``` r
make_two_group_violin <- function(df2, protein_label, unit, comp_key, show_header) {
    event_label <- setdiff(levels(df2$clinical_group), "Control")
    a <- df2 %>% filter(clinical_group == "Control")   %>% pull(conc)
    b <- df2 %>% filter(clinical_group == event_label)  %>% pull(conc)
    pw <- if (length(a) > 0 && length(b) > 0) wilcox.test(a, b) else list(p.value = NA_real_)
    sig <- !is.na(pw$p.value) && pw$p.value < 0.05
    top_expand <- if (sig) 0.32 else 0.12
    
    title_txt <- if (show_header) {
        sprintf("%s\nP=%s", comparison_headers[[comp_key]], format_p(pw$p.value))
    } else {
        sprintf("P=%s", format_p(pw$p.value))
    }
    p <- ggplot(df2, aes(x = clinical_group, y = conc, fill = clinical_group)) +
        geom_violin(alpha = 0.6, trim = FALSE, linewidth = 0.25) +
        geom_jitter(width = 0.08, alpha = 0.5, size = POINT_SIZE) +
        scale_y_log10(expand = expansion(mult = c(0.05, top_expand))) +
        scale_fill_manual(values = CB_COLORS) +
        labs(title = title_txt, x = NULL, y = sprintf("%s (%s, log10)", protein_label, unit)) +
        theme_minimal(base_size = BASE_FONT) +
        theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1),
              plot.title = element_text(size = TITLE_SIZE))
    
    if (sig) {
        p <- p + stat_compare_means(comparisons = list(c("Control", event_label)), method = "wilcox.test",
                                    label = "p.format", tip.length = 0.01, size = BRACKET_SIZE)
    }
    p
}

make_three_group_violin <- function(df3, protein_label, unit, show_header) {
    kw <- kruskal.test(conc ~ clinical_group, data = df3)
    levs <- levels(droplevels(df3$clinical_group))
    pairs <- combn(levs, 2, simplify = FALSE)
    pvals <- vapply(pairs, function(pr) {
        a <- df3$conc[df3$clinical_group == pr[1]]
        b <- df3$conc[df3$clinical_group == pr[2]]
        if (length(a) > 1 && length(b) > 1) wilcox.test(a, b)$p.value else NA_real_
    }, numeric(1))
    sig_pairs <- pairs[!is.na(pvals) & pvals < 0.05]
    # Narrowest-span comparison first (drawn lowest), widest last (drawn
    # highest) -- with only 2-3 possible pairs from 3 groups this keeps stacked
    # brackets/labels from overlapping (ggpubr stacks in list order).
    sig_pairs <- sig_pairs[order(vapply(sig_pairs, function(pr) diff(match(pr, levs)), numeric(1)))]
    top_expand <- 0.12 + 0.28 * length(sig_pairs)
    
    title_txt <- if (show_header) {
        sprintf("%s\nP=%s", comparison_headers[["Component"]], format_p(kw$p.value))
    } else {
        sprintf("P=%s", format_p(kw$p.value))
    }
    
    p <- ggplot(df3, aes(x = clinical_group, y = conc, fill = clinical_group)) +
        geom_violin(alpha = 0.6, trim = FALSE, linewidth = 0.25) +
        geom_jitter(width = 0.08, alpha = 0.5, size = POINT_SIZE) +
        scale_y_log10(expand = expansion(mult = c(0.05, top_expand))) +
        scale_fill_manual(values = CB_COLORS) +
        labs(title = title_txt, x = NULL, y = sprintf("%s (%s, log10)", protein_label, unit)) +
        theme_minimal(base_size = BASE_FONT) +
        theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1),
              plot.title = element_text(size = TITLE_SIZE))
    
    if (length(sig_pairs) > 0) {
        p <- p + stat_compare_means(comparisons = sig_pairs, method = "wilcox.test",
                                    label = "p.format", tip.length = 0.01, size = BRACKET_SIZE,
                                    step.increase = 0.18)
    }
    p
}

make_scatter <- function(df, protein_label, show_header) {
    ct <- cor.test(df$conc, df$npx, method = "spearman")
    title_txt <- if (show_header) {
        sprintf("%s\nrho=%.3f, P=%s", comparison_headers[["Scatter"]], ct$estimate, format_p(ct$p.value))
    } else {
        sprintf("rho=%.3f, P=%s", ct$estimate, format_p(ct$p.value))
    }
    ggplot(df, aes(x = npx, y = conc)) +
        geom_point(size = POINT_SIZE, alpha = 0.6) +
        labs(title = title_txt, x = "NPX", y = "Concentration") +
        theme_minimal(base_size = BASE_FONT) +
        theme(plot.title = element_text(size = TITLE_SIZE))
}
```

## Build figure

``` r
unit_map <- c(SKAP2 = "ng/ml", DAPP1 = "ng/ml", LGALS9 = "pg/ml", CCL13 = "pg/ml")
panel_order <- c("SKAP2", "DAPP1", "LGALS9", "CCL13")  # published panel order A-D

# One row of 4 panels (Cases, Cancer-vs-Cardio-vs-Control, Death, scatter) per
# protein. Column headers appear once, on the first protein's row only.
build_row <- function(protein_label, base_df, scatter_df, show_header) {
    cases_df     <- comparison_builders$Cases(base_df)
    component_df <- comparison_builders$Component(base_df)
    death_df     <- comparison_builders$Death(base_df)
    list(
        make_two_group_violin(cases_df, protein_label, unit_map[[protein_label]], "Cases", show_header),
        make_three_group_violin(component_df, protein_label, unit_map[[protein_label]], show_header),
        make_two_group_violin(death_df, protein_label, unit_map[[protein_label]], "Death", show_header),
        make_scatter(scatter_df, protein_label, show_header)
    )
}
```

``` r
figures <- list()
for (version in c("with_outliers", "no_outliers")) {
    dfs <- if (version == "with_outliers") scatter_dfs else scatter_dfs_clean
    panels <- unlist(lapply(seq_along(panel_order), function(i) {
        p <- panel_order[i]
        build_row(p, dfs[[p]], dfs[[p]], show_header = (i == 1))
    }), recursive = FALSE)
    fig <- wrap_plots(panels, ncol = 4)
    figures[[version]] <- fig
}
```

    Warning in cor.test.default(df$conc, df$npx, method = "spearman"): cannot
    compute exact p-value with ties
    Warning in cor.test.default(df$conc, df$npx, method = "spearman"): cannot
    compute exact p-value with ties
    Warning in cor.test.default(df$conc, df$npx, method = "spearman"): cannot
    compute exact p-value with ties
    Warning in cor.test.default(df$conc, df$npx, method = "spearman"): cannot
    compute exact p-value with ties

Fig. 4 (all measurements, no exclusions):

``` r
figures$with_outliers
```

![](06_fig4_elisa_panels_files/figure-commonmark/unnamed-chunk-12-1.png)

``` r
# Exported at half a DIN-A4 sheet (21 x 14.85cm) -- no further scaling
# should be applied when placing these files in a document.
# ggsave("../results/fig4_with_outliers.svg", plot = figures$with_outliers, width = 21, height = 14.85, units = "cm")
# ggsave("../results/fig4_with_outliers.pdf", plot = figures$with_outliers, width = 21, height = 14.85, units = "cm")
```

Fig. S4 (sensitivity analysis, Tukey outliers excluded):

``` r
figures$no_outliers
```

![](06_fig4_elisa_panels_files/figure-commonmark/unnamed-chunk-13-1.png)

``` r
# Exported at half a DIN-A4 sheet (21 x 14.85cm) -- no further scaling
# should be applied when placing these files in a document.
# ggsave("../results/figS4_no_outliers.svg", plot = figures$no_outliers, width = 21, height = 14.85, units = "cm")
# ggsave("../results/figS4_no_outliers.pdf", plot = figures$no_outliers, width = 21, height = 14.85, units = "cm")
```

## Full comparison table, all 4 proteins, with/without outliers

Raw/unadjusted Mann-Whitney p-values throughout (matching what is drawn
on the figure panels) – not corrected for multiple comparisons.

``` r
run_two_group <- function(base_df, comp_key) {
    df2 <- comparison_builders[[comp_key]](base_df)
    event_label <- setdiff(levels(df2$clinical_group), "Control")
    a <- df2 %>% filter(clinical_group == "Control")  %>% pull(conc)
    b <- df2 %>% filter(clinical_group == event_label) %>% pull(conc)
    tibble(comparison = paste(event_label, "vs Control"),
           n_control = length(a), n_event = length(b),
           test = "Mann-Whitney", p_value = wilcox.test(a, b)$p.value)
}

run_component <- function(base_df) {
    df3 <- comparison_builders$Component(base_df)
    kw <- kruskal.test(conc ~ clinical_group, data = df3)
    n_by_group <- table(df3$clinical_group)
    omnibus <- tibble(comparison = "Cancer vs Cardiovascular vs Control (omnibus)",
                      n_control = n_by_group[["Control"]], n_event = NA_integer_,
                      test = "Kruskal-Wallis", p_value = kw$p.value)
    levs <- levels(droplevels(df3$clinical_group))
    pairwise <- bind_rows(lapply(combn(levs, 2, simplify = FALSE), function(pr) {
        a <- df3$conc[df3$clinical_group == pr[1]]
        b <- df3$conc[df3$clinical_group == pr[2]]
        tibble(comparison = paste(pr[2], "vs", pr[1]),
               n_control = length(a), n_event = length(b),
               test = "Mann-Whitney", p_value = wilcox.test(a, b)$p.value)
    }))
    bind_rows(omnibus, pairwise)
}

run_all <- function(base_df) {
    bind_rows(
        run_two_group(base_df, "Cases"),
        run_component(base_df),
        run_two_group(base_df, "Death")
    )
}

pairwise_tab <- bind_rows(lapply(panel_order, function(p) {
    bind_rows(
        run_all(scatter_dfs[[p]])       %>% mutate(protein = p, outliers = "included"),
        run_all(scatter_dfs_clean[[p]]) %>% mutate(protein = p, outliers = "excluded")
    )
})) %>% dplyr::select(protein, outliers, comparison, test, n_control, n_event, p_value)

pairwise_tab
```

    # A tibble: 48 × 7
       protein outliers comparison                   test  n_control n_event p_value
       <chr>   <chr>    <chr>                        <chr>     <int>   <int>   <dbl>
     1 SKAP2   included Cases vs Control             Mann…        69      67   0.938
     2 SKAP2   included Cancer vs Cardiovascular vs… Krus…        69      NA   0.932
     3 SKAP2   included Cancer vs Control            Mann…        69      46   0.925
     4 SKAP2   included Cardiovascular vs Control    Mann…        69      21   0.738
     5 SKAP2   included Cardiovascular vs Cancer     Mann…        46      21   0.743
     6 SKAP2   included Death (all-cause) vs Control Mann…        69      14   0.340
     7 SKAP2   excluded Cases vs Control             Mann…        67      60   0.358
     8 SKAP2   excluded Cancer vs Cardiovascular vs… Krus…        67      NA   0.450
     9 SKAP2   excluded Cancer vs Control            Mann…        67      42   0.647
    10 SKAP2   excluded Cardiovascular vs Control    Mann…        67      18   0.210
    # ℹ 38 more rows

``` r
pairwise_wide <- pairwise_tab %>%
    mutate(cell = sprintf("%s%s", ifelse(p_value < 0.05, "**", ""), format_p(p_value))) %>%  dplyr::select(protein, outliers, comparison, cell) %>%
    pivot_wider(names_from = comparison, values_from = cell) %>%
    arrange(protein, outliers)
pairwise_wide
```

    # A tibble: 8 × 8
      protein outliers `Cases vs Control` Cancer vs Cardiovasc…¹ `Cancer vs Control`
      <chr>   <chr>    <chr>              <chr>                  <chr>              
    1 CCL13   excluded 0.1555             0.3584                 0.2193             
    2 CCL13   included 0.1395             0.3281                 0.2035             
    3 DAPP1   excluded 0.3458             **0.0423               0.9162             
    4 DAPP1   included 0.8531             0.1769                 0.3793             
    5 LGALS9  excluded 0.1805             **0.0224               0.9260             
    6 LGALS9  included 0.0696             **0.0243               0.4915             
    7 SKAP2   excluded 0.3576             0.4505                 0.6472             
    8 SKAP2   included 0.9375             0.9323                 0.9250             
    # ℹ abbreviated name: ¹​`Cancer vs Cardiovascular vs Control (omnibus)`
    # ℹ 3 more variables: `Cardiovascular vs Control` <chr>,
    #   `Cardiovascular vs Cancer` <chr>, `Death (all-cause) vs Control` <chr>

``` r
n_sig <- sum(pairwise_tab$p_value < 0.05)
sprintf("Significant (p<0.05) tests: %d / %d", n_sig, nrow(pairwise_tab))
```

    [1] "Significant (p<0.05) tests: 10 / 48"

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
     abind          1.4-8   2024-09-12 [1] CRAN (R 4.6.0)
     backports      1.5.1   2026-04-03 [1] CRAN (R 4.6.0)
     bit            4.6.0   2025-03-06 [1] CRAN (R 4.6.0)
     bit64          4.8.2   2026-05-19 [1] CRAN (R 4.6.0)
     broom          1.0.13  2026-05-14 [1] CRAN (R 4.6.0)
     car            3.1-5   2026-02-03 [1] CRAN (R 4.6.0)
     carData        3.0-6   2026-01-30 [1] CRAN (R 4.6.0)
     cli            3.6.6   2026-04-09 [1] CRAN (R 4.6.0)
     crayon         1.5.3   2024-06-20 [1] CRAN (R 4.6.0)
     digest         0.6.39  2025-11-19 [1] CRAN (R 4.6.0)
     dplyr        * 1.2.1   2026-04-03 [1] CRAN (R 4.6.0)
     evaluate       1.0.5   2025-08-27 [1] CRAN (R 4.6.0)
     farver         2.1.2   2024-05-13 [1] CRAN (R 4.6.0)
     fastmap        1.2.0   2024-05-15 [1] CRAN (R 4.6.0)
     Formula        1.2-5   2023-02-24 [1] CRAN (R 4.6.0)
     generics       0.1.4   2025-05-09 [1] CRAN (R 4.6.0)
     ggplot2      * 4.0.3   2026-04-22 [1] CRAN (R 4.6.0)
     ggpubr       * 0.6.3   2026-02-24 [1] CRAN (R 4.6.0)
     ggsignif       0.6.4   2022-10-13 [1] CRAN (R 4.6.0)
     glue           1.8.1   2026-04-17 [1] CRAN (R 4.6.0)
     gtable         0.3.6   2024-10-25 [1] CRAN (R 4.6.0)
     hms            1.1.4   2025-10-17 [1] CRAN (R 4.6.0)
     htmltools      0.5.9   2025-12-04 [1] CRAN (R 4.6.0)
     jsonlite       2.0.0   2025-03-27 [1] CRAN (R 4.6.0)
     knitr          1.51    2025-12-20 [1] CRAN (R 4.6.0)
     labeling       0.4.3   2023-08-29 [1] CRAN (R 4.6.0)
     lifecycle      1.0.5   2026-01-08 [1] CRAN (R 4.6.0)
     magrittr       2.0.5   2026-04-04 [1] CRAN (R 4.6.0)
     otel           0.2.0   2025-08-29 [1] CRAN (R 4.6.0)
     patchwork    * 1.3.2   2025-08-25 [1] CRAN (R 4.6.0)
     pillar         1.11.1  2025-09-17 [1] CRAN (R 4.6.0)
     pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.6.0)
     purrr          1.2.2   2026-04-10 [1] CRAN (R 4.6.0)
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
     tibble         3.3.1   2026-01-11 [1] CRAN (R 4.6.0)
     tidyr        * 1.3.2   2025-12-19 [1] CRAN (R 4.6.0)
     tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.6.0)
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
