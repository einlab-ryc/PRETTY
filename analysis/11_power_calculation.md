# 11. Statistical power calculation


## Purpose

This script reproduces the study’s power calculation, the sample size
justification reported in the Methods, and the power curve in Fig. S12.

## Reproducibility

This script runs entirely from the public repository. No data files are
needed, only the fixed sample sizes and protein counts already stated in
the manuscript.

``` r
library(pwr)
```

``` r
d_effect <- 0.7
```

## Scenario 1: originally available matched population (103 pairs, N=206, 368-protein panel)

``` r
n1 <- 103
k1 <- 368
alpha1 <- 0.05 / k1

power1 <- pwr.t.test(n = n1, d = d_effect, sig.level = alpha1,
                      type = "two.sample", alternative = "two.sided")$power

mde1 <- pwr.t.test(n = n1, power = 0.8, sig.level = alpha1,
                    type = "two.sample", alternative = "two.sided")$d

tibble::tibble(
  n_per_group = n1, N = 2 * n1, proteins_tested = k1,
  alpha_bonferroni = alpha1,
  power_at_d0.7 = power1, power_manuscript = 0.87,
  mde_at_80pct_power = mde1, mde_manuscript = 0.66
)
```

    # A tibble: 1 × 8
      n_per_group     N proteins_tested alpha_bonferroni power_at_d0.7
            <dbl> <dbl>           <dbl>            <dbl>         <dbl>
    1         103   206             368         0.000136         0.868
    # ℹ 3 more variables: power_manuscript <dbl>, mde_at_80pct_power <dbl>,
    #   mde_manuscript <dbl>

## Scenario 2: final analytic set after QC (69 pairs, 138 participants, 326 proteins)

``` r
n2 <- 69
k2 <- 326
alpha2 <- 0.05 / k2

power2 <- pwr.t.test(n = n2, d = d_effect, sig.level = alpha2,
                      type = "two.sample", alternative = "two.sided")$power

mde2 <- pwr.t.test(n = n2, power = 0.8, sig.level = alpha2,
                    type = "two.sample", alternative = "two.sided")$d

tibble::tibble(
  n_per_group = n2, N = 2 * n2, proteins_passing_qc = k2,
  alpha_bonferroni = alpha2,
  power_at_d0.7 = power2, power_manuscript = 0.59,
  mde_at_80pct_power = mde2, mde_manuscript = 0.81
)
```

    # A tibble: 1 × 8
      n_per_group     N proteins_passing_qc alpha_bonferroni power_at_d0.7
            <dbl> <dbl>               <dbl>            <dbl>         <dbl>
    1          69   138                 326         0.000153         0.586
    # ℹ 3 more variables: power_manuscript <dbl>, mde_at_80pct_power <dbl>,
    #   mde_manuscript <dbl>

## Sanity check against the exact figures quoted in the Methods paragraph

``` r
checks <- data.frame(
  Quantity = c("Power, Scenario 1 (368 proteins)", "MDE @80% power, Scenario 1",
               "Power, Scenario 2 (326 proteins, post-QC)", "MDE @80% power, Scenario 2"),
  R        = round(c(power1, mde1, power2, mde2), 2),
  Methods  = c(0.87, 0.66, 0.59, 0.81)
)
checks$Match <- ifelse(abs(checks$R - checks$Methods) < 0.005, "OK", "MISMATCH")
checks
```

                                       Quantity    R Methods Match
    1          Power, Scenario 1 (368 proteins) 0.87    0.87    OK
    2                MDE @80% power, Scenario 1 0.66    0.66    OK
    3 Power, Scenario 2 (326 proteins, post-QC) 0.59    0.59    OK
    4                MDE @80% power, Scenario 2 0.81    0.81    OK

## Fig. S12: power vs. total sample size

d=0.7, alpha=0.05/326 (post-QC threshold, the one carried into the
manuscript figure), balanced groups (n1 = n2).

``` r
n_seq <- seq(40, 220, by = 20)
curve_power <- sapply(n_seq / 2, function(n_per_group) {
  pwr.t.test(n = n_per_group, d = d_effect, sig.level = alpha2,
             type = "two.sample", alternative = "two.sided")$power
})
curve_df <- data.frame(N_total = n_seq, N_per_group = n_seq / 2,
                        Power = round(curve_power, 4))
curve_df
```

       N_total N_per_group  Power
    1       40          20 0.0381
    2       60          30 0.1077
    3       80          40 0.2127
    4      100          50 0.3399
    5      120          60 0.4731
    6      140          70 0.5981
    7      160          80 0.7058
    8      180          90 0.7926
    9      200         100 0.8587
    10     220         110 0.9066

``` r
plot_power_curve <- function() {
  par(cex.axis = 1.3, cex.lab = 1.5, mar = c(6, 5.5, 5, 1.5))
  plot(curve_df$N_total, curve_df$Power, type = "b", pch = 19, cex = 1.3, lwd = 2, col = "#1a3a5c",
       xlab = "Total sample size (N, balanced n1 = n2)", ylab = "Power",
       sub = sprintf("Standardized between-group difference d = %.1f; alpha = 0.05/%d (Bonferroni)",
                      d_effect, k2),
       ylim = c(0, 1), cex.main = 1.7, cex.sub = 1.15)
  abline(h = 0.8, lty = 2, lwd = 1.5, col = "grey50")
  points(2 * n2, power2, pch = 17, col = "#c0392b", cex = 1.9)
  text(2 * n2, power2, labels = sprintf("  final set\n  N=%d, power=%.0f%%", 2 * n2, 100 * power2),
       pos = 4, cex = 1.15, col = "#c0392b")
  legend("bottomright", legend = "80% power threshold", lty = 2, lwd = 1.5, col = "grey50",
         bty = "n", cex = 1.15)
}

plot_power_curve()
```

![](11_power_calculation_files/figure-commonmark/unnamed-chunk-7-1.png)

``` r
# pdf("../results/FigS12_power_calculation.pdf", width = 10.5, height = 7)
plot_power_curve()
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
     package     * version date (UTC) lib source
     cli           3.6.6   2026-04-09 [1] CRAN (R 4.6.0)
     digest        0.6.39  2025-11-19 [1] CRAN (R 4.6.0)
     evaluate      1.0.5   2025-08-27 [1] CRAN (R 4.6.0)
     fastmap       1.2.0   2024-05-15 [1] CRAN (R 4.6.0)
     glue          1.8.1   2026-04-17 [1] CRAN (R 4.6.0)
     htmltools     0.5.9   2025-12-04 [1] CRAN (R 4.6.0)
     jsonlite      2.0.0   2025-03-27 [1] CRAN (R 4.6.0)
     knitr         1.51    2025-12-20 [1] CRAN (R 4.6.0)
     lifecycle     1.0.5   2026-01-08 [1] CRAN (R 4.6.0)
     magrittr      2.0.5   2026-04-04 [1] CRAN (R 4.6.0)
     otel          0.2.0   2025-08-29 [1] CRAN (R 4.6.0)
     pillar        1.11.1  2025-09-17 [1] CRAN (R 4.6.0)
     pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.6.0)
     pwr         * 1.3-0   2020-03-17 [1] CRAN (R 4.6.0)
     rlang         1.2.0   2026-04-06 [1] CRAN (R 4.6.0)
     rmarkdown     2.31    2026-03-26 [1] CRAN (R 4.6.0)
     rstudioapi    0.18.0  2026-01-16 [1] CRAN (R 4.6.0)
     sessioninfo   1.2.3   2025-02-05 [1] CRAN (R 4.6.0)
     tibble        3.3.1   2026-01-11 [1] CRAN (R 4.6.0)
     vctrs         0.7.3   2026-04-11 [1] CRAN (R 4.6.0)
     xfun          0.60    2026-07-09 [1] CRAN (R 4.6.0)
     yaml          2.3.12  2025-12-10 [1] CRAN (R 4.6.0)

     [1] /opt/Rlibs
     [2] /usr/local/lib/R/site-library
     [3] /usr/lib/R/site-library
     [4] /usr/lib/R/library
     * ── Packages attached to the search path.

    ──────────────────────────────────────────────────────────────────────────────
