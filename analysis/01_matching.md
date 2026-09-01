# 1. Propensity-score matching


## Purpose

This script matches each case to one event-free control with an
available sample, using propensity scores based on baseline demographic,
clinical, immunovirological, and metabolic variables. This is the first
step in the pipeline and it decides who counts as a case and who counts
as a control before anything else happens.

## Reproducibility

This script needs the full, controlled-access CoRIS database and can’t
be run from the public repository. It’s included so the matching process
is auditable, not so it can be re-executed.

The `case_id`/`control_id` pairs produced here become `sid`/`sid_pair`
once pseudonymized in `02_1_preprocess_metadata.qmd`.

## Design summary

- Cases: participants with `evento_ever == 1`
- Eligible controls: participants with `evento_ever == 0` and
  `muestra_ok == 1`
- Propensity score estimated with logistic regression
- 1:1 nearest-neighbour matching
- Matching without replacement
- Caliper = 0.4
- Random ordering set before matching to avoid deterministic selection
  among equivalent candidates

## Variables included in the propensity-score model

| Variable     | Description                                        |
|--------------|----------------------------------------------------|
| `age`        | Age                                                |
| `sex`        | Sex                                                |
| `mode`       | HIV acquisition category                           |
| `CD4CD8_0`   | Baseline CD4/CD8 ratio                             |
| `cv100`      | Baseline HIV-1 RNA category / viral load indicator |
| `aids`       | Prior AIDS diagnosis                               |
| `edulvl`     | Educational level                                  |
| `paisorigen` | Country/region of origin                           |
| `art_reg`    | ART regimen                                        |
| `col_0`      | Baseline total cholesterol                         |
| `hdl_0`      | Baseline HDL cholesterol                           |

## Matching code (Stata)

It requires the private CoRIS database and Stata (`psmatch2`), which are
not part of this repository.

``` stata
clear all
set more off

* Load analysis dataset
use "Dataset", clear

* Create internal row identifier
gen id = _n

* Restrict the pool to:
*   1) all cases
*   2) controls with an available sample
keep if evento_ever == 1 | (evento_ever == 0 & muestra_ok == 1)

* Randomize row order before matching so that ties or equally suitable
* controls are not selected deterministically
set seed 12345
gen u = runiform()
sort u

* Estimate the propensity score
logistic evento_ever ///
age sex i.mode CD4CD8_0 cv100 aids ///
i.edulvl i.paisorigen i.art_reg ///
col_0 hdl_0

predict ps, pr

* Perform 1:1 nearest-neighbour propensity score matching
* - neighbor(1): one control per case
* - noreplacement: each control can be used only once
* - caliper(0.4): maximum allowed distance in propensity score
psmatch2 evento_ever, ///
pscore(ps) ///
neighbor(1) ///
noreplacement ///
caliper(0.4)

* Assess post-matching covariate balance
pstest ///
age sex i.mode CD4CD8_0 cv100 aids ///
i.edulvl i.paisorigen i.art_reg ///
col_0 hdl_0, both graph

* Create pair-level identifiers:
*   case_id    = id of treated participant (case)
*   control_id = id of matched untreated participant (control)
gen case_id    = id if _treated == 1
gen control_id = id[_n1] if _treated == 1 & !missing(_n1)

* Keep only successfully matched cases
keep if _treated == 1 & !missing(control_id)

* Inspect matched pairs
list case_id control_id ps if !missing(case_id), noobs sepby(case_id)

* Save matched pairs for downstream pseudonymization
save "results/matched_pairs.dta", replace
```

## Outputs

- `ps`: estimated propensity score
- Matching diagnostics from `pstest`
- `case_id` / `control_id`: identifiers for matched pairs, carried
  forward into `02_1_preprocess_metadata.qmd` as `sid` / `sid_pair`

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
     package     * version date (UTC) lib source
     cli           3.6.6   2026-04-09 [1] CRAN (R 4.6.0)
     digest        0.6.39  2025-11-19 [1] CRAN (R 4.6.0)
     evaluate      1.0.5   2025-08-27 [1] CRAN (R 4.6.0)
     fastmap       1.2.0   2024-05-15 [1] CRAN (R 4.6.0)
     htmltools     0.5.9   2025-12-04 [1] CRAN (R 4.6.0)
     jsonlite      2.0.0   2025-03-27 [1] CRAN (R 4.6.0)
     knitr         1.51    2025-12-20 [1] CRAN (R 4.6.0)
     otel          0.2.0   2025-08-29 [1] CRAN (R 4.6.0)
     rlang         1.2.0   2026-04-06 [1] CRAN (R 4.6.0)
     rmarkdown     2.31    2026-03-26 [1] CRAN (R 4.6.0)
     rstudioapi    0.18.0  2026-01-16 [1] CRAN (R 4.6.0)
     sessioninfo   1.2.3   2025-02-05 [1] CRAN (R 4.6.0)
     xfun          0.60    2026-07-09 [1] CRAN (R 4.6.0)
     yaml          2.3.12  2025-12-10 [1] CRAN (R 4.6.0)

     [1] /opt/Rlibs
     [2] /usr/local/lib/R/site-library
     [3] /usr/lib/R/site-library
     [4] /usr/lib/R/library

    ──────────────────────────────────────────────────────────────────────────────
