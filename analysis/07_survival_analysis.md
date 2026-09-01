# 7. Survival analysis


## Purpose

This script fits Cox proportional hazards models for the four
ELISA-validated biomarkers across all outcomes. It produces the main
models results (Table 2), the outcome-specific tables (Tables S12 to
S14), the main Kaplan-Meier figure (Fig. 5) and the outcome-specific KM
figures (Figs. S7 to S9). It also produces the non-linearity spline
models (Fig. S6), the proportional hazards diagnostics (Fig. S11, Table
S16), and the quartile forest plot (Fig. S5).

## Reproducibility

This script runs entirely from the public repository. Only
`elisa_validation.tsv` and `metadata.tsv` are needed.

``` r
library(broom)
library(ggpubr)
library(ggsurvfit)
library(gtsummary)
library(patchwork)
library(survival)
library(tidyverse)
```

## Load public data

`time` (time-to-SNAE-or-censoring) is precomputed and included in
`metadata.tsv` directly, so this script is fully reproducible from the
public repository. See `01_matching.qmd` /
`02_1_preprocess_metadata.qmd` for how the underlying dates were used to
derive it (not runnable publicly, documented there for transparency).

``` r
metadata <- read_tsv("../data/metadata.tsv", show_col_types = FALSE) |>
    mutate(CD4_0_200 = factor(ifelse(CD4_0 < 200, 1, 0)))
```

# Models

The four ELISA-validated biomarkers (LGALS9, CCL13, DAPP1, SKAP2) go
through identical data prep and Cox modelling.

## Helper functions

``` r
markers <- c("LGALS9", "CCL13", "DAPP1", "SKAP2")
outcomes <- c("group", "nadm", "mace", "death_all")

# Builds one marker's analysis dataset: ELISA concentration + metadata +
# time-to-event, quartile-coded.
build_marker_data <- function(marker) {
    read_tsv("../data/elisa_validation.tsv") |>
        filter(assay == marker) |>
        left_join(
            metadata |> dplyr::select(sid, group, nadm, mace, death_all, CD4_0, CD4_0_200, time),
            by = "sid"
        ) |>
        mutate(quantile = factor(ntile(mean_conc, 4))) |>
        filter(time > 0) |>
        mutate(group = if_else(group == "Event", 1, 0))
}

# Categorical quartile Cox model (Q2/Q3/Q4 vs Q1), for one outcome, crude
# or CD4<200-adjusted.
fit_categorical_model <- function(data, event_var, adjusted = FALSE) {
    model_data <- data |>
        mutate(
            time, event = as.numeric(as.character(.data[[event_var]])),
            quantile, CD4_0_200,
            .keep = "none"
        ) |>
        filter(!is.na(time), time > 0, !is.na(event), !adjusted | !is.na(CD4_0_200))

    rhs <- if (adjusted) "quantile + CD4_0_200" else "quantile"
    coxph(as.formula(paste0("Surv(time, event) ~ ", rhs)), data = model_data)
}

# Linear-trend Cox model (one HR per one-quartile increase), for one
# outcome, crude or CD4<200-adjusted. Used for Global P values throughout
# (Table 2/S12-S14 and the KM figure annotations).
fit_trend_model <- function(data, event_var, adjusted = FALSE) {
    model_data <- data |>
        mutate(
            time, event = as.numeric(as.character(.data[[event_var]])),
            q = as.numeric(quantile), CD4_0_200,
            .keep = "none"
        ) |>
        filter(!is.na(time), time > 0, !is.na(event), !adjusted | !is.na(CD4_0_200))

    rhs <- if (adjusted) "q + CD4_0_200" else "q"
    coxph(as.formula(paste0("Surv(time, event) ~ ", rhs)), data = model_data)
}

# Event counts per quantile, for each outcome (used to populate Tables).
count_events_by_quantile <- function(data) {
    purrr::map(outcomes, function(outcome) {
        data |> group_by(quantile) |> dplyr::count(.data[[outcome]]) |> filter(.data[[outcome]] == 1)
    }) |> setNames(outcomes)
}
```

## Fit models

``` r
marker_data <- purrr::map(markers, build_marker_data) |> setNames(markers)
```

    Rows: 552 Columns: 3
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): sid, assay
    dbl (1): mean_conc

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    Rows: 552 Columns: 3
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): sid, assay
    dbl (1): mean_conc

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    Rows: 552 Columns: 3
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): sid, assay
    dbl (1): mean_conc

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    Rows: 552 Columns: 3
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr (2): sid, assay
    dbl (1): mean_conc

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Categorical (Q2/Q3/Q4 vs Q1) models, crude and CD4\<200-adjusted, for
every marker and outcome:

``` r
categorical_models <- purrr::map(marker_data, function(data) {
    purrr::map(outcomes, function(outcome) {
        list(
            crude = fit_categorical_model(data, outcome, adjusted = FALSE),
            adjusted = fit_categorical_model(data, outcome, adjusted = TRUE)
        )
    }) |> setNames(outcomes)
})
```

    Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
    Loglik converged before variable 1,2,3 ; coefficient may be infinite.
    Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
    Loglik converged before variable 1,2,3 ; coefficient may be infinite.
    Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
    Loglik converged before variable 1,2,3 ; coefficient may be infinite.
    Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
    Loglik converged before variable 1,2,3 ; coefficient may be infinite.

    Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
    Loglik converged before variable 2 ; coefficient may be infinite.
    Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
    Loglik converged before variable 2 ; coefficient may be infinite.

Linear-trend models, crude and CD4\<200-adjusted, for the Global P
values:

``` r
trend_models <- purrr::map(marker_data, function(data) {
    purrr::map(outcomes, function(outcome) {
        list(
            crude = fit_trend_model(data, outcome, adjusted = FALSE),
            adjusted = fit_trend_model(data, outcome, adjusted = TRUE)
        )
    }) |> setNames(outcomes)
})
```

Event counts per marker/outcome/quartile:

``` r
marker_event_counts <- purrr::map(marker_data, count_events_by_quantile)
```

# Tables

## Marker-level summary tables (unadjusted + adjusted, main SNAE comparison)

``` r
print_marker_tables <- function(models, marker) {
    print_one <- function(model) {
        tidy(model, conf.int = TRUE, exponentiate = TRUE) |> print()
        model |>
            tbl_regression(exponentiate = TRUE, label = list(quantile ~ paste(marker, "Quantile"))) |>
            bold_p(t = 0.05) |>
            bold_labels() |>
            print()
    }
    print_one(models$group$crude)
    print_one(models$group$adjusted)
}

purrr::iwalk(categorical_models, ~ print_marker_tables(.x, .y))
```

    # A tibble: 3 × 7
      term      estimate std.error statistic p.value conf.low conf.high
      <chr>        <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
    1 quantile2     2.14     0.375      2.03  0.0425    1.03       4.46
    2 quantile3     2.04     0.380      1.88  0.0604    0.969      4.30
    3 quantile4     2.37     0.368      2.35  0.0188    1.15       4.88
    <div id="nulphpthco" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
      <style>#nulphpthco table {
      font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
      -webkit-font-smoothing: antialiased;
      -moz-osx-font-smoothing: grayscale;
    }

    #nulphpthco thead, #nulphpthco tbody, #nulphpthco tfoot, #nulphpthco tr, #nulphpthco td, #nulphpthco th {
      border-style: none;
    }

    #nulphpthco p {
      margin: 0;
      padding: 0;
    }

    #nulphpthco .gt_table {
      display: table;
      border-collapse: collapse;
      line-height: normal;
      margin-left: auto;
      margin-right: auto;
      color: #333333;
      font-size: 16px;
      font-weight: normal;
      font-style: normal;
      background-color: #FFFFFF;
      width: auto;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #A8A8A8;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #A8A8A8;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
    }

    #nulphpthco .gt_caption {
      padding-top: 4px;
      padding-bottom: 4px;
    }

    #nulphpthco .gt_title {
      color: #333333;
      font-size: 125%;
      font-weight: initial;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-color: #FFFFFF;
      border-bottom-width: 0;
    }

    #nulphpthco .gt_subtitle {
      color: #333333;
      font-size: 85%;
      font-weight: initial;
      padding-top: 3px;
      padding-bottom: 5px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-color: #FFFFFF;
      border-top-width: 0;
    }

    #nulphpthco .gt_heading {
      background-color: #FFFFFF;
      text-align: center;
      border-bottom-color: #FFFFFF;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #nulphpthco .gt_bottom_border {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #nulphpthco .gt_col_headings {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #nulphpthco .gt_col_heading {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 6px;
      padding-left: 5px;
      padding-right: 5px;
      overflow-x: hidden;
    }

    #nulphpthco .gt_column_spanner_outer {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      padding-top: 0;
      padding-bottom: 0;
      padding-left: 4px;
      padding-right: 4px;
    }

    #nulphpthco .gt_column_spanner_outer:first-child {
      padding-left: 0;
    }

    #nulphpthco .gt_column_spanner_outer:last-child {
      padding-right: 0;
    }

    #nulphpthco .gt_column_spanner {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 5px;
      overflow-x: hidden;
      display: inline-block;
      width: 100%;
    }

    #nulphpthco .gt_spanner_row {
      border-bottom-style: hidden;
    }

    #nulphpthco .gt_group_heading {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      text-align: left;
    }

    #nulphpthco .gt_empty_group_heading {
      padding: 0.5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: middle;
    }

    #nulphpthco .gt_from_md > :first-child {
      margin-top: 0;
    }

    #nulphpthco .gt_from_md > :last-child {
      margin-bottom: 0;
    }

    #nulphpthco .gt_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      margin: 10px;
      border-top-style: solid;
      border-top-width: 1px;
      border-top-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      overflow-x: hidden;
    }

    #nulphpthco .gt_stub {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
    }

    #nulphpthco .gt_stub_row_group {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
      vertical-align: top;
    }

    #nulphpthco .gt_row_group_first td {
      border-top-width: 2px;
    }

    #nulphpthco .gt_row_group_first th {
      border-top-width: 2px;
    }

    #nulphpthco .gt_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #nulphpthco .gt_first_summary_row {
      border-top-style: solid;
      border-top-color: #D3D3D3;
    }

    #nulphpthco .gt_first_summary_row.thick {
      border-top-width: 2px;
    }

    #nulphpthco .gt_last_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #nulphpthco .gt_grand_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #nulphpthco .gt_first_grand_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-style: double;
      border-top-width: 6px;
      border-top-color: #D3D3D3;
    }

    #nulphpthco .gt_last_grand_summary_row_top {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: double;
      border-bottom-width: 6px;
      border-bottom-color: #D3D3D3;
    }

    #nulphpthco .gt_striped {
      background-color: rgba(128, 128, 128, 0.05);
    }

    #nulphpthco .gt_table_body {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #nulphpthco .gt_footnotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #nulphpthco .gt_footnote {
      margin: 0px;
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #nulphpthco .gt_sourcenotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #nulphpthco .gt_sourcenote {
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #nulphpthco .gt_left {
      text-align: left;
    }

    #nulphpthco .gt_center {
      text-align: center;
    }

    #nulphpthco .gt_right {
      text-align: right;
      font-variant-numeric: tabular-nums;
    }

    #nulphpthco .gt_font_normal {
      font-weight: normal;
    }

    #nulphpthco .gt_font_bold {
      font-weight: bold;
    }

    #nulphpthco .gt_font_italic {
      font-style: italic;
    }

    #nulphpthco .gt_super {
      font-size: 65%;
    }

    #nulphpthco .gt_footnote_marks {
      font-size: 75%;
      vertical-align: 0.4em;
      position: initial;
    }

    #nulphpthco .gt_asterisk {
      font-size: 100%;
      vertical-align: 0;
    }

    #nulphpthco .gt_indent_1 {
      text-indent: 5px;
    }

    #nulphpthco .gt_indent_2 {
      text-indent: 10px;
    }

    #nulphpthco .gt_indent_3 {
      text-indent: 15px;
    }

    #nulphpthco .gt_indent_4 {
      text-indent: 20px;
    }

    #nulphpthco .gt_indent_5 {
      text-indent: 25px;
    }

    #nulphpthco .katex-display {
      display: inline-flex !important;
      margin-bottom: 0.75em !important;
    }

    #nulphpthco div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
      height: 0px !important;
    }
    </style>
      <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
      <thead>
        <tr class="gt_col_headings">
          <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="label"><span data-qmd-base64="KipDaGFyYWN0ZXJpc3RpYyoq"><span class='gt_from_md'><strong>Characteristic</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipIUioq"><span class='gt_from_md'><strong>HR</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="conf.low"><span data-qmd-base64="Kio5NSUgQ0kqKg=="><span class='gt_from_md'><strong>95% CI</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="p.value"><span data-qmd-base64="KipwLXZhbHVlKio="><span class='gt_from_md'><strong>p-value</strong></span></span></th>
        </tr>
      </thead>
      <tbody class="gt_table_body">
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">LGALS9 Quantile</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    2</td>
    <td headers="estimate" class="gt_row gt_center">2.14</td>
    <td headers="conf.low" class="gt_row gt_center">1.03, 4.46</td>
    <td headers="p.value" class="gt_row gt_center" style="font-weight: bold;">0.043</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    3</td>
    <td headers="estimate" class="gt_row gt_center">2.04</td>
    <td headers="conf.low" class="gt_row gt_center">0.97, 4.30</td>
    <td headers="p.value" class="gt_row gt_center">0.060</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    4</td>
    <td headers="estimate" class="gt_row gt_center">2.37</td>
    <td headers="conf.low" class="gt_row gt_center">1.15, 4.88</td>
    <td headers="p.value" class="gt_row gt_center" style="font-weight: bold;">0.019</td></tr>
      </tbody>
      <tfoot>
        <tr class="gt_sourcenotes">
          <td class="gt_sourcenote" colspan="4"><span data-qmd-base64="QWJicmV2aWF0aW9uczogQ0kgPSBDb25maWRlbmNlIEludGVydmFsLCBIUiA9IEhhemFyZCBSYXRpbw=="><span class='gt_from_md'>Abbreviations: CI = Confidence Interval, HR = Hazard Ratio</span></span></td>
        </tr>
      </tfoot>
    </table>
    </div>
    # A tibble: 4 × 7
      term       estimate std.error statistic p.value conf.low conf.high
      <chr>         <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
    1 quantile2     2.27      0.398     2.06   0.0398    1.04       4.95
    2 quantile3     2.10      0.395     1.88   0.0595    0.971      4.56
    3 quantile4     2.46      0.403     2.23   0.0256    1.12       5.42
    4 CD4_0_2001    0.939     0.275    -0.228  0.820     0.547      1.61
    <div id="yitnyimznx" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
      <style>#yitnyimznx table {
      font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
      -webkit-font-smoothing: antialiased;
      -moz-osx-font-smoothing: grayscale;
    }

    #yitnyimznx thead, #yitnyimznx tbody, #yitnyimznx tfoot, #yitnyimznx tr, #yitnyimznx td, #yitnyimznx th {
      border-style: none;
    }

    #yitnyimznx p {
      margin: 0;
      padding: 0;
    }

    #yitnyimznx .gt_table {
      display: table;
      border-collapse: collapse;
      line-height: normal;
      margin-left: auto;
      margin-right: auto;
      color: #333333;
      font-size: 16px;
      font-weight: normal;
      font-style: normal;
      background-color: #FFFFFF;
      width: auto;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #A8A8A8;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #A8A8A8;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
    }

    #yitnyimznx .gt_caption {
      padding-top: 4px;
      padding-bottom: 4px;
    }

    #yitnyimznx .gt_title {
      color: #333333;
      font-size: 125%;
      font-weight: initial;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-color: #FFFFFF;
      border-bottom-width: 0;
    }

    #yitnyimznx .gt_subtitle {
      color: #333333;
      font-size: 85%;
      font-weight: initial;
      padding-top: 3px;
      padding-bottom: 5px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-color: #FFFFFF;
      border-top-width: 0;
    }

    #yitnyimznx .gt_heading {
      background-color: #FFFFFF;
      text-align: center;
      border-bottom-color: #FFFFFF;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #yitnyimznx .gt_bottom_border {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #yitnyimznx .gt_col_headings {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #yitnyimznx .gt_col_heading {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 6px;
      padding-left: 5px;
      padding-right: 5px;
      overflow-x: hidden;
    }

    #yitnyimznx .gt_column_spanner_outer {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      padding-top: 0;
      padding-bottom: 0;
      padding-left: 4px;
      padding-right: 4px;
    }

    #yitnyimznx .gt_column_spanner_outer:first-child {
      padding-left: 0;
    }

    #yitnyimznx .gt_column_spanner_outer:last-child {
      padding-right: 0;
    }

    #yitnyimznx .gt_column_spanner {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 5px;
      overflow-x: hidden;
      display: inline-block;
      width: 100%;
    }

    #yitnyimznx .gt_spanner_row {
      border-bottom-style: hidden;
    }

    #yitnyimznx .gt_group_heading {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      text-align: left;
    }

    #yitnyimznx .gt_empty_group_heading {
      padding: 0.5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: middle;
    }

    #yitnyimznx .gt_from_md > :first-child {
      margin-top: 0;
    }

    #yitnyimznx .gt_from_md > :last-child {
      margin-bottom: 0;
    }

    #yitnyimznx .gt_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      margin: 10px;
      border-top-style: solid;
      border-top-width: 1px;
      border-top-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      overflow-x: hidden;
    }

    #yitnyimznx .gt_stub {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
    }

    #yitnyimznx .gt_stub_row_group {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
      vertical-align: top;
    }

    #yitnyimznx .gt_row_group_first td {
      border-top-width: 2px;
    }

    #yitnyimznx .gt_row_group_first th {
      border-top-width: 2px;
    }

    #yitnyimznx .gt_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #yitnyimznx .gt_first_summary_row {
      border-top-style: solid;
      border-top-color: #D3D3D3;
    }

    #yitnyimznx .gt_first_summary_row.thick {
      border-top-width: 2px;
    }

    #yitnyimznx .gt_last_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #yitnyimznx .gt_grand_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #yitnyimznx .gt_first_grand_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-style: double;
      border-top-width: 6px;
      border-top-color: #D3D3D3;
    }

    #yitnyimznx .gt_last_grand_summary_row_top {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: double;
      border-bottom-width: 6px;
      border-bottom-color: #D3D3D3;
    }

    #yitnyimznx .gt_striped {
      background-color: rgba(128, 128, 128, 0.05);
    }

    #yitnyimznx .gt_table_body {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #yitnyimznx .gt_footnotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #yitnyimznx .gt_footnote {
      margin: 0px;
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #yitnyimznx .gt_sourcenotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #yitnyimznx .gt_sourcenote {
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #yitnyimznx .gt_left {
      text-align: left;
    }

    #yitnyimznx .gt_center {
      text-align: center;
    }

    #yitnyimznx .gt_right {
      text-align: right;
      font-variant-numeric: tabular-nums;
    }

    #yitnyimznx .gt_font_normal {
      font-weight: normal;
    }

    #yitnyimznx .gt_font_bold {
      font-weight: bold;
    }

    #yitnyimznx .gt_font_italic {
      font-style: italic;
    }

    #yitnyimznx .gt_super {
      font-size: 65%;
    }

    #yitnyimznx .gt_footnote_marks {
      font-size: 75%;
      vertical-align: 0.4em;
      position: initial;
    }

    #yitnyimznx .gt_asterisk {
      font-size: 100%;
      vertical-align: 0;
    }

    #yitnyimznx .gt_indent_1 {
      text-indent: 5px;
    }

    #yitnyimznx .gt_indent_2 {
      text-indent: 10px;
    }

    #yitnyimznx .gt_indent_3 {
      text-indent: 15px;
    }

    #yitnyimznx .gt_indent_4 {
      text-indent: 20px;
    }

    #yitnyimznx .gt_indent_5 {
      text-indent: 25px;
    }

    #yitnyimznx .katex-display {
      display: inline-flex !important;
      margin-bottom: 0.75em !important;
    }

    #yitnyimznx div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
      height: 0px !important;
    }
    </style>
      <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
      <thead>
        <tr class="gt_col_headings">
          <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="label"><span data-qmd-base64="KipDaGFyYWN0ZXJpc3RpYyoq"><span class='gt_from_md'><strong>Characteristic</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipIUioq"><span class='gt_from_md'><strong>HR</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="conf.low"><span data-qmd-base64="Kio5NSUgQ0kqKg=="><span class='gt_from_md'><strong>95% CI</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="p.value"><span data-qmd-base64="KipwLXZhbHVlKio="><span class='gt_from_md'><strong>p-value</strong></span></span></th>
        </tr>
      </thead>
      <tbody class="gt_table_body">
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">LGALS9 Quantile</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    2</td>
    <td headers="estimate" class="gt_row gt_center">2.27</td>
    <td headers="conf.low" class="gt_row gt_center">1.04, 4.95</td>
    <td headers="p.value" class="gt_row gt_center" style="font-weight: bold;">0.040</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    3</td>
    <td headers="estimate" class="gt_row gt_center">2.10</td>
    <td headers="conf.low" class="gt_row gt_center">0.97, 4.56</td>
    <td headers="p.value" class="gt_row gt_center">0.059</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    4</td>
    <td headers="estimate" class="gt_row gt_center">2.46</td>
    <td headers="conf.low" class="gt_row gt_center">1.12, 5.42</td>
    <td headers="p.value" class="gt_row gt_center" style="font-weight: bold;">0.026</td></tr>
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">CD4_0_200</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    0</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">0.94</td>
    <td headers="conf.low" class="gt_row gt_center">0.55, 1.61</td>
    <td headers="p.value" class="gt_row gt_center">0.8</td></tr>
      </tbody>
      <tfoot>
        <tr class="gt_sourcenotes">
          <td class="gt_sourcenote" colspan="4"><span data-qmd-base64="QWJicmV2aWF0aW9uczogQ0kgPSBDb25maWRlbmNlIEludGVydmFsLCBIUiA9IEhhemFyZCBSYXRpbw=="><span class='gt_from_md'>Abbreviations: CI = Confidence Interval, HR = Hazard Ratio</span></span></td>
        </tr>
      </tfoot>
    </table>
    </div>
    # A tibble: 3 × 7
      term      estimate std.error statistic p.value conf.low conf.high
      <chr>        <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
    1 quantile2     1.45     0.378     0.990  0.322     0.693      3.05
    2 quantile3     2.35     0.359     2.38   0.0175    1.16       4.75
    3 quantile4     1.30     0.383     0.695  0.487     0.616      2.76
    <div id="mfjkoarhsp" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
      <style>#mfjkoarhsp table {
      font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
      -webkit-font-smoothing: antialiased;
      -moz-osx-font-smoothing: grayscale;
    }

    #mfjkoarhsp thead, #mfjkoarhsp tbody, #mfjkoarhsp tfoot, #mfjkoarhsp tr, #mfjkoarhsp td, #mfjkoarhsp th {
      border-style: none;
    }

    #mfjkoarhsp p {
      margin: 0;
      padding: 0;
    }

    #mfjkoarhsp .gt_table {
      display: table;
      border-collapse: collapse;
      line-height: normal;
      margin-left: auto;
      margin-right: auto;
      color: #333333;
      font-size: 16px;
      font-weight: normal;
      font-style: normal;
      background-color: #FFFFFF;
      width: auto;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #A8A8A8;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #A8A8A8;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
    }

    #mfjkoarhsp .gt_caption {
      padding-top: 4px;
      padding-bottom: 4px;
    }

    #mfjkoarhsp .gt_title {
      color: #333333;
      font-size: 125%;
      font-weight: initial;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-color: #FFFFFF;
      border-bottom-width: 0;
    }

    #mfjkoarhsp .gt_subtitle {
      color: #333333;
      font-size: 85%;
      font-weight: initial;
      padding-top: 3px;
      padding-bottom: 5px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-color: #FFFFFF;
      border-top-width: 0;
    }

    #mfjkoarhsp .gt_heading {
      background-color: #FFFFFF;
      text-align: center;
      border-bottom-color: #FFFFFF;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #mfjkoarhsp .gt_bottom_border {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #mfjkoarhsp .gt_col_headings {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #mfjkoarhsp .gt_col_heading {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 6px;
      padding-left: 5px;
      padding-right: 5px;
      overflow-x: hidden;
    }

    #mfjkoarhsp .gt_column_spanner_outer {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      padding-top: 0;
      padding-bottom: 0;
      padding-left: 4px;
      padding-right: 4px;
    }

    #mfjkoarhsp .gt_column_spanner_outer:first-child {
      padding-left: 0;
    }

    #mfjkoarhsp .gt_column_spanner_outer:last-child {
      padding-right: 0;
    }

    #mfjkoarhsp .gt_column_spanner {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 5px;
      overflow-x: hidden;
      display: inline-block;
      width: 100%;
    }

    #mfjkoarhsp .gt_spanner_row {
      border-bottom-style: hidden;
    }

    #mfjkoarhsp .gt_group_heading {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      text-align: left;
    }

    #mfjkoarhsp .gt_empty_group_heading {
      padding: 0.5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: middle;
    }

    #mfjkoarhsp .gt_from_md > :first-child {
      margin-top: 0;
    }

    #mfjkoarhsp .gt_from_md > :last-child {
      margin-bottom: 0;
    }

    #mfjkoarhsp .gt_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      margin: 10px;
      border-top-style: solid;
      border-top-width: 1px;
      border-top-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      overflow-x: hidden;
    }

    #mfjkoarhsp .gt_stub {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
    }

    #mfjkoarhsp .gt_stub_row_group {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
      vertical-align: top;
    }

    #mfjkoarhsp .gt_row_group_first td {
      border-top-width: 2px;
    }

    #mfjkoarhsp .gt_row_group_first th {
      border-top-width: 2px;
    }

    #mfjkoarhsp .gt_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #mfjkoarhsp .gt_first_summary_row {
      border-top-style: solid;
      border-top-color: #D3D3D3;
    }

    #mfjkoarhsp .gt_first_summary_row.thick {
      border-top-width: 2px;
    }

    #mfjkoarhsp .gt_last_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #mfjkoarhsp .gt_grand_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #mfjkoarhsp .gt_first_grand_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-style: double;
      border-top-width: 6px;
      border-top-color: #D3D3D3;
    }

    #mfjkoarhsp .gt_last_grand_summary_row_top {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: double;
      border-bottom-width: 6px;
      border-bottom-color: #D3D3D3;
    }

    #mfjkoarhsp .gt_striped {
      background-color: rgba(128, 128, 128, 0.05);
    }

    #mfjkoarhsp .gt_table_body {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #mfjkoarhsp .gt_footnotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #mfjkoarhsp .gt_footnote {
      margin: 0px;
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #mfjkoarhsp .gt_sourcenotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #mfjkoarhsp .gt_sourcenote {
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #mfjkoarhsp .gt_left {
      text-align: left;
    }

    #mfjkoarhsp .gt_center {
      text-align: center;
    }

    #mfjkoarhsp .gt_right {
      text-align: right;
      font-variant-numeric: tabular-nums;
    }

    #mfjkoarhsp .gt_font_normal {
      font-weight: normal;
    }

    #mfjkoarhsp .gt_font_bold {
      font-weight: bold;
    }

    #mfjkoarhsp .gt_font_italic {
      font-style: italic;
    }

    #mfjkoarhsp .gt_super {
      font-size: 65%;
    }

    #mfjkoarhsp .gt_footnote_marks {
      font-size: 75%;
      vertical-align: 0.4em;
      position: initial;
    }

    #mfjkoarhsp .gt_asterisk {
      font-size: 100%;
      vertical-align: 0;
    }

    #mfjkoarhsp .gt_indent_1 {
      text-indent: 5px;
    }

    #mfjkoarhsp .gt_indent_2 {
      text-indent: 10px;
    }

    #mfjkoarhsp .gt_indent_3 {
      text-indent: 15px;
    }

    #mfjkoarhsp .gt_indent_4 {
      text-indent: 20px;
    }

    #mfjkoarhsp .gt_indent_5 {
      text-indent: 25px;
    }

    #mfjkoarhsp .katex-display {
      display: inline-flex !important;
      margin-bottom: 0.75em !important;
    }

    #mfjkoarhsp div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
      height: 0px !important;
    }
    </style>
      <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
      <thead>
        <tr class="gt_col_headings">
          <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="label"><span data-qmd-base64="KipDaGFyYWN0ZXJpc3RpYyoq"><span class='gt_from_md'><strong>Characteristic</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipIUioq"><span class='gt_from_md'><strong>HR</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="conf.low"><span data-qmd-base64="Kio5NSUgQ0kqKg=="><span class='gt_from_md'><strong>95% CI</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="p.value"><span data-qmd-base64="KipwLXZhbHVlKio="><span class='gt_from_md'><strong>p-value</strong></span></span></th>
        </tr>
      </thead>
      <tbody class="gt_table_body">
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">CCL13 Quantile</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    2</td>
    <td headers="estimate" class="gt_row gt_center">1.45</td>
    <td headers="conf.low" class="gt_row gt_center">0.69, 3.05</td>
    <td headers="p.value" class="gt_row gt_center">0.3</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    3</td>
    <td headers="estimate" class="gt_row gt_center">2.35</td>
    <td headers="conf.low" class="gt_row gt_center">1.16, 4.75</td>
    <td headers="p.value" class="gt_row gt_center" style="font-weight: bold;">0.018</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    4</td>
    <td headers="estimate" class="gt_row gt_center">1.30</td>
    <td headers="conf.low" class="gt_row gt_center">0.62, 2.76</td>
    <td headers="p.value" class="gt_row gt_center">0.5</td></tr>
      </tbody>
      <tfoot>
        <tr class="gt_sourcenotes">
          <td class="gt_sourcenote" colspan="4"><span data-qmd-base64="QWJicmV2aWF0aW9uczogQ0kgPSBDb25maWRlbmNlIEludGVydmFsLCBIUiA9IEhhemFyZCBSYXRpbw=="><span class='gt_from_md'>Abbreviations: CI = Confidence Interval, HR = Hazard Ratio</span></span></td>
        </tr>
      </tfoot>
    </table>
    </div>
    # A tibble: 4 × 7
      term       estimate std.error statistic p.value conf.low conf.high
      <chr>         <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
    1 quantile2      1.36     0.406     0.749  0.454     0.612      3.00
    2 quantile3      2.39     0.373     2.34   0.0193    1.15       4.97
    3 quantile4      1.34     0.393     0.744  0.457     0.620      2.89
    4 CD4_0_2001     1.14     0.256     0.524  0.600     0.693      1.89
    <div id="boiogcsbzw" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
      <style>#boiogcsbzw table {
      font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
      -webkit-font-smoothing: antialiased;
      -moz-osx-font-smoothing: grayscale;
    }

    #boiogcsbzw thead, #boiogcsbzw tbody, #boiogcsbzw tfoot, #boiogcsbzw tr, #boiogcsbzw td, #boiogcsbzw th {
      border-style: none;
    }

    #boiogcsbzw p {
      margin: 0;
      padding: 0;
    }

    #boiogcsbzw .gt_table {
      display: table;
      border-collapse: collapse;
      line-height: normal;
      margin-left: auto;
      margin-right: auto;
      color: #333333;
      font-size: 16px;
      font-weight: normal;
      font-style: normal;
      background-color: #FFFFFF;
      width: auto;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #A8A8A8;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #A8A8A8;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
    }

    #boiogcsbzw .gt_caption {
      padding-top: 4px;
      padding-bottom: 4px;
    }

    #boiogcsbzw .gt_title {
      color: #333333;
      font-size: 125%;
      font-weight: initial;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-color: #FFFFFF;
      border-bottom-width: 0;
    }

    #boiogcsbzw .gt_subtitle {
      color: #333333;
      font-size: 85%;
      font-weight: initial;
      padding-top: 3px;
      padding-bottom: 5px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-color: #FFFFFF;
      border-top-width: 0;
    }

    #boiogcsbzw .gt_heading {
      background-color: #FFFFFF;
      text-align: center;
      border-bottom-color: #FFFFFF;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #boiogcsbzw .gt_bottom_border {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #boiogcsbzw .gt_col_headings {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #boiogcsbzw .gt_col_heading {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 6px;
      padding-left: 5px;
      padding-right: 5px;
      overflow-x: hidden;
    }

    #boiogcsbzw .gt_column_spanner_outer {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      padding-top: 0;
      padding-bottom: 0;
      padding-left: 4px;
      padding-right: 4px;
    }

    #boiogcsbzw .gt_column_spanner_outer:first-child {
      padding-left: 0;
    }

    #boiogcsbzw .gt_column_spanner_outer:last-child {
      padding-right: 0;
    }

    #boiogcsbzw .gt_column_spanner {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 5px;
      overflow-x: hidden;
      display: inline-block;
      width: 100%;
    }

    #boiogcsbzw .gt_spanner_row {
      border-bottom-style: hidden;
    }

    #boiogcsbzw .gt_group_heading {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      text-align: left;
    }

    #boiogcsbzw .gt_empty_group_heading {
      padding: 0.5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: middle;
    }

    #boiogcsbzw .gt_from_md > :first-child {
      margin-top: 0;
    }

    #boiogcsbzw .gt_from_md > :last-child {
      margin-bottom: 0;
    }

    #boiogcsbzw .gt_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      margin: 10px;
      border-top-style: solid;
      border-top-width: 1px;
      border-top-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      overflow-x: hidden;
    }

    #boiogcsbzw .gt_stub {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
    }

    #boiogcsbzw .gt_stub_row_group {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
      vertical-align: top;
    }

    #boiogcsbzw .gt_row_group_first td {
      border-top-width: 2px;
    }

    #boiogcsbzw .gt_row_group_first th {
      border-top-width: 2px;
    }

    #boiogcsbzw .gt_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #boiogcsbzw .gt_first_summary_row {
      border-top-style: solid;
      border-top-color: #D3D3D3;
    }

    #boiogcsbzw .gt_first_summary_row.thick {
      border-top-width: 2px;
    }

    #boiogcsbzw .gt_last_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #boiogcsbzw .gt_grand_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #boiogcsbzw .gt_first_grand_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-style: double;
      border-top-width: 6px;
      border-top-color: #D3D3D3;
    }

    #boiogcsbzw .gt_last_grand_summary_row_top {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: double;
      border-bottom-width: 6px;
      border-bottom-color: #D3D3D3;
    }

    #boiogcsbzw .gt_striped {
      background-color: rgba(128, 128, 128, 0.05);
    }

    #boiogcsbzw .gt_table_body {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #boiogcsbzw .gt_footnotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #boiogcsbzw .gt_footnote {
      margin: 0px;
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #boiogcsbzw .gt_sourcenotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #boiogcsbzw .gt_sourcenote {
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #boiogcsbzw .gt_left {
      text-align: left;
    }

    #boiogcsbzw .gt_center {
      text-align: center;
    }

    #boiogcsbzw .gt_right {
      text-align: right;
      font-variant-numeric: tabular-nums;
    }

    #boiogcsbzw .gt_font_normal {
      font-weight: normal;
    }

    #boiogcsbzw .gt_font_bold {
      font-weight: bold;
    }

    #boiogcsbzw .gt_font_italic {
      font-style: italic;
    }

    #boiogcsbzw .gt_super {
      font-size: 65%;
    }

    #boiogcsbzw .gt_footnote_marks {
      font-size: 75%;
      vertical-align: 0.4em;
      position: initial;
    }

    #boiogcsbzw .gt_asterisk {
      font-size: 100%;
      vertical-align: 0;
    }

    #boiogcsbzw .gt_indent_1 {
      text-indent: 5px;
    }

    #boiogcsbzw .gt_indent_2 {
      text-indent: 10px;
    }

    #boiogcsbzw .gt_indent_3 {
      text-indent: 15px;
    }

    #boiogcsbzw .gt_indent_4 {
      text-indent: 20px;
    }

    #boiogcsbzw .gt_indent_5 {
      text-indent: 25px;
    }

    #boiogcsbzw .katex-display {
      display: inline-flex !important;
      margin-bottom: 0.75em !important;
    }

    #boiogcsbzw div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
      height: 0px !important;
    }
    </style>
      <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
      <thead>
        <tr class="gt_col_headings">
          <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="label"><span data-qmd-base64="KipDaGFyYWN0ZXJpc3RpYyoq"><span class='gt_from_md'><strong>Characteristic</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipIUioq"><span class='gt_from_md'><strong>HR</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="conf.low"><span data-qmd-base64="Kio5NSUgQ0kqKg=="><span class='gt_from_md'><strong>95% CI</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="p.value"><span data-qmd-base64="KipwLXZhbHVlKio="><span class='gt_from_md'><strong>p-value</strong></span></span></th>
        </tr>
      </thead>
      <tbody class="gt_table_body">
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">CCL13 Quantile</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    2</td>
    <td headers="estimate" class="gt_row gt_center">1.36</td>
    <td headers="conf.low" class="gt_row gt_center">0.61, 3.00</td>
    <td headers="p.value" class="gt_row gt_center">0.5</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    3</td>
    <td headers="estimate" class="gt_row gt_center">2.39</td>
    <td headers="conf.low" class="gt_row gt_center">1.15, 4.97</td>
    <td headers="p.value" class="gt_row gt_center" style="font-weight: bold;">0.019</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    4</td>
    <td headers="estimate" class="gt_row gt_center">1.34</td>
    <td headers="conf.low" class="gt_row gt_center">0.62, 2.89</td>
    <td headers="p.value" class="gt_row gt_center">0.5</td></tr>
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">CD4_0_200</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    0</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">1.14</td>
    <td headers="conf.low" class="gt_row gt_center">0.69, 1.89</td>
    <td headers="p.value" class="gt_row gt_center">0.6</td></tr>
      </tbody>
      <tfoot>
        <tr class="gt_sourcenotes">
          <td class="gt_sourcenote" colspan="4"><span data-qmd-base64="QWJicmV2aWF0aW9uczogQ0kgPSBDb25maWRlbmNlIEludGVydmFsLCBIUiA9IEhhemFyZCBSYXRpbw=="><span class='gt_from_md'>Abbreviations: CI = Confidence Interval, HR = Hazard Ratio</span></span></td>
        </tr>
      </tfoot>
    </table>
    </div>
    # A tibble: 3 × 7
      term      estimate std.error statistic p.value conf.low conf.high
      <chr>        <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
    1 quantile2     1.08     0.361     0.202 0.840      0.530      2.18
    2 quantile3     2.99     0.352     3.11  0.00187    1.50       5.97
    3 quantile4     1.05     0.378     0.123 0.902      0.499      2.20
    <div id="hmgwnvuuua" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
      <style>#hmgwnvuuua table {
      font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
      -webkit-font-smoothing: antialiased;
      -moz-osx-font-smoothing: grayscale;
    }

    #hmgwnvuuua thead, #hmgwnvuuua tbody, #hmgwnvuuua tfoot, #hmgwnvuuua tr, #hmgwnvuuua td, #hmgwnvuuua th {
      border-style: none;
    }

    #hmgwnvuuua p {
      margin: 0;
      padding: 0;
    }

    #hmgwnvuuua .gt_table {
      display: table;
      border-collapse: collapse;
      line-height: normal;
      margin-left: auto;
      margin-right: auto;
      color: #333333;
      font-size: 16px;
      font-weight: normal;
      font-style: normal;
      background-color: #FFFFFF;
      width: auto;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #A8A8A8;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #A8A8A8;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
    }

    #hmgwnvuuua .gt_caption {
      padding-top: 4px;
      padding-bottom: 4px;
    }

    #hmgwnvuuua .gt_title {
      color: #333333;
      font-size: 125%;
      font-weight: initial;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-color: #FFFFFF;
      border-bottom-width: 0;
    }

    #hmgwnvuuua .gt_subtitle {
      color: #333333;
      font-size: 85%;
      font-weight: initial;
      padding-top: 3px;
      padding-bottom: 5px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-color: #FFFFFF;
      border-top-width: 0;
    }

    #hmgwnvuuua .gt_heading {
      background-color: #FFFFFF;
      text-align: center;
      border-bottom-color: #FFFFFF;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #hmgwnvuuua .gt_bottom_border {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #hmgwnvuuua .gt_col_headings {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #hmgwnvuuua .gt_col_heading {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 6px;
      padding-left: 5px;
      padding-right: 5px;
      overflow-x: hidden;
    }

    #hmgwnvuuua .gt_column_spanner_outer {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      padding-top: 0;
      padding-bottom: 0;
      padding-left: 4px;
      padding-right: 4px;
    }

    #hmgwnvuuua .gt_column_spanner_outer:first-child {
      padding-left: 0;
    }

    #hmgwnvuuua .gt_column_spanner_outer:last-child {
      padding-right: 0;
    }

    #hmgwnvuuua .gt_column_spanner {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 5px;
      overflow-x: hidden;
      display: inline-block;
      width: 100%;
    }

    #hmgwnvuuua .gt_spanner_row {
      border-bottom-style: hidden;
    }

    #hmgwnvuuua .gt_group_heading {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      text-align: left;
    }

    #hmgwnvuuua .gt_empty_group_heading {
      padding: 0.5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: middle;
    }

    #hmgwnvuuua .gt_from_md > :first-child {
      margin-top: 0;
    }

    #hmgwnvuuua .gt_from_md > :last-child {
      margin-bottom: 0;
    }

    #hmgwnvuuua .gt_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      margin: 10px;
      border-top-style: solid;
      border-top-width: 1px;
      border-top-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      overflow-x: hidden;
    }

    #hmgwnvuuua .gt_stub {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
    }

    #hmgwnvuuua .gt_stub_row_group {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
      vertical-align: top;
    }

    #hmgwnvuuua .gt_row_group_first td {
      border-top-width: 2px;
    }

    #hmgwnvuuua .gt_row_group_first th {
      border-top-width: 2px;
    }

    #hmgwnvuuua .gt_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #hmgwnvuuua .gt_first_summary_row {
      border-top-style: solid;
      border-top-color: #D3D3D3;
    }

    #hmgwnvuuua .gt_first_summary_row.thick {
      border-top-width: 2px;
    }

    #hmgwnvuuua .gt_last_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #hmgwnvuuua .gt_grand_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #hmgwnvuuua .gt_first_grand_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-style: double;
      border-top-width: 6px;
      border-top-color: #D3D3D3;
    }

    #hmgwnvuuua .gt_last_grand_summary_row_top {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: double;
      border-bottom-width: 6px;
      border-bottom-color: #D3D3D3;
    }

    #hmgwnvuuua .gt_striped {
      background-color: rgba(128, 128, 128, 0.05);
    }

    #hmgwnvuuua .gt_table_body {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #hmgwnvuuua .gt_footnotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #hmgwnvuuua .gt_footnote {
      margin: 0px;
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #hmgwnvuuua .gt_sourcenotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #hmgwnvuuua .gt_sourcenote {
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #hmgwnvuuua .gt_left {
      text-align: left;
    }

    #hmgwnvuuua .gt_center {
      text-align: center;
    }

    #hmgwnvuuua .gt_right {
      text-align: right;
      font-variant-numeric: tabular-nums;
    }

    #hmgwnvuuua .gt_font_normal {
      font-weight: normal;
    }

    #hmgwnvuuua .gt_font_bold {
      font-weight: bold;
    }

    #hmgwnvuuua .gt_font_italic {
      font-style: italic;
    }

    #hmgwnvuuua .gt_super {
      font-size: 65%;
    }

    #hmgwnvuuua .gt_footnote_marks {
      font-size: 75%;
      vertical-align: 0.4em;
      position: initial;
    }

    #hmgwnvuuua .gt_asterisk {
      font-size: 100%;
      vertical-align: 0;
    }

    #hmgwnvuuua .gt_indent_1 {
      text-indent: 5px;
    }

    #hmgwnvuuua .gt_indent_2 {
      text-indent: 10px;
    }

    #hmgwnvuuua .gt_indent_3 {
      text-indent: 15px;
    }

    #hmgwnvuuua .gt_indent_4 {
      text-indent: 20px;
    }

    #hmgwnvuuua .gt_indent_5 {
      text-indent: 25px;
    }

    #hmgwnvuuua .katex-display {
      display: inline-flex !important;
      margin-bottom: 0.75em !important;
    }

    #hmgwnvuuua div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
      height: 0px !important;
    }
    </style>
      <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
      <thead>
        <tr class="gt_col_headings">
          <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="label"><span data-qmd-base64="KipDaGFyYWN0ZXJpc3RpYyoq"><span class='gt_from_md'><strong>Characteristic</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipIUioq"><span class='gt_from_md'><strong>HR</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="conf.low"><span data-qmd-base64="Kio5NSUgQ0kqKg=="><span class='gt_from_md'><strong>95% CI</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="p.value"><span data-qmd-base64="KipwLXZhbHVlKio="><span class='gt_from_md'><strong>p-value</strong></span></span></th>
        </tr>
      </thead>
      <tbody class="gt_table_body">
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">DAPP1 Quantile</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    2</td>
    <td headers="estimate" class="gt_row gt_center">1.08</td>
    <td headers="conf.low" class="gt_row gt_center">0.53, 2.18</td>
    <td headers="p.value" class="gt_row gt_center">0.8</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    3</td>
    <td headers="estimate" class="gt_row gt_center">2.99</td>
    <td headers="conf.low" class="gt_row gt_center">1.50, 5.97</td>
    <td headers="p.value" class="gt_row gt_center" style="font-weight: bold;">0.002</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    4</td>
    <td headers="estimate" class="gt_row gt_center">1.05</td>
    <td headers="conf.low" class="gt_row gt_center">0.50, 2.20</td>
    <td headers="p.value" class="gt_row gt_center">>0.9</td></tr>
      </tbody>
      <tfoot>
        <tr class="gt_sourcenotes">
          <td class="gt_sourcenote" colspan="4"><span data-qmd-base64="QWJicmV2aWF0aW9uczogQ0kgPSBDb25maWRlbmNlIEludGVydmFsLCBIUiA9IEhhemFyZCBSYXRpbw=="><span class='gt_from_md'>Abbreviations: CI = Confidence Interval, HR = Hazard Ratio</span></span></td>
        </tr>
      </tfoot>
    </table>
    </div>
    # A tibble: 4 × 7
      term       estimate std.error statistic p.value conf.low conf.high
      <chr>         <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
    1 quantile2      1.19     0.388     0.442 0.659      0.555      2.54
    2 quantile3      3.27     0.374     3.17  0.00154    1.57       6.80
    3 quantile4      1.17     0.397     0.403 0.687      0.539      2.56
    4 CD4_0_2001     1.12     0.257     0.425 0.671      0.673      1.85
    <div id="ijkzmuahro" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
      <style>#ijkzmuahro table {
      font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
      -webkit-font-smoothing: antialiased;
      -moz-osx-font-smoothing: grayscale;
    }

    #ijkzmuahro thead, #ijkzmuahro tbody, #ijkzmuahro tfoot, #ijkzmuahro tr, #ijkzmuahro td, #ijkzmuahro th {
      border-style: none;
    }

    #ijkzmuahro p {
      margin: 0;
      padding: 0;
    }

    #ijkzmuahro .gt_table {
      display: table;
      border-collapse: collapse;
      line-height: normal;
      margin-left: auto;
      margin-right: auto;
      color: #333333;
      font-size: 16px;
      font-weight: normal;
      font-style: normal;
      background-color: #FFFFFF;
      width: auto;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #A8A8A8;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #A8A8A8;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
    }

    #ijkzmuahro .gt_caption {
      padding-top: 4px;
      padding-bottom: 4px;
    }

    #ijkzmuahro .gt_title {
      color: #333333;
      font-size: 125%;
      font-weight: initial;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-color: #FFFFFF;
      border-bottom-width: 0;
    }

    #ijkzmuahro .gt_subtitle {
      color: #333333;
      font-size: 85%;
      font-weight: initial;
      padding-top: 3px;
      padding-bottom: 5px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-color: #FFFFFF;
      border-top-width: 0;
    }

    #ijkzmuahro .gt_heading {
      background-color: #FFFFFF;
      text-align: center;
      border-bottom-color: #FFFFFF;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #ijkzmuahro .gt_bottom_border {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #ijkzmuahro .gt_col_headings {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #ijkzmuahro .gt_col_heading {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 6px;
      padding-left: 5px;
      padding-right: 5px;
      overflow-x: hidden;
    }

    #ijkzmuahro .gt_column_spanner_outer {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      padding-top: 0;
      padding-bottom: 0;
      padding-left: 4px;
      padding-right: 4px;
    }

    #ijkzmuahro .gt_column_spanner_outer:first-child {
      padding-left: 0;
    }

    #ijkzmuahro .gt_column_spanner_outer:last-child {
      padding-right: 0;
    }

    #ijkzmuahro .gt_column_spanner {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 5px;
      overflow-x: hidden;
      display: inline-block;
      width: 100%;
    }

    #ijkzmuahro .gt_spanner_row {
      border-bottom-style: hidden;
    }

    #ijkzmuahro .gt_group_heading {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      text-align: left;
    }

    #ijkzmuahro .gt_empty_group_heading {
      padding: 0.5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: middle;
    }

    #ijkzmuahro .gt_from_md > :first-child {
      margin-top: 0;
    }

    #ijkzmuahro .gt_from_md > :last-child {
      margin-bottom: 0;
    }

    #ijkzmuahro .gt_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      margin: 10px;
      border-top-style: solid;
      border-top-width: 1px;
      border-top-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      overflow-x: hidden;
    }

    #ijkzmuahro .gt_stub {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
    }

    #ijkzmuahro .gt_stub_row_group {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
      vertical-align: top;
    }

    #ijkzmuahro .gt_row_group_first td {
      border-top-width: 2px;
    }

    #ijkzmuahro .gt_row_group_first th {
      border-top-width: 2px;
    }

    #ijkzmuahro .gt_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #ijkzmuahro .gt_first_summary_row {
      border-top-style: solid;
      border-top-color: #D3D3D3;
    }

    #ijkzmuahro .gt_first_summary_row.thick {
      border-top-width: 2px;
    }

    #ijkzmuahro .gt_last_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #ijkzmuahro .gt_grand_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #ijkzmuahro .gt_first_grand_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-style: double;
      border-top-width: 6px;
      border-top-color: #D3D3D3;
    }

    #ijkzmuahro .gt_last_grand_summary_row_top {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: double;
      border-bottom-width: 6px;
      border-bottom-color: #D3D3D3;
    }

    #ijkzmuahro .gt_striped {
      background-color: rgba(128, 128, 128, 0.05);
    }

    #ijkzmuahro .gt_table_body {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #ijkzmuahro .gt_footnotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #ijkzmuahro .gt_footnote {
      margin: 0px;
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #ijkzmuahro .gt_sourcenotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #ijkzmuahro .gt_sourcenote {
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #ijkzmuahro .gt_left {
      text-align: left;
    }

    #ijkzmuahro .gt_center {
      text-align: center;
    }

    #ijkzmuahro .gt_right {
      text-align: right;
      font-variant-numeric: tabular-nums;
    }

    #ijkzmuahro .gt_font_normal {
      font-weight: normal;
    }

    #ijkzmuahro .gt_font_bold {
      font-weight: bold;
    }

    #ijkzmuahro .gt_font_italic {
      font-style: italic;
    }

    #ijkzmuahro .gt_super {
      font-size: 65%;
    }

    #ijkzmuahro .gt_footnote_marks {
      font-size: 75%;
      vertical-align: 0.4em;
      position: initial;
    }

    #ijkzmuahro .gt_asterisk {
      font-size: 100%;
      vertical-align: 0;
    }

    #ijkzmuahro .gt_indent_1 {
      text-indent: 5px;
    }

    #ijkzmuahro .gt_indent_2 {
      text-indent: 10px;
    }

    #ijkzmuahro .gt_indent_3 {
      text-indent: 15px;
    }

    #ijkzmuahro .gt_indent_4 {
      text-indent: 20px;
    }

    #ijkzmuahro .gt_indent_5 {
      text-indent: 25px;
    }

    #ijkzmuahro .katex-display {
      display: inline-flex !important;
      margin-bottom: 0.75em !important;
    }

    #ijkzmuahro div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
      height: 0px !important;
    }
    </style>
      <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
      <thead>
        <tr class="gt_col_headings">
          <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="label"><span data-qmd-base64="KipDaGFyYWN0ZXJpc3RpYyoq"><span class='gt_from_md'><strong>Characteristic</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipIUioq"><span class='gt_from_md'><strong>HR</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="conf.low"><span data-qmd-base64="Kio5NSUgQ0kqKg=="><span class='gt_from_md'><strong>95% CI</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="p.value"><span data-qmd-base64="KipwLXZhbHVlKio="><span class='gt_from_md'><strong>p-value</strong></span></span></th>
        </tr>
      </thead>
      <tbody class="gt_table_body">
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">DAPP1 Quantile</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    2</td>
    <td headers="estimate" class="gt_row gt_center">1.19</td>
    <td headers="conf.low" class="gt_row gt_center">0.55, 2.54</td>
    <td headers="p.value" class="gt_row gt_center">0.7</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    3</td>
    <td headers="estimate" class="gt_row gt_center">3.27</td>
    <td headers="conf.low" class="gt_row gt_center">1.57, 6.80</td>
    <td headers="p.value" class="gt_row gt_center" style="font-weight: bold;">0.002</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    4</td>
    <td headers="estimate" class="gt_row gt_center">1.17</td>
    <td headers="conf.low" class="gt_row gt_center">0.54, 2.56</td>
    <td headers="p.value" class="gt_row gt_center">0.7</td></tr>
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">CD4_0_200</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    0</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">1.12</td>
    <td headers="conf.low" class="gt_row gt_center">0.67, 1.85</td>
    <td headers="p.value" class="gt_row gt_center">0.7</td></tr>
      </tbody>
      <tfoot>
        <tr class="gt_sourcenotes">
          <td class="gt_sourcenote" colspan="4"><span data-qmd-base64="QWJicmV2aWF0aW9uczogQ0kgPSBDb25maWRlbmNlIEludGVydmFsLCBIUiA9IEhhemFyZCBSYXRpbw=="><span class='gt_from_md'>Abbreviations: CI = Confidence Interval, HR = Hazard Ratio</span></span></td>
        </tr>
      </tfoot>
    </table>
    </div>
    # A tibble: 3 × 7
      term      estimate std.error statistic p.value conf.low conf.high
      <chr>        <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
    1 quantile2     1.16     0.344     0.424   0.672    0.589      2.27
    2 quantile3     1.12     0.349     0.330   0.741    0.566      2.22
    3 quantile4     1.09     0.349     0.243   0.808    0.549      2.16
    <div id="tpsoutpwee" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
      <style>#tpsoutpwee table {
      font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
      -webkit-font-smoothing: antialiased;
      -moz-osx-font-smoothing: grayscale;
    }

    #tpsoutpwee thead, #tpsoutpwee tbody, #tpsoutpwee tfoot, #tpsoutpwee tr, #tpsoutpwee td, #tpsoutpwee th {
      border-style: none;
    }

    #tpsoutpwee p {
      margin: 0;
      padding: 0;
    }

    #tpsoutpwee .gt_table {
      display: table;
      border-collapse: collapse;
      line-height: normal;
      margin-left: auto;
      margin-right: auto;
      color: #333333;
      font-size: 16px;
      font-weight: normal;
      font-style: normal;
      background-color: #FFFFFF;
      width: auto;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #A8A8A8;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #A8A8A8;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
    }

    #tpsoutpwee .gt_caption {
      padding-top: 4px;
      padding-bottom: 4px;
    }

    #tpsoutpwee .gt_title {
      color: #333333;
      font-size: 125%;
      font-weight: initial;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-color: #FFFFFF;
      border-bottom-width: 0;
    }

    #tpsoutpwee .gt_subtitle {
      color: #333333;
      font-size: 85%;
      font-weight: initial;
      padding-top: 3px;
      padding-bottom: 5px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-color: #FFFFFF;
      border-top-width: 0;
    }

    #tpsoutpwee .gt_heading {
      background-color: #FFFFFF;
      text-align: center;
      border-bottom-color: #FFFFFF;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #tpsoutpwee .gt_bottom_border {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #tpsoutpwee .gt_col_headings {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #tpsoutpwee .gt_col_heading {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 6px;
      padding-left: 5px;
      padding-right: 5px;
      overflow-x: hidden;
    }

    #tpsoutpwee .gt_column_spanner_outer {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      padding-top: 0;
      padding-bottom: 0;
      padding-left: 4px;
      padding-right: 4px;
    }

    #tpsoutpwee .gt_column_spanner_outer:first-child {
      padding-left: 0;
    }

    #tpsoutpwee .gt_column_spanner_outer:last-child {
      padding-right: 0;
    }

    #tpsoutpwee .gt_column_spanner {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 5px;
      overflow-x: hidden;
      display: inline-block;
      width: 100%;
    }

    #tpsoutpwee .gt_spanner_row {
      border-bottom-style: hidden;
    }

    #tpsoutpwee .gt_group_heading {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      text-align: left;
    }

    #tpsoutpwee .gt_empty_group_heading {
      padding: 0.5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: middle;
    }

    #tpsoutpwee .gt_from_md > :first-child {
      margin-top: 0;
    }

    #tpsoutpwee .gt_from_md > :last-child {
      margin-bottom: 0;
    }

    #tpsoutpwee .gt_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      margin: 10px;
      border-top-style: solid;
      border-top-width: 1px;
      border-top-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      overflow-x: hidden;
    }

    #tpsoutpwee .gt_stub {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
    }

    #tpsoutpwee .gt_stub_row_group {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
      vertical-align: top;
    }

    #tpsoutpwee .gt_row_group_first td {
      border-top-width: 2px;
    }

    #tpsoutpwee .gt_row_group_first th {
      border-top-width: 2px;
    }

    #tpsoutpwee .gt_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #tpsoutpwee .gt_first_summary_row {
      border-top-style: solid;
      border-top-color: #D3D3D3;
    }

    #tpsoutpwee .gt_first_summary_row.thick {
      border-top-width: 2px;
    }

    #tpsoutpwee .gt_last_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #tpsoutpwee .gt_grand_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #tpsoutpwee .gt_first_grand_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-style: double;
      border-top-width: 6px;
      border-top-color: #D3D3D3;
    }

    #tpsoutpwee .gt_last_grand_summary_row_top {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: double;
      border-bottom-width: 6px;
      border-bottom-color: #D3D3D3;
    }

    #tpsoutpwee .gt_striped {
      background-color: rgba(128, 128, 128, 0.05);
    }

    #tpsoutpwee .gt_table_body {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #tpsoutpwee .gt_footnotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #tpsoutpwee .gt_footnote {
      margin: 0px;
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #tpsoutpwee .gt_sourcenotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #tpsoutpwee .gt_sourcenote {
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #tpsoutpwee .gt_left {
      text-align: left;
    }

    #tpsoutpwee .gt_center {
      text-align: center;
    }

    #tpsoutpwee .gt_right {
      text-align: right;
      font-variant-numeric: tabular-nums;
    }

    #tpsoutpwee .gt_font_normal {
      font-weight: normal;
    }

    #tpsoutpwee .gt_font_bold {
      font-weight: bold;
    }

    #tpsoutpwee .gt_font_italic {
      font-style: italic;
    }

    #tpsoutpwee .gt_super {
      font-size: 65%;
    }

    #tpsoutpwee .gt_footnote_marks {
      font-size: 75%;
      vertical-align: 0.4em;
      position: initial;
    }

    #tpsoutpwee .gt_asterisk {
      font-size: 100%;
      vertical-align: 0;
    }

    #tpsoutpwee .gt_indent_1 {
      text-indent: 5px;
    }

    #tpsoutpwee .gt_indent_2 {
      text-indent: 10px;
    }

    #tpsoutpwee .gt_indent_3 {
      text-indent: 15px;
    }

    #tpsoutpwee .gt_indent_4 {
      text-indent: 20px;
    }

    #tpsoutpwee .gt_indent_5 {
      text-indent: 25px;
    }

    #tpsoutpwee .katex-display {
      display: inline-flex !important;
      margin-bottom: 0.75em !important;
    }

    #tpsoutpwee div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
      height: 0px !important;
    }
    </style>
      <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
      <thead>
        <tr class="gt_col_headings">
          <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="label"><span data-qmd-base64="KipDaGFyYWN0ZXJpc3RpYyoq"><span class='gt_from_md'><strong>Characteristic</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipIUioq"><span class='gt_from_md'><strong>HR</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="conf.low"><span data-qmd-base64="Kio5NSUgQ0kqKg=="><span class='gt_from_md'><strong>95% CI</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="p.value"><span data-qmd-base64="KipwLXZhbHVlKio="><span class='gt_from_md'><strong>p-value</strong></span></span></th>
        </tr>
      </thead>
      <tbody class="gt_table_body">
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">SKAP2 Quantile</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    2</td>
    <td headers="estimate" class="gt_row gt_center">1.16</td>
    <td headers="conf.low" class="gt_row gt_center">0.59, 2.27</td>
    <td headers="p.value" class="gt_row gt_center">0.7</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    3</td>
    <td headers="estimate" class="gt_row gt_center">1.12</td>
    <td headers="conf.low" class="gt_row gt_center">0.57, 2.22</td>
    <td headers="p.value" class="gt_row gt_center">0.7</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    4</td>
    <td headers="estimate" class="gt_row gt_center">1.09</td>
    <td headers="conf.low" class="gt_row gt_center">0.55, 2.16</td>
    <td headers="p.value" class="gt_row gt_center">0.8</td></tr>
      </tbody>
      <tfoot>
        <tr class="gt_sourcenotes">
          <td class="gt_sourcenote" colspan="4"><span data-qmd-base64="QWJicmV2aWF0aW9uczogQ0kgPSBDb25maWRlbmNlIEludGVydmFsLCBIUiA9IEhhemFyZCBSYXRpbw=="><span class='gt_from_md'>Abbreviations: CI = Confidence Interval, HR = Hazard Ratio</span></span></td>
        </tr>
      </tfoot>
    </table>
    </div>
    # A tibble: 4 × 7
      term       estimate std.error statistic p.value conf.low conf.high
      <chr>         <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
    1 quantile2      1.21     0.368     0.527   0.598    0.591      2.49
    2 quantile3      1.28     0.367     0.670   0.503    0.623      2.62
    3 quantile4      1.18     0.361     0.456   0.649    0.581      2.39
    4 CD4_0_2001     1.14     0.260     0.518   0.604    0.687      1.91
    <div id="pmntrcricx" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
      <style>#pmntrcricx table {
      font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
      -webkit-font-smoothing: antialiased;
      -moz-osx-font-smoothing: grayscale;
    }

    #pmntrcricx thead, #pmntrcricx tbody, #pmntrcricx tfoot, #pmntrcricx tr, #pmntrcricx td, #pmntrcricx th {
      border-style: none;
    }

    #pmntrcricx p {
      margin: 0;
      padding: 0;
    }

    #pmntrcricx .gt_table {
      display: table;
      border-collapse: collapse;
      line-height: normal;
      margin-left: auto;
      margin-right: auto;
      color: #333333;
      font-size: 16px;
      font-weight: normal;
      font-style: normal;
      background-color: #FFFFFF;
      width: auto;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #A8A8A8;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #A8A8A8;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
    }

    #pmntrcricx .gt_caption {
      padding-top: 4px;
      padding-bottom: 4px;
    }

    #pmntrcricx .gt_title {
      color: #333333;
      font-size: 125%;
      font-weight: initial;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-color: #FFFFFF;
      border-bottom-width: 0;
    }

    #pmntrcricx .gt_subtitle {
      color: #333333;
      font-size: 85%;
      font-weight: initial;
      padding-top: 3px;
      padding-bottom: 5px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-color: #FFFFFF;
      border-top-width: 0;
    }

    #pmntrcricx .gt_heading {
      background-color: #FFFFFF;
      text-align: center;
      border-bottom-color: #FFFFFF;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #pmntrcricx .gt_bottom_border {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #pmntrcricx .gt_col_headings {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
    }

    #pmntrcricx .gt_col_heading {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 6px;
      padding-left: 5px;
      padding-right: 5px;
      overflow-x: hidden;
    }

    #pmntrcricx .gt_column_spanner_outer {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: normal;
      text-transform: inherit;
      padding-top: 0;
      padding-bottom: 0;
      padding-left: 4px;
      padding-right: 4px;
    }

    #pmntrcricx .gt_column_spanner_outer:first-child {
      padding-left: 0;
    }

    #pmntrcricx .gt_column_spanner_outer:last-child {
      padding-right: 0;
    }

    #pmntrcricx .gt_column_spanner {
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: bottom;
      padding-top: 5px;
      padding-bottom: 5px;
      overflow-x: hidden;
      display: inline-block;
      width: 100%;
    }

    #pmntrcricx .gt_spanner_row {
      border-bottom-style: hidden;
    }

    #pmntrcricx .gt_group_heading {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      text-align: left;
    }

    #pmntrcricx .gt_empty_group_heading {
      padding: 0.5px;
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      vertical-align: middle;
    }

    #pmntrcricx .gt_from_md > :first-child {
      margin-top: 0;
    }

    #pmntrcricx .gt_from_md > :last-child {
      margin-bottom: 0;
    }

    #pmntrcricx .gt_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      margin: 10px;
      border-top-style: solid;
      border-top-width: 1px;
      border-top-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 1px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 1px;
      border-right-color: #D3D3D3;
      vertical-align: middle;
      overflow-x: hidden;
    }

    #pmntrcricx .gt_stub {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
    }

    #pmntrcricx .gt_stub_row_group {
      color: #333333;
      background-color: #FFFFFF;
      font-size: 100%;
      font-weight: initial;
      text-transform: inherit;
      border-right-style: solid;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
      padding-left: 5px;
      padding-right: 5px;
      vertical-align: top;
    }

    #pmntrcricx .gt_row_group_first td {
      border-top-width: 2px;
    }

    #pmntrcricx .gt_row_group_first th {
      border-top-width: 2px;
    }

    #pmntrcricx .gt_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #pmntrcricx .gt_first_summary_row {
      border-top-style: solid;
      border-top-color: #D3D3D3;
    }

    #pmntrcricx .gt_first_summary_row.thick {
      border-top-width: 2px;
    }

    #pmntrcricx .gt_last_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #pmntrcricx .gt_grand_summary_row {
      color: #333333;
      background-color: #FFFFFF;
      text-transform: inherit;
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #pmntrcricx .gt_first_grand_summary_row {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-top-style: double;
      border-top-width: 6px;
      border-top-color: #D3D3D3;
    }

    #pmntrcricx .gt_last_grand_summary_row_top {
      padding-top: 8px;
      padding-bottom: 8px;
      padding-left: 5px;
      padding-right: 5px;
      border-bottom-style: double;
      border-bottom-width: 6px;
      border-bottom-color: #D3D3D3;
    }

    #pmntrcricx .gt_striped {
      background-color: rgba(128, 128, 128, 0.05);
    }

    #pmntrcricx .gt_table_body {
      border-top-style: solid;
      border-top-width: 2px;
      border-top-color: #D3D3D3;
      border-bottom-style: solid;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
    }

    #pmntrcricx .gt_footnotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #pmntrcricx .gt_footnote {
      margin: 0px;
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #pmntrcricx .gt_sourcenotes {
      color: #333333;
      background-color: #FFFFFF;
      border-bottom-style: none;
      border-bottom-width: 2px;
      border-bottom-color: #D3D3D3;
      border-left-style: none;
      border-left-width: 2px;
      border-left-color: #D3D3D3;
      border-right-style: none;
      border-right-width: 2px;
      border-right-color: #D3D3D3;
    }

    #pmntrcricx .gt_sourcenote {
      font-size: 90%;
      padding-top: 4px;
      padding-bottom: 4px;
      padding-left: 5px;
      padding-right: 5px;
    }

    #pmntrcricx .gt_left {
      text-align: left;
    }

    #pmntrcricx .gt_center {
      text-align: center;
    }

    #pmntrcricx .gt_right {
      text-align: right;
      font-variant-numeric: tabular-nums;
    }

    #pmntrcricx .gt_font_normal {
      font-weight: normal;
    }

    #pmntrcricx .gt_font_bold {
      font-weight: bold;
    }

    #pmntrcricx .gt_font_italic {
      font-style: italic;
    }

    #pmntrcricx .gt_super {
      font-size: 65%;
    }

    #pmntrcricx .gt_footnote_marks {
      font-size: 75%;
      vertical-align: 0.4em;
      position: initial;
    }

    #pmntrcricx .gt_asterisk {
      font-size: 100%;
      vertical-align: 0;
    }

    #pmntrcricx .gt_indent_1 {
      text-indent: 5px;
    }

    #pmntrcricx .gt_indent_2 {
      text-indent: 10px;
    }

    #pmntrcricx .gt_indent_3 {
      text-indent: 15px;
    }

    #pmntrcricx .gt_indent_4 {
      text-indent: 20px;
    }

    #pmntrcricx .gt_indent_5 {
      text-indent: 25px;
    }

    #pmntrcricx .katex-display {
      display: inline-flex !important;
      margin-bottom: 0.75em !important;
    }

    #pmntrcricx div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
      height: 0px !important;
    }
    </style>
      <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
      <thead>
        <tr class="gt_col_headings">
          <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="label"><span data-qmd-base64="KipDaGFyYWN0ZXJpc3RpYyoq"><span class='gt_from_md'><strong>Characteristic</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="estimate"><span data-qmd-base64="KipIUioq"><span class='gt_from_md'><strong>HR</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="conf.low"><span data-qmd-base64="Kio5NSUgQ0kqKg=="><span class='gt_from_md'><strong>95% CI</strong></span></span></th>
          <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="p.value"><span data-qmd-base64="KipwLXZhbHVlKio="><span class='gt_from_md'><strong>p-value</strong></span></span></th>
        </tr>
      </thead>
      <tbody class="gt_table_body">
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">SKAP2 Quantile</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    2</td>
    <td headers="estimate" class="gt_row gt_center">1.21</td>
    <td headers="conf.low" class="gt_row gt_center">0.59, 2.49</td>
    <td headers="p.value" class="gt_row gt_center">0.6</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    3</td>
    <td headers="estimate" class="gt_row gt_center">1.28</td>
    <td headers="conf.low" class="gt_row gt_center">0.62, 2.62</td>
    <td headers="p.value" class="gt_row gt_center">0.5</td></tr>
        <tr><td headers="label" class="gt_row gt_left">    4</td>
    <td headers="estimate" class="gt_row gt_center">1.18</td>
    <td headers="conf.low" class="gt_row gt_center">0.58, 2.39</td>
    <td headers="p.value" class="gt_row gt_center">0.6</td></tr>
        <tr><td headers="label" class="gt_row gt_left" style="font-weight: bold;">CD4_0_200</td>
    <td headers="estimate" class="gt_row gt_center"><br /></td>
    <td headers="conf.low" class="gt_row gt_center"><br /></td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    0</td>
    <td headers="estimate" class="gt_row gt_center">—</td>
    <td headers="conf.low" class="gt_row gt_center">—</td>
    <td headers="p.value" class="gt_row gt_center"><br /></td></tr>
        <tr><td headers="label" class="gt_row gt_left">    1</td>
    <td headers="estimate" class="gt_row gt_center">1.14</td>
    <td headers="conf.low" class="gt_row gt_center">0.69, 1.91</td>
    <td headers="p.value" class="gt_row gt_center">0.6</td></tr>
      </tbody>
      <tfoot>
        <tr class="gt_sourcenotes">
          <td class="gt_sourcenote" colspan="4"><span data-qmd-base64="QWJicmV2aWF0aW9uczogQ0kgPSBDb25maWRlbmNlIEludGVydmFsLCBIUiA9IEhhemFyZCBSYXRpbw=="><span class='gt_from_md'>Abbreviations: CI = Confidence Interval, HR = Hazard Ratio</span></span></td>
        </tr>
      </tfoot>
    </table>
    </div>

## Crude and adjusted hazard ratios by outcome (Table 2 / Tables S12-S14)

One table per outcome, matching the manuscript’s published layout: a
Ref. row for Q1, per-quartile crude/adjusted HR (95% CI), number of
events per quartile, and a Global P value row (from the linear-trend
models).

``` r
format_hr <- function(estimate, low, high) {
    ifelse(is.na(estimate), NA_character_,
           ifelse(is.na(low), "Ref.", sprintf("%.2f (%.2f-%.2f)", estimate, low, high)))
}

build_outcome_table <- function(outcome) {
    quartile_rows <- purrr::map_dfr(markers, function(marker) {
        crude <- tidy(categorical_models[[marker]][[outcome]]$crude, exponentiate = TRUE, conf.int = TRUE) |>
            filter(str_detect(term, "^quantile")) |>
            mutate(quartile = str_remove(term, "quantile"), crude_estimate = estimate,
                   crude_low = conf.low, crude_high = conf.high, .keep = "none")
        adj <- tidy(categorical_models[[marker]][[outcome]]$adjusted, exponentiate = TRUE, conf.int = TRUE) |>
            filter(str_detect(term, "^quantile")) |>
            mutate(quartile = str_remove(term, "quantile"), adj_estimate = estimate,
                   adj_low = conf.low, adj_high = conf.high, .keep = "none")
        ref <- tibble(quartile = "1", crude_estimate = 1, crude_low = NA, crude_high = NA,
                      adj_estimate = 1, adj_low = NA, adj_high = NA)
        n_events <- marker_event_counts[[marker]][[outcome]] |>
            mutate(quartile = as.character(quantile), n_events = n, .keep = "none")

        bind_rows(ref, crude |> left_join(adj, by = "quartile")) |>
            left_join(n_events, by = "quartile") |>
            mutate(
                biomarker = marker, quartile, n_events,
                HR = format_hr(crude_estimate, crude_low, crude_high),
                Adjusted_HR = format_hr(adj_estimate, adj_low, adj_high),
                .keep = "none"
            )
    })

    global_p <- purrr::map_dfr(markers, function(marker) {
        tibble(
            biomarker = marker, quartile = "Global P value", n_events = NA_integer_,
            HR = formatC(tidy(trend_models[[marker]][[outcome]]$crude)$p.value[1], digits = 3, format = "f"),
            Adjusted_HR = formatC(tidy(trend_models[[marker]][[outcome]]$adjusted)$p.value[1], digits = 3, format = "f")
        )
    })

    bind_rows(quartile_rows, global_p) |>
        mutate(quartile = factor(quartile, levels = c("1", "2", "3", "4", "Global P value"))) |>
        pivot_wider(
            id_cols = quartile,
            names_from = biomarker,
            values_from = c(n_events, HR, Adjusted_HR),
            names_glue = "{biomarker}_{.value}"
        ) |>
        arrange(quartile) |>
        dplyr::select(quartile, starts_with("CCL13"), starts_with("DAPP1"),
                      starts_with("LGALS9"), starts_with("SKAP2"))
}
```

``` r
table2 <- build_outcome_table("group")
table2
```

    # A tibble: 5 × 13
      quartile     CCL13_n_events CCL13_HR CCL13_Adjusted_HR DAPP1_n_events DAPP1_HR
      <fct>                 <int> <chr>    <chr>                      <int> <chr>   
    1 1                        12 Ref.     Ref.                          14 Ref.    
    2 2                        17 1.45 (0… 1.36 (0.61-3.00)              17 1.08 (0…
    3 3                        22 2.35 (1… 2.39 (1.15-4.97)              22 2.99 (1…
    4 4                        16 1.30 (0… 1.34 (0.62-2.89)              14 1.05 (0…
    5 Global P va…             NA 0.306    0.271                         NA 0.356   
    # ℹ 7 more variables: DAPP1_Adjusted_HR <chr>, LGALS9_n_events <int>,
    #   LGALS9_HR <chr>, LGALS9_Adjusted_HR <chr>, SKAP2_n_events <int>,
    #   SKAP2_HR <chr>, SKAP2_Adjusted_HR <chr>

``` r
# table2 |> write_tsv("../results/table2.tsv")
```

``` r
tableS12 <- build_outcome_table("nadm")
tableS12
```

    # A tibble: 5 × 13
      quartile     CCL13_n_events CCL13_HR CCL13_Adjusted_HR DAPP1_n_events DAPP1_HR
      <fct>                 <int> <chr>    <chr>                      <int> <chr>   
    1 1                         9 Ref.     Ref.                          12 Ref.    
    2 2                        11 1.23 (0… 1.28 (0.50-3.28)              12 0.88 (0…
    3 3                        16 2.27 (1… 2.32 (0.98-5.47)              14 2.34 (1…
    4 4                        10 1.08 (0… 1.14 (0.45-2.89)               7 0.62 (0…
    5 Global P va…             NA 0.548    0.527                         NA 0.807   
    # ℹ 7 more variables: DAPP1_Adjusted_HR <chr>, LGALS9_n_events <int>,
    #   LGALS9_HR <chr>, LGALS9_Adjusted_HR <chr>, SKAP2_n_events <int>,
    #   SKAP2_HR <chr>, SKAP2_Adjusted_HR <chr>

``` r
# tableS12 |> write_tsv("../results/tableS12_neoplasia.tsv")
```

``` r
tableS13 <- build_outcome_table("mace")
tableS13
```

    # A tibble: 5 × 13
      quartile     CCL13_n_events CCL13_HR CCL13_Adjusted_HR DAPP1_n_events DAPP1_HR
      <fct>                 <int> <chr>    <chr>                      <int> <chr>   
    1 1                         3 Ref.     Ref.                           2 Ref.    
    2 2                         6 2.13 (0… 1.54 (0.34-6.94)               5 2.22 (0…
    3 3                         6 2.59 (0… 2.58 (0.64-10.34)              8 6.74 (1…
    4 4                         6 2.01 (0… 1.90 (0.48-7.62)               7 3.53 (0…
    5 Global P va…             NA 0.347    0.297                         NA 0.055   
    # ℹ 7 more variables: DAPP1_Adjusted_HR <chr>, LGALS9_n_events <int>,
    #   LGALS9_HR <chr>, LGALS9_Adjusted_HR <chr>, SKAP2_n_events <int>,
    #   SKAP2_HR <chr>, SKAP2_Adjusted_HR <chr>

``` r
# tableS13 |> write_tsv("../results/tableS13_cardiovascular.tsv")
```

``` r
tableS14 <- build_outcome_table("death_all")
tableS14
```

    # A tibble: 5 × 13
      quartile     CCL13_n_events CCL13_HR CCL13_Adjusted_HR DAPP1_n_events DAPP1_HR
      <fct>                 <int> <chr>    <chr>                      <int> <chr>   
    1 1                         4 Ref.     Ref.                           4 Ref.    
    2 2                         1 0.28 (0… 0.36 (0.04-3.47)               3 0.66 (0…
    3 3                         5 1.55 (0… 2.02 (0.48-8.47)               6 2.53 (0…
    4 4                         5 1.26 (0… 1.54 (0.37-6.48)               2 0.52 (0…
    5 Global P va…             NA 0.402    0.275                         NA 0.942   
    # ℹ 7 more variables: DAPP1_Adjusted_HR <chr>, LGALS9_n_events <int>,
    #   LGALS9_HR <chr>, LGALS9_Adjusted_HR <chr>, SKAP2_n_events <int>,
    #   SKAP2_HR <chr>, SKAP2_Adjusted_HR <chr>

``` r
# tableS14 |> write_tsv("../results/tableS14_mortality.tsv")
```

## Schoenfeld residuals (Table S16)

``` r
model_schoenfeld <- function(fit, transform = "km") {
    ph <- cox.zph(fit, transform = transform)
    df <- data.frame(time = ph$x, resid = ph$y[, 1])
    list(ph = ph, df = df)
}

schoenfeld_models <- purrr::imap(categorical_models, ~ model_schoenfeld(.x$group$crude))

tableS16 <- purrr::imap_dfr(schoenfeld_models, function(model, marker) {
    tibble(
        biomarker = marker,
        quantile_p = model$ph$table[1, "p"],
        global_p   = model$ph$table["GLOBAL", "p"]
    )
}) |>
    pivot_longer(cols = c(quantile_p, global_p), names_to = "term", values_to = "p_value") |>
    mutate(term = recode(term, quantile_p = "quantile (Q1 to Q4)", global_p = "Global")) |>
    pivot_wider(names_from = biomarker, values_from = p_value) |>
    select(term, CCL13, DAPP1, LGALS9, SKAP2)

tableS16
```

    # A tibble: 2 × 5
      term                 CCL13 DAPP1 LGALS9 SKAP2
      <chr>                <dbl> <dbl>  <dbl> <dbl>
    1 quantile (Q1 to Q4) 0.0614 0.499  0.986 0.551
    2 Global              0.0614 0.499  0.986 0.551

``` r
# tableS16 |> write_tsv("../results/tableS16_schoenfeld.tsv")
```

# Plots

## Kaplan-Meier figures

``` r
spread_y <- function(y, min_sep = 0.03, y_min = 0, y_max = 1) {
    if (length(y) <= 1) return(y)

    ord <- order(y)
    ys  <- y[ord]

    for (i in 2:length(ys)) {
        if (ys[i] - ys[i - 1] < min_sep) ys[i] <- ys[i - 1] + min_sep
    }

    overflow <- max(ys) - y_max
    if (overflow > 0) {
        ys <- ys - overflow
        for (i in (length(ys) - 1):1) {
            if (ys[i + 1] - ys[i] < min_sep) ys[i] <- ys[i + 1] - min_sep
        }
    }

    ys <- pmin(pmax(ys, y_min), y_max)

    out <- y
    out[ord] <- ys
    out
}

plot_cuminc_pub <- function(fit,
                            title,
                            xlab = "Months",
                            ylab = "Time-to-event prob.",
                            risktable_stats = "cum.event",
                            x_nudge = 0.6,
                            min_sep = 0.03,
                            right_expand = 0.2,
                            base_size = 10,
                            p_text = NULL) {

    p <- ggcuminc(fit) +
        labs(x = xlab, y = ylab, color = "Quartile") +
        ggtitle(title) +
        theme_bw(base_size = base_size) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0),
            axis.title = element_text(face = "bold"),
            legend.position = "none"
        )

    lab_df <-
        p$data |>
        group_by(strata) |>
        slice_max(time, n = 1, with_ties = FALSE) |>
        ungroup() |>
        mutate(
            strata_chr = as.character(strata),
            qnum  = str_extract(strata_chr, "\\d+"),
            label = paste0("Q", qnum),
            x_label = time + x_nudge
        )

    lab_df <- lab_df |>
        mutate(y_label = spread_y(estimate, min_sep = min_sep, y_min = 0, y_max = 1))

    p <- p +
        geom_text(
            data = lab_df,
            aes(x = x_label, y = y_label, label = label, color = strata),
            inherit.aes = FALSE,
            fontface = "bold",
            hjust = 0,
            size = 3,
            show.legend = FALSE
        ) +
        coord_cartesian(clip = "off") +
        scale_x_continuous(expand = expansion(mult = c(0.02, right_expand)))

    if (!is.null(p_text)) {
        p <- p + annotate("text", x = Inf, y = 0.02, label = p_text,
                          hjust = 1.1, vjust = 0, size = 3.2, fontface = "plain")
    }

    p <- p + add_risktable(risktable_stats = risktable_stats)
    ggsurvfit_build(p)
}

build_4panel_km <- function(p1, p2, p3, p4) {
    w1 <- wrap_elements(full = p1); w2 <- wrap_elements(full = p2)
    w3 <- wrap_elements(full = p3); w4 <- wrap_elements(full = p4)

    (w1 | w2) / (w3 | w4) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(face = "bold", size = 12),
              plot.tag.position = c(0.01, 0.99))
}

p_text_for <- function(marker, outcome_name, adjusted = FALSE) {
    p <- tidy(trend_models[[marker]][[outcome_name]][[if (adjusted) "adjusted" else "crude"]])$p.value[1]
    paste0("Global P value = ", formatC(p, digits = 3, format = "f"))
}
```

### Main comparison (Fig. 5)

``` r
fits_group <- purrr::map(marker_data, ~ survfit2(Surv(time, as.factor(group)) ~ quantile, data = .x))

fig5 <- build_4panel_km(
    plot_cuminc_pub(fits_group$LGALS9, "LGALS9", p_text = p_text_for("LGALS9", "group")),
    plot_cuminc_pub(fits_group$CCL13,  "CCL13",  p_text = p_text_for("CCL13", "group")),
    plot_cuminc_pub(fits_group$DAPP1,  "DAPP1",  p_text = p_text_for("DAPP1", "group")),
    plot_cuminc_pub(fits_group$SKAP2,  "SKAP2",  p_text = p_text_for("SKAP2", "group"))
)
```

    Plotting outcome "1".
    Plotting outcome "1".
    Plotting outcome "1".
    Plotting outcome "1".

``` r
# svg("../results/fig5_KM_main.svg", width = 8, height = 8)
fig5
```

![](07_survival_analysis_files/figure-commonmark/unnamed-chunk-16-1.png)

``` r
# dev.off()
```

### Outcome-specific comparisons (Figs. S7-S9)

``` r
fits_nadm  <- purrr::map(marker_data, ~ survfit2(Surv(time, as.factor(nadm)) ~ quantile, data = .x))

figS7 <- build_4panel_km(
    plot_cuminc_pub(fits_nadm$LGALS9, "LGALS9", p_text = p_text_for("LGALS9", "nadm")),
    plot_cuminc_pub(fits_nadm$CCL13,  "CCL13",  p_text = p_text_for("CCL13", "nadm")),
    plot_cuminc_pub(fits_nadm$DAPP1,  "DAPP1",  p_text = p_text_for("DAPP1", "nadm")),
    plot_cuminc_pub(fits_nadm$SKAP2,  "SKAP2",  p_text = p_text_for("SKAP2", "nadm"))
)
```

    Plotting outcome "1".
    Plotting outcome "1".
    Plotting outcome "1".
    Plotting outcome "1".

``` r
# svg("../results/figS7_KM_nadm.svg", width = 8, height = 8)
figS7
```

![](07_survival_analysis_files/figure-commonmark/unnamed-chunk-17-1.png)

``` r
# dev.off()
```

``` r
fits_mace  <- purrr::map(marker_data, ~ survfit2(Surv(time, as.factor(mace)) ~ quantile, data = .x))

figS8 <- build_4panel_km(
    plot_cuminc_pub(fits_mace$LGALS9, "LGALS9", p_text = p_text_for("LGALS9", "mace")),
    plot_cuminc_pub(fits_mace$CCL13,  "CCL13",  p_text = p_text_for("CCL13", "mace")),
    plot_cuminc_pub(fits_mace$DAPP1,  "DAPP1",  p_text = p_text_for("DAPP1", "mace")),
    plot_cuminc_pub(fits_mace$SKAP2,  "SKAP2",  p_text = p_text_for("SKAP2", "mace"))
)
```

    Plotting outcome "1".
    Plotting outcome "1".
    Plotting outcome "1".
    Plotting outcome "1".

``` r
# svg("../results/figS8_KM_mace.svg", width = 8, height = 8)
figS8
```

![](07_survival_analysis_files/figure-commonmark/unnamed-chunk-18-1.png)

``` r
# dev.off()
```

``` r
fits_death <- purrr::map(marker_data, ~ survfit2(Surv(time, as.factor(death_all)) ~ quantile, data = .x))

figS9 <- build_4panel_km(
    plot_cuminc_pub(fits_death$LGALS9, "LGALS9", p_text = p_text_for("LGALS9", "death_all")),
    plot_cuminc_pub(fits_death$CCL13,  "CCL13",  p_text = p_text_for("CCL13", "death_all")),
    plot_cuminc_pub(fits_death$DAPP1,  "DAPP1",  p_text = p_text_for("DAPP1", "death_all")),
    plot_cuminc_pub(fits_death$SKAP2,  "SKAP2",  p_text = p_text_for("SKAP2", "death_all"))
)
```

    Plotting outcome "1".
    Plotting outcome "1".
    Plotting outcome "1".
    Plotting outcome "1".

``` r
# svg("../results/figS9_KM_death.svg", width = 8, height = 8)
figS9
```

![](07_survival_analysis_files/figure-commonmark/unnamed-chunk-19-1.png)

``` r
# dev.off()
```

## Non-linearity: penalized spline models (Fig. S6)

``` r
fit_spline_plot <- function(data, marker) {
    data <- data |>
        mutate(log2_conc = log2(mean_conc)) |>
        filter(!is.na(log2_conc), !is.na(time), !is.na(group))

    cox_pspline <- coxph(Surv(time, as.numeric(group)) ~ pspline(log2_conc, df = 3), data = data)

    newdat <- data.frame(log2_conc = seq(min(data$log2_conc), max(data$log2_conc), length = 200))
    pred <- predict(cox_pspline, newdata = newdat, type = "risk", se.fit = TRUE)

    plotdat <- newdat |>
        mutate(HR = exp(pred$fit), HR_low = exp(pred$fit - 1.96 * pred$se.fit),
               HR_high = exp(pred$fit + 1.96 * pred$se.fit))
    q_vals <- quantile(data$log2_conc, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)

    ggplot(plotdat, aes(x = log2_conc, y = HR)) +
        geom_ribbon(aes(ymin = HR_low, ymax = HR_high), alpha = 0.25) +
        geom_line(linewidth = 1.2) +
        geom_vline(xintercept = q_vals, linetype = "dashed", color = "black") +
        geom_rug(data = data, aes(x = log2_conc), sides = "b", alpha = 0.4, inherit.aes = FALSE) +
        scale_x_continuous(paste0("log2 ", marker)) +
        scale_y_continuous("Relative Hazard") +
        theme_minimal(base_size = 14) +
        ggtitle(marker)
}

spline_plots <- purrr::imap(marker_data, fit_spline_plot)
figS6 <- ggarrange(plotlist = spline_plots[c("LGALS9", "DAPP1", "CCL13", "SKAP2")], labels = "AUTO")

# svg("../results/figS6_non_linearity.svg", width = 7, height = 7)
figS6
```

![](07_survival_analysis_files/figure-commonmark/unnamed-chunk-20-1.png)

``` r
# dev.off()
```

## Assessment of proportional hazards assumption (Fig. S11)

``` r
plot_schoenfeld_gg <- function(model, xlab = "Time (KM scale)",
                               ylab = "Scaled Schoenfeld residuals", title = NULL) {
    pval <- model$ph$table[1, "p"]

    ggplot(model$df, aes(time, resid)) +
        geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
        geom_point(alpha = 0.4, size = 1) +
        geom_smooth(method = "loess", se = TRUE, span = 0.9, linewidth = 0.8) +
        labs(x = xlab, y = ylab, title = title,
             subtitle = paste0("Proportional hazards test: p = ", formatC(pval, digits = 3, format = "f"))) +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 10))
}

schoenfeld_plots <- purrr::imap(schoenfeld_models, ~ plot_schoenfeld_gg(.x, title = .y))

figS11 <- ggarrange(plotlist = schoenfeld_plots[c("LGALS9", "DAPP1", "CCL13", "SKAP2")], labels = "AUTO")
```

    `geom_smooth()` using formula = 'y ~ x'
    `geom_smooth()` using formula = 'y ~ x'
    `geom_smooth()` using formula = 'y ~ x'
    `geom_smooth()` using formula = 'y ~ x'

``` r
# svg("../results/figS11_schoenfeld.svg", width = 7, height = 7)
figS11
```

![](07_survival_analysis_files/figure-commonmark/unnamed-chunk-21-1.png)

``` r
# dev.off()
```

## Forest plot: hazard ratios by quartile across outcomes (Fig. S5)

``` r
forest_data <- purrr::imap_dfr(categorical_models, function(models_by_outcome, marker) {
    purrr::imap_dfr(
        list(
            "Any SNAE" = models_by_outcome$group$crude,
            "Non-AIDS-defining malignancy" = models_by_outcome$nadm$crude,
            "Major cardiovascular event" = models_by_outcome$mace$crude,
            "Mortality" = models_by_outcome$death_all$crude
        ),
        function(model, outcome) {
            tibble(
                biomarker = marker, outcome = outcome,
                n = model$n, events = model$nevent,
                broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
            )
        }
    )
}) |>
    filter(str_detect(term, "^quantile")) |>
    mutate(
        quartile = recode(term, "quantile2" = "Q2 vs Q1", "quantile3" = "Q3 vs Q1", "quantile4" = "Q4 vs Q1"),
        quartile = factor(quartile, levels = c("Q4 vs Q1", "Q3 vs Q1", "Q2 vs Q1")),
        outcome = factor(outcome, levels = c("Any SNAE", "Non-AIDS-defining malignancy",
                                             "Major cardiovascular event", "Mortality")),
        biomarker = factor(biomarker, levels = markers),
        estimable = is.finite(estimate) & is.finite(conf.low) & is.finite(conf.high) & conf.low > 0
    )

forest_data |> arrange(biomarker, outcome, quartile)
```

    # A tibble: 48 × 13
       biomarker outcome         n events term  estimate std.error statistic p.value
       <fct>     <fct>       <int>  <dbl> <chr>    <dbl>     <dbl>     <dbl>   <dbl>
     1 LGALS9    Any SNAE      136     67 quan…   2.37e0     0.368   2.35     0.0188
     2 LGALS9    Any SNAE      136     67 quan…   2.04e0     0.380   1.88     0.0604
     3 LGALS9    Any SNAE      136     67 quan…   2.14e0     0.375   2.03     0.0425
     4 LGALS9    Non-AIDS-d…   136     45 quan…   1.32e0     0.420   0.652    0.515 
     5 LGALS9    Non-AIDS-d…   136     45 quan…   1.35e0     0.422   0.715    0.474 
     6 LGALS9    Non-AIDS-d…   136     45 quan…   1.31e0     0.421   0.636    0.525 
     7 LGALS9    Major card…   136     22 quan…   4.32e8  5874.      0.00339  0.997 
     8 LGALS9    Major card…   136     22 quan…   2.82e8  5874.      0.00331  0.997 
     9 LGALS9    Major card…   136     22 quan…   3.43e8  5874.      0.00335  0.997 
    10 LGALS9    Mortality     136     15 quan…   3.42e8  7153.      0.00275  0.998 
    # ℹ 38 more rows
    # ℹ 4 more variables: conf.low <dbl>, conf.high <dbl>, quartile <fct>,
    #   estimable <lgl>

``` r
forest_plot_data <- forest_data |>
    filter(estimable)

panel_labels <- c(
    "Any SNAE" = "A   Any SNAE",
    "Non-AIDS-defining malignancy" = "B   Non-AIDS-defining malignancy",
    "Major cardiovascular event" = "C   Major cardiovascular event",
    "Mortality" = "D   Mortality"
)

forest_plot <- ggplot(
    forest_plot_data,
    aes(x = quartile, y = estimate, ymin = conf.low, ymax = conf.high,
        colour = biomarker, shape = biomarker, group = biomarker)
) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.5) +
    geom_pointrange(position = position_dodge(width = 0.65), linewidth = 0.55,
                    key_glyph = ggplot2::draw_key_path) +
    scale_y_log10() +
    scale_colour_manual(
        values = c("LGALS9" = "#F8766D", "CCL13" = "#7CAE00", "DAPP1" = "#00BFC4", "SKAP2" = "#C77CFF"),
        drop = FALSE
    ) +
    scale_shape_manual(
        values = c("LGALS9" = 16, "CCL13" = 17, "DAPP1" = 15, "SKAP2" = 18),
        drop = FALSE
    ) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE), shape = "none") +
    facet_wrap(~ outcome, ncol = 2, labeller = as_labeller(panel_labels)) +
    coord_flip() +
    labs(x = NULL, y = "Hazard ratio with 95% confidence interval", colour = NULL, shape = NULL) +
    theme_bw(base_size = 10) +
    theme(
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"),
        legend.position = "top",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.box.just = "center",
        legend.box.background = element_rect(colour = "black", fill = "white", linewidth = 0.5),
        legend.key.width = grid::unit(1, "cm"),
        legend.key.height = grid::unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0, size = 11),
        panel.spacing = grid::unit(1, "lines")
    )

# svg("../results/figS5_forest_plot.svg", width = 7, height = 7)
forest_plot
```

![](07_survival_analysis_files/figure-commonmark/unnamed-chunk-23-1.png)

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
     package       * version date (UTC) lib source
     abind           1.4-8   2024-09-12 [1] CRAN (R 4.6.0)
     backports       1.5.1   2026-04-03 [1] CRAN (R 4.6.0)
     base64enc       0.1-6   2026-02-02 [1] CRAN (R 4.6.0)
     bit             4.6.0   2025-03-06 [1] CRAN (R 4.6.0)
     bit64           4.8.2   2026-05-19 [1] CRAN (R 4.6.0)
     broom         * 1.0.13  2026-05-14 [1] CRAN (R 4.6.0)
     broom.helpers   1.22.0  2025-09-17 [1] CRAN (R 4.6.0)
     car             3.1-5   2026-02-03 [1] CRAN (R 4.6.0)
     carData         3.0-6   2026-01-30 [1] CRAN (R 4.6.0)
     cards           0.7.1   2025-12-02 [1] CRAN (R 4.6.0)
     cli             3.6.6   2026-04-09 [1] CRAN (R 4.6.0)
     commonmark      2.0.0   2025-07-07 [1] CRAN (R 4.6.0)
     cowplot         1.2.0   2025-07-07 [1] CRAN (R 4.6.0)
     crayon          1.5.3   2024-06-20 [1] CRAN (R 4.6.0)
     digest          0.6.39  2025-11-19 [1] CRAN (R 4.6.0)
     dplyr         * 1.2.1   2026-04-03 [1] CRAN (R 4.6.0)
     evaluate        1.0.5   2025-08-27 [1] CRAN (R 4.6.0)
     farver          2.1.2   2024-05-13 [1] CRAN (R 4.6.0)
     fastmap         1.2.0   2024-05-15 [1] CRAN (R 4.6.0)
     forcats       * 1.0.1   2025-09-25 [1] CRAN (R 4.6.0)
     Formula         1.2-5   2023-02-24 [1] CRAN (R 4.6.0)
     fs              2.1.0   2026-04-18 [1] CRAN (R 4.6.0)
     generics        0.1.4   2025-05-09 [1] CRAN (R 4.6.0)
     ggplot2       * 4.0.3   2026-04-22 [1] CRAN (R 4.6.0)
     ggpubr        * 0.6.3   2026-02-24 [1] CRAN (R 4.6.0)
     ggsignif        0.6.4   2022-10-13 [1] CRAN (R 4.6.0)
     ggsurvfit     * 1.2.0   2025-09-13 [1] CRAN (R 4.6.0)
     glue            1.8.1   2026-04-17 [1] CRAN (R 4.6.0)
     gt              1.3.0   2026-01-22 [1] CRAN (R 4.6.0)
     gtable          0.3.6   2024-10-25 [1] CRAN (R 4.6.0)
     gtsummary     * 2.5.0   2025-12-05 [1] CRAN (R 4.6.0)
     haven           2.5.5   2025-05-30 [1] CRAN (R 4.6.0)
     hms             1.1.4   2025-10-17 [1] CRAN (R 4.6.0)
     htmltools       0.5.9   2025-12-04 [1] CRAN (R 4.6.0)
     jsonlite        2.0.0   2025-03-27 [1] CRAN (R 4.6.0)
     knitr           1.51    2025-12-20 [1] CRAN (R 4.6.0)
     labeling        0.4.3   2023-08-29 [1] CRAN (R 4.6.0)
     labelled        2.16.0  2025-10-22 [1] CRAN (R 4.6.0)
     lattice         0.22-9  2026-02-09 [1] CRAN (R 4.6.0)
     lifecycle       1.0.5   2026-01-08 [1] CRAN (R 4.6.0)
     litedown        0.9     2025-12-18 [1] CRAN (R 4.6.0)
     lubridate     * 1.9.5   2026-02-04 [1] CRAN (R 4.6.0)
     magrittr        2.0.5   2026-04-04 [1] CRAN (R 4.6.0)
     markdown        2.0     2025-03-23 [1] CRAN (R 4.6.0)
     Matrix          1.7-5   2026-03-21 [4] CRAN (R 4.5.3)
     mgcv            1.9-4   2025-11-07 [1] CRAN (R 4.6.0)
     nlme            3.1-169 2026-03-27 [1] CRAN (R 4.6.0)
     otel            0.2.0   2025-08-29 [1] CRAN (R 4.6.0)
     patchwork     * 1.3.2   2025-08-25 [1] CRAN (R 4.6.0)
     pillar          1.11.1  2025-09-17 [1] CRAN (R 4.6.0)
     pkgconfig       2.0.3   2019-09-22 [1] CRAN (R 4.6.0)
     purrr         * 1.2.2   2026-04-10 [1] CRAN (R 4.6.0)
     R6              2.6.1   2025-02-15 [1] CRAN (R 4.6.0)
     RColorBrewer    1.1-3   2022-04-03 [1] CRAN (R 4.6.0)
     readr         * 2.2.0   2026-02-19 [1] CRAN (R 4.6.0)
     rlang           1.2.0   2026-04-06 [1] CRAN (R 4.6.0)
     rmarkdown       2.31    2026-03-26 [1] CRAN (R 4.6.0)
     rstatix         0.7.3   2025-10-18 [1] CRAN (R 4.6.0)
     rstudioapi      0.18.0  2026-01-16 [1] CRAN (R 4.6.0)
     S7              0.2.2   2026-04-22 [1] CRAN (R 4.6.0)
     sass            0.4.10  2025-04-11 [1] CRAN (R 4.6.0)
     scales          1.4.0   2025-04-24 [1] CRAN (R 4.6.0)
     sessioninfo     1.2.3   2025-02-05 [1] CRAN (R 4.6.0)
     stringi         1.8.7   2025-03-27 [1] CRAN (R 4.6.0)
     stringr       * 1.6.0   2025-11-04 [1] CRAN (R 4.6.0)
     survival      * 3.8-6   2026-01-16 [4] CRAN (R 4.5.2)
     tibble        * 3.3.1   2026-01-11 [1] CRAN (R 4.6.0)
     tidyr         * 1.3.2   2025-12-19 [1] CRAN (R 4.6.0)
     tidyselect      1.2.1   2024-03-11 [1] CRAN (R 4.6.0)
     tidyverse     * 2.0.0   2023-02-22 [1] CRAN (R 4.6.0)
     timechange      0.4.0   2026-01-29 [1] CRAN (R 4.6.0)
     tzdb            0.5.0   2025-03-15 [1] CRAN (R 4.6.0)
     utf8            1.2.6   2025-06-08 [1] CRAN (R 4.6.0)
     vctrs           0.7.3   2026-04-11 [1] CRAN (R 4.6.0)
     vroom           1.7.1   2026-03-31 [1] CRAN (R 4.6.0)
     withr           3.0.2   2024-10-28 [1] CRAN (R 4.6.0)
     xfun            0.60    2026-07-09 [1] CRAN (R 4.6.0)
     xml2            1.5.2   2026-01-17 [1] CRAN (R 4.6.0)
     yaml            2.3.12  2025-12-10 [1] CRAN (R 4.6.0)

     [1] /opt/Rlibs
     [2] /usr/local/lib/R/site-library
     [3] /usr/lib/R/site-library
     [4] /usr/lib/R/library
     * ── Packages attached to the search path.

    ──────────────────────────────────────────────────────────────────────────────
