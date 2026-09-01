# 3. Table 1


## Purpose

This script builds Table 1 using gtsummary. Participant characteristics
are stratified by case/control status, for the full cohort and for three
subgroups (the Olink analytic set, the RNA-seq analytic set, and
participants excluded after Olink QC).

## Reproducibility

This script runs directly from the public `metadata.tsv`. Four variables
with elevated re-identification risk are present in `metadata.tsv` as
`NA` placeholder columns: `age`, `mode` (HIV acquisition risk factor),
`paisorigen` (country of origin), and `artini_period` (year of ART
initiation). The table structure, sample sizes, and every other
statistic reproduce exactly, but these four columns display as entirely
missing when run from the public repository, unlike the published Table
1, which used the corresponding values from the restricted CoRIS
extract.

``` r
library(gtsummary)
library(tidyverse)
```

## Helper functions

The four tables below (full cohort, Olink subset, RNA-seq subset,
Olink-excluded subset) differ only in *which participants* are included.

``` r
# Reads metadata.tsv and optionally restricts to a set of `sid`s.
prepare_table1_data <- function(filter_expr = TRUE) {
    mdata <- read_tsv("../data/metadata.tsv") |>
        mutate(across(c(olink_included, rnaseq_included), as.logical)) |> 
        filter({{ filter_expr }})
    
    mdata |>
        select(group, nadm, mace, death_all,
               age, sex, mode, paisorigen, edulvl,
               smoke, hta, art_reg, cvmax, nadircd4,
               CD4_0, CD8_0, CD4CD8_0, CD4_24m, CD8_24m, CD4CD8_24,
               months_art_sampling, artini_period) |>
        mutate(across(
            c(group, nadm, mace, death_all, sex,
              edulvl, smoke, hta, art_reg),
            as.factor
        ))
}
```

``` r
# Builds the gtsummary Table 1 from a prepared dataset (output of
# prepare_table1_data()).
build_gtsummary_table1 <- function(mdata) {
    n_group <- table(mdata$group)
    n_display <- setNames(sprintf("%d (%.1f%%)", n_group, 100 * prop.table(n_group)), names(n_group))
    
    mdata |>
        tbl_summary(
            by = group,
            label = list(
                nadm ~ "Non-AIDS defining malignancy",
                mace ~ "Major cardiovascular events",
                death_all ~ "Non-accidental death",
                age ~ "Age",
                sex ~ "Sex at birth",
                mode ~ "Risk factor for HIV acquisition",
                paisorigen ~ "Country of origin",
                edulvl ~ "Education attainment",
                smoke ~ "Smoker",
                hta ~ "Hypertension",
                art_reg ~ "ART class",
                cvmax ~ "Maximum viral load (copies/mL)",
                nadircd4 ~ "CD4 nadir (cells/µL)",
                CD4_0 ~ "CD4 counts at ART initiation (cells/µL)",
                CD8_0 ~ "CD8 counts at ART initiation (cells/µL)",
                CD4CD8_0 ~ "CD4/CD8 ratio at ART initiation",
                CD4_24m ~ "CD4 counts at Year 2 of ART (cells/µL)",
                CD8_24m ~ "CD8 counts at Year 2 of ART (cells/µL)",
                CD4CD8_24 ~ "CD4/CD8 ratio at Year 2 of ART",
                months_art_sampling ~ "Months from ART initiation to sampling",
                artini_period ~ "Year of ART initiation"
            ),
            statistic = all_continuous() ~ "{mean} ({sd})",
            digits = list(
                c(all_continuous(), -cvmax) ~ 1,
                cvmax ~ 0,
                all_categorical() ~ c(0, 1)
            ),
            sort = list(
                all_categorical() ~ "frequency",
                edulvl ~ "alphanumeric",
                artini_period ~ "alphanumeric"
            ),
            missing = "no",
            percent = mdata,
            type = list(
                age ~ "continuous",
                mode ~ "continuous",
                paisorigen ~ "continuous",
                artini_period ~ "continuous"
            )
        ) |>
        add_overall(last = TRUE, col_label = "**Total**") |>
        add_p(
            test = list(
                all_continuous() ~ "t.test",
                all_categorical() ~ "fisher.test"),
            pvalue_fun = label_style_pvalue(digits = 3),
            include = -c(age, mode, paisorigen, artini_period)
        ) |>
        modify_header(
            all_stat_cols(stat_0 = FALSE) ~ "**{level}**"
        ) |>
        modify_table_body(
            ~ add_row(
                .x,
                variable = "N",
                row_type = "label",
                label = "N",
                stat_1 = n_display["Event"],
                stat_2 = n_display["Control"],
                stat_0 = sprintf("%d (100.0%%)", sum(n_group)),
                .before = 1
            )
        )
}
```

## Table 1 (full cohort)

``` r
table_gts <- prepare_table1_data() |> build_gtsummary_table1()
```

    Rows: 187 Columns: 27
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (6): sid, group, sex, edulvl, art_reg, sid_pair
    dbl (17): nadm, mace, death_all, smoke, hta, cvmax, nadircd4, CD4_0, CD8_0, ...
    lgl  (4): age, mode, paisorigen, artini_period

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
table_gts
```

<div id="jygjmaoxfz" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#jygjmaoxfz table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#jygjmaoxfz thead, #jygjmaoxfz tbody, #jygjmaoxfz tfoot, #jygjmaoxfz tr, #jygjmaoxfz td, #jygjmaoxfz th {
  border-style: none;
}
&#10;#jygjmaoxfz p {
  margin: 0;
  padding: 0;
}
&#10;#jygjmaoxfz .gt_table {
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
&#10;#jygjmaoxfz .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#jygjmaoxfz .gt_title {
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
&#10;#jygjmaoxfz .gt_subtitle {
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
&#10;#jygjmaoxfz .gt_heading {
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
&#10;#jygjmaoxfz .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#jygjmaoxfz .gt_col_headings {
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
&#10;#jygjmaoxfz .gt_col_heading {
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
&#10;#jygjmaoxfz .gt_column_spanner_outer {
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
&#10;#jygjmaoxfz .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#jygjmaoxfz .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#jygjmaoxfz .gt_column_spanner {
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
&#10;#jygjmaoxfz .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#jygjmaoxfz .gt_group_heading {
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
&#10;#jygjmaoxfz .gt_empty_group_heading {
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
&#10;#jygjmaoxfz .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#jygjmaoxfz .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#jygjmaoxfz .gt_row {
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
&#10;#jygjmaoxfz .gt_stub {
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
&#10;#jygjmaoxfz .gt_stub_row_group {
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
&#10;#jygjmaoxfz .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#jygjmaoxfz .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#jygjmaoxfz .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jygjmaoxfz .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#jygjmaoxfz .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#jygjmaoxfz .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#jygjmaoxfz .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jygjmaoxfz .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#jygjmaoxfz .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#jygjmaoxfz .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#jygjmaoxfz .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#jygjmaoxfz .gt_footnotes {
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
&#10;#jygjmaoxfz .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jygjmaoxfz .gt_sourcenotes {
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
&#10;#jygjmaoxfz .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jygjmaoxfz .gt_left {
  text-align: left;
}
&#10;#jygjmaoxfz .gt_center {
  text-align: center;
}
&#10;#jygjmaoxfz .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#jygjmaoxfz .gt_font_normal {
  font-weight: normal;
}
&#10;#jygjmaoxfz .gt_font_bold {
  font-weight: bold;
}
&#10;#jygjmaoxfz .gt_font_italic {
  font-style: italic;
}
&#10;#jygjmaoxfz .gt_super {
  font-size: 65%;
}
&#10;#jygjmaoxfz .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#jygjmaoxfz .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#jygjmaoxfz .gt_indent_1 {
  text-indent: 5px;
}
&#10;#jygjmaoxfz .gt_indent_2 {
  text-indent: 10px;
}
&#10;#jygjmaoxfz .gt_indent_3 {
  text-indent: 15px;
}
&#10;#jygjmaoxfz .gt_indent_4 {
  text-indent: 20px;
}
&#10;#jygjmaoxfz .gt_indent_5 {
  text-indent: 25px;
}
&#10;#jygjmaoxfz .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#jygjmaoxfz div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>

<table class="gt_table" data-quarto-postprocess="true"
data-quarto-disable-processing="false" data-quarto-bootstrap="false">
<colgroup>
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
</colgroup>
<thead>
<tr class="gt_col_headings">
<th id="label" class="gt_col_heading gt_columns_bottom_border gt_left"
data-quarto-table-cell-role="th"
scope="col"><strong>Characteristic</strong></th>
<th id="stat_1"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th"
scope="col"><strong>Control</strong><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="stat_2"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th" scope="col"><strong>Event</strong><span
class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="stat_0"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th" scope="col"><strong>Total</strong><span
class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="p.value"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th"
scope="col"><strong>p-value</strong><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span></th>
</tr>
</thead>
<tbody class="gt_table_body">
<tr>
<td class="gt_row gt_left" headers="label">N</td>
<td class="gt_row gt_center" headers="stat_1">92 (49.2%)</td>
<td class="gt_row gt_center" headers="stat_2">95 (50.8%)</td>
<td class="gt_row gt_center" headers="stat_0">187 (100.0%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Non-AIDS defining
malignancy</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&lt;0.001</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">95 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">28 (30.4%)</td>
<td class="gt_row gt_center" headers="stat_0">123 (65.8%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">64 (69.6%)</td>
<td class="gt_row gt_center" headers="stat_0">64 (34.2%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Major cardiovascular
events</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&lt;0.001</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">95 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">64 (69.6%)</td>
<td class="gt_row gt_center" headers="stat_0">159 (85.0%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">28 (30.4%)</td>
<td class="gt_row gt_center" headers="stat_0">28 (15.0%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Non-accidental death</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&lt;0.001</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">95 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">69 (75.0%)</td>
<td class="gt_row gt_center" headers="stat_0">164 (87.7%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">23 (25.0%)</td>
<td class="gt_row gt_center" headers="stat_0">23 (12.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Age</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Sex at birth</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.731</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Male</td>
<td class="gt_row gt_center" headers="stat_1">74 (77.9%)</td>
<td class="gt_row gt_center" headers="stat_2">69 (75.0%)</td>
<td class="gt_row gt_center" headers="stat_0">143 (76.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Female</td>
<td class="gt_row gt_center" headers="stat_1">21 (22.1%)</td>
<td class="gt_row gt_center" headers="stat_2">23 (25.0%)</td>
<td class="gt_row gt_center" headers="stat_0">44 (23.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Risk factor for HIV
acquisition</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Country of origin</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Education attainment</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.497</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    High school</td>
<td class="gt_row gt_center" headers="stat_1">29 (30.5%)</td>
<td class="gt_row gt_center" headers="stat_2">23 (25.0%)</td>
<td class="gt_row gt_center" headers="stat_0">52 (27.8%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    No studies</td>
<td class="gt_row gt_center" headers="stat_1">6 (6.3%)</td>
<td class="gt_row gt_center" headers="stat_2">4 (4.3%)</td>
<td class="gt_row gt_center" headers="stat_0">10 (5.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Other</td>
<td class="gt_row gt_center" headers="stat_1">1 (1.1%)</td>
<td class="gt_row gt_center" headers="stat_2">3 (3.3%)</td>
<td class="gt_row gt_center" headers="stat_0">4 (2.1%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Primary school</td>
<td class="gt_row gt_center" headers="stat_1">12 (12.6%)</td>
<td class="gt_row gt_center" headers="stat_2">17 (18.5%)</td>
<td class="gt_row gt_center" headers="stat_0">29 (15.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Secondary school</td>
<td class="gt_row gt_center" headers="stat_1">26 (27.4%)</td>
<td class="gt_row gt_center" headers="stat_2">19 (20.7%)</td>
<td class="gt_row gt_center" headers="stat_0">45 (24.1%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    University</td>
<td class="gt_row gt_center" headers="stat_1">21 (22.1%)</td>
<td class="gt_row gt_center" headers="stat_2">26 (28.3%)</td>
<td class="gt_row gt_center" headers="stat_0">47 (25.1%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Smoker</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.609</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">35 (36.8%)</td>
<td class="gt_row gt_center" headers="stat_2">34 (37.0%)</td>
<td class="gt_row gt_center" headers="stat_0">69 (36.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">39 (41.1%)</td>
<td class="gt_row gt_center" headers="stat_2">30 (32.6%)</td>
<td class="gt_row gt_center" headers="stat_0">69 (36.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Hypertension</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.129</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">24 (25.3%)</td>
<td class="gt_row gt_center" headers="stat_2">28 (30.4%)</td>
<td class="gt_row gt_center" headers="stat_0">52 (27.8%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">14 (14.7%)</td>
<td class="gt_row gt_center" headers="stat_2">7 (7.6%)</td>
<td class="gt_row gt_center" headers="stat_0">21 (11.2%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">ART class</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.994</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 NNRTI</td>
<td class="gt_row gt_center" headers="stat_1">42 (44.2%)</td>
<td class="gt_row gt_center" headers="stat_2">40 (43.5%)</td>
<td class="gt_row gt_center" headers="stat_0">82 (43.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 PI</td>
<td class="gt_row gt_center" headers="stat_1">37 (38.9%)</td>
<td class="gt_row gt_center" headers="stat_2">36 (39.1%)</td>
<td class="gt_row gt_center" headers="stat_0">73 (39.0%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 INSTI</td>
<td class="gt_row gt_center" headers="stat_1">12 (12.6%)</td>
<td class="gt_row gt_center" headers="stat_2">11 (12.0%)</td>
<td class="gt_row gt_center" headers="stat_0">23 (12.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Other</td>
<td class="gt_row gt_center" headers="stat_1">4 (4.2%)</td>
<td class="gt_row gt_center" headers="stat_2">5 (5.4%)</td>
<td class="gt_row gt_center" headers="stat_0">9 (4.8%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Maximum viral load
(copies/mL)</td>
<td class="gt_row gt_center" headers="stat_1">482,101 (1,164,967)</td>
<td class="gt_row gt_center" headers="stat_2">413,461 (1,063,940)</td>
<td class="gt_row gt_center" headers="stat_0">448,332 (1,113,942)</td>
<td class="gt_row gt_center" headers="p.value">0.674</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 nadir (cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">240.3 (178.4)</td>
<td class="gt_row gt_center" headers="stat_2">197.7 (150.8)</td>
<td class="gt_row gt_center" headers="stat_0">219.0 (166.1)</td>
<td class="gt_row gt_center" headers="p.value">0.082</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 counts at ART initiation
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">289.8 (211.9)</td>
<td class="gt_row gt_center" headers="stat_2">240.5 (179.3)</td>
<td class="gt_row gt_center" headers="stat_0">265.3 (197.3)</td>
<td class="gt_row gt_center" headers="p.value">0.100</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD8 counts at ART initiation
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">997.4 (564.0)</td>
<td class="gt_row gt_center" headers="stat_2">968.0 (530.0)</td>
<td class="gt_row gt_center" headers="stat_0">982.8 (546.0)</td>
<td class="gt_row gt_center" headers="p.value">0.724</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4/CD8 ratio at ART
initiation</td>
<td class="gt_row gt_center" headers="stat_1">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="stat_2">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="stat_0">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="p.value">0.248</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 counts at Year 2 of ART
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">556.5 (270.5)</td>
<td class="gt_row gt_center" headers="stat_2">487.6 (274.5)</td>
<td class="gt_row gt_center" headers="stat_0">521.6 (274.0)</td>
<td class="gt_row gt_center" headers="p.value">0.090</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD8 counts at Year 2 of ART
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">1,004.5 (502.3)</td>
<td class="gt_row gt_center" headers="stat_2">998.5 (575.9)</td>
<td class="gt_row gt_center" headers="stat_0">1,001.5 (539.3)</td>
<td class="gt_row gt_center" headers="p.value">0.940</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4/CD8 ratio at Year 2 of
ART</td>
<td class="gt_row gt_center" headers="stat_1">0.6 (0.4)</td>
<td class="gt_row gt_center" headers="stat_2">0.6 (0.4)</td>
<td class="gt_row gt_center" headers="stat_0">0.6 (0.4)</td>
<td class="gt_row gt_center" headers="p.value">0.217</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Months from ART initiation to
sampling</td>
<td class="gt_row gt_center" headers="stat_1">24.2 (5.7)</td>
<td class="gt_row gt_center" headers="stat_2">24.2 (7.3)</td>
<td class="gt_row gt_center" headers="stat_0">24.2 (6.5)</td>
<td class="gt_row gt_center" headers="p.value">0.985</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Year of ART initiation</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
</tbody><tfoot>
<tr class="gt_footnotes">
<td colspan="5" class="gt_footnote"><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span>
n (%); Mean (SD)</td>
</tr>
<tr class="gt_footnotes">
<td colspan="5" class="gt_footnote"><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span>
Fisher’s exact test; Welch Two Sample t-test</td>
</tr>
</tfoot>
&#10;</table>

</div>

``` r
write_tsv(table_gts |> as_tibble(), "../results/table1.tsv", na = "")
```

## Table 1 for Olink

``` r
table_gts_olink <- prepare_table1_data(olink_included) |> build_gtsummary_table1()
```

    Rows: 187 Columns: 27
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (6): sid, group, sex, edulvl, art_reg, sid_pair
    dbl (17): nadm, mace, death_all, smoke, hta, cvmax, nadircd4, CD4_0, CD8_0, ...
    lgl  (4): age, mode, paisorigen, artini_period

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
table_gts_olink
```

<div id="owfjkwhnqk" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#owfjkwhnqk table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#owfjkwhnqk thead, #owfjkwhnqk tbody, #owfjkwhnqk tfoot, #owfjkwhnqk tr, #owfjkwhnqk td, #owfjkwhnqk th {
  border-style: none;
}
&#10;#owfjkwhnqk p {
  margin: 0;
  padding: 0;
}
&#10;#owfjkwhnqk .gt_table {
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
&#10;#owfjkwhnqk .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#owfjkwhnqk .gt_title {
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
&#10;#owfjkwhnqk .gt_subtitle {
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
&#10;#owfjkwhnqk .gt_heading {
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
&#10;#owfjkwhnqk .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#owfjkwhnqk .gt_col_headings {
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
&#10;#owfjkwhnqk .gt_col_heading {
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
&#10;#owfjkwhnqk .gt_column_spanner_outer {
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
&#10;#owfjkwhnqk .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#owfjkwhnqk .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#owfjkwhnqk .gt_column_spanner {
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
&#10;#owfjkwhnqk .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#owfjkwhnqk .gt_group_heading {
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
&#10;#owfjkwhnqk .gt_empty_group_heading {
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
&#10;#owfjkwhnqk .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#owfjkwhnqk .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#owfjkwhnqk .gt_row {
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
&#10;#owfjkwhnqk .gt_stub {
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
&#10;#owfjkwhnqk .gt_stub_row_group {
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
&#10;#owfjkwhnqk .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#owfjkwhnqk .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#owfjkwhnqk .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#owfjkwhnqk .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#owfjkwhnqk .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#owfjkwhnqk .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#owfjkwhnqk .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#owfjkwhnqk .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#owfjkwhnqk .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#owfjkwhnqk .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#owfjkwhnqk .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#owfjkwhnqk .gt_footnotes {
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
&#10;#owfjkwhnqk .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#owfjkwhnqk .gt_sourcenotes {
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
&#10;#owfjkwhnqk .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#owfjkwhnqk .gt_left {
  text-align: left;
}
&#10;#owfjkwhnqk .gt_center {
  text-align: center;
}
&#10;#owfjkwhnqk .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#owfjkwhnqk .gt_font_normal {
  font-weight: normal;
}
&#10;#owfjkwhnqk .gt_font_bold {
  font-weight: bold;
}
&#10;#owfjkwhnqk .gt_font_italic {
  font-style: italic;
}
&#10;#owfjkwhnqk .gt_super {
  font-size: 65%;
}
&#10;#owfjkwhnqk .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#owfjkwhnqk .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#owfjkwhnqk .gt_indent_1 {
  text-indent: 5px;
}
&#10;#owfjkwhnqk .gt_indent_2 {
  text-indent: 10px;
}
&#10;#owfjkwhnqk .gt_indent_3 {
  text-indent: 15px;
}
&#10;#owfjkwhnqk .gt_indent_4 {
  text-indent: 20px;
}
&#10;#owfjkwhnqk .gt_indent_5 {
  text-indent: 25px;
}
&#10;#owfjkwhnqk .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#owfjkwhnqk div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>

<table class="gt_table" data-quarto-postprocess="true"
data-quarto-disable-processing="false" data-quarto-bootstrap="false">
<colgroup>
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
</colgroup>
<thead>
<tr class="gt_col_headings">
<th id="label" class="gt_col_heading gt_columns_bottom_border gt_left"
data-quarto-table-cell-role="th"
scope="col"><strong>Characteristic</strong></th>
<th id="stat_1"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th"
scope="col"><strong>Control</strong><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="stat_2"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th" scope="col"><strong>Event</strong><span
class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="stat_0"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th" scope="col"><strong>Total</strong><span
class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="p.value"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th"
scope="col"><strong>p-value</strong><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span></th>
</tr>
</thead>
<tbody class="gt_table_body">
<tr>
<td class="gt_row gt_left" headers="label">N</td>
<td class="gt_row gt_center" headers="stat_1">69 (50.0%)</td>
<td class="gt_row gt_center" headers="stat_2">69 (50.0%)</td>
<td class="gt_row gt_center" headers="stat_0">138 (100.0%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Non-AIDS defining
malignancy</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&lt;0.001</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">69 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">22 (31.9%)</td>
<td class="gt_row gt_center" headers="stat_0">91 (65.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">47 (68.1%)</td>
<td class="gt_row gt_center" headers="stat_0">47 (34.1%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Major cardiovascular
events</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&lt;0.001</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">69 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">47 (68.1%)</td>
<td class="gt_row gt_center" headers="stat_0">116 (84.1%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">22 (31.9%)</td>
<td class="gt_row gt_center" headers="stat_0">22 (15.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Non-accidental death</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&lt;0.001</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">69 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">54 (78.3%)</td>
<td class="gt_row gt_center" headers="stat_0">123 (89.1%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">15 (21.7%)</td>
<td class="gt_row gt_center" headers="stat_0">15 (10.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Age</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Sex at birth</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.844</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Male</td>
<td class="gt_row gt_center" headers="stat_1">53 (76.8%)</td>
<td class="gt_row gt_center" headers="stat_2">51 (73.9%)</td>
<td class="gt_row gt_center" headers="stat_0">104 (75.4%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Female</td>
<td class="gt_row gt_center" headers="stat_1">16 (23.2%)</td>
<td class="gt_row gt_center" headers="stat_2">18 (26.1%)</td>
<td class="gt_row gt_center" headers="stat_0">34 (24.6%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Risk factor for HIV
acquisition</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Country of origin</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Education attainment</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.886</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    High school</td>
<td class="gt_row gt_center" headers="stat_1">21 (30.4%)</td>
<td class="gt_row gt_center" headers="stat_2">18 (26.1%)</td>
<td class="gt_row gt_center" headers="stat_0">39 (28.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    No studies</td>
<td class="gt_row gt_center" headers="stat_1">4 (5.8%)</td>
<td class="gt_row gt_center" headers="stat_2">3 (4.3%)</td>
<td class="gt_row gt_center" headers="stat_0">7 (5.1%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Other</td>
<td class="gt_row gt_center" headers="stat_1">1 (1.4%)</td>
<td class="gt_row gt_center" headers="stat_2">3 (4.3%)</td>
<td class="gt_row gt_center" headers="stat_0">4 (2.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Primary school</td>
<td class="gt_row gt_center" headers="stat_1">10 (14.5%)</td>
<td class="gt_row gt_center" headers="stat_2">13 (18.8%)</td>
<td class="gt_row gt_center" headers="stat_0">23 (16.7%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Secondary school</td>
<td class="gt_row gt_center" headers="stat_1">17 (24.6%)</td>
<td class="gt_row gt_center" headers="stat_2">15 (21.7%)</td>
<td class="gt_row gt_center" headers="stat_0">32 (23.2%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    University</td>
<td class="gt_row gt_center" headers="stat_1">16 (23.2%)</td>
<td class="gt_row gt_center" headers="stat_2">17 (24.6%)</td>
<td class="gt_row gt_center" headers="stat_0">33 (23.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Smoker</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.557</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">28 (40.6%)</td>
<td class="gt_row gt_center" headers="stat_2">27 (39.1%)</td>
<td class="gt_row gt_center" headers="stat_0">55 (39.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">29 (42.0%)</td>
<td class="gt_row gt_center" headers="stat_2">21 (30.4%)</td>
<td class="gt_row gt_center" headers="stat_0">50 (36.2%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Hypertension</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.391</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">17 (24.6%)</td>
<td class="gt_row gt_center" headers="stat_2">21 (30.4%)</td>
<td class="gt_row gt_center" headers="stat_0">38 (27.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">10 (14.5%)</td>
<td class="gt_row gt_center" headers="stat_2">7 (10.1%)</td>
<td class="gt_row gt_center" headers="stat_0">17 (12.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">ART class</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&gt;0.999</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 PI</td>
<td class="gt_row gt_center" headers="stat_1">29 (42.0%)</td>
<td class="gt_row gt_center" headers="stat_2">29 (42.0%)</td>
<td class="gt_row gt_center" headers="stat_0">58 (42.0%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 NNRTI</td>
<td class="gt_row gt_center" headers="stat_1">29 (42.0%)</td>
<td class="gt_row gt_center" headers="stat_2">28 (40.6%)</td>
<td class="gt_row gt_center" headers="stat_0">57 (41.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 INSTI</td>
<td class="gt_row gt_center" headers="stat_1">8 (11.6%)</td>
<td class="gt_row gt_center" headers="stat_2">8 (11.6%)</td>
<td class="gt_row gt_center" headers="stat_0">16 (11.6%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Other</td>
<td class="gt_row gt_center" headers="stat_1">3 (4.3%)</td>
<td class="gt_row gt_center" headers="stat_2">4 (5.8%)</td>
<td class="gt_row gt_center" headers="stat_0">7 (5.1%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Maximum viral load
(copies/mL)</td>
<td class="gt_row gt_center" headers="stat_1">504,163 (1,294,304)</td>
<td class="gt_row gt_center" headers="stat_2">384,818 (1,038,571)</td>
<td class="gt_row gt_center" headers="stat_0">444,491 (1,170,668)</td>
<td class="gt_row gt_center" headers="p.value">0.551</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 nadir (cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">248.5 (187.3)</td>
<td class="gt_row gt_center" headers="stat_2">191.8 (156.3)</td>
<td class="gt_row gt_center" headers="stat_0">220.2 (174.2)</td>
<td class="gt_row gt_center" headers="p.value">0.056</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 counts at ART initiation
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">301.6 (223.0)</td>
<td class="gt_row gt_center" headers="stat_2">228.2 (178.0)</td>
<td class="gt_row gt_center" headers="stat_0">265.5 (204.7)</td>
<td class="gt_row gt_center" headers="p.value">0.040</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD8 counts at ART initiation
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">1,035.5 (566.4)</td>
<td class="gt_row gt_center" headers="stat_2">982.3 (561.3)</td>
<td class="gt_row gt_center" headers="stat_0">1,009.3 (562.3)</td>
<td class="gt_row gt_center" headers="p.value">0.592</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4/CD8 ratio at ART
initiation</td>
<td class="gt_row gt_center" headers="stat_1">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="stat_2">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="stat_0">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="p.value">0.341</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 counts at Year 2 of ART
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">562.9 (270.9)</td>
<td class="gt_row gt_center" headers="stat_2">466.2 (279.0)</td>
<td class="gt_row gt_center" headers="stat_0">514.2 (278.2)</td>
<td class="gt_row gt_center" headers="p.value">0.041</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD8 counts at Year 2 of ART
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">1,001.9 (472.4)</td>
<td class="gt_row gt_center" headers="stat_2">998.3 (627.3)</td>
<td class="gt_row gt_center" headers="stat_0">1,000.1 (553.8)</td>
<td class="gt_row gt_center" headers="p.value">0.969</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4/CD8 ratio at Year 2 of
ART</td>
<td class="gt_row gt_center" headers="stat_1">0.7 (0.4)</td>
<td class="gt_row gt_center" headers="stat_2">0.5 (0.4)</td>
<td class="gt_row gt_center" headers="stat_0">0.6 (0.4)</td>
<td class="gt_row gt_center" headers="p.value">0.106</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Months from ART initiation to
sampling</td>
<td class="gt_row gt_center" headers="stat_1">24.1 (5.4)</td>
<td class="gt_row gt_center" headers="stat_2">24.5 (5.8)</td>
<td class="gt_row gt_center" headers="stat_0">24.3 (5.6)</td>
<td class="gt_row gt_center" headers="p.value">0.663</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Year of ART initiation</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
</tbody><tfoot>
<tr class="gt_footnotes">
<td colspan="5" class="gt_footnote"><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span>
n (%); Mean (SD)</td>
</tr>
<tr class="gt_footnotes">
<td colspan="5" class="gt_footnote"><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span>
Fisher’s exact test; Welch Two Sample t-test</td>
</tr>
</tfoot>
&#10;</table>

</div>

``` r
write_tsv(table_gts_olink |> as_tibble(), "../results/tableS1_olink.tsv", na = "")
```

## Table 1 for RNA-seq

``` r
table_gts_rnaseq <- prepare_table1_data(rnaseq_included) |> build_gtsummary_table1()
```

    Rows: 187 Columns: 27
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (6): sid, group, sex, edulvl, art_reg, sid_pair
    dbl (17): nadm, mace, death_all, smoke, hta, cvmax, nadircd4, CD4_0, CD8_0, ...
    lgl  (4): age, mode, paisorigen, artini_period

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
table_gts_rnaseq
```

<div id="decljpwksl" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#decljpwksl table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#decljpwksl thead, #decljpwksl tbody, #decljpwksl tfoot, #decljpwksl tr, #decljpwksl td, #decljpwksl th {
  border-style: none;
}
&#10;#decljpwksl p {
  margin: 0;
  padding: 0;
}
&#10;#decljpwksl .gt_table {
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
&#10;#decljpwksl .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#decljpwksl .gt_title {
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
&#10;#decljpwksl .gt_subtitle {
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
&#10;#decljpwksl .gt_heading {
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
&#10;#decljpwksl .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#decljpwksl .gt_col_headings {
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
&#10;#decljpwksl .gt_col_heading {
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
&#10;#decljpwksl .gt_column_spanner_outer {
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
&#10;#decljpwksl .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#decljpwksl .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#decljpwksl .gt_column_spanner {
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
&#10;#decljpwksl .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#decljpwksl .gt_group_heading {
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
&#10;#decljpwksl .gt_empty_group_heading {
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
&#10;#decljpwksl .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#decljpwksl .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#decljpwksl .gt_row {
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
&#10;#decljpwksl .gt_stub {
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
&#10;#decljpwksl .gt_stub_row_group {
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
&#10;#decljpwksl .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#decljpwksl .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#decljpwksl .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#decljpwksl .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#decljpwksl .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#decljpwksl .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#decljpwksl .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#decljpwksl .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#decljpwksl .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#decljpwksl .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#decljpwksl .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#decljpwksl .gt_footnotes {
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
&#10;#decljpwksl .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#decljpwksl .gt_sourcenotes {
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
&#10;#decljpwksl .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#decljpwksl .gt_left {
  text-align: left;
}
&#10;#decljpwksl .gt_center {
  text-align: center;
}
&#10;#decljpwksl .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#decljpwksl .gt_font_normal {
  font-weight: normal;
}
&#10;#decljpwksl .gt_font_bold {
  font-weight: bold;
}
&#10;#decljpwksl .gt_font_italic {
  font-style: italic;
}
&#10;#decljpwksl .gt_super {
  font-size: 65%;
}
&#10;#decljpwksl .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#decljpwksl .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#decljpwksl .gt_indent_1 {
  text-indent: 5px;
}
&#10;#decljpwksl .gt_indent_2 {
  text-indent: 10px;
}
&#10;#decljpwksl .gt_indent_3 {
  text-indent: 15px;
}
&#10;#decljpwksl .gt_indent_4 {
  text-indent: 20px;
}
&#10;#decljpwksl .gt_indent_5 {
  text-indent: 25px;
}
&#10;#decljpwksl .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#decljpwksl div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>

<table class="gt_table" data-quarto-postprocess="true"
data-quarto-disable-processing="false" data-quarto-bootstrap="false">
<colgroup>
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
</colgroup>
<thead>
<tr class="gt_col_headings">
<th id="label" class="gt_col_heading gt_columns_bottom_border gt_left"
data-quarto-table-cell-role="th"
scope="col"><strong>Characteristic</strong></th>
<th id="stat_1"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th"
scope="col"><strong>Control</strong><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="stat_2"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th" scope="col"><strong>Event</strong><span
class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="stat_0"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th" scope="col"><strong>Total</strong><span
class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="p.value"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th"
scope="col"><strong>p-value</strong><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span></th>
</tr>
</thead>
<tbody class="gt_table_body">
<tr>
<td class="gt_row gt_left" headers="label">N</td>
<td class="gt_row gt_center" headers="stat_1">75 (48.4%)</td>
<td class="gt_row gt_center" headers="stat_2">80 (51.6%)</td>
<td class="gt_row gt_center" headers="stat_0">155 (100.0%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Non-AIDS defining
malignancy</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&lt;0.001</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">80 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">23 (30.7%)</td>
<td class="gt_row gt_center" headers="stat_0">103 (66.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">52 (69.3%)</td>
<td class="gt_row gt_center" headers="stat_0">52 (33.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Major cardiovascular
events</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&lt;0.001</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">80 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">52 (69.3%)</td>
<td class="gt_row gt_center" headers="stat_0">132 (85.2%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">23 (30.7%)</td>
<td class="gt_row gt_center" headers="stat_0">23 (14.8%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Non-accidental death</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&lt;0.001</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">80 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">56 (74.7%)</td>
<td class="gt_row gt_center" headers="stat_0">136 (87.7%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">19 (25.3%)</td>
<td class="gt_row gt_center" headers="stat_0">19 (12.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Age</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Sex at birth</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.579</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Male</td>
<td class="gt_row gt_center" headers="stat_1">62 (77.5%)</td>
<td class="gt_row gt_center" headers="stat_2">55 (73.3%)</td>
<td class="gt_row gt_center" headers="stat_0">117 (75.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Female</td>
<td class="gt_row gt_center" headers="stat_1">18 (22.5%)</td>
<td class="gt_row gt_center" headers="stat_2">20 (26.7%)</td>
<td class="gt_row gt_center" headers="stat_0">38 (24.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Risk factor for HIV
acquisition</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Country of origin</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Education attainment</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.578</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    High school</td>
<td class="gt_row gt_center" headers="stat_1">26 (32.5%)</td>
<td class="gt_row gt_center" headers="stat_2">18 (24.0%)</td>
<td class="gt_row gt_center" headers="stat_0">44 (28.4%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    No studies</td>
<td class="gt_row gt_center" headers="stat_1">4 (5.0%)</td>
<td class="gt_row gt_center" headers="stat_2">3 (4.0%)</td>
<td class="gt_row gt_center" headers="stat_0">7 (4.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Other</td>
<td class="gt_row gt_center" headers="stat_1">1 (1.3%)</td>
<td class="gt_row gt_center" headers="stat_2">3 (4.0%)</td>
<td class="gt_row gt_center" headers="stat_0">4 (2.6%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Primary school</td>
<td class="gt_row gt_center" headers="stat_1">8 (10.0%)</td>
<td class="gt_row gt_center" headers="stat_2">11 (14.7%)</td>
<td class="gt_row gt_center" headers="stat_0">19 (12.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Secondary school</td>
<td class="gt_row gt_center" headers="stat_1">23 (28.8%)</td>
<td class="gt_row gt_center" headers="stat_2">18 (24.0%)</td>
<td class="gt_row gt_center" headers="stat_0">41 (26.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    University</td>
<td class="gt_row gt_center" headers="stat_1">18 (22.5%)</td>
<td class="gt_row gt_center" headers="stat_2">22 (29.3%)</td>
<td class="gt_row gt_center" headers="stat_0">40 (25.8%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Smoker</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.589</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">37 (46.3%)</td>
<td class="gt_row gt_center" headers="stat_2">28 (37.3%)</td>
<td class="gt_row gt_center" headers="stat_0">65 (41.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">30 (37.5%)</td>
<td class="gt_row gt_center" headers="stat_2">29 (38.7%)</td>
<td class="gt_row gt_center" headers="stat_0">59 (38.1%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Hypertension</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.189</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">22 (27.5%)</td>
<td class="gt_row gt_center" headers="stat_2">25 (33.3%)</td>
<td class="gt_row gt_center" headers="stat_0">47 (30.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">14 (17.5%)</td>
<td class="gt_row gt_center" headers="stat_2">7 (9.3%)</td>
<td class="gt_row gt_center" headers="stat_0">21 (13.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">ART class</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.762</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 NNRTI</td>
<td class="gt_row gt_center" headers="stat_1">37 (46.3%)</td>
<td class="gt_row gt_center" headers="stat_2">31 (41.3%)</td>
<td class="gt_row gt_center" headers="stat_0">68 (43.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 PI</td>
<td class="gt_row gt_center" headers="stat_1">30 (37.5%)</td>
<td class="gt_row gt_center" headers="stat_2">31 (41.3%)</td>
<td class="gt_row gt_center" headers="stat_0">61 (39.4%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 INSTI</td>
<td class="gt_row gt_center" headers="stat_1">11 (13.8%)</td>
<td class="gt_row gt_center" headers="stat_2">9 (12.0%)</td>
<td class="gt_row gt_center" headers="stat_0">20 (12.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Other</td>
<td class="gt_row gt_center" headers="stat_1">2 (2.5%)</td>
<td class="gt_row gt_center" headers="stat_2">4 (5.3%)</td>
<td class="gt_row gt_center" headers="stat_0">6 (3.9%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Maximum viral load
(copies/mL)</td>
<td class="gt_row gt_center" headers="stat_1">482,966 (1,245,839)</td>
<td class="gt_row gt_center" headers="stat_2">465,781 (1,170,338)</td>
<td class="gt_row gt_center" headers="stat_0">474,651 (1,206,007)</td>
<td class="gt_row gt_center" headers="p.value">0.930</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 nadir (cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">260.8 (181.5)</td>
<td class="gt_row gt_center" headers="stat_2">202.9 (154.6)</td>
<td class="gt_row gt_center" headers="stat_0">232.4 (170.8)</td>
<td class="gt_row gt_center" headers="p.value">0.035</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 counts at ART initiation
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">309.9 (214.3)</td>
<td class="gt_row gt_center" headers="stat_2">247.0 (183.6)</td>
<td class="gt_row gt_center" headers="stat_0">279.3 (201.7)</td>
<td class="gt_row gt_center" headers="p.value">0.060</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD8 counts at ART initiation
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">1,035.2 (569.7)</td>
<td class="gt_row gt_center" headers="stat_2">986.7 (536.9)</td>
<td class="gt_row gt_center" headers="stat_0">1,011.6 (552.6)</td>
<td class="gt_row gt_center" headers="p.value">0.600</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4/CD8 ratio at ART
initiation</td>
<td class="gt_row gt_center" headers="stat_1">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="stat_2">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="stat_0">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="p.value">0.275</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 counts at Year 2 of ART
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">584.4 (265.0)</td>
<td class="gt_row gt_center" headers="stat_2">513.2 (268.7)</td>
<td class="gt_row gt_center" headers="stat_0">549.0 (268.3)</td>
<td class="gt_row gt_center" headers="p.value">0.104</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD8 counts at Year 2 of ART
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">1,028.9 (457.7)</td>
<td class="gt_row gt_center" headers="stat_2">1,038.4 (609.1)</td>
<td class="gt_row gt_center" headers="stat_0">1,033.6 (536.5)</td>
<td class="gt_row gt_center" headers="p.value">0.914</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4/CD8 ratio at Year 2 of
ART</td>
<td class="gt_row gt_center" headers="stat_1">0.7 (0.4)</td>
<td class="gt_row gt_center" headers="stat_2">0.6 (0.4)</td>
<td class="gt_row gt_center" headers="stat_0">0.6 (0.4)</td>
<td class="gt_row gt_center" headers="p.value">0.241</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Months from ART initiation to
sampling</td>
<td class="gt_row gt_center" headers="stat_1">23.7 (4.5)</td>
<td class="gt_row gt_center" headers="stat_2">24.6 (7.2)</td>
<td class="gt_row gt_center" headers="stat_0">24.1 (6.0)</td>
<td class="gt_row gt_center" headers="p.value">0.395</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Year of ART initiation</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
</tbody><tfoot>
<tr class="gt_footnotes">
<td colspan="5" class="gt_footnote"><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span>
n (%); Mean (SD)</td>
</tr>
<tr class="gt_footnotes">
<td colspan="5" class="gt_footnote"><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span>
Fisher’s exact test; Welch Two Sample t-test</td>
</tr>
</tfoot>
&#10;</table>

</div>

``` r
write_tsv(table_gts_rnaseq |> as_tibble(), "../results/tableS6_rnaseq.tsv", na = "")
```

## Table 1 for excluded participants in Olink after QC

``` r
table_gts_olink_excluded <- prepare_table1_data(!olink_included) |> build_gtsummary_table1()
```

    Rows: 187 Columns: 27
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (6): sid, group, sex, edulvl, art_reg, sid_pair
    dbl (17): nadm, mace, death_all, smoke, hta, cvmax, nadircd4, CD4_0, CD8_0, ...
    lgl  (4): age, mode, paisorigen, artini_period

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
table_gts_olink_excluded
```

<div id="jnluanggdf" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#jnluanggdf table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#jnluanggdf thead, #jnluanggdf tbody, #jnluanggdf tfoot, #jnluanggdf tr, #jnluanggdf td, #jnluanggdf th {
  border-style: none;
}
&#10;#jnluanggdf p {
  margin: 0;
  padding: 0;
}
&#10;#jnluanggdf .gt_table {
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
&#10;#jnluanggdf .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#jnluanggdf .gt_title {
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
&#10;#jnluanggdf .gt_subtitle {
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
&#10;#jnluanggdf .gt_heading {
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
&#10;#jnluanggdf .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#jnluanggdf .gt_col_headings {
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
&#10;#jnluanggdf .gt_col_heading {
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
&#10;#jnluanggdf .gt_column_spanner_outer {
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
&#10;#jnluanggdf .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#jnluanggdf .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#jnluanggdf .gt_column_spanner {
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
&#10;#jnluanggdf .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#jnluanggdf .gt_group_heading {
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
&#10;#jnluanggdf .gt_empty_group_heading {
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
&#10;#jnluanggdf .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#jnluanggdf .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#jnluanggdf .gt_row {
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
&#10;#jnluanggdf .gt_stub {
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
&#10;#jnluanggdf .gt_stub_row_group {
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
&#10;#jnluanggdf .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#jnluanggdf .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#jnluanggdf .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jnluanggdf .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#jnluanggdf .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#jnluanggdf .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#jnluanggdf .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jnluanggdf .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#jnluanggdf .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#jnluanggdf .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#jnluanggdf .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#jnluanggdf .gt_footnotes {
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
&#10;#jnluanggdf .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jnluanggdf .gt_sourcenotes {
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
&#10;#jnluanggdf .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jnluanggdf .gt_left {
  text-align: left;
}
&#10;#jnluanggdf .gt_center {
  text-align: center;
}
&#10;#jnluanggdf .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#jnluanggdf .gt_font_normal {
  font-weight: normal;
}
&#10;#jnluanggdf .gt_font_bold {
  font-weight: bold;
}
&#10;#jnluanggdf .gt_font_italic {
  font-style: italic;
}
&#10;#jnluanggdf .gt_super {
  font-size: 65%;
}
&#10;#jnluanggdf .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#jnluanggdf .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#jnluanggdf .gt_indent_1 {
  text-indent: 5px;
}
&#10;#jnluanggdf .gt_indent_2 {
  text-indent: 10px;
}
&#10;#jnluanggdf .gt_indent_3 {
  text-indent: 15px;
}
&#10;#jnluanggdf .gt_indent_4 {
  text-indent: 20px;
}
&#10;#jnluanggdf .gt_indent_5 {
  text-indent: 25px;
}
&#10;#jnluanggdf .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#jnluanggdf div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>

<table class="gt_table" data-quarto-postprocess="true"
data-quarto-disable-processing="false" data-quarto-bootstrap="false">
<colgroup>
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
<col style="width: 20%" />
</colgroup>
<thead>
<tr class="gt_col_headings">
<th id="label" class="gt_col_heading gt_columns_bottom_border gt_left"
data-quarto-table-cell-role="th"
scope="col"><strong>Characteristic</strong></th>
<th id="stat_1"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th"
scope="col"><strong>Control</strong><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="stat_2"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th" scope="col"><strong>Event</strong><span
class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="stat_0"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th" scope="col"><strong>Total</strong><span
class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
<th id="p.value"
class="gt_col_heading gt_columns_bottom_border gt_center"
data-quarto-table-cell-role="th"
scope="col"><strong>p-value</strong><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span></th>
</tr>
</thead>
<tbody class="gt_table_body">
<tr>
<td class="gt_row gt_left" headers="label">N</td>
<td class="gt_row gt_center" headers="stat_1">23 (46.9%)</td>
<td class="gt_row gt_center" headers="stat_2">26 (53.1%)</td>
<td class="gt_row gt_center" headers="stat_0">49 (100.0%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Non-AIDS defining
malignancy</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&lt;0.001</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">26 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">6 (26.1%)</td>
<td class="gt_row gt_center" headers="stat_0">32 (65.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">17 (73.9%)</td>
<td class="gt_row gt_center" headers="stat_0">17 (34.7%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Major cardiovascular
events</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.007</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">26 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">17 (73.9%)</td>
<td class="gt_row gt_center" headers="stat_0">43 (87.8%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">6 (26.1%)</td>
<td class="gt_row gt_center" headers="stat_0">6 (12.2%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Non-accidental death</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.001</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">26 (100.0%)</td>
<td class="gt_row gt_center" headers="stat_2">15 (65.2%)</td>
<td class="gt_row gt_center" headers="stat_0">41 (83.7%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_2">8 (34.8%)</td>
<td class="gt_row gt_center" headers="stat_0">8 (16.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Age</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Sex at birth</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&gt;0.999</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Male</td>
<td class="gt_row gt_center" headers="stat_1">21 (80.8%)</td>
<td class="gt_row gt_center" headers="stat_2">18 (78.3%)</td>
<td class="gt_row gt_center" headers="stat_0">39 (79.6%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Female</td>
<td class="gt_row gt_center" headers="stat_1">5 (19.2%)</td>
<td class="gt_row gt_center" headers="stat_2">5 (21.7%)</td>
<td class="gt_row gt_center" headers="stat_0">10 (20.4%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Risk factor for HIV
acquisition</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Country of origin</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Education attainment</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.351</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    High school</td>
<td class="gt_row gt_center" headers="stat_1">8 (30.8%)</td>
<td class="gt_row gt_center" headers="stat_2">5 (21.7%)</td>
<td class="gt_row gt_center" headers="stat_0">13 (26.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    No studies</td>
<td class="gt_row gt_center" headers="stat_1">2 (7.7%)</td>
<td class="gt_row gt_center" headers="stat_2">1 (4.3%)</td>
<td class="gt_row gt_center" headers="stat_0">3 (6.1%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Primary school</td>
<td class="gt_row gt_center" headers="stat_1">2 (7.7%)</td>
<td class="gt_row gt_center" headers="stat_2">4 (17.4%)</td>
<td class="gt_row gt_center" headers="stat_0">6 (12.2%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Secondary school</td>
<td class="gt_row gt_center" headers="stat_1">9 (34.6%)</td>
<td class="gt_row gt_center" headers="stat_2">4 (17.4%)</td>
<td class="gt_row gt_center" headers="stat_0">13 (26.5%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    University</td>
<td class="gt_row gt_center" headers="stat_1">5 (19.2%)</td>
<td class="gt_row gt_center" headers="stat_2">9 (39.1%)</td>
<td class="gt_row gt_center" headers="stat_0">14 (28.6%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Smoker</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&gt;0.999</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">10 (38.5%)</td>
<td class="gt_row gt_center" headers="stat_2">9 (39.1%)</td>
<td class="gt_row gt_center" headers="stat_0">19 (38.8%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">7 (26.9%)</td>
<td class="gt_row gt_center" headers="stat_2">7 (30.4%)</td>
<td class="gt_row gt_center" headers="stat_0">14 (28.6%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Hypertension</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">0.119</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    1</td>
<td class="gt_row gt_center" headers="stat_1">7 (26.9%)</td>
<td class="gt_row gt_center" headers="stat_2">7 (30.4%)</td>
<td class="gt_row gt_center" headers="stat_0">14 (28.6%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    0</td>
<td class="gt_row gt_center" headers="stat_1">4 (15.4%)</td>
<td class="gt_row gt_center" headers="stat_2">0 (0.0%)</td>
<td class="gt_row gt_center" headers="stat_0">4 (8.2%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">ART class</td>
<td class="gt_row gt_center" headers="stat_1"><br />
</td>
<td class="gt_row gt_center" headers="stat_2"><br />
</td>
<td class="gt_row gt_center" headers="stat_0"><br />
</td>
<td class="gt_row gt_center" headers="p.value">&gt;0.999</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 NNRTI</td>
<td class="gt_row gt_center" headers="stat_1">13 (50.0%)</td>
<td class="gt_row gt_center" headers="stat_2">12 (52.2%)</td>
<td class="gt_row gt_center" headers="stat_0">25 (51.0%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 PI</td>
<td class="gt_row gt_center" headers="stat_1">8 (30.8%)</td>
<td class="gt_row gt_center" headers="stat_2">7 (30.4%)</td>
<td class="gt_row gt_center" headers="stat_0">15 (30.6%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    2 NRTI + 1 INSTI</td>
<td class="gt_row gt_center" headers="stat_1">4 (15.4%)</td>
<td class="gt_row gt_center" headers="stat_2">3 (13.0%)</td>
<td class="gt_row gt_center" headers="stat_0">7 (14.3%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">    Other</td>
<td class="gt_row gt_center" headers="stat_1">1 (3.8%)</td>
<td class="gt_row gt_center" headers="stat_2">1 (4.3%)</td>
<td class="gt_row gt_center" headers="stat_0">2 (4.1%)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Maximum viral load
(copies/mL)</td>
<td class="gt_row gt_center" headers="stat_1">423,551 (735,770)</td>
<td class="gt_row gt_center" headers="stat_2">499,390 (1,156,716)</td>
<td class="gt_row gt_center" headers="stat_0">459,149 (946,924)</td>
<td class="gt_row gt_center" headers="p.value">0.789</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 nadir (cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">215.7 (149.3)</td>
<td class="gt_row gt_center" headers="stat_2">215.4 (134.6)</td>
<td class="gt_row gt_center" headers="stat_0">215.6 (140.6)</td>
<td class="gt_row gt_center" headers="p.value">0.994</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 counts at ART initiation
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">252.8 (172.1)</td>
<td class="gt_row gt_center" headers="stat_2">276.4 (182.1)</td>
<td class="gt_row gt_center" headers="stat_0">264.9 (175.6)</td>
<td class="gt_row gt_center" headers="p.value">0.665</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD8 counts at ART initiation
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">878.0 (552.7)</td>
<td class="gt_row gt_center" headers="stat_2">926.3 (435.3)</td>
<td class="gt_row gt_center" headers="stat_0">902.7 (490.7)</td>
<td class="gt_row gt_center" headers="p.value">0.753</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4/CD8 ratio at ART
initiation</td>
<td class="gt_row gt_center" headers="stat_1">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="stat_2">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="stat_0">0.3 (0.2)</td>
<td class="gt_row gt_center" headers="p.value">0.535</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4 counts at Year 2 of ART
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">536.5 (274.5)</td>
<td class="gt_row gt_center" headers="stat_2">551.6 (256.0)</td>
<td class="gt_row gt_center" headers="stat_0">544.2 (262.3)</td>
<td class="gt_row gt_center" headers="p.value">0.850</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD8 counts at Year 2 of ART
(cells/µL)</td>
<td class="gt_row gt_center" headers="stat_1">1,012.5 (597.8)</td>
<td class="gt_row gt_center" headers="stat_2">999.1 (394.2)</td>
<td class="gt_row gt_center" headers="stat_0">1,005.7 (498.3)</td>
<td class="gt_row gt_center" headers="p.value">0.930</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">CD4/CD8 ratio at Year 2 of
ART</td>
<td class="gt_row gt_center" headers="stat_1">0.6 (0.3)</td>
<td class="gt_row gt_center" headers="stat_2">0.6 (0.4)</td>
<td class="gt_row gt_center" headers="stat_0">0.6 (0.4)</td>
<td class="gt_row gt_center" headers="p.value">0.698</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Months from ART initiation to
sampling</td>
<td class="gt_row gt_center" headers="stat_1">24.6 (6.4)</td>
<td class="gt_row gt_center" headers="stat_2">23.5 (10.7)</td>
<td class="gt_row gt_center" headers="stat_0">24.1 (8.6)</td>
<td class="gt_row gt_center" headers="p.value">0.662</td>
</tr>
<tr>
<td class="gt_row gt_left" headers="label">Year of ART initiation</td>
<td class="gt_row gt_center" headers="stat_1">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_2">NA (NA)</td>
<td class="gt_row gt_center" headers="stat_0">NA (NA)</td>
<td class="gt_row gt_center" headers="p.value"><br />
</td>
</tr>
</tbody><tfoot>
<tr class="gt_footnotes">
<td colspan="5" class="gt_footnote"><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span>
n (%); Mean (SD)</td>
</tr>
<tr class="gt_footnotes">
<td colspan="5" class="gt_footnote"><span class="gt_footnote_marks"
style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>2</sup></span>
Fisher’s exact test; Welch Two Sample t-test</td>
</tr>
</tfoot>
&#10;</table>

</div>

``` r
write_tsv(table_gts_olink_excluded |> as_tibble(), "../results/tableS15_olink_excluded.tsv", na = "")
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
     backports      1.5.1   2026-04-03 [1] CRAN (R 4.6.0)
     base64enc      0.1-6   2026-02-02 [1] CRAN (R 4.6.0)
     bit            4.6.0   2025-03-06 [1] CRAN (R 4.6.0)
     bit64          4.8.2   2026-05-19 [1] CRAN (R 4.6.0)
     broom          1.0.13  2026-05-14 [1] CRAN (R 4.6.0)
     cards          0.7.1   2025-12-02 [1] CRAN (R 4.6.0)
     cardx          0.3.2   2026-02-05 [1] CRAN (R 4.6.0)
     cli            3.6.6   2026-04-09 [1] CRAN (R 4.6.0)
     commonmark     2.0.0   2025-07-07 [1] CRAN (R 4.6.0)
     crayon         1.5.3   2024-06-20 [1] CRAN (R 4.6.0)
     digest         0.6.39  2025-11-19 [1] CRAN (R 4.6.0)
     dplyr        * 1.2.1   2026-04-03 [1] CRAN (R 4.6.0)
     evaluate       1.0.5   2025-08-27 [1] CRAN (R 4.6.0)
     farver         2.1.2   2024-05-13 [1] CRAN (R 4.6.0)
     fastmap        1.2.0   2024-05-15 [1] CRAN (R 4.6.0)
     forcats      * 1.0.1   2025-09-25 [1] CRAN (R 4.6.0)
     fs             2.1.0   2026-04-18 [1] CRAN (R 4.6.0)
     generics       0.1.4   2025-05-09 [1] CRAN (R 4.6.0)
     ggplot2      * 4.0.3   2026-04-22 [1] CRAN (R 4.6.0)
     glue           1.8.1   2026-04-17 [1] CRAN (R 4.6.0)
     gt             1.3.0   2026-01-22 [1] CRAN (R 4.6.0)
     gtable         0.3.6   2024-10-25 [1] CRAN (R 4.6.0)
     gtsummary    * 2.5.0   2025-12-05 [1] CRAN (R 4.6.0)
     hms            1.1.4   2025-10-17 [1] CRAN (R 4.6.0)
     htmltools      0.5.9   2025-12-04 [1] CRAN (R 4.6.0)
     jsonlite       2.0.0   2025-03-27 [1] CRAN (R 4.6.0)
     knitr          1.51    2025-12-20 [1] CRAN (R 4.6.0)
     lifecycle      1.0.5   2026-01-08 [1] CRAN (R 4.6.0)
     litedown       0.9     2025-12-18 [1] CRAN (R 4.6.0)
     lubridate    * 1.9.5   2026-02-04 [1] CRAN (R 4.6.0)
     magrittr       2.0.5   2026-04-04 [1] CRAN (R 4.6.0)
     markdown       2.0     2025-03-23 [1] CRAN (R 4.6.0)
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
     sass           0.4.10  2025-04-11 [1] CRAN (R 4.6.0)
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
     vctrs          0.7.3   2026-04-11 [1] CRAN (R 4.6.0)
     vroom          1.7.1   2026-03-31 [1] CRAN (R 4.6.0)
     withr          3.0.2   2024-10-28 [1] CRAN (R 4.6.0)
     xfun           0.60    2026-07-09 [1] CRAN (R 4.6.0)
     xml2           1.5.2   2026-01-17 [1] CRAN (R 4.6.0)
     yaml           2.3.12  2025-12-10 [1] CRAN (R 4.6.0)

     [1] /opt/Rlibs
     [2] /usr/local/lib/R/site-library
     [3] /usr/lib/R/site-library
     [4] /usr/lib/R/library
     * ── Packages attached to the search path.

    ──────────────────────────────────────────────────────────────────────────────
