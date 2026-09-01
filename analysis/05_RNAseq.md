# 5. RNA-seq analysis


## Purpose

This script runs differential expression analysis on the RNA-seq data
across the four outcomes using DESeq2, producing the `deg_*.tsv` files
(Table S7 to S10) and the FDR-adjusted volcano plots (Fig. 2). It also
runs the cell-type deconvolution statistical analysis (Wilcoxon tests
per cell type, BH-corrected within each outcome), reproducing Table S11
and Fig. S2.

## Reproducibility

This script runs entirely from the public repository. Only
`rnaseq_transcriptomics_counts.tsv`, `metadata.tsv`, and
`cell_deconvolution_sid.csv` are needed.

``` r
library(AnnotationDbi)
library(biomaRt)
library(org.Hs.eg.db)
library(DESeq2)
library(ggpubr)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
```

## Data generation

RNA-seq reads were processed with nf-core/rnaseq v3.17.0, using
STAR/RSEM for alignment and quantification, against the GRCh38.112
Ensembl reference. `rnaseq_transcriptomics_counts.tsv` is derived from
the pipeline’s `rsem.merged.gene_counts.tsv` output.

``` bash
nextflow run \
    nf-core/rnaseq -r 3.17.0 \
    --input "samplesheet.csv" \
    --outdir "results" \
    --gtf "Homo_sapiens.GRCh38.112.gtf.gz" \
    --fasta "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz" \
    --aligner "star_rsem" \
    --igenomes_ignore \
    --genome null \
    -profile singularity \
    -c nextflow.config
```

Reference files (`Homo_sapiens.GRCh38.112.gtf.gz`,
`Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz`) are the standard
Ensembl release 112 GRCh38 assembly, publicly available from Ensembl.

## Load public data and build the dds object

Raw RNA-seq counts (genes x samples) and the public metadata file are
combined here to build the `DESeqDataSet` from scratch.

``` r
counts_mat <- read_tsv("../data/rnaseq_transcriptomics_counts.tsv") |>
    column_to_rownames("gene_ensembl_id") |>
    as.matrix()
```

``` r
coldata <- read_tsv("../data/metadata.tsv") |>
    filter(sid %in% colnames(counts_mat)) |>
    column_to_rownames("sid")

coldata <- coldata[colnames(counts_mat), , drop = FALSE] # match sample order
```

``` r
stopifnot(all(colnames(counts_mat) == rownames(coldata))) # Must be TRUE
```

``` r
dds <- DESeqDataSetFromMatrix(
    countData = round(counts_mat),
    colData   = coldata,
    design    = ~1
)
```

    converting counts to integer mode

## dds inspection

``` r
head(counts(dds)) # to check the raw counts
```

                    PTY_001 PTY_003 PTY_004 PTY_005 PTY_006 PTY_007 PTY_008 PTY_009
    ENSG00000000003      10       6      23       5       0      11      11       4
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     452     582    1465     410     687     559     523     611
    ENSG00000000457     235     326     477     116     168     149     160     146
    ENSG00000000460      98      46      88      45      46      39      43      49
    ENSG00000000938    6015    5741    9566    7065    7140    5568    9784   13625
                    PTY_010 PTY_011 PTY_014 PTY_015 PTY_016 PTY_017 PTY_018 PTY_019
    ENSG00000000003       6      12       0      16       0       9       1      10
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     482     644     552     482     456     558     560     361
    ENSG00000000457     178     323     229     351     171     354     272     285
    ENSG00000000460      37      87      76     106      67      80     101      91
    ENSG00000000938   12333    9113    9417    9283   18539   11700   16111    9186
                    PTY_020 PTY_024 PTY_025 PTY_026 PTY_027 PTY_028 PTY_029 PTY_030
    ENSG00000000003       3      15       5       3       1       0       2       7
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     549     541     420     477     540     536     649     546
    ENSG00000000457     219     241     214     256     212     128     104      76
    ENSG00000000460      51      23      46      67      62      68      52      43
    ENSG00000000938    9669    5204   13026   10165   13191   14017    6994    4428
                    PTY_031 PTY_032 PTY_034 PTY_035 PTY_036 PTY_037 PTY_038 PTY_039
    ENSG00000000003      11      24       0       0       1       9       3       5
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     508     648     596     224     598     902     710     411
    ENSG00000000457     120     329     163     147     112     260     119     426
    ENSG00000000460      16      15      55      50      59      21      89      23
    ENSG00000000938    2547    2890   11484   17437    5473    3117    5596    2169
                    PTY_040 PTY_041 PTY_042 PTY_044 PTY_045 PTY_046 PTY_048 PTY_049
    ENSG00000000003       2      10      31       5       1       5       0      14
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     461     433     708     282     483     715     604     629
    ENSG00000000457     100     103     311     181      56     192     185     355
    ENSG00000000460      11      29      19      80      78      63      63      86
    ENSG00000000938    6071    4979    2565   13559    9756    9936    9868    4826
                    PTY_050 PTY_051 PTY_052 PTY_053 PTY_054 PTY_056 PTY_057 PTY_058
    ENSG00000000003       3      16      15       8      12       3      17       4
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     573     554     406     523     887     300     593     692
    ENSG00000000457     236     206     383     170     375     175     252     172
    ENSG00000000460      76      65      60      24     113      33     120      39
    ENSG00000000938    9989    5710    6892    6274   18525   17078   10879    5554
                    PTY_059 PTY_060 PTY_061 PTY_062 PTY_065 PTY_066 PTY_068 PTY_069
    ENSG00000000003       4       1       1      16       8      11       3       6
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     530     455     334     865    1167     734     554     824
    ENSG00000000457     233     182      69     180     261     178     243      88
    ENSG00000000460      51      69       6      20      17      18      23      59
    ENSG00000000938    5906    9797    2018    4110    4381    2552    4006    5477
                    PTY_071 PTY_072 PTY_073 PTY_074 PTY_076 PTY_079 PTY_083 PTY_084
    ENSG00000000003       9      12      14       0       3      10       0       6
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     659     959     530     333     533     635     750     439
    ENSG00000000457     112     151     180     126      94     204     113     197
    ENSG00000000460      47      22      11      33      12      22      63      30
    ENSG00000000938    4041    2039    1033    8738    6867    3324    1372    2219
                    PTY_086 PTY_087 PTY_088 PTY_090 PTY_092 PTY_093 PTY_094 PTY_096
    ENSG00000000003       2       5       1       7      10       0       1       3
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     365     621     407     565     571     307     504     490
    ENSG00000000457     155      80     151      99     194      87     116     235
    ENSG00000000460      69       7      15      54      47      42      20      63
    ENSG00000000938    9804    4053    3489    4905    6237    5964    6573    3609
                    PTY_098 PTY_100 PTY_101 PTY_102 PTY_103 PTY_104 PTY_105 PTY_106
    ENSG00000000003       6       1       8       5       1      13      13       6
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     222     345     450     409     351     470     629     619
    ENSG00000000457     135     104     121      49      86     111     214      99
    ENSG00000000460      85      47      25      42      19      21      23      10
    ENSG00000000938   11048    9403    6151   11530   12654    4483    4222    2765
                    PTY_107 PTY_108 PTY_109 PTY_110 PTY_111 PTY_112 PTY_113 PTY_114
    ENSG00000000003       6       2       9       8      15       0       4       5
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     508     489     550     495     406     456     516     416
    ENSG00000000457     276     240     208     307     325      92     366     125
    ENSG00000000460      39      91      61      56      84      33      54      44
    ENSG00000000938    5496   13949    8784   10550   14079   15568    6820   11123
                    PTY_115 PTY_116 PTY_117 PTY_118 PTY_119 PTY_120 PTY_121 PTY_122
    ENSG00000000003      16       2      10       9      14       4       0       6
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     474     667     348     511     363     488     382     462
    ENSG00000000457     181     203     345     230     266     221     279     220
    ENSG00000000460      52      69      79      35      70      74      63      71
    ENSG00000000938    4580    7759    8987    5193    4822   14155   13215    9210
                    PTY_123 PTY_124 PTY_125 PTY_126 PTY_127 PTY_128 PTY_129 PTY_130
    ENSG00000000003       5       0       5       6      25       3       7       8
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     380     452     396     559     629     707     742     673
    ENSG00000000457     214     367      75     119     259     125     141     323
    ENSG00000000460      67      66      15      42      39      36      17      39
    ENSG00000000938   11001    6337   10758    7633    3255    4720    3027    3853
                    PTY_131 PTY_132 PTY_133 PTY_134 PTY_135 PTY_136 PTY_137 PTY_138
    ENSG00000000003      11       1       5       2       4       6       2       4
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     312     366     551     429     596     714     526     840
    ENSG00000000457     215      54     200     107     379      75     120     141
    ENSG00000000460      26      23      11      44      64      77       5      54
    ENSG00000000938    5694    4122    2513    8119    8889    3784    4699    5012
                    PTY_139 PTY_140 PTY_141 PTY_142 PTY_143 PTY_144 PTY_145 PTY_146
    ENSG00000000003       2       7      15       3       2       9       0       5
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     673     748    1146     578     296     684     441     559
    ENSG00000000457     145     142     224     148      55     307      95     142
    ENSG00000000460      87      58      16      26      34      20      66      40
    ENSG00000000938    7473    6044    1965    6537    8904    5176   11614   10668
                    PTY_147 PTY_148 PTY_149 PTY_150 PTY_151 PTY_152 PTY_154 PTY_156
    ENSG00000000003       6       0       7       4      10       7      19       2
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     357     450     586     278     349     775    1116     598
    ENSG00000000457     112     100     150     155     351     337     163     143
    ENSG00000000460      30      47      61      62      72      69      29      17
    ENSG00000000938    6937    9701    8143   12508   10966   10127    2631    3212
                    PTY_157 PTY_158 PTY_159 PTY_161 PTY_162 PTY_163 PTY_164 PTY_166
    ENSG00000000003       2       3       1       0       6       3       5       5
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     233     487     590     504     516     616     603     517
    ENSG00000000457     323     257     104     119      79     114      89     230
    ENSG00000000460     103      30      73      31      32      17      53      74
    ENSG00000000938   15412    5398    6563    2568    7254    3677    4869    6402
                    PTY_168 PTY_169 PTY_170 PTY_171 PTY_172 PTY_173 PTY_174 PTY_175
    ENSG00000000003       0      14       0       2       5       0       0       0
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     717     660     756     760    1076     580     217     546
    ENSG00000000457     116      72     132     149     134      84     115     251
    ENSG00000000460      32      32      76      42      50      77      88      22
    ENSG00000000938    5596    4305    3365    2992    3462    5051    5785    3409
                    PTY_176 PTY_177 PTY_178 PTY_179 PTY_180 PTY_181 PTY_182 PTY_183
    ENSG00000000003       1      12       0       1      13       0       2       2
    ENSG00000000005       0       0       0       0       0       0       0       0
    ENSG00000000419     378     640      23     244     755     248     365     326
    ENSG00000000457      70     362       0     107     182     310     145     129
    ENSG00000000460      32      15       0      43      24      60      48      84
    ENSG00000000938    7739    3501     191   13142    4054   12723   10644   13047
                    PTY_184 PTY_185 PTY_187
    ENSG00000000003      10       4       2
    ENSG00000000005       0       0       0
    ENSG00000000419     639     519     690
    ENSG00000000457     233     155     126
    ENSG00000000460      68      33      68
    ENSG00000000938    2519    6556    5217

``` r
colData(dds) # to check the sample info
```

    DataFrame with 155 rows and 26 columns
                  group      nadm      mace death_all       age         sex
            <character> <numeric> <numeric> <numeric> <logical> <character>
    PTY_001     Control         0         0         0        NA        Male
    PTY_003     Control         0         0         0        NA        Male
    PTY_004       Event         1         0         0        NA        Male
    PTY_005       Event         1         0         0        NA        Male
    PTY_006       Event         0         1         0        NA        Male
    ...             ...       ...       ...       ...       ...         ...
    PTY_182       Event         1         0         0        NA      Female
    PTY_183       Event         0         1         1        NA      Female
    PTY_184       Event         1         0         1        NA      Female
    PTY_185     Control         0         0         0        NA      Female
    PTY_187     Control         0         0         0        NA      Female
                 mode paisorigen           edulvl     smoke       hta
            <logical>  <logical>      <character> <numeric> <numeric>
    PTY_001        NA         NA      High school         0        NA
    PTY_003        NA         NA       University         1        NA
    PTY_004        NA         NA       University         0        NA
    PTY_005        NA         NA      High school         1         1
    PTY_006        NA         NA   Primary school         1         1
    ...           ...        ...              ...       ...       ...
    PTY_182        NA         NA      High school         1         1
    PTY_183        NA         NA Secondary school         0         1
    PTY_184        NA         NA Secondary school         1         1
    PTY_185        NA         NA Secondary school        NA        NA
    PTY_187        NA         NA      High school        NA        NA
                     art_reg     cvmax  nadircd4     CD4_0     CD8_0  CD4CD8_0
                 <character> <numeric> <numeric> <numeric> <numeric> <numeric>
    PTY_001    2 NRTI + 1 PI    467100       207       298      1144  0.260490
    PTY_003 2 NRTI + 1 NNRTI     66640       837       837      1129  0.741364
    PTY_004            Other    110972       416       421      1964  0.214358
    PTY_005    2 NRTI + 1 PI    204932       126       126       513  0.245614
    PTY_006 2 NRTI + 1 INSTI    113700        55        55       515  0.106796
    ...                  ...       ...       ...       ...       ...       ...
    PTY_182    2 NRTI + 1 PI   1000000       140       228       964  0.236515
    PTY_183 2 NRTI + 1 NNRTI    117527       200       200       700  0.285714
    PTY_184 2 NRTI + 1 NNRTI    189170       190       250      1810  0.138122
    PTY_185 2 NRTI + 1 NNRTI     19000        NA        NA        NA  0.409268
    PTY_187    2 NRTI + 1 PI   2000000        NA        NA        NA  0.409268
              CD4_24m   CD8_24m CD4CD8_24 months_art_sampling artini_period
            <numeric> <numeric> <numeric>           <numeric>     <logical>
    PTY_001       554      1478  0.374831             29.1333            NA
    PTY_003      1389      1566  0.886973             18.0667            NA
    PTY_004       840      1057  0.794702             21.2069            NA
    PTY_005       522      1158  0.450777             25.4333            NA
    PTY_006       419      2144  0.195429             13.3667            NA
    ...           ...       ...       ...                 ...           ...
    PTY_182       790      1370  0.576642             51.2333            NA
    PTY_183       396       680  0.582353             27.7000            NA
    PTY_184       600      1200  0.500000             24.6000            NA
    PTY_185        NA        NA        NA             26.5161            NA
    PTY_187        NA        NA        NA             26.8387            NA
               sid_pair      time olink_included rnaseq_included
            <character> <numeric>      <numeric>       <numeric>
    PTY_001     PTY_015        10              1               1
    PTY_003     PTY_049        41              1               1
    PTY_004     PTY_032        31              1               1
    PTY_005     PTY_075        21              1               1
    PTY_006     PTY_041         6              1               1
    ...             ...       ...            ...             ...
    PTY_182     PTY_185       111              0               1
    PTY_183     PTY_176        50              1               1
    PTY_184     PTY_132        24              1               1
    PTY_185     PTY_182         4              0               1
    PTY_187     PTY_133        29              0               1

``` r
design(dds) # to check the design formula
```

    ~1

## Pre-filtering

``` r
smallestGroupSize <- 3

keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds_filtered <- dds[keep, ]
```

## Help functions

``` r
# A function to run DESeq2
run_deseq <- function(dds, design_var) {
    dds_tmp <- DESeqDataSet(dds, design = as.formula(paste0("~ ", design_var)))
    dds_tmp <- DESeq(dds_tmp)
    res <- results(dds_tmp)
    list(dds = dds_tmp, res = res)
}
```

``` r
# A function for Ensembl annotation
annotate_genes <- function(ensembl_ids) {
    mart <- useEnsembl("genes", "hsapiens_gene_ensembl")
    
    getBM(
        attributes = c(
            "ensembl_gene_id",
            "external_gene_name",
            "description",
            "gene_biotype"
        ),
        filters = "ensembl_gene_id",
        values = unique(ensembl_ids),
        mart = mart
    )
}
```

``` r
tidy_results <- function(res, lfc = 0, alpha = 0.05, filter_sig = TRUE) {
    
    # 1. Convert DESeq2 results to tibble and classify genes
    res_tb <- as_tibble(res) |>
        mutate(
            gene = rownames(res),
            diffexpressed = case_when(
                log2FoldChange >  lfc & padj < alpha ~ "upregulated",
                log2FoldChange < -lfc & padj < alpha ~ "downregulated",
                TRUE ~ "not_de"
            )
        ) |>
        relocate(gene)
    
    # 2. Annotate ALL genes
    annot <- annotate_genes(res_tb$gene) |>
        mutate(
            description = sub("\\s*\\[.*\\]$", "", description)
        )
    
    # 3. Join annotation
    res_tb <- res_tb |>
        left_join(
            annot,
            by = c("gene" = "ensembl_gene_id")
        ) |>
        arrange(pvalue)
    
    # 4. Optionally filter significant genes
    if (filter_sig) {
        res_tb <- res_tb |> filter(padj < alpha)
    }
    
    return(res_tb)
}
```

``` r
# FDR-adjusted volcano plot
plot_volcano_fdr <- function(
        res_tb,
        title,
        lfc_cutoff = 1,
        fdr_cutoff = 0.05,
        label_pc_only = TRUE
) {
    cols <- c("upregulated" = "#ffad73", "downregulated" = "#26b3ff", "not_de" = "grey")
    alphas <- c("upregulated" = 1, "downregulated" = 1, "not_de" = 0.5)
    
    plot_tb <- res_tb |>
        mutate(
            genelabels = ifelse(
                diffexpressed != "not_de" &
                    (!label_pc_only | gene_biotype == "protein_coding"),
                external_gene_name,
                ""
            )
        )
    xmax <- max(abs(plot_tb$log2FoldChange), na.rm = TRUE)
    
    ggplot(
        plot_tb,
        aes(x = log2FoldChange, y = -log10(padj), fill = diffexpressed, alpha = diffexpressed)
    ) +
        geom_point(shape = 21, colour = "black") +
        geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed",
                   linewidth = 0.2, colour = "black") +
        geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed") +
        scale_fill_manual(values = cols, name = "Differential expression",
                          breaks = c("upregulated", "downregulated", "not_de"),
                          labels = c("Up", "Down", "n.s.")) +
        scale_alpha_manual(values = alphas) +
        scale_x_continuous(limits = c(-xmax, xmax)) +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.background = element_rect(colour = "black", linewidth = 0.5),
              axis.title = element_text(size = 9)) +
        labs(title = title, x = expression(log[2](FC)),
             y = expression(-log[10]("FDR-adjusted P value"))) +
        guides(alpha = "none")
}
```

## Run analyses

``` r
res_group      <- run_deseq(dds_filtered, "group")
```

    Warning in DESeqDataSet(dds, design = as.formula(paste0("~ ", design_var))):
    some variables in design formula are characters, converting to factors

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 262 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

``` r
res_nadm       <- run_deseq(dds_filtered, "nadm")
```

      the design formula contains one or more numeric variables with integer values,
      specifying a model with increasing fold change for higher values.
      did you mean for this to be a factor? if so, first convert
      this variable to a factor using the factor() function

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 235 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

``` r
res_mace       <- run_deseq(dds_filtered, "mace")
```

      the design formula contains one or more numeric variables with integer values,
      specifying a model with increasing fold change for higher values.
      did you mean for this to be a factor? if so, first convert
      this variable to a factor using the factor() function

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 257 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

``` r
res_death_all  <- run_deseq(dds_filtered, "death_all")
```

      the design formula contains one or more numeric variables with integer values,
      specifying a model with increasing fold change for higher values.
      did you mean for this to be a factor? if so, first convert
      this variable to a factor using the factor() function

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

    -- replacing outliers and refitting for 270 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)

    estimating dispersions

    fitting model and testing

## Tidy + plot

``` r
tb_group <- tidy_results(res_group$res, filter_sig = FALSE)
```

    Ensembl site unresponsive, trying asia mirror

``` r
tb_group
```

    # A tibble: 19,457 × 11
       gene       baseMean log2FoldChange lfcSE  stat   pvalue    padj diffexpressed
       <chr>         <dbl>          <dbl> <dbl> <dbl>    <dbl>   <dbl> <chr>        
     1 ENSG00000…   138.            1.55  0.212  7.28 3.32e-13 5.96e-9 upregulated  
     2 ENSG00000…  1775.            1.13  0.185  6.11 1.01e- 9 9.04e-6 upregulated  
     3 ENSG00000…    65.8           1.22  0.217  5.62 1.96e- 8 1.17e-4 upregulated  
     4 ENSG00000…   758.            1.41  0.285  4.94 7.96e- 7 3.57e-3 upregulated  
     5 ENSG00000…     6.15          1.47  0.304  4.83 1.34e- 6 4.80e-3 upregulated  
     6 ENSG00000…     9.68         -1.59  0.345 -4.62 3.87e- 6 9.72e-3 downregulated
     7 ENSG00000…    13.2          -1.43  0.310 -4.60 4.27e- 6 9.72e-3 downregulated
     8 ENSG00000…   120.            1.81  0.394  4.59 4.33e- 6 9.72e-3 upregulated  
     9 ENSG00000…    36.9           1.13  0.247  4.56 5.04e- 6 1.01e-2 upregulated  
    10 ENSG00000…  1728.            0.949 0.217  4.38 1.18e- 5 2.12e-2 upregulated  
    # ℹ 19,447 more rows
    # ℹ 3 more variables: external_gene_name <chr>, description <chr>,
    #   gene_biotype <chr>

``` r
tb_nadm  <- tidy_results(res_nadm$res, filter_sig = FALSE)
```

    Ensembl site unresponsive, trying asia mirror

``` r
tb_nadm
```

    # A tibble: 19,457 × 11
       gene        baseMean log2FoldChange lfcSE  stat  pvalue    padj diffexpressed
       <chr>          <dbl>          <dbl> <dbl> <dbl>   <dbl>   <dbl> <chr>        
     1 ENSG000001…    138.           1.35  0.234  5.77 7.95e-9 1.55e-4 upregulated  
     2 ENSG000000…     15.8         -1.88  0.335 -5.60 2.15e-8 2.09e-4 downregulated
     3 ENSG000001…    118.          -1.42  0.282 -5.02 5.17e-7 3.36e-3 downregulated
     4 ENSG000001…    243.          -1.70  0.349 -4.88 1.09e-6 5.29e-3 downregulated
     5 ENSG000001…     16.9         -1.42  0.299 -4.73 2.23e-6 8.67e-3 downregulated
     6 ENSG000001…    269.          -2.05  0.437 -4.68 2.83e-6 9.16e-3 downregulated
     7 ENSG000001…    110.          -1.02  0.222 -4.61 4.07e-6 1.13e-2 downregulated
     8 ENSG000002…     16.4         -1.51  0.331 -4.57 4.98e-6 1.21e-2 downregulated
     9 ENSG000001…     11.8         -2.48  0.570 -4.35 1.35e-5 2.71e-2 downregulated
    10 ENSG000001…     65.9         -0.926 0.214 -4.33 1.50e-5 2.71e-2 downregulated
    # ℹ 19,447 more rows
    # ℹ 3 more variables: external_gene_name <chr>, description <chr>,
    #   gene_biotype <chr>

``` r
tb_mace  <- tidy_results(res_mace$res, filter_sig = FALSE)
```

    Ensembl site unresponsive, trying useast mirror
    Ensembl site unresponsive, trying asia mirror

``` r
tb_mace
```

    # A tibble: 19,457 × 11
       gene        baseMean log2FoldChange lfcSE  stat  pvalue    padj diffexpressed
       <chr>          <dbl>          <dbl> <dbl> <dbl>   <dbl>   <dbl> <chr>        
     1 ENSG000001…   486.            2.27  0.428  5.29 1.21e-7 0.00216 upregulated  
     2 ENSG000001…    75.5           0.812 0.171  4.76 1.94e-6 0.0135  upregulated  
     3 ENSG000001…    72.3           1.21  0.259  4.68 2.93e-6 0.0135  upregulated  
     4 ENSG000001…   118.            1.67  0.358  4.65 3.31e-6 0.0135  upregulated  
     5 ENSG000001…   454.            0.870 0.188  4.62 3.75e-6 0.0135  upregulated  
     6 ENSG000000…    15.8           1.98  0.433  4.56 5.06e-6 0.0137  upregulated  
     7 ENSG000001…   118.            1.68  0.370  4.55 5.36e-6 0.0137  upregulated  
     8 ENSG000001…    60.3           1.17  0.262  4.49 7.24e-6 0.0153  upregulated  
     9 ENSG000001…     5.59          2.42  0.542  4.47 7.78e-6 0.0153  upregulated  
    10 ENSG000001…  2055.            1.46  0.329  4.45 8.52e-6 0.0153  upregulated  
    # ℹ 19,447 more rows
    # ℹ 3 more variables: external_gene_name <chr>, description <chr>,
    #   gene_biotype <chr>

``` r
tb_death <- tidy_results(res_death_all$res, filter_sig = FALSE)
```

    Ensembl site unresponsive, trying useast mirror
    Ensembl site unresponsive, trying asia mirror

``` r
tb_death
```

    # A tibble: 19,457 × 11
       gene          baseMean log2FoldChange lfcSE  stat  pvalue  padj diffexpressed
       <chr>            <dbl>          <dbl> <dbl> <dbl>   <dbl> <dbl> <chr>        
     1 ENSG00000173…   249.            -2.86 0.630 -4.54 5.68e-6 0.110 not_de       
     2 ENSG00000102…    10.6           -2.66 0.655 -4.06 4.82e-5 0.239 not_de       
     3 ENSG00000158…    25.5           -1.42 0.350 -4.06 4.99e-5 0.239 not_de       
     4 ENSG00000106…   228.            -2.03 0.507 -4.01 6.09e-5 0.239 not_de       
     5 ENSG00000182…    33.3           -1.68 0.418 -4.01 6.15e-5 0.239 not_de       
     6 ENSG00000285…    12.8           -1.32 0.336 -3.94 8.23e-5 0.247 not_de       
     7 ENSG00000123…    15.4           -1.14 0.291 -3.92 8.87e-5 0.247 not_de       
     8 ENSG00000256…     5.04          -5.01 1.30  -3.85 1.20e-4 0.292 not_de       
     9 ENSG00000131…     8.86          -1.75 0.459 -3.81 1.36e-4 0.295 not_de       
    10 ENSG00000100…   415.            -1.20 0.318 -3.77 1.66e-4 0.324 not_de       
    # ℹ 19,447 more rows
    # ℹ 3 more variables: external_gene_name <chr>, description <chr>,
    #   gene_biotype <chr>

## Export results

``` r
# tb_group |> write_tsv("../results/tableS7_deg_group.tsv")
# tb_nadm |> write_tsv("../results/tableS8_deg_nadm.tsv")
# tb_mace |> write_tsv("../results/tableS9_deg_mace.tsv")
# tb_death |> write_tsv("../results/tableS10_deg_death_all.tsv")
```

## Volcano plots

``` r
pA <- plot_volcano_fdr(tb_group, "PWH with vs without SNAEs")
pA
```

    Warning: Removed 1509 rows containing missing values or values outside the scale range
    (`geom_point()`).

![](05_RNAseq_files/figure-commonmark/unnamed-chunk-18-1.png)

``` r
pB <- plot_volcano_fdr(tb_nadm, "PWH with vs without neoplasias")
pB
```

![](05_RNAseq_files/figure-commonmark/unnamed-chunk-18-2.png)

``` r
pC <- plot_volcano_fdr(tb_mace, "PWH with vs without CV events")
pC
```

    Warning: Removed 1509 rows containing missing values or values outside the scale range
    (`geom_point()`).

![](05_RNAseq_files/figure-commonmark/unnamed-chunk-18-3.png)

``` r
pD <- plot_volcano_fdr(tb_death, "PWH with vs without mortality")
pD
```

![](05_RNAseq_files/figure-commonmark/unnamed-chunk-18-4.png)

All together in one figure:

``` r
# svg("../results/fig2_AtoD.svg", width = 7, height = 5)
ggarrange(pA, pB, pC, pD, labels = c("A", "B", "C", "D"),
          common.legend = TRUE, legend = "bottom")
```

    Warning: Removed 1509 rows containing missing values or values outside the scale range
    (`geom_point()`).
    Removed 1509 rows containing missing values or values outside the scale range
    (`geom_point()`).
    Removed 1509 rows containing missing values or values outside the scale range
    (`geom_point()`).

![](05_RNAseq_files/figure-commonmark/unnamed-chunk-19-1.png)

``` r
# dev.off()
```

## Cell deconvolution statistical analysis

This section reproduces the statistics reported in Table S11/Fig. S2
directly from the public deconvolution results file.

``` r
deconv <- read_csv("../data/cell_deconvolution_sid.csv")
metadata <- read_tsv("../data/metadata.tsv")

outcomes <- c("group", "nadm", "mace", "death_all")
diagnostic_columns <- c("sid", "P-value", "Correlation", "RMSE")
cell_types <- setdiff(names(deconv), diagnostic_columns)
```

``` r
fit_qc <- deconv |>
    summarise(
        samples = n(),
        cibersort_p_lt_0.001 = sum(`P-value` < 0.001),
        correlation_median = median(Correlation),
        correlation_q1 = quantile(Correlation, 0.25),
        correlation_q3 = quantile(Correlation, 0.75),
        rmse_median = median(RMSE),
        rmse_q1 = quantile(RMSE, 0.25),
        rmse_q3 = quantile(RMSE, 0.75)
    )
fit_qc
```

    # A tibble: 1 × 8
      samples cibersort_p_lt_0.001 correlation_median correlation_q1 correlation_q3
        <int>                <int>              <dbl>          <dbl>          <dbl>
    1     155                  155              0.811          0.785          0.836
    # ℹ 3 more variables: rmse_median <dbl>, rmse_q1 <dbl>, rmse_q3 <dbl>

``` r
celltype_qc <- deconv |>
    dplyr::select(all_of(cell_types)) |>
    pivot_longer(everything(), names_to = "cell_type", values_to = "cell_proportion") |>
    group_by(cell_type) |>
    summarise(
        mean_proportion = mean(cell_proportion),
        median_proportion = median(cell_proportion),
        q1_proportion = quantile(cell_proportion, 0.25),
        q3_proportion = quantile(cell_proportion, 0.75),
        nonzero_n = sum(cell_proportion > 0),
        zero_n = sum(cell_proportion == 0),
        zero_percent = 100 * mean(cell_proportion == 0),
        .groups = "drop"
    ) |>
    arrange(desc(mean_proportion))
celltype_qc
```

    # A tibble: 22 × 8
       cell_type       mean_proportion median_proportion q1_proportion q3_proportion
       <chr>                     <dbl>             <dbl>         <dbl>         <dbl>
     1 T cells CD8             0.367             0.365          0.297        0.434  
     2 Monocytes               0.290             0.290          0.227        0.349  
     3 NK cells resti…         0.0819            0.0757         0.0366       0.115  
     4 B cells naive           0.0639            0.0583         0.0371       0.0828 
     5 T cells CD4 me…         0.0434            0.0295         0            0.0661 
     6 NK cells activ…         0.0381            0.0320         0.0125       0.0566 
     7 Macrophages M2          0.0311            0.0279         0.0201       0.0415 
     8 T cells regula…         0.0256            0.0238         0.0123       0.0352 
     9 T cells CD4 me…         0.0129            0.00351        0            0.0196 
    10 Mast cells act…         0.00896           0              0            0.00542
    # ℹ 12 more rows
    # ℹ 3 more variables: nonzero_n <int>, zero_n <int>, zero_percent <dbl>

``` r
# group/nadm/mace/death_all: convert group's "Control"/"Event" labels to
# 0/1 first, so all four outcomes share the same 0/1 convention.
metadata_rnaseq <- metadata |>
    filter(sid %in% deconv$sid) |>
    mutate(group = if_else(group == "Event", 1, 0)) |>
    mutate(across(all_of(outcomes), as.integer))

analysis_data <- deconv |>
    left_join(metadata_rnaseq, by = "sid")
```

``` r
sample_counts <- analysis_data |>
    dplyr::select(sid, all_of(outcomes)) |>
    pivot_longer(cols = all_of(outcomes), names_to = "outcome", values_to = "event") |>
    filter(!is.na(event)) |>
    dplyr::count(outcome, event, name = "n")
sample_counts
```

    # A tibble: 8 × 3
      outcome   event     n
      <chr>     <int> <int>
    1 death_all     0   136
    2 death_all     1    19
    3 group         0    80
    4 group         1    75
    5 mace          0   132
    6 mace          1    23
    7 nadm          0   103
    8 nadm          1    52

``` r
deconv_long <- analysis_data |>
    dplyr::select(sid, all_of(cell_types), all_of(outcomes)) |>
    pivot_longer(cols = all_of(cell_types), names_to = "cell_type", values_to = "cell_proportion") |>
    pivot_longer(cols = all_of(outcomes), names_to = "outcome", values_to = "event") |>
    filter(!is.na(event))

compare_groups <- function(data) {
    event0 <- data$cell_proportion[data$event == 0]
    event1 <- data$cell_proportion[data$event == 1]
    
    test <- suppressWarnings(wilcox.test(event1, event0, paired = FALSE, exact = FALSE, correct = TRUE))
    
    combined_ranks <- rank(c(event1, event0), ties.method = "average")
    n1 <- length(event1); n0 <- length(event0)
    mann_whitney_u <- sum(combined_ranks[seq_len(n1)]) - n1 * (n1 + 1) / 2
    rank_biserial <- (2 * mann_whitney_u / (n1 * n0)) - 1
    
    tibble(
        n_event0 = n0, n_event1 = n1,
        median_event0 = median(event0), q1_event0 = quantile(event0, 0.25, names = FALSE),
        q3_event0 = quantile(event0, 0.75, names = FALSE),
        median_event1 = median(event1), q1_event1 = quantile(event1, 0.25, names = FALSE),
        q3_event1 = quantile(event1, 0.75, names = FALSE),
        median_difference = median(event1) - median(event0),
        rank_biserial = rank_biserial, p_value = test$p.value
    )
}

celltype_results <- deconv_long |>
    group_by(outcome, cell_type) |>
    group_modify(~ compare_groups(.x)) |>
    ungroup() |>
    group_by(outcome) |>
    mutate(p_adj = p.adjust(p_value, method = "BH")) |>
    ungroup() |>
    arrange(outcome, p_adj, p_value)
celltype_results
```

    # A tibble: 88 × 14
       outcome   cell_type       n_event0 n_event1 median_event0 q1_event0 q3_event0
       <chr>     <chr>              <int>    <int>         <dbl>     <dbl>     <dbl>
     1 death_all NK cells resti…      136       19      0.0775    0.0398     0.118  
     2 death_all T cells CD4 na…      136       19      0         0          0      
     3 death_all T cells CD8          136       19      0.359     0.296      0.425  
     4 death_all Dendritic cell…      136       19      0.000906  0.000241   0.00222
     5 death_all T cells CD4 me…      136       19      0.00378   0          0.0210 
     6 death_all Mast cells act…      136       19      0         0          0.00699
     7 death_all Monocytes            136       19      0.289     0.223      0.342  
     8 death_all Macrophages M1       136       19      0         0          0      
     9 death_all T cells follic…      136       19      0.00260   0          0.0102 
    10 death_all Mast cells res…      136       19      0.00237   0          0.0100 
    # ℹ 78 more rows
    # ℹ 7 more variables: median_event1 <dbl>, q1_event1 <dbl>, q3_event1 <dbl>,
    #   median_difference <dbl>, rank_biserial <dbl>, p_value <dbl>, p_adj <dbl>

``` r
nominal_findings <- celltype_results |> filter(p_value < 0.05)
fdr_findings <- celltype_results |> filter(p_adj < 0.05)
nominal_findings
```

    # A tibble: 7 × 14
      outcome   cell_type        n_event0 n_event1 median_event0 q1_event0 q3_event0
      <chr>     <chr>               <int>    <int>         <dbl>     <dbl>     <dbl>
    1 death_all NK cells resting      136       19       0.0775     0.0398   0.118  
    2 group     NK cells resting       80       75       0.0847     0.0536   0.124  
    3 mace      B cells memory        132       23       0          0        0.0110 
    4 mace      Neutrophils           132       23       0          0        0      
    5 mace      Macrophages M2        132       23       0.0296     0.0214   0.0427 
    6 nadm      T cells CD4 mem…      103       52       0.00462    0        0.0277 
    7 nadm      Mast cells acti…      103       52       0          0        0.00918
    # ℹ 7 more variables: median_event1 <dbl>, q1_event1 <dbl>, q3_event1 <dbl>,
    #   median_difference <dbl>, rank_biserial <dbl>, p_value <dbl>, p_adj <dbl>

``` r
fdr_findings
```

    # A tibble: 0 × 14
    # ℹ 14 variables: outcome <chr>, cell_type <chr>, n_event0 <int>,
    #   n_event1 <int>, median_event0 <dbl>, q1_event0 <dbl>, q3_event0 <dbl>,
    #   median_event1 <dbl>, q1_event1 <dbl>, q3_event1 <dbl>,
    #   median_difference <dbl>, rank_biserial <dbl>, p_value <dbl>, p_adj <dbl>

``` r
# celltype_results |> write_tsv("../results/tableS11_cell_deconvolution_results.tsv")
```

## Deconvolution figure (Fig. S2)

``` r
group_plot_data <- deconv_long |>
    filter(outcome == "group") |>
    mutate(
        event = factor(event, levels = c(0, 1), labels = c("No event", "Event")),
        cell_type = as.character(cell_type)
    )

cell_order <- group_plot_data |>
    group_by(cell_type) |>
    summarise(mean_abundance = mean(cell_proportion), .groups = "drop") |>
    arrange(desc(mean_abundance)) |>
    pull(cell_type)

group_plot_data <- group_plot_data |>
    mutate(cell_type = factor(cell_type, levels = cell_order))

label_positions <- group_plot_data |>
    mutate(cell_type = as.character(cell_type)) |>
    group_by(cell_type) |>
    summarise(y_position = max(cell_proportion) * 1.12 + 0.002, .groups = "drop")

nominal_labels <- celltype_results |>
    filter(outcome == "group", p_value < 0.05) |>
    dplyr::select(cell_type, p_value) |>
    mutate(cell_type = as.character(cell_type)) |>
    left_join(label_positions, by = "cell_type") |>
    mutate(
        cell_type = factor(cell_type, levels = cell_order),
        p_label = format.pval(p_value, digits = 1, eps = 0.001)
    )

deconvolution_figure <- ggplot(group_plot_data, aes(x = cell_type, y = cell_proportion, fill = event)) +
    geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"),
                 width = 0.68, outlier.size = 0.7, outlier.alpha = 0.5, linewidth = 0.35) +
    geom_text(data = nominal_labels, aes(x = cell_type, y = y_position, label = p_label),
              inherit.aes = FALSE, size = 3.2, fontface = "bold") +
    scale_fill_manual(values = c("No event" = "#0072B2", "Event" = "#D55E00")) +
    scale_y_sqrt(labels = scales::label_percent(accuracy = 1),
                 expand = expansion(mult = c(0.02, 0.14))) +
    labs(x = NULL, y = "Estimated cell proportion", fill = NULL) +
    theme_bw(base_size = 10) +
    theme(
        legend.position = "top",
        legend.background = element_rect(colour = "black", fill = "white", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8),
        axis.title.y = element_text(margin = margin(r = 8)),
        plot.caption = element_text(hjust = 0, size = 8)
    )

deconvolution_figure
```

![](05_RNAseq_files/figure-commonmark/unnamed-chunk-28-1.png)

``` r
# svg("../results/figS2_cell_deconvolution.svg", width = 7, height = 7)
deconvolution_figure
```

![](05_RNAseq_files/figure-commonmark/unnamed-chunk-28-2.png)

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
     package              * version   date (UTC) lib source
     abind                  1.4-8     2024-09-12 [1] CRAN (R 4.6.0)
     AnnotationDbi        * 1.74.0    2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     backports              1.5.1     2026-04-03 [1] CRAN (R 4.6.0)
     Biobase              * 2.72.0    2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     BiocFileCache          3.2.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     BiocGenerics         * 0.58.1    2026-05-14 [1] Bioconductor 3.23 (R 4.6.0)
     BiocParallel           1.46.0    2026-04-29 [1] Bioconductor 3.23 (R 4.6.0)
     biomaRt              * 2.68.0    2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     Biostrings             2.80.1    2026-05-22 [1] Bioconductor 3.23 (R 4.6.0)
     bit                    4.6.0     2025-03-06 [1] CRAN (R 4.6.0)
     bit64                  4.8.2     2026-05-19 [1] CRAN (R 4.6.0)
     blob                   1.3.0     2026-01-14 [1] CRAN (R 4.6.0)
     broom                  1.0.13    2026-05-14 [1] CRAN (R 4.6.0)
     cachem                 1.1.0     2024-05-16 [1] CRAN (R 4.6.0)
     car                    3.1-5     2026-02-03 [1] CRAN (R 4.6.0)
     carData                3.0-6     2026-01-30 [1] CRAN (R 4.6.0)
     cli                    3.6.6     2026-04-09 [1] CRAN (R 4.6.0)
     codetools              0.2-20    2024-03-31 [4] CRAN (R 4.4.0)
     cowplot                1.2.0     2025-07-07 [1] CRAN (R 4.6.0)
     crayon                 1.5.3     2024-06-20 [1] CRAN (R 4.6.0)
     curl                   7.1.0     2026-04-22 [1] CRAN (R 4.6.0)
     DBI                    1.3.0     2026-02-25 [1] CRAN (R 4.6.0)
     dbplyr                 2.5.2     2026-02-13 [1] CRAN (R 4.6.0)
     DelayedArray           0.38.1    2026-04-30 [1] Bioconductor 3.23 (R 4.6.0)
     DESeq2               * 1.52.0    2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     digest                 0.6.39    2025-11-19 [1] CRAN (R 4.6.0)
     dplyr                * 1.2.1     2026-04-03 [1] CRAN (R 4.6.0)
     evaluate               1.0.5     2025-08-27 [1] CRAN (R 4.6.0)
     farver                 2.1.2     2024-05-13 [1] CRAN (R 4.6.0)
     fastmap                1.2.0     2024-05-15 [1] CRAN (R 4.6.0)
     filelock               1.0.3     2023-12-11 [1] CRAN (R 4.6.0)
     forcats              * 1.0.1     2025-09-25 [1] CRAN (R 4.6.0)
     Formula                1.2-5     2023-02-24 [1] CRAN (R 4.6.0)
     generics             * 0.1.4     2025-05-09 [1] CRAN (R 4.6.0)
     GenomicRanges        * 1.64.0    2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     ggplot2              * 4.0.3     2026-04-22 [1] CRAN (R 4.6.0)
     ggpubr               * 0.6.3     2026-02-24 [1] CRAN (R 4.6.0)
     ggrepel              * 0.9.8     2026-03-17 [1] CRAN (R 4.6.0)
     ggsignif               0.6.4     2022-10-13 [1] CRAN (R 4.6.0)
     glue                   1.8.1     2026-04-17 [1] CRAN (R 4.6.0)
     gridExtra              2.3       2017-09-09 [1] CRAN (R 4.6.0)
     gtable                 0.3.6     2024-10-25 [1] CRAN (R 4.6.0)
     hms                    1.1.4     2025-10-17 [1] CRAN (R 4.6.0)
     htmltools              0.5.9     2025-12-04 [1] CRAN (R 4.6.0)
     httr                   1.4.8     2026-02-13 [1] CRAN (R 4.6.0)
     httr2                  1.2.2     2025-12-08 [1] CRAN (R 4.6.0)
     IRanges              * 2.46.0    2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     jsonlite               2.0.0     2025-03-27 [1] CRAN (R 4.6.0)
     KEGGREST               1.52.0    2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     knitr                  1.51      2025-12-20 [1] CRAN (R 4.6.0)
     labeling               0.4.3     2023-08-29 [1] CRAN (R 4.6.0)
     lattice                0.22-9    2026-02-09 [1] CRAN (R 4.6.0)
     lifecycle              1.0.5     2026-01-08 [1] CRAN (R 4.6.0)
     locfit                 1.5-9.12  2025-03-05 [1] CRAN (R 4.6.0)
     lubridate            * 1.9.5     2026-02-04 [1] CRAN (R 4.6.0)
     magrittr               2.0.5     2026-04-04 [1] CRAN (R 4.6.0)
     Matrix                 1.7-5     2026-03-21 [4] CRAN (R 4.5.3)
     MatrixGenerics       * 1.24.0    2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     matrixStats          * 1.5.0     2025-01-07 [1] CRAN (R 4.6.0)
     memoise                2.0.1     2021-11-26 [1] CRAN (R 4.6.0)
     org.Hs.eg.db         * 3.23.1    2026-05-22 [1] Bioconductor
     otel                   0.2.0     2025-08-29 [1] CRAN (R 4.6.0)
     pheatmap             * 1.0.13    2025-06-05 [1] CRAN (R 4.6.0)
     pillar                 1.11.1    2025-09-17 [1] CRAN (R 4.6.0)
     pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.6.0)
     png                    0.1-9     2026-03-15 [1] CRAN (R 4.6.0)
     prettyunits            1.2.0     2023-09-24 [1] CRAN (R 4.6.0)
     progress               1.2.3     2023-12-06 [1] CRAN (R 4.6.0)
     purrr                * 1.2.2     2026-04-10 [1] CRAN (R 4.6.0)
     R6                     2.6.1     2025-02-15 [1] CRAN (R 4.6.0)
     rappdirs               0.3.4     2026-01-17 [1] CRAN (R 4.6.0)
     RColorBrewer         * 1.1-3     2022-04-03 [1] CRAN (R 4.6.0)
     Rcpp                   1.1.1-1.1 2026-04-24 [1] CRAN (R 4.6.0)
     readr                * 2.2.0     2026-02-19 [1] CRAN (R 4.6.0)
     rlang                  1.2.0     2026-04-06 [1] CRAN (R 4.6.0)
     rmarkdown              2.31      2026-03-26 [1] CRAN (R 4.6.0)
     RSQLite                3.53.1    2026-05-23 [1] CRAN (R 4.6.0)
     rstatix                0.7.3     2025-10-18 [1] CRAN (R 4.6.0)
     rstudioapi             0.18.0    2026-01-16 [1] CRAN (R 4.6.0)
     S4Arrays               1.12.0    2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     S4Vectors            * 0.50.1    2026-05-13 [1] Bioconductor 3.23 (R 4.6.0)
     S7                     0.2.2     2026-04-22 [1] CRAN (R 4.6.0)
     scales                 1.4.0     2025-04-24 [1] CRAN (R 4.6.0)
     Seqinfo              * 1.2.0     2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     sessioninfo            1.2.3     2025-02-05 [1] CRAN (R 4.6.0)
     SparseArray            1.12.2    2026-05-01 [1] Bioconductor 3.23 (R 4.6.0)
     stringi                1.8.7     2025-03-27 [1] CRAN (R 4.6.0)
     stringr              * 1.6.0     2025-11-04 [1] CRAN (R 4.6.0)
     SummarizedExperiment * 1.42.0    2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     tibble               * 3.3.1     2026-01-11 [1] CRAN (R 4.6.0)
     tidyr                * 1.3.2     2025-12-19 [1] CRAN (R 4.6.0)
     tidyselect             1.2.1     2024-03-11 [1] CRAN (R 4.6.0)
     tidyverse            * 2.0.0     2023-02-22 [1] CRAN (R 4.6.0)
     timechange             0.4.0     2026-01-29 [1] CRAN (R 4.6.0)
     tzdb                   0.5.0     2025-03-15 [1] CRAN (R 4.6.0)
     utf8                   1.2.6     2025-06-08 [1] CRAN (R 4.6.0)
     vctrs                  0.7.3     2026-04-11 [1] CRAN (R 4.6.0)
     vroom                  1.7.1     2026-03-31 [1] CRAN (R 4.6.0)
     withr                  3.0.2     2024-10-28 [1] CRAN (R 4.6.0)
     xfun                   0.60      2026-07-09 [1] CRAN (R 4.6.0)
     xml2                   1.5.2     2026-01-17 [1] CRAN (R 4.6.0)
     XVector                0.52.0    2026-04-28 [1] Bioconductor 3.23 (R 4.6.0)
     yaml                   2.3.12    2025-12-10 [1] CRAN (R 4.6.0)

     [1] /opt/Rlibs
     [2] /usr/local/lib/R/site-library
     [3] /usr/lib/R/site-library
     [4] /usr/lib/R/library
     * ── Packages attached to the search path.

    ──────────────────────────────────────────────────────────────────────────────
