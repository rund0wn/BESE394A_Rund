## Experimental design

    # Read data
    urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
    path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
    GSE198256_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

    # Read Meta data
    library(GEOquery)

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
    ##     as.data.frame, basename, cbind, colnames, dirname, do.call,
    ##     duplicated, eval, evalq, get, grep, grepl, intersect, is.unsorted,
    ##     lapply, mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    ##     pmin.int, rank, rbind, rownames, sapply, setdiff, sort, table,
    ##     tapply, union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

    gds <- getGEO("GSE198256")

    ## Found 1 file(s)

    ## GSE198256_series_matrix.txt.gz

    Meta_GSE198256 <- pData(gds$GSE198256_series_matrix.txt.gz@phenoData)
    Group <- Meta_GSE198256[,c("disease state:ch1")]

    dim(GSE198256_count)

    ## [1] 39376    34

    Group

    ##  [1] "Healthy"                  "Healthy"                 
    ##  [3] "Healthy"                  "Healthy"                 
    ##  [5] "Healthy"                  "Healthy"                 
    ##  [7] "Healthy"                  "Healthy"                 
    ##  [9] "Healthy"                  "Healthy"                 
    ## [11] "Healthy"                  "Covid19: Acute infection"
    ## [13] "Covid19: Acute infection" "Covid19: Acute infection"
    ## [15] "Covid19: Acute infection" "Covid19: Acute infection"
    ## [17] "Covid19: Acute infection" "Covid19: Acute infection"
    ## [19] "Covid19: Recovery 3Mo"    "Covid19: Recovery 3Mo"   
    ## [21] "Covid19: Recovery 3Mo"    "Covid19: Recovery 3Mo"   
    ## [23] "Covid19: Recovery 3Mo"    "Covid19: Recovery 3Mo"   
    ## [25] "Covid19: Recovery 6Mo"    "Covid19: Recovery 6Mo"   
    ## [27] "Covid19: Recovery 6Mo"    "Covid19: Recovery 6Mo"   
    ## [29] "Covid19: Recovery 6Mo"    "Covid19: Recovery 6Mo"   
    ## [31] "Covid19: Recovery 6Mo"    "Covid19: Recovery 6Mo"   
    ## [33] "Covid19: Recovery 6Mo"    "Covid19: Recovery 6Mo"

## Limma: Normalize and set design

    # set DGE class
    require(limma)

    ## Loading required package: limma

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

    require(edgeR)

    ## Loading required package: edgeR

    dge <- DGEList(counts=GSE198256_count)

    # Make sure on the metadata
    rownames(Meta_GSE198256)==colnames(GSE198256_count)

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [31] TRUE TRUE TRUE TRUE

    Group[Group=="Covid19: Acute infection"] <- "Covid19AI"
    Group[Group=="Covid19: Recovery 3Mo"] <- "Covid193Mo"
    Group[Group=="Covid19: Recovery 6Mo"] <- "Covid196Mo"
    design <- model.matrix(~ Group )

    # Filter
    keep <- filterByExpr(dge, design=design)
    dge <- dge[keep,,keep.lib.sizes=FALSE]

    # Normalization
    dge <- calcNormFactors(dge)

## Limma: Voom or Trend?

-   Noticeable variability in sequencing depth among the samples. The
    read counts range from as low as about 8.8 million to as high as
    about 32 million
-   limma-voom is the recommended approach as it will identify
    differentially expressed genes with an adjustment for the
    variability associated with each read count

<!-- -->

    total_reads <- colSums(GSE198256_count)
    boxplot(log10(total_reads) ~ Group, main = "Log10 Total Reads per Group", ylab = "Log10 Total Reads", xlab = "Group", las = 2)

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/Explore%20Data-1.png)

    barplot(total_reads,
            main = "Total Reads per Sample",
            ylab = "Total Reads",
            xlab = "Samples",
            las = 2,
            col = "skyblue",
            names.arg = colnames(GSE198256_count))

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/Explore%20Data-2.png)

    # Trend (for comparison)
    logCPM <- cpm(dge, log=TRUE, prior.count=3)
    fitT <- lmFit(logCPM, design)
    fitT <- eBayes(fitT, trend=TRUE)
    topTable(fitT, coef=ncol(design))

    ##               logFC    AveExpr         t      P.Value    adj.P.Val         B
    ## 29090      3.287737  1.4550637  7.827982 2.035253e-09 3.339647e-05 11.121251
    ## 2919      -3.427478  1.9415851 -5.894407 8.189137e-07 4.194087e-03  5.690910
    ## 727684    -2.868134 -1.1139954 -5.852319 9.348748e-07 4.194087e-03  5.569920
    ## 105379599 -2.518509 -0.8079827 -5.823881 1.022387e-06 4.194087e-03  5.488153
    ## 1435      -3.960957  0.1185299 -5.567817 2.287780e-06 4.263485e-03  4.751644
    ## 7013       1.420032  3.8470271  5.556059 2.373908e-06 4.263485e-03  4.717831
    ## 120534     1.482396  3.8701683  5.506442 2.774453e-06 4.263485e-03  4.575161
    ## 92935      2.590940  0.8298893  5.486884 2.950258e-06 4.263485e-03  4.518935
    ## 112935969  2.537765  0.3804905  5.486785 2.951173e-06 4.263485e-03  4.518651
    ## 59274     -1.718493  4.6857453 -5.480403 3.010933e-06 4.263485e-03  4.500304

    resultsTrend <- topTable(fitT, coef=ncol(design), number=Inf)

    # Voom
    v <- voom(dge, design, plot=TRUE)

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/Voom%20or%20Trend-1.png)

    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    topTable(fit, coef=ncol(design))

    ##               logFC    AveExpr         t      P.Value   adj.P.Val        B
    ## 29090      4.161867  1.2216892  6.416136 1.684408e-07 0.002657002 6.030796
    ## 727684    -4.207169 -2.3402062 -6.206378 3.238469e-07 0.002657002 5.670244
    ## 2919      -3.749755  1.7717956 -5.832236 1.042316e-06 0.003642824 5.417822
    ## 5778      -2.157988  2.6630123 -5.720712 1.477187e-06 0.003642824 5.206465
    ## 101928173 -3.668782 -0.2985125 -5.808920 1.121143e-06 0.003642824 5.067868
    ## 105376281 -3.538534 -0.2886424 -5.704496 1.554011e-06 0.003642824 4.901121
    ## 105379599 -3.707908 -1.5669060 -5.924411 7.813501e-07 0.003642824 4.812158
    ## 123879    -1.672692  3.4958284 -5.523271 2.737879e-06 0.003945484 4.539344
    ## 59274     -1.727872  4.6755229 -5.527930 2.698330e-06 0.003945484 4.433349
    ## 57862     -1.339288  4.8662912 -5.506472 2.885356e-06 0.003945484 4.357617

    resultsVoom <- topTable(fit, coef=ncol(design), number=Inf)

## ACTIVITY 1:

-   How would you compare the results between voom and trend?

1.  Differences in LogFC: expected due to the different approaches to
    modeling the mean-variance relationship and normalization. For
    instance, gene 29090 shows a higher logFC in the voom results
    compared to trend.

2.  Differences in Statistical Significance: the statistical
    significance of the differential expression varies between the two
    methods. This reflects the different weighting and normalization
    procedures.

3.  Gene Rankings: while some genes appear in the top results of both
    methods (e.g. 29090, 2919, 727684), their rankings and statistical
    metrics vary. This suggests that each method may prioritize
    different aspects of the variance.

-   Is it required to run more than analysis?

No, an exploratory analysis can be run on the data (such as plotting
read counts above) to decide which method would be more appropriate for
the distribution of the data.

-   What exactly are we asking with this differential expression?

1.  Gene Expression Differences: the primary question is whether there
    are significant differences in gene expression levels between the
    compared groups or conditions in the study

2.  Biological Significance: beyond statistical significance, we’re also
    interested in the biological relevance of these differences—how
    changes in gene expression levels might relate to the biological or
    clinical states being studied

<!-- -->

    mergedResults <- merge(resultsTrend, resultsVoom, by="row.names", suffixes=c("_Trend", "_Voom"))

    # Plot difference
    logFC_range <- range(c(resultsTrend$logFC, resultsVoom$logFC))
    pVal_range <- range(c(-log10(resultsTrend$P.Value), -log10(resultsVoom$P.Value)))
    logFC_range <- logFC_range + c(-1, 1) * diff(logFC_range) * 0.05
    pVal_range <- pVal_range + c(-1, 1) * diff(pVal_range) * 0.05

    # Plot for Trend
    plot(resultsTrend$logFC, -log10(resultsTrend$P.Value),
         xlim=logFC_range, ylim=pVal_range,
         main="Volcano Plot [LogFC] - Trend",
         xlab="LogFC", ylab="-log10 P-Value",
         col=ifelse(resultsTrend$adj.P.Val < 0.05, "red", "black"), pch=20,
         las=1)
    abline(h=-log10(0.05), col="blue", lty=2)

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/ACTIVITY%201-1.png)

    # Plot for Voom
    plot(resultsVoom$logFC, -log10(resultsVoom$P.Value),
         xlim=logFC_range, ylim=pVal_range,
         main="Volcano Plot [LogFC] - Voom",
         xlab="LogFC", ylab="-log10 P-Value",
         col=ifelse(resultsVoom$adj.P.Val < 0.05, "blue", "grey"), pch=20,
         las=1)
    abline(h=-log10(0.05), col="red", lty=2)

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/ACTIVITY%201-2.png)

## ACTIVITY 2:

-   Plan the next analysis: questions, steps,…

### Objective

To explore the biological implications of gene expression changes
observed in monocytes across different stages related to COVID-19
infection and recovery.

-   ORA
-   GSEA

<!-- -->

    library(clusterProfiler)

    ## 

    ## clusterProfiler v4.10.0  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
    ## 
    ## If you use clusterProfiler in published research, please cite:
    ## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    library(msigdbr)
    library(org.Hs.eg.db)

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     rename

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     I, expand.grid, unname

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     slice

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     select

    ## 

    library(magrittr)

    # Add more contrasts
    v <- voom(dge, design, plot=TRUE)

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/ACTIVITY%202-1.png)

    colnames(design) <- c("Intercept","Covid196Mo","Covid19AI","Healthy")
    fit <- lmFit(v, design)
    contrast.matrix <- makeContrasts(Covid19AI-Healthy, Healthy, 
                                     Covid196Mo-Healthy,    
                                     levels=design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    topTable(fit2) 

    ##           Covid19AI...Healthy     Healthy Covid196Mo...Healthy  AveExpr
    ## 55022               -3.974775 -0.80182133          -0.23414173 6.588543
    ## 57118               -2.453259 -0.20757069           0.30146433 4.786936
    ## 7045                -1.612767 -0.08140918           0.15375768 8.237514
    ## 920                 -2.126645 -0.32080239           0.22487968 7.254670
    ## 4082                -2.051550 -0.44937687          -0.01858625 6.958858
    ## 8460                 4.169807 -1.30913080          -0.92296012 2.408063
    ## 9308                -3.244526 -1.68099602           0.36905444 5.994052
    ## 3553                -2.923238 -0.85005431           0.50525122 7.018174
    ## 7071                -1.093365 -0.45610552           0.56332786 8.148792
    ## 101928032           -2.905109 -0.89842287          -1.18620270 2.493322
    ##                  F      P.Value    adj.P.Val
    ## 55022     34.89562 6.646463e-11 7.307875e-07
    ## 57118     34.15000 8.907155e-11 7.307875e-07
    ## 7045      31.75166 2.359435e-10 1.290532e-06
    ## 920       30.96399 3.286371e-10 1.348152e-06
    ## 4082      28.18133 1.113262e-09 3.653503e-06
    ## 8460      27.22403 1.726628e-09 4.722039e-06
    ## 9308      26.74856 2.155542e-09 5.052899e-06
    ## 3553      25.56334 3.791849e-09 7.731782e-06
    ## 7071      25.29712 4.315026e-09 7.731782e-06
    ## 101928032 25.11693 4.711915e-09 7.731782e-06

    topTable(fit2,coef=1) 

    ##           logFC   AveExpr         t      P.Value    adj.P.Val         B
    ## 55022 -3.974775  6.588543 -8.456265 3.458570e-10 4.719191e-06 13.224973
    ## 57118 -2.453259  4.786936 -8.163766 8.189639e-10 4.719191e-06 12.402454
    ## 7045  -1.612767  8.237514 -8.146181 8.627931e-10 4.719191e-06 12.245184
    ## 920   -2.126645  7.254670 -7.800212 2.423153e-09 9.940381e-06 11.239943
    ## 4082  -2.051550  6.958858 -7.555063 5.075351e-09 1.665629e-05 10.511915
    ## 8460   4.169807  2.408063  7.401420 8.091020e-09 2.212759e-05 10.092079
    ## 241    2.033411  5.227170  7.150635 1.740173e-08 3.597256e-05  9.374371
    ## 2289   2.782224  7.148611  7.138107 1.808290e-08 3.597256e-05  9.244277
    ## 7850   5.988983  4.521831  7.044320 2.411470e-08 3.597256e-05  9.097811
    ## 6280   1.490481 11.830660  7.063257 2.275175e-08 3.597256e-05  9.015975

    topTable(fit2,coef=2) 

    ##               logFC    AveExpr         t      P.Value   adj.P.Val        B
    ## 29090      4.161867  1.2216892  6.416136 1.684408e-07 0.002657002 6.030796
    ## 727684    -4.207169 -2.3402062 -6.206378 3.238469e-07 0.002657002 5.670244
    ## 2919      -3.749755  1.7717956 -5.832236 1.042316e-06 0.003642824 5.417822
    ## 5778      -2.157988  2.6630123 -5.720712 1.477187e-06 0.003642824 5.206465
    ## 101928173 -3.668782 -0.2985125 -5.808920 1.121143e-06 0.003642824 5.067868
    ## 105376281 -3.538534 -0.2886424 -5.704496 1.554011e-06 0.003642824 4.901121
    ## 105379599 -3.707908 -1.5669060 -5.924411 7.813501e-07 0.003642824 4.812158
    ## 123879    -1.672692  3.4958284 -5.523271 2.737879e-06 0.003945484 4.539344
    ## 59274     -1.727872  4.6755229 -5.527930 2.698330e-06 0.003945484 4.433349
    ## 57862     -1.339288  4.8662912 -5.506472 2.885356e-06 0.003945484 4.357617

    topTable(fit2,coef=3)

    ##                logFC    AveExpr         t      P.Value adj.P.Val         B
    ## 3008       1.1618793  7.8027241  4.907330 1.858358e-05 0.1645265 2.5845782
    ## 8971       1.0875381  6.4957990  4.882675 2.005320e-05 0.1645265 2.5443619
    ## 1408       0.9498234  4.0373516  4.697746 3.542046e-05 0.1743946 2.1553583
    ## 654322    -1.7825313  0.4373567 -4.505206 6.376792e-05 0.1743946 1.4308976
    ## 8349       0.8532691  7.0057546  4.391223 9.009259e-05 0.1984745 1.0758274
    ## 100188987  3.2110970 -1.2120058  4.601462 4.755685e-05 0.1743946 1.0500846
    ## 85453      2.2107700  0.2883908  4.254885 1.358308e-04 0.2275806 0.8415751
    ## 80014     -0.8687124  3.5332089 -4.131811 1.961991e-04 0.2397987 0.6409434
    ## 3021       0.6324829  8.8828672  4.247933 1.386925e-04 0.2275806 0.6060103
    ## 105372599  2.8142572 -2.2719391  4.548985 5.581255e-05 0.1743946 0.5872291

## ORA

    library(org.Hs.eg.db)
    library(dplyr)

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    library(msigdbr)
    library(clusterProfiler)

    ENSEMBL_vector <- mapIds(
      org.Hs.eg.db,
      keys = rownames(GSE198256_count),
      keytype = "ENTREZID",
      column = "ENSEMBL",
      multiVals = "first")

    ## 'select()' returned 1:many mapping between keys and columns

    gene_key_df <- data.frame(
      ensembl_id = ENSEMBL_vector,
      entrez_id = names(ENSEMBL_vector),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(!is.na(ensembl_id))

    # Step 1: determine genes of interest.
    diff_table <- topTable(fit2,coef=1,p.value=0.01,number=10000) 
    genes_dif<- rownames(diff_table)

    # Step 2: determine background.
    background_set <- unique(rownames(logCPM))

    # Step 3: Determine gene sets.
    hs_msigdb_df <- msigdbr(species = "Homo sapiens")
    head(hs_msigdb_df)

    ## # A tibble: 6 × 15
    ##   gs_cat gs_subcat      gs_name        gene_symbol entrez_gene ensembl_gene   
    ##   <chr>  <chr>          <chr>          <chr>             <int> <chr>          
    ## 1 C3     MIR:MIR_Legacy AAACCAC_MIR140 ABCC4             10257 ENSG00000125257
    ## 2 C3     MIR:MIR_Legacy AAACCAC_MIR140 ABRAXAS2          23172 ENSG00000165660
    ## 3 C3     MIR:MIR_Legacy AAACCAC_MIR140 ACTN4                81 ENSG00000130402
    ## 4 C3     MIR:MIR_Legacy AAACCAC_MIR140 ACTN4                81 ENSG00000282844
    ## 5 C3     MIR:MIR_Legacy AAACCAC_MIR140 ACVR1                90 ENSG00000115170
    ## 6 C3     MIR:MIR_Legacy AAACCAC_MIR140 ADAM9              8754 ENSG00000168615
    ## # ℹ 9 more variables: human_gene_symbol <chr>, human_entrez_gene <int>,
    ## #   human_ensembl_gene <chr>, gs_id <chr>, gs_pmid <chr>, gs_geoid <chr>,
    ## #   gs_exact_source <chr>, gs_url <chr>, gs_description <chr>

    hs_kegg_df <- hs_msigdb_df %>%
      dplyr::filter(
        gs_cat == "C2", # This is to filter only to the C2 curated gene sets
        gs_subcat == "CP:KEGG" # This is because we only want KEGG pathways
      )

    # Step 4: conduct ORA.
    kegg_ora_results <- enricher(
      gene = genes_dif, # A vector of your genes of interest
      pvalueCutoff = 0.1, # Can choose a FDR cutoff
      pAdjustMethod = "BH", # Method to be used for multiple testing correction
      universe = background_set,
      TERM2GENE = dplyr::select(
        hs_kegg_df,
        gs_name,
        human_entrez_gene))


    # Step 5: Visualize / explore
    enrich_plot <- enrichplot::dotplot(kegg_ora_results)
    enrich_plot

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/Prepare%20ORA%20and%20GSEA-1.png)

    upset_plot <- enrichplot::upsetplot(kegg_ora_results)
    upset_plot

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/Prepare%20ORA%20and%20GSEA-2.png)

\##EXERCISE: alternatives to KEGG?

    go_bp_results <- enrichGO(
      gene         = genes_dif,
      OrgDb        = org.Hs.eg.db,
      keyType      = "ENTREZID",
      ont          = "BP", # Biological Processes
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1,
      universe     = background_set
    )
    dotplot(go_bp_results, showCategory=10) + 
      ggplot2::ggtitle("Top 10 Enriched GO Terms in Biological Process")

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/ORA%20GO-1.png)

## GSEA

    diff_table_all <- topTable(fit2,coef=1,p.value=1,number=nrow(logCPM)) 
    # Step 4: conduct GSEA

    list_ordered <- diff_table_all[,"B"]
    names(list_ordered) <- rownames(diff_table_all)
      
      
    gsea_results <- GSEA(
      geneList = list_ordered, # Ordered ranked gene list
      minGSSize = 25, # Minimum gene set size
      maxGSSize = 500, # Maximum gene set set
      pvalueCutoff = 0.05, # p-value cutoff
      eps = 0, # Boundary for calculating the p value
      seed = TRUE, # Set seed to make results reproducible
      pAdjustMethod = "BH", # Benjamini-Hochberg correction
      TERM2GENE = dplyr::select(
        hs_kegg_df,
        gs_name,
        human_entrez_gene
      )
    )

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (0.01% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : There were 2 pathways for which P-values were not calculated
    ## properly due to unbalanced (positive and negative) gene-level statistic values.
    ## For such pathways pval, padj, NES, log2err are set to NA. You can try to
    ## increase the value of the argument nPermSimple (for example set it nPermSimple
    ## = 10000)

    ## leading edge analysis...

    ## done...

    # Step 5: Visualize / explore
    head(gsea_results@result)

    ##                                                                          ID
    ## KEGG_SPLICEOSOME                                           KEGG_SPLICEOSOME
    ## KEGG_REGULATION_OF_ACTIN_CYTOSKELETON KEGG_REGULATION_OF_ACTIN_CYTOSKELETON
    ## KEGG_INSULIN_SIGNALING_PATHWAY               KEGG_INSULIN_SIGNALING_PATHWAY
    ## KEGG_ENDOCYTOSIS                                           KEGG_ENDOCYTOSIS
    ## KEGG_NEUROTROPHIN_SIGNALING_PATHWAY     KEGG_NEUROTROPHIN_SIGNALING_PATHWAY
    ## KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS     KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS
    ##                                                                 Description
    ## KEGG_SPLICEOSOME                                           KEGG_SPLICEOSOME
    ## KEGG_REGULATION_OF_ACTIN_CYTOSKELETON KEGG_REGULATION_OF_ACTIN_CYTOSKELETON
    ## KEGG_INSULIN_SIGNALING_PATHWAY               KEGG_INSULIN_SIGNALING_PATHWAY
    ## KEGG_ENDOCYTOSIS                                           KEGG_ENDOCYTOSIS
    ## KEGG_NEUROTROPHIN_SIGNALING_PATHWAY     KEGG_NEUROTROPHIN_SIGNALING_PATHWAY
    ## KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS     KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS
    ##                                       setSize enrichmentScore       NES
    ## KEGG_SPLICEOSOME                          126      -0.4612522 -3.250786
    ## KEGG_REGULATION_OF_ACTIN_CYTOSKELETON     147      -0.3387370 -2.453525
    ## KEGG_INSULIN_SIGNALING_PATHWAY            111      -0.3688939 -2.538368
    ## KEGG_ENDOCYTOSIS                          153      -0.3091265 -2.246320
    ## KEGG_NEUROTROPHIN_SIGNALING_PATHWAY       105      -0.3398888 -2.305316
    ## KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS       124      -0.3228639 -2.268650
    ##                                             pvalue     p.adjust       qvalue
    ## KEGG_SPLICEOSOME                      1.535955e-17 1.904585e-15 8.569014e-16
    ## KEGG_REGULATION_OF_ACTIN_CYTOSKELETON 7.312274e-09 4.533610e-07 2.039740e-07
    ## KEGG_INSULIN_SIGNALING_PATHWAY        1.503211e-08 6.213274e-07 2.795446e-07
    ## KEGG_ENDOCYTOSIS                      5.616062e-08 1.740979e-06 7.832929e-07
    ## KEGG_NEUROTROPHIN_SIGNALING_PATHWAY   6.118289e-07 1.215352e-05 5.468052e-06
    ## KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS   5.673356e-07 1.215352e-05 5.468052e-06
    ##                                       rank                   leading_edge
    ## KEGG_SPLICEOSOME                      5073 tags=71%, list=31%, signal=49%
    ## KEGG_REGULATION_OF_ACTIN_CYTOSKELETON 4481 tags=47%, list=27%, signal=34%
    ## KEGG_INSULIN_SIGNALING_PATHWAY        4340 tags=52%, list=26%, signal=39%
    ## KEGG_ENDOCYTOSIS                      4258 tags=46%, list=26%, signal=35%
    ## KEGG_NEUROTROPHIN_SIGNALING_PATHWAY   4253 tags=49%, list=26%, signal=36%
    ## KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS   4868 tags=52%, list=30%, signal=37%
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       core_enrichment
    ## KEGG_SPLICEOSOME                      84321/5093/10929/22916/6625/3310/84991/10286/51691/8559/9410/3304/51503/9939/11157/10465/10569/23451/10189/220988/51690/6432/3303/10285/10450/10713/3192/4116/58517/8896/8683/10992/10946/9416/10772/25804/6631/84844/56259/144983/10262/6635/9092/8449/22985/26121/57461/9128/6427/56949/9343/6633/6628/988/10523/27258/151903/6634/9775/6637/11325/83443/22827/11338/1659/9785/4809/22938/55696/84950/29896/9716/55119/1665/23450/153527/6429/3178/6428/25949/55660/23350/9879/57187/3312/4670/6434/3190/3183
    ## KEGG_REGULATION_OF_ACTIN_CYTOSKELETON                                                                                                                                               10095/3682/5501/3985/71/9459/23365/8396/9138/23191/4893/3984/8874/23533/5296/81873/4478/7430/5605/1398/2934/3845/23396/5829/3688/10093/387/5305/5294/3678/5594/10094/10163/4659/5962/5595/5291/28964/6654/200576/5499/7409/4627/10672/3071/26999/369/9475/7454/89846/3685/8826/929/5500/87/55845/5293/1445/5295/324/10788/1730/10451/1399/998/7114/5894/3683/6093
    ## KEGG_INSULIN_SIGNALING_PATHWAY                                                                                                                                                                                                                  3098/5501/6464/8569/1978/5565/867/5106/4893/10603/5257/5584/5567/5573/3643/23533/5296/5605/1398/3845/2002/5571/7248/5834/6198/6720/6199/2885/51763/5294/5594/122809/57521/5595/5291/6654/7249/2475/5499/9021/5256/3551/6009/369/5566/9470/23265/2932/5500/5293/2872/8660/5170/5295/1399/5836/5894/805
    ## KEGG_ENDOCYTOSIS                                                                                                                                          2869/10254/26119/3304/4734/8396/867/11311/89853/79643/7037/5584/5868/80230/11267/3303/10617/60682/8027/1785/3949/2870/382/5337/23396/409/29082/27183/5338/163/116983/22841/55738/57403/5878/161/84364/51160/26286/28964/200576/7189/1212/6455/160/3480/51652/9230/156/1173/9146/5869/30011/51028/23550/30844/8411/10193/7852/3106/2060/9135/25978/10015/3107/9101/3312/998/9815/128866/3105
    ## KEGG_NEUROTROPHIN_SIGNALING_PATHWAY                                                                                                                                                                                                                                                       10782/6464/3656/25970/596/2309/4893/10603/4793/814/23533/5296/7534/818/5605/1398/3845/7157/11108/7529/581/2885/387/5294/5594/396/5595/3654/5291/8986/6654/7189/25/10818/3551/5970/3725/2932/9261/5293/8660/1445/5170/5295/1399/57498/5663/998/5894/5908/805
    ## KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS                                                                                                                                                  27338/8945/26091/57448/672/1161/22888/9820/51366/92912/51588/4734/7322/867/89910/10075/10393/246184/996/7329/29882/8454/9063/9978/63893/25898/9690/10273/10277/118424/10055/9040/55120/8924/25847/55236/8450/29945/51529/7189/64682/9616/9021/5371/7318/51433/329/8453/331/55284/23291/9354/6923/7323/10477/65264/8554/7320/8065/6500/8451/11060/10054/1642/8452

    gsea_result_df <- data.frame(gsea_results@result)
    gsea_result_df %>%
      # This returns the 3 rows with the largest NES values
      dplyr::slice_max(NES, n = 3)

    ##                                                                                            ID
    ## KEGG_ALZHEIMERS_DISEASE                                               KEGG_ALZHEIMERS_DISEASE
    ## KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY
    ## KEGG_JAK_STAT_SIGNALING_PATHWAY                               KEGG_JAK_STAT_SIGNALING_PATHWAY
    ##                                                                                   Description
    ## KEGG_ALZHEIMERS_DISEASE                                               KEGG_ALZHEIMERS_DISEASE
    ## KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY
    ## KEGG_JAK_STAT_SIGNALING_PATHWAY                               KEGG_JAK_STAT_SIGNALING_PATHWAY
    ##                                                setSize enrichmentScore
    ## KEGG_ALZHEIMERS_DISEASE                            139      -0.2143559
    ## KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY      95      -0.2400268
    ## KEGG_JAK_STAT_SIGNALING_PATHWAY                     93      -0.2424305
    ##                                                      NES      pvalue   p.adjust
    ## KEGG_ALZHEIMERS_DISEASE                        -1.531763 0.009359863 0.02072541
    ## KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY -1.599569 0.007917800 0.01818161
    ## KEGG_JAK_STAT_SIGNALING_PATHWAY                -1.605445 0.008813781 0.01987107
    ##                                                     qvalue rank
    ## KEGG_ALZHEIMERS_DISEASE                        0.009324675 5071
    ## KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY 0.008180183 3133
    ## KEGG_JAK_STAT_SIGNALING_PATHWAY                0.008940294 4238
    ##                                                                  leading_edge
    ## KEGG_ALZHEIMERS_DISEASE                        tags=41%, list=31%, signal=29%
    ## KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY tags=34%, list=19%, signal=27%
    ## KEGG_JAK_STAT_SIGNALING_PATHWAY                tags=39%, list=26%, signal=29%
    ##                                                                                                                                                                                                                                                                                                                           core_enrichment
    ## KEGG_ALZHEIMERS_DISEASE                        7124/23621/572/1329/102/4512/6622/513/23385/4694/4698/7388/1537/4722/1349/4714/4697/11261/4700/55851/27089/4716/317/842/5331/51107/4708/489/9451/4723/4707/5533/5594/4724/1327/4695/5595/517/5330/7386/4731/54205/10975/6389/5530/4725/22926/2932/841/9167/6868/1350/514/351/5663/7381/805
    ## KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY                                                                                                                         11261/3384/23533/5296/5605/8797/3845/8793/3937/2885/3455/5294/5533/5594/5595/5291/6654/5579/10870/7409/2214/5530/369/3106/5777/5293/5295/3107/10451/5894/3683/3105
    ## KEGG_JAK_STAT_SIGNALING_PATHWAY                                                                                                               51588/10254/30837/867/1439/200734/6778/23533/3586/5296/2033/10379/8027/9063/2885/3597/3455/9655/5294/6772/122809/5291/6654/9021/3566/3716/6776/7297/8554/5777/5293/5295/3570/3717/6777/6773

    most_positive_nes_plot <- enrichplot::gseaplot(
      gsea_results,
      geneSetID = "KEGG_JAK_STAT_SIGNALING_PATHWAY",
      title = "KEGG_JAK_STAT_SIGNALING_PATHWAY",
      color.line = "#0d76ff"
    )
    most_positive_nes_plot

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GSEA-1.png)

    gsea_result_df %>%
      # Return the 3 rows with the smallest (most negative) NES values
      dplyr::slice_min(NES, n = 3)

    ##                                                            ID
    ## KEGG_SPLICEOSOME                             KEGG_SPLICEOSOME
    ## KEGG_INSULIN_SIGNALING_PATHWAY KEGG_INSULIN_SIGNALING_PATHWAY
    ## KEGG_RNA_DEGRADATION                     KEGG_RNA_DEGRADATION
    ##                                                   Description setSize
    ## KEGG_SPLICEOSOME                             KEGG_SPLICEOSOME     126
    ## KEGG_INSULIN_SIGNALING_PATHWAY KEGG_INSULIN_SIGNALING_PATHWAY     111
    ## KEGG_RNA_DEGRADATION                     KEGG_RNA_DEGRADATION      55
    ##                                enrichmentScore       NES       pvalue
    ## KEGG_SPLICEOSOME                    -0.4612522 -3.250786 1.535955e-17
    ## KEGG_INSULIN_SIGNALING_PATHWAY      -0.3688939 -2.538368 1.503211e-08
    ## KEGG_RNA_DEGRADATION                -0.4263806 -2.500411 1.630981e-06
    ##                                    p.adjust       qvalue rank
    ## KEGG_SPLICEOSOME               1.904585e-15 8.569014e-16 5073
    ## KEGG_INSULIN_SIGNALING_PATHWAY 6.213274e-07 2.795446e-07 4340
    ## KEGG_RNA_DEGRADATION           2.528021e-05 1.137395e-05 4297
    ##                                                  leading_edge
    ## KEGG_SPLICEOSOME               tags=71%, list=31%, signal=49%
    ## KEGG_INSULIN_SIGNALING_PATHWAY tags=52%, list=26%, signal=39%
    ## KEGG_RNA_DEGRADATION           tags=60%, list=26%, signal=44%
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                core_enrichment
    ## KEGG_SPLICEOSOME               84321/5093/10929/22916/6625/3310/84991/10286/51691/8559/9410/3304/51503/9939/11157/10465/10569/23451/10189/220988/51690/6432/3303/10285/10450/10713/3192/4116/58517/8896/8683/10992/10946/9416/10772/25804/6631/84844/56259/144983/10262/6635/9092/8449/22985/26121/57461/9128/6427/56949/9343/6633/6628/988/10523/27258/151903/6634/9775/6637/11325/83443/22827/11338/1659/9785/4809/22938/55696/84950/29896/9716/55119/1665/23450/153527/6429/3178/6428/25949/55660/23350/9879/57187/3312/4670/6434/3190/3183
    ## KEGG_INSULIN_SIGNALING_PATHWAY                                                                                                                                                                                                           3098/5501/6464/8569/1978/5565/867/5106/4893/10603/5257/5584/5567/5573/3643/23533/5296/5605/1398/3845/2002/5571/7248/5834/6198/6720/6199/2885/51763/5294/5594/122809/57521/5595/5291/6654/7249/2475/5499/9021/5256/3551/6009/369/5566/9470/23265/2932/5500/5293/2872/8660/5170/5295/1399/5836/5894/805
    ## KEGG_RNA_DEGRADATION                                                                                                                                                                                                                                                                                                                               51691/219988/51013/4850/23019/27257/11340/11157/4849/10200/23644/51690/23404/3313/5394/80349/25804/9652/10438/55802/6499/22894/5073/27258/87178/29883/9125/23517/4848/9337/1656/22803/54464

    most_negative_nes_plot <- enrichplot::gseaplot(
      gsea_results,
      geneSetID = "KEGG_SPLICEOSOME",
      title = "KEGG_SPLICEOSOME",
      color.line = "#0d76ff"
    )
    most_negative_nes_plot

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GSEA-2.png)

    # Step 7: EXERCISE: compare GSEA vs ORA?

### GSEA vs ORA

-   ORA is useful when there is a strong hypothesis about specific genes
    or pathways of interest, and the analysis focuses on these
    predefined categories
-   GSEA is useful for exploratory analyses where the goal is to uncover
    potentially novel associations between gene expression patterns and
    biological functions or pathways, especially when the effects are
    subtle

## GeneSetCluster

    # Healthy vs Group Covid19

    # Healthy vs Covid19AI 
    Diff_HvsAI <- topTable(fit2,coef=1,p.value=1,number=nrow(logCPM))
    # Healthy vs Covid196Mo 
    Diff_Hvs6Mo <- topTable(fit2,coef=3,p.value=1,number=nrow(logCPM))


    hs_msigdb_df <- msigdbr(species = "Homo sapiens")
    hs_kegg_df <- hs_msigdb_df %>%
      dplyr::filter(
        gs_cat == "C2", # This is to filter only to the C2 curated gene sets
        gs_subcat == "CP:KEGG" # This is because we only want KEGG pathways
      )


    doGSEA <- function(diff_table) {
      list_ordered <- diff_table[,"B"]
      names(list_ordered) <- rownames(diff_table)
      
      return(GSEA(
        geneList = list_ordered, # Ordered ranked gene list
        minGSSize = 25, # Minimum gene set size
        maxGSSize = 500, # Maximum gene set set
        pvalueCutoff = 0.05, # p-value cutoff
        eps = 0, # Boundary for calculating the p value
        seed = TRUE, # Set seed to make results reproducible
        pAdjustMethod = "BH", # Benjamini-Hochberg correction
        TERM2GENE = dplyr::select(
          hs_kegg_df,
          gs_name,
          human_entrez_gene
        )
      ))
    }

    GSEA_HvsAI <- doGSEA(Diff_HvsAI)

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (0.01% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize =
    ## minSize, : There were 2 pathways for which P-values were not calculated
    ## properly due to unbalanced (positive and negative) gene-level statistic values.
    ## For such pathways pval, padj, NES, log2err are set to NA. You can try to
    ## increase the value of the argument nPermSimple (for example set it nPermSimple
    ## = 10000)

    ## leading edge analysis...

    ## done...

    GSEA_Hvs6Mo <- doGSEA(Diff_Hvs6Mo)

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (0.01% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## leading edge analysis...

    ## done...

    path <- "~/Bioinformatics_Pipelines_Class"

    write.csv(GSEA_HvsAI, file = paste0(path, "/GSEA_HvsAI.csv"), row.names = FALSE)
    write.csv(GSEA_Hvs6Mo, file = paste0(path, "/GSEA_Hvs6Mo.csv"), row.names = FALSE)

    library(GeneSetCluster)

    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

    ## Warning: replacing previous import 'AnnotationDbi::select' by
    ## 'clusterProfiler::select' when loading 'GeneSetCluster'

    ## Warning: replacing previous import 'ComplexHeatmap::%v%' by 'network::%v%' when
    ## loading 'GeneSetCluster'

    ## 

    ## Warning: replacing previous import 'cowplot::align_plots' by
    ## 'patchwork::align_plots' when loading 'GeneSetCluster'

    ## Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
    ## 'GeneSetCluster'

    ## Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when
    ## loading 'GeneSetCluster'

    ## Warning: replacing previous import 'AnnotationDbi::tail' by 'utils::tail' when
    ## loading 'GeneSetCluster'

    ## Warning: replacing previous import 'AnnotationDbi::head' by 'utils::head' when
    ## loading 'GeneSetCluster'

    GSEA.files <- paste0(path, "/", list.files(path, pattern = ".csv"))

    # Load the data and create Pathway object
    # Automatically for GSEA, GREAT or IPA
    GSEA.Object1 <- LoadGeneSets(file_location = GSEA.files, 
                                  groupnames= c("GSEA_Hvs6Mo", "GSEA_HvsAI"), # names of the groups
                                  P.cutoff = 0.05, # cut off the p.adjust
                                  Mol.cutoff = 15, # minimum number of genes per pathway
                                  Source = "GSEA", # the analysis (GSEA, GREAT or IPA)
                                  structure = "ENTREZID", # Gene type (SYMBOL, ENTREZID, ENSEMBLID)
                                  Organism = "org.Hs.eg.db", # database: Homo Sapiens or Mus musculus
                                  seperator = "/") # the separator used for listing genes

    ## [=========================================================]

    ## [<<<<            LoadGeneSets START                  >>>>>]

    ## -----------------------------------------------------------

    ## Loading data from ~/Bioinformatics_Pipelines_Class/GSEA_Hvs6Mo.csv

    ## Loading data from ~/Bioinformatics_Pipelines_Class/GSEA_HvsAI.csv

    ## [=========================================================]

    ## [<<<<            ObjectCreator START                 >>>>>]

    ## -----------------------------------------------------------
    ## -----------------------------------------------------------

    ## [<<<<<               ObjectCreator END              >>>>>>]

    ## [=========================================================]

    ## [You may want to process CombineGeneSets next.            ]

    ## [or merge objects using MergeObjects                     ]

    ## [or select certain types from objects using ManageGeneSets]

    ## -----------------------------------------------------------

    ## [<<<<<               LoadGeneSets END               >>>>>>]

    ## [=========================================================]

    ## [You may want to process CombineGeneSets next.            ]

    ## [or merge objects using MergeObjects                     ]

    ## [or select certain types from objects using ManageGeneSets]

    # IMPORTANT when created manually, it is assumed that the pathways have been filtered by p-value and minimum number of genes per pathway
    # Make sure you have filtered your data
    GSEA.Object1Manual <- ObjectCreator(Pathways = c(GSEA_HvsAI@result$ID, 
                                                     GSEA_Hvs6Mo@result$ID),
                                        Molecules = c(GSEA_HvsAI@result$core_enrichment, 
                                                      GSEA_Hvs6Mo@result$core_enrichment),
                                        Groups = c(rep("GSEA_HvsAI", times=nrow(GSEA_HvsAI@result)), 
                                                   rep("GSEA_Hvs6Mo", times=nrow(GSEA_Hvs6Mo@result))),
                                        Pvalues = c(GSEA_HvsAI@result$p.adjust,  # optional
                                                    GSEA_Hvs6Mo@result$p.adjust),
                                        enrichmentScore = c(GSEA_HvsAI@result$NES, # optional
                                                            GSEA_Hvs6Mo@result$NES), 
                                        structure = "ENTREZID", Type = "", sep = "/",
                                        Source = "GSEA", organism = "org.Hs.eg.db")

    ## [=========================================================]

    ## [<<<<            ObjectCreator START                 >>>>>]

    ## -----------------------------------------------------------
    ## -----------------------------------------------------------

    ## [<<<<<               ObjectCreator END              >>>>>>]

    ## [=========================================================]

    ## [You may want to process CombineGeneSets next.            ]

    ## [or merge objects using MergeObjects                     ]

    ## [or select certain types from objects using ManageGeneSets]

    GSEA.Object2 <- CombineGeneSets(Object = GSEA.Object1,
                                    combineMethod = "Standard", threads = 8)

    ## [=========================================================]

    ## [<<<<            CombineGeneSets START               >>>>>]

    ## -----------------------------------------------------------

    ## Combining experiments

    ## preparing condensed display

    ## calulating RR

    ## Parallelizing processes, 8 cores have been selected.

    ## -----------------------------------------------------------

    ## [<<<<<             CombineGeneSets END              >>>>>>]

    ## [=========================================================]

    ## [You may want to process ClusterGeneSets next.            ]

    OptimalGeneSets(Object = GSEA.Object2, 
                    uniquePathway = FALSE, # consider all the pathways (also repeated) or the unique pathways
                    method = "silhouette", max_cluster= 24, cluster_method = "kmeans", main= "Kmeans for 24 clusters")

    ## [=========================================================]

    ## [<<<<            OptimalGeneSets START               >>>>>]

    ## -----------------------------------------------------------

    ## Finding optimal number of clusters

    ## Clustering method =  kmeans

    ## Optimizing method =  silhouette

    ## -----------------------------------------------------------

    ## [<<<<<             OptimalGeneSets END              >>>>>>]

    ## [=========================================================]

    ## [You may want to process ClusterGeneSets next.            ]

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-1.png)

    OptimalGeneSets(Object = GSEA.Object2, 
                    uniquePathway = TRUE, # consider all the pathways (also repeated) or the unique pathways
                    method = "silhouette", max_cluster= 24, cluster_method = "kmeans", main= "Kmeans for 24 clusters")

    ## [=========================================================]

    ## [<<<<            OptimalGeneSets START               >>>>>]

    ## -----------------------------------------------------------

    ## Finding optimal number of clusters

    ## Clustering method =  kmeans

    ## Optimizing method =  silhouette

    ## -----------------------------------------------------------

    ## [<<<<<             OptimalGeneSets END              >>>>>>]

    ## [=========================================================]

    ## [You may want to process ClusterGeneSets next.            ]

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-2.png)

    # in both cases the optimal cluster is 2

    GSEA.Object3 <- ClusterGeneSets(Object = GSEA.Object2, 
                                    clusters = 2, # consider all the pathways (also repeated) or the unique pathways
                                    method = "Hierarchical", # Hierarchical clustering or kmeans
                                    order = "cluster",
                                    molecular.signature = "All")

    ## [=========================================================]

    ## [<<<<            ClusterGeneSets START               >>>>>]

    ## -----------------------------------------------------------

    ## Using all Gene sets

    ## Running Hierarchical

    ## Ordering pathway clusters

    ## -----------------------------------------------------------

    ## [<<<<<             ClusterGeneSets END              >>>>>>]

    ## [=========================================================]

    ## [You may want to process HighlightGeneSets next.            ]

    ## [You may want to plot the results using PlotGeneSets next. ]

    # plot results for both all pathways and unique pathways
    plotnounique <- PlotGeneSets(GSEA.Object3, 
                                 uniquePathways = FALSE, 
                                 wordcloud = FALSE, # wordcloud only supported for GO terms
                                 doORA = T) # do ora per cluster

    ## Performing ORA for 2 clusters and choosing the most overrepresented pathway...

    ## No legend element is put in the last 6 rows under `nrow = 8`, maybe you
    ## should set `by_row = FALSE`? Reset `nrow` to 2.

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-3.png)

    plotunique <- PlotGeneSets(GSEA.Object3, 
                               uniquePathways = TRUE, 
                               wordcloud = FALSE, # wordcloud only supported for GO terms
                               doORA = T) # do ora per cluster

    ## Performing ORA for 2 clusters and choosing the most overrepresented pathway...
    ## No legend element is put in the last 6 rows under `nrow = 8`, maybe you
    ## should set `by_row = FALSE`? Reset `nrow` to 2.

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-4.png)

    # let's say we are interested in exploring cluster 2 in plotunique. Lets break up this cluste for further analysis 

    plotoptimalcluster2 <- OptimalGeneSets(Object = GSEA.Object3, 
                    uniquePathway = TRUE, # consider all the pathways (also repeated) or the unique pathways
                    cluster = 2, # which cluster
                    method = "silhouette", max_cluster= 24, cluster_method = "kmeans", main= "Kmeans for 24 clusters in cluster 1")

    ## [=========================================================]

    ## [<<<<            OptimalGeneSets START               >>>>>]

    ## -----------------------------------------------------------

    ## Finding optimal number of clusters

    ## Clustering method =  kmeans

    ## Optimizing method =  silhouette

    ## -----------------------------------------------------------

    ## [<<<<<             OptimalGeneSets END              >>>>>>]

    ## [=========================================================]

    ## [You may want to process ClusterGeneSets next.            ]

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-5.png)

    plotoptimalcluster2 # optimal 2 break up cluster 2 in 2 clusters

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-6.png)

    GSEA.Object3breakup <- BreakUpCluster(GSEA.Object3, 
                                          breakup.cluster = 2, # which cluster
                                          sub.cluster = 2, # in how many cluster split up
                                          uniquePathways = TRUE) # conside unique pathways

    ## [=========================================================]

    ## [<<<<            BreakUpCluster START                >>>>>]

    ## -----------------------------------------------------------

    ## Loading required package: RColorBrewer

    ## Warning in BreakUpCluster(GSEA.Object3, breakup.cluster = 2, sub.cluster = 2, :
    ## Currently only working for Kmeans and Hierarchical

    ## Ordering pathway clusters

    ## -----------------------------------------------------------

    ## [<<<<<             BreakUpCluster END               >>>>>>]

    ## [=========================================================]

    ## [You may want to process HighlightGeneSets next.           ]

    ## [You may want to plot the results using PlotGeneSets next. ]

    plotuniquebreakup <- PlotGeneSets(GSEA.Object3breakup, 
                                      uniquePathways = TRUE, 
                                      wordcloud = FALSE, # wordcloud only supported for GO terms
                                      doORA = T) # do ora per cluster

    ## Performing ORA for 3 clusters and choosing the most overrepresented pathway...

    ## No legend element is put in the last 5 rows under `nrow = 8`, maybe you
    ## should set `by_row = FALSE`? Reset `nrow` to 3.

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-7.png)

    plotuniquebreakup

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-8.png)

    # Now break up the cluster 1 
    plotoptimalcluster1 <- OptimalGeneSets(Object = GSEA.Object3, 
                    uniquePathway = TRUE, # consider all the pathways (also repeated) or the unique pathways
                    cluster = 1, # which cluster
                    method = "silhouette", max_cluster= 24, cluster_method = "kmeans", main= "Kmeans for 24 clusters in cluster 1")

    ## [=========================================================]

    ## [<<<<            OptimalGeneSets START               >>>>>]

    ## -----------------------------------------------------------

    ## Finding optimal number of clusters

    ## Clustering method =  kmeans

    ## Optimizing method =  silhouette

    ## -----------------------------------------------------------

    ## [<<<<<             OptimalGeneSets END              >>>>>>]

    ## [=========================================================]

    ## [You may want to process ClusterGeneSets next.            ]

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-9.png)

    plotoptimalcluster1 # optimal 1 break up cluster 1 in 9 clusters

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-10.png)

    GSEA.Object3breakup2 <- BreakUpCluster(GSEA.Object3breakup, 
                                          breakup.cluster = 1, # which cluster
                                          sub.cluster = 9, # in how many cluster split up
                                          uniquePathways = TRUE) # conside unique pathways

    ## [=========================================================]

    ## [<<<<            BreakUpCluster START                >>>>>]

    ## -----------------------------------------------------------

    ## Warning in BreakUpCluster(GSEA.Object3breakup, breakup.cluster = 1, sub.cluster
    ## = 9, : Currently only working for Kmeans and Hierarchical

    ## Ordering pathway clusters
    ## -----------------------------------------------------------

    ## [<<<<<             BreakUpCluster END               >>>>>>]

    ## [=========================================================]

    ## [You may want to process HighlightGeneSets next.           ]

    ## [You may want to plot the results using PlotGeneSets next. ]

    plotuniquebreakup2 <- PlotGeneSets(GSEA.Object3breakup2, 
                                       uniquePathways = TRUE, 
                                       wordcloud = FALSE, # wordcloud only supported for GO terms
                                       doORA = T) # do ora per cluster

    ## Performing ORA for 11 clusters and choosing the most overrepresented pathway...

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-11.png)

    plotuniquebreakup2

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/run%20GeneSetCluster-12.png)

    # plot results for both all pathways and unique pathways
    plotnounique <- PlotGeneSets(GSEA.Object3, 
                                 uniquePathways = FALSE, 
                                 wordcloud = FALSE, # wordcloud only supported for GO terms
                                 doORA = T) # do ora per cluster

    ## Performing ORA for 2 clusters and choosing the most overrepresented pathway...

    ## No legend element is put in the last 6 rows under `nrow = 8`, maybe you
    ## should set `by_row = FALSE`? Reset `nrow` to 2.

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/explore%20results-1.png)

    plotunique <- PlotGeneSets(GSEA.Object3, 
                               uniquePathways = TRUE, 
                               wordcloud = FALSE, # wordcloud only supported for GO terms
                               doORA = T) # do ora per cluster

    ## Performing ORA for 2 clusters and choosing the most overrepresented pathway...
    ## No legend element is put in the last 6 rows under `nrow = 8`, maybe you
    ## should set `by_row = FALSE`? Reset `nrow` to 2.

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/explore%20results-2.png)

    # let's say we are interested in exploring cluster 2 in plotunique. Lets break up this cluste for further analysis 

    plotoptimalcluster2 <- OptimalGeneSets(Object = GSEA.Object3, 
                    uniquePathway = TRUE, # consider all the pathways (also repeated) or the unique pathways
                    cluster = 2, # which cluster
                    method = "silhouette", max_cluster= 24, cluster_method = "kmeans", main= "Kmeans for 24 clusters in cluster 1")

    ## [=========================================================]

    ## [<<<<            OptimalGeneSets START               >>>>>]

    ## -----------------------------------------------------------

    ## Finding optimal number of clusters

    ## Clustering method =  kmeans

    ## Optimizing method =  silhouette

    ## -----------------------------------------------------------

    ## [<<<<<             OptimalGeneSets END              >>>>>>]

    ## [=========================================================]

    ## [You may want to process ClusterGeneSets next.            ]

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/explore%20results-3.png)

    plotoptimalcluster2 # optimal 2 break up cluster 2 in 2 clusters

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/explore%20results-4.png)

    GSEA.Object3breakup <- BreakUpCluster(GSEA.Object3, 
                                          breakup.cluster = 2, # which cluster
                                          sub.cluster = 2, # in how many cluster split up
                                          uniquePathways = TRUE) # conside unique pathways

    ## [=========================================================]

    ## [<<<<            BreakUpCluster START                >>>>>]

    ## -----------------------------------------------------------

    ## Warning in BreakUpCluster(GSEA.Object3, breakup.cluster = 2, sub.cluster = 2, :
    ## Currently only working for Kmeans and Hierarchical

    ## Ordering pathway clusters
    ## -----------------------------------------------------------

    ## [<<<<<             BreakUpCluster END               >>>>>>]

    ## [=========================================================]

    ## [You may want to process HighlightGeneSets next.           ]

    ## [You may want to plot the results using PlotGeneSets next. ]

    plotuniquebreakup <- PlotGeneSets(GSEA.Object3breakup, 
                                      uniquePathways = TRUE, 
                                      wordcloud = FALSE, # wordcloud only supported for GO terms
                                      doORA = T) # do ora per cluster

    ## Performing ORA for 3 clusters and choosing the most overrepresented pathway...

    ## No legend element is put in the last 5 rows under `nrow = 8`, maybe you
    ## should set `by_row = FALSE`? Reset `nrow` to 3.

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/explore%20results-5.png)

    plotuniquebreakup

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/explore%20results-6.png)

    # Now break up the cluster 1 
    plotoptimalcluster1 <- OptimalGeneSets(Object = GSEA.Object3, 
                    uniquePathway = TRUE, # consider all the pathways (also repeated) or the unique pathways
                    cluster = 1, # which cluster
                    method = "silhouette", max_cluster= 24, cluster_method = "kmeans", main= "Kmeans for 24 clusters in cluster 1")

    ## [=========================================================]

    ## [<<<<            OptimalGeneSets START               >>>>>]

    ## -----------------------------------------------------------

    ## Finding optimal number of clusters

    ## Clustering method =  kmeans

    ## Optimizing method =  silhouette

    ## -----------------------------------------------------------

    ## [<<<<<             OptimalGeneSets END              >>>>>>]

    ## [=========================================================]

    ## [You may want to process ClusterGeneSets next.            ]

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/explore%20results-7.png)

    plotoptimalcluster1 # optimal 1 break up cluster 1 in 9 clusters

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/explore%20results-8.png)

    GSEA.Object3breakup2 <- BreakUpCluster(GSEA.Object3breakup, 
                                          breakup.cluster = 1, # which cluster
                                          sub.cluster = 9, # in how many cluster split up
                                          uniquePathways = TRUE) # conside unique pathways

    ## [=========================================================]

    ## [<<<<            BreakUpCluster START                >>>>>]

    ## -----------------------------------------------------------

    ## Warning in BreakUpCluster(GSEA.Object3breakup, breakup.cluster = 1, sub.cluster
    ## = 9, : Currently only working for Kmeans and Hierarchical

    ## Ordering pathway clusters
    ## -----------------------------------------------------------

    ## [<<<<<             BreakUpCluster END               >>>>>>]

    ## [=========================================================]

    ## [You may want to process HighlightGeneSets next.           ]

    ## [You may want to plot the results using PlotGeneSets next. ]

    plotuniquebreakup2 <- PlotGeneSets(GSEA.Object3breakup2, 
                                       uniquePathways = TRUE, 
                                       wordcloud = FALSE, # wordcloud only supported for GO terms
                                       doORA = T) # do ora per cluster

    ## Performing ORA for 11 clusters and choosing the most overrepresented pathway...

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/explore%20results-9.png)

    plotuniquebreakup2

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/explore%20results-10.png)

    # let's say we are interested in exploring cluster 2 in plotunique. Lets break up this cluste for further analysis 

    plotoptimalcluster2 <- OptimalGeneSets(Object = GSEA.Object3, 
                    uniquePathway = TRUE, # consider all the pathways (also repeated) or the unique pathways
                    cluster = 2, # which cluster
                    method = "silhouette", max_cluster= 24, cluster_method = "kmeans", main= "Kmeans for 24 clusters in cluster 1")

    ## [=========================================================]

    ## [<<<<            OptimalGeneSets START               >>>>>]

    ## -----------------------------------------------------------

    ## Finding optimal number of clusters

    ## Clustering method =  kmeans

    ## Optimizing method =  silhouette

    ## -----------------------------------------------------------

    ## [<<<<<             OptimalGeneSets END              >>>>>>]

    ## [=========================================================]

    ## [You may want to process ClusterGeneSets next.            ]

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/breakup%20cluster%20number%201-1.png)

    plotoptimalcluster2 # optimal 2 break up cluster 2 in 2 clusters

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/breakup%20cluster%20number%201-2.png)

    GSEA.Object3breakup <- BreakUpCluster(GSEA.Object3, 
                                          breakup.cluster = 2, # which cluster
                                          sub.cluster = 2, # in how many cluster split up
                                          uniquePathways = TRUE) # conside unique pathways

    ## [=========================================================]

    ## [<<<<            BreakUpCluster START                >>>>>]

    ## -----------------------------------------------------------

    ## Warning in BreakUpCluster(GSEA.Object3, breakup.cluster = 2, sub.cluster = 2, :
    ## Currently only working for Kmeans and Hierarchical

    ## Ordering pathway clusters
    ## -----------------------------------------------------------

    ## [<<<<<             BreakUpCluster END               >>>>>>]

    ## [=========================================================]

    ## [You may want to process HighlightGeneSets next.           ]

    ## [You may want to plot the results using PlotGeneSets next. ]

    plotuniquebreakup <- PlotGeneSets(GSEA.Object3breakup, 
                                      uniquePathways = TRUE, 
                                      wordcloud = FALSE, # wordcloud only supported for GO terms
                                      doORA = T) # do ora per cluster

    ## Performing ORA for 3 clusters and choosing the most overrepresented pathway...

    ## No legend element is put in the last 5 rows under `nrow = 8`, maybe you
    ## should set `by_row = FALSE`? Reset `nrow` to 3.

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/breakup%20cluster%20number%201-3.png)

    plotuniquebreakup

![](GSE198256_ANALYSIS_DAY2_GENESETCLUSTER_files/figure-markdown_strict/breakup%20cluster%20number%201-4.png)
