## Introduction

This document presents an analysis of RNA-Seq data from the GEO
database, focusing on differential gene expression in the context of
COVID-19. The analysis includes data preprocessing, gene annotation,
quality control, normalization, and differential expression analysis.

## Load libraries

    library(R.utils)
    library(DESeq2)
    library(GEOquery)
    library(data.table)
    library(NOISeq)
    library(biomaRt)
    library("pheatmap")

## Data and metadata retrieval

Download and load the RNA-Seq count data from the GEO database.

    # Load data
    mydata <- read.table('My_data.tsv', sep = "\t", header = TRUE, row.names = 1)
    # Get metadata
    metadata <- read.table('metadata_mydata.txt', sep = "\t", header = TRUE, row.names = 1)
    Factors <- metadata[,c("infection.ch1")]

## Map genes using biomaRt

    # Map gene names using biomaRt
    ensemblIDs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    genes <- as.character(rownames(mydata))
    attributes <- c("percentage_gene_gc_content", "start_position", "end_position", "chromosome_name", "gene_biotype", "hgnc_symbol")
    annotgene <- getBM(attributes = attributes, 
                       filters = "hgnc_symbol", 
                       values = genes, 
                       mart = ensemblIDs)

    ## Batch submitting query [==>----------------------------] 11% eta: 40sBatch
    ## submitting query [======>------------------------] 22% eta: 25sBatch submitting
    ## query [=========>---------------------] 33% eta: 17sBatch submitting query
    ## [=============>-----------------] 44% eta: 12sBatch submitting query
    ## [================>--------------] 56% eta: 9sBatch submitting query
    ## [====================>----------] 67% eta: 7sBatch submitting query
    ## [=======================>-------] 78% eta: 5sBatch submitting query
    ## [===========================>---] 89% eta: 2s

    annotgene <- annotgene[annotgene$chromosome_name %in% c(as.character(1:22) ,"X","Y"),] # Only keeping genes in these chromosomes
    annotgene_filt <- annotgene[!duplicated(annotgene$hgnc_symbol),] # Removing duplicates
    annotgene_filt <- annotgene_filt[!is.na(annotgene_filt$hgnc_symbol), ]
    rownames(annotgene_filt) <- as.character(annotgene_filt$hgnc_symbol) # Checking our genes + annotation overlap
    sum(as.character(rownames(annotgene_filt)) %in% rownames(mydata))

    ## [1] 28381

## NOISeq exploration

    # Only keeping genes that overlap in both
    mydata_count_filt <- mydata[rownames(mydata) %in% rownames(annotgene_filt),]
    annotgene_ord <- annotgene_filt[rownames(mydata_count_filt),]

    # NOISeq
    Factors <- metadata[rownames(metadata) %in% colnames(mydata_count_filt), "infection.ch1", drop = FALSE]
    Factors <- Factors[match(colnames(mydata_count_filt), rownames(Factors)), , drop = FALSE]
    colnames(Factors)[1]<- "Group"

    lengthuse <- abs(annotgene_ord$end-annotgene_ord$start)
    names(lengthuse) <- rownames(annotgene_ord)
    gc <- annotgene_ord$percentage_gene_gc_content # GC%
    names(gc) <- rownames(annotgene_ord)
    biotype <-annotgene_ord$gene_biotype
    names(biotype) <- rownames(annotgene_ord)

    # Creating the object
    data_NOISEQ <- readData(data = mydata_count_filt,
                            length=lengthuse,
                            gc= gc,
                            biotype= biotype,
                            chromosome = annotgene_ord[,c("chromosome_name","start_position","end_position")],
                            factors = Factors) # After fixing

    # Check representation of RNA biotypes in your sample compares to their distribution in the genome
    myexplodata <- dat(data_NOISEQ, type = "biodetection")

    ## Biotypes detection is to be computed for:
    ##  [1] "SRR26705739" "SRR26705741" "SRR26705742" "SRR26705743" "SRR26705744"
    ##  [6] "SRR26705745" "SRR26705746" "SRR26705747" "SRR26705748" "SRR26705749"
    ## [11] "SRR26705750"

    explo.plot(myexplodata, plottype = "persample")

![](Rund_ExtraMile_files/figure-markdown_strict/unnamed-chunk-4-1.png)![](Rund_ExtraMile_files/figure-markdown_strict/unnamed-chunk-4-2.png)

    # Check sequencing depth and feature distribution
    mysaturation = dat(data_NOISEQ, k = 0, ndepth = 7, type = "saturation")
    explo.plot(mysaturation, toplot = 1, samples = 1:11, yleftlim = NULL, yrightlim = NULL)

![](Rund_ExtraMile_files/figure-markdown_strict/unnamed-chunk-4-3.png)

    # Check GC content bias
    myGCbias = dat(data_NOISEQ, factor = "Group", type = "GCbias")

    ## [1] "Warning: 9571 features with 0 counts in all samples are to be removed for this analysis."
    ## [1] "GC content bias detection is to be computed for:"
    ## [1] "COVID19" "Control"
    ## [1] "COVID19"
    ## 
    ## Call:
    ## lm(formula = datos[, i] ~ bx)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -12.5916  -4.1914  -0.0387   3.2244  18.7272 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    7.859      6.267   1.254 0.213431    
    ## bx1           10.151      8.863   1.145 0.255453    
    ## bx2          -21.932     29.936  -0.733 0.465870    
    ## bx3           38.022     10.612   3.583 0.000575 ***
    ## bx4           27.511      7.580   3.630 0.000493 ***
    ## bx5           31.414      7.250   4.333 4.14e-05 ***
    ## bx6           11.696      7.353   1.591 0.115508    
    ## bx7           15.299      7.536   2.030 0.045587 *  
    ## bx8            4.468      8.172   0.547 0.586010    
    ## bx9            6.592     10.325   0.638 0.525007    
    ## bx10          19.907     29.187   0.682 0.497128    
    ## bx11        -155.965    164.343  -0.949 0.345399    
    ## bx12              NA         NA      NA       NA    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.267 on 82 degrees of freedom
    ## Multiple R-squared:  0.7157, Adjusted R-squared:  0.6775 
    ## F-statistic: 18.76 on 11 and 82 DF,  p-value: < 2.2e-16
    ## 
    ## [1] "Control"
    ## 
    ## Call:
    ## lm(formula = datos[, i] ~ bx)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -12.4358  -4.2112  -0.6428   2.4541  20.0820 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    6.726      6.578   1.022 0.309564    
    ## bx1           13.287      9.303   1.428 0.157002    
    ## bx2          -20.438     31.420  -0.650 0.517196    
    ## bx3           38.154     11.138   3.426 0.000961 ***
    ## bx4           28.439      7.955   3.575 0.000591 ***
    ## bx5           29.833      7.609   3.921 0.000182 ***
    ## bx6           12.920      7.717   1.674 0.097902 .  
    ## bx7           13.561      7.909   1.715 0.090204 .  
    ## bx8            5.630      8.577   0.656 0.513389    
    ## bx9            6.249     10.837   0.577 0.565776    
    ## bx10          21.874     30.634   0.714 0.477217    
    ## bx11        -162.677    172.488  -0.943 0.348391    
    ## bx12              NA         NA      NA       NA    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.578 on 82 degrees of freedom
    ## Multiple R-squared:  0.6954, Adjusted R-squared:  0.6546 
    ## F-statistic: 17.02 on 11 and 82 DF,  p-value: < 2.2e-16

    explo.plot(myGCbias, samples = NULL, toplot = "global")

![](Rund_ExtraMile_files/figure-markdown_strict/unnamed-chunk-4-4.png)

    # Count depth test
    mycd = dat(data_NOISEQ, type = "cd", norm = FALSE, refColumn = 1)

    ## [1] "Warning: 9571 features with 0 counts in all samples are to be removed for this analysis."
    ## [1] "Reference sample is: SRR26705739"
    ## [1] "Confidence intervals for median of M:"
    ##             0.25%                 99.75%                Diagnostic Test
    ## SRR26705741 "0.0409692555279669"  "0.0409692555279669"  "FAILED"       
    ## SRR26705742 "-0.105007230149952"  "-0.0705861695957229" "FAILED"       
    ## SRR26705743 "0.725748785849823"   "0.736357962133176"   "FAILED"       
    ## SRR26705744 "0.672048944936015"   "0.768360134568675"   "FAILED"       
    ## SRR26705745 "0.161293353706316"   "0.280337508632329"   "FAILED"       
    ## SRR26705746 "0.555096334495822"   "0.701540611597368"   "FAILED"       
    ## SRR26705747 "-0.0735146089877083" "-0.0384522234268172" "FAILED"       
    ## SRR26705748 "0.294722660786065"   "0.374608314508942"   "FAILED"       
    ## SRR26705749 "0.0640860752958727"  "0.0640860752958727"  "FAILED"       
    ## SRR26705750 "0.19535872398136"    "0.19535872398136"    "FAILED"       
    ## [1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

    explo.plot(mycd,samples = 1:11)

![](Rund_ExtraMile_files/figure-markdown_strict/unnamed-chunk-4-5.png)

Comparing my results obtained from fastq files to directly downloaded
count tables from GEO:

1.  Sequencing depth variation: my samples (42-55%) vs their samples
    (60-69%)
2.  Length bias in gene expression: length bias present in both
3.  GC bias in gene expression: statistically significant in both
4.  Count depth: diagnostic test failed in both, recommending
    normalization

## Differential gene expression analysis

    # Clean up metadata
    pDataUSE <- pData(data_NOISEQ)
    levels(pDataUSE$Group) <- c("Control", "COVID19")
    pDataUSE[,1] <- as.factor(pDataUSE[,1])

    # Create DESeq2 object
    mydata_DESeq2 <- DESeqDataSetFromMatrix(countData = mydata_count_filt,
                                               colData = pDataUSE,
                                               design = ~ Group)

    ##   it appears that the last variable in the design formula, 'Group',
    ##   has a factor level, 'Control', which is not the reference level. we recommend
    ##   to use factor(...,levels=...) or relevel() to set this as the reference level
    ##   before proceeding. for more information, please see the 'Note on factor levels'
    ##   in vignette('DESeq2').

    keep <- rowSums(counts(mydata_DESeq2) >= 10) >= 6 # Filter out genes that do not have at least 10 counts in at least 6 samples
    mydata_DESeq2_F <- mydata_DESeq2[keep,]


    # Run DESeq2
    mydata_DESeq2_F<- DESeq(mydata_DESeq2_F)

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 22 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    mydata_res <- results(mydata_DESeq2_F)

    # Filter DEGs based on criteria
    significant_DEGs <- mydata_res$padj < 0.05 & abs(mydata_res$log2FoldChange) > 1
    num_significant_DEGs <- sum(significant_DEGs, na.rm = TRUE)
    upregulated_genes <- mydata_res$log2FoldChange > 1 & mydata_res$padj < 0.05
    downregulated_genes <- mydata_res$log2FoldChange < -1 & mydata_res$padj < 0.05
    num_upregulated_genes <- sum(upregulated_genes, na.rm = TRUE)
    num_downregulated_genes <- sum(downregulated_genes, na.rm = TRUE)

DEGs: (Comparison between my results obtained from fastq files to
directly downloaded count tables from GEO)

1.  Total significant DEGs: 121 vs 477
2.  Upregulated genes: 58 vs 168
3.  Downregulated genes: 63 vs 309

## Heatmap

Expression of the top 20 DEGs in the dataset

    vsd <- vst(mydata_DESeq2_F, blind=FALSE)
    select <- order(rowMeans(counts(mydata_DESeq2_F,normalized=TRUE)),
                    decreasing=TRUE)[1:11]
    df <- as.data.frame(colData(mydata_DESeq2_F)[,c("Group")])
    colnames(df) <- "Group"
    pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
             cluster_cols=FALSE)

![](Rund_ExtraMile_files/figure-markdown_strict/unnamed-chunk-6-1.png)
