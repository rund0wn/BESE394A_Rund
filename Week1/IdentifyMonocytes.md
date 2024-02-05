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

    urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
    path <- paste(urld, "acc=GSE247186", "file=GSE247186_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
    GSE247186_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
    gds <- getGEO("GSE247186")
    Meta_GSE247186 <- pData(gds$GSE247186_series_matrix.txt.gz@phenoData)
    Factors_GSE247186 <- Meta_GSE247186[,c("infection:ch1")]

## Map genes using biomaRt

    ensemblIDs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    genes <- as.character(rownames(GSE247186_count))
    attributes <- c("percentage_gene_gc_content", "start_position", "end_position", "chromosome_name", "gene_biotype", "entrezgene_id")
    annotgene <- getBM(attributes = attributes, 
                       filters = "entrezgene_id", 
                       values = genes, 
                       mart = ensemblIDs)
    annotgene <- annotgene[annotgene$chromosome_name %in% c(as.character(1:22) ,"X","Y"),] # Only keeping genes in these chromosomes
    annotgene_filt <- annotgene[!duplicated(annotgene$entrezgene_id),] # Removing duplicates
    rownames(annotgene_filt) <- as.character(annotgene_filt$entrezgene_id)
    GSE247186_count_filt <- GSE247186_count[rownames(GSE247186_count) %in% rownames(annotgene_filt),]
    annotgene_ord <- annotgene_filt[rownames(GSE247186_count_filt),]

## NOISeq exploration

    # Setting up variables we need
    Factors_GSE247186 <- data.frame(Meta_GSE247186 [ colnames(GSE247186_count_filt),c("infection:ch1")])
    colnames(Factors_GSE247186)[1]<- "Group"
    lengthuse <- abs(annotgene_ord$end-annotgene_ord$start)
    names(lengthuse) <- rownames(annotgene_ord)
    gc <- annotgene_ord$percentage_gene_gc_content # GC%
    names(gc) <- rownames(annotgene_ord)
    biotype <-annotgene_ord$gene_biotype
    names(biotype) <- rownames(annotgene_ord)

    # Creating the object
    data_NOISEQ <- readData(data = GSE247186_count_filt,
                            length=lengthuse,
                            gc= gc,
                            biotype= biotype,
                            chromosome = annotgene_ord[,c("chromosome_name","start_position","end_position")],
                            factors = Factors_GSE247186)

    # Check representation of RNA biotypes in your sample compares to their distribution in the genome
    myexplodata <- dat(data_NOISEQ, type = "biodetection")

    ## Biotypes detection is to be computed for:
    ##  [1] "GSM7884936" "GSM7884937" "GSM7884938" "GSM7884939" "GSM7884940"
    ##  [6] "GSM7884942" "GSM7884943" "GSM7884944" "GSM7884945" "GSM7884946"
    ## [11] "GSM7884947"

    explo.plot(myexplodata, plottype = "persample") 

![](RNA_Seq_Data_Analysis_files/figure-markdown_strict/unnamed-chunk-4-1.png)![](RNA_Seq_Data_Analysis_files/figure-markdown_strict/unnamed-chunk-4-2.png)

    # Check sequencing depth and feature distribution
    mysaturation = dat(data_NOISEQ, k = 0, ndepth = 7, type = "saturation")
    explo.plot(mysaturation, toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)

![](RNA_Seq_Data_Analysis_files/figure-markdown_strict/unnamed-chunk-4-3.png)

    # Check GC content bias
    myGCbias = dat(data_NOISEQ, factor = "Group", type = "GCbias")

    ## [1] "Warning: 7876 features with 0 counts in all samples are to be removed for this analysis."
    ## [1] "GC content bias detection is to be computed for:"
    ## [1] "COVID-19" "control" 
    ## [1] "COVID-19"
    ## 
    ## Call:
    ## lm(formula = datos[, i] ~ bx)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -14.459  -2.913   0.000   2.671  15.080 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    7.827      6.026   1.299 0.197864    
    ## bx1           12.681      8.522   1.488 0.140818    
    ## bx2         -121.716     99.480  -1.224 0.224865    
    ## bx3           40.614     14.737   2.756 0.007302 ** 
    ## bx4           27.424      7.799   3.516 0.000737 ***
    ## bx5           34.844      7.071   4.928 4.65e-06 ***
    ## bx6           18.320      7.095   2.582 0.011722 *  
    ## bx7           15.219      7.249   2.099 0.039053 *  
    ## bx8           10.486      7.822   1.341 0.183985    
    ## bx9            2.112      9.801   0.215 0.829960    
    ## bx10          31.660     27.614   1.147 0.255136    
    ## bx11        -183.953    152.206  -1.209 0.230522    
    ## bx12              NA         NA      NA       NA    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.026 on 77 degrees of freedom
    ## Multiple R-squared:  0.7263, Adjusted R-squared:  0.6872 
    ## F-statistic: 18.57 on 11 and 77 DF,  p-value: < 2.2e-16
    ## 
    ## [1] "control"
    ## 
    ## Call:
    ## lm(formula = datos[, i] ~ bx)
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -11.86  -3.69   0.00   1.92  19.23 
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    4.969      5.996   0.829 0.409787    
    ## bx1           18.144      8.480   2.140 0.035542 *  
    ## bx2         -119.090     98.987  -1.203 0.232628    
    ## bx3           42.033     14.664   2.867 0.005350 ** 
    ## bx4           30.415      7.760   3.919 0.000191 ***
    ## bx5           34.432      7.036   4.894 5.31e-06 ***
    ## bx6           20.007      7.060   2.834 0.005871 ** 
    ## bx7           13.692      7.213   1.898 0.061411 .  
    ## bx8            9.774      7.783   1.256 0.213000    
    ## bx9            5.019      9.752   0.515 0.608279    
    ## bx10          21.827     27.477   0.794 0.429426    
    ## bx11        -123.756    151.451  -0.817 0.416373    
    ## bx12              NA         NA      NA       NA    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 5.996 on 77 degrees of freedom
    ## Multiple R-squared:  0.7529, Adjusted R-squared:  0.7176 
    ## F-statistic: 21.33 on 11 and 77 DF,  p-value: < 2.2e-16

    explo.plot(myGCbias, samples = NULL, toplot = "global")

![](RNA_Seq_Data_Analysis_files/figure-markdown_strict/unnamed-chunk-4-4.png)

    # Count depth test
    mycd = dat(data_NOISEQ, type = "cd", norm = FALSE, refColumn = 1)

    ## [1] "Warning: 7876 features with 0 counts in all samples are to be removed for this analysis."
    ## [1] "Reference sample is: GSM7884936"
    ## [1] "Confidence intervals for median of M:"
    ##            0.25%                99.75%               Diagnostic Test
    ## GSM7884937 "-0.155624998483786" "-0.155624998483786" "FAILED"       
    ## GSM7884938 "0.0224454344595077" "0.0856428215624738" "FAILED"       
    ## GSM7884939 "-0.233861498673988" "-0.233861498673988" "FAILED"       
    ## GSM7884940 "0.0381012019025336" "0.131210606294015"  "FAILED"       
    ## GSM7884942 "0.171301784377284"  "0.268216814312323"  "FAILED"       
    ## GSM7884943 "0.301396847025916"  "0.409440771210217"  "FAILED"       
    ## GSM7884944 "-0.274179880493065" "-0.274179880493065" "FAILED"       
    ## GSM7884945 "-0.130374295832555" "-0.130374295832555" "FAILED"       
    ## GSM7884946 "-0.340270518621757" "-0.340270518621757" "FAILED"       
    ## GSM7884947 "-0.208564640914304" "-0.208564640914304" "FAILED"       
    ## [1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

    explo.plot(mycd,samples = 1:6)

![](RNA_Seq_Data_Analysis_files/figure-markdown_strict/unnamed-chunk-4-5.png)

In summary, the analysis suggests the following:

1.  Sequencing depth variation: sequencing depth varies across samples
    (60-69%), indicating the importance of normalization
2.  Length bias in gene expression: longer or shorter genes are more
    likely to be expressed depending on the condition (control
    vs. COVID-19)
3.  GC bias in gene expression: statistically significant, indicating
    there may be GC content bias within the groups
4.  Count depth: diagnostic test failed, recommending normalization

## Differential gene expression analysis

    # Clean up metadata
    pDataUSE <- pData(data_NOISEQ)
    levels(pDataUSE$Group) <- c("Control", "COVID19")
    pDataUSE[pDataUSE=="COVID-19"] <- "COVID19"
    pDataUSE[pDataUSE=="control"] <- "Control"
    pDataUSE[,1] <- as.factor(pDataUSE[,1])

    # Create DESeq2 object
    GSE247186_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE247186_count_filt,
                                               colData = pDataUSE,
                                               design = ~ Group)

    ##   it appears that the last variable in the design formula, 'Group',
    ##   has a factor level, 'Control', which is not the reference level. we recommend
    ##   to use factor(...,levels=...) or relevel() to set this as the reference level
    ##   before proceeding. for more information, please see the 'Note on factor levels'
    ##   in vignette('DESeq2').

    keep <- rowSums(counts(GSE247186_DESeq2) >= 10) >= 6 # Filter out genes that do not have at least 10 counts in at least 6 samples
    GSE247186_DESeq2_F <- GSE247186_DESeq2[keep,]

    # Run DESeq2
    GSE247186_DESeq2_F<- DESeq(GSE247186_DESeq2_F)

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    GSE247186_res <- results(GSE247186_DESeq2_F)
    head(GSE247186_res)

    ## log2 fold change (MLE): Group Control vs COVID19 
    ## Wald test p-value: Group Control vs COVID19 
    ## DataFrame with 6 rows and 6 columns
    ##            baseMean log2FoldChange     lfcSE      stat      pvalue        padj
    ##           <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
    ## 643837      35.4033      0.0995002  0.341033  0.291762 7.70469e-01 8.97696e-01
    ## 148398      10.2287     -0.4182504  0.476755 -0.877286 3.80331e-01 6.35799e-01
    ## 26155       10.9216     -0.3785550  0.477084 -0.793476 4.27500e-01 6.76052e-01
    ## 9636       761.3109      4.2712748  0.825268  5.175625 2.27150e-07 2.82054e-05
    ## 100288175   67.0114     -0.4019391  0.490876 -0.818821 4.12889e-01 6.64947e-01
    ## 8784        19.4787     -0.2855731  0.665114 -0.429360 6.67661e-01 8.42082e-01

    # Filter DEGs based on criteria
    significant_DEGs <- GSE247186_res$padj < 0.05 & abs(GSE247186_res$log2FoldChange) > 1
    num_significant_DEGs <- sum(significant_DEGs, na.rm = TRUE)
    upregulated_genes <- GSE247186_res$log2FoldChange > 1 & GSE247186_res$padj < 0.05
    downregulated_genes <- GSE247186_res$log2FoldChange < -1 & GSE247186_res$padj < 0.05
    num_upregulated_genes <- sum(upregulated_genes, na.rm = TRUE)
    num_downregulated_genes <- sum(downregulated_genes, na.rm = TRUE)

DEGs:

1.  Total significant DEGs: 878
2.  Upregulated genes: 492
3.  Downregulated genes: 386

## Heatmap

Expression of the top 20 DEGs in the dataset

    vsd <- vst(GSE247186_DESeq2_F, blind=FALSE)
    select <- order(rowMeans(counts(GSE247186_DESeq2_F,normalized=TRUE)),
                    decreasing=TRUE)[1:20]
    df <- as.data.frame(colData(GSE247186_DESeq2_F)[,c("Group")])
    colnames(df) <- "Group"
    pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
             cluster_cols=FALSE)

![](RNA_Seq_Data_Analysis_files/figure-markdown_strict/unnamed-chunk-6-1.png)
