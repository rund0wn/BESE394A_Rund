# Libraries
library(R.utils)
library(DESeq2)
library(GEOquery)
library(data.table)
library(NOISeq)
library(biomaRt)
library(ggplot2)

# Load data
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
GSE198256_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# Get metadata
gds <- getGEO("GSE198256")
Meta_GSE198256 <- pData(gds$GSE198256_series_matrix.txt.gz@phenoData)
Meta_GSE198256 <- Meta_GSE198256[,c("title","source_name_ch1","characteristics_ch1","characteristics_ch1.1","description","cell type:ch1","disease state:ch1")]
Factors_GSE198256 <- Meta_GSE198256[,c("disease state:ch1")]
GSE198256_count

# Map gene names using biomaRt
ensemblIDs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- rownames(GSE198256_count)
genes <- as.character(genes)
attributes <- c("percentage_gene_gc_content", "start_position", "end_position", "chromosome_name", "gene_biotype", "entrezgene_id")
annotgene <- getBM(attributes = attributes, 
                         filters = "entrezgene_id", 
                         values = genes, 
                         mart = ensemblIDs)
sum(rownames(GSE198256_count) %in% annotgene$entrezgene_id) # Total number of annotated genes
annotgene <- annotgene[annotgene$chromosome_name %in% c(as.character(1:22) ,"X","Y"),] # Only keeping genes in these chromosomes
annotgene_filt <- annotgene[!duplicated(annotgene$entrezgene_id),] # Removing duplicates
rownames(annotgene_filt) <- as.character(annotgene_filt$entrezgene_id) # Checking our genes + annotation overlap
sum(as.character(rownames(annotgene_filt)) %in% rownames(GSE198256_count))

####### QC #######
# Only keeping genes that overlap in both
GSE198256_count_filt <- GSE198256_count[rownames(GSE198256_count) %in% rownames(annotgene_filt),]
GSE198256_count_exc <- GSE198256_count[!(rownames(GSE198256_count) %in% rownames(annotgene_filt)),]
annotgene_ord <- annotgene_filt[rownames(GSE198256_count_filt),]

# NOISeq
Factors_GSE198256 <- data.frame(Meta_GSE198256 [ colnames(GSE198256_count_filt),c("disease state:ch1")])
colnames(Factors_GSE198256)[1]<- "Group"
data_NOISEQ <- readData(data = GSE198256_count_filt,
                        length=abs(annotgene_ord$end-annotgene_ord$start),
                        gc=annotgene_ord$GC,
                        biotype= annotgene_ord$type ,
                        chromosome = annotgene_ord[,c("chromosome_name","start_position","end_position")],
                        factors = Factors_GSE198256)
myexplodata <- dat(data_NOISEQ, type = "countsbio")
explo.plot(myexplodata, plottype = "boxplot")
mynicedata <- dat2save(myexplodata)
mybiodetection <- dat(data_NOISEQ, k = 0, type = "countsbio", factor = NULL)
lengthuse <- abs(annotgene_ord$end-annotgene_ord$start)
names(lengthuse) <- rownames(annotgene_ord)
gc <- annotgene_ord$percentage_gene_gc_content # GC%
names(gc) <- rownames(annotgene_ord)
biotype <-annotgene_ord$gene_biotype
names(biotype) <- rownames(annotgene_ord)
chromosome <- annotgene_ord[,c("chromosome_name","start_position","end_position")]
data_NOISEQ <- readData(data = GSE198256_count_filt,
                        length=lengthuse,
                        gc=gc,
                        biotype= biotype ,
                        chromosome = annotgene_ord[,c("chromosome_name","start_position","end_position")],
                        factors = Factors_GSE198256) # After fixing
myexplodata <- dat(data_NOISEQ, type = "biodetection")
explo.plot(myexplodata, plottype = "persample")
par(mfrow = c(1, 2))
explo.plot(myexplodata, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")
mycountsbio = dat(data_NOISEQ, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")
mysaturation = dat(data_NOISEQ, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:2, yleftlim = NULL, yrightlim = NULL)
explo.plot(mysaturation, toplot = "protein_coding", samples = 1:4)
explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")
mylengthbias = dat(data_NOISEQ, factor = "Group", type = "lengthbias")
explo.plot(mylengthbias, samples = NULL, toplot = "global")

# Check GC bias
myGCbias = dat(data_NOISEQ, factor = "Group", type = "GCbias")
explo.plot(myGCbias, samples = NULL, toplot = "global")
mycd = dat(data_NOISEQ, type = "cd", norm = FALSE, refColumn = 1)
explo.plot(mycd,samples = 1:12)
myPCA = dat(data_NOISEQ, type = "PCA")
explo.plot(myPCA, factor = "Group")
QCreport(data_NOISEQ, samples = NULL, factor = "Group", norm = FALSE)
save(data_NOISEQ,GSE198256_count_filt,annotgene_ord,file="GSE198256_step1.Rda")

####### Normalization and Differential Expression #######
myRPKM = rpkm(assayData(data_NOISEQ)$exprs, long = lengthuse, k = 0, lc = 1)
myUQUA = uqua(assayData(data_NOISEQ)$exprs, long = lengthuse, lc = 0.5, k = 0)
myTMM = tmm(assayData(data_NOISEQ)$exprs, long = 1000, lc = 0)

# DESeq2
pDataUSE <- pData(data_NOISEQ)
pDataUSE[pDataUSE=="Covid19: Acute infection"] <- "Acute"
pDataUSE[pDataUSE=="Covid19: Recovery 3Mo"] <- "EarlyRecovery"
pDataUSE[pDataUSE=="Covid19: Recovery 6Mo"] <- "LateRecovery"
pDataUSE[,1] <- as.factor(pDataUSE[,1])
resultsNames(GSE198256_DESeq2)
GSE198256_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt,
                                           colData = pDataUSE,
                                           design = ~ Group)

keep <- rowSums(counts(GSE198256_DESeq2) >= 20) >= 1 # Filter out genes that do not have at least 10 counts in at least 6 samples
GSE198256_DESeq2_F <- GSE198256_DESeq2[keep,] # DE analysis
GSE198256_DESeq2_F<- DESeq(GSE198256_DESeq2_F)
GSE198256_res <- results(GSE198256_DESeq2_F)
GSE198256_res
resultsNames(GSE198256_DESeq2_F)


####### Compare with Publication #######
res_acute_vs_control <- results(GSE198256_DESeq2_F, contrast = c("Group", "Acute", "Healthy"))
res_early_recovery_vs_control <- results(GSE198256_DESeq2_F, contrast = c("Group", "EarlyRecovery", "Healthy"))
res_late_recovery_vs_control <- results(GSE198256_DESeq2_F, contrast = c("Group", "LateRecovery", "Healthy"))

# Function to count significant DEGs
count_sig_DEGs <- function(res) {
  sig_DEGs <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > log2(2)), ]
  upregulated <- sum(sig_DEGs$log2FoldChange > 0)
  downregulated <- sum(sig_DEGs$log2FoldChange < 0)
  return(list("upregulated" = upregulated, "downregulated" = downregulated))
}

# Counting DEGs for each comparison
DEGs_acute_vs_control <- count_sig_DEGs(res_acute_vs_control)
DEGs_early_recovery_vs_control <- count_sig_DEGs(res_early_recovery_vs_control)
DEGs_late_recovery_vs_control <- count_sig_DEGs(res_late_recovery_vs_control)

# Print the results
print(DEGs_acute_vs_control)
print(DEGs_early_recovery_vs_control)
print(DEGs_late_recovery_vs_control)

# Compare specific genes
gene_symbols <- c("JUNB", "ATF3", "NFkB2")
genes_mapped <- getBM(attributes = c('hgnc_symbol', 'entrezgene_id'),
                      filters = 'hgnc_symbol',
                      values = gene_symbols,
                      mart = ensemblIDs)
specific_genes_data <- res_acute_vs_control[rownames(res_acute_vs_control) %in% genes_mapped$entrezgene_id, ]
print(specific_genes_data)
print(genes_mapped)

####### ORA with ClusterProfiler #######
GSE198256_res <- na.omit(GSE198256_res)
sig_genes <- rownames(GSE198256_res[which(GSE198256_res$padj < 0.05), ])
sig_genes <- unique(sig_genes)
nonsig_genes <- unique(rownames(GSE198256_res))

ego <- enrichGO(gene         = sig_genes,
                universe = nonsig_genes,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable     = TRUE)
head(ego)[, 1:7]
dotplot(ego)
barplot(ego, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)

ekg <- enrichKEGG(gene         = sig_genes,
                  universe = nonsig_genes,
                  organism     = 'hsa', # For human
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

head(ekg)
barplot(ekg, showCategory=20)

####### GSEA #######
GSE198256_res_sig <- GSE198256_res[which(GSE198256_res$padj < 0.05), ]
geneList <- GSE198256_res_sig$log2FoldChange
names(geneList) <- rownames(GSE198256_res_sig)
geneList = sort(geneList, decreasing = TRUE)

gse <- gseGO(geneList = geneList,
             ont = "BP",
             OrgDb = org.Hs.eg.db,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH")

dotplot(gse, showCategory = 10)
