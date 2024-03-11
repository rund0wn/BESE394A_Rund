
################################################################################
#                            Libraries and Data                                #
################################################################################
options(width = 210)
library(Seurat)
library(TFBSTools)
library(motifmatchr)
library(Signac)
library(ggplot2)
library(Azimuth)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)

# Plotting and saving function
plot_QC_features <- function(seurat_object, pdf_path){
  pdf(pdf_path, width =10) 
  print(VlnPlot(
    object = seurat_object,
    features = c("nCount_RNA", "nCount_ATAC", 'nFeature_RNA',  'mitoPct', 
                 'FRiP', 'log10GenesPerUMI', "TSS.enrichment", 
                 "nucleosome_signal"),
    pt.size = 0, ncol =4))
  DefaultAssay(seurat_object) <- "ATAC"
  print(FragmentHistogram(object = seurat_object, region = 'chr1-1-10000000', group.by = 'nucleosome_group'))
  print(TSSPlot(seurat_object, group.by = 'high.tss') + NoLegend())
  print(DensityScatter(seurat_object, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE))
  print(ggplot(seurat_object@meta.data, aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          scale_x_log10() + 
          geom_vline(xintercept = 300) + ggtitle('GENESpercell'))
  print(ggplot(seurat_object@meta.data, aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
          geom_boxplot() + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          theme(plot.title = element_text(hjust=0.5, face="bold")) +
          ggtitle("NCells vs NGenes"))
  print(ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, colour=mitoPct, group=orig.ident)) + 
          geom_point() + 
          scale_colour_gradient(low = "gray90", high = "black", limits=c(0,100)) +
          stat_smooth(method=lm) +
          scale_x_log10() + 
          scale_y_log10() + 
          theme_classic() +
          geom_vline(xintercept = 500) +
          geom_hline(yintercept = 6000) +
          geom_hline(yintercept = 250) + ggtitle(paste0('UMIvsGENESpercell  Ncell:: ', ncol(seurat_object))))
  dev.off()
  
}

# Load data tables and 10X Genomics data
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.tar.gz
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz.tbi

counts <- Read10X("/home/tawfiqre/Downloads/Class_Homework11-3/filtered_feature_bc_matrix")
fragpath <- '/home/tawfiqre/Downloads/Class_Homework11-3/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz'

# Download gene annotations
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# Create Seurat object for expression data
mouse_brain <- CreateSeuratObject(
  counts = counts[['Gene Expression']],
  assay = "RNA"
)

# Prepare ATAC-seq data for analysis
atac_counts <- counts$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Merge ATAC-seq data and expression data
mouse_brain[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

################################################################################
#                         Calculate metrics for QC                             #
################################################################################

# Calculate QC features, such as counts and signals
DefaultAssay(mouse_brain) <- "ATAC"

mouse_brain <- NucleosomeSignal(mouse_brain)
mouse_brain <- TSSEnrichment(mouse_brain, fast=FALSE)
total_fragments <- CountFragments(fragments = fragpath)
rownames(total_fragments) <- total_fragments$CB
mouse_brain$fragments <- total_fragments[colnames(mouse_brain), "frequency_count"]
mouse_brain <- FRiP(object = mouse_brain, assay = 'ATAC', total.fragments = 'fragments')
mouse_brain$nucleosome_group <- ifelse(mouse_brain$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
mouse_brain$high.tss <- ifelse(mouse_brain$TSS.enrichment > 3, 'High', 'Low')

DefaultAssay(mouse_brain) <- "RNA"
mouse_brain$mitoPct <- PercentageFeatureSet(mouse_brain, pattern = "^mt-")
mouse_brain$RPSPct  <- PercentageFeatureSet(object = mouse_brain, pattern = "^Rp[sl]")
mouse_brain$log10GenesPerUMI <- log10(mouse_brain$nFeature_RNA) / log10(mouse_brain$nCount_RNA)

# Plot calculated features and save to PDF, save data as .rds
plot_QC_features(mouse_brain, '/home/kurowsaa/5k_mouse_brain_GEX_QC_Pre.pdf')
saveRDS(mouse_brain, './mouse_brain_multiome_Pre_QC.rds')

################################################################################
#                   Filtering based on calculated QC metrics                   #
################################################################################

# Criteria: <100k read counts, <25k RNA counts, <1k ATAC/RNA counts, high nucleosome, <1 TSS enrichment, FRiP <0.2

mouse_brain <- subset(
  x = mouse_brain,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 & 
    FRiP > 0.2
)

# Save filtered object and QC plots after filtering
saveRDS(mouse_brain, '/home/tawfiqre/Downloads/Class_Homework11-3/Data/mouse_brain_multiome_Post_QC.rds')
plot_QC_features(mouse_brain, '/home/tawfiqre/Downloads/Class_Homework11-3/Plots/5k_mouse_brain_GEX_QC_Post.pdf')

################################################################################
#    Process the data (normalization, scaling, PCA, clustering, UMAP, ....)    #
################################################################################

## Transcriptomic analysis
DefaultAssay(mouse_brain) <- "RNA" # Set assay to RNA
mouse_brain <- NormalizeData(mouse_brain) # Normalize read counts to make them comparable
mouse_brain <- FindVariableFeatures(mouse_brain, 
            selection.method = "vst", 
            nfeatures = 2000) # Identify top 2000 variable features across cells

# For PCA: scale, find neighbors, cluster cells and perform PCA
mouse_brain <- ScaleData(mouse_brain)
mouse_brain <- FindNeighbors(mouse_brain, dims = 1:30)
mouse_brain <- FindClusters(mouse_brain, 
            resolution = 0.4, 
            algorithm = 3, 
            cluster.name="RNA_clusters_03")
mouse_brain <- RunPCA(mouse_brain)

# Run UMAP based on PCA results
mouse_brain <- RunUMAP(mouse_brain, dims=1:30, 
            reduction = "pca", 
            reduction.name = "rna_umap")

# Map cells to reference data
mouse_brain <- Azimuth::RunAzimuth(mouse_brain, 
            reference = "mousecortexref")


pdf('/home/tawfiqre/Downloads/Class_Homework11-3/Plots/UMAP_RNA.pdf') # Save results to PDF
DimPlot(mouse_brain, reduction = "rna_umap", group.by='RNA_clusters_03', 
        label = TRUE, repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "rna_umap", group.by = "predicted.subclass", 
        label = TRUE, label.size = 3) + NoLegend()
dev.off()

## Epigenomic analysis

DefaultAssay(mouse_brain) <- "ATAC" # Set assay to ATAC
mouse_brain <- FindTopFeatures(mouse_brain, min.cutoff = 5) # Identify top features
mouse_brain <- RunTFIDF(mouse_brain) # Calculate scores for peak regions
mouse_brain <- RunSVD(mouse_brain) # SVD for dimensionality reduction

pdf('/home/tawfiqre/Downloads/Class_Homework11-3/Plots/Correlation.pdf') # Save results to pdf
DepthCor(mouse_brain)
dev.off()

mouse_brain <- RunUMAP(mouse_brain, dims=1:30, 
                reduction = "lsi", 
                reduction.name = "atac_umap") # UMAP for ATAC based on LSI results

pdf('/home/tawfiqre/Downloads/Class_Homework11-3/Plots/UMAP_ATAC.pdf')
DimPlot(mouse_brain, reduction = "atac_umap", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "atac_umap", group.by = "predicted.subclass",
label = TRUE, label.size = 3) + NoLegend()
dev.off()

saveRDS(mouse_brain, '/home/tawfiqre/Downloads/Class_Homework11-3/Data/mouse_brain_multiome.rds')

################################################################################
#                        Integrate the RNA and ATAC data                       #
################################################################################

# Use RNA and ATAC data to find multimodal neighbors
mouse_brain <- FindMultiModalNeighbors(
  object = mouse_brain,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 1:40),
  verbose = TRUE
)

# Build a joint UMAP visualization
mouse_brain <- RunUMAP(
  object = mouse_brain, 
  reduction.name = "wnn.umap", dims = 1:30,
  assay = "RNA",
  verbose = TRUE
)

saveRDS(mouse_brain, '/home/tawfiqre/Downloads/Class_Homework11-3/Data/mouse_brain_multiome_Integrated.rds')

pdf('/home/tawfiqre/Downloads/Class_Homework11-3/Plots/Umap_integrated.pdf')
DimPlot(mouse_brain, reduction = "rna_umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "atac_umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "wnn.umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
dev.off()

pdf('/home/tawfiqre/Downloads/Class_Homework11-3/Plots/Umap_integrated_weights.pdf')
Idents(mouse_brain) <- "predicted.subclass"
VlnPlot(mouse_brain, features = c("RNA.weight", "ATAC.weight"),  pt.size = 0, ncol = 1)
dev.off()

################################################################################
#                  Perform a Differential accessibility analysis               #
################################################################################

# Perform differential accessibility (DA) analysis between Oligodendrocytes and Astrocytes using a logistic regression (LR) model, adjusting for total ATAC read count
DefaultAssay(mouse_brain) <- 'ATAC'
Idents(mouse_brain) <- "predicted.subclass" # Set the cell identity classes to 'predicted.subclass' for differential analysis
da_peaks <- FindMarkers(
  object = mouse_brain,
  ident.1 = c("Oligo"), 
  ident.2 = c("Astro"),
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)
saveRDS(da_peaks, '/home/tawfiqre/Downloads/Class_Homework11-3/Data/DA_peaks.rds')

da_peaks <- da_peaks[order(da_peaks$avg_log2FC, decreasing = TRUE), ] # Sort DA peaks by average log2 fold change in descending order for visualization priority

pdf('/home/tawfiqre/Downloads/Class_Homework11-3/Plots/Peaks_DA.pdf', width=12, height=20)
cowplot::plot_grid(
  VlnPlot(
    object = mouse_brain,
    assay = 'ATAC',
    features = rownames(da_peaks)[1:3],
    pt.size = 0.1,
    group.by='predicted.subclass', 
    ncol=3
  ),
  VlnPlot(   # Violin plots of the expression of genes closest to the top 3 DA peaks in RNA data
    object = mouse_brain,
    assay = 'RNA',
    features = ClosestFeature(mouse_brain, rownames(da_peaks)[1:3])$gene_name,
    pt.size = 0.1,
    group.by='predicted.subclass', 
    ncol=3
  ),
  FeaturePlot(   # Feature plots showing spatial distribution of the top 3 DA peaks in UMAP space
    object = mouse_brain,
    reduction = 'wnn.umap', 
    order=TRUE,
    features = rownames(da_peaks)[1:3],
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol=3
  ) & NoLegend(),
  FeaturePlot(   # Feature plots showing spatial distribution of genes closest to the top 3 DA peaks in UMAP space
    object = mouse_brain,
    reduction = 'wnn.umap', 
    order=TRUE,
    features = ClosestFeature(mouse_brain, rownames(da_peaks)[1:3])$gene_name,
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol=3
  ) & NoLegend(),
nrow=4)
dev.off()


# Annotate peaks with his closest feature
open_Oligo <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_Astro <- rownames(da_peaks[da_peaks$avg_log2FC < 3, ])
closest_Oligo <- ClosestFeature(mouse_brain, open_Oligo)
closest_Astro <- ClosestFeature(mouse_brain, open_Astro)

# https://www.cellsignal.com/pathways/neuronal-and-glial-cell-markers
# Visualize the coverage of the peaks
pdf('/home/tawfiqre/Downloads/Class_Homework11-3/Plots/Coverage_selected.pdf', height=12)
CoveragePlot(
  object = mouse_brain,
  region = c('Olig1', 
             'Gfap'),
  extend.upstream = 1000,
  extend.downstream = 1000,
  ncol = 1
)
dev.off()


################################################################################
#                              Motif analysis                                  #
################################################################################

# Motif analysis with the DA peaks
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Add motif information
mouse_brain <- AddMotifs(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
DefaultAssay(mouse_brain) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 'Mus musculus', all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(mouse_brain), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
mouse_brain <- SetAssayData(mouse_brain, assay = 'ATAC', layer = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
mouse_brain <- RunChromVAR(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

enriched.motifs <- FindMotifs(
  object = mouse_brain,
  features = top.da.peak
)

pdf('/home/tawfiqre/Downloads/Class_Homework11-3/Plots/Motif_enrichment.pdf')
MotifPlot(
  object = mouse_brain,
  motifs = head(rownames(enriched.motifs))
)
dev.off()

################################################################################
#           Generate a RNA activity matrix based on the ATAC-seq data          #
################################################################################

# We can create a proxy of the gene expression from the ATAC-seq data using  GeneActivity function. 
# We can also create a proxy of the gene expression from the ATAC-seq data using the chromVAR package. 
# This package uses the motif accessibility to infer the gene expression. 
# We can use the motif models from JASPAR2020 to perform this analysis.

gene.activities <- GeneActivity(mouse_brain)
mouse_brain[['RNA_ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
mouse_brain <- NormalizeData(
  object = mouse_brain,
  assay = 'RNA_ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(mouse_brain$nCount_RNA)
)

DefaultAssay(mouse_brain) <- 'RNA'
RNA_plot <- FeaturePlot(
    object = mouse_brain,
    order=TRUE,
    features =c("Sema5a","Dennd4a","Nkain1"),
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol = 3
  )& NoLegend()

DefaultAssay(mouse_brain) <- 'RNA_ACTIVITY'
RNA_Activity_plot <- FeaturePlot(
    object = mouse_brain,
    order=TRUE,
    features =c("Sema5a","Dennd4a","Nkain1"),
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol = 3
  )& NoLegend()

pdf('/home/tawfiqre/Downloads/Class_Homework11-3/Plots/RNA_comparison.pdf', height=12)
cowplot::plot_grid(
  RNA_plot,
  RNA_Activity_plot,
nrow=2)
dev.off()


################################################################################
#                             Group 3 Assignment                               #
################################################################################
# Generate motifs for the top 5 differentially accessible peaks per cell type. 
# Use the JASPAR database to obtain the motifs.

mouse_brain <- readRDS('/home/tawfiqre/Downloads/Class_Homework11-3/mouse_brain_multiome_Integrated.rds')
DefaultAssay(mouse_brain) <- 'ATAC'
Idents(mouse_brain) <- "predicted.subclass"

# Cell line analysis
cell_lines <- c('Pvalb', 'Lamp5', 'Endo', 'OPC')
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
) # JASPAR motifs
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 'Mus musculus', all_versions = FALSE))

# Loop through cell lines
for (cell_line in cell_lines) {
  da_peaks <- FindMarkers(object = mouse_brain, ident.1 = cell_line, test.use = 'LR', latent.vars = 'nCount_ATAC')

  da_peaks <- da_peaks[order(da_peaks$avg_log2FC, decreasing = TRUE), ]
  top_da_peaks <- rownames(da_peaks[da_peaks$p_val < 0.005, ])[1:5]
  saveRDS(da_peaks, paste0('/home/tawfiqre/Downloads/Class_Homework11-3/DA_peaks_', cell_line, '.rds'))
  generate_and_save_plots(mouse_brain, da_peaks, cell_line)
  mouse_brain <- AddMotifs(object = mouse_brain, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)
  
  # Update the mouse_brain object with motif data for ATAC assay
  motif.matrix <- CreateMotifMatrix(features = granges(mouse_brain), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
  motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
  mouse_brain <- SetAssayData(mouse_brain, assay = 'ATAC', layer = 'motifs', new.data = motif.object)

  # Run ChromVAR and find enriched motifs for the top differentially accessible peaks
  mouse_brain <- RunChromVAR(object = mouse_brain, genome = BSgenome.Mmusculus.UCSC.mm10)
  enriched.motifs <- FindMotifs(object = mouse_brain, features = top_da_peaks)
  saveRDS(enriched.motifs, paste0('./Motif_enrichment_', cell_line, '.rds'))

  # Plot the enriched motifs
  pdf(paste0('./Motif_enrichment_', cell_line, '.pdf'))
  MotifPlot(object = mouse_brain, motifs = head(rownames(enriched.motifs)))
  dev.off()
}

# This function encapsulates the plotting logic for DA peaks and associated genes
generate_and_save_plots <- function(mouse_brain, da_peaks, cell_line) {
  pdf(paste0('./Peaks_DA_', cell_line, '.pdf'), width=12, height=20)
  cowplot::plot_grid(
    VlnPlot(object = mouse_brain, assay = 'ATAC', features = rownames(da_peaks)[1:3], pt.size = 0.1, group.by = 'predicted.subclass', ncol = 3),
    VlnPlot(object = mouse_brain, assay = 'RNA', features = ClosestFeature(mouse_brain, rownames(da_peaks)[1:3])$gene_name, pt.size = 0.1, group.by = 'predicted.subclass', ncol = 3),
    FeaturePlot(object = mouse_brain, reduction = 'wnn.umap', order = TRUE, features = rownames(da_peaks)[1:3], pt.size = 0.1, max.cutoff = 'q95', ncol = 3) & NoLegend(),
    FeaturePlot(object = mouse_brain, reduction = 'wnn.umap', order = TRUE, features = ClosestFeature(mouse_brain, rownames(da_peaks)[1:3])$gene_name, pt.size = 0.1, max.cutoff = 'q95', ncol = 3) & NoLegend(),
    nrow = 4
  )
  dev.off()
}
