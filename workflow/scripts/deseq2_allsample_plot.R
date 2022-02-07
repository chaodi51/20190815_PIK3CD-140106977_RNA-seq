# Overall view of all the samples by PCA, Clustering for expression values, clustering by sample distance
library(DESeq2)
library(pheatmap)
library(genefilter)
library(ggplot2)
library(RColorBrewer)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["count_table"]], header=TRUE, row.names="gene", check.names=FALSE)
coldata <- read.table(snakemake@params[["sample_table"]], header=TRUE, row.names="sample", check.names=FALSE)

dds <- DESeqDataSetFromMatrix(countData=cts, colData=coldata, design = ~ cell_type + condition)
dds$condition <- relevel(dds$condition, "Pre") # use "Pre" as the reference
# Using a grouping variable as contrast 
dds$group <- factor(paste0(dds$cell_type, dds$condition))
design(dds) <- ~ group

# remove uninformative columns
dds <- dds[rowSums(counts(dds)) > 1, ]
# normalization and pre-processing
dds <- DESeq(dds, parallel=parallel)

# raw count normalization
norm_counts <- counts(dds, normalized=TRUE) 
# count transformation, log2 scale, either rlog or vst
vsd <- vst(dds, blind=FALSE)

## visualizations ##
# pca plot
pdf(snakemake@output[["pca_plot"]])
pcaData <- plotPCA(vsd, intgroup=c("condition", "cell_type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=cell_type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

# heatmap of the count matrix
#select <- order(rowMeans(norm_counts), decreasing=TRUE)[1:30] # most highly expressed
select <- head(order(-rowVars(assay(vsd))),35) # most variable genes
pdata <- assay(vsd)[select,]
df <- as.data.frame(colData(dds)[,c("condition","cell_type")])
rownames(df) <- rownames(colData(dds))
colnames(df) <- c("condition","cell_type")
pdf(snakemake@output[["heatmap_exp"]])
pheatmap(pdata, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

# heatmap of sample-sample distances
pdf(snakemake@output[["heatmap_dist"]])
sampleDists = dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$cell_type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

