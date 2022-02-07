## Setup contrast as Post vs Pre Leniolisib treatment (PI3Kδ inhibitor) in four cell types: 
## CD4Tn, CD4Tnn, CD8Tn, CD8Tnn, oupout tables with the standard DESeq2 result format including " 
## baseMean log2FoldChange lfcSE stat pvalue padj" plus the “Log fold change shrinked” normalized readcounts 

library(DESeq2)
library(pheatmap)
library(genefilter)
library(ggplot2)
library(RColorBrewer)
library(fdrtool)
library(EnhancedVolcano)

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

# get the current contrast/cell_type from snakemake output, e.g., "CD8Tnn_Post_vs_Pre"
output_file <- snakemake@output[["table"]]
comp = gsub(".diffexp.tsv", "", tail(unlist(strsplit(output_file, "/")),1))
cell = gsub("_Post_vs_Pre", "", comp)
# get Post_vs_Pre by groups(i.e., cell_types)
res <- results(dds, contrast = c("group", paste0(cell,"Post"), paste0(cell,"Pre")), parallel = parallel)

# use fdrtool to correct the overestimated p-value,
# https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html
res <- res[!is.na(res$pvalue),]
res <- res[!is.na(res$padj),]
res <- res[,-which(names(res)=="padj")]
FDR.res <- fdrtool(res$stat, statistic="normal", plot=F)
res[,"padj"]  <- p.adjust(FDR.res$pval, method = "BH")
message(comp, paste0(" : # Up = ", length(res[which(res$padj<=0.1 & res$log2FoldChange>0),]$padj)),
        paste0(" # Down = ", length(res[which(res$padj<=0.1 & res$log2FoldChange<0),]$padj)))

# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast = c("group", paste0(cell,"Post"), paste0(cell,"Pre")), res=res, type="ashr")
# extract the current cell_type samples
df_vsd = as.data.frame(assay(vsd))
df_vsd_cell = df_vsd[,grep(paste0(cell,'$'),colnames(df_vsd))]
# merge with normalized count data and output the table
resdata <- merge(df_vsd_cell, as.data.frame(res), by="row.names",sort=FALSE)
names(resdata)[1] <- "Gene"
#print(head(resdata))
write.table(resdata, file=snakemake@output[["table"]], sep="\t", quote=FALSE, row.names=FALSE)

## basic plots for Data quality assessment
# M-A plot, points are red with padj < 0.1, points fall out of the window are open triangles 
pdf(snakemake@output[["ma_plot"]])
plotMA(res, main=comp, colLine="red")
dev.off()

