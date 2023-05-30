## Find NHA9 loops

## Load required libraries
library(mariner)
library(SummarizedExperiment)
library(DESeq2)

## Load loop counts
loopCounts <- readRDS("data/ControlMergedLoopCounts.rds")

## Update path to count data
# path(loopCounts) <- "data/ControlMergedLoopCounts.h5"

## Build colData from filenames
coldata <- 
    colData(loopCounts)$fileNames |>
    strsplit(split = "_") |>
    do.call(rbind, args=_) |>
    { \(x) data.frame(condition = x[, 3])} ()

## Set condition & replicates
coldata$condition <- gsub(
    pattern = ".*(NHA9|parental)",
    replacement = "\\1",
    x = coldata$condition
) |>
factor() |>
relevel("parental")
coldata$replicate <- factor(c(1,2,1,2,3,4,3,4))
colData(loopCounts) <- cbind(colData(loopCounts), coldata)

## Build DESeq dataset & run DESeq analysis
counts(loopCounts) <- as.matrix(counts(loopCounts))
dds <- DESeqDataSet(
    se = loopCounts,
    design = ~ replicate + condition
)
# sizeFactors(dds) <- rep(1, ncol(dds))
dds <- DESeq(dds)

## Get shrunken results
res <- lfcShrink(
    dds,
    coef = "condition_NHA9_vs_parental",
    type = "apeglm"
)

## Visualize results
pdf("plots/ParentalVsNha9_MAplot.pdf")
plotMA(res)
dev.off()

## Add to rowData of InteractionMatrix
rowData(loopCounts) <- cbind(rowData(loopCounts), res)

which(res$padj <= 0.1 & res$log2FoldChange > 0) |> length()
which(res$padj <= 0.1 & res$log2FoldChange < 0) |> length()

## Save object with differential results
saveRDS(loopCounts, file = "data/ControlDiffLoopCounts.rds")
