## Find differential loops specific to 
## rad21 degron cells

## Load required libraries
library(mariner)
library(SummarizedExperiment)
library(DESeq2)

## Load loop counts
loopCounts <- readRDS("data/rad21Nha9MergedLoopCounts.rds")

## Update path to count data
# path(loopCounts) <- "data/rad21Nha9MergedLoopCounts.h5"

## Build colData from filenames
coldata <- colData(loopCounts)$fileNames |>
    strsplit(split = "_") |>
    do.call(rbind, args = _) |>
    {
        \(x) x[, c(3,4,6)]
    }() |> # subset
    `colnames<-`(c("genotype", "condition", "biorep")) |>
    as.data.frame()
colData(loopCounts)$replicate <- factor(coldata$biorep)
colData(loopCounts)$condition <- relevel(factor(coldata$condition), "Control")

## Build DESeq dataset & run DESeq analysis
counts(loopCounts) <- as.matrix(counts(loopCounts))
dds <- DESeqDataSet(se = loopCounts, design = ~ replicate + condition)
# sizeFactors(dds) <- rep(1, ncol(dds))
dds <- DESeq(dds)


## Get shrunken results
res <- lfcShrink(dds, coef = "condition_5PhIAA_vs_Control", type = "apeglm")

## Add to rowData of InteractionMatrix
rowData(loopCounts) <- cbind(rowData(loopCounts), res)

which(res$padj <= 0.1 & res$log2FoldChange > 0) |> length()
which(res$padj <= 0.1 & res$log2FoldChange < 0) |> length()

## Save object with differential results
saveRDS(loopCounts, file = "data/rad21Nha9DiffLoopCounts.rds")