## Survey plots for the differential
## loops called between parental and
## NHA9 cells (no auxin).
## These are the NHA9 loops in HCT
## cells.

## Load required packages
library(InteractionSet)
library(plotgardener)
library(mariner)
library(grid)


## Load differential looping data
loopCounts <- readRDS("data/ControlDiffLoopCounts.rds")

## Isolate diff loop results
loops <- interactions(loopCounts)
diffLoops <- loops[which(loops$padj <= 0.1)]

## Define paths to .hic files
ctrlFile <- "data/raw/hic/condition/MOPS_HCT_CTCFparental_Control_0h_inter_30.hic"
nha9File <- "data/raw/hic/condition/MOPS_HCT_CTCFNHA9_Control_0h_inter_30.hic"
degrFile <- "data/raw/hic/condition/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30.hic"

pdf("plots/surveyParentalVsNha9_CTCF.pdf", 5, 7)
for (i in seq_along(diffLoops)) {
    ## Define parameters
    p <- pgParams(
        assembly = "hg38",

        ## Region
        chrom = seqnames1(diffLoops)[i],
        chromstart = start1(diffLoops)[i] - 100e3,
        chromend = end2(diffLoops)[i] + 100e3,

        ## Hi-C parameters
        resolution = 5e3,
        zrange = c(0, 50),
        norm = "KR",

        ## Default position
        x = 0.5,
        y = 0.5,
        width = 4,
        height = 1.5,
        space = 0.1
    )

    ## Begin visualization
    pageCreate(width = 5, height = 7, showGuides=FALSE)
    ctrlHicPlot <- plotHicRectangle(
        params = p,
        data = ctrlFile
    )
    nha9HicPlot <- plotHicRectangle(
        params = p,
        data = nha9File,
        y = "0.1b"
    )
    degrHicPlot <- plotHicRectangle(
        params = p,
        data = degrFile,
        y = '0.1b'
    )
    plotGenes(
        params = p,
        height = 0.75,
        y = "0.1b"
    )
    plotGenomeLabel(
        params = p,
        length = p$width,
        y = "0b"
    )

    annoHeatmapLegend(
        plot = ctrlHicPlot,
        x = p$x + p$width + p$space,
        y = p$y,
        width = p$space,
        height = p$height * 0.75,
        fontcolor = "black"
    )
    annoHeatmapLegend(
        plot = nha9HicPlot,
        x = p$x + p$width + p$space,
        y = nha9HicPlot$y,
        width = p$space,
        height = p$height * 0.75,
        fontcolor = "black"
    )
    annoHeatmapLegend(
        plot = degrHicPlot,
        x = p$x + p$width + p$space,
        y = degrHicPlot$y,
        width = p$space,
        height = p$height * 0.75,
        fontcolor = "black"
    )

    plotText(
        label = "HCT",
        x = ctrlHicPlot$x + unit(p$space / 2, "inches"),
        y = ctrlHicPlot$y + unit(p$space / 2, "inches"),
        just = c("left", "top")
    )
    plotText(
        label = "HCT +NHA9",
        x = nha9HicPlot$x + unit(p$space / 2, "inches"),
        y = nha9HicPlot$y + unit(p$space / 2, "inches"),
        just = c("left", "top")
    )
    plotText(
        label = "HCT +NHA9 -CTCF",
        x = degrHicPlot$x + unit(p$space / 2, "inches"),
        y = degrHicPlot$y + unit(p$space / 2, "inches"),
        just = c("left", "top")
    )
}
dev.off()