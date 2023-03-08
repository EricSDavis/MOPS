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

## Define paths to .hic files
ctrlFile <- "data/raw/hic/condition/MOPS_HCT_CTCFparental_Control_0h_inter_30.hic"
nha9File <- "data/raw/hic/condition/MOPS_HCT_CTCFNHA9_Control_0h_inter_30.hic"
degrFile <- "data/raw/hic/condition/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30.hic"

## Define parameters
p <- pgParams(
    assembly = "hg38",

    ## Region (originally the second differential loop)
    chrom = "chr9", #seqnames1(diffLoops)[2],
    chromstart = 125650000, #start1(diffLoops)[2] - 100e3,
    chromend = 126010000, #end2(diffLoops)[2] + 100e3,

    ## Hi-C parameters
    resolution = 5e3,
    zrange = c(0, 50),
    norm = "KR",

    ## Default position
    x = 0.5,
    y = 0.5,
    width = 11,
    height = 5,
    space = 0.1
)

## HCT parental (-NHA9)
png(
    "plots/exampleParentalVsNha9_CTCF_-NHA9.png",
    width=12, height=5.65, unit='in', res=300
)
pageCreate(12.5, 6, showGuides = FALSE)
ctrlHicPlot <- plotHicRectangle(
    params = p,
    data = ctrlFile
)
plotGenomeLabel(
    params = p,
    length = p$width,
    y = "0b",
    fontsize = 20
)
annoHeatmapLegend(
    plot = ctrlHicPlot,
    x = p$x + p$width + p$space*2,
    y = p$y,
    width = p$space * 2,
    height = p$height * 0.75,
    fontcolor = "black",
    fontsize = 18
)
dev.off()

## HCT +NHA9
png(
    "plots/exampleParentalVsNha9_CTCF_+NHA9.png",
    width = 12, height = 5.65, unit = "in", res = 300
)
pageCreate(12.5, 6, showGuides = FALSE)
nha9HicPlot <- plotHicRectangle(
    params = p,
    data = nha9File
)
plotGenomeLabel(
    params = p,
    length = p$width,
    y = "0b",
    fontsize = 20
)
annoHeatmapLegend(
    plot = nha9HicPlot,
    x = p$x + p$width + p$space * 2,
    y = p$y,
    width = p$space * 2,
    height = p$height * 0.75,
    fontcolor = "black",
    fontsize = 18
)
dev.off()

## HCT +NHA9 -CTCF
png(
    "plots/exampleParentalVsNha9_CTCF_+NHA9-CTCF.png",
    width = 12, height = 5.65, unit = "in", res = 300
)
pageCreate(12.5, 6, showGuides = FALSE)
degrHicPlot <- plotHicRectangle(
    params = p,
    data = degrFile
)
plotGenomeLabel(
    params = p,
    length = p$width,
    y = "0b",
    fontsize = 20
)
annoHeatmapLegend(
    plot = degrHicPlot,
    x = p$x + p$width + p$space * 2,
    y = p$y,
    width = p$space * 2,
    height = p$height * 0.75,
    fontcolor = "black",
    fontsize = 18
)
dev.off()