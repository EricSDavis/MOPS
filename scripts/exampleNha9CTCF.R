## Example plots for the differential
## loops after degradation of CTCF
## with NHA9 present.
## Specifically looking for cases
## where loss of CTCF insulation
## leads to the merging of TADs.
## The method to find these candidates
## has been commented out to avoid
## changes in code affecting the
## example. But uncomment to see
## other potential candidates.

## Load required packages
library(InteractionSet)
library(plotgardener)
library(mariner)
library(grid)

## Load differential looping data
# loopCounts <- readRDS("data/Nha9CtcfDiffLoopCounts.rds")

## Isolate strongly gained loops
# loops <- interactions(loopCounts)
# lostLoops <-
#     loops[which(loops$padj <= 0.1 &
#         loops$log2FoldChange <= 0)]
# gainedLoops <-
#     loops[which(loops$padj <= 0.1 &
#         loops$log2FoldChange >= 2)]

## Define paths to .hic files
nha9File <- "data/raw/hic/condition/MOPS_HCT_CTCFNHA9_Control_0h_inter_30.hic"
degrFile <- "data/raw/hic/condition/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30.hic"

## Identify putative TAD boundaries
## by looking around lost loop ends
# boundaries <- c(
#     GRanges(
#         seqnames = seqnames1(lostLoops),
#         ranges = IRanges(start = start1(lostLoops))
#     ),
#     GRanges(
#         seqnames = seqnames1(lostLoops),
#         ranges = IRanges(start = end2(lostLoops))
#     )
# )
# boundaries <- unique(boundaries)

# gainedRanges <- GRanges(
#     seqnames = seqnames1(gainedLoops),
#     ranges = IRanges(
#         start = start1(gainedLoops),
#         end = end2(gainedLoops)
#     )
# )

## Define regions of interest as 
## gained regions that contain 
## lost loop boundaries
# roi <- subsetByOverlaps(gainedRanges, boundaries)

## Define parameters (maybe 4, 7)
i = 7
p <- pgParams(
    assembly = "hg38",

    ## Region (originally the 7th roi loop)
    chrom = "chrX", #as.character(seqnames(roi)[i]),
    chromstart = 23780000 - 50e3, #start(roi)[i] - 50e3,
    chromend = 24060000 + 50e3, #end(roi)[i] + 50e3,

    ## Hi-C parameters
    resolution = 10e3,
    zrange = c(0, 50),
    norm = "KR",

    ## Default position
    x = 0.5,
    y = 0.5,
    width = 11,
    height = 5,
    space = 0.1
)

png(
    "plots/exampleNha9CtcfLostBoundary.png",
    width = 14, height = 14.65, unit = "in", res = 300
)
pageCreate(12.5, 12, showGuides = TRUE)
plotHicRectangle(
    params = p,
    data = nha9File
)
plotHicRectangle(
    params = p,
    data = degrFile,
    y = '0.1b'
)
dev.off()

## Lost boundary (+CTCF)
png(
    "plots/exampleNha9CtcfLostBoundary.png",
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


## Gained (NHA9) loop (-CTCF)
png(
    "plots/exampleNha9CtcfGainedLoop.png",
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
