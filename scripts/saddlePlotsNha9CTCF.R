## Example plots for the differential
## loops after degradation of CTCF
## with NHA9 present.
## Specifically looking for cases
## where loss of CTCF insulation
## leads to the merging of TADs.

## Load required packages
library(InteractionSet)
library(plotgardener)
library(mariner)
library(grid)

## Load differential looping data
loopCounts <- readRDS("data/Nha9CtcfDiffLoopCounts.rds")

## Isolate strongly gained loops
loops <- interactions(loopCounts)
lostLoops <-
    loops[which(loops$padj <= 0.1 &
        loops$log2FoldChange <= 0)]
gainedLoops <-
    loops[which(loops$padj <= 0.1 &
        loops$log2FoldChange >= 2)]

## Define paths to .hic files
nha9File <- "data/raw/hic/condition/MOPS_HCT_CTCFNHA9_Control_0h_inter_30.hic"
degrFile <- "data/raw/hic/condition/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30.hic"


## Identify putative TAD boundaries
## by looking around lost loop ends
boundaries <- c(
    GRanges(
        seqnames = seqnames1(lostLoops),
        ranges = IRanges(start = start1(lostLoops))
    ),
    GRanges(
        seqnames = seqnames1(lostLoops),
        ranges = IRanges(start = end2(lostLoops))
    )
)
boundaries <- unique(boundaries)

gainedRanges <- GRanges(
    seqnames = seqnames1(gainedLoops),
    ranges = IRanges(
        start = start1(gainedLoops),
        end = end2(gainedLoops)
    )
)

## Define regions of interest as 
## gained regions that contain 
## lost loop boundaries
roi <- subsetByOverlaps(gainedRanges, boundaries)

## Define parameters (maybe 4, 7)
i = 7
p <- pgParams(
    assembly = "hg38",

    ## Region (originally the second differential loop)
    chrom = as.character(seqnames(roi)[i]),
    chromstart = start(roi)[i] - 50e3,
    chromend = end(roi)[i] + 50e3,

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
    "plots/temp.png",
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


## Testing something
## Pull matrices in a square
## around the TAD boundaries
## and compare between Hi-C
## files
gr1 <- GRanges(
    seqnames = seqnames(boundaries),
    ranges = IRanges(
        start = start(boundaries) - 500e3,
        end = start(boundaries) + 500e3
    )
)
library(InteractionSet)
gi1 <- GInteractions(gr1, gr1)

tmp1 <- pullHicMatrices(
    x = gi1,
    files = c(nha9File, degrFile),
    binSize = 10e3,
    blockSize = 6e6
)

mats1 <- aggHicMatrices(tmp1, by = "files")

png(
    "plots/temp.png",
    width=5, height= 5, units='in', res = 300
)
plotMatrix(data = mats[, , 2] / mats[, , 1],
palette = colorRampPalette(c("blue", "white", "red")))
dev.off()

## Control using all loop boundaries
staticLoops <- loops[which(loops$padj > 0.1)]

boundaries2 <- c(
    GRanges(
        seqnames = seqnames1(staticLoops),
        ranges = IRanges(start = start1(staticLoops))
    ),
    GRanges(
        seqnames = seqnames1(staticLoops),
        ranges = IRanges(start = end2(staticLoops))
    )
)
boundaries2 <- unique(boundaries2)

gr2 <- GRanges(
    seqnames = seqnames(boundaries2),
    ranges = IRanges(
        start = start(boundaries2) - 500e3,
        end = start(boundaries2) + 500e3
    )
)
library(InteractionSet)
gi2 <- GInteractions(gr2, gr2)

tmp2 <- pullHicMatrices(
    x = gi2,
    files = c(nha9File, degrFile),
    binSize = 10e3,
    blockSize = 6e6
)

mats2 <- aggHicMatrices(tmp2, by = "files")

png(
    "plots/temp.png",
    width = 5, height = 5, units = "in", res = 300
)
plotMatrix(
    data = mats2[, , 2] / mats2[, , 1],
    palette = colorRampPalette(c("blue", "white", "red"))
)
dev.off()


## Set scale the same
rng <- range(c(mats[, , 2] / mats[, , 1],
 mats2[, , 2] / mats2[, , 1]))
png(
    "plots/temp.png",
    width=5, height= 5, units='in', res = 300
)
plotMatrix(
    data = mats[, , 2] / mats[, , 1],
    palette = colorRampPalette(c("blue", "white", "red")),
    zrange = rng
)
dev.off()
png(
    "plots/temp.png",
    width = 5, height = 5, units = "in", res = 300
)
plotMatrix(
    data = mats2[, , 2] / mats2[, , 1],
    palette = colorRampPalette(c("blue", "white", "red")),
    zrange = rng
)
dev.off()

