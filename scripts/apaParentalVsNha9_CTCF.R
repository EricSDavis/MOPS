## APA plots of NHA9 loops before
## and after degradation of CTCF

## Load required packages
library(data.table)
library(InteractionSet)
library(mariner)
library(plotgardener)

## Load differential looping data
loopCounts <- readRDS("data/ControlDiffLoopCounts.rds")

## Isolate diff loop results
loops <- interactions(loopCounts)
diffLoops <- loops[which(loops$padj <= 0.1)]

## Load hic file paths
hicFiles <- c(
    "data/raw/hic/condition/MOPS_HCT_CTCFNHA9_Control_0h_inter_30.hic",
    "data/raw/hic/condition/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30.hic"
)

## Set the binSize & buffer
binSize <- 10e3
buffer <- 10

## Ensure loops are the correct binSize
if (unique(unlist(width(diffLoops))) != (binSize+1)) {
    diffLoops <- binPairs(diffLoops, binSize)
}

## Temporary work-around to work with
## hictoolsr function
## Filter out short loops
diffLoops <- diffLoops |>
    as.data.table() |>
    as_ginteractions() |>
    hictoolsr::filterBedpe(
        res = binSize,
        buffer = buffer
    )

## Extract pixels for each Hi-C file
iset <- diffLoops |>
    pixelsToMatrices(buffer = buffer) |>
    pullHicMatrices(
        files = hicFiles,
        binSize = binSize,
        h5File = "data/apaParentalVsNha9_CTCF.h5",
        norm = "KR",
        blockSize = 6e6
    )

## Save iset
saveRDS(iset, file = "data/apaParentalVsNha9_CTCF.rds")
