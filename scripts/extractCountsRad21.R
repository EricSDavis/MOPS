## Load required libraries
library(mariner)
library(SummarizedExperiment)

## List Hi-C files
hicFiles <- list.files(
    path = "data/raw/hic/bioreps",
    pattern = "RAD21NHA9",
    full.names = TRUE
)

## Load merged loops
loops <- readRDS("data/mergedLoops/rad21Nha9MergedLoops.rds")

## Bin to 10Kb
loops <- binPairs(loops, 10e3)

## Extract replicate Hi-C counts for each .hic file
loopCounts <- pullHicPixels(
    x = loops,
    binSize = 10e3,
    files = hicFiles,
    h5File = "data/rad21Nha9MergedLoopCounts.h5",
    norm = "NONE"
)

## Add md5sums (optional)
colData(loopCounts)$md5 <- tools::md5sum(colData(loopCounts)$files)

## Save results
saveRDS(loopCounts, file = "data/rad21Nha9MergedLoopCounts.rds")
