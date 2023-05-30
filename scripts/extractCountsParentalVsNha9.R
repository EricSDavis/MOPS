## Load required libraries
library(mariner)
library(SummarizedExperiment)

## List Hi-C files
hicFiles <- list.files(
    path = "data/raw/hic/bioreps",
    pattern = "Control",
    full.names = TRUE
)

## Load merged loops
loops <- readRDS("data/mergedLoops/mergedLoopsControl.rds")

## Bin to 10Kb
loops <- binPairs(loops, 10e3)

## Extract replicate Hi-C counts for each .hic file
loopCounts <- pullHicPixels(
    x = loops,
    binSize = 10e3,
    files = hicFiles,
    h5File = "data/ControlMergedLoopCounts.h5",
    norm = "NONE"
)

## Add md5sums (optional)
colData(loopCounts)$md5 <- tools::md5sum(colData(loopCounts)$files)

## Save results
saveRDS(loopCounts, file = "data/ControlMergedLoopCounts.rds")
