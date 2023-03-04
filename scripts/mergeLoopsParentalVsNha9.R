## Merge loop calls among all control
## (non auxin treated cells)

## Load required packages
library(data.table)
library(mariner)

## File paths for all loops called
allLoopFiles <- list.files(
    path = "data/raw/sip_loops",
    pattern = "5kbLoops.txt",
    full.names = TRUE,
    recursive = TRUE
)

## Isolate loops not treated with auxin
## (i.e. Control)
loopFiles <- grep("Control", allLoopFiles, value=TRUE)

## Read in loops
loops <- 
  lapply(loopFiles, fread) |>
  lapply(as_ginteractions) |>
  setNames(gsub("^.*loops/(.*)/5kb.*$", "\\1", loopFiles))

## Count total number of loops in each set:
lapply(loops, length)

## Merge loops
mergedLoops <- mergePairs(
    x = loops,
    radius = 10e3,
    column = "APScoreAvg"
)

saveRDS(
    object = mergedLoops,
    file = "data/mergedLoops/mergedLoopsControl.rds"
)