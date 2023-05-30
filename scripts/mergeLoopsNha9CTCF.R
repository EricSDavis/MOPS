## Merge loops of NHA9 containing
## CTCF degron cells (+/- auxin)

## Load required packages
library(data.table)
library(mariner)

## File paths for NHA9 CTCF degron loops
loopFiles <- c(
    "data/raw/sip_loops/MOPS_HCT_CTCFNHA9_Control_0h_inter_30/5kbLoops.txt",
    "data/raw/sip_loops/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30/5kbLoops.txt"
)

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

## Save results
saveRDS(
    object = mergedLoops,
    file = "data/mergedLoops/mergedLoopsNha9Ctcf.rds"
)