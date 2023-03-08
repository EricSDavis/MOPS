## Merge loops between parental lines
## expressing CTCF degron +/- auxin

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

## Isolate loops containin CTCF degron
## Excluding NHA9 to avoid confounding
## effects of NHA9 binding.
loopFiles <- grep("CTCFparen", allLoopFiles, value = TRUE)

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

## Interesting stats about the loops:
## Number of loops that are shared +/- auxin
## treatment...
subsetBySource(
    x = mergedLoops,
    include = sources(mergedLoops)
) |> length()

## Save result
saveRDS(
    object = mergedLoops,
    file = "data/mergedLoops/mergedLoopsParentalCTCF.rds"
)