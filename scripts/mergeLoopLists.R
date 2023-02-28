## Merge loop calls from various conditions

## Load required packages
library(data.table)
library(mariner)

## TODO: Consider merging all Hi-C data and calling
##       loops in the "omega" dataset.
## File paths for loops
loopFiles <- list.files(path="data/raw/sip_loops",
                        pattern="5kbLoops.txt",
                        full.names=TRUE,
                        recursive=TRUE)

## Read in loops
loops <- 
  lapply(loopFiles, fread) |>
  lapply(as_ginteractions) |>
  setNames(gsub("^.*loops/(.*)/5kb.*$", "\\1", loopFiles))


## Count total number of loops in each set:
lapply(loops, length)

## Merge loops by condition ----------------------------------------------------

## Define sets to merge
sets <- 
  list(
    "ctcfNhad9" = loops[1:2],
    "ctcfParental" = loops[3:4],
    "rad21Nha9" = loops[5:6],
    "rad21Parental" = loops[7:8]
    )

## Apply loop merging across sets
mergedLoops <- lapply(sets, mergePairs, radius=10e3, column="APScoreAvg")

## Save merged loop lists
lapply(seq_along(mergedLoops), \(i) {
  saveRDS(object=mergedLoops[[i]],
          file=paste0("data/mergedLoops/",
                      names(mergedLoops[i]),
                      "MergedLoops.rds"))
}) |> invisible()


## Merge loops from all conditions ---------------------------------------------

allMergedLoops <- mergePairs(x=loops,
                             radius=10e3,
                             column="APScoreAvg")

saveRDS(object=allMergedLoops, file="data/mergedLoops/allMergedLoops.rds")
