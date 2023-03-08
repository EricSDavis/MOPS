## APA plots of all loops in the untreated
## condition in parental lines expressing
## CTCF degron +/- auxin to show what
## happens to existing CTCF loops

## Load required packages
library(data.table)
library(mariner)
library(plotgardener)

## Load loop file paths
loopFiles <- file.path(
    "data/raw/sip_loops",
    "MOPS_HCT_CTCFparental_Control_0h_inter_30",
    "5kbLoops.txt"
)

## Load hic file paths
hicFiles <- list.files(
    path = "data/raw/hic/condition",
    pattern = "CTCFparental",
    recursive = TRUE,
    full.names = TRUE
) |> rev() # reverse so that control is first

## Set the binSize & buffer
binSize <- 10e3
buffer <- 10

## Read in loops & expand to 10Kb
loops <-
    fread(loopFiles) |>
    binPairs(binSize = binSize)

## Filter out short loops
loops <- hictoolsr::filterBedpe(
    bedpe = loops,
    res = binSize,
    buffer = buffer
)

## Extract pixels for each Hi-C file
iset <- loops |>
    pixelsToMatrices(buffer = buffer) |>
    pullHicMatrices(
        files = hicFiles,
        binSize = binSize,
        h5File = "data/ctcfParentalCounts.h5",
        norm = "KR"
    )

## Save iset
saveRDS(iset, file="data/ctcfParentalCounts.rds")
