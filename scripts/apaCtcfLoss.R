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

## Aggregate loops per hic files
mats <- aggHicMatrices(iset, by = "files")
mats <- mats / length(iset)

## Define plotting parameters
p <- pgParams(
    x = 0.65,
    y = 0.5,
    width = 5,
    height = 5,
    space = 0.1,
    zrange = c(0, round(max(mats)))
)

## Visualize results
## +CTCF
png(
    "plots/apaCtcfCtl.png", 
    width = 5.75, 
    height = 5.5, 
    units = "in", 
    res = 300
)
pageCreate(6.5, 6, showGuides = FALSE)
apa <- plotMatrix(
    params = p,
    data = mats[, , 1]
)
annoHeatmapLegend(
    plot = apa,
    params = p,
    x = p$x + p$width + p$space*2,
    width = p$space * 1.2,
    height = p$height * 0.5,
    fontcolor = 'black',
    fontsize = 14
)
dev.off()

## -CTCF
png(
    "plots/apaCtcfAux.png", 
    width = 5.75, 
    height = 5.5, 
    units = "in", 
    res = 300
)
pageCreate(6.5, 6, showGuides = FALSE)
apa <- plotMatrix(
    params = p,
    data = mats[, , 2]
)
annoHeatmapLegend(
    plot = apa,
    params = p,
    x = p$x + p$width + p$space*2,
    width = p$space * 1.2,
    height = p$height * 0.5,
    fontcolor = 'black',
    fontsize = 14
)
dev.off()
