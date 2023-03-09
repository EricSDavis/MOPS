## APA plots of all loops in the untreated
## condition in parental lines expressing
## CTCF degron +/- auxin to show what
## happens to existing CTCF loops

## Load required packages
library(data.table)
library(mariner)
library(plotgardener)

## Load data
iset <- readRDS("data/apaParentalVsNha9_CTCF.rds")

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
## Parental Control
png(
    "plots/apaNha9PareCtl.png",
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
    x = p$x + p$width + p$space * 2,
    width = p$space * 2,
    height = p$height * 0.5,
    fontcolor = "black",
    fontsize = 18
)
dev.off()

## Visualize results
## +CTCF
png(
    "plots/apaNha9CtcfCtl.png", 
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
    width = p$space * 2,
    height = p$height * 0.5,
    fontcolor = 'black',
    fontsize = 18
)
dev.off()

## -CTCF
png(
    "plots/apaNha9CtcfAux.png", 
    width = 5.75, 
    height = 5.5, 
    units = "in", 
    res = 300
)
pageCreate(6.5, 6, showGuides = FALSE)
apa <- plotMatrix(
    params = p,
    data = mats[, , 3]
)
annoHeatmapLegend(
    plot = apa,
    params = p,
    x = p$x + p$width + p$space*2,
    width = p$space * 2,
    height = p$height * 0.5,
    fontcolor = 'black',
    fontsize = 18
)
dev.off()
