## Example of a loop that is lost
## after degredation of CTCF

## Load required libraries
library(data.table)
library(mariner)
library(plotgardener)
library(InteractionSet)

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

## Read in loops & expand to 10Kb
loops <- 
    fread(loopFiles) |>
    as_ginteractions()

## Finding candidate loops
## Pick a loop with the best change in
## enrichment. (It happens to be 1)
# enrich <- calcLoopEnrichment(
#     x = loops,
#     files = hicFiles
# )
# which(rank(enrich[,1] / enrich[,2]) == 1)


## Visualize loops
## Define shared parameters
i <- 1
p <- pgParams(

    ## Hi-C params
    chrom = seqnames1(loops)[i],
    chromstart = start1(loops)[i] - 150e3,
    chromend = end2(loops)[i] + 150e3,

    ## Default position
    x = 0.5,
    y = 0.5,
    width = 11,
    height = 5,
    space = 0.1,

    ## Hi-C params
    resolution = 25e3,
    zrange = c(0, 150)
)

png("plots/ctcfLossExampleCtl.png", width=12, height=5.5, unit='in', res=300)
pageCreate(12.5, 6, showGuides = FALSE)
hicPlot <- plotHicRectangle(
    params = p,
    data = hicFiles[1]
)
annoHeatmapLegend(
    params = p,
    plot = hicPlot,
    x = p$x + p$width + p$space*2,
    width = p$space * 1.2,
    height = p$height * 0.5,
    fontcolor = 'black',
    fontsize = 14
)
annoGenomeLabel(
    plot = p,
    x = p$x,
    y = p$y + p$height,
    fontsize = 16
)
dev.off()


png("plots/ctcfLossExampleAux.png", width=12, height=5.5, unit='in', res=300)
pageCreate(12.5, 6, showGuides = FALSE)
hicPlot <- plotHicRectangle(
    params = p,
    data = hicFiles[2]
)
annoHeatmapLegend(
    params = p,
    plot = hicPlot,
    x = p$x + p$width + p$space*2,
    width = p$space * 1.2,
    height = p$height * 0.5,
    fontcolor = 'black',
    fontsize = 14
)
annoGenomeLabel(
    plot = p,
    x = p$x,
    y = p$y + p$height,
    fontsize = 16
)
dev.off()
