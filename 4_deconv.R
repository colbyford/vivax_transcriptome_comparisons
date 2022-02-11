# devtools::install_github("ctlab/linseed")

library(linseed)
lo <- LinseedObject$new("GSE19830", samples=10:42, topGenes=10000)

lo$calculatePairwiseLinearity()
lo$calculateSpearmanCorrelation()
lo$calculateSignificanceLevel(100)
lo$significancePlot(0.01)

lo$filterDatasetByPval(0.01)

lo$svdPlot()

lo$setCellTypeNumber(3)
lo$project("full") # projecting full dataset
lo$projectionPlot(color="filtered")


lo$project("filtered")
lo$smartSearchCorners(dataset="filtered", error="norm")


lo$deconvolveByEndpoints()
plotProportions(lo$proportions)


# lets select 100 genes closest to the simplex corners 
lo$selectGenes(100)
lo$tsnePlot()



data("proportionsLiverBrainLung")
dotPlotPropotions(lo$proportions, proportionsLiverBrainLung[, 10:42], guess=TRUE)