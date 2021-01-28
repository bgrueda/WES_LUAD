setwd("~/Desktop/Bioinformatic_tools/hg38")

library(maftools)
library("viridis")

xann = read.maf(maf = "xannotation.maf")
xann
getSampleSummary(xann)
getGeneSummary(xann)
getFields(xann)
write.mafSummary(maf = xann, basename = 'xann')

vc_cols = scale_color_viridis(discrete = TRUE, option = "D")

plotmafSummary(maf = xann, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE + scale_color_viridis(discrete = TRUE, option = "D"))

#oncoplot for top ten mutated genes.
oncoplot(maf = xann, top = 10)

xann.titv = titv(maf = xann, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = xann.titv)

#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
lollipopPlot(maf = laml, gene = 'DNMT3A', AACol = 'Protein_Change', showMutationRate = TRUE)
lollipopPlot(maf = laml, gene = 'DNMT3A', showDomainLabel = FALSE, labelPos = 882)

brca <- system.file("extdata", "brca.maf.gz", package = "maftools")
brca = read.maf(maf = brca, verbose = FALSE)

rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.4)