# ENCODE ATAC peaks (https://rockefelleruniversity.github.io/RU_ATACseq/presentations/singlepage/RU_ATAC_part3.html#ENCODE_ATAC_data)
library(dplyr)
library(MotifDb)
library(BSgenome.Mmusculus.UCSC.mm39)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(SummarizedExperiment)
library(chromVAR)

# Select motifs to scan
opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["collection"]] <- "CORE"
opts[["all_versions"]] <- FALSE
motifsToScan <- getMatrixSet(JASPAR2020, opts)
#test <- getMatrixByName(JASPAR2020, "SREBF1")

myCounts <- readRDS("data/005.1c_Peaks_counts_non_redundant.rds")
peakRanges <- rowRanges(myCounts)
seqlevelsStyle(peakRanges) <- "UCSC"
peakRangesCentered <- resize(peakRanges, fix = "center", width = 100)
peakSeqs <- getSeq(BSgenome.Mmusculus.UCSC.mm39, peakRangesCentered)
names(peakSeqs) <- as.character(peakRangesCentered)
peakSeqs

motifHits <- matchMotifs(motifsToScan, peakSeqs, out = "matches")

mmMatrix <- motifMatches(motifHits)
dim(mmMatrix)

totalMotifOccurence <- colSums(mmMatrix)
totalMotifOccurence[1:4]

# Summarizing ATAC signal to motifs
myCounts <- myCounts[rowSums(assay(myCounts)) > 5, ]
seqlevelsStyle(myCounts) <- "UCSC"
myCounts <- addGCBias(myCounts, genome = BSgenome.Mmusculus.UCSC.mm39)

motif_ix <- matchMotifs(motifsToScan, myCounts, genome = BSgenome.Mmusculus.UCSC.mm39)

dir.create("data/005.1d")
fname <- "data/005.1d/deviations_in_motifs.rds"
if(!file.exists(fname)){
  deviations <- computeDeviations(object = myCounts, annotations = motif_ix)
  saveRDS(deviations, fname)
}else{deviations <- readRDS(fname)}
fname <- "data/005.1d/Peaks_variability_in_motifs.rds"
if(!file.exists(fname)){
  variability_Known <- computeVariability(deviations)
  saveRDS(variability_Known, fname)
}else{variability_Known <- readRDS(fname)}

devZscores <- deviationScores(deviations)
devZscores[1:2, ]

variability_Known <- variability_Known[order(variability_Known$p_value), ]
variability_Known[1:10, ]

topVariable <- variability_Known[1:50, ]
devTop <- merge(topVariable[, 1, drop = FALSE], devZscores, by = 0)
devTop[1:50, ]

dev <- merge(variability_Known, devZscores, by = 0)
dev <- dev[dev$p_value_adj < 0.05,]
dev <- dev[order(dev$p_value_adj),]

png("plots/005.1d/000_Variability_of_motifs_in_peaks.png", h = 1200, w = 2400, res = 300)
plotVariability(dev, n = 0, use_plotly = T)
dev.off()

# Plot chrmovar results
library(pheatmap)
dir.create("plots/005.1d")
n <- which(dev$name == "BACH2")
devToPlot <- as.matrix(dev[1:n, -c(1:7)])
rownames(devToPlot) <- dev$name[1:n]
colnames(devToPlot) <- c("wt 1", "siAgc1 1", "wt 2", "wt 3", "wt 4", "siAgc1 2", "siAgc1 3", "siAgc1 4")
png("plots/005.1d/001_Binding_motifs_variability.png", h = 3000, w = 3000, res = 300)
pheatmap(devToPlot, main = "Variability of binding motifs across samples (Z-scores)", angle_col = 45)
dev.off()

# De novo motif discovery
library(DESeq2)
atacDDS <- readRDS("data/005.1c_ATAC_deseq_obj.rds")
myRes <- results(atacDDS, contrast = c("Group", "kd", "wt"), format = "GRanges")
seqlevelsStyle(myRes) <- "UCSC"
myRes <- myRes[order(myRes$padj), ]
upRegions <- myRes[myRes$log2FoldChange > 0 & myRes$padj < 0.05 & !is.na(myRes$padj),]
downRegions <- myRes[myRes$log2FoldChange < 0 & myRes$padj < 0.05 & !is.na(myRes$padj), ]
upRegions

upRegions <- resize(upRegions, fix = "center", width = 100)
downRegions <- resize(downRegions, fix = "center", width = 100)

upStrings <- getSeq(BSgenome.Mmusculus.UCSC.mm39, upRegions)
downStrings <- getSeq(BSgenome.Mmusculus.UCSC.mm39, downRegions)
names(upStrings) <- as.character(upRegions)
names(downStrings) <- as.character(downRegions)
writeXStringSet(upStrings, file = "data/005.1d_UpRegions.fa")
writeXStringSet(downStrings, file = "data/005.1d_DownRegions.fa")

# submit to https://meme-suite.org/meme/tools/meme-chip and use differential enrichment mode, use up as primary sequence
# input motifs = Eukaryote DNA, Vertebrates. Download tar and extract
library(universalmotif)
memeMotifs <- read_meme("data/appMEMECHIP_5.4.116502600188291740058233/combined.meme")
memeMotifs

memeMotifsTFBStools <- convert_motifs(memeMotifs, "TFBSTools-PWMatrix")
memeMotifsTFBStools

