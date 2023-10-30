library(ATACseqQC)
library(Rsamtools)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(ggplot2)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm39)
library(ChIPpeakAnno)
library(MotifDb)
library(GenomicAlignments)

# https://bioconductor.org/packages/release/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html
bamlist <- list.files("data/new_bam/", ".sorted.bam$", full.names = TRUE)
txdb_mm39 <- TxDb.Mmusculus.UCSC.mm39.refGene
seqlevelsStyle(txdb_mm39) <- "NCBI"
bsgenome_mm39 <- BSgenome.Mmusculus.UCSC.mm39
seqlevelsStyle(bsgenome_mm39) <- "NCBI"

for(bamfile in bamlist[4:8]){
  
  bamfile.labels <- gsub(".sorted.bam", "", basename(bamfile))
  dir.create(paste0("plots/005.1e/", bamfile.labels))

  print(Sys.time())
  if(!file.exists(paste0("plots/005.1e/", bamfile.labels, "/001_Library_complexity_", bamfile.labels, ".png"))){
    print(Sys.time())
    print(paste0("[1] Estimating Library complexity for sample ", bamfile.labels))
    png(paste0("plots/005.1e/", bamfile.labels, "/001_Library_complexity_", bamfile.labels, ".png"), h = 1500, w = 2000, res = 300)
    estimateLibComplexity(readsDupFreq(bamfile))
    dev.off()
  }

  if(!file.exists(paste0("plots/005.1e/", bamfile.labels, "/002_Fragment_size_distribution_", bamfile.labels, ".png"))){
    print(Sys.time())
    print(paste0("Plotting fragment size distribution for sample ", bamfile.labels))
    png(paste0("plots/005.1e/", bamfile.labels, "/002_Fragment_size_distribution_", bamfile.labels, ".png"), h = 1800, w = 2500, res = 300)
    fragSize <- fragSizeDist(bamfile, bamfile.labels)
    dev.off()
  }

  possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2",
                                  "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                  "TC", "UQ"),
                      "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                                    "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                                    "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                                    "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                                    "U2"))

  bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 1000),
                       param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
  tags <- names(bamTop100)[lengths(bamTop100)>0]
  tags

  ## files will be output into outPath
  outPath <- paste0("data/005.1e/", bamfile.labels)
  dir.create(outPath, showWarnings = F)
  ## shift the coordinates of 5'ends of alignments in the bam file

  ## if you don't have an available TxDb, please refer
  ## GenomicFeatures::makeTxDbFromGFF to create one from gff3 or gtf file.
  seqlev <- c(1:19, "X", "Y")

  seqinformation <- seqinfo(txdb_mm39)
  print(paste0("[2] Shifting coordinates of bam file for sample ", bamfile.labels))
  fname <- paste0(outPath, "shifted_alignment.rds")
  shiftedBamfile <- file.path(outPath, paste0(bamfile.labels, "_shifted.bam"))
  if(!file.exists(fname)){
    which <- as(seqinformation[seqlev], "GRanges")
    gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
    gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
    saveRDS(gal1, fname)
  }else{gal1 <- readRDS(fname)}

  # Promoters/Transcript score
  print(paste0("[2.1] Plotting PT score for sample ", bamfile.labels))
  txs <- transcripts(txdb_mm39)
  fname <- paste0(outPath, "pt_score.rds")
  if(!file.exists(fname)){
    pt <- PTscore(gal1, txs)
    saveRDS(pt, fname)
  }else{pt <- readRDS(fname)}
  
  pt_plot <- paste0("plots/005.1e/", bamfile.labels, "/003_PT_score_", bamfile.labels, ".png")
  if(!file.exists(pt_plot)){
    pt |>
      as.data.frame() |>
      filter(log2meanCoverage > -35) |>
    ggplot(aes(log2meanCoverage, PT_score)) +
      ggpointdensity::geom_pointdensity(alpha = 0.7) +
      labs(title = paste0("Promoter vs Transcript coverage ", bamfile.labels),
           x = "log2 mean coverage",
           y = "Promoter vs Transcripts score") +
      scale_color_viridis_c(option = "C")
    ggsave(pt_plot, h = 1500, w = 2500, units = "px")
  }

  if(!file.exists(paste0("plots/005.1e/", bamfile.labels, "/004_NFR_score_", bamfile.labels, ".png"))){
    print(Sys.time())
    print(paste0("[3] Plotting NFR score for sample ", bamfile.labels))
    nfr <- NFRscore(gal1, txs)
    nfr |> as.data.frame() |>
       filter(log2meanCoverage > -15) |>
    ggplot(aes(log2meanCoverage, NFR_score)) +
      ggpointdensity::geom_pointdensity(alpha = 0.7) +
      labs(title = paste0("NFRscore for 200bp flanking TSSs ", bamfile.labels),
           x = "log2 mean coverage",
           y = "Nucleosome Free Regions score") +
      scale_color_viridis_c(option = "C", end = 0.9)
    ggsave(paste0("plots/005.1e/", bamfile.labels, "/004_NFR_score_", bamfile.labels, ".png"), h = 1500, w = 2500, units = "px")
  }
  # Skipped TSS ES
  txs <- txs[seqnames(txs) %in% seqlev]
  ## split the reads into NucleosomeFree, mononucleosome,
  ## dinucleosome and trinucleosome.
  ## and save the binned alignments into bam files.

  fname <- paste0(outPath, "/split_obj.rds")
  if(!file.exists(fname)){
  print(Sys.time())
  print(paste0("[4] Splitting BAMs for sample ", bamfile.labels))
  objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=bsgenome_mm39, outPath = outPath)
  saveRDS(objs, fname)
  }else{objs <- readRDS(fname)}
  # Heatmap and coverage curve for nucleosome positions
  bamfiles <- file.path(outPath,
                        c("NucleosomeFree.bam",
                          "mononucleosome.bam",
                          "dinucleosome.bam",
                          "trinucleosome.bam"))

  TSS <- promoters(txs, upstream=0, downstream=1)
  TSS <- unique(TSS)
  ## estimate the library size for normalization
  librarySize <- estLibSize(bamfiles)

  ## calculate the signals around TSSs.
  print(paste0("[5] Plotting heatmap and coverage of NFR and mononucleosome-bound regions for sample ", bamfile.labels))
  NTILE <- 101
  dws <- ups <- 1010
  feature_hm <- paste0("plots/005.1e/", bamfile.labels, "/005_Heatmap_NFR_and_nucleosome_positions_", bamfile.labels, ".png")
  if(!file.exists(feature_hm)){
    sigs <- enrichedFragments(gal=objs[c("NucleosomeFree",
                                         "mononucleosome",
                                         "dinucleosome",
                                         "trinucleosome")],
                              TSS=TSS,
                              librarySize=librarySize,
                              seqlev=seqlev,
                              TSS.filter=0.5,
                              n.tile = NTILE,
                              upstream = ups,
                              downstream = dws)
    rm(objs)
    gc()
  
    ## log2 transformed signals
    sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
    #plot heatmap
    png(feature_hm, h = 2500, w = 1000, res = 300)
    featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                          zeroAt=.5, n.tile=NTILE, gap = 0.05, margin = c(0.1,0.05,0.25,0.05))
    dev.off()
  }
  ## get signals normalized for nucleosome-free and nucleosome-bound regions.
  fad_plot <- paste0("plots/005.1e/", bamfile.labels, "/006_Coverage_NFR_and_nucleosome_positions_", bamfile.labels, ".png")
  if(!file.exists(fad_plot)){
    png(fad_plot, h = 2000, w = 2000, res = 300)
    out <- featureAlignedDistribution(sigs,
                                      reCenterPeaks(TSS, width=ups+dws),
                                      zeroAt=.5, n.tile=NTILE, type="l",
                                      ylab="Averaged coverage",
                                      main = paste0("NFR and nucleosome-bound\nnormalized signal ", bamfile.labels))
    dev.off()

  }
  
  
  mat_plot <- paste0("plots/005.1e/", bamfile.labels, "/007_Scaled_coverage_NFR_and_nucleosome_positions_", bamfile.labels, ".png")
  if(!file.exists(mat_plot)){
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    out <- apply(out, 2, range01)
    png(mat_plot, h = 2000, w = 2000, res = 300)
    matplot(out, type="l", xaxt="n",
            xlab="Position (bp)",
            ylab="Fraction of signal",
            main = paste0("NFR and nucleosome-bound\nnormalized rescaled signal ", bamfile.labels))
    axis(1, at=seq(0, 100, by=10)+1,
         labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
    abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
    dev.off()
}
  # plot Footprints
  outPath <- paste0("data/005.1e/", bamfile.labels)
  shiftedBamfile <- file.path(outPath, paste0(bamfile.labels, "_shifted.bam"))
  seqlev <- c(1:19, "X", "Y")
  for (i in c("SREBF1", "SREBF2", "LEF1", "HLX")) {
    
    dir.create(paste0("plots/005.1e/", bamfile.labels), showWarnings = F)
    print(Sys.time())
    print(paste0("[6] Footprinting ", i, " for sample ", bamfile.labels))
    TF <- query(MotifDb, i)
    TF <- as.list(TF)
    print(TF[[1]], digits=2)
    plot_file <- paste0("plots/005.1e/", bamfile.labels, "/008_Footprint_coverage_", i, "_", bamfile.labels, ".png")
    if(!file.exists(plot_file)){
    png(plot_file, h = 2000, w = 2000, res = 300)
    sigs <- factorFootprints(shiftedBamfile, pfm=TF[[1]], 
                             genome=bsgenome_mm39, ## Don't have a genome? ask ?factorFootprints for help
                             min.score="90%", seqlev=seqlev,
                             upstream=100, downstream=100)
    dev.off()
    
    print(Sys.time())
    print(paste0("[7] Vplot ", i, " for sample ", bamfile.labels))
    png(paste0("plots/005.1e/", bamfile.labels, "/009_vplot_", i, "_", bamfile.labels, ".png"), h = 2000, w = 2000, res = 300)
    vp <- vPlot(shiftedBamfile, pfm=TF[[1]], 
                genome=bsgenome_mm39, min.score="90%", seqlev=seqlev,
                upstream=200, downstream=200, 
                ylim=c(30, 250), bandwidth=c(2, 1))
    dev.off()
    }else{next}
  }
  
  rm(list = grep("bamlist|txdb_mm39|bsgenome_mm39", ls(), value = T, invert = T))
  gc()
}

