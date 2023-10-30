library(tidyverse)
library(Rsubread)
library(Rsamtools)
library(GenomicAlignments)
library(soGGi)
library(rtracklayer)
library(ChIPQC)
library(DT)
library(ChIPseeker)
library(rGREAT)
library(BSgenome.Mmusculus.UCSC.mm39)
library(TxDb.Mmusculus.UCSC.mm39.refGene) # Guide uses knownGene, but as of 03Feb22 only refGene is available for mm39

# Guide followed: https://rockefelleruniversity.github.io/RU_ATAC_Workshop.html
# Check also https://rockefelleruniversity.github.io/RU_ATACseq/
# Download reference genome from (https://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz)
genome <- "data/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
rsub_ind <- gsub("\\.fa.gz$", "", genome)

# Build index from reference genome, sort and align reads
buildindex(rsub_ind, genome)

fq <- c(list.files("data/new_fastq", pattern = "R1_001.fastq.gz", full.names = TRUE),
        list.files("data/1Fastq", pattern = "R1_001.fastq.gz", full.names = TRUE))

# We used a reference genome that has NCBI annotations, so IT IS FUNDAMENTAL to pay attention when loading
# objects that MIGHT have different annotations. If they are different, they MUST be changed to match those of 
# the reference genome. In order to change the annotation type just set the seqlevelsStyle like shown below. 
# It is possible to check the levels using seqlevels(obj) and assess if they match across the different objects.
# For reference if the chromsomea are only numbers ("1", "2", etc) they are NCBI formattted
# "chr1", "chr2" etc are UCSC
TSSs <- resize(genes(TxDb.Mmusculus.UCSC.mm39.refGene), fix = "start", 1) # Need this inside the loop
seqlevelsStyle(TSSs) <- c("NCBI")

for(rname in fq){
  
  outBAM <- paste0("data/new_bam/", 
                   basename(rname) %>% str_extract("^[[:alnum:]]*_[[:alnum:]]*"),
                   ".bam")
  align(rsub_ind,
        readfile1 = rname,
        readfile2 = gsub("_R1_", "_R2_", rname),
        output_file = outBAM,
        nthreads = 5,
        type = 1,
        unique = TRUE,
        maxFragLength = 2000)

  sortedBAM <- str_replace(outBAM, "\\.bam", ".sorted.bam")
  sortBam(outBAM, gsub("\\.bam", "", sortedBAM))

  indexBam(sortedBAM)

   sample <- basename(sortedBAM) %>% 
     str_extract("[:alnum:]{2,3}_[:alnum:]{2,3}")
   plt_title <- paste0("ATAC-Seq Mapped Reads for sample ", sample)

  idxstatsBam(sortedBAM) %>%
    filter(seqnames %in% c(1:19, "X", "Y")) %>%
    ggplot(aes(seqnames, mapped, fill = seqnames)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    coord_flip() +
    labs(title = plt_title,
         x = "Chromosome",
         y = "Mapped Reads") +
    scale_fill_viridis_d() +
    theme_light()

  ggsave(paste0("plots/005.1_ATACSeq/005.1_1_Mapped_reads/005.1_1_Mapped_reads_", sample, ".png"),
         dpi = 600,
         height = 3.5, width = 6)

# Read in mapped reads

# bampar <- scanBam(sortedBAM)
bampar[[1]]$mapq %>%
  as_data_frame() %>%
  ggplot(aes(x = value)) +
  geom_bar() +
  stat_bin(aes(label=..count..), vjust=0, geom="text", position="identity")
sum(is.na(bampar[[1]]$mapq)) # 3362019

unique(bampar[[1]]$mapq)  # 40 20 13  8 10  6 NA

bampar[[1]]$mapq %>%
  as_data_frame() %>%
  group_by(value) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(perc = n/sum(n))

# names(atacReads[[1]])

  atacReads <- readGAlignmentPairs(sortedBAM, 
                                   param = ScanBamParam(mapqFilter = 20, 
                                                        what = c("qname", "mapq", "isize"),
                                                        flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE)))
  
  atacReads <- atacReads[seqnames(atacReads) %in% c(as.character(seq(1:19)), "X", "Y")]
  atacReads_read1 <- GenomicAlignments::first(atacReads)
  insertSizes <- abs(elementMetadata(atacReads_read1)$isize)

  # Inserts and higlight on nucleosome zones ----
  insert_plot <- table(insertSizes) %>%
    data.frame %>%
    rename(InsertSize = insertSizes, Count = Freq) %>%
    mutate(InsertSize = as.numeric(as.vector(InsertSize)),
           Count = as.numeric(as.vector(Count))) %>%
    filter(Count > 1)

  ggplot() +
    geom_line(data = insert_plot, aes(x = InsertSize, y = Count), size = 0.6) +
    geom_rect(aes(xmin = 0, xmax = 100, ymin = 0, ymax = Inf, fill = "Nucleosome-free\n(0-100)"), colour = NA, alpha = 0.3) +
    geom_rect(aes(xmin = 180, xmax = 247, ymin = 0, ymax = Inf, fill = "Mono-Nucleosome\n(180-247)"), colour = NA, alpha = 0.3) +
    geom_rect(aes(xmin = 315, xmax = 437, ymin = 0, ymax = Inf, fill = "Di-nucleosome\n(315-437)"), colour = NA, alpha = 0.3) +
    labs(title = paste0("Inserts size distribution in sample ", sample),
         color = "Legend") +
    scale_y_log10() +
    scale_fill_manual('Inserts region',
                      values = c("blue", "red", "green"),
                      guide = guide_legend(override.aes = list(alpha = 0.25))) +
    theme_light() +
    theme(legend.spacing.y = unit(0.2, 'cm'),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8)) +
    guides(fill = guide_legend(byrow = TRUE))

  ggsave(paste0("plots/005.1_ATACSeq/005.1_2_Insert_size_distribution/005.1_2_Insert_size_distribution_highlight_", sample, ".png"),
         h = 2500,
         w = 4200,
         units = "px",
         dpi = 600)

  
  # Get coverage + GC content and plots ----
  ## Run this first part on Linux
  ## cd $agc1v2/data/
  ## for rname in Sorted*.bam
  ## do
  ## output=`basename $rname .bam`
  ## bamCoverage -b $rname -o ${output}.bedgraph -of bedgraph --numberOfProcessors 4 --binSize 50000 --minMappingQuality 20 --maxFragmentLength 2000
  ## done

  # for(file in list.files("data/", pattern = "*.bedgraph", full.names = TRUE)){
  # 
  #   bed <- read_delim(file, col_names = c("Chromosome", "Start", "End", "Count"))
  #   sample <- str_extract(file, "_[:alnum:]{2}_[:alnum:]{3}.") %>%
  #     str_replace_all("[:punct:]", " ")
  # 
  #   bedfil <- bed %>%
  #     mutate(Region = rowMeans(select(bed, Start, End), na.rm = T), .after = Chromosome) %>%
  #     filter(Chromosome %in% c(seq(1:19), "X"), Count != 0) %>%
  #     rowwise() %>%
  #     mutate(Chromosome = paste0("chr", Chromosome),
  #            Chromosome = factor(Chromosome, levels = paste0("chr", c(seq(1:19), "X"))))
  # 
  #   bedgr <- bedfil %>%
  #     select(Chromosome, Start, End) %>%
  #     makeGRangesFromDataFrame()
  # 
  #   freqs <- alphabetFrequency(getSeq(BSgenome.Mmusculus.UCSC.mm39, bedgr))
  #   gc <- (freqs[,'C'] + freqs[,'G'])/rowSums(freqs)
  # 
  #   count_me <- mean(bedfil$Count)
  #   count_sd <- sd(bedfil$Count)
  #   gc_me <- mean(gc)
  #   gc_sd <- sd(gc)
  # 
  #   bedgc <- bedfil %>%
  #     bind_cols(gc) %>%
  #     rename(GC_content = last_col()) %>%
  #     mutate(count_zscore = (Count - count_me)/count_sd,
  #            gc_zscore = (GC_content - gc_me)/gc_sd) %>% # Compute Z score of Counts to remove outliers
  #     filter(abs(count_zscore) < 3, abs(gc_zscore) < 3) %>%
  #     group_by(Chromosome) %>%
  #     mutate(p_val = cor.test(log10(Count), GC_content)$p.value,
  #            p_val = replace(p_val, p_val == 0e+00, 0 + .Machine$double.eps),
  #            p_adj = p.adjust(p_val, method = "bonferroni") %>% scales::scientific(),
  #            adj_rsq = summary(lm(Count ~ GC_content))$adj.r.squared %>% signif(3),
  #            chr_label = paste0(Chromosome, " - ",
  #                           "p_adj = ", p_adj, " - ",
  #                           "adj_rsq = ", adj_rsq)) %>%
  #     ungroup()
  # 
  #   ggplot(bedgc, aes(x = Region)) +
  #     geom_point(aes(y = Count/mean(Count)), size = 1, color = "black", alpha = 0.2) +
  #     geom_point(aes(y = GC_content/mean(GC_content)), size = 1, color = "#ED254E", alpha = 0.1)  +
  #     scale_y_continuous(name = "Normalized Counts",
  #                        sec.axis = sec_axis(trans = ~., name = "Normalized GC content"),
  #                        trans = "log10") +
  #     facet_wrap(~Chromosome, scales = "free") +
  #     theme_light() +
  #     labs(title = paste0("Counts and GC content distribution in sample", sample)) +
  #     theme(plot.title = element_text(size = 25),
  #           axis.title.x = element_text(size = 18),
  #           axis.title.y.left = element_text(color = "black", size = 18),
  #           axis.text.y.left = element_text(color = "black"),
  #           axis.title.y.right = element_text(color = "#ED254E", size = 18),
  #           axis.text.y.right = element_text(color = "#ED254E"),
  #           strip.text.x = element_text(color = "black", size = 12),
  #           strip.background = element_rect(fill = "#ED254E")
  #     )
  # 
  #   ggsave(paste0("plots/005_ATACSeq/005.1_3_Counts_and_GC_content_over_region",
  #                 str_replace_all(sample, "\\s", "_"),
  #                 ".png"),
  #          h = 35,
  #          w = 70,
  #          units = "cm",
  #          dpi = 600)
  # 
  #   ggplot(bedgc, aes(GC_content, Count)) +
  #     geom_point(alpha = 0.3) +
  #     geom_smooth(method = "lm", formula = y ~ x, color = "#4775FF") +
  #     scale_y_log10() +
  #     facet_wrap(~chr_label, scales = "free") +
  #     theme_light() +
  #     labs(title = paste0("Counts distribution over GC content in sample", sample)) +
  #     theme(plot.title = element_text(size = 25),
  #           axis.title.x = element_text(size = 18),
  #           axis.title.y = element_text(size = 18),
  #           strip.background = element_rect(fill = "#4775FF"),
  #           strip.text = element_text(size = 12))
  # 
  # 
  #     ggsave(paste0("plots/005_ATACSeq/005.1_4_Counts_over_GC_content",
  #                   str_replace_all(sample, "\\s", "_"),
  #                   ".png"),
  #            h = 35,
  #            w = 70,
  #            units = "cm",
  #            dpi = 600)
  # }
  
  # ATACSeq signal of TSS ----
  nucFree <- regionPlot(bamFile = sortedBAM, testRanges = TSSs, style = "point",
                        format = "bam", paired = TRUE, minFragmentLength = 0, maxFragmentLength = 100,
                        forceFragment = 50)

  # Mononucleosome
  monoNuc <- regionPlot(bamFile = sortedBAM, testRanges = TSSs, style = "point",
                        format = "bam", paired = TRUE, minFragmentLength = 180, maxFragmentLength = 240,
                        forceFragment = 80)

  # Dinucleosome
  diNuc <- regionPlot(bamFile = sortedBAM, testRanges = TSSs, style = "point",
                      format = "bam", paired = TRUE, minFragmentLength = 315, maxFragmentLength = 437,
                      forceFragment = 160)

  # Plot regions
  plotRegion(nucFree) +
    labs(title = paste0("Signal from nucleosome free regions in sample ", sample),
                             x = "Relative position from TSS") +
    theme_light() +
    theme(plot.title = element_text(size = 17),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          strip.text.y = element_blank(),
          strip.text.x = element_blank(),
          strip.background = element_blank())
  ggsave(paste0("plots/005.1_ATACSeq/005.1_5_NFR_signal/005.1_5_NFR_signal_", sample, ".png"), h = 3500, w = 4100, units = "px" ,dpi = 600)
  plotRegion(monoNuc) + labs(title = paste0("Signal from mononucleosome regions in sample ", sample),
                             x = "Relative position from TSS") +
    theme_light() +
    theme(plot.title = element_text(size = 17),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          strip.text.y = element_blank(),
          strip.text.x = element_blank(),
          strip.background = element_blank())
  ggsave(paste0("plots/005.1_ATACSeq/005.1_6_mononucleosome_signal/005.1_6_mononucleosome_signal_", sample, ".png"), h = 3500, w = 4000, units = "px",dpi = 600)
  plotRegion(diNuc) + labs(title = paste0("Signal from dinucleosome regions in sample ", sample),
                           x = "Relative position from TSS") +
    theme_light() +
    theme(plot.title = element_text(size = 17),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          strip.text.y = element_blank(),
          strip.text.x = element_blank(),
          strip.background = element_blank())
  ggsave(paste0("plots/005.1_ATACSeq/005.1_7_dinucleosome_signal/005.1_7_dinucleosome_signal_", sample, ".png"), h = 3500, w = 4000, units = "px",dpi = 600)

  # Subsetting ATAC-seq reads files by insert sizes
  atacReads_Open <- atacReads[insertSizes < 100, ]
  saveRDS(atacReads_Open, gsub("\\.bam", "_openRegions\\.rds", sortedBAM))
  atacReads_MonoNuc <- atacReads[insertSizes > 180 & insertSizes < 240, ]
  monoNucBam <- gsub("\\.bam", "_monoNuc\\.rds", sortedBAM)
  atacReads_diNuc <- atacReads[insertSizes > 315 & insertSizes < 437, ]
  diNucBam <- gsub("\\.bam", "_diNuc\\.rds", sortedBAM)

  # Create Bam files split by insert sizes
  openRegionBam <- gsub("\\.bam", "_openRegions\\.bam", sortedBAM)
  monoNucBam <- gsub("\\.bam", "_monoNuc\\.bam", sortedBAM)
  diNucBam <- gsub("\\.bam", "_diNuc\\.bam", sortedBAM)

  export(atacReads_Open, openRegionBam, format = "bam")
  export(atacReads_MonoNuc, monoNucBam, format = "bam")
  export(atacReads_diNuc, diNucBam, format = "bam")

  # Creating files split by insert sizes
  openRegionBigWig <- gsub("\\.bam", "_openRegions\\.bw", sortedBAM)
  atacFragments_Open <- granges(atacReads_Open)
  export.bw(coverage(atacFragments_Open), openRegionBigWig)

  monoNucbw <- gsub("\\.bam", "_monoNuc\\.bw", sortedBAM)
  atacFragments_monoNuc <- granges(atacReads_MonoNuc)
  export.bw(coverage(atacFragments_monoNuc), monoNucbw)

  diNucbw <- gsub("\\.bam", "_diNuc\\.bw", sortedBAM)
  atacFragments_diNuc <- granges(atacReads_diNuc)
  export.bw(coverage(atacFragments_diNuc), diNucbw)
  
  rm(list = ls(pattern = "atacReads*")) # Remove elements in order not to occupy too much RAM
  gc()
}
  
# Run macs2 (005.1_macs2_callpeak.sh)
# #Download blacklisted regions from (https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz)
blkList <- import.bed("data/mm10_blacklist_ENCFF547MET.bed.gz") 
# Since we aligned on mm39, we need to liftover this bed (downlaod chain file from https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz)
# liftover guide: https://genviz.org/module-01-intro/0001/06/02/liftoverTools/
# custom annotation guide: https://bioconductor.org/help/course-materials/2014/BioC2014/Bioc2014_ChIPQC_Practical.pdf
chain_39to10 <- import.chain("data/mm39ToMm10.over.chain")
chain_10to39 <- import.chain("data/mm10ToMm39.over.chain")
blklist_lift <- liftOver(blkList, chain_10to39)
seqlevelsStyle(blklist_lift) <- "NCBI"
export(blklist_lift, "data/blacklist_mm39.bed", format = "bed")
bl_bed <- import.bed("data/blacklist_mm39.bed")

peaks <- list.files("data/new_bam", pattern = ".narrowPeak", full.names = T)
bams <- list.files("data/new_bam", pattern = "openRegions.bam$", full.names = T)
blacklisted <- data.frame(Sample = character(), Blacklisted = double(), Not_Blacklisted = double())

mm39_ncbi <- TxDb.Mmusculus.UCSC.mm39.refGene
seqlevelsStyle(mm39_ncbi) <- "NCBI" 
genes_mm39 <- genes(mm39_ncbi)
mm39_ann <- list(version = "custom_mm39", genes_mm39)
chain_39to10 <- import.chain("data/mm39ToMm10.over.chain")

macs_stat <- vector(mode = "list", length = 8L)
great_list <- vector(mode = "list", length = 8L)

for(i in 1:length(peaks)){
  
  peak <- peaks[i]
  bam <- bams[i]
  sample <- str_extract(basename(peak), "[:alnum:]{2,3}_[:alnum:]{2,3}")
  
  qcRes <- ChIPQCsample(bam, peaks = peak, annotation = mm39_ann,
                        blacklist = bl_bed, chromosomes = c(as.character(1:19)))

  QCmetrics(qcRes) %>% t %>% data.frame %>%
    dplyr:::select(Reads, starts_with(c("Filt")), starts_with(c("RiP")), starts_with(c("RiBL"))) %>%
    datatable(rownames = NULL)

  flagtagcounts(qcRes) %>% t %>% data.frame %>%
    mutate(Dup_Percent = (DuplicateByChIPQC/Mapped) * 100) %>%
    dplyr:::select(Mapped, Dup_Percent) %>% datatable(rownames = NULL)

  MacsCalls <- granges(qcRes)

  blacklisted <- rbind(blacklisted,
                       data.frame(
                         Sample = sample,
                         Blacklisted = sum(MacsCalls %over% bl_bed),
                         Not_Blacklisted = sum(!MacsCalls %over% bl_bed)))

  macs_filt <- MacsCalls[!MacsCalls %over% bl_bed]
  macs_anno <- annotatePeak(macs_filt, TxDb = mm39_ncbi) #Errors are only regarding scaffolds
  
  macs_gr <- macs_filt
  seqlevelsStyle(macs_gr) <- "UCSC"
  macs_lift <- liftOver(macs_gr, chain_39to10)
  macs_lift <- unlist(macs_lift)
  great_Job <- submitGreatJob(macs_lift, species = "mm10", version = "3.0.0")
  great_ResultTable = getEnrichmentTables(great_Job, category = c("GO", "Pathway Data"))
  great_list[[i]] <- map_dfr(names(great_ResultTable), function(x){
    great_ResultTable[[x]] %>%
      filter(Binom_Adjp_BH < 0.05) %>% 
      mutate(database = x,
             sample = sample) %>%
      dplyr::select("ID", "name", "Binom_Fold_Enrichment", "Binom_Adjp_BH", "Binom_Region_Set_Coverage", "sample", "database") 
  })
  
  png(paste0("plots/005.1_ATACSeq/005.1_8_Feature_distribution/005.1_8_Feature_distribution", sample,".png"), h = 3000, w = 5000, res = 400)
  plotAnnoPie(macs_anno)
  dev.off()

  macs_gr <- as.GRanges(macs_anno)
  tss_macs_gr <- macs_gr[abs(macs_gr$distanceToTSS) < 500]
  
  macs_gr <- readRDS(paste0("data/Filtered_peaks/", sample, "filtered_peaks.rds"))
  export(macs_gr, paste0("data/Filtered_peaks/", sample, "filtered_peaks.bam"), format = "bam")
  macs_stat[[i]] <- macs_anno@annoStat %>% 
    as_tibble() %>% 
    mutate("Sample" = sample)
}

saveRDS(great_list, "data/005.1_GREAT_list.rds")
saveRDS(macs_stat, "data/Filtered_peaks/Annotation_stats.rds")

sample_labels <- paste(rep(c("control", "siAgc1"), each = 4), 1:4)
sample_levels <- c("C3_S14", "sC1_S1", "sC2_S2", "sC3_S5", "S2_S13",  "sS1_S3", "sS2_S4", "sS3_S6")
macs_stat_plot <- map_dfr(macs_stat, `[`) %>% 
  mutate(Sample_labels = factor(Sample, levels = sample_levels, labels = sample_labels) |> 
           fct_relevel(rev))

macs_stat_plot %>% 
  ggplot(aes(Sample_labels, Frequency, fill = Feature)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Peaks feature distribution",
       x = "Sample",
       y = "Feature (%)") +
  scale_fill_viridis_d(option = "C", end = 0.9)
 
ggsave("plots/005.1_ATACSeq/005.1_8_Feature_distribution/005.1_9_Total_Feature_distribution.png",
       h = 2500,
       w = 5000,
       units = "px",
       dpi = 600)

# Plot GREAT results (to do)
# great_df <- map_dfr(great_list, `[`) %>% 
#   as_tibble() %>% 
#   filter(!sample %in% c("sC3_S5", "sS3_S6"), 
#          database %in% c("GO Biological Process", "MSigDB Pathway", "BioCyc Pathway", "PANTHER Pathway")) %>% 
#   mutate(genotype = ifelse(grepl(pattern = "C", .$sample), "Control", "siAGC1"))
# 
# best_pathway_ctr <- great_df %>% 
#   filter(genotype == "Control") %>% 
#   group_by(name) %>% 
#   summarise(database, score = mean(Binom_Fold_Enrichment*Binom_Region_Set_Coverage*-log10(Binom_Adjp_BH))) %>%
#   slice_head() %>% 
#   ungroup() %>% 
#   group_by(database) %>% 
#   slice_max(order_by = score) %>% 
#   ungroup()
# 
# best_pathway_si <- great_df %>% 
#   filter(genotype == "siAGC1") %>% 
#   group_by(name) %>% 
#   summarise(database, score = mean(Binom_Fold_Enrichment*Binom_Region_Set_Coverage*-log10(Binom_Adjp_BH))) %>%
#   slice_head() %>% 
#   ungroup() %>% 
#   group_by(database) %>% 
#   slice_max(order_by = score) %>% 
#   ungroup()
