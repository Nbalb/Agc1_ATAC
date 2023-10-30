# Test with working dataset
library(Rsubread)
library(Rsamtools)
library(tidyverse)

# Download fastq
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891269/SRR891269_1.fastq.gz
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR891/SRR891269/SRR891269_2.fastq.gz

fq <- list.files("data/Positive_test/", pattern = ".fastq.gz", full.names = TRUE)

memory.limit(size = 500000)
genome <- "data/Positive_test/hg19_Genome.fa"
rsub_ind <- gsub("\\.fa$", "", genome)
buildindex(rsub_ind, genome)

outBAM <- paste0("data/Positive_test/", 
                 basename(fq[1]) |> str_extract("^[[:alnum:]]*"),
                 ".bam")
align(rsub_ind, 
      readfile1 = fq[1],
      readfile2 = gsub("_1.", "_2.", fq[1]),
      output_file = outBAM,
      nthreads = 4,
      type = 1,
      unique = TRUE,
      maxFragLength = 2000)

sortedBAM <- file.path(dirname(outBAM), paste0("Sorted_", basename(outBAM)))

sortBam(outBAM,
        gsub("\\.bam", "", sortedBAM))

indexBam(sortedBAM, gsub("\\.bam", "", sortedBAM))

idxstatsBam(sortedBAM) |> ggplot(aes(seqnames, mapped, fill = seqnames)) + 
  geom_bar(stat = "identity") + coord_flip()

# Reads quality
library(GenomicAlignments)

bampar <- scanBam(sortedBAM)
bampar[[1]]$mapq |> 
  as_data_frame() |> 
  ggplot(aes(x = value)) +
  geom_bar() + 
  stat_bin(aes(label=..count..), vjust = 0, geom="text", position="identity")
sum(is.na(bampar[[1]]$mapq)) # 3362019
unique(bampar[[1]]$mapq)  # 40 20 13  8 10  6 NA
bampar[[1]]$mapq |> 
  as_data_frame() |>
  group_by(value) |> 
  summarise(n = n()) |> 
  ungroup() |> 
  mutate(perc = n/sum(n))

# Plot fragments
atacReads <- readGAlignmentPairs(
  sortedBAM, 
  param = ScanBamParam(mapqFilter = 20, 
                       what = c("qname", "mapq", "isize"),
                       flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE)))

chromosomes <- levels(seqnames(atacReads))[levels(seqnames(atacReads)) |> str_starts("NC_.*")]
atacReads_fil <- atacReads[seqnames(atacReads) %in% paste0("chr", c(1:22, "X", "Y"))]
atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)

table(insertSizes) |> 
  data.frame |> 
  rename(InsertSize = insertSizes, Count = Freq) |> 
  mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
         Count = as.numeric(as.vector(Count))) |> 
  ggplot(aes(x = InsertSize, y = Count)) + 
  geom_line() +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0,2000,100)) +
  labs(title = "ATAC-Seq counts distribution from Buenrostro et al(2013)")

ggsave("plots/005_ATACSeq/005.5_Counts distribution positive test.png",
       h = 3000,
       w = 5000,
       units = "px",
       dpi = 600)
