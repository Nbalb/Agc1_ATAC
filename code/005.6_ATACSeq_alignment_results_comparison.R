library(ATACseqQC)

# Import bamfile
bamfile_subr <- list.files(recursive = TRUE, pattern = "^Sorted_ATAC.*.bam$", include.dirs = TRUE)
bamfile_bt <- list.files(recursive = TRUE, pattern = "NOMT.sorted.bam$", include.dirs = TRUE)
bamfile_bwa <- list.files(recursive = TRUE, pattern = "pruned.bam$", include.dirs = TRUE)
bamfile_pos <- list.files(recursive = TRUE, pattern = "^Sorted_S.*.bam$", include.dirs = TRUE)

# Metrics

bamfile.labels <- gsub(".pruned.bam", "", basename(bamfile_bwa))
estimateLibComplexity(readsDupFreq(bamfile))

pdf("plots/005_ATACSeq/005.1_ATACSeq_Pipeline comparison.pdf", w = 8, h = 5)
fragSizeDist(bamfile_pos, paste0("Positive control (SRR891269) ATAC-Seq (Subread)"))
fragSizeDist(bamfile_subr[1], paste0(bamfile.labels[1], " ATAC-Seq (Subread)"))
fragSizeDist(bamfile_subr[2], paste0(bamfile.labels[2], " ATAC-Seq (Subread)"))
fragSizeDist(bamfile_bt[1], paste0(bamfile.labels[1], " ATAC-Seq (Bowtie2)"))
fragSizeDist(bamfile_bt[2], paste0(bamfile.labels[2], " ATAC-Seq (Bowtie2)"))
fragSizeDist(bamfile_bwa[1], paste0(bamfile.labels[1], " ATAC-Seq (Bwa-mem2)"))
fragSizeDist(bamfile_bwa[2], paste0(bamfile.labels[2], " ATAC-Seq (Bwa-mem2)"))
dev.off()
