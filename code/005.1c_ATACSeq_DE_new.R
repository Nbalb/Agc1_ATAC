library(limma)
library(soGGi)
library(ChIPQC)
library(Rsubread)
library(GenomicAlignments)
library(DESeq2)
library(tracktables)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(clusterProfiler)
library(ChIPseeker)
library(msigdbr)
library(tidyr)
library(stringr)
library(forcats)
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr)
source("code/functions/p_star.R")

# Follows this guide (https://rockefelleruniversity.github.io/RU_ATACseq/presentations/singlepage/RU_ATAC_part2.html)
# Identifying non redundant peaks
peaks <- list.files("data/new_bam", pattern = ".narrowPeak", full.names = T)

myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)
allPeaksSet_nR <- reduce(unlist(GRangesList(myPeaks)))
overlap <- list()
for (i in 1:length(myPeaks)) {
  overlap[[i]] <- allPeaksSet_nR %over% myPeaks[[i]]
}
overlapMatrix <- do.call(cbind, overlap)
colnames(overlapMatrix) <- basename(peaks)
mcols(allPeaksSet_nR) <- overlapMatrix

blklist <- rtracklayer::import.bed("data/blacklist_mm39.bed")
seqlevels(blklist)
nrToCount <- allPeaksSet_nR[!allPeaksSet_nR %over% blklist & seqnames(allPeaksSet_nR) %in% c(1:19, "X", "Y")]
nr_df <- elementMetadata((nrToCount)) %>%
  as.data.frame() %>% 
  dplyr::rename_with(~ gsub("(.*)_(.*)_.*", "\\1_\\2", .x))
colnames(nr_df) <- c("control")

wt_col <- brewer.pal(n = 6, name = "Blues")[3:6]
kd_col <- brewer.pal(n = 6, name = "Reds")[3:6]
dir.create("plots/005.1c", showWarnings = F)
png("plots/005.1c/005.1c_1_Overlap_open_regions_control.png", h = 4000, w = 4000, res = 600)
nr_df %>% 
  dplyr::select(1, 3:5) %>% 
  vennDiagram(main = "Overlap for Control samples' peaks",
              circle.col = wt_col,
              cex = c(1.8, 1.4, 1.2),
              names = c(paste("control", 1:4)),
              lwd = 2)
dev.off()

png("plots/005.1c/005.1c_1_Overlap_open_regions_kd.png", h = 4000, w = 4000, res = 600)
nr_df %>% 
  dplyr::select(2 ,6:8) %>% 
  vennDiagram(main = "Overlap for kd samples' peaks",
              circle.col = kd_col,
              cex = c(1.8, 1.4, 1.2),
              names = c(paste("siAgc1", 1:4)),
              lwd = 2)
dev.off()

nr_occ <- nr_df %>% 
  mutate(occurrence = rowSums(across(where(is.logical))),
         control_sum = rowSums(across(all_of(c("C3_S14", "sC1_S1", "sC2_S2", "sC3_S5")))),
         siAgc1_sum = rowSums(across(all_of(c("S2_S13", "sS1_S3", "sS2_S4", "sS3_S6"))))) %>% 
  rowwise() %>%  
  mutate(control_pres = ifelse(sum(C3_S14, sC1_S1, sC2_S2, sC3_S5) > 0 & occurrence >= 2, T, F),
         siAgc1_pres = ifelse(sum(S2_S13, sS1_S3, sS2_S4, sS3_S6) > 0 & occurrence >= 2, T, F))

png("plots/005.1c/005.1c_1_Overlap_open_regions_kd_vs_wt.png", h = 4000, w = 4000, res = 600)
nr_occ %>% 
  select(control_pres, siAgc1_pres) %>% 
  vennDiagram(main = "Overlap for wt and kd samples",
              circle.col = c("cornflowerblue", "salmon"),
              lwd = 2)
dev.off()

# Counting for differential ATACseq
occurrences <- rowSums(nr_df)

nrToCount <- nrToCount[occurrences >= 2, ]
nr <- nrToCount %>% 
  as_tibble() %>% 
  select(ends_with("Peak")) 
nr <- nr*1 # Turns logical into numeric

# pca <- nr %>% 
#   as.matrix() %>% 
#   t() %>% 
#   prcomp()
# pca_vars <- summary(pca)
# pca$x %>%  
#   as_tibble(rownames = "Sample") %>%
#   mutate(Sample = str_remove(Sample, "_peaks.narrowPeak")) %>% 
#   ggplot(aes(x = PC1, y = PC2, colour = Sample)) +
#   geom_point(size = 5) +
#   labs(title = "PCA ATAC Peaks",
#        x = paste0("PC1 (", pca_vars$importance[2,1]*100, "%)"),
#        y = paste0("PC1 (", pca_vars$importance[2,2]*100, "%)"))
# ggsave("plots/005.1c/005.1c_2_PCA_open_regions.png",
#        h = 1500,
#        w = 2000,
#        units = "px")

bamsToCount <- dir("data/new_bam", full.names = TRUE, pattern = "sorted.bam$")

if(!file.exists("data/005.1c_Peaks_counts_non_redundant.rds")){
  myCounts <- summarizeOverlaps(nrToCount, bamsToCount, singleEnd = FALSE)
  saveRDS(myCounts, "data/005.1c_Peaks_counts_non_redundant.rds")
}else(myCounts <- readRDS("data/005.1c_Peaks_counts_non_redundant.rds"))

colnames(myCounts) <- basename(peaks) %>% stringr::str_extract("[:alnum:]{2,3}_[:alnum:]{2,3}")

# Run DESeq2
Group <- factor(c("control", "siAgc1", rep("control", 3), rep("siAgc1", 3)))
metaData <- data.frame(Group, row.names = colnames(myCounts))

fname <- "data/005.1c_ATAC_deseq_obj.rds"
if(!file.exists(fname)){
  atacDDS <- DESeqDataSetFromMatrix(assay(myCounts), metaData, ~Group, rowRanges = rowRanges(myCounts))
  atacDDS <- DESeq(atacDDS)
  saveRDS(atacDDS, fname)
}else{atacDDS <- readRDS(fname)}

atac_rlog <- rlog(atacDDS)
colData(atac_rlog)$sample <- gsub(".sorted.bam", "", rownames(colData(atac_rlog)))

sample_labels <- paste(rep(c("control", "siAgc1"), each = 4), 1:4)
sample_levels <- c("C3_S14", "sC1_S1", "sC2_S2", "sC3_S5", "S2_S13",  "sS1_S3", "sS2_S4", "sS3_S6")
plot_pca <- atac_rlog
plot_pca@colData$Sample <- rownames(colData(atac_rlog)) %>% 
  factor(levels = sample_levels, labels = sample_labels)
plotPCA(plot_pca, intgroup = "Sample", ntop = nrow(atac_rlog)) + 
  ggtitle("ATACSeq Peaks PCA rlog") +
  geom_point(aes(shape = Group, fill = Sample), size = 5, color = "black") +
  scale_fill_manual(values = c(wt_col, kd_col), name = "Sample") +
  scale_shape_manual(values = c(22, 21)) +
  theme_bw() +
  guides(fill=guide_legend(override.aes=list(alpha=1, size=3, 
                                             color = c(wt_col, kd_col),
                                             shape = c(rep(21, 4), rep(22, 4)))),
         color = "none")
ggsave("plots/005.1c/005.1c_3_ATACSeq_Peaks_PCA_rlog.png", h = 1500, w = 1500, units = "px")

res <- results(atacDDS, c("Group", "kd", "wt"))
res_gr <- results(atacDDS, c("Group", "kd", "wt"), format = "GRanges")
seqlevelsStyle(res_gr) <- "UCSC"
res_gr <- res_gr[seqnames(res_gr) %in% paste0("chr", c(1:19, "X", "Y"))]
res_gr <- res_gr[order(res_gr$pvalue)]

mm39 <- TxDb.Mmusculus.UCSC.mm39.refGene
seqlevelsStyle(mm39) <- "UCSC"
toOverLap <- promoters(mm39, 1000, 500)  # Ask if this range is good
res_gr <- res_gr[(!is.na(res_gr$padj)) & res_gr %over% toOverLap,]
myReport <- makebedtable(res_gr, "data/005.1c_res_gr.html", getwd())
browseURL(myReport)

# Annotate differential ATACSeq regions + pathway analysis. results are not significant, need to check which samples to use for DE ----
anno_res <- annotatePeak(res_gr, TxDb = mm39, verbose = FALSE)

mmdb <- msigdbr(species = "Mus musculus") %>% 
  filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "CP:WIKIPATHWAYS", "GO:BP"))

atac_de <- as.data.frame(anno_res) %>% 
  as_tibble() %>% 
  mutate(geneId = as.character(geneId)) %>% 
  left_join(mmdb %>% 
              dplyr::select(entrez_gene, 
                            gene_symbol, 
                            human_gene_symbol, 
                            num_ortholog_sources, 
                            ensembl_gene) %>% 
              mutate(entrez_gene = as.character(entrez_gene)),
            by = c("geneId" = "entrez_gene"),
            relationship = "many-to-many") %>% 
  relocate(gene_symbol, .before = 1) %>% 
  group_by(gene_symbol) %>%   # There is more than 1 entry per gene_name (transcripts), pick the most recognized one
  arrange(num_ortholog_sources, padj) %>%
  slice_head() %>% 
  ungroup() %>% 
  mutate(gene_symbol = ifelse(is.na(gene_symbol), "Flcn", gene_symbol)) # Checked this transcript on IGV and actually found it is Flcn, (only NA in the column)
write.csv(atac_de %>% arrange(padj), "data/005.1c_ATAC_de.csv")

# Plot most differentialy expressed peaks
de_plot <- atac_de %>%
  mutate(gene_label = paste0(gene_symbol, " - chr", geneChr)) %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  slice_max(n = 10, order_by = log2FoldChange) %>%
  arrange(desc(log2FoldChange)) %>% 
  bind_rows(atac_de %>% 
              mutate(gene_label = paste0(gene_symbol, " - chr", geneChr)) %>% 
              filter(padj < 0.05 & log2FoldChange < 0) %>% 
              slice_min(n = 10, order_by = log2FoldChange) %>%
              arrange(desc(log2FoldChange))) %>% 
  mutate(gene_label = fct_reorder(gene_label, log2FoldChange)) 

de_plot %>% 
  ggplot(aes(gene_label, log2FoldChange, 
             fill = as_factor(sign(log2FoldChange)),
             color = as_factor(sign(log2FoldChange)))) +
  geom_bar(stat='identity', alpha = 0.7) +
  theme_light() +
  theme(legend.position = "none") +
  labs(title = "Best differentially called peaks in aggregated siAgc1 vs ctr", y = "log2FoldChange", x="") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10)) +
  coord_flip() +
  geom_text(aes(label = gene_symbol, 
                y = ifelse(log2FoldChange < 0, 0.05, -0.05),
                hjust = ifelse(log2FoldChange < 0, 0, 1)),
            position = position_dodge(width = 0),
            size = 3,
            lineheight = 0.85,
            color = "black") +
  geom_text(aes(label = p_star(padj), 
                hjust = ifelse(sign(log2FoldChange) < 0, 1.3, -0.3),
                vjust = 0.75),
            color = "black") +
  scale_fill_viridis_d(option = "C", end = 0.7) +
  scale_color_viridis_d(option = "C", end = 0.7)

ggsave("plots/005.1c/005.1c_4_Best_differentially_called_peaks.png",
       h = 1500,
       w = 2500,
       units = "px")

# Overlap chromsome regions with RNASeq
res_chr <- readRDS("data/004_results_chromosome.rds")

for(chr in str_subset(unique(res_chr$mm_chrom), "^chr")){
  
  subchr <- res_chr %>% 
    filter(mm_chrom == chr) %>% 
    filter(!is.na(log2FoldChange) & padj < 0.05) %>% 
    mutate(start = factor(start)) %>% 
    left_join(atac_de %>% filter(seqnames == chr), by = c("symbol" = "gene_symbol"), suffix = c("_rna", "_atac")) %>% 
    mutate(log2FoldChange_atac = ifelse(is.na(log2FoldChange_atac), 0, log2FoldChange_atac)) %>% 
    select(symbol, starts_with("log2"), start_rna) %>% 
    pivot_longer(starts_with("log2"), names_to = "assay", values_to = "log2FoldChange") %>% 
    mutate(assay = str_remove(assay, "log2FoldChange_"))
  
  top_labels <- subchr %>% 
    arrange(desc(abs(log2FoldChange))) %>% 
    slice_head(n = 15) %>% 
    pull(symbol)
  
  if(nrow(subchr) > 10){
    
    png(paste0("plots/005.1c/005.1c_chromsome_comparison/005.1c_", chr, "_comparison.png"), h = 3000, w = 5000, res = 600)
    p <- ggplot(subchr, aes(x = start_rna, y = log2FoldChange, fill = assay)) +
      geom_bar(stat = "identity", position = "identity", alpha = 0.6) +
      labs(title = paste0(chr, "Localized RNA DE and differential ATAC peak calls"),
           subtitle = "padj < 0.05",) + 
      xlab("Chromosomal coordinate") +
      ggrepel::geom_label_repel(aes(label = ifelse(symbol %in% top_labels & assay == "rna", symbol, "")),
                       fontface = "bold",
                       size = 2,
                       # segment.size = 1,
                       # label.size = 1,
                       box.padding = 0.5,
                       max.overlaps = Inf,
                       fill = "#00000020",
                       seed = 41) +
      scale_x_discrete(breaks = br <- levels(subchr$start_rna)[floor(seq(1, nlevels(subchr$start_rna), length.out = 10))],
                       labels = scales::scientific(as.numeric(br))) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
      theme_bw()
    print(p)
    dev.off()
    
  }
}

# Get signatures ----
for(db in c("CP:KEGG", "CP:REACTOME", "CP:WIKIPATHWAYS", "GO:BP")){
  
  fname <- paste0("data/005.1c_GSEA/GSEA_ATAC_", str_replace(db, ":", "_"), ".rds")
  
  if(!file.exists(fname)){
      library(fgsea)
      library(data.table)
      
      mdb_sub <- mmdb %>% 
        filter(gs_cat == db)
      
      mlist <- split(x = mmdb$gene_symbol, f = mmdb$gs_name)
      
      df <- atac_de %>%
        filter(!is.na(padj) & !is.na(log2FoldChange) & !is.na(gene_symbol)) %>% 
        mutate(padj = replace(padj, padj == 0, 2.225074e-308)) 
      
      sig <- setNames(df$stat, df$gene_symbol)
      #Use fgseaMultilevel for better accuracy than fgseaSimple (https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/fgsea/inst/doc/fgseaMultilevel-tutorial.html)
      fgseaRes <- fgseaMultilevel(pathways = mlist,
                                  stats = sig,
                                  eps = 0,
                                  nPermSimple = 10000
      )  
      
      saveRDS(fgseaRes, fname)
      res_csv <- as_tibble(fgseaRes) %>% 
        select(-leadingEdge)
      write.csv(res_csv, paste0("data/005.1c_GSEA/GSEA_ATAC_", db, ".csv"))
      
    }else{fgseaRes <- readRDS(fname)}
    
    # Plots
    if(nrow(fgseaRes[padj < 0.05]) > 0){
      collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                            mlist, sig)
      
      mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
        order(-NES), pathway]
    
      ## Barplot gsea
      source("code/functions/p_star.R")
      source("code/functions/pretty_path_label.R")
      source("code/functions/plotEnrichment2.R")
      
      ends <- fgseaRes %>%
        filter(pathway %in% mainPathways) %>% 
        filter(NES > 0 & padj < 0.05) %>% 
        slice_min(padj, n = 10) %>% 
        arrange(-NES) %>% 
        bind_rows(fgseaRes %>%
                    filter(NES < 0 & padj < 0.05) %>% 
                    slice_min(padj, n = 10) %>% 
                    arrange(-NES)) %>% 
        arrange(NES)
      
      
      png(paste0("plots/006_RNAseq_patients/006_08_agc1_patients_GSEA_barplot.png"), h = 3500, w = 4500, res = 600)
      
      plt <- ggplot(ends, aes(x = c(1:nrow(ends)), y = NES, fill = as.factor(sign(-NES)))) +
        geom_bar(stat='identity') +
        theme_light() +
        theme(legend.position = "none") +
        labs(title = "Best scoring pathways in siAgc1 OliNeu cells", y = "Combined Score", x = "") +
        theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
              axis.title.y = element_text(size = 12),
              axis.text.y = element_text(size = 10),
              axis.text.x = element_text(size = 10)) +
        coord_flip() +
        geom_text(aes(label = pretty_path_label(pathway), 
                      y = ifelse(NES < 0, 0.05, -0.05),
                      hjust = ifelse(NES < 0, 0, 1)),
                  position = position_dodge(width = 0),
                  size = 3,
                  lineheight = 0.85) +
        geom_text(aes(label = p_star(padj)), 
                  hjust = ifelse(sign(ends$NES) < 0, 1.3, -0.3),
                  vjust = 0.75)
      
      
      print(plt)
      dev.off()
      } else {next}
  }

# Volcano plot
EnhancedVolcano(atac_de,
                lab = atac_de$gene_symbol,
                x = "log2FoldChange",
                y = "padj",
                title = "Differential accessibility LFCs in siAgc1 vs control",
                subtitle = NULL,
                selectLab = c("Flcn", "Cops3", "H2ac25", "Tom1l2", "Arf1", "Guk1", "Mrpl55", "Rai1",
                              "Atpaf2", "Cenpv", "Zkscan17", "Zswim7", "Zfp867", "Srebf1", "Trim11", "Mprip", 
                              "Tspan17", "S100a6", "Zfp39"),
                caption = paste0("sig. upregulated = ", atac_de %>% 
                                   filter(log2FoldChange > 0, padj < 0.05) %>% 
                                   nrow(),
                                 "       ",
                                 "sig. downregulated = ", atac_de %>% 
                                   filter(log2FoldChange < 0, padj < 0.05) %>% 
                                   nrow()),
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 1,
                axisLabSize = 14,
                titleLabSize = 14,
                subtitleLabSize = 10,
                captionLabSize = 14,
                labSize = 4.5,
                legendLabSize = 10,
                # legendPosition = "right",
                drawConnectors = T,
                max.overlaps = 20
                ) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("plots/005.1c/005.1c_5_Volcano_plot.png", h = 1500, w = 2500, units = "px")
