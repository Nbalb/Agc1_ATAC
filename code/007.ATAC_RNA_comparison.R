library(fgsea)
library(msigdbr)
library(ggrepel)
library(ggpointdensity)
library(ggpubr)
library(msigdbr)
library(fgsea)
library(foreach)
library(tidyverse)

atac_de <- read_csv("data/005.1c_ATAC_de.csv")
dge <- readRDS("data/002_res_tbl.rds")
tfl <- read_tsv("data/TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0.tsv.gz") # TFLink data

merged <- atac_de |> 
  select(-1) |> 
  inner_join(dge, 
             by = c("gene_symbol" = "symbol"), suffix = c("_atac", "_rna"))

s_cor <- cor.test(merged$log2FoldChange_atac, merged$log2FoldChange_rna, method = "spearman")
thr_atac <- 0.7
thr_rna <- 0.5

point_plot <- merged |> 
  mutate(labels = ifelse(gene_symbol %in% unique(tfl$Name.TF) & 
                           abs(log2FoldChange_rna) > thr_rna & 
                           abs(log2FoldChange_atac) > thr_atac, gene_symbol, ""))
  # mutate(labels = ifelse((log2FoldChange_rna > 2 & log2FoldChange_atac > 0) |
  #                          (log2FoldChange_rna > 2 & log2FoldChange_atac < 0) |
  #                          (log2FoldChange_rna < -2 & log2FoldChange_atac > 0) |
  #                          (log2FoldChange_rna < -3 & log2FoldChange_atac < -1) |
  #                          gene_symbol %in% c("Hlx", "Lef1"),
  #                        gene_symbol, ""))
point_plot |> 
  na.omit() |> 
  ggplot(aes(log2FoldChange_atac, log2FoldChange_rna, label = labels)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(thr_atac, -thr_atac), color = "grey", linewidth = 1) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(thr_rna, -thr_rna), color = "grey", linewidth = 1) +
  geom_pointdensity() +
  geom_smooth(method = "lm", 
              formula = y ~ x) +
  geom_text_repel(segment.curvature = -0.1,
                   segment.ncp = 3,
                   segment.angle = 20,
                   # nudge_x = 0.5, 
                   max.overlaps = Inf,
                   min.segment.length = 0) +
  labs(title = "ATAC/RNA LFC comparison",
       y = "RNASeq log2FoldChange",
       x = "ATACSeq log2FoldChange",
       subtitle = paste0(nrow(point_plot), " genes compared      ",
       "Spearman r: ", round(s_cor$estimate, 3), "     ",
       "p-value: ", scales::scientific(s_cor$p.value, 3))) + 
  scale_color_viridis_c(option = "C") +
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_bw()
ggsave("plots/007/001_logFC_ATAC_vs_RNA_thr.png", h = 1200, w = 2000, units = "px")

# Barplot log2FC best genes ATAC
merged |> 
  filter(padj_atac < 0.05) |> 
  arrange(-log2FoldChange_atac) |> 
  slice(-(11:(n()-10))) |> 
  select(gene_symbol, matches("log2|padj")) |> 
  pivot_longer(!gene_symbol, 
               names_to = c(".value", "contrast"), 
               names_sep = "_",
               values_to = ".value",
               names_transform = list(contrast = toupper)) |> 
  mutate(color = ifelse(padj < 0.05, "black", NA)) |> 
  ggplot(aes(reorder(gene_symbol, log2FoldChange), 
             log2FoldChange, 
             fill = contrast,
             color = color)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9), linewidth = 1) +
  theme_bw() +
  labs(title = "Largest FC in differential accessibility vs corresponding differential gene expression",
       x = "Gene symbol",
       color = "padj") + 
  scale_color_identity(guide = "legend", labels = c("signif", "ns")) +
  guides(color = guide_legend(override.aes = list(fill = "grey"))) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/007/ATAC_FC_comparison.png", h = 1200, w = 2400, units = "px")

# Barplot log2FC best genes ATAC
merged |> 
  filter(padj_rna < 0.05) |> 
  arrange(-log2FoldChange_rna) |> 
  slice(-(11:(n()-10))) |> 
  select(gene_symbol, matches("log2|padj")) |> 
  pivot_longer(!gene_symbol, 
               names_to = c(".value", "contrast"), 
               names_sep = "_",
               values_to = ".value",
               names_transform = list(contrast = toupper)) |> 
  mutate(color = ifelse(padj < 0.05, "black", NA)) |> 
  ggplot(aes(reorder(gene_symbol, log2FoldChange), 
             log2FoldChange, 
             fill = contrast,
             color = color)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9), linewidth = 1) +
  theme_bw() +
  labs(title = "Largest FC in differential expression vs corresponding differential accessibility",
       x = "Gene symbol",
       color = "padj") + 
  scale_color_identity(guide = "legend", labels = c("signif", "ns")) +
  guides(color = guide_legend(override.aes = list(fill = "grey"))) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("plots/007/RNA_FC_comparison.png", h = 1200, w = 2400, units = "px")

# Density plot 
dens_plot <- point_plot |> 
  select(gene_symbol, starts_with("Log")) |> 
  pivot_longer(-gene_symbol, names_to = "experiment", names_prefix = "log2FoldChange_", values_to = "log2FC") |> 
  mutate(experiment = str_to_upper(experiment))

mw_res <- wilcox.test(dens_plot |> filter(experiment == "ATAC") |> pull(log2FC),
                      dens_plot |> filter(experiment == "RNA") |> pull(log2FC),
                      alternative = "two.sided")

dens_plot |> 
  ggplot(aes(log2FC, color = experiment, fill = experiment)) +
  geom_density(lwd = 1, alpha = 0.3) +
  theme_bw() +
  labs(title = "Density distributions of log2FoldChanges",
       subtitle = paste0("MWU test p-value: ", mw_res$p.value |> round(4))) +
  scale_color_discrete(name = "Experiment", aesthetics = c("color", "fill"))
ggsave("plots/007/002_logFC_density_comparison.png", h = 1200, w = 2000, units = "px")

mw_res_pos <- wilcox.test(dens_plot |> filter(experiment == "ATAC", log2FC > 0) |> pull(log2FC),
                          dens_plot |> filter(experiment == "RNA", log2FC > 0) |> pull(log2FC),
                          alternative = "two.sided")

dens_plot |> 
  filter(log2FC > 0) |> 
  ggplot(aes(log2FC, color = experiment, fill = experiment)) +
  geom_density(lwd = 1, alpha = 0.3) +
  theme_bw() +
  labs(title = "Density distributions of positive log2FoldChanges",
       subtitle = paste0("MWU test p-value: ", mw_res_pos$p.value |> round(4))) +
  scale_color_discrete(name = "Experiment", aesthetics = c("color", "fill"))
ggsave("plots/007/002_logFC_density_comparison_positive.png", h = 1200, w = 2000, units = "px")

mw_res_neg <- wilcox.test(dens_plot |> filter(experiment == "ATAC", log2FC < 0) |> pull(log2FC),
                          dens_plot |> filter(experiment == "RNA", log2FC < 0) |> pull(log2FC),
                          alternative = "two.sided")

dens_plot |> 
  filter(log2FC < 0) |> 
  ggplot(aes(log2FC, color = experiment, fill = experiment)) +
  geom_density(lwd = 1, alpha = 0.3) +
  theme_bw() +
  labs(title = "Density distributions of negative log2FoldChanges",
       subtitle = paste0("MWU test p-value: ", mw_res_neg$p.value |> round(4))) +
  scale_color_discrete(name = "Experiment", aesthetics = c("color", "fill"))
ggsave("plots/007/002_logFC_density_comparison_negative.png", h = 1200, w = 2000, units = "px")

# Enrichment in ATAC for TFTs 
tft_db <- msigdbr(species = "Mus musculus", category = "C3") |> 
  filter(gs_subcat %in% c("TFT:GTRD", "TFT:TFT_Legacy"))

mlist <- split(x = tft_db$gene_symbol, f = tft_db$gs_name)
df <- atac_de |> 
  as.data.frame() |> 
  filter(!is.na(padj) & !is.na(log2FoldChange) & !is.na(gene_symbol)) |> 
  mutate(padj = replace(padj, padj == 0, 2.225074e-308)) 
sig <- setNames(df$log2FoldChange, df$gene_symbol)
fgseaRes <- fgseaMultilevel(pathways = mlist,
                            stats = sig,
                            eps = 0
) # Only significant results is FOXN3_TARGET_GENES

# find targets of interesting TFs
test_targets <- function(TF, save_name = NULL){
  targets <- tfl |> 
    filter(Name.TF == TF) |> 
    pull(Name.Target)
  targets_lfc <- merged |> filter(gene_symbol %in% targets) |> pull(log2FoldChange_rna)
  not_targets_lfc <- merged |> filter(!(gene_symbol %in% targets)) |> pull(log2FoldChange_rna)
  ks <- ks.test(targets_lfc, not_targets_lfc)
  p <- merged |> 
    mutate(Target = ifelse(gene_symbol %in% targets, "Yes", "No")) |> 
    na.omit() |> 
    ggdensity(x = "log2FoldChange_rna", 
              color = "Target", 
              fill = "Target", 
              size = 1.5, alpha = 0.25,
              add = "mean") + 
    coord_cartesian(xlim = quantile(merged$log2FoldChange_rna, probs = seq(0, 1, 0.01), na.rm = T)[c(3,99)]) +
    labs(title = paste("KS test for", TF, "targets"),
         x = "log2FC RNA",
         subtitle = paste("KS p-value:", scales::scientific(ks$p.value))) +
    scale_fill_discrete(breaks = c("Yes", "No"),
                        labels = c(paste0("Yes"),
                                   paste0("No")), 
                        aesthetics = c("color", "fill"),) 
  print(p)
  if(!is.null(save_name)){
    ggsave(paste0("plots/007/", save_name), plot = p, h = 1200, w = 2000, units = "px")
  }
    
}

test_targets("Srebf1")
test_targets("Srebf1", "Srebf1_ks_test.png")
test_targets("Zeb1", "Zeb1_ks_test.png")
# test_targets("Hlx") # ns
test_targets("Lef1", "Lef1_ks_test.png")
test_targets("Ncor1", "Ncor1_ks_test.png")
test_targets("Rai1", "Rai1_ks_test.png")

# Check which genes have strong downregulation in both RNA-Seq and ATAC-Seq using gsea and fisher test
merge_dn <- atac_de |> 
  filter(log2FoldChange < -0.5, padj < 0.05) |> 
  select(gene_symbol, log2FoldChange, padj) |> 
  inner_join(dge |> 
               filter(log2FoldChange < -0.5, padj < 0.05) |> 
               select(symbol, log2FoldChange, padj),
             by = c("gene_symbol" = "symbol"),
             suffix = c("_atac", "_rna"))
saveRDS(merge_dn, "data/007_OliNeu_signature.rds")
mdf <- msigdbr(species = "Mus musculus") |> 
  filter(gs_cat %in% c("C2", "C5", "H"), 
         gs_subcat %in% c("", "CP:KEGG", "CP:WIKIPATHWAYS", "CP:REACTOME"))
mlist <- split(x = mdf$gene_symbol, f = mdf$gs_name)
sig <- pull(merge_dn, log2FoldChange_rna, gene_symbol)
fgsea_res <- fgseaMultilevel(pathways = mlist, 
                             stats = sig)

cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
ft_results <- foreach(i = unique(mdf$gs_name),
                      .packages = "dplyr",
                   .combine = "c",
                   .inorder = T) %dopar% {
             mdf_sub <- mdf |> 
               filter(gs_name == i)
             y_both <- sum(mdf_sub$gene_symbol %in% merge_dn$gene_symbol)
             y_list <- nrow(merge_dn)-y_both
             y_path <- nrow(mdf_sub)-y_both
             n_both <- n_distinct(mdf$gene_symbol)-y_both-y_list-y_path
             
             mat <- matrix(c(y_both, y_list, y_path, n_both),
                           nrow = 2,
                           dimnames = list(in_list = c("Yes", "No"),
                                           in_pathway = c("Yes", "No")))
             ft_res <- fisher.test(mat, alternative = "two.sided")
             ft_res$p.value
                   }
parallel::stopCluster(cl)
p_adj <- p.adjust(ft_results, method = "BH")
