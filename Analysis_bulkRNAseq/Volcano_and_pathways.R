library(tidyverse)
library(ggrepel)
library(gghighlight)
library(msigdbr)
library(ggpubr)

# Set working directory
setwd("/DATA/users/k.bresser/QPCTL_bulk/")

DE.results <- read_tsv("./tumor/DE_results_tumor.tsv" )
cpm_table <- read_tsv("./cpm_table_CD45neg.tsv")


##### Volcano plot #######

# Tidy up, add column for highlighting and labeling
DE.results %>% 
  rownames_to_column("genes") %>% 
  mutate(sig =  FDR < 0.05, lab = FDR < 0.01 ) %>% 
  arrange(sig) -> to_plot

# Make the plot
ggplot(to_plot, aes( x = logFC, y = log10(FDR), color = sig, label = genes) ) + 
  geom_point()+
  scale_y_reverse()+
  scale_color_manual(values=c("#c0c0c0", "red"))+
  geom_text_repel(data = subset(to_plot, lab == TRUE ),box.padding = 1, max.overlaps = 15 )+
  labs(title = "")+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  #  xlim(-0.4,0.6)+
  geom_vline(xintercept = 0, linetype ="dotted")+
  geom_hline(yintercept = log10(0.05), linetype ="dotted")
ggsave("./tumor/VolcanoRed_tumor_samples.pdf", device = "pdf", width = 4, height = 6, useDingbats = F)


##### Melanogenisis boxplots #######

# Transpose and select relevant genes
cpm_table %>%
  column_to_rownames("genes") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::select(one_of(c("Tyrp1", "Tyr", "Dct", "Mitf", "Trpm1","Met", "Mlana","Pmel", "Gpnmb","Mc1r"))) %>% 
  # Add group identities
  mutate(group = factor(c("WT","WT","WT","WT","WT","KO","KO" ,"KO","KO","KO","KO"),levels = c("WT", "KO"))) %>%
  # Switch to long data and add median centered cpm
  gather("gene", "cpm", -group) %>% 
  group_by(gene) %>% 
  mutate(norm_cpm = cpm / median(cpm)) %>% 
### plot function
ggplot(aes(x = group, y = norm_cpm))+
  geom_boxplot(aes(color = group))+
  geom_jitter(width = 0.05, color = "black", size = 1)+
  facet_wrap(~gene, nrow = 2)
ggsave("./tumor/Boxplot_Melanogenisis.pdf", device = "pdf", width = 4, height = 3, useDingbats = F)

# Same as above, but calculate P values
cpm_table %>%
  column_to_rownames("genes") %>% 
  t() %>% 
  as.data.frame() %>% 
  select(one_of(c("Tyrp1", "Tyr", "Dct", "Mitf", "Trpm1","Met", "Mlana","Pmel", "Gpnmb","Mc1r"))) %>% 
  mutate(group = factor(c("WT","WT","WT","WT","WT","KO","KO" ,"KO","KO","KO","KO"),levels = c("WT", "KO"))) %>% 
  gather("gene", "cpm", -group) %>% 
  group_by(gene) %>% 
  summarise(pval = t.test(cpm[group=="WT"], cpm[group=="KO"])$p.value) %>% 
  mutate(adj_pval = p.adjust(pval, method = "fdr")) %>% 
  write_tsv(file = "./tumor/table_pval_melanogenisis.tsv")



##### Cell-cycle plots #######

pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "H"))
pathways <- split(pathways.hallmark[, 5], pathways.hallmark[, 3])


### cell cycle waterfalls

# Set pathways to plot
paths <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS")
# Make empty plot list
plot_list <- list()

# For-loop making the plots
for(path in paths){
  p <- DE.results %>% 
    # filter for genes in the pathway (and significant ones)
    filter(Gene.Symbol %in% pathways[[path]]) %>% 
    filter(PValue < 0.05) %>% 
    # reorder gene names by logFC
    mutate(Gene.Symbol = reorder(Gene.Symbol, logFC) ) %>% 
  # Make the plot
  ggplot(aes(x = Gene.Symbol, y = logFC, fill=logFC > 0))+
    geom_bar(stat="identity")+ 
    coord_flip()+
    theme_minimal()+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
    scale_fill_manual(values=c("blue", "red"))+
    labs(title = path)
  plot_list[[path]] <- p
}
ggarrange(plotlist = plot_list)
ggsave("./tumor/CellCycle_waterfalls.pdf", width = 12, height = 8)


### cell cycle Sigs Boxplots

# Set pathways to plot
paths <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS")
# Make empty plot list
plot_list <- list()

# For-loop making the plots
for(path in paths){
  p <- cpm_table %>% 
    # filter for genes in the pathway
    filter(genes %in% pathways[[path]]) %>% 
    # switch to long data
    pivot_longer(cols = !genes, names_to = "samples", values_to = "cpm") %>% 
    # Add grouping factor
    mutate(group = factor(case_when(samples %in% colnames(cpm_table)[1:5] ~ "WT",
                                    TRUE ~ "KO"),levels = c("WT", "KO"))) %>% 
    # Get signature sums
    group_by(group, samples) %>% 
    summarise(signature_sum = sum(cpm)) %>% 
  # Make plots  
  ggplot(aes(x = group, y = signature_sum, color = group))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.1, color = "black", size = 2)+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
    labs(title = path)
  plot_list[[path]] <- p
}
ggarrange(plotlist = plot_list)
ggsave("./tumor/CellCycle_Box_Sum.pdf", width = 12, height = 8)


