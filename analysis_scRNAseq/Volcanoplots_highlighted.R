library(msigdbr)
library(tidyverse)

# set working directory
setwd("/DATA/users/k.bresser/scRNAseq_QPCTL/Analysis/")


# Import tumor DE data
marks <- read_rds("./misc_data/Marks_tumor.rds")


## Build TGFb signature
pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "H"))
pathways.hallmark %>%
  dplyr::filter(str_detect(gs_name, "TGF")) %>% 
  pull(gene_symbol) %>% 
  unique -> genes.TGFb
pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2"))
pathways.hallmark %>%
  dplyr::filter(str_detect(gs_name, "TGF")) %>% 
  dplyr::filter(str_detect(gs_name, "UP")) %>% 
  pull(gene_symbol) %>% 
  union(genes.TGFb) %>% 
  unique -> genes.TGFb

# Save Tgfb signature
write_rds(x = genes.TGFb, file =  "./misc_data/TGFb_signature.rds")


# Build IFN signature
pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2"))
pathways.hallmark %>%
  dplyr::filter(str_detect(gs_name, "IFN")) %>% 
  dplyr::filter(!str_detect(gs_name, "DN")) %>% 
  pull(gene_symbol) %>% 
  unique -> genes.IFN
pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "H"))
pathways.hallmark %>%
  dplyr::filter(str_detect(gs_name, "INTE")) %>% 
  pull(gene_symbol) %>% 
  unique %>% 
  union( genes.IFN) %>% 
  union(grep("^H2\\-K|^H2\\-D|Ifi|Gbp", to_plot$genes, value = T)) -> genes.IFN

# Save IFN signature
write_rds(x = genes.IFN, file =  "./misc_data/IFN_signature.rds")

# Add column to highlight significant genes in TGFb signature
marks %>% 
  rownames_to_column("genes") %>% 
  mutate(sig = ((genes %in% genes.TGFb) & (p_val_adj < 0.05)) ) %>% 
  arrange(sig) -> to_plot
# Plot volcano plot
ggplot(to_plot, aes( y= log10(p_val) , x=avg_log2FC, color = sig, label = genes) ) + 
  geom_point()+
  scale_y_reverse()+
  scale_color_manual(values=c("#c0c0c0", "red"))+
  labs(title = "DE tumor | TGFb signature")+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlim(-.6,1)+
  geom_vline(xintercept = 0, linetype ="dotted")+
  geom_hline(yintercept = log10(1.731084e-06), linetype ="dotted")
ggsave(filename = "./metacell_tumor_figs/Volcano_tumor_TGF.pdf", device = "pdf", width = 3, height = 4, useDingbats = F)

# Check amount of genes 
to_plot %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down") )%>% 
  group_by(sig, direction) %>% 
  tally

# Add column to highlight significant genes in IFN signature
marks %>% 
  rownames_to_column("genes") %>% 
  mutate(sig = ((genes %in% genes.IFN) & (p_val_adj < 0.05)) ) %>% 
  arrange(sig) -> to_plot
# Make volcanoplot
ggplot(to_plot, aes( y= log10(p_val) , x=avg_log2FC, color = sig, label = genes) ) + 
  geom_point()+
  scale_y_reverse()+
  scale_color_manual(values=c("#c0c0c0", "red"))+
  #  geom_text_repel(data = subset(to_plot, sig == TRUE ),box.padding = 1, max.overlaps = 15 )+
  labs(title = "DE tumor | IFN signature")+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  xlim(-0.6,1.0)+
  geom_vline(xintercept = 0, linetype ="dotted")+
  geom_hline(yintercept = log10(1.731084e-06), linetype ="dotted")
ggsave(filename = "./metacell_tumor_figs/Volcano_tumor_IFN.pdf", device = "pdf", width = 3, height = 4, useDingbats = F)

# Check amount of genes 
to_plot %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down") )%>% 
  group_by(sig, direction) %>% 
  tally
