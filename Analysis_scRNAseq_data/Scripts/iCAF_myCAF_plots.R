library(metacell)
library(Seurat)
library(tidyverse)

# Set working directory 
setwd("/DATA/users/k.bresser/scRNAseq_QPCTL/Analysis/")

# Initiate MC database and get MC object
scdb_init("./metacell_fibroblast_db/", force_reinit=T) 
mc <- scdb_mc("QPCTL_fibro_MC")

# Get log2 enrichment values
lfp <- as.data.frame(log2(mc@mc_fp))

# Import iCAF and myCAF signatures
iCAF <- as.character(read.table("./misc_data/iCAF.txt", header = F)$V1)
myCAF <- as.character(read.table("./misc_data/myCAF.txt", header = F)$V1)

# Calculate CAF signature scores
lfp %>% 
  rownames_to_column("genes") %>% 
  # Assign signatures to genes
  mutate(signature = case_when(genes %in% iCAF ~ "iCAF",
                               genes %in% myCAF ~ "myCAF",
                               TRUE ~ "none")) %>% 
  # Switch to long data and summarise by signature
  pivot_longer(cols = -c(genes, signature), names_to = "MC",values_to = "lfp" ) %>% 
  group_by(signature, MC) %>% 
  summarise(sig_score = sum(lfp)) %>% 
  ungroup %>% 
  # switch to long for dotplot
  pivot_wider( names_from = signature, values_from = sig_score) -> CAFsigs
  
# Plot signature scores
ggplot(CAFsigs, aes(x = iCAF, y = myCAF, label = MC))+
  geom_point(size = 7,color = "dodgerblue")+
  geom_point(size = 7, pch=21,color = "black")+
  geom_text()+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave("./metacell_fibroblast_figs/sig_exprCAFs.pdf", device = "pdf", width = 4, height = 4, scale = .75, useDingbats = F)


############################ Plot the ratio's of myCAFs/iCAFS in each sample

# Import Seurat object
full.seurat <- read_rds("./misc_data/seurat.Fibr.rds")

# Combine sample-hashtags and MCs
enframe(full.seurat$hash.ID, "cellcode", "Hashtag") %>% 
  right_join(enframe(as.factor(mc@mc), "cellcode", "MC")) %>% 
  # Count cells per Hashtag-MC combination
  count(Hashtag, MC) %>% 
  # Normalize per hashtag
  group_by(Hashtag)%>%
  mutate(normalized.count = (n/sum(n))*1000 ) %>% 
  ungroup() %>% 
  # filter to the iCAF and myCAF MCs
  filter(MC %in% c("2","4")) %>% 
  # put counts in separate columns and calculate ratio
  pivot_wider(id_cols = -n, names_from = MC, values_from = normalized.count) %>% 
  mutate(ratio = `2`/`4`) %>% 
  # Add genotype identity
  mutate(genotype = fct_collapse(Hashtag,WT = c("HTO1","HTO2","HTO3"), 
                                 KO = c("HTO4", "HTO5", "HTO6")  )) %>% 
  # Plot the data
  ggplot(aes(x = genotype, y = log2(ratio))) +
  geom_jitter(position=position_jitter(0.1))+
  geom_hline( yintercept = 0,  linetype="dotted")+ 
  stat_summary(fun.data=mean_sdl,  geom="pointrange", color="blue", cex = 0.2)+
  labs(y = "log2 myCAF/iCAF ratio")
ggsave("./metacell_fibroblast_figs/sig_exprCAFsRatio.pdf", device = "pdf", width = 3, height = 4, scale = .75, useDingbats = F)

