library(tidyverse)

# Set working directory
setwd("/DATA/users/k.bresser/QPCTL_bulk/")

# Import DE results
DE.results <- read_tsv("./tumor/DE_results_tumor.tsv" )

# Import TGFb responsive signature defined during scRNAseq analysis
TGFb_sig <- read_rds("/DATA/users/k.bresser/scRNAseq_QPCTL/Analysis/misc_data/TGFb_signature.rds")


### Make TGFb plots
DE.results %>% 
  filter(Gene.Symbol %in% TGFb_sig) %>% 
    filter(PValue < 0.05) %>% 
  mutate(Gene.Symbol = reorder(Gene.Symbol, logFC) ) %>% 
  ggplot(aes(x = Gene.Symbol, y = logFC, fill=logFC > 0))+
  geom_bar(stat="identity")+ 
  coord_flip()+
  theme_classic()+
  xlab("")+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c("blue", "red"))
ggsave("./tumor/waterfall_Tgfb.pdf", device = "pdf", width = 3, height = 6, useDingbats = F)

