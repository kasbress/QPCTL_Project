library(fgsea)
library(msigdbr)
library(tidyverse)
library(metacell)

# set working directory
setwd("/DATA/users/k.bresser/scRNAseq_QPCTL/Analysis/")

# initiate MC database and get MC object
scdb_init("./metacell_tumor_db/", force_reinit=T) 
mc <- scdb_mc("QPCTL_Tum_MC2")

# Get log2 enrichment values
lfp <- as.data.frame(log2(mc@mc_fp))

## Load Hallmark pathways from MSigDB
pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "H"))

# Create list containing pathways
pathways <- split(pathways.hallmark[, 5], pathways.hallmark[, 3])

# Create gene-stats object to perform enrichment analysis on
lfp %>% 
  rownames_to_column("genes") %>% 
  # select MC12 - metacell enriched in KO tumors
  dplyr::select("genes", "12") %>% 
  # store as named list and sort by lfp value
  deframe %>% 
  sort(decreasing = T) -> stats

# subset to 200 highest and 200 lowest values
stats <- stats[dense_rank(stats) <= 200 | dense_rank(-stats) <= 200]

# Run FGSEA
fgseaRes <- fgsea(pathways=pathways, stats=stats, minSize = 10)

# Plot NES
fgseaRes %>% 
  dplyr::filter(padj < 0.05) %>% 
  mutate(pathway = reorder(pathway, NES)) %>% 
  ggplot(aes(x = pathway, y = NES,  fill=NES > 0))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=c("blue", "red"))+
  ggtitle("Pathways adjusted pval < 0.05")+
  theme(legend.position = "none", plot.title = element_text(hjust = -1))+
  coord_flip()
ggsave(filename = "./metacell_tumor_figs/Pathways_hallmark.pdf",device = "pdf", width = 6, height = 4,useDingbats=FALSE )


### Plot specific pathway
pw <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
nes <- fgseaRes %>%  dplyr::filter(pathway == pw) %>%  pull(NES) %>% round(digits = 2)
plotEnrichment(pathway = pathways[[pw]], stats)+
  ggtitle(pw)+
  geom_text(aes(label = paste0(c("NES = ", as.character(nes)), collapse = ""), x = Inf, y = Inf  ), vjust = "inward", hjust = "inward")
ggsave(filename = "./metacell_tumor_figs/Pathways_IFNy.pdf",device = "pdf", width = 5, height = 3,useDingbats=FALSE )

### Plot specific pathway
pw <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
nes <- fgseaRes %>%  dplyr::filter(pathway == pw) %>%  pull(NES) %>% round(digits = 2)
plotEnrichment(pathway = pathways[[pw]], stats)+
  ggtitle(pw)+
  geom_text(aes(label = paste0(c("NES = ", as.character(nes)), collapse = ""), x = Inf, y = Inf  ), vjust = "inward", hjust = "inward")
ggsave(filename = "./metacell_tumor_figs/Pathways_IFNa.pdf",device = "pdf", width = 5, height = 3,useDingbats=FALSE )

