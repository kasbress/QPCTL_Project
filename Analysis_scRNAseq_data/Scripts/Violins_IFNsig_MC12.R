library(fgsea)
library(msigdbr)
library(tidyverse)
library(metacell)

# Set working directory
setwd("/DATA/users/k.bresser/scRNAseq_QPCTL/Analysis/")

# Import seurat object containing tumor cells
seurat.Object <- readRDS("./misc_data/seurat.Tum.rds")

scdb_init("./metacell_tumor_db/", force_reinit=T) # tell R that this is a "workdirectory"  - RUN THIS EVERYTIME YOU START
mc <- scdb_mc("QPCTL_Tum_MC2")

# Normalize and scale gene expression data
seurat.Object <- NormalizeData(seurat.Object, assay = "RNA", normalization.method = "CLR")
seurat.Object <- ScaleData(seurat.Object, assay = "RNA")

# Load Hallmark pathways
pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "H"))
pathways <- split(pathways.hallmark[, 5], pathways.hallmark[, 3])

# Get scaled data and subset on genes from the IFN signature
seurat.Object %>% 
  GetAssayData( slot = "scale.data", assay = "RNA")%>% 
  as.matrix %>%
  t %>% 
  as.data.frame %>% 
  rownames_to_column("cellcode") %>% 
  as_tibble %>% 
  dplyr::select(one_of(c("cellcode", pathways[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]))) -> full.seurat.data

# Set genes for plotting, selected these from the leading edge from the GSEA performed on MC12
IFNgenes <- c("Ifitm3", "B2m", "Bst2", "Stat1", "Cxcl10", "Psmb9", "Psmb8", "Psmb10", "Psme1", 
              "Isg15", "Psme2", "Tap1", "Rnf213", "Zbp1", "Rsad2", "Ifit3", "Lgals3bp", "Tapbp", "Xaf1", 
              "Gbp3", "Irf1", "Ifi35", "Herc6", "Parp12", "Nmi", "Parp14", "Rtp4", "Lap3", "Eif2ak2", "Cd74", 
              "Samd9l", "Ddx58", "Ogfr", "Psma2")

#Prepare for plotting
full.seurat.data %>% 
  # Switch to long data and join with MC info
  pivot_longer(!cellcode, names_to = "gene", values_to =  "umi") %>% 
  left_join(enframe(mc@mc, "cellcode", "MC")) %>% 
  na.omit %>% 
  mutate(MC = as.factor(MC)) %>% 
  # subset on genes to plot
  filter(gene %in% IFNgenes) %>%
  # Plot as Violins
  ggplot(aes(x = MC, y = umi, fill = MC))+
  geom_violin(scale = "width")+
  facet_wrap(~gene, scales = "free")+
  theme( strip.background = element_blank())
ggsave("./metacell_tumor_figs/Violins_LeadingEdge_MC9_IFN.pdf", width= 12, height = 8)


