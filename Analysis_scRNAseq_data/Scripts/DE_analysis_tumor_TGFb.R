library(metacell)
library(Seurat)
library(limma)
library(tidyverse)
library(fgsea)
library(msigdbr)


# Set working directory
setwd("/DATA/users/k.bresser/scRNAseq_QPCTL/Analysis/")

# Initiate MC database and import MC object
scdb_init("./metacell_tumor_db/", force_reinit=T) 
mc <- scdb_mc("QPCTL_Tum_MC2")

# Get seuratobject containing tumor cells
seurat.Object <- readRDS("./misc_data/seurat.Tum.rds")

# Add MC identities as Metadata
seurat.Object <- AddMetaData(object = seurat.Object, metadata = as.factor(mc@mc), col.name = "MC")

# Collapse groups together for DE comparisons, add as metadata
seurat.Object@meta.data %>% 
  rownames_to_column("cellcode") %>% 
  select(cellcode, hash.ID) %>% 
  deframe %>% 
  as.factor %>% 
  fct_collapse(WT = c("HTO1", "HTO2", "HTO3"), KO = c("HTO4", "HTO5", "HTO6")) -> genotype
seurat.Object <- AddMetaData(object = seurat.Object, metadata = genotype, col.name = "genotype")

# Set idents to genotype and perform DE-analysis 
Idents(seurat.Object) <- "genotype"
marks <- FindMarkers(object = seurat.Object, ident.1 = 'KO', ident.2 = 'WT', logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1, slot = "counts")
write_rds(x = marks, file = "./misc_data/Marks_tumor.rds")


##### Pathway analysis
marks <- read_rds(file = "./misc_data/Marks_tumor.rds")


# Load pathways
pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2"))
pathways <- split(pathways.hallmark[, 5], pathways.hallmark[, 3])

# get stats
marks %>% 
  rownames_to_column("genes") %>% 
  dplyr::filter(p_val < 0.05) %>%  
  dplyr::select(genes, avg_log2FC) %>% 
  deframe %>% 
  sort(decreasing = T) -> stats

## perform analysis
fgseaRes <- fgsea(pathways=pathways, stats=stats, minSize = 10)

BIOCARTA_TGFB_PATHWAY
fgseaRes %>% 
  dplyr::filter(pval < 0.05) %>% 
  dplyr::filter(str_detect(pathway, "TGF"))%>% 
  mutate(pathway = reorder(pathway, NES)) %>% 
  ggplot(aes(x = pathway, y = NES,  fill=NES > 0))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=c("blue", "red"))+
  ggtitle("Pathways adjusted pval < 0.05")+
  theme(legend.position = "none", plot.title = element_text(hjust = 5))+
  coord_flip()+
  gghighlight(padj < 0.05)
ggsave(filename = "./metacell_tumor_figs/Pathways_C2_TGFB.pdf",device = "pdf", width = 6, height = 4,useDingbats=FALSE )


