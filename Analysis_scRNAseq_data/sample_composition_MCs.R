library(metacell)
library(Seurat)
library(viridis)
library(tidyverse)

setwd("/DATA/users/k.bresser/scRNAseq_QPCTL/Analysis/")

############################################################ TUMOR

# Import data
scdb_init("./metacell_tumor_db/", force_reinit=T) 
mc <- scdb_mc("QPCTL_Tum_MC2")
seurat.Object = readRDS(file = "./misc_data/seurat.Tum.rds")

# Set groups
seurat.Object@meta.data$Hashtags <- fct_recode(seurat.Object@meta.data$hash.ID, WT_1 = "HTO1", WT_2 = "HTO2", WT_3 = "HTO3",
                                               KO_1 = "HTO4", KO_2 = "HTO5", KO_3 = "HTO6")


## Get metacell annotations and sample hastags and merge to dataframe
seurat.Object@meta.data["Hashtags"] %>% 
  rownames_to_column("cellcode") %>% 
  right_join(enframe(mc@mc, "cellcode", "MC")) %>% 
  filter(Hashtags %in% c("WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3")) -> df

# Get normalized counts
norm.counts <- df %>%
  # first filter unused cells remove NAs
  dplyr::filter( !is.na(MC) ) %>%
  # count MCs by HTO
  dplyr::count(MC, Hashtags)%>%
  # perform the normalization within hash.IDs
  group_by(Hashtags)%>%
  mutate(normalized.count = (n/sum(n))*1000 ) %>% 
  mutate(MC = as.factor(MC) ) 

# Plot the barchart
ggplot(norm.counts, aes(fill=Hashtags, y=normalized.count, x=MC)) + 
  geom_bar( stat="identity", position = "fill")+
  scale_fill_manual(values=c("#2F4F4F", "#528B8B", "#B4CDCD", "#CD5555","#A52A2A", "#8B2323"))+
  theme(legend.title = element_blank())
ggsave(filename = "./metacell_tumor_figs/HTO_fractions.pdf",device = "pdf", width = 4, height = 3,useDingbats=FALSE )

############################################################ FIBROBLAST

# Import data
scdb_init("./metacell_fibroblast_db/", force_reinit=T) 
mc <- scdb_mc("QPCTL_fibro_MC")
seurat.Object = readRDS(file = "./misc_data/seurat.Fibr.rds")

# Set groups
seurat.Object@meta.data$Hashtags <- fct_recode(seurat.Object@meta.data$hash.ID, WT_1 = "HTO1", WT_2 = "HTO2", WT_3 = "HTO3",
                                               KO_1 = "HTO4", KO_2 = "HTO5", KO_3 = "HTO6")


## Get metacell annotations and sample hastags and merge to dataframe
seurat.Object@meta.data["Hashtags"] %>% 
  rownames_to_column("cellcode") %>% 
  right_join(enframe(mc@mc, "cellcode", "MC")) %>% 
  filter(Hashtags %in% c("WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3")) -> df

# Get normalized counts
norm.counts <- df %>%
  # first filter unused cells remove NAs
  dplyr::filter( !is.na(MC) ) %>%
  # count MCs by HTO
  dplyr::count(MC, Hashtags)%>%
  # perform the normalization within hash.IDs
  group_by(Hashtags)%>%
  mutate(normalized.count = (n/sum(n))*1000 ) %>% 
  mutate(MC = as.factor(MC) ) 

# Plot the barchart
ggplot(norm.counts, aes(fill=Hashtags, y=normalized.count, x=MC)) + 
  geom_bar( stat="identity", position = "fill")+
  scale_fill_manual(values=c("#2F4F4F", "#528B8B", "#B4CDCD", "#CD5555","#A52A2A", "#8B2323"))+
  theme(legend.title = element_blank())
ggsave(filename = "./metacell_fibroblast_figs/HTO_fractions.pdf",device = "pdf", width = 4, height = 3,useDingbats=FALSE )


############################################################ IMMUNE

# Import data
scdb_init("./metacell_Imm_db2/", force_reinit=T) 
mc <- scdb_mc("QPCTL_Imm_MC")
seurat.Object = readRDS(file = "./misc_data/seurat.immune.rds")


# Set groups
seurat.Object@meta.data$Hashtags <- fct_recode(seurat.Object@meta.data$hash.ID, WT_1 = "HTO1", WT_2 = "HTO2", WT_3 = "HTO3",
                                               KO_1 = "HTO4", KO_2 = "HTO5", KO_3 = "HTO6")


## Get metacell annotations and sample hastags and merge to dataframe
seurat.Object@meta.data["Hashtags"] %>% 
  rownames_to_column("cellcode") %>% 
  right_join(enframe(mc@mc, "cellcode", "MC")) %>% 
  filter(Hashtags %in% c("WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3")) -> df

# Get normalized counts
norm.counts <- df %>%
  # first filter unused cells remove NAs
  dplyr::filter( !is.na(MC) ) %>%
  # count MCs by HTO
  dplyr::count(MC, Hashtags)%>%
  # perform the normalization within hash.IDs
  group_by(Hashtags)%>%
  mutate(normalized.count = (n/sum(n))*1000 ) %>% 
  mutate(MC = as.factor(MC) ) 

# Plot the barchart
ggplot(norm.counts, aes(fill=Hashtags, y=normalized.count, x=MC)) + 
  geom_bar( stat="identity", position = "fill")+
  scale_fill_manual(values=c("#2F4F4F", "#528B8B", "#B4CDCD", "#CD5555","#A52A2A", "#8B2323"))+
  theme(legend.title = element_blank())
ggsave(filename = "./metacell_Imm_figs2/HTO_fractions.pdf",device = "pdf", width = 4, height = 3,useDingbats=FALSE )

