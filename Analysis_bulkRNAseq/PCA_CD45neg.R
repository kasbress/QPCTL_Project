library(edgeR)
library(tidyverse)
library(ggfortify)
library(pheatmap)

setwd("/DATA/users/k.bresser/QPCTL_bulk/")

all.counts <- read_tsv(file = "./all_gene_counts_bulkRNA_QPCTL.tsv", col_types = cols(chromosome_name = "c"))


## Filter data (Note that CD45 negative cells are annotated as 'Tum')
all.counts %>% 
  na.omit() %>% 
  # select tumor samples and geneIDs
  dplyr::select(contains(c("Tum", "external"))) %>% 
  # Trash rows without reads
  mutate(sum_counts = rowSums(dplyr::select(.,1:12)) ) %>% 
  filter(sum_counts > 0 ) %>% 
  # In case of double gene_ids, select the one with the highest read-count
  group_by(external_gene_id) %>%
  top_n(1, sum_counts) %>% 
  # Simplify names
  setNames(gsub(pattern = "Sample_6362_", replacement = "S", x = names(.))) %>% 
  # Tumor sample 4 was removed due to insufficient quality of sample
  dplyr::select(-S40_Tum04) %>% 
  column_to_rownames("external_gene_id")-> counts

# Create grouping factor for DGE object
group <- factor(c("WT","WT","WT","WT","WT","KO","KO" ,"KO","KO","KO","KO"),levels = c("WT", "KO"))

# Make the DGE object
DGE.object <- DGEList(counts=counts[1:11], 
                      group = group,
                      genes = rownames(counts))


# Get counts per Million reads
countsPerMillion <- as.data.frame(cpm(DGE.object))

# We estimate the variance for each row in the  matrix
var_genes <- apply(countsPerMillion, 1, var)
head(var_genes)

# select top 1000 variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:1000]
countsPerMillion <- countsPerMillion[select_var, ]

# perform PCA on remaining genes
countsPerMillion %>% 
  t() %>% 
  as.data.frame() %>% 
  prcomp(center = T, scale = T) -> pca

# Plot the PCA
as.data.frame(pca$x) %>% 
  mutate( genotype = factor(c("WT","WT","WT","WT","WT","KO","KO" ,"KO","KO","KO","KO"),levels = c("WT", "KO"))) %>% 
  ggplot(aes(x = PC1, y = PC2, group = genotype))+
  geom_point( aes(color=genotype), size = 2 )+
  ggtitle("PCA 1000 most variable genes")
ggsave("./tumor/PCA_samples.pdf", device = "pdf", width = 2.5, height = 2, useDingbats = F, scale = 1.5)
