library(edgeR)
library(tidyverse)
library(ggfortify)
library(pheatmap)

# Set working directory
setwd("/DATA/users/k.bresser/QPCTL_bulk/")

# Import count data
all.counts <- read_tsv(file = "./all_gene_counts_bulkRNA_QPCTL.tsv", col_types = cols(chromosome_name = "c"))

## Select samples
all.counts %>% 
  na.omit() %>% 
  # Select columns containing baseline tissue data
  dplyr::select( c(2:25, contains(c( "external"))) ) %>% 
  # trash rows that have no reads
  mutate(sum_counts = rowSums(dplyr::select(.,1:24)) ) %>% 
  filter(sum_counts > 0 ) %>% 
  # In case of double gene_ids, select the one with the highest read-count
  group_by(external_gene_id) %>%
  top_n(1, sum_counts) %>% 
  # Simplify names
  setNames(gsub(pattern = "Sample_6362_", replacement = "S", x = names(.))) %>% 
  column_to_rownames("external_gene_id")-> counts

## Make the DGE object
DGE.object <- DGEList(counts=counts[1:24], 
                     # group = group,
                      genes = rownames(counts))

# Get counts per million
countsPerMillion <- as.data.frame(cpm(DGE.object) )

# We estimate the variance for each row in the  matrix
var_genes <- apply(countsPerMillion, 1, var)
head(var_genes)

# select the 1000 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:1000]
countsPerMillion <- countsPerMillion[select_var, ]

# perform PCA on remaining genes
countsPerMillion %>% 
  t() %>% 
  as.data.frame() %>% 
  prcomp(center = T, scale = T) -> pca

# Plot the PCA 
as.data.frame(pca$x) %>% 
  mutate(organ = as.factor(c( rep("LN", 8), rep("Spleen", 8), rep("BM", 8) )),
         genotype = as.factor(rep( c("WT","WT","WT","WT","KO","KO","KO","KO"), 3 ))) %>% 
  ggplot(aes(x = PC1, y = PC2, group = organ))+
  geom_point( aes(shape=organ, color=genotype), size = 2 )+
  ggtitle("PCA 1000 most variable genes")
ggsave("./baseline/PCA_samples.pdf", device = "pdf", width = 2.5, height = 2, useDingbats = F, scale = 1.5)


####################### Below for heatmap

cpm_table <- cpm(DGE.object, log = T)

# We estimate the variance for each row in the  matrix
var_genes <- apply(cpm_table, 1, var)
head(var_genes)

# Select the 250 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:250]
highly_variable_genes <- cpm_table[select_var,]

## Set annotations for heatmap
annotation <- data.frame( organ = as.factor(c( rep("LN", 8), rep("Spleen", 8), rep("BM", 8) )),
                          genotype = as.factor(rep( c("WT","WT","WT","WT","KO","KO","KO","KO"), 3 )) )
rownames(annotation) <- colnames(highly_variable_genes) 



pdf(file = "./baseline/heatmap_top250.pdf", width = 6, height = 4)
pheatmap(highly_variable_genes, scale = "row",breaks = seq(-1.5,1.5,by=0.03), 
         clustering_distance_cols = "euclidean" , annotation_col = annotation, show_rownames = F,
         border_color = NA, clustering_distance_rows = "manhattan", main = "250 most variable genes")
dev.off()



