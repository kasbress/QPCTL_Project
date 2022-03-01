library(edgeR)
library(tidyverse)
library(pheatmap)

# Set working directory
setwd("/DATA/users/k.bresser/QPCTL_bulk/")

# Import count table
all.counts <- read_tsv(file = "./genecounts.txt", col_types = cols(chromosome_name = "c"))

## Filter data
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

# Set grouping vector for DGE object
group <- factor(c("WT","WT","WT","WT","WT","KO","KO" ,"KO","KO","KO","KO"),levels = c("WT", "KO"))



# Make the DGE object
DGE.object <- DGEList(counts=counts[1:11], 
                      group = group,
                      genes = rownames(counts))

## First we will filter the DGE object for noisyness

# get counts per million reads
countsPerMillion <- as.data.frame(cpm(DGE.object))

# make logicle matrix to check cpm threshold at 1.5
thresh <- countsPerMillion > 1.5
head(thresh)

# Summary of how many TRUEs there are in each row
table(rowSums(thresh))

# we would like to keep genes for which at least 7 samples have > 1.5 cpm
keep <- rowSums(thresh) >= 7
summary(keep)

DGE.object <- DGE.object[keep, , keep.lib.sizes=FALSE]

# Calculate sample normalization factors
DGE.object <- calcNormFactors(DGE.object)



##### Now we can start prepping for the heatmap

# Get cpm values
cpm_table <- as.data.frame(cpm(DGE.object, log = T))

# We estimate the variance for each row in the  matrix
var_genes <- apply(cpm_table, 1, var)

# Select top 250 variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:250]
head(select_var)
highly_variable_genes <- cpm_table[select_var,]



## Set column annotations for heatmap
annotation <- data.frame( genotype = c(rep("WT", 5),rep("KO", 6)) )
rownames(annotation) <- colnames(highly_variable_genes) 

# Plot initial heatmap to inspect gene clusters
pheatmap(highly_variable_genes, scale = "row",breaks = seq(-1.5,1.5,by=0.03), 
         clustering_distance_cols = "euclidean" , annotation_col = annotation, show_rownames = F,
         border_color = NA, cutree_rows = 4)


## perform hierachrical clustering
gene_hclust <- highly_variable_genes %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t() %>% 
  dist %>% 
  hclust(method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 5, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

# Get row annotations
gene_cluster <- cutree(gene_hclust, k = 4) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  transmute(gene = name, cluster = as.factor(value)) %>% 
  column_to_rownames("gene")

# Make the final heatmap
pdf(file = "./tumor/heatmap_tumor_samples.pdf", width = 4, height = 6)
pheatmap(highly_variable_genes, scale = "row",breaks = seq(-1.5,1.5,by=0.03), 
         clustering_distance_cols = "euclidean" , annotation_col = annotation, show_rownames = F,
         border_color = NA, cutree_rows = 4, annotation_row = gene_cluster)
dev.off()

# Write out gene-clusters
gene_cluster %>% 
  rownames_to_column("gene") %>% 
  write_tsv("./gene_clusters_CD45neg.tsv")
