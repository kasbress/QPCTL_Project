library(edgeR)
library(tidyverse)
library(ggfortify)

# Set working directory
setwd("/DATA/users/k.bresser/QPCTL_bulk/")

# Import count table
all.counts <- read_tsv(file = "./all_gene_counts_bulkRNA_QPCTL.tsv", col_types = cols(chromosome_name = "c"))

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

# Create grouping factor for DGE object
group <- factor(c("WT","WT","WT","WT","WT","KO","KO" ,"KO","KO","KO","KO"),levels = c("WT", "KO"))

# Make the DGE object
DGE.object <- DGEList(counts=counts[1:11], 
                      group = group,
                      genes = rownames(counts))


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

## Calculate sample normalization factors
DGE.object <- calcNormFactors(DGE.object)
DGE.object$samples

# Set design for comparison
design <- model.matrix(~ 1 + group)

# Estimate dispersions
DGE.object <- estimateDisp(DGE.object, design, robust=TRUE)

# General plots
plotBCV(DGE.object)
plotMDS(DGE.object)

# Perform the DE test
fit <- glmQLFit(DGE.object, design)
qlf <- glmQLFTest(fit, coef=2)

# output MD plot
plotMD(qlf)

# Get results and calculate FDR
qlf$table %>% 
  mutate(Gene.Symbol = qlf$genes$genes)%>% 
  mutate(FDR = p.adjust(PValue, method = "fdr"))  -> DE.results

# Write out DE results
write_tsv(x = DE.results, file = "./tumor/DE_results_tumor.tsv" )

