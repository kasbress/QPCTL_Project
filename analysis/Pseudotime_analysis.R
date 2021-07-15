library(RColorBrewer)
library(Seurat)
library(metacell)
library(slingshot)
library(tradeSeq)
library(tidyverse)
library(SummarizedExperiment)
library(clusterExperiment)

# Set working directory
setwd("/DATA/users/k.bresser/scRNAseq_QPCTL/Analysis/")

# Import MC object
scdb_init("./metacell_Imm_db2/", force_reinit=T) 
mc <- scdb_mc("QPCTL_Imm_MC")

## 
MCs <- mc@mc[mc@mc %in% c(1,2,3)]

## Load seurat data and add MC identity of relevant MCs
countsfil <- readRDS("./misc_data/seurat.immune.rds")
countsfil <- AddMetaData(object = countsfil, metadata = MCs, col.name = "MCs")

## get MC-2d coordinates
projection  <- scdb_mc2d("QPCTL_Imm_MC")
coo = merge(projection@sc_x, projection@sc_y, by=0) #by=0 merges on rownames - makes a column from rownames
coo <- merge(coo, countsfil@meta.data[,c(13,14)], by.x = "Row.names", by.y = 0, all = T)
coo <- as.matrix(data.frame(row.names = coo[,1], MC_1 = coo[,2], MC_2 = coo[,3]))

# We will now store this as a custom dimensional reduction called 'mcs'
countsfil[["mcs"]] <- CreateDimReducObject(embeddings = coo, key = "MC_", assay = DefaultAssay(countsfil))

# We can now use this as you would any other dimensional reduction in all downstream functions
DimPlot(countsfil, reduction = "mcs", pt.size = 0.5)

## keep only relevant cells by removing NAs for the MC metadata
countsfil.subset <- subset(countsfil, cells = row.names(subset(countsfil@meta.data, !is.na(MCs) ) ))

# Check dimplot, should only include MC1,2,3
DimPlot(countsfil.subset, reduction = "mcs", pt.size = 0.5)

# cleanup
rm(countsfil)
rm(coo)
rm(projection)

# add normalized counts
countsfil.subset <- NormalizeData(object = countsfil.subset, normalization.method = "LogNormalize", scale.factor = 10000)

# create sce and perform slingshot
countsfil.sce <- as.SingleCellExperiment(countsfil.subset)
countsfil.sce <- countsfil.sce[,names(MCs)]
countsfil.sce <- slingshot(data = countsfil.sce, clusterLabels = "MCs", reducedDim = "MCS")



## get vargenes used for MC generation
vargenes <- scdb_gset("QPCTL_Imm_feats_filt")
vargenes <- names(vargenes@gene_set)

# get the log transformed counts from the sce object, subset for vargenes and cellcodes
logexpr <- as.matrix(countsfil.sce@assays@data$logcounts)[vargenes,names(MCs)]
# same, but to get non-transformed data
expr <- as.matrix(countsfil.sce@assays@data$counts)[vargenes,names(MCs)]

# Extract pseudotime values and cellWeights
lineages <- slingPseudotime(countsfil.sce)
cellWeights <- slingCurveWeights(countsfil.sce)

# Evaluate the number of knots needed for GAM fitting
icMat2 <- evaluateK(counts = expr, pseudotime = lineages, cellWeights = cellWeights,
                    k=3:7, nGenes = 100, verbose = FALSE, plot = TRUE)

## fit GAMs using knots = 5
countsfil.sce <- fitGAM(counts = expr, sds = SlingshotDataSet(countsfil.sce), nknots = 5)

# test for dynamic expression, output will have genes associated with pseudotime
assoRes <- associationTest(countsfil.sce)
head(assoRes)

## get significantly associated genes
assoRes %>% 
  rownames_to_column("genes") %>% 
  mutate(padj = p.adjust(pvalue, method = "bonferroni")) %>% 
  filter(padj < 0.05) %>% 
  # slice_max(order_by = waldStat, n = 30) %>% 
  pull(genes) -> genes


###################### Use this to plot clusters

# Subset logexpression matrix to genes associated with pseudotime
hclust_matrix <- logexpr[genes, ]

# scale for clustering
hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

# perform hierchical clustering
gene_dist <- dist(hclust_matrix)
gene_hclust <- hclust(gene_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 30, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

# Decided on 6 clusters, use cuttree to set clusters
gene_cluster <- cutree(gene_hclust, k = 6) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  transmute(gene = name, cluster = value)
head(gene_cluster)

# define function for sd scaling of count data
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# Make dataframe containing pseudotimes per cellcode
pseudotimes <- lineages %>%
  as.data.frame %>% 
  rownames_to_column("cellcode") %>% 
  gather(key = "lineage", value = "pseudotime", -cellcode) %>% 
  drop_na()

# make dataframe containing scaled logcounts of genes per cellcode
logcounts <- (logexpr[genes,]) %>% 
  as.data.frame %>% 
  rownames_to_column("gene") %>%
  gather(key = "cellcode", value = "logcount", -gene) %>% 
  filter(cellcode %in% pseudotimes$cellcode) %>% 
  group_by(gene) %>% 
  mutate(scaled_counts = scale_this(logcount))

# make combined table with pseudotimes, GAM-cluster and gene-expression
to.plot <- pseudotimes %>% 
  left_join(logcounts, by = "cellcode") %>% 
  left_join(gene_cluster, by = "gene")
write_tsv(x = to.plot, "./misc_data/Imm_gene_clusters_pseudotime.tsv")


# Plot the gene(GAM)-clusters
ggplot(to.plot, aes(x = pseudotime, y = scaled_counts))+
  #    geom_point(size = 0.2)+ 
  geom_smooth(aes(group =gene), se = F, color = "lightgrey", method = "gam", size = 0.75)+
  geom_smooth(method = "gam")+
  scale_x_reverse()+
  facet_wrap(~cluster, scales = "free_y")
ggsave(filename = "./metacell_Imm_figs2/Pseudotime_clusters.pdf", device = "pdf", width = 7,height = 4, useDingbats = F )

# To inspect gene-clusters
gene_cluster %>% 
  filter(cluster == 2) %>% 
  pull(gene)




################################### use this for single gene plots
plot_genes <- c("Ly6c2", "Ccr2","Ms4a4c", "Plac8", "Irf8", "Ifitm3")
plot_genes <- c("Eps8", "Creg1", "Pld3","Serpinb6a","Ecm1","Cd5l","Ly6c2", "Ms4a4c","Lyz1", "Lamp1", "Cstb", "Ctsl", "Ccr2", "H2-Ab1", "Cd74", "Mrc1", "Retnla", "Fn1", "Lipa", "Ctsd")
plot_genes <- c("H2-DMa",  "H2-DMb1", "H2-Ab1" , "H2-Aa"  , "H2-Eb1",  "Cd74")

# get logcounts of desired genes
logcounts <- (logexpr[plot_genes,]) %>% 
  as.data.frame %>% 
  rownames_to_column("gene") %>%
  gather(key = "cellcode", value = "logcount", -gene) %>% 
  filter(cellcode %in% pseudotimes$cellcode) 

# merge with pseudotimes
to.plot <- merge(pseudotimes,logcounts, by = "cellcode")

## Make the plot
ggplot(to.plot, aes(x = pseudotime, y = logcount))+
  geom_point(size = 0.2, color = "grey")+ 
  geom_smooth(method = "gam")+
  scale_x_reverse()+
  facet_wrap(~gene, scales = "free_y", ncol = 4)
ggsave(filename = "./metacell_Imm_figs2/Pseudotime_genes.pdf", device = "pdf", width = 7,height = 5, useDingbats = F )
ggsave(filename = "./metacell_Imm_figs2/Pseudotime_genes_points.pdf", device = "pdf", width = 4,height = 7, useDingbats = F )


############################# Plot the rolling average of the sample quantities across pseudotime
library(RcppRoll)

### Plot rolling average cell counts across pseudotime
countsfil.subset@meta.data %>% 
  rownames_to_column("cellcode") %>% 
  dplyr::select(cellcode, hash.ID) %>% 
  # Add count column and normalize to total hashtag size
  mutate(count = 1) %>% 
  group_by(hash.ID)%>%
  mutate(normalized.count = (count/sum(count))*1000 ) %>% 
  ungroup() %>% 
  # collapse genotypes
  mutate(genotype = fct_collapse(hash.ID,WT = c("HTO1", "HTO2", "HTO3"), 
                                 KO = c("HTO4", "HTO5", "HTO6"))) %>% 
  # Add pseudotime data
  left_join(pseudotimes, by = "cellcode") %>% 
  # Make spread WT and KO counts to separate columns
  spread(genotype, normalized.count, fill = 0) %>% 
  # Sort by Pseudotime
  arrange(desc(pseudotime))  %>%
  # Get rolling sums
  mutate(WT_sums = roll_sum(x = WT, n = 60, by = 1, align = "center", fill = NA), 
         KO_sums = roll_sum(x = KO, n = 60, by = 1, align = "center" , fill = NA)) %>% 
  # Prepare for plot
  gather("hash", "count", ends_with("sums")) %>% 
  # Plot data
  ggplot(aes(x = pseudotime, y = count, color = hash))+
  scale_x_reverse()+
  geom_line(size = 1)
ggsave(filename = "./metacell_Imm_figs2/moving_window_Pseudotime.pdf", device = "pdf", width = 6,height = 3, useDingbats = F )
  


########################################################## Analysis looking at terminals instead of continious

# Extract data and plot
countsfil.subset@meta.data %>% 
  rownames_to_column("cellcode") %>% 
  dplyr::select(cellcode, hash.ID) %>% 
  filter(hash.ID != "Doublet") %>% 
  # Add count column and normalize to total hashtag size
  mutate(count = 1) %>% 
  group_by(hash.ID)%>%
  mutate(normalized.count = (count/sum(count))*1000 ) %>% 
  ungroup() %>% 
  # collapse genotypes
  mutate(genotype = fct_collapse(hash.ID,WT = c("HTO1", "HTO2", "HTO3"), 
                                 KO = c("HTO4", "HTO5", "HTO6")))%>% 
  # Add pseudotime data
  left_join(pseudotimes, by = "cellcode") %>% 
  ## subset on the right and left most terminals
  filter(dense_rank(pseudotime) <= 120 | dense_rank(dplyr::desc(pseudotime)) <= 120) %>% 
  mutate(terminal = case_when(pseudotime < 75 ~ "right",
                              pseudotime > 75 ~ "left")) %>% 
  group_by( genotype,terminal, hash.ID) %>% 
  summarise(sum.count = sum(normalized.count)) %>% 
  ## plotting function
  ggplot(aes(x = genotype, y = sum.count))+
  geom_jitter(width = 0.2, size = 2)+
  facet_wrap(~terminal, scales = "free_y")+
  stat_summary(fun.data=mean_se,  geom="pointrange", color="blue", cex = 0.4 )
ggsave("./metacell_Imm_figs2/Pseudotime_terminals_bySide.pdf", device = "pdf", width = 3, height = 4, scale = .75, useDingbats = F)



## Get data
pseudotime.table <- read_tsv("./misc_data/Imm_gene_clusters_pseudotime.tsv")
# Subset on terminals
pseudotime.table %>% 
  dplyr::filter(dense_rank(pseudotime) <= 120 | dense_rank(dplyr::desc(pseudotime)) <= 120) %>% 
  mutate(terminal = case_when(pseudotime < 75 ~ "right",
                              pseudotime > 75 ~ "left")) -> for.plot

# plot expression of genes of interest
for.plot %>% 
  dplyr::filter(gene %in% c("Ly6c2", "Mrc1", "Lyz1", "Lamp1", "Ms4a4c","Eps8", "Ctsd", "H2-Aa")) %>% 
  mutate(polarization = fct_recode(.$terminal, Mono = "left", Macro = "right")) %>% 
  ggplot( aes(x = polarization, y = logcount))+
  geom_violin(aes(fill = polarization))+
  geom_jitter(color = "black", size = 0.1, width = 0.25)+
  facet_wrap(~fct_relevel(gene, c("Ly6c2", "Mrc1", "Lyz1", "Lamp1", "Ms4a4c","Eps8", "Ctsd", "H2-Aa")), nrow = 2)+
  theme_classic()+
  theme(strip.background = element_blank())
ggsave("./metacell_Imm_figs2/Vln_plot_Phenotypes_terminals.pdf", device = "pdf", width = 5, height = 2.5, useDingbats = F)



