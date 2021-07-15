library(metacell)
library(Seurat)
library(viridis)
library(tidyverse)
library(clustree)
library(factoextra)
library(cluster)

setwd("/DATA/users/k.bresser/scRNAseq_QPCTL/Analysis/")

# Set up MetaCell folders
if(!dir.exists("./metacell_full_db/")) dir.create("./metacell_full_db/") 
scdb_init("./metacell_full_db/", force_reinit=T) 

if(!dir.exists("./metacell_full_figs")) dir.create("./metacell_full_figs") 
scfigs_init("./metacell_full_figs") 

# Import seurat object
seurat_data <- readRDS("./misc_data/seurat_object.rds")

# Include extra cells
seurat_data@meta.data$new_hash.ID <- ifelse(seurat_data@meta.data$hash.ID == "Negative" & seurat_data@meta.data$HTO_margin > 0.25,
                                            as.character(seurat_data@meta.data$HTO_maxID), as.character(seurat_data@meta.data$hash.ID))


table(seurat_data@meta.data$hash.ID)
table(seurat_data@meta.data$new_hash.ID)
HTOHeatmap(seurat_data, assay = "HTO", ncells = 5000)


## Generate metacell mat object
sce = as.SingleCellExperiment(seurat_data)
mat = scm_import_sce_to_mat(sce)
scdb_add_mat(id = "QPCTL_exp",  mat = mat)

## filter small cells
mcell_plot_umis_per_cell("QPCTL_exp", min_umis_cutoff = 450)  
mcell_mat_ignore_small_cells("QPCTL_exp", "QPCTL_exp", 450)

## Clean-up
remove(sce)
remove(seurat_data)


## Get Mito genes 
mito_genes <- grep(pattern = "^mt-", x = mat@genes, value = TRUE)


## Get Mitofractions
uc = Matrix::colSums(mat@mat)
mito_f = Matrix::colSums(mat@mat[mito_genes, ]) / uc

# Check mito fraction and set threshold to exclude cells
plot(density(mito_f))
abline(v = .2, col = "red", lwd = 2)



## CLEAN UP MAT OBJECT

#filter high mito fractions
mcell_mat_ignore_cells("QPCTL_exp_clean", "QPCTL_exp", names(uc)[mito_f >= 0.20] )
clean_mat = scdb_mat("QPCTL_exp_clean")

#remove doublets and negative (hashtags)
mcell_mat_ignore_cells("QPCTL_exp_clean", "QPCTL_exp_clean", 
                       union(clean_mat@ignore_cells, names(uc)[ !(clean_mat@cell_metadata$new_hash.ID %in% c("HTO1","HTO2","HTO3","HTO4","HTO5","HTO6" )) ] ))
clean_mat = scdb_mat("QPCTL_exp_clean")



# Make gstat object
mcell_add_gene_stat(mat_id = "QPCTL_exp_clean", gstat_id = "QPCTL_exp_gs", force = T)
gstat <- scdb_gstat("QPCTL_exp_gs")
dim(gstat)

# generate feats_gset and plot stats
mcell_gset_filter_varmean(gstat_id = "QPCTL_exp_gs", gset_id = "QPCTL_exp_feats", T_vm=0.10, force_new=T)
mcell_gset_filter_cov(gstat_id = "QPCTL_exp_gs", gset_id = "QPCTL_exp_feats", T_tot=80, T_top3=2)

# Check length feature genes and plot statistics
feats_gset <- scdb_gset("QPCTL_exp_feats")
length(names(feats_gset@gene_set))
mcell_plot_gstats(gstat_id = "QPCTL_exp_gs", "QPCTL_exp_feats")


# caclulate gene-gene corrections from feature gene list
gene.anchors <- names(feats_gset@gene_set)
mcell_mat_rpt_cor_anchors(mat_id = "QPCTL_exp_clean", gene_anchors = gene.anchors, cor_thresh = 0.1, gene_anti = c(),
                          tab_fn = "./misc_data/g2g_correlations.txt", sz_cor_thresh = 0.2)


### Read correlation matrix and generate  gene-set to use for gene-clustering
gcor.mat <- read.table("./misc_data/g2g_correlations.txt", header = T)
foc.genes <- apply(gcor.mat[,-1], 1, which.max)
gset <- gset_new_gset(sets = foc.genes, desc = "Diff Expressed Genes")
scdb_add_gset("All_corr_diff_genes", gset)

# Check amount of clusters with elbow method
fviz_nbclust(t(gcor.mat), kmeans, method = "wss", k.max = 50) + theme_minimal() + ggtitle("the Elbow Method")
ggsave(filename = "./metacell_full_figs/k_means_corrs/elbow.pdf", width = 8, height = 4)

# Check amount of clusters with gap statistic
gap_stat <- clusGap(t(gcor.mat), FUN = kmeans, nstart = 20, K.max = 20, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")



### generate mat object containing only correlated genes
mcell_mat_ignore_genes(new_mat_id = "QPCTL_clean_diff", mat_id = "QPCTL_exp_clean",
                       ig_genes = names(foc.genes), reverse = T)
### cluster  genes and identify modules | Set here to 50 clusters
mcell_gset_split_by_dsmat(gset_id = "All_corr_diff_genes" , mat_id = "QPCTL_clean_diff", K = 50)

feats_gset <- scdb_gset("All_corr_diff_genes")
mat = scdb_mat("QPCTL_clean_diff")

# Plot heatmaps of gene-gene correlations per module
mcell_plot_gset_cor_mats(gset_id = "All_corr_diff_genes", scmat_id = "QPCTL_clean_diff")

# Manually inspect and annotate

## 2 = Ribosomal
## 3 = Unknown good cluster
## 4 = Loose cluster, TFs
## 5 = Some HSPs and slicing factors, translation stuff
## 6 = Unknown, loose cluster
## 7 = DNA replication
## 8 = Large loose cluster
## 9 = Bunch of enzymes, Tyr & Mitf are here! Also cell cycle stuff
## 10 = Immune genes !!!!!
## 11 = Looks to be immune associated
## 12 = respiratory chain components
## 13 = Loose cluster
## 14 = huge loose cluster
## 15 = Loose immune cluster? contains QPCT, some ILs, CXCL13
## 16 = loose cluster
## 17 = Tight cluster Collagens !!!!!
## 18 = loose cluster
## 19 = Antigen processing, Tap, Class II, Proteasome !!!!!!
## 20 = Cell Cycle, Ki67, Top2a !!!!!
## 21 = Unknown, contains Fibronectin, S100, Cd47!!!!!
## 22 = Pretty cluster, unknown, contains cd55
## 23 = loose cluster, contains Vim, Lgals
## 24 = aSMA and collagens!!!!
## 25 = loose cluster
## 26 = CD11b and varion chemokines!!!!!
## 27 = loose Immune cluster, Sell, Cd11a
## 28 = IFN induced cluster!!!!
## 29 = chemokines and chemokine receptors!!!!
## 30 = chemokines, Ccr7, Immune associated !!!!!
## 31 = Small cluster, Cxcl5, cxcl1
## 32 = small, immune or fibroblast cluster, itga2
## 33 = Immune cluster!!!!!
## 34 = Pretty cluster, unknown
## 35 = Immune cluster? Scavenger !!!!
## 36 = Pretty cluster, unknown
## 37 = loose cluster, contains Cxcl14 !!!!
## 38 = small immune, Ig-heavy, Ly6c, ifn activated genes
## 39 = small immune cluster, Lyz1, Retlna !!!!!!!
## 40 = Immune, class II associated !!!!
## 41 = Some tRNA loaders, loose
## 42 = Fos, Jun, HSP
## 43 = Pretty cluster, unknown, contains Gimaps
## 44 = Interesting, contains glycolosis, Melanocyte antigen and Mif!!!!!!
## 45 = HSPs and endoplasmatic reticulum, protein folding?
## 46 = loose cluster
## 47 = small tight immune cluster!!!!
## 48 = Small T cell cluster !!!!!!
## 49 = Hemoglobin !!!!!
## 50 = mitochondrium!!!!!!

# write out general annotations
feats_gset <- scdb_gset("All_corr_diff_genes")
enframe(feats_gset@gene_set, "gene", "cluster") %>% 
  mutate(annotation = case_when( cluster %in% c(2) ~ "ribosome", 
                                 cluster %in% c(5) ~ "translation",
                                 cluster %in% c(7, 20) ~ "Cell cycle",
                                 cluster %in% c(10, 19, 26,29,30,33,39,40,47,48) ~ "Immune genes",
                                 cluster %in% c(12) ~ "respiration",
                                 cluster %in% c(17,21,24) ~ "Fibroblast genes",
                                 cluster %in% c(28) ~ "IFN response",
                                 cluster %in% c(44) ~ "melanocyte",
                                 cluster %in% c(49) ~ "Hemoglobin",
                                 cluster %in% c(50) ~ "mitochondrium",
                                 TRUE ~ "ambigious")) %>% 
  write_tsv("./misc_data/gene_modules.tsv")


# set gene to remove for clustering to lateral
to.lateral <- names(feats_gset@gene_set[feats_gset@gene_set %in% c(2, 5, 12,20, 28)])
test <- rep(1, length(to.lateral))
names(test) <- to.lateral
scdb_add_gset("lateral",gset_new_gset(test, "lateral"))

# make set of feature genes to find cell types
genes <- names(foc.genes)
try.out <- names(feats_gset@gene_set[feats_gset@gene_set %in% c(10,17,19,21, 24,26,29,30,33,39,40,44,47,48,49)])
test <- rep(1, length(try.out))
names(test) <- try.out
scdb_add_gset("find_types",gset_new_gset(test, "find_types"))
feats_gset <- scdb_gset("find_types")



# Run MetaCell clustering
mcell_add_cgraph_from_mat_bknn(mat_id = "QPCTL_exp_all", gset_id = "find_types", graph_id = "QPCTL_exp_graph_all", K=150, dsamp=T)
mcell_coclust_from_graph_resamp(coc_id = "QPCTL_exp_coc_all", graph_id = "QPCTL_exp_graph_all", min_mc_size=50, p_resamp=0.75, n_resamp=400)
mcell_mc_from_coclust_balanced(mc_id = "QPCTL_exp_MC_all", coc_id =  "QPCTL_exp_coc_all", mat_id = "QPCTL_exp_all", K=40, min_mc_size=50, alpha=2)

# Import clustering and add colors
mc <- scdb_mc("QPCTL_exp_MC_all")
length(names(mc@mc))
length(names(mc@annots))
mc@colors <- viridis(length(mc@annots))
scdb_add_mc("QPCTL_exp_MC_all",mc)

# Plot MetaCell 2D projection
mcell_mc2d_force_knn(mc2d_id = "QPCTL_exp_MC_all" ,mc_id =  "QPCTL_exp_MC_all", "QPCTL_exp_graph_all", ignore_mismatch = T)
tgconfig::set_param("mcell_mc2d_height",800, "metacell")
tgconfig::set_param("mcell_mc2d_width",800, "metacell")
mcell_mc2d_plot("QPCTL_exp_MC_all")
