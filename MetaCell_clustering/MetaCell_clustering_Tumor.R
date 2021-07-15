library(metacell)
library(Seurat)
library(viridis)
library(tidyverse)
library(clustree)
library(factoextra)
library(cluster)

setwd("/DATA/users/k.bresser/scRNAseq_QPCTL/Analysis/")

# Initiate metacell database containing all cells
scdb_init("./metacell_full_db/", force_reinit=T) # tell R that this is a "workdirectory"  - RUN THIS EVERYTIME YOU START

# Import MetaCell object
mc <- scdb_mc("QPCTL_exp_MC_all")
table(mc@mc)

# Select all cells from Tumor MetaCells
Tum <- names(mc@mc[mc@mc %in% c(1:72, 90)])

# Import seurat object containing all data
seurat_data <- readRDS("./misc_data/seurat_object.rds")

# Subset on Tumor cells and save object
seurat_data <- subset(seurat_data, cells = Tum)
saveRDS(object = seurat_data, "./misc_data/seurat.Tum.rds")

# Import mat object
mat <- scdb_mat("QPCTL_exp_all")

# Import gene-modules gene-set
gene.clusters <- scdb_gset("All_corr_diff_genes")

# Switch to new directory devoted to Tumor cell analysis
if(!dir.exists("./metacell_tumor_db/")) dir.create("./metacell_tumor_db/")
scdb_init("./metacell_tumor_db/", force_reinit=T)

if(!dir.exists("./metacell_tumor_figs")) dir.create("./metacell_tumor_figs") 
scfigs_init("./metacell_tumor_figs") 

# Add mat object to new Database
scdb_add_mat("QPCTL_exp_all", mat)

# assign hashtag sample IDs
seurat_data <- HTODemux(seurat_data, positive.quantile = .95)
table(seurat_data@meta.data$hash.ID)

# Check classification with some informative plots
table(seurat_data@meta.data$HTO_classification.global)
ggplot(seurat_data@meta.data, aes(x = log10(nCount_HTO), color = HTO_classification.global))+
  geom_density()
ggplot(seurat_data@meta.data, aes(x = HTO_margin, color = HTO_classification.global))+
  geom_density()
ggplot(seurat_data@meta.data, aes(x = nCount_RNA, color = HTO_classification.global))+
  geom_density()
ggplot(seurat_data@meta.data, aes(x = nFeature_RNA, color = HTO_classification.global))+
  geom_density()


# Make doublet-to-singlet edit in the seurat metadata and re-save object
seurat_data@meta.data %>% 
  mutate(HTO_classification.global = ifelse(HTO_margin > 1, "Singlet",as.character(HTO_classification.global))) %>% 
  mutate(hash.ID = ifelse(HTO_margin > 1, as.character(HTO_maxID) ,as.character(hash.ID))) -> seurat_data@meta.data
saveRDS(object = seurat_data, "./misc_data/seurat.Tum.rds")

# Use plots to set tresholds, select Singlets based on gene and UMI counts
seurat_data@meta.data %>% 
  rownames_to_column("cellcodes") %>% 
  dplyr::filter(HTO_classification.global == "Singlet") %>% 
  dplyr::filter(nCount_RNA < 30000 & nFeature_RNA < 5800) %>% 
  pull(cellcodes) -> to.keep

# Select single fibroblasts and save as new mat object
mcell_mat_ignore_cells(new_mat_id = "QPCTL_exp_Tum", mat_id = "QPCTL_exp_all",ig_cells = to.keep, reverse = T )



### get true
### MAKE GSTAT OBJECT
mcell_add_gene_stat(mat_id = "QPCTL_exp_Tum", gstat_id = "QPCTL_Tum_gs", force = T)
gstat <- scdb_gstat("QPCTL_Tum_gs")
dim(gstat)

# generate feats_gset and plot stats
mcell_gset_filter_varmean(gstat_id = "QPCTL_Tum_gs", gset_id = "QPCTL_Tum_feats", T_vm=0.1, force_new=T)
mcell_gset_filter_cov(gstat_id = "QPCTL_Tum_gs", gset_id = "QPCTL_Tum_feats", T_tot=75, T_top3=2)

# Check features genes and plot statistics
feats_gset <- scdb_gset("QPCTL_Tum_feats")
length(names(feats_gset@gene_set))
mcell_plot_gstats(gstat_id = "QPCTL_Tum_gs", "QPCTL_Tum_feats")

# Set cell cycle and ribosomal proteins to lateral
to.lateral <- names(gene.clusters@gene_set[gene.clusters@gene_set %in% c(20, 2,7)])
test <- rep(1, length(to.lateral))
names(test) <- to.lateral
scdb_add_gset("lateral",gset_new_gset(test, "lateral"))

# filter feats gset for lateral genes
lateral_gset_id = "lateral"
lateral_gset = scdb_gset("lateral")
feats_gset = gset_new_restrict_gset(feats_gset, lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
scdb_add_gset(id = "QPCTL_Tum_feats_filt", gset = feats_gset)
feats_gset <- scdb_gset("QPCTL_Tum_feats_filt")


# metacell generation 
mcell_add_cgraph_from_mat_bknn(mat_id = "QPCTL_exp_Tum", gset_id = "QPCTL_Tum_feats_filt", graph_id = "QPCTL_Tum_graph", K=120, dsamp=T)
mcell_coclust_from_graph_resamp(coc_id = "QPCTL_Tum_coc", graph_id = "QPCTL_Tum_graph", min_mc_size= 100, p_resamp=0.75, n_resamp=300)
mcell_mc_from_coclust_balanced(mc_id = "QPCTL_Tum_MC2", coc_id = "QPCTL_Tum_coc",mat_id = "QPCTL_exp_Tum" , K=70, min_mc_size=200,alpha = 2 )

# Add colors to metacells and plot 2D projection
mc <- scdb_mc("QPCTL_Tum_MC2")
length(names(mc@mc))
length(names(mc@annots))
mc@colors <- viridis(length(mc@annots))
scdb_add_mc("QPCTL_Tum_MC2",mc)
mc <- scdb_mc("QPCTL_Tum_MC2")
mcell_mc2d_force_knn(mc2d_id = "QPCTL_Tum_MC2" ,mc_id =  "QPCTL_Tum_MC2", "QPCTL_Tum_graph", ignore_mismatch = T)
tgconfig::set_param("mcell_mc2d_height",500, "metacell")
tgconfig::set_param("mcell_mc2d_width",500, "metacell")
mcell_mc2d_plot("QPCTL_Tum_MC2")
