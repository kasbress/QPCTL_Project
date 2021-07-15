library(metacell)
library(Seurat)
library(viridis)
library(tidyverse)
library(clustree)
library(factoextra)
library(cluster)

setwd("/DATA/users/k.bresser/scRNAseq_QPCTL/Analysis/")

# Initiate metacell database containing all cells
scdb_init("./metacell_full_db/", force_reinit=T) 
# Import MetaCell object
mc <- scdb_mc("QPCTL_exp_MC_all")
table(mc@mc)

# Select all cells from Immune MetaCells
Imm <- names(mc@mc[mc@mc %in% c(73:76, 77:83, 89)])

# Import seurat object containing all data
seurat_data <- readRDS("./misc_data/seurat_object.rds")

# Subset on Immune cells and save object
seurat_data <- subset(seurat_data, cells = Imm)
saveRDS(object = seurat_data, "./misc_data/seurat.immune.rds")

# Import mat object
mat <- scdb_mat("QPCTL_exp_all")

# Import gene-modules gene-set
gene.clusters <- scdb_gset("All_corr_diff_genes")


# Switch to new directory devoted to Immune cell analysis
if(!dir.exists("./metacell_Imm_db2/")) dir.create("./metacell_Imm_db2/") 
scdb_init("./metacell_Imm_db2/", force_reinit=T)
if(!dir.exists("./metacell_Imm_figs2")) dir.create("./metacell_Imm_figs2") 
scfigs_init("./metacell_Imm_figs2") 

# Add mat object to new Database
scdb_add_mat("QPCTL_exp_all", mat)

# assign hashtag sample IDs
seurat_data <- HTODemux(seurat_data, positive.quantile = .9999)

# Check classification with some informative plots
ggplot(seurat_data@meta.data, aes(x = HTO_margin, color = HTO_classification.global))+
  geom_density()
ggplot(seurat_data@meta.data, aes(x = nCount_RNA, color = HTO_classification.global))+
  geom_density()
ggplot(seurat_data@meta.data, aes(x = nFeature_RNA, color = HTO_classification.global))+
  geom_density()

# Use plots to set tresholds, select Singlets based on gene and UMI counts
seurat_data@meta.data %>% 
  rownames_to_column("cellcodes") %>% 
  filter(HTO_classification.global == "Singlet") %>% 
  filter(nCount_RNA < 10000 & nFeature_RNA < 3000) %>% 
  pull(cellcodes) -> to.keep

# Filter mat object to single immune cells
mcell_mat_ignore_cells(new_mat_id = "QPCTL_exp_Imm", mat_id = "QPCTL_exp_all",ig_cells = to.keep, reverse = T )

mat <- scdb_mat("QPCTL_exp_Imm")

# MAKE GSTAT OBJECT
mcell_add_gene_stat(mat_id = "QPCTL_exp_Imm", gstat_id = "QPCTL_Imm_gs", force = T)
gstat <- scdb_gstat("QPCTL_Imm_gs")
dim(gstat)

# generate feats_gset and plot stats
mcell_gset_filter_varmean(gstat_id = "QPCTL_Imm_gs", gset_id = "QPCTL_Imm_feats", T_vm=0.02, force_new=T)
mcell_gset_filter_cov(gstat_id = "QPCTL_Imm_gs", gset_id = "QPCTL_Imm_feats", T_tot=75, T_top3=2)

# Check features genes and plot statistics
feats_gset <- scdb_gset("QPCTL_Imm_feats")
length(names(feats_gset@gene_set))
mcell_plot_gstats(gstat_id = "QPCTL_Imm_gs", "QPCTL_Imm_feats")

# Set cell cycle and ribosomal proteins to lateral
to.lateral <- names(gene.clusters@gene_set[gene.clusters@gene_set %in% c(20, 2,7)])
test <- rep(1, length(to.lateral))
names(test) <- to.lateral
scdb_add_gset("lateral",gset_new_gset(test, "lateral"))

# filter feats gset for lateral genes
lateral_gset_id = "lateral"
lateral_gset = scdb_gset("lateral")
feats_gset = gset_new_restrict_gset(feats_gset, lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
scdb_add_gset(id = "QPCTL_Imm_feats_filt", gset = feats_gset)
feats_gset <- scdb_gset("QPCTL_Imm_feats_filt")


# metacell generation 
mcell_add_cgraph_from_mat_bknn(mat_id = "QPCTL_exp_Imm", gset_id = "QPCTL_Imm_feats_filt", graph_id = "QPCTL_Imm_graph", K=120, dsamp=T)
mcell_coclust_from_graph_resamp(coc_id = "QPCTL_Imm_coc", graph_id = "QPCTL_Imm_graph", min_mc_size= 40, p_resamp=0.75, n_resamp=300)
mcell_mc_from_coclust_balanced(mc_id = "QPCTL_Imm_MC", coc_id = "QPCTL_Imm_coc",mat_id = "QPCTL_exp_Imm" , K=50, min_mc_size=80,alpha = 2 )

# Add colors to metacells and plot 2D projection
mc <- scdb_mc("QPCTL_Imm_MC")
length(names(mc@mc))
length(names(mc@annots))
mc@colors <- viridis(length(mc@annots))
scdb_add_mc("QPCTL_Imm_MC",mc)
mc <- scdb_mc("QPCTL_Imm_MC")
mcell_mc2d_force_knn(mc2d_id = "QPCTL_Imm_MC" ,mc_id =  "QPCTL_Imm_MC", "QPCTL_Imm_graph", ignore_mismatch = T)
tgconfig::set_param("mcell_mc2d_height",500, "metacell")
tgconfig::set_param("mcell_mc2d_width",500, "metacell")
mcell_mc2d_plot("QPCTL_Imm_MC")
