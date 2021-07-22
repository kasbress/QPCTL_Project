library(flowCore)
library(flowSpecs)
library(tidyverse)
library(umap)
library(rstatix)
library(ggpubr)

# Set working directory
setwd("/DATA/users/k.bresser/QPCTL_Phenotype/")

# Import Flow cytometry data as flowset
# Direct read.flowSet function to location of fcs files
fset <- read.flowSet( path = "./2019_12 WT vs KO 9466/" )

# Down-sample to 5000 cells from each file
dsFilt <- sampleFilter(size = 5000, filterId="dsFilter")
result <- flowCore::filter(fset, dsFilt)
fset <- flowCore::Subset(fset, result)

# calculate logicle transformation based on 1 of the file
toTrans <- as.character(grep(pattern = "FJ", x = colnames(fset), value = T ))
lgcl <- estimateLogicle(fset@frames$`export_20191216 TME DAY2_Tum_M1_CD45+.fcs`, channels =  toTrans)
# Apply transformation
trans.fset <- transform(fset, lgcl)

# Export data as long dataframe
# idInfo argument extracts shorthand sample name from fcs file name
flow.data <- flowSet2LongDf(flowObj = trans.fset, idInfo = list("Sample" = "export_20191216 TME DAY2_Tum_|_CD45\\+\\.fcs"))

# Will perform this analysis in triplicate for each group
flow.data %>% 
  rownames_to_column("cell_id") %>% 
  dplyr::filter(grepl('5|7|10|4|6|8', Sample)) %>% 
  group_by(Sample) -> flow.data


# Get names of the markers in the experiment
trans.fset@frames$`export_20191216 TME DAY2_Tum_M1_CD45+.fcs`@parameters@data %>% 
  na.omit %>% 
  pull(desc) -> markers

# Select markers, drop other columns, and rename according to proteins
flow.data %>% 
  dplyr::select(starts_with("FJ"), Sample, cell_id ) %>% 
  dplyr::select(-c(5,11)) %>% # drop channels that were not used
  setNames(., c(markers, "Sample", "cell_id")) %>% # rename columns 
  dplyr::select(!one_of("CD45", "7AAD")) -> flow.data # drop irrelevant markers

# Perform UMAP on flow data
reducU <- umap::umap(flow.data[1:11])

# Extract UMAP as dataframe
umap_plot_df <- data.frame(reducU$layout)

# H-clustering, ward.D usually looks nice. 
hc.norm = hclust(dist(umap_plot_df), method = "ward.D")

# Cut the tree at 10 clusters
flow.data$cluster = factor(cutree(hc.norm, 10))


# Add the coordinates to the flow data and plot the clusters
flow.data %>% 
  bind_cols(umap_plot_df) %>% 
  group_by(Sample) %>% 
  slice_sample(n = 2000) %>% 
ggplot( aes( x = X1, y = X2, color = cluster)) +
  geom_point(size = 0.2) + 
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave("./analysis_tum/UMAP_clust.pdf", width = 4.5, height = 4)

# Plot phenotype plots (Violins)
flow.data %>% 
  gather("marker", "FI", -c(Sample, cluster, cell_id)) %>% 
ggplot( aes(x = as.factor(cluster), y = FI, fill = as.factor(cluster)))+
  geom_violin() +
  facet_wrap(~marker, scales = "free_y", nrow = 2)+
  theme_grey()+
  theme( strip.background = element_blank())
ggsave("./analysis_tum/Violin_phenotypes.pdf", width = 160, height = 40, units = "mm", scale = 2.5)



## Calculate statistics for the WT and KO cells in each cluster
flow.data %>% 
  count(cluster, Sample ) %>% # count amount of cells from each sample per cluster
  mutate(genotype = factor(case_when(Sample %in% c( 'M4', 'M6', 'M8') ~ "WT",
                                     TRUE ~ "KO"), levels = c('WT', 'KO') ) ) %>% # assign genotype info
  group_by(cluster) %>% # we want to test per cluster
  t_test(n ~ genotype ) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  add_xy_position(fun = "mean_se", x = "cluster", dodge = 0.9) %>% # x-position has to dodge around cluster
  mutate(p.adj = round(p.adj, 4))  -> stat.test

# Prep data for plotting
flow.data %>%
  count(cluster, Sample ) %>% 
  mutate(genotype = factor(case_when(Sample %in% c( 'M4', 'M6', 'M8') ~ "WT",
                                     TRUE ~ "KO"), levels = c('WT', 'KO') ) ) %>%
  # make barplot
  ggplot(aes(x = cluster, y = n))+ 
  geom_bar(position = "dodge", stat = "summary", fun = "mean", aes(fill = genotype))+
  stat_summary(color = "black", geom = "errorbar", fun.data = "mean_se", 
               size = 0.4,width = 0.4, position = position_dodge(.9), 
               aes(group = genotype) )+
  stat_pvalue_manual(data = stat.test,  label = "p", 
                     tip.length = 0.01,hide.ns = F, label.size = 4 )+
  scale_fill_manual(values=c("#1B76BD", "#A51E23"))+
  theme_classic()
ggsave("./analysis_spleen/barchart_counts.pdf", width = 5, height = 4)


## Write data
flow.data %>% 
  bind_cols(umap_plot_df) %>% 
  write_tsv( file = "./analysis_tum/tum_data.tsv")
