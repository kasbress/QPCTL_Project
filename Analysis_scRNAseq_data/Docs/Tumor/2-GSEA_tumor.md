Gene set enrichment analysis tumor cells
================
Kaspar Bresser
20/10/2021

-   [Import and normalize](#import-and-normalize)
-   [DE testing](#de-testing)
-   [Pathway analysis](#pathway-analysis)
    -   [Hallmark pathways](#hallmark-pathways)
    -   [TGFb pathways](#tgfb-pathways)
-   [Volcano plots](#volcano-plots)
    -   [Build signatures](#build-signatures)
        -   [TGFb signature](#tgfb-signature)
        -   [IFN signature](#ifn-signature)
        -   [Cell cycle](#cell-cycle)
    -   [Plot volcano’s](#plot-volcanos)


# Import and normalize

Import the Seurat object containing the tumor cells.

``` r
(seurat.object <- read_rds(here("Data", "seurat_Tum.rds")))
```

    ## An object of class Seurat 
    ## 31075 features across 10153 samples within 3 assays 
    ## Active assay: RNA (31053 features, 0 variable features)
    ##  2 other assays present: HTO, ADT

``` r
scdb_init(here("Data", "Metacell_files_tumor"))
mc <- scdb_mc("QPCTL_Tum_MC2")

seurat.object <- subset(seurat.object, cells = names(mc@mc))

seurat.object@meta.data %>% 
  mutate(hash.ID = case_when(hash.ID != "Singlet" ~ HTO_maxID,
                             TRUE ~ hash.ID)) -> seurat.object@meta.data 

seurat.object
```

    ## An object of class Seurat 
    ## 31075 features across 5013 samples within 3 assays 
    ## Active assay: RNA (31053 features, 0 variable features)
    ##  2 other assays present: HTO, ADT

``` r
seurat.object <- NormalizeData(seurat.object, normalization.method = "CLR", assay = "RNA")
```

# DE testing

Add genotype into the metadata table, first three hashtags are WT cells,
the others are KO cells.

``` r
seurat.object@meta.data %>% 
  mutate(genotype = fct_collapse(hash.ID, 
                                 WT = c("HTO1", "HTO2", "HTO3"), 
                                 KO = c("HTO4", "HTO5", "HTO6"))) -> seurat.object@meta.data 
```

Calculate differentially expressed genes between KO and WT.

``` r
Idents(seurat.object) <- "genotype"
marks <- FindMarkers(object = seurat.object, ident.1 = 'KO', ident.2 = 'WT', logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1)

marks <- as_tibble(marks, rownames = "gene")

write_tsv(marks, here("Output", "tumor_KOvsWT_marks.tsv"))

marks
```

    ## # A tibble: 8,752 × 6
    ##    gene     p_val avg_log2FC pct.1 pct.2 p_val_adj
    ##    <chr>    <dbl>      <dbl> <dbl> <dbl>     <dbl>
    ##  1 B2m   1.68e-77      0.266 0.935 0.899  5.22e-73
    ##  2 Psmb9 4.35e-66      0.262 0.246 0.054  1.35e-61
    ##  3 H2-D1 2.27e-61      0.264 0.769 0.653  7.05e-57
    ##  4 Gbp7  3.58e-52      0.258 0.293 0.111  1.11e-47
    ##  5 Psmb8 8.03e-52      0.250 0.206 0.049  2.49e-47
    ##  6 Gbp2  5.77e-51      0.269 0.227 0.064  1.79e-46
    ##  7 Igtp  2.80e-45      0.230 0.2   0.053  8.71e-41
    ##  8 Cyba  3.52e-45     -0.173 0.459 0.673  1.09e-40
    ##  9 Bst2  5.31e-45      0.238 0.736 0.609  1.65e-40
    ## 10 Iigp1 1.26e-43      0.272 0.171 0.038  3.93e-39
    ## # … with 8,742 more rows

# Pathway analysis

## Hallmark pathways

Lets first import the hallmark pathways

``` r
pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "H"))
pathways <- split(pathways.hallmark$gene_symbol, pathways.hallmark$gs_name)

pathways %>% 
  head(10) %>% 
  str()
```

    ## List of 10
    ##  $ HALLMARK_ADIPOGENESIS           : chr [1:200] "Abca1" "Abcb8" "Acaa2" "Acadl" ...
    ##  $ HALLMARK_ALLOGRAFT_REJECTION    : chr [1:207] "Aars" "Abce1" "Abi1" "Ache" ...
    ##  $ HALLMARK_ANDROGEN_RESPONSE      : chr [1:125] "Abcc4" "Abhd2" "Acsl3" "Actn1" ...
    ##  $ HALLMARK_ANGIOGENESIS           : chr [1:36] "Apoh" "App" "Ccnd2" "Col3a1" ...
    ##  $ HALLMARK_APICAL_JUNCTION        : chr [1:199] "Acta1" "Actb" "Actc1" "Actg1" ...
    ##  $ HALLMARK_APICAL_SURFACE         : chr [1:44] "Adam10" "Adipor2" "Afap1l2" "Akap7" ...
    ##  $ HALLMARK_APOPTOSIS              : chr [1:161] "Add1" "Aifm3" "Ank" "Anxa1" ...
    ##  $ HALLMARK_BILE_ACID_METABOLISM   : chr [1:115] "Abca1" "Abca2" "Abca3" "Abca4" ...
    ##  $ HALLMARK_CHOLESTEROL_HOMEOSTASIS: chr [1:75] "Abca2" "Acat3" "Acss2" "Actg1" ...
    ##  $ HALLMARK_COAGULATION            : chr [1:139] "A2m" "Acox2" "Adam9" "Ang" ...

Get the stats to feed to `fgsea`

``` r
marks %>% 
  filter(p_val < 0.05) %>%  
  select(gene, avg_log2FC) %>% 
  deframe() %>% 
  sort(decreasing = T) -> stats
```

``` r
fgseaRes <- fgseaMultilevel(pathways, stats, eps = 0, minSize = 10)

(fgseaRes <- as_tibble(fgseaRes))
```

    ## # A tibble: 45 × 8
    ##    pathway                pval      padj log2err     ES    NES  size leadingEdge
    ##    <chr>                 <dbl>     <dbl>   <dbl>  <dbl>  <dbl> <int> <list>     
    ##  1 HALLMARK_ADIPOGEN… 3.83e- 1   5.07e-1  0.0838 -0.199 -1.06     62 <chr [8]>  
    ##  2 HALLMARK_ALLOGRAF… 2.00e-10   2.25e-9  0.827   0.728  3.27     30 <chr [15]> 
    ##  3 HALLMARK_ANDROGEN… 4.22e- 1   5.27e-1  0.0796 -0.237 -1.03     29 <chr [14]> 
    ##  4 HALLMARK_APICAL_J… 1.01e- 3   5.68e-3  0.455  -0.425 -2.01     38 <chr [19]> 
    ##  5 HALLMARK_APOPTOSIS 8.34e- 2   1.71e-1  0.288   0.294  1.50     46 <chr [14]> 
    ##  6 HALLMARK_BILE_ACI… 6.77e- 1   7.25e-1  0.0552 -0.243 -0.832    14 <chr [5]>  
    ##  7 HALLMARK_CHOLESTE… 4.75e- 2   1.19e-1  0.322  -0.396 -1.59     22 <chr [13]> 
    ##  8 HALLMARK_COAGULAT… 1.35e- 2   4.67e-2  0.381  -0.483 -1.78     17 <chr [6]>  
    ##  9 HALLMARK_COMPLEME… 7.99e- 2   1.71e-1  0.249   0.335  1.50     30 <chr [8]>  
    ## 10 HALLMARK_DNA_REPA… 8.42e- 1   8.42e-1  0.0644  0.148  0.712    38 <chr [27]> 
    ## # … with 35 more rows

``` r
fgseaRes %>% 
  filter(pval < 0.05)  %>% 
  mutate(pathway = gsub("HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway),
         pathway = reorder(pathway, NES)) %>% 
  ggplot( aes(x = NES, y = pathway))+
    geom_point(aes(size = size, color  = padj))+
    grids(axis = "y", linetype = "dashed")+
    geom_vline(xintercept = 0, color = "blue", linetype = "dotted")
```

<img src="1-GSEA_tumor_files/figure-gfm/plot_pathways-1.png" style="display: block; margin: auto;" />

``` r
ggsave(filename = here("Figs", "tumor", "Pathways_Hallmark.pdf"), device = "pdf", width = 6, height = 4, useDingbats=FALSE )

fgseaRes %>% 
  filter(pval < 0.05)  %>% 
  write_tsv(here("Output", "tumor_GSEA_Hallmark.tsv"))
```

## TGFb pathways

Lets first import the immune related pathways (C2) pathways

``` r
pathways.hallmark <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2"))
pathways <- split(pathways.hallmark$gene_symbol, pathways.hallmark$gs_name)
```

Perform GSEA

``` r
fgseaRes <- fgseaMultilevel(pathways, stats, eps = 0, minSize = 10)

(fgseaRes <- as_tibble(fgseaRes))
```

    ## # A tibble: 2,530 × 8
    ##    pathway                  pval    padj log2err     ES    NES  size leadingEdge
    ##    <chr>                   <dbl>   <dbl>   <dbl>  <dbl>  <dbl> <int> <list>     
    ##  1 ABBUD_LIF_SIGNALING_… 3.08e-4 0.00280  0.498   0.745  2.25     10 <chr [5]>  
    ##  2 ABRAMSON_INTERACT_WI… 8.07e-1 0.887    0.0483 -0.237 -0.705    10 <chr [7]>  
    ##  3 ACEVEDO_LIVER_CANCER… 1.48e-1 0.301    0.143  -0.205 -1.25    120 <chr [28]> 
    ##  4 ACEVEDO_LIVER_CANCER… 2.49e-3 0.0156   0.432   0.231  1.69    225 <chr [27]> 
    ##  5 ACEVEDO_LIVER_CANCER… 3.41e-1 0.507    0.114   0.290  1.10     20 <chr [6]>  
    ##  6 ACEVEDO_LIVER_CANCER… 4.56e-1 0.608    0.0738 -0.236 -1.01     30 <chr [7]>  
    ##  7 ACEVEDO_LIVER_CANCER… 1.13e-1 0.253    0.171  -0.442 -1.42     12 <chr [6]>  
    ##  8 ACEVEDO_LIVER_CANCER… 9.18e-1 0.957    0.0417 -0.171 -0.590    15 <chr [5]>  
    ##  9 ACEVEDO_LIVER_TUMOR_… 3.18e-1 0.488    0.0931 -0.231 -1.11     44 <chr [8]>  
    ## 10 ACEVEDO_LIVER_TUMOR_… 1.25e-1 0.269    0.209   0.166  1.20    206 <chr [29]> 
    ## # … with 2,520 more rows

Select the pathways that contain TGF

``` r
fgseaRes %>% 
  filter(pval < 0.05) %>% 
  filter(str_detect(pathway, "TGF"))%>% 
  mutate(pathway = gsub("_", " ", pathway),
         pathway = reorder(pathway, NES)) %>% 
  ggplot( aes(x = NES, y = pathway))+
    geom_point(aes(size = size, color  = padj))+
    grids(axis = "y", linetype = "dashed")+
    geom_vline(xintercept = 0, color = "blue", linetype = "dotted")
```

<img src="1-GSEA_tumor_files/figure-gfm/plot_TGFb-1.png" style="display: block; margin: auto;" />

``` r
ggsave(filename = here("Figs", "tumor", "Pathways_C2_TGFB.pdf"),device = "pdf", width = 4, height = 3, useDingbats=FALSE )

fgseaRes %>% 
  filter(pval < 0.05) %>% 
  filter(str_detect(pathway, "TGF"))%>% 
  write_tsv(here("Output", "tumor_GSEA_C2_TGFb.tsv"))
```

# Volcano plots

## Build signatures

We’d like to superimpose a few gene-signatures on top of the volcano
plot. Start off with building those from the MsigDB signatures.

### TGFb signature

Will collect TGFb responsive genes by combining the Hallmark TGFb
response from the HALLMARK collection with all signatures from the C2
collection that contain the strings `"TGF"` and `"UP"`.

``` r
msigdbr(species = "Mus musculus", category = "H") %>% 
  as_tibble() %>% 
  filter(str_detect(gs_name, "TGF")) %>% 
  pull(gene_symbol) %>% 
  unique -> genes.TGFb

msigdbr(species = "Mus musculus", category = "C2") %>% 
  as_tibble() %>% 
  filter(str_detect(gs_name, "TGF")) %>% 
  filter(str_detect(gs_name, "UP")) %>% 
  pull(gene_symbol) %>% 
  union(genes.TGFb) %>% 
  unique -> genes.TGFb

write_lines(genes.TGFb, here("Output", "TGFb_signature.txt"))
```

### IFN signature

``` r
msigdbr(species = "Mus musculus", category = "C2") %>% 
  as_tibble() %>% 
  dplyr::filter(str_detect(gs_name, "IFN")) %>% 
  dplyr::filter(!str_detect(gs_name, "DN")) %>% 
  pull(gene_symbol) %>% 
  unique -> genes.IFN

msigdbr(species = "Mus musculus", category = "H") %>% 
  as_tibble() %>% 
  dplyr::filter(str_detect(gs_name, "INTE")) %>% 
  pull(gene_symbol) %>% 
  unique %>% 
  union( genes.IFN) %>% 
  union(grep("^H2\\-K|^H2\\-D|Ifi|Gbp", marks$gene, value = T)) -> genes.IFN

write_lines(genes.IFN, here("Output", "IFN_signature.txt"))
```

### Cell cycle

``` r
paths <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_GLYCOLYSIS")

msigdbr(species = "Mus musculus", category = "H") %>% 
  as_tibble() %>% 
  filter(gs_name %in% paths) %>% 
  pull(gene_symbol) %>% 
  unique() -> genes.cellcycle

write_lines(genes.cellcycle, here("Output", "CellCycle_signature.txt"))
```

## Plot volcano’s

``` r
marks %>%
  mutate(sig = ((gene %in% genes.TGFb) & (p_val_adj < 0.05)) ) %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down") )%>% 
  group_by(sig, direction) %>% 
  tally() %>% 
  pivot_wider(names_from = direction, values_from = n) %>% 
  arrange(desc(sig))-> numbers

marks %>% 
  mutate(sig = ((gene %in% genes.TGFb) & (p_val_adj < 0.05)) ) %>% 
  arrange(sig) %>% 
  ggplot( aes( y= log10(p_val) , x=avg_log2FC, color = sig, label = gene) ) + 
    geom_point()+
    scale_y_reverse()+
    scale_color_manual(values=c("#c0c0c0", "red"))+
    labs(title = "DE tumor | TGFb signature")+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
    xlim(-.3,.3)+
    geom_vline(xintercept = 0, linetype ="dotted")+
    geom_hline(yintercept = log10(1.731084e-06), linetype ="dotted")+
    annotate(geom = "table", x = -.28, y = log10(1e-84), label = list(numbers))
```

<img src="1-GSEA_tumor_files/figure-gfm/TGFb_volcano-1.png" style="display: block; margin: auto;" />

``` r
ggsave(here("Figs", "tumor", "Volcano_TGFb.pdf"), width = 3, height = 4, scale = 1.2)
```

``` r
marks %>%
  mutate(sig = ((gene %in% genes.IFN) & (p_val_adj < 0.05)) ) %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down") )%>% 
  group_by(sig, direction) %>% 
  tally() %>% 
  pivot_wider(names_from = direction, values_from = n) %>% 
  arrange(desc(sig))-> numbers

marks %>% 
  mutate(sig = ((gene %in% genes.IFN) & (p_val_adj < 0.05)) ) %>% 
  arrange(sig) %>% 
  ggplot( aes( y= log10(p_val) , x=avg_log2FC, color = sig, label = gene) ) + 
    geom_point()+
    scale_y_reverse()+
    scale_color_manual(values=c("#c0c0c0", "red"))+
    labs(title = "DE tumor | IFN signature")+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
    xlim(-.3,.3)+
    geom_vline(xintercept = 0, linetype ="dotted")+
    geom_hline(yintercept = log10(1.731084e-06), linetype ="dotted")+
    annotate(geom = "table", x = -.28, y = log10(1e-84), label = list(numbers))
```

<img src="1-GSEA_tumor_files/figure-gfm/IFN_volcano-1.png" style="display: block; margin: auto;" />

``` r
ggsave(here("Figs", "tumor", "Volcano_IFN.pdf"), width = 3, height = 4, scale = 1.2)
```

``` r
melan <- c("Tyrp1", "Tyr", "Dct", "Mitf", "Trpm1","Met", "Mlana","Pmel", "Gpnmb","Mc1r")

marks %>% 
  mutate(sig = gene %in% melan ) %>% 
  arrange(sig) %>% 
  ggplot( aes( y= log10(p_val) , x=avg_log2FC, color = sig, label = gene) ) + 
    geom_point()+
    scale_y_reverse()+
    geom_text_repel(data = . %>% filter(sig == TRUE), box.padding = 2, max.overlaps = 15 )+
    scale_color_manual(values=c("#c0c0c0", "red"))+
    labs(title = "DE tumor | Melanogenisis genes")+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
    xlim(-.3,.3)+
    geom_vline(xintercept = 0, linetype ="dotted")+
    geom_hline(yintercept = log10(1.731084e-06), linetype ="dotted")
```

<img src="1-GSEA_tumor_files/figure-gfm/melan_volcano-1.png" style="display: block; margin: auto;" />

``` r
ggsave(here("Figs", "tumor", "Volcano_melan.pdf"), width = 3, height = 4, scale = 1.2)
```

    ## Warning: ggrepel: 2 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
marks %>%
  mutate(sig = ((gene %in% genes.cellcycle) & (p_val_adj < 0.05)) ) %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down") )%>% 
  group_by(sig, direction) %>% 
  tally() %>% 
  pivot_wider(names_from = direction, values_from = n) %>% 
  arrange(desc(sig))-> numbers

marks %>% 
  mutate(sig = ((gene %in% genes.cellcycle) & (p_val_adj < 0.05)) ) %>% 
  arrange(sig) %>% 
  ggplot( aes( y= log10(p_val) , x=avg_log2FC, color = sig, label = gene) ) + 
    geom_point()+
    scale_y_reverse()+
    scale_color_manual(values=c("#c0c0c0", "red"))+
    labs(title = "DE tumor | Cell cycle related genes")+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
    xlim(-.3,.3)+
    geom_vline(xintercept = 0, linetype ="dotted")+
    geom_hline(yintercept = log10(1.731084e-06), linetype ="dotted")+
    annotate(geom = "table", x = -.28, y = log10(1e-84), label = list(numbers))
```

<img src="1-GSEA_tumor_files/figure-gfm/cell_cycle_volcano-1.png" style="display: block; margin: auto;" />

``` r
ggsave(here("Figs", "tumor", "Volcano_CellCycle.pdf"), width = 3, height = 4, scale = 1.2)
```
