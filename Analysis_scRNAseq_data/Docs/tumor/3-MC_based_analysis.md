MetaCell based analyses - Tumor compartment
================
Kaspar Bresser
20/10/2021

-   [Import data](#import-data)
-   [Plot marker genes](#plot-marker-genes)
-   [GSEA MC12](#gsea-mc12)
-   [Violin plots](#violin-plots)

Used the analysis below to investigate genes of interest in Tumor
MetaCell 12, which was specifically enriched in QPCTL deficient TMEs

## Import data

Get the metacell object

``` r
scdb_init(here("Data", "Metacell_files_tumor"))
mc <- scdb_mc("QPCTL_Tum_MC2")
```

## Plot marker genes

Check the top and bottom 5 genes of MC12.

``` r
lfp <- as.data.frame(log2(mc@mc_fp))

lfp %>% 
  rownames_to_column("genes") %>% 
  dplyr::select("genes", "12") %>% 
  filter(dense_rank(`12`) <= 5 | dense_rank(-`12`) <= 5) -> TopBot.MC12

TopBot.MC12
```

    ##     genes         12
    ## 1     B2m  2.0402133
    ## 2    Gbp7  1.6783662
    ## 3   Bnip3 -0.4924163
    ## 4  Ifitm3  2.1572273
    ## 5    Bst2  1.8541797
    ## 6   Ero1l -0.4452218
    ## 7    Ugp2 -0.4680698
    ## 8     Fos -0.5096215
    ## 9   Iigp1  1.9753261
    ## 10 Malat1 -0.5892173

and plot as waterfall.

``` r
TopBot.MC12 %>% 
  mutate(genes = reorder(genes, `12`)) %>% 
  ggplot(aes(x = genes, y = `12`, fill = `12` > 0))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values=c("blue", "red"))+
    labs(title = "top and bottom MC12", x = "gene", y = "log2 enrichment")+
    theme(legend.position = "none", plot.title = element_text(hjust = -1))+
    coord_flip()
```

<img src="2-MC_based_analysis_files/figure-gfm/top_bottom_plot-1.png" style="display: block; margin: auto;" />

``` r
ggsave(here("Figs", "tumor", "MC12_top_bottom.pdf"), width = 1, height = 2, scale = 2.5)
```

## GSEA MC12

First get the top and bottom 200 genes of MC12. Will use as input for
fgsea analysis.

``` r
lfp %>% 
  rownames_to_column("genes") %>% 
  dplyr::select("genes", "12") %>% 
  deframe() %>% 
  sort(decreasing = T) -> stats

# subset to 200 highest and 200 lowest values
stats <- stats[dense_rank(stats) <= 200 | dense_rank(-stats) <= 200]

str(stats)
```

    ##  Named num [1:400] 2.16 2.04 1.98 1.85 1.68 ...
    ##  - attr(*, "names")= chr [1:400] "Ifitm3" "B2m" "Iigp1" "Bst2" ...

Now import the pathways

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

And perform the analysis.

``` r
fgseaRes <- fgsea(pathways=pathways, stats=stats, minSize = 10)
```

    ## Warning in fgseaMultilevel(...): For some pathways, in reality P-values are less
    ## than 1e-10. You can set the `eps` argument to zero for better estimation.

Now we can plot the enrichment of the pathways foudn by fgsea.

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

<img src="2-MC_based_analysis_files/figure-gfm/plot_enrichment-1.png" style="display: block; margin: auto;" />

``` r
ggsave(filename = here("Figs", "tumor", "MC12_Pathways_Hallmark.pdf"), device = "pdf", width = 5, height = 4, useDingbats=FALSE )

fgseaRes %>% 
  filter(pval < 0.05)  %>% 
  write_tsv(here("Output", "tumor_MC12_GSEA_Hallmark.tsv"))
```

And plot the top scoring pathways.

``` r
pw <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"

fgseaRes %>%  
  dplyr::filter(pathway == pw) %>%  
  pull(NES) %>% round(digits = 2) -> nes

plotEnrichment(pathway = pathways[[pw]], stats)+
  ggtitle(pw)+
  geom_text(aes(label = paste0(c("NES = ", as.character(nes)), collapse = ""), x = Inf, y = Inf  ), 
            vjust = "inward", hjust = "inward")
```

<img src="2-MC_based_analysis_files/figure-gfm/example_1-1.png" style="display: block; margin: auto;" />

``` r
ggsave(filename = here("Figs", "tumor", "MC12_Pathway_IFNy.pdf"), device = "pdf", width = 5, height = 3, useDingbats=FALSE )
```

``` r
pw <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"

fgseaRes %>%  
  dplyr::filter(pathway == pw) %>%  
  pull(NES) %>% round(digits = 2) -> nes

plotEnrichment(pathway = pathways[[pw]], stats)+
  ggtitle(pw)+
  geom_text(aes(label = paste0(c("NES = ", as.character(nes)), collapse = ""), x = Inf, y = Inf  ), 
            vjust = "inward", hjust = "inward")
```

<img src="2-MC_based_analysis_files/figure-gfm/example_2-1.png" style="display: block; margin: auto;" />

``` r
ggsave(filename = here("Figs", "tumor", "MC12_Pathway_IFNa.pdf"), device = "pdf", width = 5, height = 3, useDingbats=FALSE )
```

## Violin plots

Next, we’ll check how various IFN responsive genes are expressed across
the different MCs.

Import the seurat object.

``` r
seurat.object <- read_rds(here("Data", "seurat_Tum.rds"))
```

Normalize and scale the data.

``` r
seurat.object <- NormalizeData(seurat.object, assay = "RNA", normalization.method = "CLR")
seurat.object <- ScaleData(seurat.object, assay = "RNA")
```

Extract the expression data for all genes in the IFNg geneset.

``` r
IFN.genes <- pathways[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]

seurat.object %>% 
  GetAssayData( slot = "scale.data", assay = "RNA")%>% 
  as.matrix() %>%
  t() %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "cellcode") %>% 
  dplyr::select(one_of(c("cellcode", IFN.genes))) -> IFN.gene.expression
```

    ## Warning: Unknown columns: `Cxcl11`, `H2-Bl`, `Ifi30`, `Marchf1`, `Mx2`

``` r
IFN.gene.expression
```

    ## # A tibble: 10,153 × 206
    ##    cellcode     Adar   Apol6 Arid5b  Arl4a   Auts2     B2m   Bank1  Batf2   Bpgm
    ##    <chr>       <dbl>   <dbl>  <dbl>  <dbl>   <dbl>   <dbl>   <dbl>  <dbl>  <dbl>
    ##  1 AAACCCAAG… -0.493 -0.0887 -0.729 -0.527 -0.0516  0.0164 -0.0199 -0.172  2.11 
    ##  2 AAACCCAAG… -0.493 -0.0887  0.716 -0.527 -0.0516 -0.860  -0.0199 -0.172 -0.643
    ##  3 AAACCCAAG… -0.493 -0.0887 -0.729  2.79  -0.0516 -0.234  -0.0199 -0.172 -0.643
    ##  4 AAACCCAAG… -0.493 -0.0887  0.716 -0.527 -0.0516  0.238  -0.0199 -0.172  2.11 
    ##  5 AAACCCAAG… -0.493 -0.0887 -0.729 -0.527 -0.0516 -0.522  -0.0199 -0.172  2.11 
    ##  6 AAACCCACA… -0.493 -0.0887  0.716 -0.527 -0.0516 -0.860  -0.0199 -0.172 -0.643
    ##  7 AAACCCACA… -0.493 -0.0887 -0.729  1.53  -0.0516  1.33   -0.0199 -0.172 -0.643
    ##  8 AAACCCACA…  1.69  -0.0887 -0.729 -0.527 -0.0516 -0.234  -0.0199 -0.172 -0.643
    ##  9 AAACCCAGT… -0.493 -0.0887  0.716  2.79  -0.0516 -0.234  -0.0199 -0.172 -0.643
    ## 10 AAACCCAGT…  1.69  -0.0887 -0.729 -0.527 -0.0516  0.238  -0.0199 -0.172 -0.643
    ## # … with 10,143 more rows, and 196 more variables: Bst2 <dbl>, Btg1 <dbl>,
    ## #   C1rb <dbl>, C1s2 <dbl>, Casp1 <dbl>, Casp3 <dbl>, Casp4 <dbl>, Casp7 <dbl>,
    ## #   Casp8 <dbl>, Ccl2 <dbl>, Ccl5 <dbl>, Ccl7 <dbl>, Cd274 <dbl>, Cd38 <dbl>,
    ## #   Cd40 <dbl>, Cd69 <dbl>, Cd74 <dbl>, Cd86 <dbl>, Cdkn1a <dbl>, Cfb <dbl>,
    ## #   Cfh <dbl>, Ciita <dbl>, Cmklr1 <dbl>, Cmpk2 <dbl>, Cmtr1 <dbl>,
    ## #   Csf2rb2 <dbl>, Cxcl10 <dbl>, Cxcl9 <dbl>, Ddx58 <dbl>, Ddx60 <dbl>,
    ## #   Dhx58 <dbl>, Eif2ak2 <dbl>, Eif4e3 <dbl>, Epsti1 <dbl>, Fas <dbl>, …

I’ve selected a number of genes of interest, picked out of the leading
edge of the GSEA.

``` r
select.genes <- c("Ifitm3", "B2m", "Bst2", "Stat1", "Cxcl10", "Psmb9", "Psmb8", "Psmb10", "Psme1", 
              "Isg15", "Psme2", "Tap1", "Rnf213", "Zbp1", "Rsad2", "Ifit3", "Lgals3bp", "Tapbp", "Xaf1", 
              "Gbp3", "Irf1", "Ifi35", "Herc6", "Parp12", "Nmi", "Parp14", "Rtp4", "Lap3", "Eif2ak2", "Cd74", 
              "Samd9l", "Ddx58", "Ogfr", "Psma2")
```

prep for plotting.

``` r
IFN.gene.expression %>% 
  pivot_longer(!cellcode, names_to = "gene", values_to =  "umi") %>% 
  left_join(enframe(mc@mc, "cellcode", "MC")) %>% 
  na.omit() %>% 
  mutate(MC = as.factor(MC)) %>% 
  filter(gene %in% select.genes) -> violin.data

violin.data
```

    ## # A tibble: 170,442 × 4
    ##    cellcode           gene       umi MC   
    ##    <chr>              <chr>    <dbl> <fct>
    ##  1 AAACCCACAATCAAGA-1 B2m     -0.860 2    
    ##  2 AAACCCACAATCAAGA-1 Bst2    -0.929 2    
    ##  3 AAACCCACAATCAAGA-1 Cd74    -0.297 2    
    ##  4 AAACCCACAATCAAGA-1 Cxcl10  -0.427 2    
    ##  5 AAACCCACAATCAAGA-1 Ddx58    1.73  2    
    ##  6 AAACCCACAATCAAGA-1 Eif2ak2 -0.691 2    
    ##  7 AAACCCACAATCAAGA-1 Gbp3    -0.275 2    
    ##  8 AAACCCACAATCAAGA-1 Herc6   -0.581 2    
    ##  9 AAACCCACAATCAAGA-1 Ifi35   -0.494 2    
    ## 10 AAACCCACAATCAAGA-1 Ifit3   -0.249 2    
    ## # … with 170,432 more rows

``` r
violin.data %>% 
  ggplot(aes(x = MC, y = umi, fill = MC))+
    geom_violin(scale = "width")+
    facet_wrap(~gene, scales = "free")+
    theme( strip.background = element_blank())
```

<img src="2-MC_based_analysis_files/figure-gfm/violins-1.png" style="display: block; margin: auto;" />

``` r
ggsave(here("Figs", "tumor", "MC_IFN_Violins.pdf"), width = 12, height = 6, scale = 1.2)
```
