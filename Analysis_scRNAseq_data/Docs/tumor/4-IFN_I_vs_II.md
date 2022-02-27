Compare IFN type I vs type II
================
Kaspar Bresser
11/10/2021

-   [CytoSig analysis (from dataset)](#cytosig-analysis-from-dataset)
    -   [Import and tidy](#import-and-tidy)
    -   [Calculate scores](#calculate-scores)
    -   [Map mouse to human gene
        symbols](#map-mouse-to-human-gene-symbols)
    -   [Compare signatures](#compare-signatures)
-   [CytoSig analysis (using webtool)](#cytosig-analysis-using-webtool)

In the bulk and scRNAseq analysis of QPCTL WT and KO mice we detected a
enhanced IFN response signature in tumors in a completely KO setting. In
the analyses below we’ll investigate if this is likely due to type I or
type II IFN sensing. To figure out whether we can demonstrate this I
made use of the [CytoSig](https://cytosig.ccr.cancer.gov/) database.

``` r
library(here)
library(metacell)
library(tidytext)
library(rstatix)
library(GGally)
library(biomaRt)
library(tidyverse)
```

# CytoSig analysis (from dataset)

## Import and tidy

Import the cytosig dataset. Each column represents a cytokine treatment
experiment on a cell model under certain conditions. Each row represents
a target gene. Each value is a log2 fold change between treatment and
control conditions.

``` r
cytosig.data <- read_tsv(here( "Data", "Cytosig", "diff.merge"))

cytosig.data
```

    ## # A tibble: 19,918 × 2,057
    ##    Gene     `TGFB1@Condition… `TGFB1@Conditio… `IFNA&Duration:… `IFNA&Duration:…
    ##    <chr>                <dbl>            <dbl>            <dbl>            <dbl>
    ##  1 A1BG              -0.0995           0.114             0.0238       -0.0261   
    ##  2 A1BG-AS1          -0.120            0.0128           -0.0326       -0.0000583
    ##  3 A1CF               0.00718          0.00718           0.117        -0.110    
    ##  4 A2M                1.05            -0.242             0.0392       -0.00126  
    ##  5 A4GALT             0.391           -0.632            -0.0101        0.182    
    ##  6 A4GNT              0               -0.0244            0.128        -0.0402   
    ##  7 AAAS              -0.808           -0.132            -0.0533        0.0782   
    ##  8 AACS              -0.488           -0.0374           -0.0364        0.0176   
    ##  9 AADAC              0.191            0.0488            0.100        -0.0996   
    ## 10 AADAT              0.0553           0.330             0.0463        0.0731   
    ## # … with 19,908 more rows, and 2,052 more variables:
    ## #   IFNA&Duration:1h@Condition:Huh7@GSE48400.MicroArray.GPL10558 <dbl>,
    ## #   IFNA&Duration:24h@Condition:Huh7@GSE48400.MicroArray.GPL10558 <dbl>,
    ## #   IFNA&Duration:2h@Condition:Huh7@GSE48400.MicroArray.GPL10558 <dbl>,
    ## #   IFNA&Duration:4h@Condition:Huh7@GSE48400.MicroArray.GPL10558 <dbl>,
    ## #   IFNA&Duration:6h@Condition:Huh7@GSE48400.MicroArray.GPL10558 <dbl>,
    ## #   IFNB&Duration:0.5h@Condition:Huh7@GSE48400.MicroArray.GPL10558 <dbl>, …

Select only to columns that relate to our cytokines of interest, switch
to long form, and extract cytokine identity into a separate column.

``` r
cytosig.data %>% 
  select(Gene, contains("IFNA"), contains("IFNB"), contains("IFNG")) %>% 
  pivot_longer(cols = -Gene, names_to = "assay", values_to = "FC") %>% 
  mutate(cytokine = str_extract(assay, "^IFN\\w{1}")) -> cytosig.IFN

cytosig.IFN
```

    ## # A tibble: 8,126,544 × 4
    ##    Gene  assay                                                       FC cytokine
    ##    <chr> <chr>                                                    <dbl> <chr>   
    ##  1 A1BG  IFNA&Duration:0.5h@Condition:Huh7@GSE48400.MicroArra…  0.0238  IFNA    
    ##  2 A1BG  IFNA&Duration:12h@Condition:Huh7@GSE48400.MicroArray… -0.0261  IFNA    
    ##  3 A1BG  IFNA&Duration:1h@Condition:Huh7@GSE48400.MicroArray.…  0.0232  IFNA    
    ##  4 A1BG  IFNA&Duration:24h@Condition:Huh7@GSE48400.MicroArray… -0.0114  IFNA    
    ##  5 A1BG  IFNA&Duration:2h@Condition:Huh7@GSE48400.MicroArray.… -0.0377  IFNA    
    ##  6 A1BG  IFNA&Duration:4h@Condition:Huh7@GSE48400.MicroArray.… -0.0328  IFNA    
    ##  7 A1BG  IFNA&Duration:6h@Condition:Huh7@GSE48400.MicroArray.… -0.0535  IFNA    
    ##  8 A1BG  IFNA@Condition:SLE keratinocytes@GSE124939.RNASeq.SR… -0.00718 IFNA    
    ##  9 A1BG  IFNA@Condition:normal keratinocytes@GSE124939.RNASeq…  0.00718 IFNA    
    ## 10 A1BG  IFNA&Duration:12h@Condition:BE(2)-C@GSE16450.MicroAr…  0.183   IFNA    
    ## # … with 8,126,534 more rows

## Calculate scores

Remove any `NA`’s that were formed, as some genes were probably not
detected in each assay. Then calculate the median log2 FC. Probably not
the ideal metric, but should give us a rough idea of the association
between gene-induction and IFN subtype.

``` r
cytosig.IFN %>% 
  na.omit() %>% 
  group_by(Gene, cytokine) %>% 
  summarise(score = median(FC, na.rm = T)) %>% 
  pivot_wider(names_from = cytokine, values_from = score) -> cytosig.IFN

cytosig.IFN
```

    ## # A tibble: 18,292 × 4
    ## # Groups:   Gene [18,292]
    ##    Gene         IFNA     IFNB     IFNG
    ##    <chr>       <dbl>    <dbl>    <dbl>
    ##  1 A1BG     -0.0157  -0.00898 -0.0135 
    ##  2 A1BG-AS1 -0.0436  -0.0275   0      
    ##  3 A1CF     -0.00677 -0.0130   0.00104
    ##  4 A2M      -0.00775  0        0.0640 
    ##  5 A2ML1     0.00521  0        0      
    ##  6 A4GALT   -0.0210  -0.00614 -0.0310 
    ##  7 A4GNT     0       -0.00796 -0.0148 
    ##  8 AAAS     -0.0315  -0.0709  -0.0566 
    ##  9 AACS     -0.0165  -0.0272  -0.0676 
    ## 10 AACSP1    0        0        0      
    ## # … with 18,282 more rows

Lets have a look at the correlations between these signatures.

``` r
ggpairs(cytosig.IFN, columns = c("IFNA", "IFNB", "IFNG"))
```

<img src="3-IFN_I_vs_II_files/figure-gfm/correlations-1.png" style="display: block; margin: auto;" />

## Map mouse to human gene symbols

The CytoSig database was constructed based on human data, our data is a
mouse data-set, can use `biomaRt` to map these to each other, af far as
possible.

Retrieve the mouse gene symbols.

``` r
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

conversion <- getLDS(attributes = c("hgnc_symbol"), 
                     filters = "hgnc_symbol", 
                     values = cytosig.IFN$Gene , 
                     mart = human, 
                     attributesL = c("mgi_symbol"), 
                     martL = mouse, 
                     uniqueRows=T)

tail(conversion, 10)
```

    ##       HGNC.symbol MGI.symbol
    ## 18424       SPTA1      Spta1
    ## 18425      CHRNA9     Chrna9
    ## 18426       FOLR2      Folr2
    ## 18427      RAD23B     Rad23b
    ## 18428       RNF25      Rnf25
    ## 18429      ADARB1     Adarb1
    ## 18430      ZNF582     Gm3854
    ## 18431       EIF5B      Eif5b
    ## 18432    CDC42SE2   Cdc42se2
    ## 18433     ALDH1L2    Aldh1l2

Now import the tumor data.

``` r
scdb_init(here( "Data", "Metacell_files_tumor"), force_reinit=T)
mc.tumor <- scdb_mc("QPCTL_Tum_MC2")

lfp <- as_tibble(log2(mc.tumor@mc_fp), rownames = "Gene" )
```

MetaCell 12 was the cell cluster that exhibited an enhanced IFN
signature.

Join with conversion table and with cytosig data.

``` r
lfp %>% 
  rename(MC12 = `12`) %>% 
  dplyr::select(Gene, MC12) %>% 
  inner_join(conversion, c("Gene" = "MGI.symbol")) %>% 
  inner_join(cytosig.IFN, c("HGNC.symbol" = "Gene")) -> combined.data
```

## Compare signatures

Check the correlations.

``` r
ggpairs(combined.data, columns = c("IFNA", "IFNB", "IFNG", "MC12"), )
```

<img src="3-IFN_I_vs_II_files/figure-gfm/correlations2-1.png" style="display: block; margin: auto;" />

Gene-enrichment in MC12 does not seem to have a stronger correlation
with any of the signatures.

Another way of looking at the data could be to zoom in on the top 100
enriched genes in MC12 and plot the score distribution for each of the
signatures.

``` r
combined.data %>% 
  slice_max(order_by = MC12, n = 100) %>% 
  pivot_longer(cols = contains("IFN"), names_to = "type", values_to = "score") %>% 
  t_test(score~type)
```

    ## # A tibble: 3 × 10
    ##   .y.   group1 group2    n1    n2 statistic    df     p p.adj p.adj.signif
    ## * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <dbl> <chr>       
    ## 1 score IFNA   IFNB     100   100    -0.699  196. 0.485 0.485 ns          
    ## 2 score IFNA   IFNG     100   100    -2.13   182. 0.035 0.104 ns          
    ## 3 score IFNB   IFNG     100   100    -1.45   191. 0.149 0.298 ns

``` r
combined.data %>% 
  slice_max(order_by = MC12, n = 100) %>% 
  pivot_longer(cols = contains("IFN"), names_to = "type", values_to = "score") %>% 
  ggplot(aes(type, score, fill = type))+
  geom_violin(scale = "width")+
  geom_boxplot(color = "black", width = .1, fill = NA)
```

<img src="3-IFN_I_vs_II_files/figure-gfm/plot_sig-1.png" style="display: block; margin: auto;" />

# CytoSig analysis (using webtool)

Another implementation of CytoSig, is providing the webtool with
gene-expression data. The tool will attempt to find gene-expression
patterns that match specific cytokine profiles in their dataset.

To perform this analysis, I’ll first write out the gene-enrichment
values calculated by MetaCell. Should be able to use this as input for
CytoSig. Again, CytoSig is build on human data, so I’ll make the
conversion and provide human gene-symbols.

``` r
lfp %>% 
  inner_join(conversion, c("Gene" = "MGI.symbol")) %>% 
  select(HGNC.symbol, matches("\\d")) %>% 
  write_tsv( here( "Output", "lfp_for_cytosig.tsv"))
```

Used the file generated above as input for
[CytoSig](https://cytosig.ccr.cancer.gov/), downloaded the results files
and stored them in `Data/Cytosig`.

Import the results

``` r
cytosig.files <- list.files(here( "Data", "Cytosig"), pattern = "Cytosig")

here( "Data", "Cytosig", cytosig.files) %>% 
  map(read_csv) %>% 
  set_names(str_extract(cytosig.files, "MC\\d+")) %>% 
  enframe("MetaCell", "data") %>% 
  unnest(cols = data) -> cytosig.results

cytosig.results
```

    ## # A tibble: 516 × 6
    ##    MetaCell ID       Coef  StdErr Zscore Pvalue
    ##    <chr>    <chr>   <dbl>   <dbl>  <dbl>  <dbl>
    ##  1 MC1      IL12   0.0214 0.00408   5.26      0
    ##  2 MC1      LIF   -0.0215 0.00445  -4.80      0
    ##  3 MC1      IL6   -0.0425 0.00428  -9.91      0
    ##  4 MC1      IL4   -0.0235 0.00455  -5.20      0
    ##  5 MC1      IL36  -0.0138 0.00383  -3.60      0
    ##  6 MC1      IL3    0.0163 0.00430   3.81      0
    ##  7 MC1      IL2    0.0216 0.00411   5.25      0
    ##  8 MC1      TGFB1 -0.0150 0.00408  -3.67      0
    ##  9 MC1      TGFB3 -0.0666 0.00401 -16.6       0
    ## 10 MC1      IL15   0.0193 0.00405   4.80      0
    ## # … with 506 more rows

Lets look at the results from all MetaCells, using a stringent P value
cutoff.

``` r
mc.order <- paste0("MC", 1:12)

cytosig.results %>% 
  filter(Pvalue < 0.001) %>% 
  mutate(ID = reorder_within(ID, Coef, MetaCell),
         MetaCell = factor(MetaCell, levels = mc.order)) %>% 
  ggplot(aes(ID, Coef, fill = Coef > 0))+
  geom_bar(stat = "identity", color = "black")+
  facet_wrap(~MetaCell, scales = "free_y")+
  scale_x_reordered()+
  coord_flip()
```

<img src="3-IFN_I_vs_II_files/figure-gfm/plot_cytosig-1.png" style="display: block; margin: auto;" />

``` r
ggsave(here( "Figs", "tumor", "cytosig_all.pdf"), height = 10, width = 12)
```

Or focus on the MetaCell of interest, MC12, plotting all results.

``` r
cytosig.results %>% 
  filter(MetaCell == "MC12") %>% 
  filter(Pvalue < 0.05) %>% 
  mutate(ID = reorder(ID, -Coef)) %>% 
    ggplot(aes(ID, Coef, fill = Coef > 0))+
    geom_bar(stat = "identity", color = "black")+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

<img src="3-IFN_I_vs_II_files/figure-gfm/plot_MC12-1.png" style="display: block; margin: auto;" />

``` r
ggsave(here( "Figs", "tumor", "cytosig_MC12.pdf"), height = 4, width = 6, scale = 1.5)
```
