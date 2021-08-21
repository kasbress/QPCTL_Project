library(igraph)
library(RColorBrewer)
library(tidyverse)
library(httr)

# Set working directory
setwd("/DATA/users/k.bresser/QPCTL_bulk/")

# import RNA data 
DE.results <- read_tsv("./tumor/DE_results_tumor.tsv" )

#Species ID, mouse is 10090, human is 9606
species <- "10090"

# set-up genes to search
DE.results %>% 
  filter(FDR < 0.05) %>% 
  pull(Gene.Symbol) -> genes

# number of maximum nodes to return
n.nodes <- 25

# set output filename
out.file <- "./test.png"


## paste genes to search together and make API
genes.search <- paste0(genes, collapse = "%0d")
api <- paste0("https://string-db.org/api/tsv/network?identifiers=", genes.search,"&species=",species , collapse = "")

## Store search in vector and extract content
x <- GET(api)
x <- content(x)

## get a unique linkage DF for mapping the stringIDs to the gene symbols
x %>% 
  dplyr::select("stringId_A", "preferredName_A") %>% 
  rename(stringID = "stringId_A", GeneID = "preferredName_A") %>% 
  group_by(stringID) %>% 
  dplyr::slice(1) -> link.short1
x %>% 
  dplyr::select("stringId_B", "preferredName_B") %>% 
  rename(stringID = "stringId_B", GeneID = "preferredName_B") %>% 
  group_by(stringID) %>% 
  dplyr::slice(1) -> link.short2
link.short1 %>% 
  bind_rows(link.short2) %>% 
  group_by(stringID) %>% 
  dplyr::slice(1) -> link.short


### extract edges that I found IDs for
x %>% dplyr::select(stringId_A, stringId_B, score) -> edges.plot

## create net
net <- graph_from_data_frame(d=edges.plot, vertices=link.short, directed=T)
class(net)


#### Next bit of code is to set the colors of the nodes

# get logFC values
DE.results %>% 
  pull(logFC) -> values
  
## set colors
## Use n equally spaced breaks to assign each value to n-1 equal sized bins 
ii <- cut(values, breaks = seq(min(values), max(values), len = 100), 
          include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(99)[ii]

## Set zero values to darkgrey
zeros <- which(values == 0)
colors[zeros] <- "darkgrey"

## add gene names to color vector
names(colors) <- DE.results$Gene.Symbol

# Only keep colors that are in the linkage, i.e. in the Net
colors <- na.omit(colors[names(colors) %in% link.short$GeneID])

# This bit below fixes any genes that were dropped in the above code for some reason
not.present <- link.short$GeneID[!(link.short$GeneID %in% names(colors))]
not.present <- setNames(rep("darkgrey", length(not.present)),not.present)
colors <- c(colors, not.present)
colors <- colors[link.short$GeneID]


############ Now we can plot the network

# set colors to the Net
V(net)$color <- colors

# scale edges
E(net)$width <- E(net)$score * 6


# Set outfile and plot the network
out.file <- "./tumor/Net_DE_tumor.pdf"
pdf(file = out.file, width = 10, height = 10)
plot(net, edge.arrow.size=0,vertex.label = V(net)$GeneID, vertex.size = 15,
     vertex.label.color = "black", vertex.label.font = 2,vertex.label.cex = .9,
     layout=layout_with_fr, edge.curved=0, )

dev.off()
