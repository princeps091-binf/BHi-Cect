#make BPT schematic 
library(readr)
library(Matrix)
library(data.tree)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(GenomicRanges)
library(IRanges)
library(igraph)
library(dplyr)
library(tidyr)
#########################################
#load the chr_spec_res
load(file="path/to/BHi-Cect/out_file")

#build the bpt
chr_bpt<-FromListSimple(chr_spec_res$part_tree)
plot(chr_bpt)
g_bpt<-as.igraph(chr_bpt,directed = T, direction = 'climb')
################################

