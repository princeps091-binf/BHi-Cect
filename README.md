# BHi-Cect
***A top-down algorithm for identifying the multi-scale hierarchical structure of chromosomes***

---

## Table of Contents

- [Installation](#installation)
- [Features](#features)
- [Usage](#usage)
- [License](#license)

---

BHi-Cect aims to find preferentially self-interacting chromosome clusters using Hi-C data. It only assumes that chromosome architecture consists of a nested hierarchy where we expect to find smaller-sized clusters nested within larger-sized clusters.

---

## Installation

- Copy-paste the different functions present in the BHi-Cect script within your R session
- Implementation as a stand-alone package ongoing

### Setup

#### Required Packages

```r
install.packages(c('readr','caret','Matrix','igraph','RSpectra','dplyr','data.tree'))
```
---

## Features

- Provides a complete hierarchical description of chromosome structure as a tree

---

## Usage

### Input
- Expects a 3-column sparse contact matrix as input
```r
  #load the 3-column sparse matrix Hi-C data
  chr_dat<- read_delim(./path/to/hic_file.txt)
  #convert into sparse matrix format
  chr_mat<-Rao_matrix(chr_dat)
  #power transform
  chr_pow<- pow_trans(chr_mat,1)
  #build the Hi-C interaction graph
  g_chr1<- graph_from_adjacency_matrix(chr_pow,weighted = T)
  #Spectral Clustering
  chr_spec_res<- spec_bipart(chr_pow ,g_chr1)
  
```
---

### Output
- Ouputs a list with 2 elements 
  - a nested list defining the tree representing the hierarchical structure of the input matrix
  - a namesd list specifying the loci contained by all the clusters found by BHi-Cect
  
### Visualisation
In the course of developing the algorithm I have developped a series of visualisation routines to aid the interpretation of the obtained clustering. These are just suggested approaches and the algorithm's output is minimally sufficient to be adapted to any visualisation strategies.

Here is a sample of the routines I use to illustrate the clusters found

- Tree (chr17 IMR90)
![alt text](https://github.com/princeps091-binf/BHi-Cect/blob/master/images/tree.png)

- Cluster map (chr22 IMR90)
![alt text](https://github.com/princeps091-binf/BHi-Cect/blob/master/images/50kb_chr19_bencl.png)

- Nestedness map (chr18 IMR90)
![alt text](https://github.com/princeps091-binf/BHi-Cect/blob/master/images/top_dom_nest_chr18.png)
---

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**


