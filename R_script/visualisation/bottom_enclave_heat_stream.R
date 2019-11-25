#Enclave pixelmap
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
#create simple edgelist fn knowing that the considered partition
#are mutually exclusive
chr_edgelist<-function(f_clust,g_chr1,chr_spec_res,test){
  fclust_sdf<-data.frame(ego=NA,alter=NA,w=NA,fclust=NA)
  #order f_clust in order of depth
  temp_ix<-sort(unlist(test$Get('level')[f_clust]),index.return=T)$ix
  f_clust_ix<-f_clust[temp_ix]
  #loop through considered partitions
  for (i in f_clust_ix){
    print(i)
    #create subgraph for each f_clust
    s_g1<-induced_subgraph(g_chr1,chr_spec_res$cl_member[[i]])
    #Making sunset plots with ggraph
    #build edge lablelling based on partition depth
    #reorder edges so as top put deeper edges last(front on plot)
    temp_edge<-get.edgelist(s_g1)
    
    sg1_df<-data.frame(ego=temp_edge[,1],alter=temp_edge[,2],w=E(s_g1)$weight)
    
    
    sg1_df$fclust<-rep(i,dim(sg1_df)[1])
    fclust_sdf<-rbind(fclust_sdf,sg1_df)
  }
  fclust_sdf<-fclust_sdf[-1,]
  fclust_sdf$fsize<-as.numeric(unlist(lapply(strsplit(fclust_sdf$fclust,split="_"),"[",2)))
  return(fclust_sdf)
}
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#########################################
#plot params
options(scipen=999999999)
#########################################
#examine tree of choice
chr_set<- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')
res_file<-"~/../../media/vipin/DISQUEDUR/PhD_jap/HiC_net/rebuttal/spec_res/sexton2012/10K/"
res_file<-"~/../../media/vipin/DISQUEDUR/PhD_jap/HiC_net/rebuttal/spec_res/resolution/5kb/"
res_file<-"~/../../media/vipin/DISQUEDUR/PhD_jap/HiC_net/data/spec_res/50kb/"

res<-5000
for(chr_idx in chr_set[-1]){
  chr<-chr_idx
  print(chr)
  #load the contact matrix
  load(paste0(res_file,chr,"/",chr,"_pow_mat.rda"))
  chr_50_mat<-chr_pow
  rm(chr_pow)
  #load the chr_spec_res
  load(file=paste(res_file,chr,"/",chr,"_spec_res.rda",sep=''))
  chr_spec_res_50kb<-chr_spec_res
  rm(chr_spec_res)
  
  #build the bpt
  chr_bpt<-FromListSimple(chr_spec_res_50kb$part_tree)
  #plot(as.dendrogram(chr_bpt),nodePar=list(cex=2,pch=19),center=T,leaflab='none',yaxt='n')
  #plot(chr_bpt)
  #load the chr graph
  load(paste0(res_file,chr,"/",chr,"_gchr.rda"))
  g50_chr<-g_chr1
  rm(g_chr1)
  #include gaps in the heatmap
  full_range<-seq(range(as.numeric(V(g50_chr)$name))[1],range(as.numeric(V(g50_chr)$name))[2],by=res)
  full_g<-g50_chr+vertices(full_range[which(!(full_range %in% V(g50_chr)$name))])
  full_g<-full_g + edges(rep(as.character(full_range[which(!(full_range %in% V(g50_chr)$name))]),each=2))
  E(full_g)$weight[is.na(E(full_g)$weight)]<-0
  
  full_chr_edge<-get.edgelist(full_g)
  full_chr_edge<- data.frame(ego=as.character(full_chr_edge[,1]),alter=as.character(full_chr_edge[,2]),w=E(full_g)$weight,stringsAsFactors = F)
  #update the heatmap to include gaps
  id_conv<-(1:length(unique(c(full_chr_edge$ego,full_chr_edge$alter))))
  names(id_conv)<-as.character(sort(as.numeric(unique(c(full_chr_edge$ego,full_chr_edge$alter)))))
  
  full_chr_edge$ego_id<-id_conv[full_chr_edge$ego]
  full_chr_edge$alter_id<-id_conv[full_chr_edge$alter]
  
  chr_f_mat<-sparseMatrix(i=full_chr_edge$ego_id,full_chr_edge$alter_id,x=full_chr_edge$w,dimnames = list(names(id_conv),names(id_conv)))
  #extract bpt partitions of interest -> domains, subdomains,
  #enclaves
  
  # domains are the largest set of mutually exclusive partitions
  # covering chromosome completely
  #such partition would yield strictly less than 2 children and wouldn't have
  #any ancestor with less than 2 children.
  #get node ancestors
  node_ancestors<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  #eleminate first element of every list member
  node_ancestors<-lapply(node_ancestors,'[',-1)
  #get node children number
  node_children_count<-chr_bpt$Get(function(x) {x$Get('count')})[['Root']]
  
  dom_check<-lapply(node_ancestors[names(which(node_children_count<2))],function(x) all(x %in% names(which(node_children_count==2))))
  dom_part<-names(which(unlist(dom_check)))
  #domains are the 'highest hierarchical level.  Sub-domains are
  #the following split partitions with enclaves representing
  #the most insulated level.
  #find all split partitions
  core_50<-names(which(chr_bpt$Get('count')>1))
  bi_child<-names(which(unlist(lapply(chr_bpt$Get('path'),function(x) x[length(x)-1] %in% core_50))))
  enclave<-bi_child
  if(any(unlist(lapply(chr_bpt$Get(function(x) x$Get('name'))[bi_child],function(x)sum(x %in%chr_bpt$Get('name',filterFun = isLeaf))>1)))){
    enclave<-bi_child[-which(unlist(lapply(chr_bpt$Get(function(x) x$Get('name'))[bi_child],function(x)sum(x %in%chr_bpt$Get('name',filterFun = isLeaf))>1)))]
  }
  
  bi_child_count<-lapply(node_ancestors[bi_child],function(x) sum(x %in% bi_child))
  bcc_df<-data.frame(node=names(bi_child_count),count=unlist(bi_child_count))
  bcc_df<-bcc_df%>%mutate(count_b=ifelse(node %in% enclave,max(bcc_df$count)+1,count))
  bcc_df<- bcc_df%>%arrange(count_b)
  #specific partitions like enclaves
  temp_dom<-dom_part[1]
  dom_bi<-bi_child[which(bi_child %in% unlist(chr_bpt$Get(function(x) x$Get('name'))[temp_dom]))]
  #temporary subset for considered domain
  temp_bcc<-bcc_df%>%filter(node %in% c(temp_dom,dom_bi))
  #check itteratively each domain split level
    print(idx)
    h_level<-idx
    chr_h_level<-as.character(unlist(temp_bcc%>%filter(count==h_level)%>%dplyr::select(node)))
    ext_enclave<-as.character(unlist(temp_bcc%>%filter(count<=h_level & count_b==max(temp_bcc$count_b))%>%dplyr::select(node)))
    chr_h_level<-unique(c(chr_h_level,ext_enclave))
    
  chr_h_level<-enclave
  edge_l_50<- chr_edgelist(chr_h_level,full_g,chr_spec_res_50kb,chr_bpt)
  
  full_edge_list<-full_join(full_chr_edge,edge_l_50,by=c('ego','alter'))
  #or for specific enclave set
  enclave_id_conv<- seq(2,length(chr_h_level)+1)
  #order the f_clust by genome coordinate
  temp_sort<-sort(as.numeric(unlist(lapply(strsplit(chr_h_level,split='_'),'[',3))),index.return=T)$ix
  names(enclave_id_conv)<-chr_h_level[temp_sort]
  #run for every itteration of domain to enclave progression
  test_f<-full_edge_list$fclust
  test_f<-as.factor(test_f)
  levels(test_f)<-enclave_id_conv[levels(test_f)]
  test_f<-as.numeric(as.character(test_f))
  test_f[is.na(test_f)]<-0
  full_edge_list$enclave_id<-test_f
  
  ######################################################
  id_conv<-(1:length(unique(c(full_edge_list$ego,full_edge_list$alter))))
  names(id_conv)<-as.character(sort(as.numeric(unique(c(full_edge_list$ego,full_edge_list$alter)))))
  
  full_edge_list$ego_id<-id_conv[full_edge_list$ego]
  full_edge_list$alter_id<-id_conv[full_edge_list$alter]
  #with individual enclave colours
  sparse_emap<-sparseMatrix(i=full_edge_list$ego_id,j=full_edge_list$alter_id,x=full_edge_list$enclave_id,dimnames = list(names(id_conv),names(id_conv)))
  image(as.matrix(sparse_emap),xaxt= "n", yaxt= "n" ,breaks=c(0,1,enclave_id_conv),col=c('black',gg_color_hue(length(enclave_id_conv))))
    
  
}
