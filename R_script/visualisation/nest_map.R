#Nestedness maps
#single Domain depth explore
library(readr)
library(Matrix)
library(data.tree)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(igraph)
library(dplyr)
library(tidyr)
#########################################
#chr domains edge listing
chr_dom_edgelist<-function(t_dom,g_chr1,chr_spec_res_50kb,chr_bpt){
  fclust_sdf<-data.frame(ego=NA,alter=NA,w=NA,depth=NA,part=NA)
  for(j in t_dom){
  #create subgraph for considered chromosome domain
  s_g1<-induced_subgraph(g_chr1,chr_spec_res_50kb$cl_member[[j]])
  #build edge lablelling based on partition depth
  core<-names(which(chr_bpt$Get('count')>1))
  bi_child<-names(which(unlist(lapply(chr_bpt$Get('path'),function(x) x[length(x)-1] %in% core))))
  f_clust<-bi_child
  if(any(unlist(lapply(chr_bpt$Get(function(x) x$Get('name'))[bi_child],function(x)sum(x %in%chr_bpt$Get('name',filterFun = isLeaf))>1)))){
    f_clust<-bi_child[-which(unlist(lapply(chr_bpt$Get(function(x) x$Get('name'))[bi_child],function(x)sum(x %in%chr_bpt$Get('name',filterFun = isLeaf))>1)))]
  }
  #extract all Domain subtree partitions
  temp_set<-unlist(chr_bpt$Get(function(x) x$Get('name'))[j])
  #order the partitions by height(lower to higher)
  s_part<-names(sort(chr_bpt$Get('level')[temp_set]))
  
  #assign depth for each partition
  if(j %in% f_clust){
    #level based on tree level of considered partitions
    #core should include every leaf of domain at least
    #garantees all domain partitions will have at least 
    #one core
    bpt_level<-sort(chr_bpt$Get('level')[temp_set])
    #assign to each possible level a normalised value between
    #0(top domain layer) to 1(lowest leaf level)
    temp_depth<-((unique(bpt_level))-min(bpt_level))/(max(bpt_level)-min(bpt_level))
    #assign each of those normalised values to their 
    #corresponding original BPT level
    names(temp_depth)<-as.character(unique(bpt_level))
    #get the vector of normalised BPT level for every
    #domain partition
    depth_conv<-temp_depth[as.character(bpt_level)]
    names(depth_conv)<-names(bpt_level)
    
    
  }
  ##################################################
  if(!(j %in% f_clust)){
    #find the domain part that are not in enclaves
    non_enclave_dom<-s_part[which(!(s_part %in% unlist(chr_bpt$Get(function(x) x$Get('name'))[f_clust])))]
    non_enclave_dom_nl<-(chr_bpt$Get('level')[non_enclave_dom]-chr_bpt$Get('level')[j])/(max((chr_bpt$Get('level')[non_enclave_dom]-chr_bpt$Get('level')[j]))+1)
    #get the domain enclaves
    dom_enclaves<-f_clust[which(f_clust %in% s_part)]
    #get the level of each enclave's split away from the 
    #domain
    enclave_split_level<-lapply(node_ancestors[dom_enclaves],function(x) max(non_enclave_dom_nl[x[which(x %in% non_enclave_dom)]]))
    
    #loop through enclave partitions
    depth_conv<-c()
    for (k in dom_enclaves){
      #assign a normalised level that guarantees the eclave
      #leaf to be equal to 1 and it's shell to be at the 
      #level of the enclave's split away from the domain
      #with every enclave layer evenly distributed between
      #the two
      temp_part<-unlist(chr_bpt$Get(function(x) x$Get('name'))[k])
      if(length(temp_part)==1){depth_conv<-c(depth_conv,1);names(depth_conv)[length(depth_conv)]<-k}
      if(length(temp_part)>1){depth_conv<-c(depth_conv,(1-enclave_split_level[[k]])*(chr_bpt$Get('level')[temp_part]-min(chr_bpt$Get('level')[temp_part]))/(max(chr_bpt$Get('level')[temp_part])-min(chr_bpt$Get('level')[temp_part]))+enclave_split_level[[k]])}
    }
    depth_conv<-c(depth_conv,non_enclave_dom_nl)
    
  }
  
  #get the edgelist of considered partition
  temp_edge<-get.edgelist(s_g1)
  #edge partition vector
  temp_depth<-rep(NA,dim(temp_edge)[1])
  #edge partition vector
  temp_part<-rep(NA,dim(temp_edge)[1])
  
  print(length(s_part))
  for (i in s_part){
    temp_part[which(apply(temp_edge,1,function(x) x[1] %in% chr_spec_res_50kb$cl_member[[i]] & x[2] %in% chr_spec_res_50kb$cl_member[[i]]))]<-i
    
    temp_depth[which(apply(temp_edge,1,function(x) x[1] %in% chr_spec_res_50kb$cl_member[[i]] & x[2] %in% chr_spec_res_50kb$cl_member[[i]]))]<-depth_conv[i]
  }
  
  fclust_sdf<-rbind(fclust_sdf,data.frame(ego=temp_edge[,1],alter=temp_edge[,2],w=E(s_g1)$weight,depth=temp_depth,part=temp_part))
  #sort dataframe so as to have the core edges at the bottom
  fclust_sdf<-fclust_sdf%>%arrange(depth)
  }
  return(fclust_sdf)
}
image.nan <- function(z,  zlim, col, na.color='gray', outside.below.color='black', outside.above.color='white',...)
{
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.below.outside <- zlim[1] - zstep # new z for values below zlim
  newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA
  
  z[which(z<zlim[1])] <- newz.below.outside # we affect newz.below.outside
  z[which(z>zlim[2])] <- newz.above.outside # we affect newz.above.outside
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  
  zlim[1] <- zlim[1] - zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na
  
  col <- c(outside.below.color, col, outside.above.color, na.color) # we construct the new color range by including: na.color and na.outside
  
  image(z=z,  zlim=zlim, col=col, ...) # we finally call image(...)
}


options(scipen=999999999)

#################################################
#load chr of interest
chr='chr18'
#load the contact matrix
load(paste0("~/../../media/vipin/DISQUEDUR/PhD_jap/HiC_net/data/spec_res/50kb/",chr,"/",chr,"_pow_mat.rda"))
chr_50_mat<-chr_pow
rm(chr_pow)
#load the chr_spec_res
load(file=paste("~/../../media/vipin/DISQUEDUR/PhD_jap/HiC_net/data/spec_res/50kb/",chr,"/",chr,"_spec_res.rda",sep=''))
chr_spec_res_50kb<-chr_spec_res
rm(chr_spec_res)

#build the bpt
chr_bpt<-FromListSimple(chr_spec_res_50kb$part_tree)
#load the chr graph
load(paste0("~/../../media/vipin/DISQUEDUR/PhD_jap/HiC_net/data/spec_res/50kb/",chr,"/",chr,"_gchr.rda"))
g50_chr<-g_chr1
rm(g_chr1)
#include gaps in the heatmap
full_range<-seq(range(as.numeric(V(g50_chr)$name))[1],range(as.numeric(V(g50_chr)$name))[2],by=50000)
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

#get node ancestors
node_ancestors<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
#eleminate first element of every list member
node_ancestors<-lapply(node_ancestors,'[',-1)
#get node children number
node_children_count<-chr_bpt$Get(function(x) {x$Get('count')})[['Root']]

dom_check<-lapply(node_ancestors[names(which(node_children_count<2))],function(x) all(x %in% names(which(node_children_count==2))))
dom_part<-names(which(unlist(dom_check)))

test_dom_edge<-chr_dom_edgelist(dom_part,full_g,chr_spec_res_50kb,chr_bpt)
full_edge_list<-left_join(full_chr_edge,test_dom_edge,by=c('ego','alter'))

id_conv<-(1:length(unique(c(full_edge_list$ego,full_edge_list$alter))))
names(id_conv)<-as.character(sort(as.numeric(unique(c(full_edge_list$ego,full_edge_list$alter)))))

full_edge_list$ego_id<-id_conv[full_edge_list$ego]
full_edge_list$alter_id<-id_conv[full_edge_list$alter]

sparse_emap<-sparseMatrix(i=full_edge_list$ego_id,j=full_edge_list$alter_id,x=full_edge_list$depth,dimnames = list(names(id_conv),names(id_conv)))
image.nan(as.matrix(sparse_emap),zlim = c(0,range(sparse_emap@x,na.rm = T)[2]),col = c('black',magma(100)),na.color = 'black',xaxt= "n", yaxt= "n")
png(paste0('~/Documents/github_readme/top_dom_nest','_',chr,'.png'), width =50,height = 50,units = 'mm',pointsize = 2,type='cairo',res=1000)
par(mar=c(0,0,0,0))
image.nan(as.matrix(sparse_emap),zlim = c(0,range(sparse_emap@x,na.rm = T)[2]),col = c('black',magma(100)),na.color = 'black',xaxt= "n", yaxt= "n")
dev.off()