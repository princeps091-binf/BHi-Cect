#Spectral clustering
library(readr)
library(caret)
library(Matrix)
library(igraph)
library(RSpectra)
library(dplyr)
library(data.tree)
options(scipen=999999999)
################################################################################
#Building adjacency matrix from the three column output of juicer_tool
Rao_matrix<-function(x){
  x2<-x%>%filter(!(is.na(X3)))
  chr1_in_idx<-seq(1,length(unique(c(x2$X1,x2$X2))))
  names(chr1_in_idx)<-as.character(sort(unique(c(x2$X1,x2$X2))))
  HiCO_in_chr1_mat <- sparseMatrix(i = chr1_in_idx[c(as.character(x2$X1),as.character(x2$X2))],
                                   j = chr1_in_idx[c(as.character(x2$X2),as.character(x2$X1))],
                                   x = c(x2$X3,x2$X3),
                                   dimnames = list(names(chr1_in_idx),names(chr1_in_idx)))
  diag(HiCO_in_chr1_mat)<-NA
  return(HiCO_in_chr1_mat)
}

#power transform
pow_trans<-function(chr1_mat,coef){
  temp_df<-summary(chr1_mat)
  temp_df<-as.data.frame(temp_df)
  preprocessParams <- BoxCoxTrans(temp_df$x,na.rm = T)
  
  transformed <- predict(preprocessParams, temp_df$x)
  check<-chr1_mat
  #only allow positive values and re-scale the difference to order coef
  check@x<-coef*(transformed+(1-min(transformed,na.rm = T)))
  #essential for Laplacian that the Diagonal be set to zero
  diag(check)<-0
  return(check)
}
#############################################################
#Recursive Bi-partitioning
#Laplacian function
lp_fn<-function(x){
  
  Dinv=Diagonal(nrow(x),1/Matrix::rowSums(x))
  
  lp_chr1=Diagonal(nrow(x),1)-Dinv %*% x
  if(dim(lp_chr1)[1] > 10000){
    return(eigs_sym(lp_chr1,k=2,sigma = 0, which='LM',maxitr=10000))
  }
  else{
    temp<-eigen(lp_chr1)
    return(list(vectors=temp[['vectors']][,c(length(temp$values)-1,length(temp$values))],values=temp[['values']][c(length(temp$values)-1,length(temp$values))]))
    
    
  }
}
#Bipartition function
#delineat bi-partitions using kmeans and compute expansion
part_cond_calc<-function(x,reff_g,g_chr1){
  #perform kmeans with 2 clusters in first 2 smallest eigen vector space
  res<-kmeans(x$vectors,2,nstart=5)
  #calculate the expansion of the resulting k clusters
  l_temp_exp<-c()
  sub_g_list<-list()
  for (j in 1:2){
    #create the subnetwork
    sub_g_temp<- induced_subgraph(reff_g,which(res$cluster==j))
    #create cluster label in considered subnetwork
     
    g3<-make_clusters(reff_g,ifelse(V(reff_g)$name %in% V(sub_g_temp)$name,1,2),modularity = F)
    
    #expansion
     l_temp_exp<- c(l_temp_exp,sum(E(reff_g)$weight[which(igraph::crossing(g3,reff_g))])/(sum(graph.strength(sub_g_temp)))) 
    
    #save the members of considered cluster
    temp_name<-paste(length(V(sub_g_temp)$name),length(E(sub_g_temp)),min(as.numeric(unlist(lapply(strsplit(V(sub_g_temp)$name,','),'[',1)))),max(as.numeric(unlist(lapply(strsplit(V(sub_g_temp)$name,','),'[',1)))),sep='_')
    sub_g_list[[temp_name]]<-V(sub_g_temp)$name
    rm(sub_g_temp,g3)
    
  }
  return(list(sub_g_list,l_temp_exp))
}

#loop through all elements in nested list
ff = function(x){ 
  if (class(x) == "list" & length(x)>0) 
    lapply(x, ff) 
  else 
    TRUE
}

#recursive bi-partitioning
spec_bipart<-function(chr1_mat,g_chr1){
  options(scipen=999999999)
  #initialisation
  #container to save cluster hierarchy as list of lists
  chr1_tree_list<-list()
  #containers for cluster member list, cluster conductance/expansion
  chr1_tree_cl<- list()
  
  #whole chromosome laplacian
  lpe_chr1<-lp_fn(chr1_mat)
  #spectral clusters
  res_chr1<-part_cond_calc(lpe_chr1,g_chr1,g_chr1)
  
  #save cluster membership and expansion
  chr1_tree_cl<-list(chr1_tree_cl, res_chr1[1])
  chr1_tree_cl<-do.call(c, unlist(chr1_tree_cl, recursive=FALSE))
  
  #temporary list of candidate cluster to further partition
  ok_part<-names(chr1_tree_cl)
  
  #initiate the tree
  for (i in ok_part){chr1_tree_list[[i]]<-list()}  
  
  #save path to all considered leaves
  lnames <- names(unlist(lapply(chr1_tree_list, ff)))
  names(lnames)<-names(chr1_tree_cl)
  
  #recursive looping
  while(length(ok_part)>0){
    temp_part<-c() 
    for(i in ok_part){
      print(paste(which(ok_part==i),'out of',length(ok_part)))
      #create the subnetwork of considered cluster
      sub_g1<- induced.subgraph(g_chr1,V(g_chr1)$name %in% chr1_tree_cl[[i]])
      #create the corresponding adjacency matrix
      sub_g1_adj<-chr1_mat[rownames(chr1_mat) %in% chr1_tree_cl[[i]],colnames(chr1_mat) %in% chr1_tree_cl[[i]]]
      if(any(colSums(sub_g1_adj)==0)){
        out<-which(colSums(sub_g1_adj)==0)
        sub_g1_adj<-sub_g1_adj[-out,]
        sub_g1_adj<-sub_g1_adj[,-out]
      }
      
      #only consider clusters of at least 3 loci 
      if(length(V(sub_g1)$name)<3){next}
      
      print(paste('eigen decomposition:',i))
      lpe_sub_g1<- lp_fn(sub_g1_adj)
      if(dim(lpe_sub_g1$vectors)[2]<2){next}
      #if resulting partition have theoritical upper bound expansion >1 within subnetwork skip
      print(paste('cheeger:',sqrt(2*Re(lpe_sub_g1$values[1]))))
      #if(sqrt(2*Re(lpe_sub_g1$values[1]))>1){next}
      
      print('kmeans')
      #find actual sub-structures using kmeans on fiedler vector
      res_subg1<-part_cond_calc(lpe_sub_g1,sub_g1,g_chr1)
      
      print('cl append')
      for(k in names(res_subg1[[1]])){chr1_tree_cl[[k]]<-res_subg1[[1]][[k]]}
      
         
      #Only consider future cluster partition if their expansion is majoritarily inside the cluster
      ok_part_temp<-names(res_subg1[[1]])[which(res_subg1[[2]]<1)]
      
      print('tree growth')
      for (j in ok_part_temp){chr1_tree_list[[c(unlist(strsplit(lnames[i],split='\\.')),j)]]<-list()} 
      
      temp_part<-c(temp_part,ok_part_temp)
      rm(sub_g1,sub_g1_adj,lpe_sub_g1,res_subg1)
    }
    print(length(temp_part))
    ok_part<-temp_part 
    #gather updated paths
    lnames <- names(unlist(lapply(chr1_tree_list, ff)))
    #name each path according to leaf of considered path
    names(lnames)<-lapply(lnames,function(x)unlist(strsplit(x,split='\\.'))[length(unlist(strsplit(x,split='\\.')))])
    
  }
  
  return(list(cl_member=chr1_tree_cl,part_tree=chr1_tree_list))  
  
  
  
  
}

################################################################################
chr_set<- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX')
timing<-list()
for(i in chr_set){
  #load chromosome data
  #data input
  chr_dat<- read_delim(paste("D:/K562/5kb/",i,".txt",sep=''),"\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
  
  #process
  print(paste(i,':','sparse matrix'))
  chr_mat<-Rao_matrix(chr_dat)
  print(paste(i,':','power tranform'))
  ptm <- proc.time()
  chr_pow<- pow_trans(chr_mat,1)
  #filter zero sum bins
  if(any(colSums(chr_pow)==0)){
    out<-which(colSums(chr_pow)==0)
    chr_pow<-chr_pow[-out,]
    chr_pow<-chr_pow[,-out]
  }
  
  #save(chr_pow,file=paste("D:/spec_res/robustness/resolution/5kb/",i,"/",i,"_pow_mat.rda",sep=''))
  save(chr_pow,file=paste("D:/spec_res/K562/5kb/",i,"/",i,"_pow_mat.rda",sep=''))
  #save(chr_pow,file=paste("~/../../media/vipin/DISQUEDUR/PhD_jap/HiC_net/data/spec_res/100kb/",i,"/",i,"_pow_mat.rda",sep=''))
  #create the corresponding graph
  print(paste(i,':','graph building'))
  g_chr1<- graph_from_adjacency_matrix(chr_pow,weighted = T)
  #eleminate self loop 
  g_chr1<-delete.edges(g_chr1,E(g_chr1)[which(which_loop(g_chr1))])
  #save(g_chr1,file=paste("D:/spec_res/robustness/resolution/5kb/",i,"/",i,"_gchr.rda",sep=''))
  save(g_chr1,file=paste("D:/spec_res/K562/5kb/",i,"/",i,"_gchr.rda",sep=''))
  
  #save(g_chr1,file=paste("~/../../media/vipin/DISQUEDUR/PhD_jap/HiC_net/data/spec_res/100kb/",i,"/",i,"_gchr.rda",sep=''))
  print(paste(i,':','spectral clustering'))
  chr_spec_res<- spec_bipart(chr_pow ,g_chr1)
  timing[[i]]<-proc.time() - ptm
  #save(chr_spec_res,file=paste("~/../../media/vipin/DISQUEDUR/PhD_jap/HiC_net/data/spec_res/100kb/",i,"/",i,"_spec_res.rda",sep=''))
  #save(chr_spec_res,file=paste("D:/spec_res/robustness/resolution/5kb/",i,"/",i,"_spec_res.rda",sep=''))
  save(chr_spec_res,file=paste("D:/spec_res/K562/5kb/",i,"/",i,"_spec_res.rda",sep=''))
  
  
  
  
}
