library('phangorn')
source('scripts/ASR_pipeline/asr_v2.R')
source('scripts/ASR_pipeline/mapmutations_v2.R')
source('scripts/ASR_pipeline/tree2regions.R')
source('scripts/ASR_pipeline/get_mutations_by_region.R')
source('scripts/ASR_pipeline/tree2distances.R')
tree_folder<-'rooted/' #folder with trees
ids_tree<-dir(tree_folder,pattern="*treefile") #tree extension file
ids_tree<-paste0(tree_folder,ids_tree)
t<-lapply(ids_tree,read.tree)
aln_folder<- 'output/' #aln folder
ids<-dir(aln_folder,pattern="*aln_w_cons.fasta", recursive = TRUE) #aln extension file
ids<-paste0(aln_folder,ids)
aln<-lapply(ids,function(x) {read.phyDat(x,format = "fasta",type="DNA")})
ntrees<-length(ids)
print(paste0('Number of trees: ', ntrees))
mean_mutation<-vector(mode = "list",length = ntrees)
std_mutation<-vector(mode = "list",length = ntrees)
br_sum_forward<-vector(mode="list",length = ntrees)
br_sum_backward<-vector(mode="list",length = ntrees)
fjb<-vector(mode="list",length = ntrees)
n_anc_sample<-10 # how many random samples for the ancestral reconstructions

for (i in 1:ntrees){
  print(paste0('analysing tree ',as.character(i)))
  t[[i]]$tip.label <- gsub('_1295_', '(1295)', t[[i]]$tip.label)
  cid<-grep('^CON_B.+',t[[i]]$tip.label) #get consensus idx
  if(length(cid)!=1){
    stop('Could not find Consensus')
  }
  fjb[[i]]<-tail(tree2distances(t[[i]],cid),n=3)
  r<-tree2regions(t[[i]],cid)
  anc<-asr_v2(t[[i]],aln[[i]],n_anc_sample)
  at<-attributes(anc[[1]])
  r_state<-anc[[1]][[cid]]
  mut_per_sample<-matrix(NA,ncol=30,nrow = n_anc_sample)
  cn<-NULL
  rn<-NULL
  br_sum_forward[[i]]<-matrix(NA,nrow=dim(t[[i]]$edge)[1],ncol=n_anc_sample)
  br_sum_backward[[i]]<-matrix(NA,nrow=dim(t[[i]]$edge)[1],ncol=n_anc_sample)
  for (j in 1:n_anc_sample){
    print(paste0('sample ',as.character(j)))
    mp<-mapmutations_v2(anc[[j]],t[[i]],TRUE,TRUE)
    mut_by_region<-get_mutations_by_region(mp[[1]],at$weight,mp[[2]],mp[[5]],r,r_state)
    df<-head(mut_by_region,n=1)
    cn<-colnames(df[[1]])
    rn<-rownames(df[[1]])
    mut_per_sample[j,]<-unlist(df)
    br_sum_forward[[i]][,j]<-mut_by_region[[2]]
    br_sum_backward[[i]][,j]<-mut_by_region[[3]]
  }
  mean_mutation[[i]]<-as.data.frame(matrix(apply(mut_per_sample,2,mean),nrow = 5,ncol=6),row.names = rn)
  colnames(mean_mutation[[i]])<-cn
  std_mutation[[i]]<-as.data.frame(matrix(apply(mut_per_sample,2,sd),nrow = 5,ncol=6),row.names = rn)
  colnames(std_mutation[[i]])<-cn
}
save.image(file="mutmap.RData") # output file
