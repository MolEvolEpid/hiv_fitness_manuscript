tree2distances <- function(t,cid){
  library('phytools')
  ntip<-length(t$tip.label)
  nnode<-Nnode(t)
  root_node<-ntip+1
  intnodes<-root_node:(nnode+ntip)
  d<-dist.nodes(t)
  tf_c<-t$edge[,2]==cid
  if (sum(tf_c)!=1){
    stop('Could not find an unambiguous reference (consensus) branch')
  }
  length_c<-t$edge.length[tf_c]
  anc_c<-t$edge[tf_c,1]
  dis_root2c<-d[anc_c,root_node]
  dis_c2tip<-max(d[setdiff(intnodes,anc_c),root_node])
  tf<-d[,root_node]==dis_c2tip
  b<-which(tf)
  b<-setdiff(b,anc_c)
  if(length(b)==0){
    stop('Could not find B node')
  }
  if(length(b)>1){
    warning('More than one B node possible. Ties were broken randomly')
    b<-b[1]
  }
  dis_c2tip<-d[anc_c,b]
  #desc_nodes<-getDescendants(t,anc_c)
  #desc_nodes<-setdiff(desc_nodes,cid)
  #dis_c2tip<-max(d[anc_c,desc_nodes])
  #b<-which.max(d[anc_c,desc_nodes])
  #b<-desc_nodes[b]
  return(c(length_c,dis_root2c,dis_c2tip,root_node,anc_c,b)) # dists: [JR, FJ, JB]; nodes: [F, J, B]
}