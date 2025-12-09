tree2regions <- function(t,cid){
  library('ape')
  source('scripts/ASR_pipeline/tree2distances.R')
  nbr<-dim(t$edge)[1]
  regions<-vector(mode = "character",length = nbr)
  fjb<-tail(tree2distances(t,cid),n=3)
  f_j<-nodepath(t,from=fjb[1],to=fjb[2])
  j_b<-nodepath(t,from=fjb[2],to=fjb[3])
  for (i in 1:nbr){
    ad<-t$edge[i,]
    if(ad[2]==cid){
      regions[i]<-'JR'
    }
    else if(all(ad %in% f_j)){
      regions[i]<-'FJ'
    }
    else if(all(ad %in% j_b)){
      regions[i]<-'JB'
    }
    else if(sum(nodepath(t,from=ad[2],to=fjb[2]) %in% f_j)>1){
      regions[i]<-'FJD'
    }
    else {
      regions[i]<-'JBD'
    }
  }
  return(regions)
}
