library("phangorn")
asr_v2<- function(t,aln,nsamples){
########################################
  single_asr <- function(t,aln){
  fit<-pml(t,aln,k=4)
  fit<-optim.pml(fit,optBf=T,optNni = F,optEdge = F,optGamma = T,optRate = T,model = "GTR")
  anc<-ancestral.pml(fit, "ml",return="prob")
  return(anc)
  }
########################################
  trees<-NULL
  if(inherits(t,'multiPhylo')){
    nsamples<-1
    trees<-t
  }
  else if (inherits(t,'phylo')){
    trees<-list(t)
  }
  else{
    stop('Tree must be a \'phylo\' or \'multiPhylo\' object')
  }
  if(!inherits(aln,'phyDat')){
    stop('Alignment must be a \'phyDat\' object')
  }
  ntrees<-length(trees)
  list_anc<-vector(mode="list",length = ntrees*nsamples)
  for (i in 1:ntrees){
    t<-trees[[i]]
    anc<-single_asr(t,aln)
    ntips<-length(t$tip.label)
    nnodes<-length(anc)
    at_anc<-attributes(anc)
    ii<-match(names(aln),at_anc$names)
    if(any(is.na(ii))){
      stop('Could not find name(s)')
    }
    n_ii<-length(ii)
    for (idx in 1:n_ii){
      anc[[ii[idx]]]<-aln[[idx]]
    }
    inodes<-setdiff(1:nnodes,ii)
    ninodes<-length(inodes)
    for (idx in 1:nsamples){
      anc2<-anc
      for (j in 1:ninodes){
        prob<-anc2[[inodes[j]]]
        r<-matrix(rep(runif(at_anc$nr,0,1),at_anc$nc),ncol=at_anc$nc)
        cs<-t(apply(prob,1,cumsum))
        tf<-cs>=r
        anc2[[inodes[j]]]<-apply(tf,1,function (x) {match(TRUE,x)})
      }
      list_anc[[i*idx]]<-anc2
    }
  }
  return(list_anc)
}

