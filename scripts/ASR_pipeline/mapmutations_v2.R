# map mutations onto branches of a tree. Input ancestral sequence reconstruction ('phyDat' object obtained with phangorn), 
# a tree object (class='phylo'), a boolean value indicating whether the program should give the full output, and
# a boolean value indicating whether sites with non-canonical characters should be ignored (highly recommended)
# Output: if 'full' is false then the output is a list with four elements: (1st) 'residue_map' is a mutation matrix 
# for canonical characters (i.e. ,non-gaps,non-ambiguous), where rows correspond to branches as in the tree object,
# the columns are the unique site-patterns in the alignment (as in the 'phyDat' object)
# the values in the mutation matrix are codes (integer-values) that indicate the type of mutation.
# (2nd) 'residue_transition' is a named matrix with the translation of the mutation codes.
# (3rd) 'homoplasy' is a list with length equals to the number of site-patterns, homoplasy[[i]] is a list with the values indicating the index of  
# homoplastic branches for the i-th site-pattern, and names(homoplasy[[i]]) indicates the mutation code. If no homoplasy occurs in the i-th site then homoplasy[[i]] is null. 
# (4th) 'stats' is a data frame with summary statistics
# If 'full' is true then three additional elements are given: (5th) 'is_homoplasy' is a boolean matrix with same dimensions as 'residue_map', indicating mutations 
# that are homoplasies. (6th) 'map' is a mutation matrix that includes
# mutations from/to non-canonical characters. (7th) 'transition_vector' is a translation vector that includes the non-canonical mutations
mapmutations_v2 <- function(asr,t,full,remove_ambiguous){
  nbr<-dim(t$edge)[1]
  nnodes<-max(t$edge)
  ntips<-length(t$tip.label)
  edge<-t$edge
  at<-attributes(asr)
  nstates<-length(at$allLevels)
  transition_vector<-matrix(0:(nstates^2-1),ncol=nstates)
  tfu<-upper.tri(transition_vector,diag = F)
  transition_vector[tfu]<-0
  transition_vector<-t(transition_vector)+transition_vector*-1
  dimnames(transition_vector)<-list(at$allLevels,at$allLevels) # this matrix translates indexes into state transitions
  residue_transitions<-transition_vector[1:at$nc,1:at$nc]
  residue_transitions_index<-as.vector(residue_transitions)
  m<-matrix(unlist(asr),ncol=at$nr,nrow=nnodes,byrow = T)
  nodenames<-append(t$tip.label,paste0('Node',as.character(1:(nnodes-ntips)))) # CHANGED THIS 
  ii<-match(nodenames,at$names)
  if (any(is.na(ii))){
    stop('node names mismatch')
  }
  m<-m[ii,] #Ensure node order is the same in the tree and the 'asr' data
  map<-matrix(0,nrow = nbr,ncol = at$nr)
  residue_map<-matrix(0,nrow = nbr,ncol = at$nr)
  homoplasy<-vector(mode = "list",length = at$nr)
  is_homoplasy<-matrix(FALSE,nrow=nbr,ncol=at$nr)
  for (i in 1:at$nr){
    is_canonical<-m[,i] %in% 1:at$nc
    if(any(!is_canonical)&&remove_ambiguous){
      next
    }
    a<-m[edge[,1],i]
    d<-m[edge[,2],i]
    tf<-a!=d
    # transform matrix-index (two elements: [row,col]) into linear index ('lind')
    # the linear index is the code of the transition as in 'transition_vector'
    # lind = (c-1)*m+r, 
    # where c is the column index, r is the row index, and m is the total number of rows in the matrix
    lind <- (d-1)*nstates+a
    map[tf,i]<-transition_vector[lind[tf]]
    tf<-map[,i] %in% residue_transitions_index
    residue_map[tf,i]<-map[tf,i]
    uv<-unique(residue_map[,i])
    uv<-setdiff(uv,0)
    nuv<-length(uv)
    if(nuv<1){
      next
    }
    homoplasy[[i]]<-vector(mode = "list",length = nuv)
    names(homoplasy[[i]])<-as.character(uv)
    uvm<-matrix(rep(uv,nbr),ncol = nuv,nrow=nbr,byrow = T)
    rmm<-matrix(rep(residue_map[,i],nuv),ncol=nuv,nrow = nbr)
    tf<-uvm==rmm
    s<-apply(tf,1,sum)
    if(any(s>1)){
      stop('More than one mutation assigned to the same branch (!)')
    }
    for (j in 1:nuv){
      ii<-which(tf[,j])
      if(length(ii)<2){
        next
      }
      is_homoplasy[ii,i]<-TRUE
      homoplasy[[i]][[j]]<-ii
    }
    bol<-sapply(homoplasy[[i]], is.null)
    homoplasy[[i]]<-homoplasy[[i]][!bol]
    #if(length(homoplasy[[i]])==0){
    #  homoplasy[[i]]<-NULL
    #}
  }
  stats<-data.frame(n_mutations=sum(residue_map!=0),n_homplasy=sum(is_homoplasy))
  if (full){
    return(list(residue_map,residue_transitions,homoplasy,stats,is_homoplasy,map,transition_vector))
  }
  else{
    return(list(residue_map,residue_transitions,homoplasy,stats))
  }
}
