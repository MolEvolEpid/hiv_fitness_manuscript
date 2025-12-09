get_mutations_by_region <-function(residue_map,weight,residue_transition,is_homoplasy,r,r_state){
  get_sum_by_region<-function(m,r){
    ss1<-apply(m,1,sum)
    z<-tapply(ss1, r, sum)
    ur<-c('FJ','JB','JR','JBD','FJD')
    nur<-length(ur)
    ii<-match(names(z),ur)
    if(any(is.na(ii))){
      stop('region names mismatch')
    }
    ss2<-rep(0,nur)
    ss2[ii]<-z
    names(ss2)<-ur
    return(list(ss2,ss1))
  }
  ncharacters<-dim(residue_transition)[1]
  nsites<-dim(residue_map)[2]
  nbr<-length(r)
  tf<-residue_map==0
  mweight<-matrix(rep(weight,nbr),ncol=nsites,nrow = nbr,byrow = T)
  z<-mweight
  z[tf]<-0
  total_mutations<-get_sum_by_region(z,r)[[1]]
  z[!is_homoplasy]<-0
  total_homoplasy<-get_sum_by_region(z,r)[[1]]
  forward_weight<-matrix(0,nrow = nbr,ncol=nsites)
  backward_weight<-matrix(0,nrow = nbr,ncol=nsites)
  for (i in 1:nsites){
    bol<-r_state[i] %in% 1:ncharacters
    if(!bol){ #skip sites with non-canonical characters in the reference
      next
    }
    forward_mutations<-setdiff(residue_transition[,r_state[i]],0)
    tf<-residue_map[,i] %in% forward_mutations
    if(sum(tf)>0){
      forward_weight[tf,i]<-mweight[tf,i]
    }
    backward_mutations<-setdiff(residue_transition[r_state[i],],0)
    tf<-residue_map[,i] %in% backward_mutations
    if(sum(tf)>0){
      backward_weight[tf,i]<-mweight[tf,i]
    }
  }
  is_forward_mutation<-forward_weight>0
  is_backward_mutation<-backward_weight>0
  s_b_r<-get_sum_by_region(forward_weight,r)
  forward_mutations<-s_b_r[[1]]
  br_forward_mutations<-s_b_r[[2]]
  forward_weight[!is_homoplasy]<-0
  forward_homoplasy<-get_sum_by_region(forward_weight,r)[[1]]
  s_b_r<-get_sum_by_region(backward_weight,r)
  backward_mutations<-s_b_r[[1]]
  br_backward_mutations<-s_b_r[[2]]
  backward_weight[!is_homoplasy]<-0
  backward_homoplasy<-get_sum_by_region(backward_weight,r)[[1]]
  df<-data.frame(total_mutations,total_homoplasy,forward_mutations,forward_homoplasy,backward_mutations,backward_homoplasy,row.names = c('FJ','JB','JR','JBD','FJD'))
  return(list(df,br_forward_mutations,br_backward_mutations))
}