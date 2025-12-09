library('ape')
load('13_env_mutmap_iqtree.RData')
m<-regexec("output/(.+?)/.+?_maxepi(.+?)_repcost(.+?)_conscost(.+?)_run(\\d+)_",ids,perl=T)
zz<-regmatches(ids,m)
tid2par<-do.call(rbind.data.frame, zz)
tid2par<-tid2par[,2:6]
nn<-dim(tid2par)
for(i in 2:nn[2]){
  tid2par[,i]<-as.numeric(tid2par[,i])
}
colnames(tid2par)<-c("tid","maxepi","rfexp","conscost","run")
rm(zz,nn,m)
vars2keep<-c('t','cid','br_sum_backward','br_sum_forward','fjb','tid2par','vars2keep')
rm(list=setdiff(ls(),vars2keep))
n<-length(vars2keep)
listvars<-vector(mode="list",length=n-1)
for (i in 1:(n-1)){
  listvars[[i]]<-get(vars2keep[i])
}

load("13_env_mutmap_iqtree.RData")
utid<-unique(tid2par$tid)
nutid<-length(utid)
df2plot<-NULL
get_mut<-function(br2mut,fjb,t){
  d<-dist.nodes(t)
  nids_fb<-tail(nodepath(t,from=fjb[1],to=fjb[3]),n=-1)
  bids_fb<-match(nids_fb,t$edge[,2])
  cs_brlen<-cumsum(t$edge.length[bids_fb])
  s_br2mut<-br2mut[bids_fb,]
  s_br2mut<-apply(s_br2mut,2,cumsum)
  mn_s_br2mut<-apply(s_br2mut,1,median)
  cs_brlen<-c(0,cs_brlen)
  mn_s_br2mut<-c(0,mn_s_br2mut)
  npoints<-length(cs_brlen)
  return(list(cs_brlen,mn_s_br2mut,npoints))
}

for (i in 1:nutid)
{
  tf<-tid2par$tid==utid[i]
  ntrees<-sum(tf)
  idx<-which(tf)
  for(j in 1:ntrees){
    zz<-get_mut(listvars[[4]][[idx[j]]],listvars[[5]][[idx[j]]],listvars[[1]][[idx[j]]]) #fw
    df2plot<-rbind(df2plot,data.frame(dist_from_F=zz[[1]],n_mut=zz[[2]],tid=rep(tid2par$tid[idx[j]],zz[[3]]),max_epi=rep(tid2par$maxepi[idx[j]],zz[[3]]),rfexp=rep(tid2par$rfexp[idx[j]],zz[[3]]),conscost=rep(tid2par$conscost[idx[j]],zz[[3]]),data_type=rep("simulation",zz[[3]]),mut_type=rep('forward',zz[[3]]),run=rep(tid2par$run[idx[j]],zz[[3]])))
    zz<-get_mut(listvars[[3]][[idx[j]]],listvars[[5]][[idx[j]]],listvars[[1]][[idx[j]]]) #bw
    df2plot<-rbind(df2plot,data.frame(dist_from_F=zz[[1]],n_mut=zz[[2]],tid=rep(tid2par$tid[idx[j]],zz[[3]]),max_epi=rep(tid2par$maxepi[idx[j]],zz[[3]]),rfexp=rep(tid2par$rfexp[idx[j]],zz[[3]]),conscost=rep(tid2par$conscost[idx[j]],zz[[3]]),data_type=rep("simulation",zz[[3]]),mut_type=rep('backward',zz[[3]]),run=rep(tid2par$run[idx[j]],zz[[3]])))
  }
  idx2<-match(utid[i],s_df$aln_ids)
  zz<-get_mut(br_sum_forward[[idx2]],fjb[[idx2]],t[[idx2]])
  df2plot<-rbind(df2plot,data.frame(dist_from_F=zz[[1]],n_mut=zz[[2]],tid=rep(tid2par$tid[idx[j]],zz[[3]]),max_epi=NA,rfexp=NA,conscost=NA,data_type=rep("patient",zz[[3]]),mut_type=rep('forward',zz[[3]]),run=rep(1,zz[[3]])))
  zz<-get_mut(br_sum_backward[[idx2]],fjb[[idx2]],t[[idx2]])
  df2plot<-rbind(df2plot,data.frame(dist_from_F=zz[[1]],n_mut=zz[[2]],tid=rep(tid2par$tid[idx[j]],zz[[3]]),max_epi=NA,rfexp=NA,conscost=NA,data_type=rep("patient",zz[[3]]),mut_type=rep('backward',zz[[3]]),run=rep(1,zz[[3]])))
}
save(df2plot,file="13_env_abm_df2plot.RData")

do_lm<-function(x,y){
  mm <- lm(y ~ 0+x)
  tr_len<-tail(x,n=1)
  if(tr_len!=max(x)){
    stop('error in trunk length assignment')
  }
  return(c(mm$coefficients[1],tr_len))
}
fc<-paste0(df2plot$tid,'_',df2plot$max_epi,'_',df2plot$rfexp,'_',df2plot$conscost,'_',df2plot$data_type,'_',df2plot$mut_type,'_',df2plot$run)
ufc<-unique(fc)
nufc<-length(ufc)
fields<-list('pid','ic','rc','cc','data_type','mutation_type','run','slope','trunk_length')
df<-data.frame(matrix(ncol=length(fields), nrow=nufc, dimnames=list(NULL, fields)))

for (i in 1:nufc){
  tf<-fc==ufc[i]
  df[i,8:9]<-do_lm(df2plot[tf,1],df2plot[tf,2])
  df[i,1:7]<-unlist(strsplit(ufc[i],'_'))
}
ii<-c(2,3,4,7)
for(i in ii){
  print(i)
  df[,i]<-as.numeric(df[,i])
}
save(df,file = "fw_bw_slopes_simulations.RData")
