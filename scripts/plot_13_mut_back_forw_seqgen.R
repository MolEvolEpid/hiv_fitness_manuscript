library('phytools')
load('~/Documents/research/replicative_fitness/data_code/plot_seqgen_figure/25_mullins_seqgen_simulations_mutation_map.RData')
vars2keep<-c('t','cid','br_sum_other','br_sum_backward','br_sum_forward','fjb','prefix','vars2keep')
rm(list=setdiff(ls(),vars2keep))
n<-length(vars2keep)
listvars<-vector(mode="list",length=n-1)
for (i in 1:(n-1)){
  listvars[[i]]<-get(vars2keep[i])
}
load('~/Documents/research/replicative_fitness/data_code/plot_seqgen_figure/13_env_mutmap_iqtree.RData')
utid<-s_df$aln_ids
nutid<-length(utid)
tid<-regmatches(prefix,regexpr('[^\\.]+',listvars[[7]]))
foo<-1
########
plot_mut<-function(br2mut,fjb,t,linecolor,lwd){
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
  lines(cs_brlen,mn_s_br2mut,col=linecolor,lwd = lwd)
}

########


pdf('13_fw_bw_seqgen.pdf')
par(mfrow=c(5,3),oma=c(0,0,0,0),mar=c(2,1,1,1))
for (i in 1:nutid)
{
    tf<-tid==utid[i]
    plot(1, type = "n", xlab = "",
         main="",
         ylab = "", xlim = c(0,0.4), 
         ylim = c(0,250))

    ntrees<-sum(tf)
    idx<-which(tf)
    for(j in 1:ntrees){
      plot_mut(listvars[[5]][[idx[j]]],listvars[[6]][[idx[j]]],listvars[[1]][[idx[j]]],"lightskyblue",0.5)
      plot_mut(listvars[[4]][[idx[j]]],listvars[[6]][[idx[j]]],listvars[[1]][[idx[j]]],"lightpink",0.5)
    }
    plot_mut(br_sum_forward[[i]],fjb[[i]],t[[i]],"blue",1)
    plot_mut(br_sum_backward[[i]],fjb[[i]],t[[i]],"red",1)
    legend('topright',utid[i],border=NULL,bty='n')
    if(i==1){
      legend("topleft", legend = c("pat:to","pat:away" ,"sim:to", "sim:away"), bty = "n",
             lwd = c(1,1,1,1), cex = 0.6, col = c( 'blue', 'red',"lightskyblue","lightpink"), lty = 1)
    }
}
dev.off()