# Target directory, output directory and file extension are all hard-coded
# See comments for possible necessary changes 

library('phytools')

target_dir<-'output/' # input folder
files<-dir(target_dir,pattern = ".+?\\.treefile$", recursive = TRUE) #change file extension here
n<-length(files)
for (i in 1:n){
  t<-read.tree(file=paste0(target_dir,files[i]))
  l<-regmatches(t$tip.label,regexec("gen(\\d+)",t$tip.label))
  l<-as.numeric(sapply(l,FUN=function(x) x[2]))
  if(sum(is.na(l))!=1){
    stop('Some timepoints were not found')
  }
  ms<-min(l,na.rm = T)
  tf<-l==ms
  tf[is.na(tf)]<-FALSE
  n_out<-sum(tf)
  ii<-which(tf)
  m<-NA
  if (length(ii)==1){
    m<-ii
  }
  else{
    m<-getMRCA(t,ii)
  }
  rt<-reroot(t,m)
  write.tree(rt,file=paste0('rooted/',gsub('.*simulations/', '', files[i]))) #output directory
}
