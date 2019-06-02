rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)

d_h <- function(dat,group_list){
  cg=names(tail(sort(apply(dat,1,sd)),1000))
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F)
  n=t(scale(t(dat[cg,])))
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac)
  
}
d_p <- function(dat,group_list){
  ## 下面是画PCA的必须操作，需要看说明书。
  dat=t(dat)
  dat=as.data.frame(dat)
  dat=cbind(dat,group_list)
  dim(dat)
  library("FactoMineR")
  library("factoextra") 
  # The variable group_list (index = 54676) is removed
  # before PCA analysis
  dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
  fviz_pca_ind(dat.pca,
               geom.ind = c("point", "text"), # show points only (nbut not "text")
               col.ind = dat$group_list, # color by groups
               palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
}

load(file = 'step1-output.Rdata')
dat[1:4,1:4]
dim(dat)
dat=log2(dat)
dat[1:4,1:4]
d_h(dat,group_list)
d_p(dat,group_list)

dat=dat[,-3]
group_list=group_list[-3]    
d_h(dat,group_list)

