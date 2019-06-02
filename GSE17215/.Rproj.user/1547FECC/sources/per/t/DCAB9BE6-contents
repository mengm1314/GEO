rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step1-output.Rdata')
# 每次都要检测数据
dat[1:4,1:4]
## 下面是画PCA的必须操作，需要看说明书。
dat=t(dat)
dat=as.data.frame(dat)
dat=cbind(dat,group_list)
library("FactoMineR")
library("factoextra") 
# The variable group_list (index = ) is removed
# before PCA analysis
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca,repel =T,
             #geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('all_samples_PCA.png')
## 从主成分分析图看出，我们要删掉一个样本，它被混入另外一个组了。
load(file = 'step1-output.Rdata')
colnames(dat)
dat=dat[,-3]
group_list=group_list[-3]
dat[1:4,1:4]
## 下面是画PCA的必须操作，需要看说明书。
dat=t(dat)
dat=as.data.frame(dat)
dat=cbind(dat,group_list)
library("FactoMineR")
library("factoextra") 
# The variable group_list (index = ) is removed
# before PCA analysis
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca,repel =T,
             #geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('choose_samples_PCA.png')





rm(list = ls())  ## 魔幻操作，一键清空~
load(file = 'step1-output.Rdata')
dat[1:4,1:4]
dat[1,]
system.time(apply(dat,1,sd))
system.time(for(i in 1:nrow(dat)){
  sd(dat[i,])
})

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
         annotation_col=ac,filename = 'heatmap_top1000_sd_all.png')

load(file = 'step1-output.Rdata')
colnames(dat)
dat=dat[,-3]
group_list=group_list[-3]
dat[1:4,1:4]
dim(dat)
library(hgu133a.db)
ids=toTable(hgu133aSYMBOL)
head(ids)
dat=dat[ids$probe_id,]
dat[1:4,1:4] 
ids$median=apply(dat,1,median)
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat=dat[ids$probe_id,]
rownames(dat)=ids$symbol
dat[1:4,1:4]  
dim(dat)
save(dat,group_list,file = 'step2-output.Rdata')

