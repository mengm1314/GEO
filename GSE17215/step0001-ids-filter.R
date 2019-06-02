rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
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

boxplot(apply(dat,1,mad))
pheatmap::pheatmap(dat[order(apply(dat,1,mad), decreasing = T)[1:50],])
pheatmap::pheatmap(dat[order(apply(dat,1,mad), decreasing = F)[1:50],])


library(corrplot) 
M <- cor( dat ) 
ac=data.frame(g=group_list)
rownames(ac)=colnames(M) 
pheatmap::pheatmap(M,show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = 'heatmap_cor.png')


