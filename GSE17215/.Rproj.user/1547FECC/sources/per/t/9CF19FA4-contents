rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
# 注意查看下载文件的大小，检查数据 
f='GSE17215_eSet.Rdata'

library(GEOquery)
# 这个包需要注意两个配置，一般来说自动化的配置是足够的。
#Setting options('download.file.method.GEOquery'='auto')
#Setting options('GEOquery.inmemory.gpl'=FALSE)
if(!file.exists(f)){
  gset <- getGEO('GSE17215', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
load('GSE17215_eSet.Rdata')  ## 载入数据
class(gset)
length(gset)
class(gset[[1]])
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
a=gset[[1]]
dat=exprs(a)
dim(dat)
#  GPL3921	[HT_HG-U133A] Affymetrix HT Human Genome U133A Array
# Bioconductor - hgu133a.db
dat[1:4,1:4]
dat=log2(dat)
pd=pData(a) 
library(stringr)
group_list =  trimws(str_split(pd$title,',',simplify = T)[,2])
table(group_list)
save(dat,group_list,file = 'step1-output.Rdata')


## 下面的代码是更为详细的注释
######### 导入数据后查看数据情况，行、列名，数据维度、目标数据的大致情况，比如我们想看不同分组的各个基因的表达量，那么你就要查看分组在哪一行/列，怎么分组的，基因是什么情况，只要在掌握已有数据架构的情况下，才能随心所欲的调整表达矩阵（以下代码可以自行增减，达到查看数据情况的目的即可）
 