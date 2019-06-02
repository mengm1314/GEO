rm(list = ls())  ## 魔幻操作，一键清空~
load(file = 'deg.Rdata')
head(deg)
## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
logFC_t=1
deg$symbol=rownames(deg)
deg$g=ifelse(deg$P.Value>0.01,'stable',
            ifelse( deg$logFC > logFC_t,'UP',
                    ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)

DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)
save(DEG,file = 'anno_DEG.Rdata')


gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

geneList=DEG$logFC
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = T)


## KEGG pathway analysis
### 做KEGG数据集超几何分布检验分析，重点在结果的可视化及生物学意义的理解。
if(T){
  ###   over-representation test
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
  source('functions.R')
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  
  ggsave(g_kegg,filename = 'kegg_up_down.png')
  
  ###  GSEA 
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 120,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  
  down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
  up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
  
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = 'kegg_up_down_gsea.png')
  
  
}

### GO database analysis 
### 做GO数据集超几何分布检验分析，重点在结果的可视化及生物学意义的理解。
{
  
  g_list=list(gene_up=gene_up,
              gene_down=gene_down,
              gene_diff=gene_diff)
  
  if(F){
    go_enrich_results <- lapply( g_list , function(gene) {
      lapply( c('BP','MF','CC') , function(ont) {
        cat(paste('Now process ',ont ))
        ego <- enrichGO(gene          = gene,
                        universe      = gene_all,
                        OrgDb         = org.Hs.eg.db,
                        ont           = ont ,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.99,
                        qvalueCutoff  = 0.99,
                        readable      = TRUE)
        
        print( head(ego) )
        return(ego)
      })
    })
    save(go_enrich_results,file = 'go_enrich_results.Rdata')
    
  }
  
  
  load(file = 'go_enrich_results.Rdata')
  
  n1= c('gene_up','gene_down','gene_diff')
  n2= c('BP','MF','CC') 
  for (i in 1:3){
    for (j in 1:3){
      fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
      cat(paste0(fn,'\n'))
      png(fn,res=150,width = 1080)
      print( dotplot(go_enrich_results[[i]][[j]] ))
      dev.off()
    }
  }
  
  
}
### 对 MigDB中的全部基因集 做GSEA分析。
{
  load(file = 'anno_DEG.Rdata')
  geneList=DEG$logFC
  names(geneList)=DEG$symbol
  geneList=sort(geneList,decreasing = T)
  #选择gmt文件（MigDB中的全部基因集）
  d='./MsigDB/symbols'
  gmts <- list.files(d,pattern = 'all')
  gmts
  #GSEA分析
  library(GSEABase) # BiocManager::install('GSEABase')
  gsea_results <- lapply(gmts, function(gmtfile){
    geneset <- read.gmt(file.path(d,gmtfile)) 
    print(paste0('Now process the ',gmtfile))
    egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE)
    head(egmt)
    return(egmt)
  })
  #提取gsea结果
  gsea_results_df <- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_df)
  gseaplot(gsea_results[[2]],'KEGG_CELL_CYCLE') 
  
  
}

### 对 MigDB中的全部基因集 做GSVA分析。
{
  load(file = 'step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]  
  library(hgu133plus2.db)
  ids=toTable(hgu133plus2SYMBOL)
  head(ids)
  dat=dat[ids$probe_id,]
  dat[1:4,1:4] 
  ids$median=apply(dat,1,median)
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]
  ids=ids[!duplicated(ids$symbol),]
  dat=dat[ids$probe_id,]
  rownames(dat)=ids$symbol
  dat[1:4,1:4]  
  
  X=dat
  table(group_list)
  ## Molecular Signatures Database (MSigDb) 
  d='./MSigDB/symbols/'
  gmts=list.files(d,pattern = 'all')
  gmts
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GSVA) # BiocManager::install('GSVA')
  
  if(T){
    es_max <- lapply(gmts, function(gmtfile){ 
      geneset <- getGmt(file.path(d,gmtfile))  
      es.max <- gsva(X, geneset, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
      return(es.max)
    })
    adjPvalueCutoff <- 0.001
    logFCcutoff <- log2(2)
    es_deg <- lapply(es_max, function(es.max){
      table(group_list)
      dim(es.max)
      design <- model.matrix(~0+factor(group_list))
      colnames(design)=levels(factor(group_list))
      rownames(design)=colnames(es.max)
      design
      library(limma)
      contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),
                                     levels = design)
      contrast.matrix<-makeContrasts("TNBC-noTNBC",
                                     levels = design)
      
      contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
      
      deg = function(es.max,design,contrast.matrix){
        ##step1
        fit <- lmFit(es.max,design)
        ##step2
        fit2 <- contrasts.fit(fit, contrast.matrix) 
        ##这一步很重要，大家可以自行看看效果
        
        fit2 <- eBayes(fit2)  ## default no trend !!!
        ##eBayes() with trend=TRUE
        ##step3
        res <- decideTests(fit2, p.value=adjPvalueCutoff)
        summary(res)
        tempOutput = topTable(fit2, coef=1, n=Inf)
        nrDEG = na.omit(tempOutput) 
        #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
        head(nrDEG)
        return(nrDEG)
      }
      
      re = deg(es.max,design,contrast.matrix)
      nrDEG=re
      head(nrDEG) 
      return(nrDEG)
    })
  } 
  
  gmts
  
  save(es_max,es_deg,file='gsva_msigdb.Rdata')
  
  load(file='gsva_msigdb.Rdata')
  
  library(pheatmap)
  lapply(1:length(es_deg), function(i){
    # i=2
    print(i)
    dat=es_max[[i]]
    df=es_deg[[i]]
    df=df[df$P.Value<0.01 & abs(df$logFC) > 0.5,]
    n=rownames(df)
    dat=dat[match(n,rownames(dat)),]
    rownames(dat)=substring(rownames(dat),1,50)
    pheatmap::pheatmap(dat, fontsize_row = 8,height = 11,
                       filename = paste0('gsva_',strsplit(gmts[i],'[.]')[[1]][1],'.pdf'))
  })
  
  adjPvalueCutoff <- 0.001
  logFCcutoff <- log2(2)
  df=do.call(rbind ,es_deg)
  es_matrix=do.call(rbind ,es_max)
  df=df[df$P.Value<0.01 & abs(df$logFC) > 0.5,]
  write.csv(df,file = 'GSVA_DEG.csv')
}





