rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step2-output.Rdata')
# 每次都要检测数据
dat[1:4,1:4]
table(group_list)

boxplot(dat[1,]~group_list)
t.test(dat[1,]~group_list)
bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
bp(dat[1,])
bp(dat[2,])
bp(dat[3,])
bp(dat[4,])


library(limma)
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
#topTable(fit,coef=2,adjust='BH') 
topTable(fit,coef=2,adjust='BH')
bp(dat['KLK10',])
bp(dat['FXYD3',])
deg=topTable(fit,coef=2,adjust='BH',number = Inf)

head(deg) 
save(deg,file = 'deg.Rdata')

## for volcano 
if(T){
  nrDEG=deg
  head(nrDEG)
  attach(nrDEG)
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df=nrDEG
  df$v= -log10(P.Value)
  ggscatter(df, x = "logFC", y = "v",size=0.5)
  
  df$g=ifelse(df$P.Value>0.01,'stable',
              ifelse( df$logFC >1.5,'up',
                      ifelse( df$logFC < -1.5,'down','stable') )
  )
  table(df$g)
  df$symbol=rownames(df)
  ggscatter(df, x = "logFC", y = "v",size=0.5,color = 'g')
  ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
            label = "symbol", repel = T,
            #label.select = rownames(df)[df$g != 'stable'] ,
            label.select = rownames(head(head(deg))),
            palette = c("#00AFBB", "#E7B800", "#FC4E07") )
  
  ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
  df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                  ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
  table(df$p_c )
  ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
            palette = c("green", "red", "black") )
  
  
}

## for heatmap 
if(T){ 
  load(file = 'step2-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]
  table(group_list)
  load(file = 'deg.Rdata')
  x=deg$logFC
  names(x)=rownames(deg)
  cg=c(names(head(sort(x),100)),
       names(tail(sort(x),100)))
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F)
  n=t(scale(t(dat[cg,])))
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  pheatmap(n,show_colnames =F,
           show_rownames = F,
           cluster_cols = F,
           annotation_col=ac)
  
  
}


### 高级火山图

if(F){
  # 得到两两差异表达的结果 
  x =topTable(fit,coef =2,n=Inf,adjust.method ="BH",sort.by="P")
  save(x,file = 'ggplot_volcano.Rdata')
  head(x)
  dim(x)
  # 选取adj.p.value<0.05且|logFC|>1的基因
  sum(x$adj.P.Val<0.05)
  re.adj = x[x$adj.P.Val <0.05 & (x$logFC > 1 | x$logFC < -1),] 
  
  # 选取p.value<0.05且|logFC|>1的基因
  sum(x$P.Val<0.05)
  re.p = x[x$P.Val <0.05 & (x$logFC > 1 | x$logFC < -1),]   
  
  # 选取adj.p.value and p.value<0.05且|logFC|>1的基因
  # # write.csv(re.adj,"b-a_DEG_limma.re.adj.p.csv",quote =F)
  #  # write.csv(re.p,"b-a_DEG_limma.re.p.csv",quote =F)
  
  # plot
  x$geneID<- rownames(x)
  library(ggplot2)
  library(ggrepel)
  library(ggthemes)
  
  # 簡單的setting for color
  x[,6] <- ifelse((x$P.Value < 0.05 & x$logFC >1.5), "red", ifelse((x$P.Value < 0.05 & x$logFC < -1.5), "blue","grey30"))
  # 複雜的的setting for color
  n1 <- length(x[,1])
  cols <- rep("grey30",n1)
  names(cols)<- rownames(x)
  cols[x$P.Value < 0.05 & x$logFC >1.5]<- "#FB9A99"
  cols[x$P.Value < 0.0001 & x$logFC > 2.5]<- "#ED4F4F"
  cols[x$P.Value < 0.05 & x$logFC < -1.5]<- "#B2DF8A"
  cols[x$P.Value < 0.0001 & x$logFC < -2.5]<- "#329E3F"
  #cols[names(cols)==genelist]<- "red"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)  # set color transparence
  x[,6] <- color_transparent
  
  # 簡單的setting for size
  size <- ifelse((x$P.Value < 0.05 & abs(x$logFC) > 1.5), 4, 2)
  # 複雜的的setting for size
  n1 <- length(x[,1])
  size <- rep(1,n1)
  size[x$P.Value < 0.05 & x$logFC >1.5]<- 2
  size[x$P.Value < 0.0001 & x$logFC > 2.5]<- 4
  size[x$P.Value < 0.00001 & x$logFC > 5]<- 6
  size[x$P.Value < 0.05 & x$logFC < -1.5]<- 2
  size[x$P.Value < 0.0001 & x$logFC < -2.5]<- 4
  size[x$P.Value < 0.00001 & x$logFC < -5]<- 6
  
  # 畫圖位置x，y軸的最大最小位置
  xmin=(range(x$logFC)[1]-(range(x$logFC)[1]+10))
  xmax=(range(x$logFC)[1]+(10-range(x$logFC)[1]))
  ymin<- 0
  ymax<- 12
  #ymax<- range(-log(x$P.Val))[2]+1
  
  ##Construct the plot object
  g <- ggplot(data=x, aes(x=x[,1], y=-log10(x$P.Value))) +
    geom_point(size=size, colour=x[,6]) +
    xlim(c(xmin, xmax)) + ylim(c(ymin,ymax)) +
    xlab("log2 fold change") + ylab("-log10  p-value") +
    geom_vline(xintercept = 1.5, color="grey40", linetype="longdash", size=0.5) +
    geom_vline(xintercept = -1.5, color="grey40", linetype="longdash", size=0.5) +
    geom_vline(xintercept = 2.5, color="grey40", linetype="longdash", size=0.5) +
    geom_vline(xintercept = -2.5, color="grey40", linetype="longdash", size=0.5) +
    geom_hline(yintercept = -log10(0.05), color="grey40", linetype="longdash", size=0.5) +
    geom_hline(yintercept = -log10(0.0001), color="grey40", linetype="longdash", size=0.5) +
    guides(colour = guide_legend(override.aes = list(shape=16)))+
    geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = ifelse(logFC >8, rownames(x),"")),
                    colour="darkred",size = 5,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines")) +
    theme_bw(base_size = 12, base_family = "Times") +
    theme(legend.position="right",
          panel.grid=element_blank(),
          legend.title = element_blank(),
          legend.text= element_text(face="bold", color="black",family = "Times", size=8),
          plot.title = element_text(hjust = 0.8),
          axis.text.x = element_text(face="bold", color="black", size=16),
          axis.text.y = element_text(face="bold",  color="black", size=16),
          axis.title.x = element_text(face="bold", color="black", size=16),
          axis.title.y = element_text(face="bold",color="black", size=16))+
    labs(x="log2 (fold change)",y="-log10 (p-value)",title="Volcano plot")
  
  g
  
  
}





