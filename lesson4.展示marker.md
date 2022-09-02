在每个cluster注释好细胞类型后，需要找出每个细胞种类的marker基因并进行可视化
# 读取数据、创建seurat对象并进行锚点整合
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

##改变工作路径，将相关结果储存至对应文件夹
setwd("C:/Users/86269/Desktop/shun.C/single_cell")

#加载所需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)

#读取文件,创建Seurat对象
dir=c("CAF1/","CAF2/","dapi1/","dapi2/")
names(dir)=c("CAF1","CAF2","dapi1","dapi2")

scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts,min.cells=3,min.features=300)
}

#归一化处理与寻找高变基因
for(i in 1:length(scRNAlist)){
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]],selection.method='vst',nfeatures=3000)
}

#进行锚点整合
library(future)
scRNA.anchor <- FindIntegrationAnchors(object.list = scRNAlist,anchor.features = 3000)
scRNA1 <- IntegrateData(anchorset = scRNA.anchor)

DefaultAssay(scRNA1) <- "integrated"
```
# 降维聚类 进行SingleR注释
```
#寻找高变基因，降维聚类
scRNA1 <- ScaleData(scRNA1)
scRNA1 <- FindVariableFeatures(scRNA1,Selection.method='vst',nFeatures=3000) 
scRNA1 <- RunPCA(scRNA1,npcs=20,verbose=T)
scRNA1 <- FindNeighbors(scRNA1,reduction='pca',dims=1:7) 
scRNA1 <- FindClusters(scRNA1,resolution=0.4) 
scRNA1 <- RunUMAP(scRNA1,resolution='pca',dims=1:7)
scRNA1<- RunTSNE(scRNA1,dims=1:7)
DimPlot(scRNA1,reduction = 'umap')
#运用SingleR进行注释
library(SingleR)
#准备参考数据与注释数据
library(celldex)
refdata <- MouseRNAseqData() 
testdata <- GetAssayData(scRNA1,slot="data")


clusters <- scRNA1@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata,ref = refdata,labels = refdata$label.main,
                    method='cluster',clusters = clusters,
                    assay.type.test = "logcounts",assay.type.ref = 'logcounts')
celltype <- data.frame(ClusterID=rownames(cellpred),celltype=cellpred$labels,stringsAsFactors = T)



scRNA1@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA1@meta.data[which(scRNA1@meta.data$seurat_clusters==celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}



DimPlot(scRNA1,reduction = 'umap',group.by = 'celltype',label=T)
save(scRNA1,file="scRNA1.Rdata")
```


![image.png](https://upload-images.jianshu.io/upload_images/28382212-e0a5fc47b5fe4e9a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

该图展示了5种细胞种类的聚类情况

# 展示marker
## 1.heatmap
### ROC方式寻找marker基因
```
#更改分组信息
Idents(scRNA1)='celltype'
#找到不同细胞类型的marker
degs <- FindAllMarkers(scRNA1,logfc.threshold = 0.5,
                       test.use = 'roc',
                       return.thresh = 0.25,
                       min.pct = 0.3,
                       only.pos = T)


#对marker进行筛选
degs_sig <- degs %>% filter(pct.1>0.3& power>0.25) %>% filter(cluster != 'others') %>% arrange(cluster,-power)

#筛选出top50
 degs_top50 <- degs_sig %>% group_by(cluster) %>% top_n(50,power) %>% top_n(50,avg_diff) %>% arrange(cluster,-power)
 
 #取平均值,对每个基因在每种细胞内的表达量取平均值，例如A基因在巨噬细胞内的表达平均值与上皮细胞内的表达平均值
 avgData <- scRNA1@assays$RNA@data[degs_top50$gene,] %>% apply(1,function(x){
   tapply(x,scRNA1$celltype,mean)
 }) %>% t

 #对最大值与最小值进行scale与限定,最大值与最小值偏离过大会对画图效果造成影响
 phData <- MinMax(scale(avgData),-2,2)
 
 #可视化
 library(pheatmap)
 rownames(phData)<- 1:nrow(phData)
 phres <- pheatmap(
   phData,
   color = colorRampPalette(c("darkblue","white","red3"))(99),
   scale='row',
   cluster_rows=F,
   cluster_cols=T,#不按行聚类，按列聚类
   clustering_method='complete',
   show_rownames=F,
   annotation_row = data.frame(cluster=degs_top50$cluster)
 )
 ```
![image.png](https://upload-images.jianshu.io/upload_images/28382212-ba871020777c2213.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

### Wilcoxon方法寻找marker基因
```
Idents(scRNA1) ="celltype"
 marker = FindAllMarkers(scRNA1,logfc.threshold = 0.5)
 top10 <- marker %>% group_by(cluster) %>% top_n(n=10,wt = avg_log2FC)
 
 DimPlot(scRNA1,reduction="tsne",group.by = 'celltype',pt.size = 1)
 DoHeatmap(scRNA1,features = top10$gene,label = F,slot="data")
```
![image.png](https://upload-images.jianshu.io/upload_images/28382212-12bdd74f4f2b004e.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![heatmap](https://upload-images.jianshu.io/upload_images/28382212-4cd17d5ae39f248a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


可以观察到热图各个列的宽度差别很大，即各个种类的细胞数量差别很大
接下来每种细胞随机抽取150个再进行热图绘制
```
scRNA1 <- subset(scRNA1,downsample=150)
DoHeatmap(scRNA1,features = top10$gene,label = F,slot="data")
```
![heatmap](https://upload-images.jianshu.io/upload_images/28382212-73af5704aa4ed168.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

### 自定义方法画热图
```
#将barcode加入metadata列中
colanno <- scRNA1@meta.data
#按照细胞类型排序
colanno=colanno%>%arrange(celltype)

#将细胞类型与barcode一一对应
colanno$celltype =factor(colanno$celltype,levels=unique(colanno$celltype))
colanno1 <- colanno[,6]
colanno1<- as.data.frame(colanno1)
rownames(colanno1)=rownames(colanno)
colnames(colanno1)<-'celltype'
```
![heatmap](https://upload-images.jianshu.io/upload_images/28382212-35e7cbb9f594ad23.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


## 2.violin plot
```
library(reshape2)
#每个细胞种类挑选两个marker基因
top2 <- marker %>% group_by(cluster) %>% top_n(n=2,wt = avg_log2FC)
#提取10个marker基因的表达矩阵
Vln.df <- as.data.frame(scRNA1[["RNA"]]@data[top2$gene[1:10],])
Vln.df$gene=rownames(Vln.df)
Vln.df=melt(Vln.df,id='gene')
colnames(Vln.df)[c(2:3)]=c('barcode','exp')

#加入metadata中的信息
anno=scRNA1@meta.data
anno$barcode=rownames(anno)
Vln.df = inner_join(Vln.df,anno,by='barcode')

#可视化
Vln.df %>% ggplot(aes(celltype,exp))+geom_violin(aes(fill=celltype),scale="width")+
  facet_grid(Vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle=45,hjust=1,vjust=1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),#移除网格线
    legend.position = "none"
  )
```
![violin plot](https://upload-images.jianshu.io/upload_images/28382212-afbe5941935ab7ce.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
如果想将颜色调整为与基因的表达量相关，则代码如下
```
#将细胞类型与基因连接组成新的列
Vln.df$celltype_gene=paste(Vln.df$celltype,Vln.df$gene,sep = "_")
#计算每个marker基因在每种细胞内的平均表达量
stat.df=as.data.frame(Vln.df%>%  dplyr::group_by(celltype,gene)%>%dplyr::summarize(mean=mean(exp)))
stat.df$celltype_gene=paste(stat.df$celltype,stat.df$gene,sep = "_")
stat.df=stat.df[,c("mean","celltype_gene")]
#加上其他信息
Vln.df=inner_join(Vln.df,stat.df,by="celltype_gene")
#控制平均数范围防止出现极值
Vln.df$mean=ifelse(Vln.df$mean > 3, 3, Vln.df$mean)

Vln.df%>%ggplot(aes(celltype,exp))+geom_violin(aes(fill=mean),scale = "width")+
  facet_grid(Vln.df$gene~.,scales = "free_y")+
  scale_fill_gradient(limits=c(0,3),low = "lightgrey",high = "red")+
  scale_x_discrete("")+scale_y_continuous("",expand = c(0.02,0))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1)
  )
```
![violin plot](https://upload-images.jianshu.io/upload_images/28382212-45151d94fad9014b.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


## 3.bubble plot
```
#提取10个marker基因的表达矩阵
bubble.df=as.matrix(scRNA1[['RNA']]@data[top2$gene[1:10],])
bubble.df = t(bubble.df)
```
![bubble.df](https://upload-images.jianshu.io/upload_images/28382212-fb0ebc61bf7653a3.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


```
#添加上其他信息
bubble.df=as.data.frame(scale(bubble.df))
bubble.df $CB=rownames(bubble.df)
bubble.df=merge(bubble.df,scRNA1@meta.data[,c(1,6)],by.x="CB",by.y=0)
bubble.df$CB=NULL
```
![bubble.df](https://upload-images.jianshu.io/upload_images/28382212-a2f15fd413dbe199.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
```
celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
#两层循环，第一层针对celltype，第二层针对该celltype内的marker基因
for(i in unique(bubble.df$celltype)){
  bubble.df_small=bubble.df %>% filter(celltype==i)
  for(j in top2$gene[1:10]){
    exp_mean=mean(bubble.df_small[,j])
    exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j]))/length(bubble.df_small[,j])
    celltype_v=append(celltype_v,i)
    gene_v=append(gene_v,j)
    mean_v=append(mean_v,exp_mean)
    ratio_v=append(ratio_v,exp_ratio )
    }

  plotdf=data.frame(
  celltype=celltype_v,
  gene=gene_v,
  exp=mean_v,
  ratio=ratio_v
)
 #可视化
library(RColorBrewer)
plotdf$celltype=factor(plotdf$celltype,levels=sort(unique(plotdf$celltype)))
plotdf$gene=factor(plotdf$gene,levels=rev(as.character(top10$gene[1:10])))
#限定基因表达量最高值
plotdf$exp=ifelse(plotdf$exp>3,3,plotdf$exp)

plotdf %>% ggplot(aes(celltype,gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c('#FFD92F','#FEE391',brewer.pal(11,"Spectral")[7:11])))+
  scale_size_continuous(limits=c(0,1))+theme_classic()+
  theme(
    axis.text.x.bottom = element_text(hjust=1,vjust=1,angle=90)
  )
```
![image.png](https://upload-images.jianshu.io/upload_images/28382212-690fb4f4026b0228.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

## 4.feature plot
以基因Rbp1为例，画TSNE降维结果的聚类图
```
#提取Rbp1基因在每个细胞中的表达量
mat1 = as.data.frame(scRNA1[['RNA']]@data['Rbp1',])
colnames(mat1)='exp'
#提取降维结果并与Rbp1表达量合并
mat2=Embeddings(scRNA1,'tsne')
mat3=merge(mat2,mat1,by='row.names')
#可视化
mat3 %>% ggplot(aes(tSNE_1,tSNE_2))+geom_point(aes(color=exp))+
  scale_color_gradient(low="grey",high="purple")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.ticks=element_blank(),axis.text = element_blank(),
    legend.position = 'none',
    plot.title=element_text(hjust=0.5,size=14)
  )+scale_x_continuous("")+scale_y_continuous("")
```
![feature plot](https://upload-images.jianshu.io/upload_images/28382212-0ab96aa3c79526a3.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

点的颜色越深，基因表达量越高

## 5.H形图
```
#加载数据并进行聚类与注释
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(data.table)
library(tidyverse)

setwd("C:/Users/86269/Desktop/shun.C/single_cell")
load("scRNA_harmony.Rdata")
#进行SingleR注释
library(SingleR)
refdata <- HumanPrimaryCellAtlasData()
testdata <- GetAssayData(scRNA_harmony, slot="data")
clusters <- scRNA_harmony@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
scRNA_harmony@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

DimPlot(scRNA_harmony,group.by = "celltype")
```
![](https://upload-images.jianshu.io/upload_images/28382212-a2e4ec0dc774a63f.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
```
#得到细胞种类
sc.combined=scRNA_harmony
table(sc.combined@meta.data$celltype)
```
![](https://upload-images.jianshu.io/upload_images/28382212-fa399e3e01addc25.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
```
#找出差异基因
type=c("Chondrocytes","Endothelial_cells","Fibroblasts",
       "Macrophage","Monocyte","T_cells",
       "Tissue_stem_cells")
library(future)
r.deg=data.frame()
for (i in 1:7) {
  Idents(sc.combined)="celltype"
  #每个细胞种类在sample11相对于sample3的差异基因
  deg=FindMarkers(sc.combined,ident.1 = "sample_11",ident.2 = "sample_3",
                  group.by = "orig.ident",subset.ident =type[i]   )
  #在deg表中加入细胞类型
  deg$celltype=type[i]
  #在deg表中加入cluster数
  deg$unm=i-1
  #合并到一个表
  r.deg=rbind(deg,r.deg) 
}

```
![r.deg](https://upload-images.jianshu.io/upload_images/28382212-ea31d08a68ccb7cf.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
筛选差异基因
```
#作初步筛选，筛选出差异明显的基因
r.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
#找出上调和下调基因
r.deg$threshold <- as.factor(ifelse(r.deg$avg_log2FC > 0 , 'Up', 'Down'))
#找出显著变化与不显著变化基因
r.deg$adj_p_signi <- as.factor(ifelse(r.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
#合并得到基因的调整情况
r.deg$thr_signi <- paste0(r.deg$threshold, "_", r.deg$adj_p_signi)
r.deg$unm %<>% as.vector(.) %>% as.numeric(.)
```
![r.deg](https://upload-images.jianshu.io/upload_images/28382212-008e5f31cb16033a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
挑选log2FC为top5的基因进行展示
```
#挑选上调基因(每个cluster上调最明显的5个基因)
top_up_label <- r.deg %>% 
  subset(., threshold%in%"Up") %>% 
  group_by(unm) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  as.data.frame()
#挑选下调基因（每个cluster下调最明显的5个基因）
top_down_label <- r.deg %>% 
  subset(., threshold %in% "Down") %>% 
  group_by(unm) %>% 
  top_n(n = -5, wt = avg_log2FC) %>% 
  as.data.frame()
#合并得到差异明显的基因
top_label <- rbind(top_up_label,top_down_label)
top_label$thr_signi %<>% factor(., levels = c("Up_Highly","Down_Highly","Up_Lowly","Down_Lowly"))
```
![](https://upload-images.jianshu.io/upload_images/28382212-c52f7d61ca3f3a8d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
画图准备
```
### 准备绘制暗灰色背景所需数据 
#灰色背景的高度
background_position <- r.deg %>% group_by(unm) %>%summarise(Min = min(avg_log2FC) - 0.2, Max=max(avg_log2FC) + 0.2) %>%as.data.frame()
#灰色背景的宽度
background_position$unm %<>% as.vector(.) %>% as.numeric(.)
background_position$start <- background_position$unm - 0.4
background_position$end <- background_position$unm + 0.4

### 准备绘制中间区域cluster彩色bar所需数据
cluster_bar_position <- background_position
cluster_bar_position$start <- cluster_bar_position$unm - 0.5
cluster_bar_position$end <- cluster_bar_position$unm + 0.5
cluster_bar_position$unm %<>% 
  factor(., levels = c(0:max(as.vector(.))))

## 设置填充颜色
cols_thr_signi <- c("Up_Highly" = "#d7301f",
                    "Down_Highly" = "#225ea8",
                    "Up_Lowly" = "black",
                    "Down_Lowly" = "black")
cols_cluster <- c("0" = "#35978f",
                  "1" = "#8dd3c7",
                  "2" = "#ffffb3",
                  "3" = "#bebada",
                  "4" = "#fb8072",
                  "5" = "#80b1d3",
                  "6" = "#fdb462",
                  "7" = "#b3de69")
p= ggplot() +
  geom_rect(data = background_position, aes(xmin = start, xmax = end, ymin = Min,ymax = Max),
            fill = "#525252", alpha = 0.1) + ###添加灰色背景色
  geom_jitter(data = r.deg, aes(x =unm, y = avg_log2FC, colour = thr_signi),
              size = 1,position = position_jitter(seed = 1)) +
  scale_color_manual(values = cols_thr_signi)+#改变点的颜色
  scale_x_continuous(limits = c(-0.5, max(r.deg$unm) + 0.5),
                     breaks = seq(0, max(r.deg$unm), 1),
                     label = seq(0, max(r.deg$unm),1)) + #修改坐标轴显示刻度

  # 根据top_label标注基因名
  geom_text_repel(data = top_label, aes(x =unm, y = avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), show.legend = F, size = 2.5,
                  box.padding = unit(0, "lines")) +
#绘制彩色bar
  geom_rect(data = cluster_bar_position, aes(xmin = start, xmax = end, ymin = -0.4,ymax = 0.4, fill = unm), color = "black", alpha = 1, show.legend = F) +
 scale_fill_manual(values = cols_cluster) +
  labs(x = "Cluster", y = "average log2FC") +
  theme_bw()

plot1 <- p + theme(panel.grid.minor = element_blank(), ##去除网格线
                   panel.grid.major = element_blank(),
                   axis.text.y = element_text(colour = 'black', size = 14),
                   axis.text.x = element_text(colour = 'black', size = 14, vjust = 67), #调整x轴坐标,vjust的值按照最终结果稍加调整
                   panel.border = element_blank(), ## 去掉坐标轴
                   axis.ticks.x = element_blank(), ## 去掉的坐标刻度线
                   axis.line.y = element_line(colour = "black")) #添加y轴坐标轴

plot1
```
![H形图](https://upload-images.jianshu.io/upload_images/28382212-4656bbcbfcf02c8c.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

## 6.FoldChange比较图
加载package与数据
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

###加载所需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(ggrepel)
library(reshape2)
setwd("C:/Users/86269/Desktop/shun.C/single_cell")
load("scRNA_harmony.Rdata")
library(ggplot2)
library(cowplot)
```
提取数据并寻找差异基因
```
#取cluster数为1的细胞
t.cells <- subset(scRNA_harmony, idents = "1")
#按sample分组
Idents(t.cells) <- "orig.ident"
t.cells=subset(t.cells,ident=c("sample_11","sample_3"))
#放缩
avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))#log1p(x)表示log2(x+1)
avg.t.cells$gene <- rownames(avg.t.cells)
#得出两个sample在cluster1的差异表达基因
deg <- FindMarkers(t.cells,ident.1 = "sample_11",ident.2 = "sample_3", logfc.threshold = 2,group.by = "orig.ident")

```
抽取一部分差异表达基因在图中展示
```
#随机抽取15个基因
genes.to.label = sample(rownames(deg),15)
p1 <- ggplot(avg.t.cells, aes(sample_11,sample_3)) + geom_point() + ggtitle("T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1 
```
![p1](https://upload-images.jianshu.io/upload_images/28382212-5763bd27b89313a7.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
p1中基因label的指向不明确，下面要将被label的基因和其他基因区分开
```
degdf=avg.t.cells
degdf.l=degdf[genes.to.label,]
#提取label基因之外的点
degdf.p=degdf[-match(genes.to.label,rownames(degdf)),]
#增加一列提示基因是否被label
degdf.p$Sig="NO_DIFF"
degdf.l$Sig="DIFF"
degdf=rbind(degdf.p,degdf.l)
```
![avg.t.cells](https://upload-images.jianshu.io/upload_images/28382212-830a171516f98c0a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
![degdf](https://upload-images.jianshu.io/upload_images/28382212-e0c0b052668e58aa.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
接下来进行可视化
```
#用不同颜色将标记的基因与未标记基因区分开
p=ggplot(degdf, aes(sample_11,sample_3)) +
  geom_point( aes(color=Sig)) +
  
  scale_color_manual(values=c("red","grey"))+
  ggtitle("Stromal Cells")+
  theme(plot.title = element_text(size =18,hjust = 0.5, face = "bold")) +
  
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14)) + 
  
  theme_bw()
```
![p](https://upload-images.jianshu.io/upload_images/28382212-62eb409c35400d96.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
```
#去除边框，加深坐标轴
p2 <- p+theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))

#对选择的基因进行标记
data_selected <- degdf[genes.to.label,]
p2 <- p2 + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))

```
![p2](https://upload-images.jianshu.io/upload_images/28382212-b5a46cd9c6ee0115.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
还可以选择展示某个特定的基因,以基因TIMP1为例
```
p +
  geom_point(size = 3, shape =2, data = data_selected[c("TIMP1"),],colour="green") +
  ggrepel::geom_label_repel(
    aes(label =genes.to.label),size =2.5,
    data = data_selected,
    color="black")+theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
```
![展示TIMP1](https://upload-images.jianshu.io/upload_images/28382212-ba59ea0aef37e056.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

## 7.火山图
寻找差异基因
```
cd4.naive=subset(scRNA_harmony,idents ="1")
degdf <- FindMarkers(cd4.naive,ident.1 = "sample_11",ident.2 = "sample_3", logfc.threshold = 0.01,group.by = "orig.ident")
degdf$symbol <- rownames(degdf)
```
![degdf](https://upload-images.jianshu.io/upload_images/28382212-d5aeda0ef75490be.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
设定阈值并筛选差异表达明显的基因
```
logFC_t=0
P.Value_t = 1e-28
#log2FC>0且padj<p则是上调基因,log2FC<0且padj<p则是下调基因，padj>p则是差异表达不明显的基因
degdf$change = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < 0,"down",
                      ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 0,"up","stable"))
```
![degdf](https://upload-images.jianshu.io/upload_images/28382212-4a40531d4d514df2.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
可视化
```
p=ggplot(degdf, aes(avg_log2FC,  -log10(p_val_adj))) +
  geom_point(alpha=0.4, size=2.8, aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("green", "grey","red"))+
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p+theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
```
![火山图
越靠近两侧的基因差异表达越明显](https://upload-images.jianshu.io/upload_images/28382212-42504fc416e006e1.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
```
#展示感兴趣的基因
data_selected <- degdf["TIMP1",]
p + geom_label_repel(data=data_selected,
                     aes(label=rownames(data_selected)))
```


![展示TIMP1](https://upload-images.jianshu.io/upload_images/28382212-0a26303b213db891.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```
#更改点的形状
p +
  geom_point(size = 4, shape =2, data = data_selected,colour="black") +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = data_selected,
    color="green")

```
![](https://upload-images.jianshu.io/upload_images/28382212-5370f5a8fba4acc4.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
