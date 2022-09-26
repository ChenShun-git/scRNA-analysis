
# 1.AUCell package
假设存在一通路A有20个相关基因，要评价通路的表达情况，可以用AUCell对这20个相关基因进行打分，得分高代表通路表达较高。
#### ①加载package与数据
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())


setwd("C:/Users/86269/Desktop/shun.C/single_cell")
###加载所需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(ggplot2)
load("scRNA1.rdata")
library(AUCell)
library(ggplot2)
library(Seurat)
library(clusterProfiler)
```
#### ②对细胞进行排序
```
cells_rankings <- AUCell_buildRankings(scRNA1@assays$RNA@data,nCores=1, plotStats=TRUE)
```
![cell_rankings](https://upload-images.jianshu.io/upload_images/28382212-4bc0b2483694c777.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
这一步是对每单个细胞中每个基因的表达量进行排序，每个细胞的每个基因都得到一个排序序号
#### ③下载用于参考的细胞通路对应的基因
在GSEA broad上进行下载
![GSEA broad](https://upload-images.jianshu.io/upload_images/28382212-1c5398c0191a0db2.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

**通路数据库简介**

1.GO（gene ontology，即基因本体论）：基因本体论是对基因在不同维度和不同层次上的描述

a. Cellular Component(CC):解释基因在细胞中的位置

b. Biological Process(BP):说明该基因参与了哪些生物学过程

c. Molecular Function(MF):说明该基因在分子层面的功能

2.KEGG:对细胞所有pathway进行富集分析,KEGG数据库共有186条pathway

**二者的区别在于数据库不同**



下载KEGG pathway数据后进行读取
```
c2 <- read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt") 
```
![c2](https://upload-images.jianshu.io/upload_images/28382212-11eeb797fe3f38b1.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

提取每个通路对应的基因
```
geneSets <- lapply(unique(c2$term), function(x){print(x);c2$gene[c2$term == x]})
names(geneSets) <- unique(c2$term)
```
![geneSets](https://upload-images.jianshu.io/upload_images/28382212-4199017209e3635b.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
#### ④基因转换
cell_rankings 是鼠的基因，而geneSets是人的pathway数据，需要用msigdbr包将geneSets转换成对应的鼠的通路数据
```
library(msigdbr)
m_df<- msigdbr(species = "Mus musculus",  category = "C2", subcategory = "KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
```
![fgsea_sets](https://upload-images.jianshu.io/upload_images/28382212-ba0eff5c32483e85.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
#### ⑤基因富集计算
```
cells_AUC <- AUCell_calcAUC(fgsea_sets, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)
#提取氧化磷酸化通路相关的基因
geneSet <- "KEGG_OXIDATIVE_PHOSPHORYLATION"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
#合并metadata信息与降维结果
df<- data.frame(scRNA1@meta.data, scRNA1@reductions$umap@cell.embeddings)
#将seurat cluster信息注释到图上
class_avg <- df %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
```
#### ⑥可视化
```
p <-ggplot(df, aes(UMAP_1, UMAP_2))  +
  geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="E") +
  ggrepel::geom_label_repel(aes(label = seurat_clusters),
                            data = class_avg,
                            size = 5,
                            label.size = 1,
                            segment.color = NA
  )+   theme(legend.position = "none") + theme_bw()


```
![氧化磷酸化通路活性](https://upload-images.jianshu.io/upload_images/28382212-16805f49ba2631bf.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
分组展示
```
p<- p+facet_grid(.~orig.ident)
```
![p](https://upload-images.jianshu.io/upload_images/28382212-fdc2f5968b9ab6bb.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


# 2.过表征分析
## ①基于GO数据库的分析
接下来将对DepiNeg1、DepiNeg2两个样本的差异基因进行分析，探究它们的差异基因主要与哪些功能有关
**找到差异基因**
```

degdf <- FindMarkers(scRNA1,ident.1 = "DapiNeg1",ident.2 = "DapiNeg2", 
                     logfc.threshold = 0.3,group.by = "orig.ident",subset.ident=1)

```
**分析BP类通路**
```
library(org.Mm.eg.db)
degs.list=rownames(degdf)
erich.go.BP = enrichGO(gene =degs.list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
dotplot(erich.go.BP,showCategory = 15 )
barplot(erich.go.BP,showCategory = 8)
erich.go.BP=erich.go.BP@result 

```
![BP dotplot](https://upload-images.jianshu.io/upload_images/28382212-d81404e1ca16508b.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
![BP barplot](https://upload-images.jianshu.io/upload_images/28382212-aeb1340077d4886a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


**分析CC类通路**
```
erich.go.CC = enrichGO(gene =degs.list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.CC,showCategory = 8)
barplot(erich.go.CC,showCategory = 8)
```
![CC dotplot](https://upload-images.jianshu.io/upload_images/28382212-ced0858ddf41de05.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
![CC barplot](https://upload-images.jianshu.io/upload_images/28382212-946556019222be0f.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


**分析MF类通路**
```
erich.go.MF = enrichGO(gene =degs.list,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.MF,showCategory = 8)
barplot(erich.go.MF,showCategory = 8)
```
![MF dotplot](https://upload-images.jianshu.io/upload_images/28382212-5c0ba058371cd569.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
![MF barplot](https://upload-images.jianshu.io/upload_images/28382212-f9e38fd6c13f58e6.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

## ②基于KEGG数据库分析
**转换基因名**
将差异基因的基因名从symble名转化为entrezID
```
DEG.entrez_id = mapIds(x = org.Mm.eg.db,
                       keys =  degs.list,
                       keytype = "SYMBOL",
                       column = "ENTREZID")
```
**进行基于KEGG的富集**
```
enrich.kegg.res <- enrichKEGG(gene = DEG.entrez_id,
                             organism = "mmu",
                             keyType = "kegg")
dotplot(enrich.kegg.res)
```
如果出现报错则运行以下代码
```
library(R.utils)
R.utils::setOption('clusterProfiler.download.method','auto')
```
![KEGG富集 dotplot](https://upload-images.jianshu.io/upload_images/28382212-954243e179bfd812.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
**将结果进行保存**
```
result <- enrich.kegg.res@result
write.table(result,"enrich.kegg.res.txt",sep = "\t",col.names = NA)
```
![enrich.kegg.res](https://upload-images.jianshu.io/upload_images/28382212-d39b0672293c1c96.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

**挑选其中感兴趣的通路进行作图**
以第12-16，21-23，27-29个通路为例，将这些通路制成一个txt文件
，并读取该文件进行作图
```
kegg=read.table("KEGG.8.20.txt",header = T,sep = "\t")
```
![kegg](https://upload-images.jianshu.io/upload_images/28382212-9f672f56f89dcec5.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


**可视化**
```
k=data.frame(kegg)
library(ggplot2)
library(dplyr)
#提取GeneRatio的前与后的值
before <- as.numeric(sub("/\\d+$", "", k$GeneRatio))
after <- as.numeric(sub("^\\d+/", "", k$GeneRatio))
k$GeneRatio = before /after
font.size =10
k %>% 
  ## 按p值排序
  arrange(p.adjust) %>% 
  ## 开始ggplot2 作图，其中fct_reorder调整因子level的顺序
  ggplot(aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  ## 画出点图
  geom_point(aes(color=p.adjust, size = Count)) +
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 调整泡泡的大小
  scale_size_continuous(range=c(3, 8))+
  ## 如果用ylab("")或出现左侧空白
  labs(y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

```
![](https://upload-images.jianshu.io/upload_images/28382212-9647bc050b31e0f2.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

### ③可视化展示富集结果
**读取数据，完成降维聚类、注释**
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

setwd("C:/Users/86269/Desktop/shun.C/single_cell")
load("scRNA_harmony.Rdata")
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
load("ref_Human_all.RData")


refdata <- ref_Human_all

testdata <- GetAssayData(scRNA_harmony, slot="data")
clusters <- scRNA_harmony@meta.data$seurat_clusters
###开始用singler分析
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
scRNA_harmony@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
```
**以cluster0为例，找到cluster0的marker基因并进行GO富集**
```
#对cluster0寻找marker基因
deg =FindMarkers(scRNA_harmony,ident.1 = "sample_11",ident.2 = "sample_3",group.by = "orig.ident",
                 subset.ident = "0")
degs.list=rownames(deg)
#GO BP富集
library(org.Hs.eg.db)
library(clusterProfiler)
erich.go.BP = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

dotplot(erich.go.BP,showCategory = 30)
```
![富集结果dotplot](https://upload-images.jianshu.io/upload_images/28382212-f8bbb71d76038539.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

**挑选感兴趣的通路进行作图**
```
erich.go.BP=erich.go.BP@result
```
![enrich.go.BP](https://upload-images.jianshu.io/upload_images/28382212-078c5c44aad15468.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

挑选3个感兴趣的通路信息制成txt文件
![goplot.txt](https://upload-images.jianshu.io/upload_images/28382212-571a90974da8b8ac.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
```
#提取用于作图的信息
library(GOplot)
go1= read.table("goplot.txt",sep = "\t",header = T)
go1=erich.go.BP[1:3,c(1,2,8,6)]
rownames(go1)=go1$ID
go=go1
```
![go](https://upload-images.jianshu.io/upload_images/28382212-49994e7432189f80.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```
#调整go数据框
library(stringr)
go$geneID=str_replace_all(go$geneID,"/",",")
names(go)=c('ID','Term','Genes','adj_pval')
go$Category="BP"
```
![go](https://upload-images.jianshu.io/upload_images/28382212-b9ff47e739248237.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```
#提取每个通路对应的基因
x1=strsplit(go$Genes[1],split=",",fixed=T)
x2=strsplit(go$Genes[2],split=",",fixed=T)
x3=strsplit(go$Genes[3],split=",",fixed=T)
g1=c(x1[[1]],x2[[1]],x3[[1]])
#从差异分析结果deg中提取这些基因的信息并整合
genedata1=deg[g1,]   
genedata1$ID=rownames(genedata1)
genedata2=data.frame(ID=genedata1$ID,logFC=genedata1$avg_log2FC)
```
![genedata1](https://upload-images.jianshu.io/upload_images/28382212-0f1871dffce7cd0b.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
![genedata2](https://upload-images.jianshu.io/upload_images/28382212-7fe36017d74c099a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

**各种形式可视化**
```
circ <- circle_dat(go,genedata2)
##条形图
GOBar(subset(circ, category == 'BP'))
#气泡图
GOBubble(circ, labels = 3)
#圈图
GOCircle(circ, nsub = 3)
```
![Bar Plot](https://upload-images.jianshu.io/upload_images/28382212-ad6380fbe3150a72.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
![bubble plot](https://upload-images.jianshu.io/upload_images/28382212-b897710f3389d180.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


![circ plot](https://upload-images.jianshu.io/upload_images/28382212-fc68212917263495.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```
#和弦图
chord<-chord_dat(circ, genedata2)
GOChord(chord, gene.order = 'logFC')
#基因与GO Term的热图
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))
```
![chord plot
](https://upload-images.jianshu.io/upload_images/28382212-e48d32a2a613e64d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
![heatmap](https://upload-images.jianshu.io/upload_images/28382212-2f46bbc7bfd106cd.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

# 3.GSVA富集
GSVA富集不仅可以得出差异基因富集于哪个通路，还能分析出该通路相关的基因综合表达是上调还是下调。GSVA先将基因富集到各通路上，再通过差异分析找出差异表达的通路。
![GSVA原理](https://upload-images.jianshu.io/upload_images/28382212-15ef8d11bc7cc446.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

### ①pipeline
输入表达矩阵，得到GSVA评分矩阵
**准备分析数据与package**
```
library(GSEABase)
library(GSVA)
library(msigdbr)

###准备参考数据集###
m_df<- msigdbr(species = "human",  category = "H" )
#将gene_symbol按gs_name分组
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
###准备分析数据集###
exp=AverageExpression(scRNA_harmony,add.ident = "orig.ident") 
exp=exp[["RNA"]]
#挑选sample11与sample3的cluster0、1、2进行分析
counts2=exp[,c(1:6)]
```
![geneSets](https://upload-images.jianshu.io/upload_images/28382212-a60d5c98edfa497c.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![counts2](https://upload-images.jianshu.io/upload_images/28382212-91d1e13950f78be6.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

**进行GSVA富集**
```
GSVA_hall <- gsva(expr=as.matrix(counts2), 
                  gset.idx.list=geneSets, 
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Poisson", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"
                   )
           
```

![GSVA_hall](https://upload-images.jianshu.io/upload_images/28382212-61c4ecf6dadaaefe.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

**limma进行差异通路分析**
```
# 设置分组矩阵
group <- factor(c(rep("s.11",3), rep("s.3", 3)), levels = c( 's.11','s.3'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall)
```
![design](https://upload-images.jianshu.io/upload_images/28382212-cf72b977a572bea4.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```
###进行差异分析###
#创建比较对象，此处是s.3与s.11比较
compare <- makeContrasts(s.3 - s.11, levels=design)
#线性拟合
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
#贝叶斯
fit3 <- eBayes(fit2)
#提取差异分析结果
Diff <- topTable(fit3, coef=1, number=200)
```
![Diff](https://upload-images.jianshu.io/upload_images/28382212-e9405de49a46742e.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

**可视化**
```
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
# 去掉"HALLMARK_"
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$t  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 按t值排序
dat_plot <- dat_plot %>% arrange(t)
```
![dat_plot](https://upload-images.jianshu.io/upload_images/28382212-56fea16f8c3eb8d8.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```
library(ggplot2)
library(ggthemes)
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +#横纵坐标对调
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, S.11 VS S.3') + 
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
###添加标签###
# 小于-1的数量
low1 <- dat_plot %>% filter(t < -1) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于1总数量
high0 <- dat_plot %>% filter(t < 1) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black',size=2) + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey',size=2) + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey',size=2) + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black',size=2) # 大于1的为黑色标签
p
```
![p](https://upload-images.jianshu.io/upload_images/28382212-aef6b8b27f1f581c.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```
pheatmap::pheatmap(GSVA_hall, #热图的数据
                   cluster_rows = T,#行聚类
                   cluster_cols =T,#列聚类，可以看出样本之间的区分度
                   fontsize = 5,
                   show_colnames=T,
                   scale = "row", #以行来标准化
                   color =colorRampPalette(c("#FF7744", "white","#AAAAAA","#0044BB"))(100))

```
![heatmap](https://upload-images.jianshu.io/upload_images/28382212-f664b262d9f79746.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

# 4.自定义通路
某些情况下直接用GO、KEGG、hallmark等数据库已有的基因集进行富集分析得不到理想结果，就可以自定义通路与通路对应的基因进行GSVA分析（通过文献或已有资料删除对通路起负反馈效果的基因等）。
### ①GSVA自定义分析
**准备自定义的通路数据集**
需要准备1.基因名称的txt文件,2.通路名称文件
先在excel表中进行操作，再复制粘贴到txt中
在excel表中准备数据时需要多加一个自定义的xyz通路，并保证该通路基因数最多，因为在后续去掉表中的空值时最长的一列（即没有空值的一列）会被全部去除

```
#读取文件
gene.set=read.table("gene.txt",
                    header =F,sep = '\t',quote = '',fill=TRUE)
kegg.123=read.table("genesetname.txt",
                    header =F,sep = '\t',quote = '')
#将行转为基因，列转为通路
gene.set1=as.matrix(gene.set)
gene.set2=t(gene.set1)
#去掉空值
gmt=list()
for (i in 1:11) {
  y=as.character(gene.set2[,i]) 
  b=which(y=="")
  gmt[[i]]=y[-b]
}
#对通路命名
names(gmt)=kegg.123$V1
#去掉xyz通路
gmt=gmt[-1]
```
![gmt](https://upload-images.jianshu.io/upload_images/28382212-7e27a7d4d0e2484f.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

### ②GO富集自定义分析
**找出差异基因**
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

setwd("C:/Users/86269/Desktop/shun.C/single_cell")
load("scRNA_harmony.Rdata")
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
load("ref_Human_all.RData")


refdata <- ref_Human_all

testdata <- GetAssayData(scRNA_harmony, slot="data")
###把scRNA数据中的seurat_clusters提取出来，注意这里是因子类型的
clusters <- scRNA_harmony@meta.data$seurat_clusters
###开始用singler分析
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
scRNA_harmony@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

sc.t=scRNA_harmony[,rownames(subset(scRNA_harmony@meta.data,celltype=="T_cells"))]  
deg =FindMarkers(scRNA_harmony,ident.1 = "sample_11",ident.2 = "sample_3",group.by = "orig.ident",
                 subset.ident = "0")

```
**进行自定义的GO分析**
用enricher()进行自定义分析，enrichGO()的分析数据都来源于网络上的数据库，无法修改，但enricher()依靠自定义的通路数据，方便修改。
enricher的通路数据TERM2GENE是包含两列的数据框，一列为通路，一列为基因
```
#准备通路数据
m_df <-msigdbr(species = "human",category = "H")
m_df <- m_df[,c(3,4)]
#进行富集分析
res<-enricher(gene = rownames(deg),TERM2GENE = m_df)
```

# 5.代谢活性通路富集
仅关注与代谢有关的通路
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

setwd("C:/Users/86269/Desktop/shun.C/single_cell")

install.packages("wesanderson")
install.packages("C:/Users/86269/Desktop/shun.C/single_cell/phangorn_2.8.1.zip", repos = NULL, type = "win.binary")
install.packages("C:/Users/86269/Desktop/shun.C/single_cell/quadprog_1.5-8.zip", repos = NULL, type = "win.binary")
devtools::install_github("YosefLab/VISION")
devtools::install_github("wu-yc/scMetabolism")


library(scMetabolism)
library(ggplot2)
library(rsvd)

load(file = "pbmc_demo.rda")

countexp.Seurat<-sc.metabolism.Seurat(obj = countexp.Seurat, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")



DimPlot.metabolism(obj = countexp.Seurat, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1)

##Dot plot
input.pathway<-c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
DotPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "ident", norm = "y")



## Box plot
BoxPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "ident", ncol = 1)

metabolism.matrix<-sc.metabolism(countexp = countexp, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")

```
![dimplot](https://user-images.githubusercontent.com/112565216/188271408-b9fc4abe-5da6-424f-9d01-08efc3155a77.png)
![dotplot](https://user-images.githubusercontent.com/112565216/188271428-116ee4ba-2f64-4427-98ca-5c875630f29d.png)
![boxplot](https://user-images.githubusercontent.com/112565216/188271489-5752f98d-bcaa-4f67-b9e7-2d0fc2a1b7a3.png)

# 6.GSEA富集分析

传统的KEGG或GO富集分析的缺点是规定阈值进行筛选，但一条通路上不同的基因表达量可能会有很大差异，上游基因的微小变化可能造成下游基因明显变化，因此规定阈值进行筛选会损失很多基因信息。而GSEA方法不需要规定阈值进行筛选，规避了传统富集方法的缺点。
# ①Pipeline
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

setwd("C:/Users/86269/Desktop/shun.C/single_cell")
###加载所需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(ggplot2)
load("scRNA1.rdata")
#按细胞类型分组，寻找差异基因
Idents(scRNA1)="celltype"
deg=FindMarkers(scRNA1,ident.1 = "Adipocytes",ident.2 = "Granulocytes",
                min.pct = 0.01,logfc.threshold = 0.01)
#将avg_log2FC按从大到小提取（GSEA函数的要求）
genelist=deg$avg_log2FC
#把鼠的基因名改为人的基因名，二者的差别就在于大小写不同
names(genelist)=toupper(rownames(deg))
genelist=sort(genelist,decreasing = T)
###进行富集分析与可视化###
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(enrichplot)
#读取通路信息
geneset=read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
#进行GSEA富集
egmt=GSEA(genelist,TERM2GENE = geneset,
          minGSSize = 1,pvalueCutoff = 0.5)#TERM2GENE可以自行筛选
#查看前两条通路的富集情况
gseaplot2(egmt,geneSetID =c(1,2),pvalue_table = T,base_size = 5)
```
![富集结果](https://upload-images.jianshu.io/upload_images/28382212-26aacc659ebf6871.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

第一部分为基因Enrichment Score的折线图，大概意思是对差异基因排序列表进行扫描，每遇到一个参与该通路的差异基因，得分就会上升，如果遇到一个差异基因不在该通路里，就罚分；score纵轴为对应的Running ES, 在折线图中有个峰值，该峰值就是这个基因集的Enrichemnt score，峰值之前的基因就是该基因集下的核心基因，也称为领头亚集基因，表明对该条通路最终得分贡献最大的基因。
第二部分用线条标记位于该基因集下的基因，从左到右按照log2FC从高到低进行排序。
第三部分的纵坐标代表FoldChange值。
从图中可看出一个通路的基因富集在FoldChange>0的区域，说明该通路上调，另一个则富集在FoldChange<0的区域，说明该通路下调。
**将上调与下调的通路分别进行批量可视化并储存图片**
```
kegg.res=egmt@result

down.kegg.res<-kegg.res[(kegg.res$p_adjust<0.01 & kegg.res$enrichmentScore < -1.5),]
down.kegg.res$group=1

up.kegg.res<-kegg.res[(kegg.res$p_adjust<0.01 & kegg.res$enrichmentScore > 1.5),]
up.kegg.res$group=1

lapply(1:nrow(down.kegg.res), function(i){
  gseaplot2(egmt,down.kegg.res$ID[i],title = down.kegg.res$Description[i],pvalue_table = T)
  ggsave(paste0(gsub("/","_",down.kegg.res$Description[i]),".down.pdf"),width = 11,height =5)
}
)

down.kegg.res$ID[1]
lapply(1:nrow(up.kegg.res), function(i){
  gseaplot2(egmt,up.kegg.res$ID[i],title = up.kegg.res$Description[i],pvalue_table = T)
  ggsave(paste0(gsub("/","-",up.kegg.res$Description[i]),".up.pdf"),width = 11,height =5)
}
)
```
### ②自定义通路
fgsea能够快速对预选基因集进行GSEA富集分析，预选基因集可以是自己设定，一般使用MSigdb数据库（同样由提出GSEA方法团队提供）。
```
library(msigdbr)
library(fgsea)
#提取预选基因集
m_df<- msigdbr(species = "Mus musculus",  category = "C2", subcategory = "KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#准备用于富集分析的基因数据
genelist=deg$avg_log2FC
names(genelist)= rownames(deg)
genelist=sort(genelist,decreasing = T)

#进行富集
fgseaRes<- fgsea(fgsea_sets, stats = genelist )

#进行可视化,筛选p值<0.005的20个通路，NES<1.5的视为下调基因，反之则是上调
ggplot(fgseaRes %>% filter(padj < 0.005) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 1.5)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +
  theme_minimal() ####以1.5进行绘图填色
```
![fgseaRes](https://upload-images.jianshu.io/upload_images/28382212-f254d1a665ab0489.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
![](https://upload-images.jianshu.io/upload_images/28382212-80c4974a1a633b4d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

### ③自定义画图
**加载数据，寻找差异基因**
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

setwd("C:/Users/86269/Desktop/shun.C/single_cell")
###加载所需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(plyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
load("scRNA1.rdata")

Idents(scRNA1)="celltype"

deg=FindMarkers(scRNA1,ident.1 = "Adipocytes",ident.2 = "Granulocytes",
                min.pct = 0.01,logfc.threshold = 0.01)

```
**进行富集分析**
```
#按照FoldChange排序
gsym.fc.id.sorted <- deg[order(deg$avg_log2FC,decreasing = T),]

#将基因ID从symbol转化为entrez
expm.id <- bitr(rownames(gsym.fc.id.sorted), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
exp.fc.id <- merge(expm.id, gsym.fc.id.sorted,by.x=1, by.y=0 )

#按FoldChange排序
gsym.fc.id.sorted <- exp.fc.id[order(exp.fc.id$avg_log2FC,decreasing = T),]
#获得ENTREZID、foldchange列表，做为GSEA的输入
id.fc <- gsym.fc.id.sorted$avg_log2FC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID

#进行gseKEGG富集分析
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
#这一条语句就做完了KEGG的GSEA分析
kk <- gseKEGG(id.fc, organism = "mmu")
#选择两条通路作图
gseaplot2(kk, geneSetID = c(1,2), 
          color = c("darkgreen","chocolate4")) #自定义颜色

```
![](https://upload-images.jianshu.io/upload_images/28382212-0f2e424a9ec396cc.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

**自定义作图**
作图前数据准备
```
#把ENTREZ ID转为gene symbol，便于查看通路里的基因
kk.gsym <- setReadable(kk, 'org.Mm.eg.db', #物种
                       'ENTREZID')
#按照enrichment score从高到低排序，便于查看富集的通路
sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T),]
# 要画的通路 
geneSetID <- c("mmu03010", "mmu00982")
# 突出显示感兴趣的基因
selectedGeneID <- c("Rps28","Rps29","Rps26","Rps2","Rpl22l1")



# 自定义足够多的颜色
mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
# 字号
base_size = 12
# rank mode
rankmode <- "comb" #合起来画
#rankmode <- "sep" #如果通路多，分开画更好

#合并多条通路的数据
gsdata <- do.call(rbind, lapply(geneSetID,enrichplot:::gsInfo, object = x))#对两条通路分别用gsInfo提取富集结果
gsdata$gsym <- rep(gsym.fc.id.sorted$SYMBOL,2)
```
分别画出GSEA富集图的三个部分
```
# 画running score
p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
  geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
  scale_color_manual(values = mycol) +
  
  scale_x_continuous(expand=c(0,0)) + #两侧不留空
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
  ylab("Enrichment\n Score") +
  
  theme_bw(base_size) +
  theme(panel.grid = element_blank()) + #不画网格
  
  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +
  
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
```
![running score](https://upload-images.jianshu.io/upload_images/28382212-f66859b436d480e3.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```
# 画rank
if (rankmode == "comb") {
  rel_heights <- c(1.5, .2, 1.5) # 整个GSEA图上中下三个部分的比例
  
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_bar(position = "dodge",stat = "identity",
             aes(x, y = position*(-0.1), fill=Description))+#图中所有竖线即y值相等，都是0
    xlab(NULL) + ylab(NULL) + 
    scale_fill_manual(values = mycol) + #用自定义的颜色
    
    theme_bw(base_size) +
    theme(panel.grid = element_blank()) + #不画网格
    
    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
} else if (rankmode == "sep") {
  rel_heights <- c(1.5, .5, 1.5) # 上中下三个部分的比例
  
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
  }
  
  head(gsdata)
  
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) + 
    scale_color_manual(values = mycol) + #用自定义的颜色
    
    theme_bw(base_size) +
    theme(panel.grid = element_blank()) + #不画网格
    
    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
} else {stop("Unsupport mode")}

```

![rank](https://upload-images.jianshu.io/upload_images/28382212-7da8cc099de4aab8.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```
# 画变化倍数
df2 <- p.res$data
df2$y <- p.res$data$geneList[df2$x]


# 提取感兴趣的基因的变化倍数
selectgenes <- data.frame(gsym = selectedGeneID)
selectgenes <- merge(selectgenes, df2, by = "gsym")
selectgenes <- selectgenes[selectgenes$position == 1,]

p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) + 
  geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), 
               color = "grey") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = mycol, guide=FALSE) + #用自定义的颜色
  scale_color_manual(values = mycol, guide=FALSE) + #用自定义的颜色
  
  scale_x_continuous(expand=c(0,0)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
  ylab("Ranked list\n metric") +
  xlab("Rank in ordered dataset") +
  
  theme_bw(base_size) +
  theme(panel.grid = element_blank()) +
  
  # 显示感兴趣的基因的基因名
  geom_text_repel(data = selectgenes, 
                  show.legend = FALSE, #不显示图例
                  direction = "x", #基因名横向排列在x轴方向
                  ylim = c(2, NA), #基因名画在-2下方
                  angle = 90, #基因名竖着写
                  size = 2.5, box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.3, "lines")) +
  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
```
![p.pos](https://upload-images.jianshu.io/upload_images/28382212-b98239420e53fdd7.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```
# 组图
plotlist <- list(p.res, p2, p.pos)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text())

plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)
```
![](https://upload-images.jianshu.io/upload_images/28382212-8a8efdf99f203aba.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
