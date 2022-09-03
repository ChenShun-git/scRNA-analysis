在完成对样本的整合、聚类与可视化后，样本细胞形成了许多cluster，需要得出这些cluster代表的细胞种类并将其注释到图上。

注释方法分为手动注释与软件注释两大类

# 手动注释
## 1. 寻找各cluster特异的高表达基因（marker）
```
markers <- FindAllMarkers(object=scRNA_harmony,test.use = "wilcox",only.pos = TRUE,
                          logfc.threshold = 0.25)
#筛选每个cluster特异表达基因表达量最高的前十个
all.markers = markers %>% select(gene,everything()) %>% subset(p_val<0.05)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
```
再查阅文献资料找出top10 的marker基因对应的细胞种类

## 2.寻找conserved marker
单细胞测序分析时经常会出现不同分组，例如药物处理与没有药物处理的上皮细胞基因表达也许存在差异，而这些差异基因不能作为marker，FindConservedMarker可以解决这个问题
```
#寻找两个sample中共有的marker基因
#FindConservedMarke()每次只能处理一个cluster，下面以cluster1 为例
marker2 <- FindConservedMarkers(scRNA_harmony,ident.1 = 1,grouping.var = "orig.ident",
                                only.pos = TRUE,
                                min.diff.pct = 0.25,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
```
## 3.利用气泡图
```
#例如T细胞的marker基因IL7R，CD3D，NKG7
DotPlot(scRNA_harmony,features=c("IL7R","CD3D","NKG7"))
```
![Dotplot](https://user-images.githubusercontent.com/112565216/188261795-a916c861-999b-4347-97ff-10e503140f74.png)


可以根据气泡图看出，cluster2的T细胞marker基因表达量较高，所以可以被归类为T细胞

## 4.重命名
用上述方法得出细胞种类后对cluster进行重命名
```
#以cluster7为例
scRNA_harmony <- RenameIdents(scRNA_harmony,"7"="Osteoblastic")
```



# 软件注释
##1.SingleR package
用已知的数据注释未知数据
### 第一步 加载需要的package并load前文已经处理好的scRNA.harmony数据
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(harmony)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(org.Hs.eg.db)
library(SingleR)
library(scRNAseq)
load("../scRNA_harmony")
```
### 第二步 准备参考数据与需要被注释的数据
利用下载的已知的人类细胞图谱的数据注释未知的scRNA_harmony的cluster
```
#下载参考数据
refdata = HumanPrimaryCellAtlasData()

#提取出scRNA的转录表达数据
testdata <- GetAssayData(scRNA_harmony,slot="data")
#提取每个细胞的cluster信息
clusters <- scRNA_harmony@meta.data$seurat_clusters

#开始使用SingleR进行分析
cellpred <- SingleR(test=testdata,ref=refdata,labels=refdata$label.main,
                    method="cluster",cluster=clusters,
                   assay.type.test = "logcounts",assay.type.ref = "logcounts")
#制作细胞类型的注释文件
celltype = data.frame(ClusterID = rownames(cellpred),celltype=cellpred$labels,stringsAsFactors = FALSE)
```
### 第三步 将注释结果加入metadata中
#### 方法一
```
scRNA_harmony@meta.data$celltype="NA"
for(i in 1:nrow(celltype)){
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters==celltype$ClusterID[i]),'celltype']
  <- celltype$celltype[i]
}
```
- 先在metadata中加入一列名为celltype的新列，本列所有数据都设置为NA
- 利用for循环，将metadata中每个cluster所含的细胞所在的行找出，并将该细胞对应的celltype列赋上SingleR的标记结果

#### 方法二
```
scRNA_harmony@meta.data$SingleR = celltype[match(clusters,celltype$ClusterID),'celltype']
```
- 在metadata中加入名为SingleR的新列
- match(clusters,celltype$ClusterID)代表celltype中各个clusters所在的行数
- 再将对应行数的celltype赋给SingleR列

## 2 Garnett方法
SingleR package注释法中参考数据来源与package自带的数据，因此有较大的局限性，如果要自定义参考数据则需要采用Garnett法。
Garnett法是利用基于机器学习训练分类器，利用训练好的分类器对细胞类型进行分类。其本质类似于气泡图注释法，利用分类器对气泡图进行更准确的判断。
###第一步 制作marker file（txt文件）
格式如下
**> B cells
  expressed:CD19,MS4A1,CD79A,ACTN,ACTB
  not expressed:
  subtype of:**
- 第一行是">"后接细胞种类
- 第二行是选定的marker基因
- 第三行是选定的阴性marker
- 第四行是细胞亚型
- 前两行是必填，但是后两行可以不填
本次示例采用下载的hsPBMC marker file
### 第二步 创建CDS对象 优化marker file
优化的目的是使marker file更适合所分析的单细胞数据
```
#安装garnett之前要先安装monocle3
devtools::install_github("cole-trapnell-lab/monocle3")
library(monocle3)
devtools::install_github('cole-trapnell-lab/garnett',ref='monocle3')
library(garnett)
#加载scRNA_harmony数据
load("scRNA_harmony.rdata")
pbmc <- scRNA_harmony
##创建CDS对象
data <- GetAssayData(pbmc,assay='RNA',slot="counts")
cell_metadata <- pbmc@meta.data
#gene_annotation实际上没有任何实际意义，仅仅是为了满足输入要求
gene_annotation <- data.frame(gene_short_name=rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata=gene_annotation)
#对CDS对象进行归一化与降维处理，preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds,num_dim=30)
#优化marker file
marker_check <- check_markers(cds,"hsPBMC_markers.txt",
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")
plot_markers(marker_check)
```
得到的结果如图所示
![示例.png](https://upload-images.jianshu.io/upload_images/28382212-1d38ae927ab832c9.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

根据图片，不在database中的marker和high overlap（即特异性不高）的marker都需要删除替换
### 第三步 得到分类器
分类器一般有两种获得途径，可以上Garne官网下载已经训练好的分类器（资源较少，一般找不到所需分类器），另一种方法就是自行训练
```
#使用marker file和cds对象训练分类器
pbmc_classifier <-train_cell_classifier(cds=cds,
                                        marker_file = "hsPBMS_markers.txt",
                                        db=org.Hs.eg.db,
                                        cds_gene_id_type = "SYMBOL",
                                        num_unknown = 50,
                                        marker_file_gene_id_type = "SYMBOL")
saveRDS(pbmc_classifier,"my_classifier.rds")
```
### 第四步 使用训练器进行注释
```
#读取marker file
hsPBMC <- readRDS("hsPBMC.rds")
pData(cds)$garnett_cluster <- pData(cds)$seurat_clusters
cds <- classify_cells(cds,
                      hsPBMC,
                      db = org.Hs.eg.db,
                      cluster_extend = TRUE,
                      cds_gene_id_type = "SYMBOL")
#提取分类结果
cds.meta <- subset(pData(cds),select=c('cell_type','cluster_ext_type')) %>% as.data.frame()
#将结果返回给seurat对象
pbmc <- AddMetaData(pbmc,metadata = cds.meta)
```

## 3 nnls（非负最小二乘回归）法
计算需要被注释的数据与已知的参考数据中哪些类型相关性更大
![公式.png](https://upload-images.jianshu.io/upload_images/28382212-cadb3b5fe1b06804.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
其中Ta表示需要被注释的数据集A中的细胞的基因表达，Mb表示参考数据集中的细胞基因表达

```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
```
### 读取数据并进行归一化降维处理
```
dir=c("BC2/","BC21/")
names(dir)=c("sample2","sample21")
#读取数据并进行归一化降维处理
scRNA.list <- list()
for(i in 1:length(dir)){
  counts=Read10X(data.dir = dir[i])
  scRNA.list[[i]]<- CreateSeuratObject(counts,min.cells = 3,min.features = 300)
  
}

seu.obj <- scRNA.list[[1]]
seu.obj <- seu.obj %>% NormalizeData(verbose=FALSE) %>% FindVariableFeatures(selection.method='vst') %>% ScaleData(verbose=FALSE)%>% RunPCA(pc.genes=seu.obj@var.genes,npcs=100,verbose=FALSE)%>% FindNeighbors(dims=1:10) %>% FindClusters(resolution=0.5) %>% RunUMAP(dims=1:10)

seu.obj2 = scRNA.list[[2]]
seu.obj2 <- seu.obj2 %>% NormalizeData(verbose=FALSE) %>% FindVariableFeatures(selection.method='vst')%>%ScaleData(verbose=FALSE) %>% RunPCA(pc.genes=seu.obj2@var.genes,npcs=100,verbose=FALSE) %>%FindNeighbors(dims=1:10) %>% FindClusters(resolution=0.5) %>% RunUMAP(dims=1:10)
```

```
#找到两个数据集共有的基因
shared_gene <- intersect(rownames(seu.obj),row.names(seu.obj2))

#计算所有基因在每个cluster中的表达量之和
seu.obj.mat <- AggregateExpression(seu.obj,assays="RNA",features=shared_gene,slot="count")$RNA
#除以测序深度进行normalization
seu.obj.mat <- seu.obj.mat / rep(colSums(seu.obj.mat),each=nrow(seu.obj.mat))
seu.obj.mat <- log10(seu.obj.mat * 100000+1)

seu.obj.mat2 <- AggregateExpression(seu.obj2,assays="RNA",features=shared_gene,slot="count")$RNA
seu.obj.mat2 <- seu.obj.mat2 / rep(colSums(seu.obj.mat2),each=nrow(seu.obj.mat2))
seu.obj.mat2 <- log10(seu.obj.mat2 * 100000+1)
```
### 筛选基因
根据目标数据集（A数据集中选定的cluster）中给定的细胞类型与整体细胞类型（数据集A）中位数的倍数变化，选择前200个，得到list1；然后根据目标数据集中给定的细胞类型与其他细胞类型（去除选定的cluster的数据集A）最大值的倍数变化，选择前200个，得到list2，然后将两个list进行合并
```
#以cluster3为为例
cluster <- 3
seu.obj.gene <- seu.obj.mat[,cluster+1]
#得出list1
gene_fc <- seu.obj.gene / apply(seu.obj.mat,1,median)
gene_list1 <- names(sort(gene_fc,decreasing=TRUE)[1:200])
#得出list2
gene_fc <- seu.obj.gene/apply(seu.obj.mat[,-(cluster+1)],1,max)
gene_list2 <- names(sort(gene_fc,decreasing=TRUE)[1:200])
#合并两个list
gene_list <- unique(c(gene_list1,gene_list2))
```
### 进行相关性计算
```
#A数据集中的某个细胞类型
Ta <- seu.obj.mat[gene_list,cluster+1]
#B数据集中的所有细胞类型
Mb <- seu.obj.mat2[gene_list,]
library(lsei)
solv <- nnls(Mb,Ta)
corr <- solv$x
corr
```

![corr](https://user-images.githubusercontent.com/112565216/188267215-b2f47fe1-7eb4-4a3f-9bed-2c4299b668c8.png)

运行结果中B数据集里与cluster3相关性最强的cluster（根据运行结果可知是cluster6）可以视作与cluster同一种类的细胞

### 进行验证
由结果可知，A数据集中的cluster2与B数据集中的cluster6相关性最高，接下来尝试找出B数据集中的cluster6与A数据集中的哪个cluster相关性最高（即将上述过程反向操作）
```
cluster=6
seu.obj.gene <- seu.obj.mat2[,cluster+1]

gene_fc <- seu.obj.gene / apply(seu.obj.mat2,1,median)
gene_list1 <- names(sort(gene_fc,decreasing=TRUE)[1:200])

gene_fc <- seu.obj.gene/apply(seu.obj.mat2[,-(cluster+1)],1,max)
gene_list2 <- names(sort(gene_fc,decreasing=TRUE)[1:200])

gene_list <- unique(c(gene_list1,gene_list2))   

Ta <- seu.obj.mat2[gene_list,cluster+1]
Mb <- seu.obj.mat[gene_list,]
solv <- nnls(Mb,Ta)
corr2 <- solv$x
corr2
```
![image](https://user-images.githubusercontent.com/112565216/188269455-7ffc8b57-1a54-475e-9d44-562fce6fd8c7.png)

根据结果可知cluster3 与cluster 6相关性最大
### 批量计算相关系数
上述是计算单个cluster相关系数的pipeline
接下来对数据集中多个cluster进行批量计算
```
library(lsei)
list1 <- list()
for (cluster in seq(1,ncol(seu.obj.mat))) {
  seu.obj.gene <- seu.obj.mat[,cluster]
  
  gene_fc <- seu.obj.gene / apply(seu.obj.mat,1,median)
  gene_list1 <- names(sort(gene_fc,decreasing=TRUE)[1:200])
  
  gene_fc <- seu.obj.gene/apply(seu.obj.mat[,-cluster],1,max)
  gene_list2 <- names(sort(gene_fc,decreasing=TRUE)[1:200])

  gene_list <- unique(c(gene_list1,gene_list2))   

  Ta <- seu.obj.mat[gene_list,cluster]
  Mb <- seu.obj.mat2[gene_list,]
  solv <- nnls(Mb,Ta)
  corr <- solv$x
  corr
  list1[[cluster]] <- corr
}
#用a预测b
list2 <- list()
for (cluster in seq(1,ncol(seu.obj.mat2))) {
  seu.obj.gene <- seu.obj.mat2[,cluster]
  
  gene_fc <- seu.obj.gene / apply(seu.obj.mat2,1,median)
  gene_list1 <- names(sort(gene_fc,decreasing=TRUE)[1:200])
  gene_fc <- seu.obj.gene/apply(seu.obj.mat2[,-cluster],1,max)
  gene_list2 <- names(sort(gene_fc,decreasing=TRUE)[1:200])

  gene_list <- unique(c(gene_list1,gene_list2))   
  Ta <- seu.obj.mat2[gene_list,cluster]
  Mb <- seu.obj.mat[gene_list,]
  solv <- nnls(Mb,Ta)
  corr <- solv$x
  corr
  list2[[cluster]] <- corr
} 
#合并结果
mat1 <- do.call(rbind,list1)
mat2 <- do.call(rbind,list2)

#计算系数
beta <- 2 * (mat1+0.01) * t(mat2+0.01)
row.names(beta) <- paste0("C",1:nrow(beta)-1)
colnames(beta)<- paste0("C",1:ncol(beta)-1)
```
### 进行可视化
### 根据热图来优化聚类结果
```
library(psych)
Idents(seu.obj)="seurat_clusters"
#计算表达矩阵每个cluster的平均值
exp = AverageExpression(seu.obj)
#画出相关性热图
corrda <- corr.test(exp$RNA,exp$RNA,method="spearman")
pheatmap::pheatmap(corrda$r)
#画出聚类图
DimPlot(seu.obj,label=T)
```
根据相关性热图可以判断聚类图上距离较近的cluster是否可以归于同一个cluster
示例：
![聚类图](https://user-images.githubusercontent.com/112565216/188270368-102663ce-3551-4b49-afc7-2e7391177683.png)
![热图](https://user-images.githubusercontent.com/112565216/188270322-5313ffd3-6766-488d-b52b-c78572a9774e.png)

根据热图聚类情况可知，cluster0,6,9的相关性很强，在聚类图中距离也很相近，可以被归为一个cluster
