当多个样本一同进行分析时，需要将多个样本进行整合，目的是让不同样本的细胞均匀地分布到各个cluster中。整合方法主要有锚点法和harmony法。
# 锚点法
>锚点法整合速度很慢，且常常过度整合，因此实际操作中跨物种整合或不同的数据类型如ATAC、蛋白组数据与单细胞数据整合时可用锚点整合，单纯的单细胞数据整合用harmony法
### 第一步 设置工作路径，加载需要的package
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

##改变工作路径，将相关结果储存至对应文件夹
setwd("C:/Users/86269/Desktop/shun.C/single_cell")
##加载需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(future)
```
### 第二步 批量读取数据并转化为Seurat对象
采用这种读取方式是因为锚点整合在归一化处理时需要每个样本单独进行
```
dir=c("BC2/","BC21/")
names(dir)=c("sample2","sample21")
#读取数据并批量转化
scRNA.list <- list()
for(i in 1:length(dir)){
  counts=Read10X(data.dir = dir[i])
  scRNA.list[[i]]<- CreateSeuratObject(counts,min.cells = 3,min.features = 2000)

}
#此处省略过滤线粒体与红细胞的质控过程
```
### 第三步 进行归一化处理并寻找高变基因
```
for(i in 1:length(dir)){
  scRNA.list[[i]]<- NormalizeData(scRNA.list[[i]])
  scRNA.list[[i]]<- FindVAriableFeatures(scRNA.list[[i]],selection.method='vst')

}
#找到不同细胞共有的高变基因
features <- SelectIntegrationFeatures(object.list=ifnb.list)
#找到锚定点
scRNA.anchors <- FindIntegrationAnchors(object.list=scRNA.list,anchor.features=features)
#进行数据整合
scRNA1 <- IntegrateData(anchorset = scRNA.anchors)
#接下来进行ScaleData、PCA、TSNE、UMAP等处理（略）
```

# harmony法
与锚点整合法不同，锚点整合法在对数据进行归一化处理时，需要将不同样本的数据分开处理，但harmony法不需要拆分样本（harmony算法不建议与SCTransform一起使用）
## 第一步 读取数据并创建Seurat对象
```
scRNA2 <- Read10X("BC2/")
scRNA21 <- Read10X("BC21")
scRNA2 <- CreateSeuratObject(scRNA2,project = "sample2", min.cells=3,min.features=2000)
scRNA21 <- CreateSeuratObject(scRNA21,project = "sample21",min.cells = 3,min.features = 2000)

scRNA_harmony <- merge(scRNA2,y=c(scRNA21))
```
## 第二步 对Seurat对象进行归一化处理与PCA降维
harmony是根据PCA的数据进行整合，因此在整合前必须要对数据进行归一化处理与PCA降维
```
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
```
## 第三步 进行harmony整合
```
#group.by.vars指定分组依据
scRNA_harmony <- RunHarmony(scRNA_harmony,group.by.vars = "orig.ident")
```
##第四步 按照正常流程进行umap降维、tsne降维与可视化与后续步骤
```
#需指定reduction参数为harmony
scRNA_harmony <- RunUMAP(scRNA_harmony,reduction = "harmony",dims=1:15)
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction = "harmony",dims=1:15) %>% FindClusters(resolution = 0.5)

plot1=DimPlot(scRNA_harmony,reduction = "umap",label=T)
plot1

