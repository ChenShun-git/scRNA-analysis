# SCTransform
SCTransform一个函数包括NormalizeData，ScaleData,FindVariable三个函数的功能，因此使用SCTransform的pipeline如下
```
#读取数据，计算线粒体基因占比
scRNA_counts = Read10X("C:/Users/86269/Desktop/shun.C/single_cell/BC21")
scRNA = CreateSeuratObject(scRNA_counts,project = "sample_21",
                           min.cells = 3, min.features = 300)
                           min.cells = 3, min.features = 300)
scRNA = PercentageFeatureSet(scRNA,pattern="^MT-",col.name="percent.mt
#使用SCTransform进行归一化
scRNA =SCTransform(scRNA,vars.to.regress = "percent.mt",verbose = FALSE)

#基于SCTransform继续进行降维聚类
#PCA降维
scRNA <- RunPCA(scRNA,verbose = FALSE)

#umap聚类
scRNA <- RunUMAP(scRNA,dims=1:30,verbose=FALSE)

#tsne聚类
scRNA <- FindNeighbors(scRNA,dims=1:30,verbose=FALSE)
scRNA <- FindClusters(scRNA,verbose=FALSE)
#可视化
DimPlot(scRNA,label=TRUE)
```


# Seurat对象的结构
![Seurat对象.png](https://upload-images.jianshu.io/upload_images/28382212-41f01e1d7e8b2af3.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
Seurat对象最重要的成分是assays与meta.data
- assays主要保存了所有细胞的UMI矩阵
- meta.data
![meta.data.png](https://upl
oad-images.jianshu.io/upload_images/28382212-d431117716e36e9d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
- orig.ident：储存细胞的样本来源，也可以自行修改
- ncount_RNA：每个细胞的UMI数量
- nFeature_RNA：每个细胞检测到的基因数

# UMI与barcodes
![barcode与UMI.png](https://upload-images.jianshu.io/upload_images/28382212-55cfc0d7ad22c0a0.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
图中绿色片段是barcode，右边紧接着就是UMI。每一种barcode代表一个细胞，每一种UMI代表一个count

# marker基因的寻找
1.FindMarkers：用于寻找两个不同组别的细胞的marker基因
2.FindAllMarkers：用于寻找一个细胞集中每个cluster的marker基因
3.FindConservedMarkers：用于寻找在两个不同组别中的保守的marker基因（即这两个组别中表达相似，但与别的组别有明显差异的marker基因）

# 数据标准化
数据标准化是指将数据按照比例放缩，使其落入特定区间
数据归一化时最简单的数据标准化方法，其目的是将数据变为（0，1）之间的小数并将有量纲的数据转化为无量纲的纯数量
