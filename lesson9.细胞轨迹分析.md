# RNA velocity
RNA velocity是指未经剪切的前体mRNA与已剪切的成熟mRNA含量的比例，来预测基因的动态变化。一般来说，如果未经剪切的前体mRNA较多，则该基因的表达量增加，反之则表达量降低。利用前体mRNA与成熟mRNA丰度的比值可以
## velocyto
### 1.velocyto生成loom文件
loom文件中包含了基因的成熟与未成熟的转录本的数量信息，需要利用cellranger结果中的filtered_feature_bc_matrix生成

![image](https://user-images.githubusercontent.com/112565216/189914551-7a87b561-182d-44c3-bc74-3c1cb83953d4.png)

但要区分成熟的mRNA与未成熟mRNA还需要利用bam文件，bam文件是fastq文件与参考基因组比对的结果，而10x matrix文件则是定量的结果

**生成loom文件**
```
velocyto run10x -m lessons/lesson2/hg38_rmsk.gtf lessons/lesson2/CellrangerResult refdata-gex-GRCh38-2020-A/genes/genes.gtf
```
run10x是针对cellranger结果进行分析，给出的三个路径分别是：重复区域注释文件（需要在UCSC→tools→table browser中下载），cellranger结果，参考基因组的索引
![image](https://user-images.githubusercontent.com/112565216/190146633-30422b3a-4831-471d-8587-37e9b04783ba.png)

运行结束后在储存cellranger结果的文件夹中生成一个loom文件

### 2.安装并打开rstudio-server
根据rstudio官网指示下载R与rstudio-server后用"sudo rstudio-server start"运行
如果需要更换不同版本的R，则在终端内输入：
```
sudo vi /etc/rstudio/rserver.conf
```
并将目标R版本所在的路径输入其中
```
rsession-which-r=/home/chen/anaconda3/envs/velocyto.r/bin/R
```

### 3.进行RNA velocity分析
**通过loom文件进行RNA velocity分析**
```
library(Seurat)
library(devtools)
library(velocyto.R)
library(SeuratWrappers)
setwd("/mnt/c/lesson")

#读取loom文件
ldat <- ReadVelocity(file="CellrangerResult.loom")
```
![1663417927(1)](https://user-images.githubusercontent.com/112565216/190856878-9d1eff84-4705-4c80-9cf7-8460aa00d9b5.png)


```
#转化为seurat对象
bm <- as.Seurat(x=ldat)
#对成熟的RNA归一化处理，降维聚类
bm <- SCTransform(object = bm,assay="spliced")
bm <- RunPCA(object = bm,verbose=FALSE)
bm <- FindNeighbors(object = bm,dims=1:20)
bm <- FindClusters(object = bm)
bm <- RunUMAP(object = bm,dims=1:20)

#挑选前1000个细胞加快分析速度
bm <- bm[,1:1000]

#进行RNA velocity分析
bm <-RunVelocity(object = bm,deltaT=1,kCells=25,fit.quantile=0.02)
#根据聚类得到的cluster数量选择颜色数目
ident.colors <- (scales::hue_pal())(n = length(levels(bm)))
#给每个cluster分配一种颜色
names(ident.colors) <- levels(x = bm)
#给每个细胞分配一种颜色（同色系渐变色）
cell.colors <- ident.colors[Idents(bm)]
names(cell.colors) <- colnames(bm)
#RNA velocity结果可视化
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)
```
![image](https://user-images.githubusercontent.com/112565216/190856163-6c615fef-e5ed-4ffd-ad10-a48259e6af8b.png)

但使用loom文件的聚类结果与使用seurat对象聚类得到的结果明显不同，为了将两者统一，接下来需要用seurat对象进行RNA velocity分析

**通过seurat对象进行RNA velocity分析**

准备已经降维聚类的seurat对象
```
#在cellranger结果的filter_feature_bc_matrix中含有可被10x读取的barcode,matrix和features等信息，利用其生成seurat对象
scRNA.counts <- Read10X("filter_feature_bc_matrix/")
scRNA = CreateSeuratObject(scRNA.counts ,min.cells = 3,project="os", min.features = 300)
###质控（线粒体基因与红细胞基因）
#计算细胞中线粒体基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)
#过滤
scRNA1 <- subset(scRNA, subset = nFeature_RNA > 300& nFeature_RNA < 7000 & percent.mt < 10 & percent.HB < 3 & nCount_RNA < 100000)
###降维聚类
scRNA1 <- NormalizeData(scRNA1, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst", nfeatures = 3000) 
scRNA1 <- ScaleData(scRNA1)
#查看细胞周期对聚类结果的影响
#细胞周期评分
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA1))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA1))
scRNA1 <- CellCycleScoring(object=scRNA1,g2m.features=g2m_genes,s.features=s_genes)
#查看细胞周期基因对细胞聚类的影响
scRNAa <- RunPCA(scRNA1, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
p
```
![image](https://user-images.githubusercontent.com/112565216/190856486-eadd7191-2458-4e67-9e9d-ee5274482a2e.png)

根据p可知细胞周期对细胞聚类的影响不大
```
###进行pca,umap与tsne降维聚类
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1))
#选择dims数
ElbowPlot(scRNA1, ndims=30, reduction="pca")
```
![image](https://user-images.githubusercontent.com/112565216/190856539-39efd134-f436-4cfe-a82b-8a8437189921.png)

```
#根据ElbowPlot选择大约20个pc数
pc.num=1:20
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num) 
scRNA1 <- FindClusters(scRNA1, resolution = 1.0)
scRNA1 = RunTSNE(scRNA1, dims = pc.num)
scRNA1 <- RunUMAP(scRNA1, dims = pc.num)
plot1 = DimPlot(scRNA1, reduction = "tsne") 
plot2 = DimPlot(scRNA1, reduction = "umap")
```
![image](https://user-images.githubusercontent.com/112565216/190856691-715d0fbe-7368-4a49-8bff-a29a0901cceb.png)
![image](https://user-images.githubusercontent.com/112565216/190856701-1c4d6281-4a84-42af-84f6-3d599d567833.png)

```
###随机抽取500个细胞进行RNA velocity分析（主要是为了缩短分析时间）
seurat.object=scRNA1
seurat.object=seurat.object[,sample(colnames(seurat.object),500)]  
#读取loom文件
ldat <- read.loom.matrices("CellrangerResult.loom")
emat <- ldat$spliced
nmat <- ldat$unspliced
```
```
#将loom文件列名与seurat对象列名统一
colnames(emat)[1:10]
colnames(seurat.object)[1:10]
```
![1663418355(1)](https://user-images.githubusercontent.com/112565216/190857242-bfa9f7e4-9d49-497f-8624-739eaf49da6a.png)

```
colnames(emat) <- paste(substring(colnames(emat),18,33),"-1",sep="")
colnames(nmat) <- paste(substring(colnames(nmat),18,33),"-1",sep="")
```

```
#提取随机抽取的500个细胞的loom结果
emat = emat[,colnames(seurat.object)]
nmat = nmat[,colnames(seurat.object)]
#通过umap降维的两个维度计算细胞距离
emb <- seurat.object@reductions$umap@cell.embeddings
cell.dist <- as.dist(1-armaCor(t(seurat.object@reductions$umap@cell.embeddings)))
#计算RNA velocity
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=2,
                                            kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=0.02,
                                            n.cores=12)
```

```
#可视化
library(ggplot2)
gg <- UMAPPlot(seurat.object)
gg
```
![image](https://user-images.githubusercontent.com/112565216/190858121-6866ae68-9e22-46ab-9ce6-1a393afe5d94.png)

```
#给每个细胞设定一种颜色
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb)

p1 <- show.velocity.on.embedding.cor(emb,rvel.cd,n=30,scale='sqrt',
                                     cell.colors=ac(colors,alpha=0.5),
                                     cex=0.8,arrow.scale=2,show.grid.flow=T,
                                     min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                     do.par=F,cell.border.alpha = 0.1,
                                     n.cores=24,main="Cell Velocity")

p1
```
![image](https://user-images.githubusercontent.com/112565216/190858163-bebc8bd6-3fa3-4b7a-a25f-b7e6585461d4.png)

# scVelo
在jupyter notebook中运行scVelo
```
#加载需要的模块
import scvelo as scv
import scanpy as sc
#设置图片规格
scv.settings.verbosity = 3  
scv.settings.presenter_view = True  
scv.set_figure_params('scvelo')
#读取loom文件
adata = scv.read("CellrangerResult.loom",cache=False)
#去除基因名重复
adata.var_names_make_unique
#查看成熟mRNA与前体mRNA比例
scv.pl.proportions(adata)
```
![image](https://user-images.githubusercontent.com/112565216/190894147-5f9e7e8c-4a22-4fac-820f-0c5931ca9720.png)

```
#过滤基因（如果该基因只被很少的细胞表达则被过滤）
scv.pp.filter_genes(adata,min_shared_counts=30)
#对每个细胞进行归一化处理
scv.pp.normalize_per_cell(adata)
#计算离散度,得到高变基因
scv.pp.filter_genes_dispersion(adata,n_top_genes=2000)
scv.pp.log1p(adata)
```
以上4步可被一个函数代替 **scv.pp.filter_and_normalize**

```
#PCA降维
scv.pp.moments(adata,n_pcs=30,n_neighbors=30)
#计算RNA velocity
scv.tl.velocity(adata)#构建velocity graph
scv.tl.velocity_graph(adata)
```
```
###进行umap降维
#构建无向有权图
scv.pp.neighbors(adata,n_neighbors=30,n_pcs=40)
#聚类
sc.tl.leiden(adata)
#降维
scv.tl.umap(adata)
#可视化
sc.pl.umap(adata,color="leiden")
```
![image](https://user-images.githubusercontent.com/112565216/190895115-cc66f9be-c3f9-4dc1-9268-777950db47b0.png)

```
#RNA velocity可视化
scv.pl.velocity_embedding_stream(adata,basis="umap",color=["leiden","initial_size_spliced"])
```
![image](https://user-images.githubusercontent.com/112565216/190895340-284e970d-9903-4c9c-8776-4b67bb96b365.png)

左图表示不同的cluster之间存在的分化顺序关系，右图则表示经过剪切的成熟mRNA含量在各个cluster之间的丰度情况

```
#识别分化过程中重要基因
scv.tl.rank_velocity_genes(adata,groupby="leiden",min_corr=0.3)
df=scv.DataFrame(adata.uns["rank_velocity_genes"]["names"])
df.head()
```
![image](https://user-images.githubusercontent.com/112565216/190896242-56e97c73-b46a-4dd8-af1c-35ae20a80789.png)

```
#对特定基因可视化
scv.pl.velocity(adata,["DOCK4"],dpi=120,color=["leiden"])
```
![image](https://user-images.githubusercontent.com/112565216/190896272-002156d2-a735-4dc9-9ed9-d9e3a2533ab4.png)

从图2和图3可以看出基因表达量和RNA velocity有较强的相关性

```
#连贯性分析(代表细胞间相互转换的概率)
scv.pl.velocity_graph(adata,threshold=1,color=["leiden"])
```
![image](https://user-images.githubusercontent.com/112565216/190896453-e49c4c18-11d4-4a5b-874a-d2eae239f285.png)

线条的多少代表了细胞间的连接度，连接度越高，细胞间发生转换的概率越高

```
#pseudotime分析
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata,color="velocity_pseudotime",cmap="gnuplot")
```
![image](https://user-images.githubusercontent.com/112565216/190897230-aad9d7b9-bca3-48ca-b2c3-55ab2f5f8385.png)

深色细胞向浅色细胞方向发育分化

```
#PAGA分析
adata.uns["neighbors"]["distances"]=adata.obsp["distances"]
adata.uns["neighbors"]["connectivities"]=adata.obsp["connectivities"]
scv.tl.paga(adata,groups="leiden")
scv.pl.paga(adata,basis="umap",size=50,alpha=1,
           min_edge_width=2,node_size_scale=1.5)
```
![image](https://user-images.githubusercontent.com/112565216/190897388-75f69701-10f6-45a5-847b-b1b2ce48c6ae.png)

箭头方向代表细胞分化方向

```
#查看单个基因的velocity结果
scv.pl.velocity_graph(adata,color="GNLY")
```
![image](https://user-images.githubusercontent.com/112565216/190897772-ed6fc264-7647-4612-849d-40f0680f8076.png)

# PHATE
PHATE是Seurat pipeline中使用的tsne 与umap降维以外的另一种降维可视化方式，PHATE的优势在于可兼顾全局特征（cluster之间的联系）与局部特征（单个cluster内部的关系），且高效利用计算资源，可计算大批量数据，在可视化结果中平滑过渡代表着细胞从一种状态到另一种状态的过渡。PHATE常用语研究物种进化与细胞分化。

## 利用R进行PHATE降维
```
#安装需要的package
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
library(reticulate)
reticulate::py_install("phate", pip=TRUE)
devtools::install_github("KrishnaswamyLab/phateR")
if (!require(viridis)) install.packages("viridis")
if (!require(readr)) install.packages("readr")
if (!require(Rmagic)) install.packages("Rmagic")
library(phateR)
library(ggplot2)
library(readr)
library(viridis)
library(Rmagic)

#读取数据
setwd("C:/Users/86269/Desktop/shun.C/single_cell")
bmmsc <- read_csv("BMMC_myeloid.csv")
bmmsc <- bmmsc[,2:ncol(bmmsc)]
```
![image](https://user-images.githubusercontent.com/112565216/190995663-fa8bcc18-006b-4aeb-8923-503af7dce602.png)

```
#过滤低表达基因
keep_cols <- colSums(bmmsc > 0) > 10
bmmsc <- bmmsc[,keep_cols]
#过滤基因表达量过少的细胞
keep_rows <- rowSums(bmmsc) > 1000
bmmsc <- bmmsc[keep_rows,]
#消除测序深度影响
bmmsc <- library.size.normalize(bmmsc)
bmmsc <- sqrt(bmmsc)
```

```
#进行PCA降维
m <-prcomp(bmmsc)
bmmsc_PCA <- as.data.frame(m$x)
```
![image](https://user-images.githubusercontent.com/112565216/191009155-2549f69c-fd9e-46db-b491-b1024311c674.png)
![image](https://user-images.githubusercontent.com/112565216/191009214-9f279e3d-a522-4ef6-9d3a-2743bef5c504.png)

```
#PHATE可视化(PHATE是直接运算，而不是基于PCA结果运算)
bmmsc_PHATE <- phate(bmmsc, knn=4, decay=100, t=10, init=bmmsc_PHATE)
ggplot(bmmsc_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color=bmmsc$Mpo)) +
  labs(color="Mpo") +
  scale_color_viridis(option="B")
```
![image](https://user-images.githubusercontent.com/112565216/191016201-36ad71cc-0b1a-4a21-b164-6f569befecf0.png)

```
#利用MAGIC还原基因表达
bmmsc_MAGIC <- magic(bmmsc, t=4)
ggplot(bmmsc_PHATE) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=bmmsc_MAGIC$result$Ifitm1)) +
  scale_color_viridis(option="B") +
  labs(color="Ifitm1")
```
![image](https://user-images.githubusercontent.com/112565216/191020531-cd051ded-4710-4fd6-9a21-0c9a2cf882af.png)

## 用python进行PHATE分析
```
#加载需要的module
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import phate
import scprep
import magic
import os
sparse=True
#读取数据
T1=scprep.io.load_10X("scRNAseq/scRNAseq/T0_1A")
T2=scprep.io.load_10X("scRNAseq/scRNAseq/T2_3B")
T3=scprep.io.load_10X("scRNAseq/scRNAseq/T4_5C")
T4=scprep.io.load_10X("scRNAseq/scRNAseq/T6_7D")
T5=scprep.io.load_10X("scRNAseq/scRNAseq/T8_9E")
```
```
T1.head()
```
![image](https://user-images.githubusercontent.com/112565216/191179458-d8860525-8961-41f2-9b05-41e569a605f6.png)

```
#过滤library size过大或过小的细胞
filtered_batches=[]
for batch in [T1,T2,T3,T4,T5]:
    batch =scprep.filter.filter_library_size(batch,percentile=20,keep_cells="above")
    batch=scprep.filter.filter_library_size(batch,percentile=75,keep_cells="below")
    filtered_batches.append(batch)
#合并
EBT_counts, sample_labels=scprep.utils.combine_batches(
    filtered_batches,
    ["Day 00-03","Day 06-09","Day 12-15","Day 18-21","Day 24-27"],
    append_to_cell_names=True
)
EBT_counts.head()
```
![image](https://user-images.githubusercontent.com/112565216/191179634-bf12e392-8355-4793-89a6-db29ee70cb75.png)

```
#过滤表达量过少的基因
EBT_counts=scprep.filter.filter_rare_genes(EBT_counts,min_cells=10)
#归一化处理
EBT_counts=scprep.normalize.library_size_normalize(EBT_counts)
#过滤线粒体含量过高的细胞
mito_genes=scprep.select.get_gene_set(EBT_counts,starts_with="MT-")
scprep.plot.plot_gene_set_expression(EBT_counts,genes=mito_genes,percentile=90)
plt.show()
EBT_counts, sample_labels = scprep.filter.filter_gene_set_expression(
    EBT_counts, sample_labels, genes=mito_genes, 
    percentile=90, keep_cells='below')
```
![image](https://user-images.githubusercontent.com/112565216/191182104-75290ce6-f18c-4a5f-9964-69ca727e7a22.png)

```
#进行log transformation
EBT_counts=scprep.transform.sqrt(EBT_counts)
```

```
#进行phate可视化
phate_operator=phate.PHATE(n_jobs=2)
Y_phate=phate_operator.fit_transform(EBT_counts)
scprep.plot.scatter2d(Y_phate, c=sample_labels, figsize=(12,8), cmap="Spectral",
                      ticks=False, label_prefix="PHATE")
plt.show()
```
![image](https://user-images.githubusercontent.com/112565216/191183853-2e3ec951-643a-420f-a1c8-9d843553f520.png)

```
#进行3d的phate可视化
phate_operator.set_params(n_components=3)
Y_phate_3d=phate_operator.transform()
scprep.plot.scatter3d(Y_phate_3d,c=sample_labels,cmap="Spectral",
                      ticks=False,label_prefix="PHATE")
```
![image](https://user-images.githubusercontent.com/112565216/191184375-9e0f8b49-dba0-4172-9e1b-887cdc31586f.png)

```
#利用MAGIC算法还原数据缺失值
magic_op=magic.MAGIC()
chronic_magic=magic_op.fit_transform(EBT_counts,genes="all_genes")
```
![image](https://user-images.githubusercontent.com/112565216/191189612-067ffdba-a0cc-4bfd-88d0-2d232f64db46.png)


```
#对列名进行修改
chronic_magic.columns=[i.split(" ")[0] for i in chronic_magic.columns.tolist()]
```
![image](https://user-images.githubusercontent.com/112565216/191193086-5d6178d0-9154-42a6-848c-d8a56db79c04.png)


```
scprep.plot.scatter_2d(Y_phate,c=chronic_magic["A2M"],cmap="Reds",
                      ticks=False,label_prefix="PHATE",title="A2M"+"magic expression")
```
![image](https://user-images.githubusercontent.com/112565216/191193158-36c8b090-4151-4719-81d5-34c7fda1cca5.png)


# Cytotrace
```
library(reticulate)
conda_create("cytoTRACE",python_version = '3.7')
use_condaenv("cytoTRACE")
conda_install("cytoTRACE", "numpy") 
conda_install("cytoTRACE", "scanoconcoramaCT") 
library(CytoTRACE)
```
利用自带的数据库进行分析
```
results <- CytoTRACE(marrow_10x_expr) 
#可视化
plotCytoTRACE(results, phenotype = marrow_10x_pheno)
#将某个基因表达情况可视化
plotCytoTRACE(results, phenotype = marrow_10x_pheno,
              gene = "Gapdh" )

```
![image](https://user-images.githubusercontent.com/112565216/195977809-e8d279c7-5f82-4ac0-b3aa-0757f54dc351.png)

