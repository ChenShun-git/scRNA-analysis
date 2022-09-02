## 第一步 数据获取与前期准备
从NCBI上下载GEO数据后解压，在Rstudio中加载需要的package
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

##安装Seurat包
install.packages('Seurat')
##改变工作路径，将相关结果储存至对应文件夹
setwd("C:/Users/86269/Desktop/shun.C/single_cell")
##加载需要的包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
```
# 第二步 读取数据并进行质控
## 读取数据
```
##读取10×数据，文件夹中需要包括features.tsv,barcodes.tsv,matrix.mtx,分别代表基因名，标记序列，表达量矩阵
scRNA_counts = Read10X("C:/Users/86269/Desktop/shun.C/single_cell/BC21")
class(scRNA_counts)

##创建Seurat对象
##min.cells表示基因至少要在多少个细胞中被检测到表达才算有效
##min.features至少要检测到多少个基因才算有效
scRNA = CreateSeuratObject(scRNA_counts,project = "sample_21",
                           min.cells = 3, min.features = 300)
```
## 进行质控
过滤线粒体DNA含量过高的基因，因为线粒体基因含量过高表示细胞即将凋亡或状态不佳，过滤标准需要根据实际情况考虑，例如肝脏细胞与肌肉细胞本身含线粒体量较多，过滤阈值也应当适当增大。
过滤红细胞，将红细胞表达的标志性基因与样本基因匹配，如果一个细胞的这些基因表达量过高，就将该细胞视作红细胞。
```
#计算线粒体基因比例
scRNA[["percent.mt"]]=PercentageFeatureSet(scRNA,pattern = "^MT-")

#计算红细胞比例
HB.genes = c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes,rownames(scRNA@assays$RNA))
HB.genes = rownames(scRNA@assays$RNA)[HB_m]
HB.genes = HB.genes[!is.na(HB.genes)]
scRNA[["percent.HB"]] <- PercentageFeatureSet(scRNA,features = HB.genes)

#过滤低质量细胞
scRNA1 <- subset(scRNA,subset=nFeature_RNA>500 &
                   nCount_RNA>1000 & percent.mt<20 & percent.HB<1)
```
# 第三步 归一化，降维与聚类
1.归一化，排除测序深度的影响
2.降维与聚类
- 寻找高变基因，即在各细胞中表达量差异最大的基因，一般挑选3000个，这一步将维度降到3000
- 中心化处理，PCA降维前的必要步骤
- 对高变基因通过PCA降维，一般降至50个维度
- 制无向有权图并切图，对细胞进行聚类
- 进行降维，降至二维进行展示，有tsne与umap法
```
#归一化处理
scRNA1 <- NormalizeData(scRNA1,normalization.method="LogNormalize",scale.factor=10000)
#挑选高变基因
scRNA1 <- FindVariableFeatures(scRNA1,selection.method="vst",nfeatures=3000)
#对基因进行中心化处理，将每个细胞中高变基因的表达量改变，使细胞的平均表达为0，细胞间差异为1，成为标准数据，赋予每个基因相同的权重，因此高表达基因不占优势
scale_gene <- VariableFeatures(scRNA1)
scRNA1 <- ScaleData(scRNA1,features=scale_gene)
#将高变基因通过PCA降维,RunPCA默认选取50个PC
scRNA1 <- RunPCA(scRNA1,features=VariableFeatures(scRNA1))
##聚类处理
pc.num=1:20
scRNA1 <- FindNeighbors(scRNA1,dims=pc.num)
scRNA1 <- FindClusters(scRNA1,resolution=1.0)
##可视化聚类，将高维数据降低到二维or三维进行展示
scRNA1 <- RunTSNE(scRNA1,dims=pc.num)
embed_tsne <- Embeddings(scRNA1,"tsne")
#tsne法
DimPlot(scRNA1,reduction="tsne",label=TRUE)
#umap法
DimPlot(scRNA1,reduction="umap",label=TRUE)
```
# 第四步 细胞周期评分
scRNA测序中，单个样品里有许多细胞，每个细胞的表达的基因种类与表达量都不相同，根据基因表达的不同将细胞进行分类，但细胞基因表达的差距除了是因为细胞种类的差别，也有可能是因为细胞所处周期的差别，细胞聚类可能会被细胞周期影响，所以需要检查细胞周期的影响，如果细胞周期影响较大则需要排除细胞周期因素
```
#细胞周期评分
#cc.genes是Seurat包自带的一个list，含不同细胞周期的marker基因
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA1))
g2m_genes=cc.genes$g2m.genes
g2m_genes=CaseMatch(search = g2m_genes,match=rownames(scRNA1))
s_genes=cc.genes$s.genes
s_genes=CaseMatch(search=s_genes,match=rownames(scRNA1))
scRNA1 <- CellCycleScoring(object=scRNA1,g2m.features = g2m_genes,s.features = s_genes)
#查看细胞周期基因对细胞聚类的影响
scRNAa <- RunPCA(scRNA1,features=c(s_genes,g2m_genes))
p <- DimPlot(scRNAa,reduction = "pca",group.by="Phase")
p
```
```
#如果细胞周期有影响则消除细胞周期的影响
#scRNAb <- ScaleData(scRNA1,vars.to.regress="S.Score","G2M.Score")
```
细胞周期影响的判定需要在中心化处理后，排除细胞周期影响后就可以继续进行PCA降维等后续步骤
