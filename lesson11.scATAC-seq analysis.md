# 单个样本的scATAC-seq分析
ATAC-seq是用来检测染色质可及性的的测序方法，也是一种确定基因表达调控机制的方法。该方法可以帮助识别启动子区域、潜在的增强子或抑制子。
为了找到开放染色质区，基因组被TN5转座酶处理。在ATAC-Seq中，修饰后的TN5将与NextEra接头相对应的DNA序列插入到基因组的开放区域，同时，DNA被转座酶的活性剪切。
![image](https://user-images.githubusercontent.com/112565216/191719399-311aabbd-9c9f-41e4-954d-5dcf7d8aa478.png)

scATAC-seq的分析一般可用signac包进行

在Signac官网上下载样本数据后开始进行单个样本分析pipeline
## 1.加载需要的package
```
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(EnsDb.Mmusculus.v79)
setwd("/mnt/c/lesson/atac")
```

## 2.创建seurat对象进行后续分析
```
counts <- Read10X_h5("filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "singlecell.csv",
  header = TRUE,
  row.names = 1
)
brain_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = 'fragments.tsv.gz',
  min.cells = 1
)
brain <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)
```



```
#提取annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
#将annotation改为UCSC形式（因为参考基因组hg19来自UCSC）
seqlevelsStyle(annotations) <- 'UCSC'
#将annotation信息加入基因组
Annotation(brain) <- annotations
```

## 3.进行质控

在ATAC-seq中常使用以下指标来评估数据质量：

1.Nucleosome banding pattern（核小体条带模式）：绘制片段大小（fragment sizes）的直方图（由配对末端测序读数确定）来展示核小体的条带模式。我们计算每个单细胞，并量化单核小体与无核小体片段（存储为nucleosome_signal）的近似比率。

2.Transcriptional start site (TSS) enrichment score（转录起始位点富集得分）：ENCODE项目已根据以TSS为中心的片段与TSS侧翼区域中的片段的比率定义了ATAC-seq定位得分。较差的ATAC-seq实验通常将具有较低的TSS富集得分。我们可以使用TSSEnrichment函数为每个细胞计算该指标，并将结果存储在元数据中，列名称为TSS.enrichment。

3.Total number of fragments in peaks（peaks中片段的总数）：细胞测序深度/复杂性的量度。由于测序深度低，可能需要排除reads数很少的细胞，reads数极高的细胞可能代表了doublets或核团块。

4.Fraction of fragments in peaks（peaks中片段的比例）：表示落在ATAC-seq peaks内的总片段比例。比例较低的细胞（即<15-20％）通常代表应删除的低质量细胞或技术误差。

5.Ratio reads in ‘blacklist’ sites（比对到“blacklist”区域的reads比率）：ENCODE项目提供了一个“blacklist“”区域列表，这些区域表示通常与人为信号关联的读取。reads比对到这些区域的比例较高的细胞（与reads比对到peaks区域比例相比）通常代表了技术误差，应将其删除。Signac包中包含了针对人类（hg19和GRCh38），小鼠（mm10），果蝇（dm3）和秀丽隐杆线虫（ce10）的ENCODE“blacklist”区域。

```
#计算核小体比例
brain <- NucleosomeSignal(object = brain)
brain$nucleosome_group <- ifelse(brain$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = brain, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
```
![36fb6c3cf1f9005f5f0f0ab107947e2](https://user-images.githubusercontent.com/112565216/191746499-ef1bb91b-c92a-4cea-8ac7-44a7ae6c6a9b.png)

```
#计算TSS
brain <- TSSEnrichment(brain, fast = FALSE)
brain$high.tss <- ifelse(brain$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(brain, group.by = 'high.tss') + NoLegend()
```
![2a8618a6b62b8f958bc6c59e5ff84e5](https://user-images.githubusercontent.com/112565216/191746586-5d9574b6-70e8-4dab-a1df-578f408729b9.png)

```
#计算Fraction of fragments in peaks
brain$pct_reads_in_peaks <- brain$peak_region_fragments / brain$passed_filters * 100
#计算Ratio reads in ‘blacklist’ sites
brain$blacklist_ratio <- brain$blacklist_region_fragments / brain$peak_region_fragments
#将质控指标可视化
VlnPlot(
  object = brain,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
```
![cc3461ff8575502a0fc670bfa602b38](https://user-images.githubusercontent.com/112565216/191747295-d31bfcd4-f3d3-4a44-a021-bc6f2c0556af.png)

```
#筛选细胞
brain <- subset(
  x = brain,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

```
## 4.降维聚类
```
#归一化
brain <- RunTFIDF(brain)
#寻找高变基因
brain <- FindTopFeatures(brain, min.cutoff = 'q0')
#降维聚类
brain <- RunSVD(object = brain)
brain <- RunUMAP(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)
brain <- FindNeighbors(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)
brain <- FindClusters(
  object = brain,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)
#降维聚类结果可视化
DimPlot(object = brain, label = TRUE) + NoLegend()
```
![d7156dce16966c28ae199f9bd01a310](https://user-images.githubusercontent.com/112565216/191746832-e38a950e-4293-412a-90ea-2f6c4e6b8d54.png)

## 5.创建基因活性矩阵
```
gene.activities <- GeneActivity(brain)
brain[["RNA"]] <- CreateAssayObject(counts=gene.activities)
brain <- NormalizeData(
  brain,
  assay="RNA",
  normalization.method="LogNormalize",
  scale.factor=median(brain$nCount_RNA)
)
#将特定基因的表达量可视化
DefaultAssay(brain)<-"RNA"
FeaturePlot(
  brain,
  features=c("Sst","Pvalb","Gad2","Neurod6","Rorb","Syt6"),
  pt.size=0.1,
  max.cutoff = "q95",
  ncol=3
)
```
![e22a7cfc1c28f8d1e9eeedd1664166f](https://user-images.githubusercontent.com/112565216/191747850-c35ae1a1-f235-4bbb-868d-238bb6d5a969.png)

```
#将特定基因的不同区域DNA accessibility可视化
DefaultAssay(brain) <- "peaks"
CoveragePlot(
  brain,
  region=c("Gapdh"),
  extend.upstream = 1000,
  extend.downstream = 1000,
  ncol=1
)
```
![6228f110a66c9bf0ccb25de3be3577f](https://user-images.githubusercontent.com/112565216/191747236-fa415a07-de85-401c-a1f3-ba114fe5e1a9.png)

# 多个样本合并
当多个样品单独处理时，得到的peak一般不会完全一致，因此如果要合并多个样本，需要create a common peak set
合并多个peak可以用GenomicRanges::reduce或GenomicRanges::disjoin,二者的区别如下
```
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(20, 70, 300), end = c(120, 200, 400)))
```
