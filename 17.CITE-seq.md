CITE-seq技术是将特异性的寡核苷酸偶联到不同的抗体上，使其可以把对蛋白的测量转化为对抗体相连的DNA标签（ADTs）的测量。因此CITE-seq可以通过测序手段获得同意细胞中RNA和细胞表面蛋白的丰度（主要是膜蛋白）
```
library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)

InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
save(bm,file = "bm.rdata")

#对RNA基因降维聚类
DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

#对ADT基因降维聚类
DefaultAssay(bm) <- 'ADT'
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')
  
 #整合RNA与ADT
bm <- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)


bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
bm <- FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

# 可视化
#RNA与ADT的umap降维结果
p1 <- DimPlot(bm, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
```
![image](https://user-images.githubusercontent.com/112565216/195991208-07e70912-19be-4b22-a112-64b37b811d9f.png)

```
#RNA与ADT的pca降维结果
bm <- RunUMAP(bm, reduction = 'pca', dims = 1:30, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
bm <- RunUMAP(bm, reduction = 'apca', dims = 1:18, assay = 'ADT', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
p3 <- DimPlot(bm, reduction = 'rna.umap', group.by = 'celltype.l2', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(bm, reduction = 'adt.umap', group.by = 'celltype.l2', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4

```
![image](https://user-images.githubusercontent.com/112565216/195991240-bc8a9c72-1b7a-4d11-96d2-33921182afad.png)

```
#RNA与ADT部分基因
p5 <- FeaturePlot(bm, features = c("adt_CD45RA","adt_CD16","adt_CD161"),
                  reduction = 'wnn.umap', max.cutoff = 2, 
                  cols = c("lightgrey","darkgreen"), ncol = 3)
p6 <- FeaturePlot(bm, features = c("rna_TRDC","rna_MPO","rna_AVP"), 
                  reduction = 'wnn.umap', max.cutoff = 3, ncol = 3)
p5 / p6
```
![image](https://user-images.githubusercontent.com/112565216/195991365-0816ae94-175a-4230-82b5-ca3d1e564638.png)

**iCellR**
```
library("iCellR")

# Read RNA file
rna.data <- read.delim("CITE-Seq_sample_RNA.tsv.gz",header=TRUE)
# Read ADT file
adt.data <- read.delim("CITE-Seq_sample_ADT.tsv.gz",header=TRUE)
# make iCellR object
my.obj <- make.obj(rna.data)
my.obj <- add.adt(my.obj, adt.data = adt.data)

#为质控提供信息
my.obj <- qc.stats(my.obj,
                   s.phase.genes = s.phase, 
                   g2m.phase.genes = g2m.phase)


# filter 
my.obj <- cell.filter(my.obj,
                      min.mito = 0,
                      max.mito = 0.07 ,
                      min.genes = 500,
                      max.genes = 4000,
                      min.umis = 0,
                      max.umis = Inf)

# normalize RNA
my.obj <- norm.data(my.obj, norm.method = "ranked.glsf", top.rank = 500) 

# normalize ADT
my.obj <- norm.adt(my.obj)

#选择gene stats
my.obj <- gene.stats(my.obj, which.data = "main.data")

# find genes for PCA
my.obj <- make.gene.model(my.obj, my.out.put = "data",
                          dispersion.limit = 1.5, 
                          base.mean.rank = 500, 
                          no.mito.model = T, 
                          mark.mito = T, 
                          interactive = F,
                          no.cell.cycle = T,
                          out.name = "gene.model")

# merge RNA and ADT data
my.obj <- adt.rna.merge(my.obj, adt.data = "main")

# run PCA and the rest is as above
my.obj <- run.pca(my.obj, method = "gene.model", gene.list = my.obj@gene.model,data.type = "main")

# find the model genes to run a second round of PCA.
my.obj <- find.dim.genes(my.obj, dims = 1:20,top.pos = 20, top.neg = 20)
# second round PC
my.obj <- run.pca(my.obj, method = "gene.model", gene.list = my.obj@gene.model,data.type = "main")

#run umap
my.obj <- run.umap(my.obj, dims = 1:10)

#可视化部分ADT基因
A = gene.plot(my.obj, 
              gene = "ADT_CD3",
              plot.data.type = "umap",
              interactive = F,
              cell.transparency = 0.5)

B = gene.plot(my.obj, 
              gene = "CD3E",
              plot.data.type = "umap",
              interactive = F,
              cell.transparency = 0.5)

C = gene.plot(my.obj, 
              gene = "ADT_CD16",
              plot.data.type = "umap",
              interactive = F,
              cell.transparency = 0.5)

D = gene.plot(my.obj, 
              gene = "FCGR3A",
              plot.data.type = "umap",
              interactive = F,
              cell.transparency = 0.5)

library(gridExtra)
grid.arrange(A,B,C,D)

```
![image](https://user-images.githubusercontent.com/112565216/195993763-59897460-0890-41dc-89dc-b99319ad6bd5.png)


