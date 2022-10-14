从稀疏矩阵开始分析
```
library(data.table)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(tidyHeatmap)
library(Seurat)
library(ggthemes)
library(plyr)
library(tidyverse)
setwd("C:/Users/86269/Desktop/shun.C/single_cell/spatial")

dataA <- fread("GSM3036911_PDAC-A-ST1-filtered.txt",data.table = F, header = T,check.names = F)
```
<img width="416" alt="image" src="https://user-images.githubusercontent.com/112565216/195783841-2d31d835-6877-40a1-bc8e-cbc1dceb0ff1.png">

行名是gene symbol名，列名是位置信息，类似于坐标（空缺的部分是被过滤掉的低质量细胞）
```
#去重复
dataA <- dataA %>%
  distinct(Genes,.keep_all = T) %>% 
  column_to_rownames("Genes")

#seurat分析
stRNA <- CreateSeuratObject(counts = dataA)
stRNA <- NormalizeData(stRNA, normalization.method = "LogNormalize", scale.factor = 10000)
stRNA <- FindVariableFeatures(stRNA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(stRNA)
stRNA <- ScaleData(stRNA, features = all.genes)

stRNA <- RunPCA(stRNA, features = VariableFeatures(object = stRNA))

stRNA <- FindNeighbors(stRNA, dims = 1:5)

stRNA <- FindClusters(stRNA, resolution = 0.3)

stRNA <- RunUMAP(stRNA, dims = 1:5)
plot <- DimPlot(stRNA)
data <- plot$data
```
![image](https://user-images.githubusercontent.com/112565216/195789938-325f21af-9120-4c98-a049-492ff7320cf7.png)

```
sp_1 <- colnames(dataA)
sp_2 <- str_split(sp_1, "x")
index_get <- function(x){
  x <- as.numeric(x)
  x <- as.vector(x)
}
ind <- as.data.frame(t(sapply(sp_2, index_get)))
names(ind) <- c("row_ind", "col_ind")
rownames(ind) <- paste0(ind$row_ind,"x",ind$col_ind)
#横坐标从0开始
ind$col_ind <- ind$col_ind - 1  
ind$row_ind <- ind$row_ind - 1
ind_mat <- dplyr::arrange(ind, col_ind, row_ind)
data <- data[rownames(ind_mat),]
ind_mat$ident <- data$ident
ind_mat$shape <- "A"
```
![image](https://user-images.githubusercontent.com/112565216/195800740-23e3b417-65fe-4449-a67d-127228c3732d.png)

```
#归一化处理
data <- dataA[,rownames(ind_mat)]
tmp2 <- 1e6 / colSums(data)
ST_norm <- sweep(data, 2, tmp2, '*')
ST_c <- log(ST_norm + 1)
```
![image](https://user-images.githubusercontent.com/112565216/195811167-f537045d-d966-4c41-81bc-a011238fdd41.png)

```
##   PRSS1 基因的分布图
data.gene <- ind_mat
data.gene$gene <- as.numeric(ST_c["PRSS1",])

data.gene <- data.gene %>% 
  mutate(zscore = (gene - mean(gene))/sd(gene))

gg <- ggplot(data.gene, aes(x = row_ind, y = col_ind))

pl <- gg + geom_point(data = data.gene, aes(size = 3, 
                                            shape = shape,
                                            color=zscore)) +
  scale_shape_manual(values = c(15))+
  scale_colour_gradientn(colors=c("#D0EBF8",
                                  "#9BCCDF",
                                  "#9BCCDF",
                                  "#9BCCDF",
                                  "#D0EBF8",
                                  "#9BCCDF",
                                  "#9BCCDF",
                                  "#9BCCDF",
                                  "#9BCCDF",
                                  "#FFEE9F",
                                  "#A11D2A",
                                  "red"))+
  theme_void()+
  guides(size=F,shape=F,color=F)
pl
```
![image](https://user-images.githubusercontent.com/112565216/195817147-3b38a1bc-3d5d-4aab-a569-88f6ecfcd895.png)

```
#降维聚类图
gg <- ggplot(ind_mat, aes(x = row_ind, y = col_ind))

region <- gg + geom_point(data = ind_mat, aes(size = 3, 
                                              shape = shape,
                                              color=ident)) +
  scale_shape_manual(values = c(15)) +
  scale_colour_manual(values = c("#189D77",
                                 "#666666",
                                 "#E42A88",
                                 "#766FB1"))+
  theme_void() +
  guides(size=F,shape=F)

region

```
![image](https://user-images.githubusercontent.com/112565216/195819572-c74bd868-a150-4d23-ad29-6722261401ec.png)

```
#单细胞转录组分析
scdataA <- fread("GSE111672_PDAC-A-indrop-filtered-expMat.txt/GSE111672_PDAC-A-indrop-filtered-expMat.txt",header=T)
```
![image](https://user-images.githubusercontent.com/112565216/195837022-ceaf1827-8680-428e-8b4b-9c93264b1f3c.png)

```
cell <- data.frame(cell=c("Acinar cells", "Ductal - terminal ductal like", "Ductal - CRISP3 high/centroacinar like", 
                          "Cancer clone A", "Ductal - MHC Class II", "Cancer clone B", 
                          "mDCs A", "Ductal - APOL1 high/hypoxic", "Tuft cells", "mDCs B", 
                          "pDCs", "Endocrine cells", "Endothelial cells", "Macrophages A", 
                          "Mast cells", "Macrophages B", "T cells & NK cells", "Monocytes", 
                          "RBCs", "Fibroblasts"),
                   celltype=c("Acinar cells", "Ductal", "Ductal", 
                              "Cancer clone A", "Ductal", "Cancer clone B", 
                              "mDCs", "Ductal", "Tuft cells", "mDCs", 
                              "pDCs", "Endocrine cells", "Endothelial cells", "Macrophages", 
                              "Mast cells", "Macrophages", "T cells & NK cells", "Monocytes", 
                              "RBCs", "Fibroblasts"))

scdataA <- as.data.frame(scdataA)
names <- colnames(scdataA)[2:ncol(scdataA)]
colnames(scdataA) <- paste0(colnames(scdataA),"_",1:1927)

scdataA <- scdataA %>%
  distinct(Genes_1,.keep_all = T) %>% 
  column_to_rownames("Genes_1")

scRNA <- CreateSeuratObject(counts = scdataA)

scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)

scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))

scRNA <- FindNeighbors(scRNA, dims = 1:15)
scRNA <- FindClusters(scRNA, resolution = 0.5)

scRNA <- RunUMAP(scRNA, dims = 1:15)

scRNA@meta.data$cell <- names
scRNA@meta.data$id <- NA


for(i in 1:nrow(cell)){
  scRNA@meta.data[which(scRNA@meta.data$cell == cell$cell[i]),'id'] <- cell$celltype[i]
}
Idents(scRNA) <- 'id'
DimPlot(scRNA)

```
![image](https://user-images.githubusercontent.com/112565216/195837773-a50bee08-4265-43e2-b380-2c921f9dfe6d.png)

![image](https://user-images.githubusercontent.com/112565216/195837660-db4efb09-6ef0-4320-8934-5e0c09ac7bb1.png)


```
scRNA_marker <- FindAllMarkers(scRNA,test.use = "t",
                               only.pos = T,
                               logfc.threshold = 0.05)

scRNA_marker2 <- subset(scRNA_marker)


#对空间转录组聚类分析
newid <- c("Stroma","Cancer","Duct","Pancreatic")
names(newid) <- levels(stRNA)
stRNA <- RenameIdents(stRNA, newid)

stRNA_marker <- FindAllMarkers(stRNA,
                               only.pos = T,
                               test.use="t",
                               logfc.threshold = 0.1)

stRNA_marker2 <- subset(stRNA_marker)




#将spatial与scRNA数据整合
p.data <- data.frame(sc=NA,st=NA,p=NA)
for (i in unique(scRNA_marker2$cluster)) {
  sub_cluster = subset(scRNA_marker2,cluster==i)
  gene = sub_cluster$gene
  gene_l = length(gene)
  for (j in unique(stRNA_marker2$cluster)) {
    sub_cluster_st = subset(stRNA_marker2,cluster==j)
    gene_st = sub_cluster_st$gene
    gene_st_l = length(gene_st)
    gene_inter_l = length(intersect(gene_l,gene_st_l))
    p = -log10(1-phyper(gene_inter_l,
                        gene_st_l, 
                        19738-gene_st_l, 
                        gene_l))
    p.demo.data = data.frame(sc = paste0(i,"(",gene_l,")"),
                             st = paste0(j,"(",gene_st_l,")"),p = p)
    p.data = rbind(p.demo.data,p.data)
    p.data = na.omit(p.data)
  }
}
p.data$group <- p.data$st

## 画图
p.data <- as_tibble(p.data)
p.data %>%
  tidyHeatmap::heatmap(
    .column = st,
    .row = sc,
    .value = p,
    .scale = "both",
    show_column_names=F,
    palette_value = c("#1B4584","white","#660621"),
    rect_gp = grid::gpar(col = "black", lwd = 1),
    row_names_side = "left",
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )%>%
  add_tile(group,palette = c("#189D77",
                             "#666666",
                             "#E42A88",
                             "#766FB1")) 



```
![image](https://user-images.githubusercontent.com/112565216/195845143-f4ec921c-bcdf-4b52-86ad-bdb7c8769e95.png)
