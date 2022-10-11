# cellphonedb
cellphonedb是包含受体、配体与其相互作用的数据库，可对细胞间的通讯分子进行全面系统的分析，研究不同类型的细胞间的交流与通信网络

**1.准备cellphonedb分析使用的数据**
cellphonedb分析的核心数据是counts.txt与meta.txt
```
#################注释细胞################
Sys.setenv(LANGUAGE="en")
library(data.table)
library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(SingleR)
library(celldex)
setwd("C:/Users/86269/Desktop/shun.C/single_cell/")
load("scRNA_harmony.rData")

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
scRNA_harmony@meta.data$celltype="NA"
for(i in 1:nrow(celltype)){
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$seurat_clusters==celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
  }
  
#############制作meta.txt与count.txt
#选取1000个细胞进行分析
bc <- sample(colnames(scRNA_harmony),1000)
scRNA <- scRNA_harmony[,bc]
meta=scRNA@meta.data
meta$Cell=rownames(meta)
#准备meta.txt数据
meta = meta[,c("Cell","celltype")]
fwrite(meta,file="meta.txt",sep="\t")
#准备counts.txt数据
counts=scRNA@assays$RNA@data
counts=as.data.frame(counts)
counts1=as.data.frame(counts)
counts$Gene=rownames(counts)
counts3=counts[,c("Gene",colnames(counts1))]
#将基因名转换为ENSEMBL 
library(org.Hs.eg.db)
x2=AnnotationDbi::select(org.Hs.eg.db,keys=counts3$Gene,columns = c("ENSEMBL","SYMBOL"),keytype="SYMBOL")
x3=merge(x2,counts3,by.x="SYMBOL",by.y="Gene")
#去除SYMBOL基因行
x4=x3[,-1]
#去除重复值
x5=distinct(x4,ENSEMBL,.keep_all = T)
fwrite(x5,file="counts.txt",sep="\t")
```
![image](https://user-images.githubusercontent.com/112565216/194984385-6f2f4242-cc02-41f5-a9fb-7d0327cebada.png)
![image](https://user-images.githubusercontent.com/112565216/194984732-ea28ef5e-a10b-4736-a94d-a7b334ca8602.png)

[meta(上图)count(下图)]

**2.进行cellphonedb分析**
在ubuntu中安装cellphonedb
```
cellphonedb method statistical_analysis meta.txt counts.txt  --iterations=1000 --threads=6
```
输出结果如下

![image](https://user-images.githubusercontent.com/112565216/194985926-69986373-faf8-4f7c-be34-886faa2ed5da.png)

**3.对cellphonedb的输出结果进行分析**
```
library(igraph)
library(qgraph)
library(psych)
cellnet <- read.delim("out/count_network.txt",check.names=FALSE)
#去除0值
cellnet <- cellnet %>% filter(count>0)
net <- graph_from_data_frame(cellnet)
plot(net)
```
![image](https://user-images.githubusercontent.com/112565216/194989019-7a7692ea-09ba-4711-8345-7b319cb57465.png)

```
#######################自定义互作通路与细胞类型#######################
all_pval <- read.table("out/pvalues.txt",header=T,stringsAsFactors = F,sep='\t',comment.char ='',check.names = F )
all_means <- read.table("out/means.txt",header=T,stringsAsFactors = F,sep='\t',comment.char ='',check.names = F )
#####提取作图的受体配体对和互作的细胞类型
intr_pairs =all_pval$interacting_pair
#去除不在分析里用到的列
all_pval <- all_pval[,-c(1:11)]
all_means <- all_means[,-c(1:11)]
```
![image](https://user-images.githubusercontent.com/112565216/194989525-e64f772e-83f5-4168-b6e6-547a7499b3f3.png)
![image](https://user-images.githubusercontent.com/112565216/194989610-3e6ab8a8-5b07-4e6b-bb33-b8f6c25365bf.png)

[all_pval(上图)all_means(下图)]

```
#随机挑选10个互作通路
selected_rows <- read.table("out/means.txt",header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
selected_rows <- sample(selected_rows$interacting_pair,10)
#随机挑选互作的细胞类型
selected_columns <- sample(colnames(all_means),8)
df_names=expand.grid(selected_rows,selected_columns)
```
![image](https://user-images.githubusercontent.com/112565216/194990835-a565af70-45b8-42d9-b75e-5bf7d7cedb9f.png)

[df_names]
```
#加入对应pvalue与mean值
sel_pval =all_pval[match(selected_rows,intr_pairs),selected_columns]
sel_means =all_means[match(selected_rows,intr_pairs),selected_columns]
pval=unlist(sel_pval)
pval[pval==0]=0.0009
plot.data=cbind(df_names,pval)
pr=unlist(as.data.frame(sel_means))
pr[pr==0]=1
plot.data=cbind(plot.data,log2(pr))
colnames(plot.data)=c("pair","cluster","pvalue","mean")
```
![image](https://user-images.githubusercontent.com/112565216/194990955-8ffa649b-026c-4ef7-b88e-d221dfb83ef2.png)

[plot.data]

```
my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

ggplot(plot.data,aes(x=cluster,y=pair))+
  geom_point(aes(size=-log10(pvalue),color=mean))+
  scale_color_gradientn("Log2 mean (Molecule 1, Molecule 2)",colors=my_palette)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text=element_text(size=14,colour="black"),
        axis.text.x=element_text(angle=90,hjust=1),
        axis.text.y=element_text(size=12,colour="black"),
        axis.title=element_blank(),
        panel.border=element_rect(size=0.7,linetype="solid",colour="black")) 


```
![image](https://user-images.githubusercontent.com/112565216/194993807-8170f4f4-aa87-4d2f-b232-a8d098e12a18.png)

# cellchat
```
library(ggalluvial)
library(ggsci)
library(CellChat)
#创建cellchat对象
data.input <- GetAssayData(scRNA_harmony,assay="RNA",slot="data")
identity<- subset(scRNA_harmony@meta.data,select="celltype")
cellchat <- CreateCellChat(object=data.input,meta=identity,group.by="celltype")

#加载参考数据库（可加载CellChatDB.human/CellChatDB.mouse）
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
```
![image](https://user-images.githubusercontent.com/112565216/195010500-b2b53be3-96e0-4809-87a1-0213e38bb73a.png)

```
#选择旁分泌通路作为参考数据库
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use 
#将信号基因子集化
cellchat <- subsetData(cellchat)
#识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
#识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
#将配体受体投射到PPI(protein-protein interaction)网络
cellchat <- projectData(cellchat,PPI.human)
#针对配体对水平计算互作强度
cellchat <- computeCommunProb(cellchat,raw.use = TRUE)
#过滤细胞数过少的通讯
cellchat <- filterCommunication(cellchat,min.cells = 3)
#汇总通信概率计算细胞间聚合通信网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#计算聚合细胞互作通信网络
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow=c(1,2),xpd=TRUE)
netVisual_circle(cellchat@net$count,vertex.weight = groupSize,weight.scale = T,label.edge=F,title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight,vertex.weight = groupSize,weight.scale = T,label.edge=F,title.name = "interaction weight")
```
![image](https://user-images.githubusercontent.com/112565216/195022807-23312f2b-e5e6-4b76-a725-8431435dc30b.png)

**1.细胞层面**
```
#可视化每种细胞与其他细胞的互作结果
mat <- cellchat@net$weight
par(mfrow = c(2,4),xpd=TRUE)
for(i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow=nrow(mat),ncol=ncol(mat),dimnames=dimnames(mat))
  mat2[,i] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale=T,edge.weight.max=max(mat),title.name=rownames(mat)[i])
}
```
![image](https://user-images.githubusercontent.com/112565216/195029820-3a069fd0-1805-4fa5-b917-c121363c1315.png)

**2.通路层面**
```
#####可视化每个信号通路
#查看细胞顺序
levels(cellchat@idents)
```
![image](https://user-images.githubusercontent.com/112565216/195032703-d05fd8d7-62a8-42b4-b6ca-9895c9a1a9d2.png)

```
#指定靶细胞索引
vertex.receiver=seq(1,2)
#指定需要展示的通路
netVisual_aggregate(cellchat,signaling = "MK",vertex.receiver=vertex.receiver,layout="hierarchy")
```
![image](https://user-images.githubusercontent.com/112565216/195032505-63f3d59e-5251-453e-9513-e6936779d7ca.png)

```
#和弦图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat,signaling = "MK",layout="chord",vertex.size = groupSize)

#圈图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat,signaling = "MK",layout="circle",vertex.size = groupSize)

#热图
par(mfrow=c(1,1))
netVisual_heatmap(cellchat,signaling ="MK",color.heatmap="Reds")
```
![image](https://user-images.githubusercontent.com/112565216/195034706-f9a6985f-d67f-426b-a6db-158c56ea1e28.png)
![image](https://user-images.githubusercontent.com/112565216/195034762-b0d16fe6-37c2-417e-9825-b378c84686d1.png)
![image](https://user-images.githubusercontent.com/112565216/195034859-2bfafa99-762e-4f3e-8eba-12495e7d7b23.png)

**3.配体-受体层面**
```
#配体-受体层级的可视化
#计算各ligand-receptor+pair对信号通路的贡献程度
#以CCL通路为例
pathways.show <- "CCL"
netAnalysis_contribution(cellchat,signaling=pathways.show)
#气泡图
netVisual_bubble(cellchat,sources.use = 2,targets.use = c(1:5),remove.isolate=FALSE)
#用小提琴图绘制信号基因的表达分布
plotGeneExpression(cellchat,signaling = "CCL",enriched.only=FALSE)
```
![image](https://user-images.githubusercontent.com/112565216/195035393-8f324280-5666-4b2e-8add-aa1ed1e7634a.png)
![image](https://user-images.githubusercontent.com/112565216/195036665-0c1f3ef8-5ed4-4707-a039-917c2475f834.png)
![image](https://user-images.githubusercontent.com/112565216/195037332-3ea8e32c-4132-4e65-9415-f30500f2b7e3.png)

```

**4.两个分组进行差异分析**
```
#依据分组进行切割并创造cellchat对象
sc.sp=SplitObject(scRNA_harmony,split.by = "orig.ident")
sc.11=scRNA_harmony[,sample(colnames(sc.sp[["sample_11"]]),1000)]
sc.3=scRNA_harmony[,sample(colnames(sc.sp[["sample_3"]]),1000)]

cellchat.sc11 <- createCellChat(object =sc.11@assays$RNA@data, meta =sc.11@meta.data,  group.by ="celltype")
cellchat.sc3 <- createCellChat(object =sc.3@assays$RNA@data, meta =sc.3@meta.data,  group.by ="celltype")

#对两组数据分别进行分析
#sample11
cellchat=cellchat.sc11 
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc11 = cellchat

#sample3
cellchat=cellchat.sc3
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cc.sc3 = cellchat

#将分别处理后的结果合并
cc.list=list(SC11=cc.sc11,SC3=cc.sc3)
cellchat=mergeCellChat(cc.list,cell.prefix = T,add.names = names(cc.list))
```

```
#可视化
##所有细胞群总体观：通讯数量与强度对比
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "count")
compareInteractions(cellchat,show.legend = F,group = c(1,3),measure = "weight")
##第一个图展示通讯数量之间的差异，第二个图展示通讯强度之间的差异。 
```
![image](https://user-images.githubusercontent.com/112565216/195054095-8a258e70-d420-457d-8f35-a591925a7710.png)
![image](https://user-images.githubusercontent.com/112565216/195054139-10aa5b52-3746-4002-94e3-b31393684a22.png)

```
##数量与强度差异网络图
netVisual_diffInteraction(cellchat,weight.scale = T)
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight")
##红色是case相对于control上调的，蓝色是下调的

```
![image](https://user-images.githubusercontent.com/112565216/195054298-a716dbde-9b74-42ac-8a08-84edde954f32.png)
![image](https://user-images.githubusercontent.com/112565216/195054380-0038e133-51dc-4079-af4b-55b63769b9f3.png)

```
#数量与强度差异热图
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat,measure = "weight")
#case和control对比，红色是上调，蓝色是下调
```
![image](https://user-images.githubusercontent.com/112565216/195055039-ee32cdf5-86c2-464a-b572-8ca317f656cf.png)
![image](https://user-images.githubusercontent.com/112565216/195055079-189b586f-3242-4124-9afe-30c139115980.png)

```
#保守和特异性信号通路的识别与可视化
rankNet(cellchat,mode = "comparison",stacked = T,do.stat = T)
rankNet(cellchat,mode = "comparison",stacked =F,do.stat = T)
```
![image](https://user-images.githubusercontent.com/112565216/195056278-80b1a149-478e-4aae-83ec-f27262b4b9a9.png)
![image](https://user-images.githubusercontent.com/112565216/195056324-1a7ea2b7-12ae-4136-b4e8-819485ab4ef6.png)

```
##细胞互作数量对比网络图
weight.max=getMaxWeight(cc.list,attribute = c("idents","count"))
netVisual_circle(cc.list[[1]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc11" )

netVisual_circle(cc.list[[2]]@net$count,weight.scale = T,label.edge = F,
                 edge.weight.max =weight.max[2],edge.width.max = 12,title.name = "sc3" )
```
![image](https://user-images.githubusercontent.com/112565216/195057354-e788ef86-99cc-4fcb-9e59-060dac483820.png)
![image](https://user-images.githubusercontent.com/112565216/195057430-c50cfa3d-cba6-4efc-af2c-8648893d08eb.png)

```
#选定receptor，target与通路，并比较特定细胞种类与通路情况下不同样本中的全部L-R配对情况
netVisual_bubble(cellchat,sources.use = 4,targets.use = 5:11,comparison = c(1:2),angle.x = 45)
```
![image](https://user-images.githubusercontent.com/112565216/195059582-2bc296db-01f0-45ed-bedd-4d70fe116d4a.png)



