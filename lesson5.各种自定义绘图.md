# 1.小提琴图标注p值（接lesson4）
仍以基因Rbp1的表达为例
```
#提取Rbp1的表达量并与metadata中其他信息合并，使metadata中增加一列表示每个细胞的Rbp1表达量
scFURIN =subset(scRNA1,Rbp1>0)
FURIN_exp=scFURIN@assays$RNA@data["Rbp1",]
FURIN_exp.m =scFURIN@meta.data
FURIN_exp.m$exp=FURIN_exp
FURIN_exp.m$Cohort=factor(FURIN_exp.m$orig.ident,levels=c("dapi1","dapi2","CAF1","CAF2"))
#计算p值
library(ggsignif)
#两两计算的分组
my_comparisons=list(c("DapiNeg1","CAF1"),c("DapiNeg2",'CAF2'),
                    c("DapiNeg2" ,"CAF1"))
#可视化
ggplot(FURIN_exp.m,aes(Cohort,exp,fill=Cohort))+
  geom_violin()+theme_classic()+
  theme(plot.title = element_text(size=15,face='bold'),
        axis.text.x = element_text(size=15,face='bold',angle=25,hjust=1,vjust=1),
        axis.text.y = element_text(size=18,face='bold'),
        axis.title.x =element_text(size=0,vjust=1,face='bold'),
        axis.title.y = element_text(size=15,face='bold')) +
  labs(y='Expression',x="Group")+
  geom_signif(comparisons=my_comparisons,step_increase=0.1,map_signif_level=F,test=t.test,size=1,textsize=6)+
  geom_jitter(shape=16,position=position_jitter(0.2))+
  ggtitle("Rbp1")+theme(plot.title=element_text(size=18,hjust=0.5))+
  guides(fill=guide_legend(title=' '))+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1))+
  theme(axis.line.y=element_line(linetype=1,color="black",size=1))+
  theme(legend.text=element_text(size=15))
```
![violin plot with p value](https://upload-images.jianshu.io/upload_images/28382212-e3248d3412e98368.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

# 2.比例图标注p值
## 画比例图
加载数据与r package
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
setwd("C:/Users/86269/Desktop/shun.C/single_cell")
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)

load("sc.combined.rdata")
```
为画图挑选颜色
```
library(RColorBrewer)
qual_col_pals=brewer.pal.info[brewer.pal.info$category=='qual',]
col_vector = unlist(mapply(brewer.pal,qual_col_pals$maxcolors,rownames(qual_col_pals)))
```
提取细胞的样品来源、cluster以及数量
```
library(reshape2)
pB2_df <- table(sc.combined@meta.data$new.celltype.2,sc.combined@meta.data$type) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")

sample_color <- col_vector[1:9]
```
可视化
```
pB2 <- ggplot(data = pB2_df,aes(x=Cluster,y=Number,fill=Sample))+
  geom_bar(stat="identity",width=0.8,position="fill")+
  scale_fill_manual(values=sample_color)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x="",y="Ratio")+
  theme(axis.text.y=element_text(size=12,colour="black"))+
  theme(axis.text.x=element_text(hjust=1,vjust=1,angle=45))
pB2
```
![比例图](https://upload-images.jianshu.io/upload_images/28382212-62f03c2991c25180.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

该图存在的问题是，三种sample的细胞总数不同，因此正常情况下必须从不同sample中随机抽取同等数量的细胞作比例图才有意义。但是此处三种sample中细胞总数分别为

![image.png](https://upload-images.jianshu.io/upload_images/28382212-8148fa4b1e13f4f3.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

数量差别不大，因此可以省去随机抽样这一步






对换横纵坐标
```
pB3 <- ggplot(data=pB2_df,aes(x=Number,y=Cluster,fill=Sample))+
  geom_bar(stat='identity',width=0.8,position='fill')+
  scale_fill_manual(values=sample_color)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x="",y="Ratio")+
  theme(axis.text.y=element_text(size=12,colour="black"))+
  theme(axis.text.x=element_text(size=12,colour="black"))+
  theme(
    axis.text.x.bottom =element_text(hjust=1,vjust=1,angle=45)
    
  )
pB3
```
![比例图](https://upload-images.jianshu.io/upload_images/28382212-7cd085109495d1d5.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
将填充颜色对应为细胞种类（即cluster）
```
pB4 <- ggplot(data=pB2_df,aes(x=Number,y=Sample,fill=Cluster))+
  geom_bar(stat='identity',width=0.8,position='fill')+
  scale_fill_manual(values=col_vector[1:20])+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x="",y="Ratio")+
  theme(axis.text.y=element_text(size=12,colour="black"))+
  theme(axis.text.x=element_text(size=12,colour="black"))+
  theme(
    axis.text.x.bottom =element_text(hjust=1,vjust=1,angle=45)
    
  )
pB4
```
![比例图](https://upload-images.jianshu.io/upload_images/28382212-f4d515f190731ca4.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

## 在图中标注误差线
```
#计算细胞比例
x <- table(sc.combined@meta.data$new.celltype.2,sc.combined@meta.data$orig.ident)
#计算每个sample中每种细胞的占比
x3<-t(t(x)/rowSums(t(x)))
x4<-as.data.frame(as.table(t(x3)))
colnames(x4)=c("sample","celltype","Freq")
x4$group =x4$sample %>%str_replace("_.","")

#计算误差线上下值
top<-function(x){
  return(mean(x)+sd(x)/sqrt(length(x)))
}
bottom <- function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}

#每个细胞类型在各分组中的情况，以上皮细胞为例
dose_epi <- x4[which(x4$celltype=="Epithelial"),]
 
ggplot(data=dose_epi,aes(x=group,y=Freq,fill=group))+
  stat_summary(geom="bar",fun='mean',
               position=position_dodge(0.9))+
  stat_summary(geom='errorbar',fun.min = bottom,fun.max = top,position=position_dodge(0.9),width=0.2)+
  scale_y_continuous(expand=expansion(mult = c(0.01)))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x="Group",y="Proportion")+
  geom_point(data=dose_epi,aes(group,Freq),size=1,pch=19)+
  theme(
    axis.text.x.bottom=element_text(hjust=1,vjust=1,angle=45)
    
  )+ggtitle("Epithelial")

```

![带误差线的柱状图，points表示初始数据中单个样本上皮细胞的占比](https://upload-images.jianshu.io/upload_images/28382212-163f70ffed9f780a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
![dose_epi](https://upload-images.jianshu.io/upload_images/28382212-501b3f2acf65be70.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


# 3.聚类图上加标签
加载数据
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
setwd("C:/Users/86269/Desktop/shun.C/single_cell")
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)

load("scRNA1.Rdata")
```
将聚类结果与降维结果整合
```
x=as.data.frame(scRNA1@active.ident)
df<-scRNA1@reductions$tsne@cell.embeddings
df <- cbind(df,x)
colnames(df)<-c("x","y","ident")
```
![df](https://upload-images.jianshu.io/upload_images/28382212-3e925167e89652c5.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
初步画出聚类图
```
p<- ggplot(df,aes(x=x,y=y))+geom_point(aes(colour=ident),size=1)+
  labs(x="t-SNE 1",y="t-SNE 2")
p
```
![p](https://upload-images.jianshu.io/upload_images/28382212-340a4ea56fec6609.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

删去分组信息
```
p1 <- p+NoLegend()
p1
```
![p1](https://upload-images.jianshu.io/upload_images/28382212-fa86a8a92bba5370.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
```
#加上metadata中的其他信息
meta=scRNA1@meta.data
meta=meta[,c(1:6)]
df.m=merge(df,meta,by=0)
```
![df.m](https://upload-images.jianshu.io/upload_images/28382212-3a9a38713bb746ea.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


在图中加上celltype信息
```
df.m <- df.m %>% group_by(celltype) %>%summarise(x=median(x),y=median(y))
ggplot(df,aes(x,y))+
  geom_point(aes(colour=ident))+
  ggrepel::geom_text_repel(aes(label=celltype),
                           data=df.m,
                           size=5,
                           label.size=1,
                           segment.color=NA)+
  theme(legend.position = "none")+theme_bw()
```
![](https://upload-images.jianshu.io/upload_images/28382212-92c9cfe413660722.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
在图中加上cluster信息
```
df.m <- df.m %>% group_by(seurat_cluster) %>%summarise(x=median(x),y=median(y))
ggplot(df,aes(x,y))+
  geom_point(aes(colour=ident))+
  ggrepel::geom_text_repel(aes(label=seurat_clusters),
                           data=df.m,
                           size=5,
                           label.size=1,
                           segment.color=NA)+
  theme(legend.position = "none")+theme_bw()
```
![](https://upload-images.jianshu.io/upload_images/28382212-2cc03b35f4208538.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

# 4.在聚类图上画椭圆线（※）
以p为例
![p](https://upload-images.jianshu.io/upload_images/28382212-253073334a3029bd.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
```
obs <- p$data[,c("x","y","ident")]
ellipse_pro <- 0.98
theta <- c(seq(-pi,pi,length=50),seq(pi,-pi,length=50))
circle <- cbind(cos(theta),sin(theta))

ell <- ddply(obs,'ident',function(x){
  if(nrow(x)<=2){
    return(NULL)
  }
  sigma<-var(cbind(x$x,x$y))#求方差
  mu<-c(mean(x$x),mean(x$y))#求平均数
  ed <-sqrt(qchisq(ellipse_pro,df=2))
  data.frame(sweep(circle %*% chol(sigma)*ed,2,mu,FUN="+"))#chol是求一个正定矩阵的上三角矩阵,%*%是矩阵相乘
})
names(ell)[2:3]<- c('x','y')
ell <-ddply(ell,.(ident),function(x)x[chull(x$x,x$y),])
p1 <- p+geom_polygon(data=ell,aes(group=ident),
                     colour= 'black',
                       alpha=1,fill=NA,
                       linetype="dashed")
p1


```
![](https://upload-images.jianshu.io/upload_images/28382212-5b11cd7126e298d7.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


# 5.在聚类图中将某些点进行highlight
取某些细胞在聚类图中用不同的颜色标注
```
#以scRNA1前100个细胞为例
mutlist <- colnames(scRNA1)[1:100]
DimPlot(scRNA1,reduction="tsne",
        cells.highlight = mutlist,
        cols.highlight='#11AA4D',
        sizes.highlight=1.5,
        group.by='orig.ident')
```
![](https://upload-images.jianshu.io/upload_images/28382212-1776b4e876fa7367.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

# 6.3d PCA聚类图
加载数据与package
```
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(Seurat)
library(scatterplot3d)
setwd("C:/Users/86269/Desktop/shun.C/single_cell")
load("sc.combined.rdata")
```
数据预处理
```
#提取三种不同类型的细胞
sc.1=subset(sc.combined,idents = c("Epithelial","Endothelial","Stem_cell"))
#将分组依据换为sample
Idents(sc.1)="orig.ident"
#用取平均值的方式将每个sample的每种细胞合并成一个表达矩阵(将每个sample的每个细胞种类的所有细胞的基因表达量加和，再除以细胞数)
x=AverageExpression(sc.1 ,add.ident ="new.celltype.2")
x=x$RNA
```
![x](https://upload-images.jianshu.io/upload_images/28382212-bca793f436cf69ff.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
进行PCA降维
```
#PCA降维
pca <- prcomp(t(x))
#提取每个细胞类型的降维结果
pca.data=pca$x
#提取画图需要的PC1,PC2,PC3
pca.data=as.data.frame(pca.data)
pca.data=pca.data[,1:3]

s1 <- strsplit(rownames(pca.data),split = "_")
s2 <- sapply(s1,function(x){x[1]})
s3 <- sapply(s1,function(x){x[3]})
pca.data$Type <- s2
pca.data$celltype = s3

```
![pca.data](https://upload-images.jianshu.io/upload_images/28382212-c3ccb89aaf0ba336.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
可视化
```
###不同type用形状区分 不同celltype用颜色区分
Type.p=c(rep(11,9),rep(16,9),rep(17,9))
cell.type.p = c(rep(c("blue","red","orange"),9 ))

colors.lib <- c("blue","red","orange")
shapes.lib = c(11,16,17)

###画图
source('huage.3D.plot.R')
s3d <- scatterplot3d(pca.data[,c("PC1","PC2","PC3")],
                     pch = Type.p,color = cell.type.p,                   
                     cex.symbols = 1,grid=FALSE, box=FALSE,                     
                     main = "3D PCA plot")
#加入celltype信息
legend("topright",
       c("Epithelial","Endothelial","Stem_cell"),
       fill=c('blue',"red","orange"),
       box.col=NA,
       cex=0.6)

#加入type信息
legend("topleft",
       c('GF','SPF','FMT'),
       pch = shapes.lib,
       box.col=NA,
       cex =0.6,
       y.intersp = 0.5)
#加上平面
addgrids3d(pca.data[,c("PC1","PC2","PC3")], grid = c("xy", "xz", "yz"))
```
![3D PCA
plot](https://upload-images.jianshu.io/upload_images/28382212-2ea8248bce0d8c91.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

从图中可看出不同分组里哪些细胞类型差异较大，分散越明显差异越大
