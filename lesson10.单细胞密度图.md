# 针对某个基因的表达量的密度图
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
setwd("C:/Users/86269/Desktop/shun.C/single_cell")
load("scRNA_harmony.rData")
scRNA <- scRNA_harmony
#准备画图使用的数据集
x=as.data.frame(scRNA@active.ident)
df <- scRNA@reductions$umap@cell.embeddings
df<- cbind(df,x)
colnames(df) <- c("x","y","ident")
```
![image](https://user-images.githubusercontent.com/112565216/188883038-1cef0e04-c440-4322-912d-e8e43b7e8d88.png)

```
#挑选两个基因进行密度图绘制
library(Nebulosa)
plot_density(scRNA,c("CD4","CD8A"),joint=TRUE)
```
![image](https://user-images.githubusercontent.com/112565216/188883276-f0841757-beb2-4b67-811f-9a323dbf86c3.png)

最右图是将两个基因的表达量密度合并，体现两个基因的共表达情况。对于marker基因，将表达量密度合并后作密度图，如果在图中出现明显的密度较高的区域，也可以作为判定cluster中细胞类型的依据

# 针对细胞数量的密度图
```
library(ggpointdensity)
ggplot(df,aes(x,y)) + geom_pointdensity(adjust=4)
```
![image](https://user-images.githubusercontent.com/112565216/188883872-84bea81f-7401-4189-9073-4dac723bd77e.png)

**自定义方法作密度图**
```
x=as.data.frame(scRNA@active.ident)
df1 <- scRNA@reductions$umap@cell.embeddings
df1<- cbind(df1,as.data.frame(scRNA@meta.data))
```
![image](https://user-images.githubusercontent.com/112565216/189041820-59ac0cf4-98b1-4e46-8c23-92dcc558bc37.png)

```
#自定义一个计算density的函数（※）
get_density <-function(x,y,...){
  dens <- MASS::kde2d(x,y,...)
  ix <- findInterval(x,dens$x)
  iy <- findInterval(y,dens$y)
  ii <- cbind(ix,iy)
  return(dens$z[ii])
}
#计算各sample的细胞密度
for(i in unique(df1$orig.ident)){
  id <- which(df1$orig.ident==i)
  gd <- get_density(df1$UMAP_1[id],df1$UMAP_2[id],n=100)
  df1[id,"Density"] <- gd/max(gd)
}
#绘图
library(ggplot2)
ggplot(df1,aes(UMAP_1,UMAP_2,col=Density))+
  geom_point(size=0.1,shape=".",alpha=0.5)+
  theme_classic(base_size=15)+
  facet_grid(.~orig.ident)+
  scale_color_viridis_c(option="C",direction=1,breaks=seq(0,1,0.25))+
  theme(legend.position="bottom",
        strip.text.y.right=element_text(angle=0),
        strip.text=element_text(size=15))+
  guides(color=guide_colorbar(title.vjust=1,
                              barwidth=12,
                              barheight=1.5))
```
![image](https://user-images.githubusercontent.com/112565216/189042087-2d106806-edc3-4d1b-9309-886919c38b95.png)
