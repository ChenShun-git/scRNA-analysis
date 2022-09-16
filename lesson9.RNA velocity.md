RNA velocity是指未经剪切的前体mRNA与已剪切的成熟mRNA含量的比例，来预测基因的动态变化。一般来说，如果未经剪切的前体mRNA较多，则该基因的表达量增加，反之则表达量降低
# velocyto
## 1.velocyto生成loom文件
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

## 2.安装并打开rstudio-server
根据rstudio官网指示下载R与rstudio-server后用"sudo rstudio-server start"运行
如果需要更换不同版本的R，则在终端内输入：
```
sudo vi /etc/rstudio/rserver.conf
```
并将目标R版本所在的路径输入其中
```
rsession-which-r=/home/chen/anaconda3/envs/velocyto.r/bin/R
```

