scanpy是基于python的单细胞数据分析包,简要介绍使用scanpy分析的pipeline
### 1.加载分析需要的库
```
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity 的取值表示测试结果显示的详细程度，数字越大越详细
sc.logging.print_header()             #输出版本型号

sc.settings.set_figure_params(dpi=80, facecolor='white')#设置图片分辨率大小及样式
```

### 2.读取10x数据
```
adata = sc.read_10x_mtx(
    "pbmc3k",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True) 
```

### 3.质控
```
#将表达量最高的20个基因可视化
sc.pl.highest_expr_genes(adata, n_top=20, )

```
![image](https://user-images.githubusercontent.com/112565216/188346714-9c67c9d5-0ed0-4628-9ac7-a7195c683281.png)

```
#初步过滤基因与细胞
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
```
![image](https://user-images.githubusercontent.com/112565216/188347033-9cb71448-84ba-4281-b241-1f3239abb172.png)

过滤结果显示过滤了19024个被少于3个细胞表达的基因

```
#计算线粒体基因占比
adata.var["mt"] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
#将线粒体基因比例可视化
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
```

![image](https://user-images.githubusercontent.com/112565216/188347355-c5588d35-b254-4f20-a03a-f018752879f0.png)

```
#将测序深度与线粒体基因占比和检测到的基因数的关系可视化
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
```

![image](https://user-images.githubusercontent.com/112565216/188347632-66e1a4fe-8bf4-4354-a048-fc8922803c48.png)
![image](https://user-images.githubusercontent.com/112565216/188347727-a86c8802-5c68-45d4-a05f-7e1ba0bea7fa.png)

从图中可看出检测到的线粒体基因的含量不随测序深度的升高而增大
```
#进一步筛选细胞，要求线粒体基因含量小于5%，基因数小于2500
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

```

### 4.降维
```
#寻找高变基因并可视化
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
```

![image](https://user-images.githubusercontent.com/112565216/188348202-4d2f272c-8af0-466e-8ac9-45db057e3340.png)

```
#将归一化后的数据进行保存
adata.raw = adata
#提取高变基因
adata = adata[:, adata.var.highly_variable]
#将测序深度和线粒体基因比例进行归一化处理
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
#对数据放缩变换，使数据方差为1，平均值为0
sc.pp.scale(adata, max_value=10)
#PCA降维
sc.tl.pca(adata, svd_solver='arpack')
#将某个基因的PCA降维结果可视化
sc.pl.pca(adata, color='CST3')
```
![image](https://user-images.githubusercontent.com/112565216/188348766-ae14985b-80ff-4107-ae66-8a90c07eb87e.png)

```
#将每个PC的贡献度可视化
sc.pl.pca_variance_ratio(adata, log=True)

```
![image](https://user-images.githubusercontent.com/112565216/188348839-735cf99b-da41-4f53-b616-0d8623f14038.png)

```
#进行UMAP降维
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
#进行聚类分簇
sc.tl.leiden(adata)
sc.pl.umap(adata, color='leiden')
```
![image](https://user-images.githubusercontent.com/112565216/188350479-22b237cf-f7c8-4eb6-ba80-e6e5fb8746ee.png)

### 5.标记cluster
确定marker基因（※）

根据marker基因确定cluster的细胞种类后对细胞进行标记
```
new_cluster_names = [
    'CD4 T', 'CD14 Monocytes',
    'B', 'CD8 T',
    'NK', 'FCGR3A Monocytes',
    'Dendritic', 'Megakaryocytes']
adata.rename_categories('leiden', new_cluster_names)
sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False)
```
![image](https://user-images.githubusercontent.com/112565216/188464735-d76909ec-8f02-441d-b853-c6eee1149193.png)

