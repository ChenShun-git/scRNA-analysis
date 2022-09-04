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
    "BC2",  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True) 
```

### 3.质控
```
#将表达量最高的20个基因可视化
sc.pl.highest_expr_genes(adata, n_top=20, )

```
![highest_expr_genes](https://user-images.githubusercontent.com/112565216/188308530-8a2cd3c3-bc78-40f0-a02c-32597928ed4a.png)

```
#初步过滤基因与细胞
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
```
![过滤结果](https://user-images.githubusercontent.com/112565216/188308593-5ef958ba-8a61-474c-b6d4-7306d1af0770.png)

过滤结果显示过滤了两个基因表达量少于200的细胞和16213个被少于3个细胞表达的基因

```
#计算线粒体基因占比
adata.var["mt"] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
#将线粒体基因比例可视化
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
```

![image](https://user-images.githubusercontent.com/112565216/188308631-d44f694b-fb56-4d76-965c-4076b38be6fb.png)

```
#将测序深度与线粒体基因占比和检测到的基因数的关系可视化
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
```

![image](https://user-images.githubusercontent.com/112565216/188308789-ef499ae4-56f8-4a27-a653-aa27226d590d.png)
![image](https://user-images.githubusercontent.com/112565216/188308810-16341b21-5d0b-4d8c-abd6-1954bc7bf11f.png)

从图中可看出检测到的线粒体基因的含量不随测序深度的升高而增大
```
#进一步筛选细胞，要求线粒体基因含量小于5%，基因数小于2500
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

```
